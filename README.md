# nf-demultiplex-matfet

*matfet*: maternal/fetal. A Nextflow pipeline that tells maternal and fetal cells apart
by genotype in single-cell RNA-seq. It is a separate pipeline built on
[`nf-demultiplex`](https://github.com/cellgeni/nf-demultiplex), whose Souporcell steps it
reuses through a pinned git submodule. It has two modes:

| Mode | For | How |
|---|---|---|
| `--mode pregnancy` (default) | Pregnancy samples where mother and fetus both contribute many cells: placenta, decidua, blood, or fallopian tube for ectopic pregnancies | Souporcell discovers both genotypes per pregnancy group; a composition gate and a marker panel name them |
| `--mode fetal` | Fetal tissue (embryo, fetal organs) where maternal cells are rare and there is no maternal sample | The fetus's genotype is learned from the pooled cells and every cell is scored against it; see [Fetal mode](#fetal-mode) |

Which to use depends on the expected maternal share: fetal mode is the right tool below
about 2% maternal cells and goes blind above about 8-10%; pregnancy mode is the right tool
from about 3-5%. See [`docs/fetal_mode_design.md`](docs/fetal_mode_design.md#when-to-use-fetal-mode-and-when-pregnancy-mode).

## Pregnancy mode

Every library holds cells of at most two individuals, mother and fetus, and libraries
from the same pregnancy share both. Genotypes are therefore inferred **once per pregnancy
group**, not per library:

1. **Per library**: Souporcell remapping and allele counting.
2. **Per group**: the libraries' allele matrices are merged, Souporcell discovers the two
   genotypes on the whole group, and a **composition gate** decides whether the second
   genotype is a real second person.
3. **Per group**: a pregnancy marker panel orients the two genotypes as maternal / fetal
   via an auditable two-tier decision engine, and the labels are split back into libraries.

Souporcell separates cells *by genotype*, not expression; the marker panel only names the
two genotypes. The pipeline **preserves uncertainty**: a group with no detectable second
genotype, or with weak or conflicting marker evidence, is reported `unresolved` rather than
forced into a call.

## Quick start

Clone with submodules (the Souporcell steps live in the `nf-demultiplex` submodule):

```bash
git clone --recurse-submodules https://github.com/cellgeni/nf-demultiplex-matfet.git
# or, in an existing clone: git submodule update --init
```

Run from the directory where results should go, e.g. a ticket directory. Nextflow submits
each task to LSF; run the launcher itself in a persistent session (`screen`/`tmux`).

```bash
cd /path/to/run_dir
module load cellgen/nextflow

nextflow run /path/to/nf-demultiplex-matfet/main.nf \
    --input actions/samplesheet.csv \
    --outdir results \
    -profile sanger \
    -resume
```

Task files go to `nf-work/` in the launch directory (`-w` overrides it); keep it for
`-resume` until the results are checked. The per-library Souporcell step is the heavy one:
16 CPUs and 35 GB per library, 2.4-11.4 h (median 5.6 h) on REQ-68805's 64 libraries (24 h limit), and about 10 GB of re-aligned reads per library
under `souporcell/`. Reuse an earlier run's output with `souporcell_dir` to skip it.
This runs pregnancy mode; for fetal tissue add `--mode fetal` (see [Fetal mode](#fetal-mode)).

### Reference files

These have no default: give them on the command line, or use `-profile sanger`, which sets
the Sanger farm copies (`conf/sanger.config`). Use a FASTA of the reference the BAMs were
aligned to: Souporcell's remap and allele counts are only correct against it.

| Parameter | Default | Description |
|---|---|---|
| `--soc_fasta` | none; `-profile sanger` sets it | Pregnancy mode, and fetal mode's Souporcell check. Genome FASTA of the BAMs' reference (GRCh38, 10x 2020-A). |
| `--soc_vcf` | none; `-profile sanger` sets it | Both modes. Common-variant panel (1000 Genomes SNPs, allele frequency ≥ 2%): Souporcell genotypes cells at these sites; fetal mode piles them up and takes population frequencies from it. |
| `--container` | none; `-profile sanger` sets it | Singularity image with Souporcell, cellsnp-lite, samtools and Python; always used. |

Off the farm, build the image from the repository's [`Dockerfile`](Dockerfile), which
reproduces the farm image's tool versions, convert it with `singularity build`, and pass the
`.sif` as `--container`.

### Samplesheet

CSV with a header; relative paths resolve against the samplesheet's directory.

```csv
group_id,sample_id,bam,barcodes,expression,souporcell_dir
CW5,Embryo-Imp15350362,/path/362/possorted_genome_bam.bam,/path/362/barcodes.tsv.gz,/path/362/filtered_feature_bc_matrix,
CW5,Embryo-Imp15350363,/path/363/possorted_genome_bam.bam,/path/363/barcodes.tsv.gz,/path/363/filtered_feature_bc_matrix,/path/results/souporcell/Embryo-Imp15350363
```

| Column | Required | Description |
|---|---|---|
| `group_id` | yes | Pregnancy mode: pregnancy group (one mother, one fetus). Fetal mode: one fetus. A group of one library is allowed. |
| `sample_id` | yes | Unique library identifier (letters, digits, `.`, `_`, `-`). |
| `bam` | yes | CellRanger or STARsolo BAM with `CB`/`UB` tags (indexed, or the index is created). |
| `barcodes` | yes | Filtered cell-barcode file (optionally `.gz`). |
| `expression` | yes | Filtered expression matrix directory: CellRanger `filtered_feature_bc_matrix` or STARsolo `GeneFull/filtered`. |
| `sample_role` | no | Pregnancy mode only. Default `mixture`; `maternal_reference` is not supported yet. |
| `expected_genotypes` | no | Pregnancy mode only. Default `2`; only 2 is supported. |
| `souporcell_dir` | no | An existing per-library Souporcell output (`souporcell/<sample_id>/` of an earlier run). Pregnancy mode: the library's whole Souporcell step, remap and allele counting, is skipped and its allele counts reused; must contain `alt.mtx`, `ref.mtx`, `barcodes.tsv` and `common_variants_covered.vcf`. Fetal mode: the run joins the Souporcell cross-library check without a new run, even if the library would not be selected; must contain `alt.mtx`, `ref.mtx`, `clusters.tsv` and `common_variants_covered.vcf`. |
| `cellsnp_dir` | no | Fetal mode only. An existing cellsnp-lite output for the library, made with the same panel and flags; the pileup is skipped. Must contain `cellSNP.base.vcf.gz`, `cellSNP.tag.DP.mtx`, `cellSNP.tag.AD.mtx` and `cellSNP.samples.tsv`. |

`sample_role` and `expected_genotypes` are refused in fetal mode, and `cellsnp_dir` in
pregnancy mode.

## Composition gate

| Verdict | Rule (defaults) | Outcome |
|---|---|---|
| `two_profiles` | the two genotype profiles differ at ≥ 0.40 of jointly called loci and the minority cluster has ≥ 20 singlets | maternal/fetal orientation runs |
| `one_profile` | profiles differ, but minority < 20 singlets | no labels |
| `unresolved` | profiles differ at < 0.40 of loci (two halves of one person), or share no called loci | no labels |

In validation, fake splits differed at 0.17-0.30 of loci and real mother/fetus splits
at 0.53-0.66. See [`docs/pipeline_design.md`](docs/pipeline_design.md#composition-gate).

## The decision engine (two tiers)

| Tier | Rule | Result |
|---|---|---|
| 1 | **RPS4Y1** fraction-positive ≥ floor **and** ≥ N-fold enriched in one cluster | `fetal_high_confidence_male` (HIGH, male) |
| 2 | RPS4Y1 uninformative; ≥ `min_vote_margin` secondary markers enriched in one cluster | `fetal_marker_supported` (PROVISIONAL, likely female) |
| 3 | Insufficient / conflicting evidence, or the gate did not accept two profiles | `unresolved` |

*Absence* of RPS4Y1 is never treated as maternal evidence; a cluster is `maternal_likely`
only once the other is positively identified as fetal. Panel: `RPS4Y1` (Y, decisive) +
`HAPLN1, PITX2, GAPLINC, HPGDS, XAGE2, XAGE3, HLA-G` (secondary), in
[`assets/marker_panel.tsv`](assets/marker_panel.tsv). The engine's thresholds
(`rps4y1_*`, `min_marker_*`, `min_vote_margin`, `min_cluster_cells` in
`nextflow.config`) are starting points that have not been calibrated against known
answers, which is one reason Tier 2 calls are provisional.

## Pregnancy-mode outputs (under `--outdir`)

| Path | Description |
|---|---|
| `<sample_id>/cell_assignments.tsv.gz` | Per barcode: genotype cluster, maternal/fetal label, confidence, reason. |
| `<sample_id>/library_summary.tsv` | Per library: gate verdict and cells per label. |
| `reports/<group_id>.html` | Human-readable report per group; read it before using the labels. |
| `groups/<group_id>/` | Group matrix QC, Souporcell discovery, `gate.tsv`, marker summary, cluster labels. |
| `souporcell/<sample_id>/` | Per-library Souporcell run (allele counts; its own clusters are diagnostic only). |
| `pipeline_info/` | Nextflow run records (timeline, trace, report, DAG). |

Every column and warning is explained in the overview page, appendix B.

## Fetal mode

```bash
nextflow run /path/to/nf-demultiplex-matfet/main.nf \
    --mode fetal \
    --input actions/samplesheet.csv \
    --outdir results \
    -profile sanger \
    -resume
```

One `group_id` per fetus. First, Souporcell (k = 2) runs on the first library of each batch
(`--fetal_souporcell per_batch`; `all` or `off`), alongside the pileups. Each library is piled up with `cellsnp-lite` at the common-variant
panel (no remapping; minutes to about 5 h per library, 16 CPUs). Each fetus's libraries, in
samplesheet order, are then scored in batches of six (`--fetal_batch_size`): the pooled
genotype of a batch is taken as the fetus's, and each cell is scored for whether its reads fit
the fetus, a related mother, or a fetal-maternal pair (up to 48 GB and 12 h per batch).

| `fetal_mode` | Meaning |
|---|---|
| `maternal` | Flagged as a maternal cell. Remove. |
| `maternal_in_pair` | Flagged as a fetal-maternal pair. Remove. |
| `fetal` | At least 1,000 informative sites and not flagged. |
| `too_little_data` | Fewer sites: kept, not confirmed fetal. |

Per fetus, the one-sided 95% upper bound on maternal cells among the cells cleared as fetal is
3 ÷ cleared cells, pooled over the fetus's batches. **The result is valid only if maternal cells
are a small share of each batch** (up to about 2% in the benchmark); above about 8-10% fetal mode
flags nothing and the bound is wrong. The pipeline therefore reports in-range evidence: the sex
call, the genotype-free XIST/Y check for male fetuses, and warnings when a batch flags 2% or more
of its cells or the XIST/Y bound reaches 2%. Female fetuses have no genotype-free check.

**Souporcell cross-library check.** It covers the shares where fetal mode is blind. Each checked
library splits into two Souporcell clusters even without a mother; the pipeline measures the
genotype distance per SNP between every pair of clusters of the fetus (and of all fetuses in the
run) and draws heatmaps in the report. The fetus's own clusters sit at ≤ 0.017; a mother sits at
about 0.04 from her fetus's clusters (0.041–0.047 in the benchmark), unrelated people at
0.068–0.089. A cluster at ≥ 0.03 (`--fetal_souporcell_warn_msd`) from the fetus's other
libraries is marked for review and warned about. The check never changes the calls.

The batch size, the 1,000-site minimum and the fetal background shares (`--fetal_rhos`) are the
validated values; changing them is allowed but warned about.

| Path | Description |
|---|---|
| `<sample_id>/fetal_mode_cells.tsv.gz` | Per library: barcode, `fetal_mode`, informative sites, scores. |
| `reports/<group_id>.html` | Per fetus: result, bound, valid range, warnings, in-range evidence, batches, flagged cells. Read it before using the calls. |
| `fetal/<group_id>/fetal_summary.tsv` | The fetus and each batch: cell counts per `fetal_mode`, flagged share, bound. |
| `fetal/<group_id>/sex_check.tsv`, `warnings.txt` | Sex call, XIST/Y check, in-range warnings. |
| `fetal/<group_id>/cells.tsv.gz`, `flagged_cells.tsv` | Every cell with its evidence; the removal candidates. |
| `fetal/<group_id>/batches/` | The scorer's raw outputs per batch. |
| `fetal/souporcell_check/` | Souporcell check for the whole run: per-cluster distances and status (`libraries.tsv`), every pair (`pairs.tsv`), fetus × fetus (`groups.tsv`), `warnings.tsv`, SVG heatmaps. |
| `souporcell/<sample_id>/` | Souporcell output of each checked library. |

Method, validation (development tests, a public benchmark on Vento-Tormo 2018 donor D8, four
embryos) and limitations: [`docs/fetal_mode_design.md`](docs/fetal_mode_design.md).

## Testing

```bash
python3 tests/test_composition_gate.py      # gate unit tests (Python 3 only)
python3 tests/test_decision_engine.py       # decision-engine unit tests (Python 3 only)
bash tests/test_souporcell_wrapper.sh       # wrapper argument regression test
python3 tests/test_fetal_summary.py         # fetal-mode summary unit tests (needs pandas)
python3 tests/test_souporcell_distance.py   # fetal-mode Souporcell check unit tests (needs pandas)
python3 tests/make_fixture.py               # synthetic two-group fixture (needs numpy + scipy)
nextflow run . -profile sanger,test         # end-to-end on the fixture (or -profile test --container <image>)
```

The end-to-end test runs on the local executor in the container: run it inside an LSF
job with `module load cellgen/nextflow cellgen/singularity`, not on a head node. See
[`docs/testing.md`](docs/testing.md).

## Documentation

- [`docs/pipeline_design.md`](docs/pipeline_design.md): pregnancy-mode design, evidence and
  the rationale for every threshold.
- [`docs/fetal_mode_design.md`](docs/fetal_mode_design.md): fetal mode: method, when to use
  it, valid range, validation and outputs.
- [`docs/pipeline_overview.html`](docs/pipeline_overview.html): illustrated walkthrough for
  a mixed lab audience (steps, spike-in tests, the gate, REQ-68805 results, how to read the
  outputs). A single self-contained file that can be sent on its own; after a new run,
  refresh its data with `python3 docs/build_overview_data.py`.
- [`docs/testing.md`](docs/testing.md): the test layers and the fixture's expected verdicts.

## Status

`v0.2.0-dev`. The group-aware path is implemented and validated: on REQ-68805 (64 libraries,
27 pregnancies) it reproduces the known CW5 split, rejects every no-fetus test case, and
identifies a fetus in 19 pregnancies (10 male at high confidence, 9 likely female,
provisional). Fetal mode (`--mode fetal`) is implemented
and reproduces exactly the prototype that screened F147, F158 and F218 (2.6 million cells).
The fourth embryo, F137 (315,335 cells), was screened with the frozen prototype only, not
through the pipeline. The Souporcell cross-library check (added 2026-11-16) ports the paper's
distance heatmap and reproduces it exactly on the 30 embryo Souporcell runs, with no warnings.
