# nf-demultiplex-matfet

*matfet*: maternal/fetal. A Nextflow pipeline that tells maternal and fetal cells apart
by genotype in single-cell RNA-seq from pregnancies (placenta, decidua, blood, or
fallopian tube for ectopic pregnancies). It is a separate pipeline built on
[`nf-demultiplex`](https://github.com/cellgeni/nf-demultiplex), whose Souporcell steps it
reuses through a pinned git submodule.

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
cd /lustre/scratch127/cellgen/cellgeni/tickets/REQ-xxxxx
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

### Reference files

The defaults are Sanger farm paths. Override them elsewhere, and whenever the BAMs were
aligned to a different reference: Souporcell's remap and allele counts are only correct
against the reference the BAMs were aligned to.

| Parameter | Default | Description |
|---|---|---|
| `--soc_fasta` | `/nfs/cellgeni/STAR/human/2020A/GRCh38_v32_modified.fa` | Genome FASTA of the BAMs' reference (GRCh38, 10x 2020-A). |
| `--soc_vcf` | `/nfs/cellgeni/pipeline-files/demultiplex/reference_vcfs/chr_filtered_2p_1kgenomes_GRCh38.vcf` | Common-variant panel (1000 Genomes SNPs, allele frequency ≥ 2%); Souporcell genotypes cells at these sites. |
| `--container` | `/nfs/cellgeni/singularity/images/souporcell_v2.5-python_3.10-vireo_v0.5.8-cellsnp-lite_1.2.3-nextflow-compatible.sif` | Singularity image with Souporcell, samtools and Python; always used. |

### Samplesheet

CSV with a header; relative paths resolve against the samplesheet's directory.

```csv
group_id,sample_id,bam,barcodes,expression,souporcell_dir
CW5,Embryo-Imp15350362,/path/362/possorted_genome_bam.bam,/path/362/barcodes.tsv.gz,/path/362/filtered_feature_bc_matrix,
CW5,Embryo-Imp15350363,/path/363/possorted_genome_bam.bam,/path/363/barcodes.tsv.gz,/path/363/filtered_feature_bc_matrix,/path/results/souporcell/Embryo-Imp15350363
```

| Column | Required | Description |
|---|---|---|
| `group_id` | yes | Pregnancy group (one mother, one fetus). A group of one library is allowed. |
| `sample_id` | yes | Unique library identifier (letters, digits, `.`, `_`, `-`). |
| `bam` | yes | CellRanger BAM with cell barcodes (indexed, or the index is created). |
| `barcodes` | yes | Filtered cell-barcode file (optionally `.gz`). |
| `expression` | yes | CellRanger `filtered_feature_bc_matrix` directory. |
| `sample_role` | no | Default `mixture`; `maternal_reference` is not supported yet. |
| `expected_genotypes` | no | Default `2`; only 2 is supported. |
| `souporcell_dir` | no | An existing per-library Souporcell output (`souporcell/<sample_id>/` of an earlier run). The library's whole Souporcell step, remap and allele counting, is skipped and its allele counts reused. Must contain `alt.mtx`, `ref.mtx`, `barcodes.tsv` and `common_variants_covered.vcf`. |

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

## Outputs (under `--outdir`)

| Path | Description |
|---|---|
| `<sample_id>/cell_assignments.tsv.gz` | Per barcode: genotype cluster, maternal/fetal label, confidence, reason. |
| `<sample_id>/library_summary.tsv` | Per library: gate verdict and cells per label. |
| `reports/<group_id>.html` | Human-readable report per group; read it before using the labels. |
| `groups/<group_id>/` | Group matrix QC, Souporcell discovery, `gate.tsv`, marker summary, cluster labels. |
| `souporcell/<sample_id>/` | Per-library Souporcell run (allele counts; its own clusters are diagnostic only). |
| `pipeline_info/` | Nextflow run records (timeline, trace, report, DAG). |

Every column and warning is explained in the overview page, appendix B.

## Testing

```bash
python3 tests/test_composition_gate.py      # gate unit tests (Python 3 only)
python3 tests/test_decision_engine.py       # decision-engine unit tests (Python 3 only)
bash tests/test_souporcell_wrapper.sh       # wrapper argument regression test
python3 tests/make_fixture.py               # synthetic two-group fixture (needs numpy + scipy)
nextflow run . -profile test                # end-to-end on the fixture
```

The end-to-end test runs on the local executor in the container: run it inside an LSF
job with `module load cellgen/nextflow cellgen/singularity`, not on a head node. See
[`docs/testing.md`](docs/testing.md).

## Documentation

- [`docs/pipeline_design.md`](docs/pipeline_design.md): design, evidence and the rationale
  for every threshold.
- [`docs/pipeline_overview.html`](docs/pipeline_overview.html): illustrated walkthrough for
  a mixed lab audience (steps, spike-in tests, the gate, REQ-68805 results, how to read the
  outputs). A single self-contained file that can be sent on its own; after a new run,
  refresh its data with `python3 docs/build_overview_data.py`.
- [`docs/testing.md`](docs/testing.md): the test layers and the fixture's expected verdicts.

## Status

`v0.2.0-dev`. The group-aware path is implemented and validated: on REQ-68805 (64 libraries,
27 pregnancies) it reproduces the known CW5 split, rejects every no-fetus test case, and
identifies a fetus in 19 pregnancies (10 male at high confidence, 9 likely female,
provisional). Next: a fetal-tissue mode that flags maternal cells in embryo samples with no
maternal reference.
