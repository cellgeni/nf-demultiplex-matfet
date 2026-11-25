# Fetal mode (`--mode fetal`)

How fetal mode flags maternal cells in libraries from fetal tissue, where it is valid, how it
was validated, and how the pipeline implements it. Status: implemented
(2026-11-03); the scorer is the frozen prototype v3 and reproduces its results exactly
(see [Implementation checks](#implementation-checks)). Since 2026-11-16 fetal mode first runs
Souporcell k=2 on one library per batch and compares the clusters across libraries (the paper's
range check C′), as a reported check that never changes the calls
([Souporcell cross-library check](#souporcell-cross-library-check)).

`<dev>` below is the project's development folder, which is not public. The
benchmark and validation scripts live there and go to the HDCA analysis repository; this
repository holds the method and its regression tests.

## The problem

Single-cell libraries from fetal tissue can hold a few maternal cells (microchimerism, or
maternal blood or decidua picked up during collection). There is no maternal sample, and the
maternal cells are a tiny share. Pregnancy mode cannot handle that: Souporcell discovers two
genotypes only when both contribute many cells, and at a few maternal cells among thousands it
either misses the mother or splits the fetus.

Fetal mode does not discover two genotypes. It learns the fetus's genotype from the pooled
cells and asks, cell by cell, whether the reads fit it. Flagged cells are removed from the
atlas, so the method is tuned for **precision**: removing a fetal cell is worse than missing a
maternal one.

## When to use fetal mode and when pregnancy mode

The full guidance, with the benchmark numbers behind it, is
`<dev>/soc_ventotormo18/eval/choosing_a_method.md`. In short:

| True maternal share | Fetal mode | Pregnancy mode (Souporcell k=2) |
|---|---|---|
| < 2% | **Right tool.** Catches most maternal cells with ≥ 1,000 positions; no fetal cell lost | Loses 5–6% of fetal cells (mostly as pairs) and overestimates the share; usually misses the mother at 0.1% |
| 2–3% | Crossover: catches 65–89% at full depth | Catches 98–99%; loses 1–5% of fetal cells |
| > 3% | Declining: 55% caught at 4%, 26% at 6%, 3% at 8%; blind from about 8–10% (flags nothing) | **Right tool for removal.** Catches 96–100%; fetal losses fall to ≤ 1% from 5% |

- **Default for embryo and fetal-organ libraries: fetal mode.** It never lost a fetal cell in
  the benchmark.
- **Use pregnancy mode, or both,** for tissue at the mother–fetus boundary (placenta, decidua,
  membranes, cord, cord blood, amniotic fluid), for samples that may have picked up maternal
  blood or decidua, and when fetal mode flags many cells.
- **Pools of several embryos, or twins:** split by genotype first; fetal mode assumes one
  fetus per group.
- **Zero flags are not proof on their own.** Above about 10% fetal mode also reports zero, so
  the result needs in-range evidence (below). The Souporcell cross-library check, which fetal
  mode now runs by default, covers that band.

The two modes are separate runs of the pipeline; the practical rule in `choosing_a_method.md`
uses signals from both.

## Method

Per batch of up to six libraries of one fetus (`bin/fetal_score.py`):

1. **Allele counts.** `cellsnp-lite` counts reads at a fixed panel of common SNPs (the
   1000 Genomes panel the pipeline already uses, allele frequency ≥ 2%, which also supplies the
   population frequencies), straight from the CellRanger or STARsolo BAM, with no remapping.
2. **Fetal genotype.** All cells of the batch are pooled and each site is called hom-ref,
   het or hom-alt from the pooled allele fraction (≤ 0.03, 0.3–0.7, ≥ 0.97; at least 20 UMIs;
   autosomes only). This is the method's central assumption: the pool is the fetus.
3. **Site filters.** Fetal-het sites whose cells spread far wider than binomial (paralog
   pseudo-SNPs, imprinted or monoallelic genes) are dropped; sites within 10 bp are one event
   and only the deepest is kept.
4. **Read model.** Per cell and site, which alleles appear among the cell's UMIs (ref only,
   alt only, both), with the per-cell allele fraction drawn from a Beta around the expected
   value. Its concentration is fitted on the data, separately for fetal-het and homozygous
   sites, so a burst of one allele cannot run up a score.
5. **Maternal hypothesis.** The mother's genotype is unknown, but she shares one allele with
   the fetus at every site; her possible genotypes are weighted by the population frequency.
   A maternal cell may carry fetal background RNA: the expected fraction is mixed with a fetal
   share of 0, 0.1 or 0.25 (best kept). A fetal–maternal pair in one droplet is a share of 0.5.
6. **Score and cut-off.** Per cell, the log-likelihood ratio of maternal against fetal, and of
   pair against fetal. The cut-off is the larger of 10 and the quantile at which simulated fetal
   cells, at each cell's own coverage, are flagged at most 1 in 10,000 times.
7. **Status.**

   | `status` (scorer) | `fetal_mode` (outputs) | Meaning |
   |---|---|---|
   | `maternal` | `maternal` | Above the maternal cut-off. Remove. |
   | `fetal_maternal_doublet` | `maternal_in_pair` | Above the pair cut-off. Remove. |
   | `fetal_consistent` | `fetal` | ≥ 1,000 informative sites and not flagged. |
   | `insufficient_evidence` | `too_little_data` | Fewer sites. Kept, **not** confirmed fetal. |

8. **Bound, per fetus** (`bin/fetal_summary.py`). The one-sided 95% upper bound on the share
   of maternal cells among the cells cleared as fetal is **3 ÷ cleared cells**, the exact
   zero-event bound, pooled over all batches of the fetus. It replaces the scorer's own
   likelihood-ratio bound (`maternal_share_upper95` in the batch summaries, kept for
   comparison), which was too optimistic at zero events (decided 2026-10-29).

### Why 1,000 sites

The minimum number of informative sites a cell needs before it is cleared was chosen on real
maternal cells hidden in fetal-dominated libraries (D8, D9, D10 from REQ-68805): the smallest
of 200–1,600 at which at most 5% of hidden maternal cells above it were cleared as fetal in
every condition (full, thinned and with added fetal RNA). Worst condition by cut-off: 400
11.7%, 800 6.7%, **1,000 3.3%**. Leave-one-pregnancy-out showed 1,000 is the minimum, not a
margin: chosen on D9 and D10 alone the cut-off is 400, which lets 25% of D8's shallow maternal
cells through. It was then checked on the target kits (Hrv99 3′ v3; Hrv305 and CW16 5′ v3).

### Why batches of six

The six-library batch came from the pilot (six F147 libraries) and was kept for the final run
of all 144 embryo libraries (`<dev>/proto/final_batches.txt`), for practical reasons: one batch
needs 13–33 GB and 2–7 h. Batch size was never tested as a variable. The pipeline keeps six
because that is how the reported results were produced, and warns if `--fetal_batch_size` is
changed. The D8 benchmark scored one library per sample (6,400–10,600 cells), embryo batches
hold 14,000–240,000.

## Valid range and in-range checks

The pooled genotype is the fetus's only while maternal cells are a small share of the batch.
As the share grows, the sites where mother and fetus differ drift out of the genotype classes
and are dropped, which removes exactly the evidence the score needs. On D8 this is gradual
from 2% and complete by about 8–10%: nothing is flagged, the scorer's share estimate reads 0,
and the "at most" bound is wrong. A result is therefore only as good as the evidence that the
fetus is in range. The pipeline reports it (`fetal_summary.py`, `sex_check.tsv`,
`warnings.txt` and the report):

| Check | What it does | Warning |
|---|---|---|
| Sex call | Share of cells with Y-gene reads: male ≥ 50%, female ≤ 5%, otherwise unclear | Unclear sex: the group may mix individuals |
| XIST/Y check (male fetus) | A maternal cell expresses XIST and no Y gene. Cells with XIST ≥ 2 UMIs and no Y read are counted without genotypes; with the chance that a female cell shows XIST ≥ 2 at each depth, (candidates + 3) ÷ (chance × cells) bounds the maternal share at **any** share | Bound ≥ `--fetal_max_maternal_share` (2%) |
| Female fetus | No genotype-free check exists; stated in the outputs | – |
| Flagged share per batch | Many flags mean the batch is outside the safe range | ≥ `--fetal_warn_flagged_fraction` (2%) |
| Souporcell cross-library check | Souporcell k=2 on one library per batch; genotype distance between every pair of clusters; a mother cluster stands apart from the fetus's other libraries ([below](#souporcell-cross-library-check)) | A cluster's median distance to the other libraries ≥ `--fetal_souporcell_warn_msd` (0.03); or only one library checked |

Caveats:
- **XIST detection rates** come from the female embryo F158 (5′ v2), by total UMIs per cell
  (0.07 below 2,000 UMIs to 0.87 above 20,000). For other kits they are an approximation.
- **The flag count is not a share.** It rises to about 2.3% of cells at a true 3–4% maternal
  share and then falls (1.6% at 6%, 0.25% at 8%, about 0 at 10%). A high count shows a batch is
  out of range, not by how much; a low count alone does not show it is in range. The 2% line is
  the one `choosing_a_method.md` uses, tested at full depth only.

### Souporcell cross-library check

**Why.** Fetal mode goes blind above about 8–10% maternal, exactly where Souporcell works. So
the pipeline runs Souporcell (k = 2) **before** fetal mode can be trusted, on the first library
of each batch (`--fetal_souporcell per_batch`, the default; `all` or `off` also exist). In the
paper this was range check C/C′ on one library per batch of each embryo.

**Reading k = 2 without a mother.** Souporcell always returns two clusters. With no mother both
are the fetus, split along a technical axis. A single library cannot tell the two cases apart
(the composition gate and the genotype-discordance gap both gave false "mother found" calls:
F218 library 1754, F137 library 2348). Comparing clusters **across libraries of the same fetus**
can: every fetal cluster is the same person, so a mother cluster is the odd one out.

**Distance.** For two clusters, the mean squared difference of the reference-allele fraction
per site, over sites where every cluster of both runs has ≥ 8 reads (singlets only). It is the
team's `shared_samples.py` loss divided by the number of sites, so pairs with different site
counts can be compared. Code: `bin/souporcell_distance.py`, a port of
`<dev>/followups/loss_heatmap/compute.py`, which drew the paper's heatmaps.

| Pair of clusters | Distance per site |
|---|---|
| Same person, benchmark (fetus split in two; mother across replicates) | 0.0000–0.0065 |
| Any two clusters within one embryo (F137, F147, F158, F218; 30 runs) | 0.0048–0.0169 |
| Mother vs her own fetus (D8, D9, D10) | **0.041–0.047** |
| Different pregnancies | 0.068–0.089 |

**Warning.** For each cluster with ≥ 20 singlets (`--fetal_souporcell_min_singlets`): the median
distance to the clusters of the fetus's other checked libraries. At ≥ 0.03
(`--fetal_souporcell_warn_msd`, between the two bold-face levels above) the cluster is marked
`review` and the report warns. It is not a decision: a flagged
cluster gets the review F218's library 1754 got (cell types, cell sizes, Y-positive share,
fetal-mode status).

**Outputs.** `fetal/souporcell_check/`: `libraries.tsv` (per cluster: singlets, distance to its
sister cluster and to the other libraries, status), `pairs.tsv` (every pair), `groups.tsv`
(fetus × fetus median), `warnings.tsv` and SVG heatmaps. Each fetus's report shows its cluster
heatmap (rows and columns = each checked library's clusters) and, with two fetuses or more,
the fetus × fetus panel, as in the paper figure but without the benchmark rows.

**Limits.**
- **Sensitivity:** Souporcell gives a mother her own cluster from about 300 maternal cells
  (`<dev>/soc_ventotormo18_rep/eval/report.txt`), 1–7% of an embryo library. Below that, fetal
  mode's own flags and the flag count cover the range.
- **Coverage:** one library per batch. A mother at 30–50% in several unchecked libraries of a
  female fetus's batch is not covered (`<dev>/paper/range_checks.md` §6); use `all` when that
  matters.
- **One checked library** (a fetus with one batch): nothing to compare against; the report says
  so. Use `--fetal_souporcell all`.
- **Library choice:** the paper picked the library with the most flags in each batch, or the
  largest; the pipeline takes the first, because Souporcell runs before fetal mode's flags exist.

## Validation

### Development tests (2026-10-08 to 10-19)

Prototype versions v0–v3 were tested on real maternal cells hidden in fetal-dominated
libraries of REQ-68805 pregnancies (D8, D9, D10; then Hrv99, Hrv305, CW16 for the target kits)
and on the F147 Y-gene answer key. Every version answered a failure found in testing:

| Version | What changed | Why |
|---|---|---|
| v0 | One library; fetal fraction from that library's pool; mother shares one allele with the fetus; standardised score | First check on a real library |
| v1 | Genotype pooled over libraries; raw log-likelihood ratio with a simulated-null cut-off (1 in 10,000, floor 10); pair hypothesis; four statuses | One library is too little; a z-score has no clear cut-off |
| v2 | Over-dispersed het sites dropped; adjacent sites merged; one UMI per allele | v1 flagged erythroid "pairs" driven by four adjacent HBG1 positions (a paralog pseudo-SNP) |
| v3 | Beta read model; fetal background share in maternal cells; mother's background from the pool | Real maternal cells scored far less maternal than simulated ones and were caught mostly as pairs |
| **v3 + 1,000 sites (frozen 2026-10-19)** | Cleared only with ≥ 1,000 informative sites | Simulated detection power over-promised below a few hundred sites |
| v3.1 (dropped) | Fetal background share 0.4 added | No gain on maternal cells with 30% added fetal RNA |

At the frozen settings, on D8–D10 at full depth: 97.0% of hidden maternal cells caught (95% CI
94.4–98.4), 0% cleared as fetal; no fetal cell flagged in any set; F147 Y-positive cells
flagged 0 of 57,471. Details and scripts: `<dev>/spikein/`, `<dev>/proto/`.

### Public benchmark (2026-10-30 to 11-03)

Vento-Tormo et al. 2018 (E-MTAB-6701), donor D8: the placenta library (mostly fetal, male
fetus) plus two decidua libraries and one blood library (maternal). All 6,397 fetal placenta
cells with maternal cells added at a known share (0.1–40%), at full, half and quarter depth,
228 samples plus 90 follow-up samples at 3–8%. Labels were checked independently from each
cell's own reads (0 of 15,783 dropped). Pass marks were fixed before any score was read.

| Mark | Verdict |
|---|---|
| No fetal cell flagged | **PASS**: 0 in 1,458,516 cell scorings, and 0 in the follow-up |
| ≥ 95% of maternal cells with ≥ 1,000 positions caught | FAIL pooled; 94–100% at ≤ 2% maternal |
| ≥ 80% of caught cells called maternal, not pair | FAIL (38–82%); pairs are removed too |
| ≤ 5% of maternal cells cleared as fetal, every condition | FAIL; at most 4.7% at ≤ 2% |
| The "at most" bound covers the truth in ≥ 95% of samples | FAIL (180 of 225); holds in every sample up to 5% |

The failures all come from the shares above the valid range. Working range (all marks hold):
1–2% at full depth, up to 3–4% at quarter depth. Against Souporcell on the same samples, the
methods cross at about 3% true maternal share. Results, tables and caveats:
`<dev>/soc_ventotormo18/eval/results.md`.

### Embryo results

All 188 libraries of four embryos: F147, F158 and F218 in 24 batches of six (2026-10-19), F137
in 8 batches of five or six (2026-11-11, `<dev>/proto/final_batches_f137.txt`):

| Embryo | Kit | Cells | Cleared as fetal | Too little data | Maternal | Pairs | Bound per 100,000 cleared cells |
|---|---|---|---|---|---|---|---|
| F137 (female) | 5′ (v1/v2) | 315,335 | 244,696 | 70,600 | 39 | 0 | ≤ 1.23 |
| F147 (male) | 5′ v2 | 496,173 | 362,678 | 133,495 | 0 | 0 | ≤ 0.83 |
| F158 (female) | 5′ v2 | 444,681 | 302,581 | 142,095 | 4 | 1 | ≤ 0.99 |
| F218 (female) | 3′ v3 | 1,694,857 | 391,773 | 1,303,076 | 5 | 3 | ≤ 0.77 |

- **F218's bound:** 3 ÷ 391,773 = 7.66 per million, which rounds to **0.77** per 100,000.
  Earlier documents (the handoff and the benchmark results) give 0.76.
- **Every flag was checked** (relatedness, mtDNA, atlas QC; `<dev>/proto/flag_summary.tsv`).
  All 8 F218 flags fit the mother: **8 candidate maternal cells**, about 5 per million. None of
  F158's 5 is the mother: three are one barcode seen in three libraries that does not fit her,
  one is unclear, and the pair fits the fetus; all are tiny droplets. F137's 39 flags are 3
  antibody-only barcodes recurring in 41–44 of its libraries; pooled across libraries, all
  carry F137's own genotype (`<dev>/f137/barcode_review.tsv`): 0 maternal. After review:
  F137 0, F147 0, F158 0, F218 8.
- **In range, all four** (`<dev>/paper/range_checks.md`, `<dev>/f137/PLAN.md`): the flag count
  (A), the blind-zone statistic (B, all 32 batches silent), Souporcell k=2 across libraries
  (C: gap, pairing ratio, and the distance per site above, C′) are silent in every embryo; the XIST/Y check (D) applies to F147 only (0 candidates, maternal share
  ≤ 1.32 per 100,000 at any share). No cell has been removed so far.
- **F218 is mostly too little data:** 77% of its cells have fewer than 1,000 informative
  sites, so its bound covers only the other 23%.

### Implementation checks

The pipeline was checked against the prototype on real data (`<dev>/pipeline_port/`):

| Check | Result |
|---|---|
| Scorer (`bin/fetal_score.py`) on F147 batch 5, F218 batch 2 and three D8 samples (clean, 2%, 40%), in the container | Identical: cells, order, statuses, sites, Y counts, removed sites; scores within 3 parts in 10¹⁵ |
| End to end, F147 batch 5 from the STARsolo BAMs through Nextflow | All six pileups byte-identical; scores identical |
| Summary on the 24 batches of F147, F158, F218 | Every number in the table above for those three, their 13 flags, and F147's XIST/Y bound reproduced; no warnings |
| Souporcell check (`bin/souporcell_distance.py`) on the paper's 30 embryo runs, against `loss_heatmap/pairs.tsv` (`tests/regression_souporcell_distance.py`, `<dev>/pipeline_port/souporcell_distance/`) | Identical: all 1,770 cluster pairs, same sites matched and used, msd difference 0 (job 674883, 2026-11-16). All 60 clusters `ok`, 0 warnings; medians to other libraries 0.0012–0.0083 (rule: ≥ 0.03). F137 library 2348, a false "mother found" under the old gap rule, reads 0.0054 |

**Not checked through the pipeline:** F137 was screened later with the frozen prototype only
(`<dev>/proto/out/final/final_f137_b*`); the pipeline has not been run on it and will not be.
Its batches were 4 of six and 4 of five; the pipeline cuts batches of `--fetal_batch_size` in
samplesheet order, so it would give 7 batches of six and one of two for F137's 44 libraries,
not the same batches. The Souporcell check has not run inside Nextflow end to end; its DAG
was checked with `-preview` on real F218 inputs.

`tests/regression_fetal_score.py` is the comparison used for the scorer. The port changed only
the interface: arguments instead of hard-coded paths, Y-gene counts from each library's
expression matrix, no atlas cell-group column, and one development-only printout removed.

## Pipeline

```
SOUPORCELL k=2 (first library of each batch)  ->  SOUPORCELL_DISTANCE (whole run) --+
CELLSNP_PILEUP (per library)  ->  FETAL_SCORE (per batch of up to 6 libraries)      |
SEX_COUNTS (per library)      ->  FETAL_SUMMARY (per fetus)  ->  FETAL_REPORT (per fetus)
```

Souporcell and the pileups start together; the reports wait for the Souporcell check. The
check never feeds into the scores.

- **Samplesheet:** the same columns as pregnancy mode, with `group_id` meaning one fetus.
  `cellsnp_dir` (optional) reuses an existing pileup; `souporcell_dir` (optional) reuses an
  existing Souporcell run for the check (`alt.mtx`, `ref.mtx`, `clusters.tsv`,
  `common_variants_covered.vcf`), and that library is checked even if it would not be selected;
  `sample_role` and `expected_genotypes` are refused. Batches are each fetus's libraries in
  samplesheet order.
- **Settings:**

  | Parameter | Default | Note |
  |---|---|---|
  | `--fetal_batch_size` | 6 | Validated value; warned if changed |
  | `--fetal_min_sites` | 1000 | Validated value; warned if changed |
  | `--fetal_rhos` | `0,0.1,0.25` | Validated value; warned if changed |
  | `--fetal_warn_flagged_fraction` | 0.02 | In-range warning |
  | `--fetal_max_maternal_share` | 0.02 | In-range warning (male XIST/Y bound) |
  | `--fetal_souporcell` | `per_batch` | Libraries for the Souporcell check: `per_batch`, `all` or `off` |
  | `--fetal_souporcell_warn_msd` | 0.03 | Souporcell check warning |
  | `--fetal_souporcell_min_singlets` | 20 | Smaller clusters are shown, not judged |
  | `--soc_vcf` | 1000 Genomes panel | Sites to pile up and population frequencies |

- **Outputs:**

  | Path | Contents |
  |---|---|
  | `<sample_id>/fetal_mode_cells.tsv.gz` | Per library: barcode, `fetal_mode`, sites, scores |
  | `fetal/<group_id>/fetal_summary.tsv` | The fetus and each batch: counts, flagged share, bound |
  | `fetal/<group_id>/sex_check.tsv` | Sex call and the XIST/Y check |
  | `fetal/<group_id>/warnings.txt` | In-range warnings, one per line |
  | `fetal/<group_id>/cells.tsv.gz`, `flagged_cells.tsv` | Every cell with its evidence; the removal candidates |
  | `fetal/<group_id>/batches/` | The scorer's raw outputs per batch |
  | `fetal/souporcell_check/` | Souporcell check: `libraries.tsv`, `pairs.tsv`, `groups.tsv`, `warnings.tsv`, `<group_id>.heatmap.svg`, `groups.heatmap.svg` |
  | `souporcell/<sample_id>/` | Souporcell output of each checked library |
  | `reports/<group_id>.html` | Result, bound, valid range, warnings, in-range evidence with the Souporcell heatmaps, batches, flagged cells |

- **Resources** (from the prototype runs): pileup 16 CPUs, 1–1.3 GB, minutes to about 5 h per
  library; scoring 48 GB and 12 h per batch, doubled on retry. Souporcell 16 CPUs, 35 GB, 24 h
  per library (the 30 embryo runs peaked at 21.6 GB and 17.7 h); the distance step 8 GB (30
  runs: 2.7 GB, 5 min).

## Limitations

- **Valid only for mostly fetal batches** (above). Fetal tissue normally is; the in-range
  evidence is weaker for female fetuses.
- **Shallow cells are not cleared.** In 3′ v3 data most cells stay `too_little_data`.
- **Recall at low depth is limited by information**, not the model: maternal cells with 100–200
  informative sites were caught 24% of the time in development.
- **Cells from another sample** (cross-contamination, index hopping) can be flagged; they are
  not the mother. The checks on the embryo flags separated these by relatedness, which the
  pipeline does not yet run.
- **The batch size** was not varied.

## Possible next steps

- Switching a fetus to pregnancy mode automatically when the checks fire
  (`<dev>/PLAN_followups_2026-11-04.md` §5). Designed, not built; the Souporcell check stays a
  report until a person has reviewed a real positive.
- A relatedness check on flagged cells in the pipeline, to separate the mother from other
  contaminants (prototype: `<dev>/proto/flag_check.py`), and a `flagged_not_maternal` value.
- XIST detection rates measured on 3′ kits.
