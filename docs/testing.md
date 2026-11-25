# Testing

Four layers for pregnancy mode, from fastest to most integrated, then the fetal-mode tests.
The pregnancy-mode acceptance tests on real data are listed in
[`pipeline_design.md`](pipeline_design.md#acceptance-tests).

## 1. Unit tests (Python 3 only)

```bash
python3 tests/test_composition_gate.py    # gate verdicts on every validation case
python3 tests/test_decision_engine.py     # two-tier engine + gate override
# or, if pytest is available:
pytest tests/
```

`test_composition_gate.py` feeds the gate the separating loci, jointly called loci
and minority sizes of real runs (CW5, the spike-in curve, the highest validation
fake split, CW4 `trio_early` and its collapsed subsets, and the production groups
the threshold was checked on: CW15, Hrv43, D8, Hrv46) and checks each lands on its
documented verdict, plus the boundaries and overrides.

`test_decision_engine.py` covers the decision engine's invariants:
- RPS4Y1 identifies the fetal cluster (either cluster id) at HIGH confidence, and takes
  precedence over conflicting secondary markers.
- **Absent RPS4Y1 never yields a maternal label**; with no other evidence the result is
  `unresolved`.
- Female pattern (secondary markers only) → PROVISIONAL, likely female.
- Vote margins, ambiguous expression, single-cluster and small-cluster QC gates,
  determinism.
- Any gate verdict other than `two_profiles` withholds every label, even a clear RPS4Y1 call.

## 2. Souporcell wrapper regression test

Uses a fake `souporcell_pipeline.py`, so it neither reads a BAM nor runs Souporcell:

```bash
bash tests/test_souporcell_wrapper.sh
```

## 3. Synthetic two-group fixture

`tests/make_fixture.py` (numpy + scipy) writes `test_data/fixture/`: allele matrices in
Vartrix layout with a real genotype signal, expression matrices, a reference-panel
header and a samplesheet whose rows all give `souporcell_dir` (so no BAM is read).

| Group | Libraries | Truth | Expected |
|---|---|---|---|
| `G1` | `libA` (180 maternal + 70 fetal), `libB` (200 maternal + 8 fetal) | one mother, one male fetus | gate `two_profiles`; RPS4Y1 identifies the fetal cluster; `libB` gets its handful of fetal cells |
| `G2` | `libC` (240 cells) | one genotype | gate `unresolved` |

Each library has a `truth.tsv` (barcode → maternal / fetal) for checking assignments.

```bash
python3 tests/make_fixture.py              # or: --out <dir> to write elsewhere
```

The scripts can be driven directly on it (`build_group_matrix.py` →
`check_group_matrix.py` → ... → `make_report.py`), following the commands in
`modules/local/*.nf`.

## 4. End-to-end on the fixture (Nextflow + container)

```bash
module load cellgen/nextflow cellgen/singularity
nextflow run . -profile sanger,test   # elsewhere: -profile test --container <image>
```

Runs everything from `BUILD_GROUP_MATRIX` onwards for real, in the pinned container,
on the local executor; outputs land in `results_test/`. Run it inside an LSF job, not
on a head node. Souporcell's behaviour on a 458-cell synthetic group is not
validation evidence: the fixture checks that the pipeline runs and wires correctly,
and the verdicts above are what a working pipeline should give, but a different
outcome on the fixture should be investigated rather than assumed to be a bug.

## Real data

Acceptance runs on the farm use `-profile sanger` with a real samplesheet
(`test_data/samplesheet_test.csv` is a template). For REQ-68805, `souporcell_dir`
can point at the existing `results/souporcell/<sample_id>` directories, which skips
remapping.

## Fetal mode

Fetal mode is checked against the frozen prototype that produced its results
(`<dev>/proto/score_v3.py`, `<dev>` being the development folder; see
[`fetal_mode_design.md`](fetal_mode_design.md#implementation-checks)). It has no synthetic
fixture: its scoring needs real genotype structure, and its evidence is that the pipeline
reproduces the prototype exactly on real data.

### Unit tests

```bash
python3 tests/test_fetal_summary.py       # needs pandas; or run inside the pipeline container
python3 tests/test_souporcell_distance.py # needs pandas; seconds, head-node safe
```

The per-fetus summary (`bin/fetal_summary.py`): the zero-event bound reproduces the three
embryos' cleared-cell counts, pairs count as flagged, every scorer status maps to one
`fetal_mode` value, and the sex call and XIST/Y check warn as documented (male in and out of
range, female with no check, unclear sex).

The Souporcell cross-library check (`bin/souporcell_distance.py`), on tiny Souporcell-shaped runs
written to a temporary folder: a fetus whose libraries split into two halves of one genotype is
silent; a library whose second cluster carries another genotype is marked `review` and warned,
and only that cluster; a fetus with one checked library is warned as not checkable; clusters
below the singlet minimum are shown, not judged; the heatmaps are written and the report embeds
them with the check's warnings.

### Souporcell check regression

```bash
python3 tests/regression_souporcell_distance.py <loss_heatmap/pairs.tsv> <port outdir> [--tol 1e-9]
```

Runs on the output of `bin/souporcell_distance.py` for the paper's 30 embryo Souporcell runs
(F218, F147, F158 and F137, 8 + 8 + 6 + 8). Every pair of clusters must match the paper's
`followups/loss_heatmap/pairs.tsv` (compute.py): the same sites matched and used, distance within
`--tol`; and no cluster may be marked `review`, since the paper found no mother. Under 3 GB, 5 min:
`<dev>/pipeline_port/souporcell_distance/regression.bsub`.

### Scorer regression

```bash
python3 tests/regression_fetal_score.py <reference prefix> <port prefix> [--tol 1e-6]
```

Compares `bin/fetal_score.py`'s output with the prototype's for the same inputs: cells, order,
sites, UMIs, chosen fetal share, status and Y counts must be equal; scores within the relative
tolerance (the prototype ran with scipy 1.16, the container has 1.15); summary counts equal;
removed sites identical. Run with real pileups on LSF (one batch needs up to 48 GB and a few
hours), in the container:

```bash
singularity exec -B <data dirs> <container> python3 bin/fetal_score.py \
    --library LIB1=<pileup>/LIB1/cellsnp [--library LIB2=...] \
    [--expression LIB1=<matrix dir> ...] --panel <soc_vcf> --prefix port
singularity exec -B <data dirs> <container> python3 tests/regression_fetal_score.py <reference prefix> port
```

### Checks run on real data (2026-11-02/03)

| Check | Inputs | Result |
|---|---|---|
| Scorer regression | F147 batch 5 (5′ v2), F218 batch 2 (3′ v3), D8 clean, D8 2% and D8 40% (decidua/blood) | All five identical; largest score difference 3 parts in 10¹⁵ |
| End to end through Nextflow (`-profile sanger`) | The six libraries of F147 batch 5, from the STARsolo BAMs | Six pileups byte-identical to the prototype's; scores identical |
| Summary on all embryo batches | The prototype's 24 batch outputs for F147, F158, F218 and their 144 expression matrices | Cleared cells, flags (0 / 5 / 8), bounds, the 13 flagged cells and F147's XIST/Y bound reproduced; no warnings |
| Souporcell check regression (2026-11-16, job 674883) | The paper's 30 embryo Souporcell runs | All 1,770 cluster pairs identical to `loss_heatmap/pairs.tsv` (msd difference 0); 60/60 clusters `ok`, 0 warnings |
| Fetal-mode DAG with the Souporcell check (`-preview`, 2026-11-16) | 7 F218 libraries (2 batches), with and without a reused `souporcell_dir`; `--fetal_souporcell per_batch`, `all`, `off`, and an invalid value | Selects 3 / 7 / 0 libraries (2 / 6 / 0 to run); the invalid value is refused |

**F137 was not run through the pipeline.** It was screened (2026-11-11) with the frozen
prototype only, in batches of five and six (`proto/final_batches_f137.txt`); the pipeline would
cut its 44 libraries into seven batches of six and one of two. The checks above cover F147,
F158, F218 and D8.

The scripts and logs are in `<dev>/pipeline_port/`
(`regression.bsub`, `e2e/run_e2e.sh` and `compare_e2e.sh`, `summary_check/summary_check.bsub`).
Rerun the end-to-end check after any change to the fetal-mode modules, the scorer regression
after any change to `bin/fetal_score.py`, and the Souporcell check regression after any change
to `bin/souporcell_distance.py`; a change to the scoring itself is a new
method version and needs the benchmark, not just a regression.
