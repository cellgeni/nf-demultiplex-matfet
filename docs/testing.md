# Testing

Four layers, from fastest to most integrated. The acceptance tests on real data are
listed in [`pipeline_design.md`](pipeline_design.md#acceptance-tests).

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
nextflow run . -profile test
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
