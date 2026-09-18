# Testing

Three layers, from fastest to most integrated.

The Souporcell shell wrapper also has a lightweight argument regression test. It uses a fake
`souporcell_pipeline.py`, so it neither reads a BAM nor runs Souporcell:

```bash
bash tests/test_souporcell_wrapper.sh
```

## 1. Decision-engine unit tests (no dependencies beyond Python 3)

The two-tier decision logic lives in pure, importable functions in
`bin/label_clusters.py` (`build_cluster_stats`, `decide`). Tests cover the cases in
ROADMAP §12:

```bash
python3 tests/test_decision_engine.py     # built-in runner, prints PASS/FAIL per test
# or, if pytest is available:
pytest tests/test_decision_engine.py
```

Key invariants under test:
- RPS4Y1 identifies the fetal cluster (either cluster id) at HIGH confidence.
- **RPS4Y1 takes precedence** over conflicting secondary markers.
- **Absent RPS4Y1 never yields a maternal label** — it only moves to Tier 2; with no other
  evidence the result is `unresolved`, not maternal.
- Female pattern (secondary markers only) → PROVISIONAL, likely female.
- Vote margins, ambiguous/equal expression, single-cluster and small-cluster QC gates.
- Determinism.

## 2. Synthetic fixture + full CLI path

`tests/make_fixture.py` writes a tiny, deterministic 10x matrix and Souporcell outputs for a
male-fetus sample under `test_data/sampleM/` (cluster 0 = maternal, cluster 1 = fetal with
RPS4Y1). Regenerate with:

```bash
python3 tests/make_fixture.py
```

You can drive the three `bin/` scripts directly (`aggregate_markers.py` →
`label_clusters.py` → `make_report.py`) against this fixture; see the commands baked into the
`LABEL_ONLY` Nextflow path below.

## 3. End-to-end labeling smoke test (Nextflow, no BAM/samtools needed)

The `LABEL_ONLY` entry runs `MARKER_SUMMARY → LABEL_CLUSTERS → BUILD_REPORT` against the
fixture Souporcell output — exercising the whole labeling subworkflow without a Souporcell run:

```bash
module load cellgen/nextflow/24.10.1
nextflow run . -entry LABEL_ONLY -profile test
```

Expected verdict for `sampleM`: `rps4y1` / HIGH / male, fetal cluster `1`, and doublets /
unassigned barcodes preserved (never coerced into a maternal/fetal label). Outputs land in
`results_test/`.

## Full Souporcell run (farm)

Requires the Souporcell container, reference FASTA/VCF, and a real BAM. Use
`test_data/samplesheet_test.csv` as a template and run with `-profile sanger,singularity`.
This is validated as part of milestone M1 against pilot data, not in CI.
