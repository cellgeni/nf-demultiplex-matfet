# nf-demux-preg

Nextflow pipeline for **maternal/fetal demultiplexing** of single-cell samples from
pregnancies (placenta, or fallopian tube for ectopic pregnancies). Each sample contains
exactly two individuals — mother and fetus.

Two stages, kept deliberately separate (see [`ROADMAP.md`](ROADMAP.md)):

1. **`SOUPORCELL_DEMUX`** — genotype demultiplexing with [Souporcell](https://github.com/wheaton5/souporcell)
   at **K=2**. Wraps the proven path from
   [`nf-demultiplex`](https://github.com/cellgeni/nf-demultiplex) (pinned as a git submodule).
2. **`LABEL_MATERNAL_FETAL`** — labels the two genotype clusters as maternal / fetal /
   unresolved using a pregnancy marker panel, via an auditable **two-tier decision engine**.

Souporcell separates cells *by genotype*, not expression. The marker panel does **not** do
the separation — it only assigns a label to each of the two clusters. The pipeline
**preserves uncertainty**: weak, conflicting, or biologically implausible evidence is
reported as `unresolved` rather than forced into a call, and the original Souporcell calls
are always retained.

## Quick start

```bash
module load cellgen/nextflow/24.10.1   # or any Nextflow >= 22.10

# Full pipeline (BAM -> labeled cells) on the farm:
nextflow run . \
    --input samplesheet.csv \
    --outdir results \
    -profile sanger,singularity
```

### Samplesheet (full pipeline)

CSV with a header. K is fixed at 2 internally and is **not** a column.

```csv
sample_id,bam,barcodes,expression
P1_placenta,/path/P1/possorted_genome_bam.bam,/path/P1/barcodes.tsv.gz,/path/P1/filtered_feature_bc_matrix
```

| Column | Description |
|---|---|
| `sample_id` | Unique library identifier (no spaces). |
| `bam` | CellRanger BAM with cell barcodes (indexed, or the index is created). |
| `barcodes` | Filtered cell-barcode file (optionally `.gz`). |
| `expression` | CellRanger `filtered_feature_bc_matrix` directory. |

### Run only the labeling step (existing Souporcell output)

`LABEL_MATERNAL_FETAL` is independently runnable — its interface is the Souporcell output
files plus the expression matrix, so you can re-label without re-running Souporcell:

```bash
nextflow run . -entry LABEL_ONLY --input label_samplesheet.csv --outdir results
```

```csv
sample_id,souporcell_dir,expression
P1_placenta,/path/results/souporcell/P1_placenta,/path/P1/filtered_feature_bc_matrix
```

## The decision engine (two tiers)

| Tier | Rule | Result |
|---|---|---|
| 1 | **RPS4Y1** fraction-positive ≥ floor **and** ≥ N-fold enriched in one cluster | `fetal_high_confidence_male` (HIGH, male) |
| 2 | RPS4Y1 uninformative; ≥ `min_vote_margin` secondary markers enriched in one cluster | `fetal_marker_supported` (PROVISIONAL, likely female) |
| 3 | Insufficient / conflicting evidence, or a QC gate failed | `unresolved` |

**Key safeguard:** *absence* of RPS4Y1 is never treated as maternal evidence — it only moves
the decision to Tier 2. A cluster is labeled `maternal_likely` **only** once the other cluster
is positively identified as fetal. RPS4Y1 precedence and this safeguard are covered by unit
tests. Female-fetus (Tier 2) calls are internally consistent but **not yet independently
validated** — they are capped at PROVISIONAL confidence.

Panel: `RPS4Y1` (Y, decisive) + `HAPLN1, PITX2, GAPLINC, HPGDS, XAGE2, XAGE3, HLA-G`
(secondary). Edit [`assets/marker_panel.tsv`](assets/marker_panel.tsv); genes are resolved to
the matrix at runtime (symbol → alias → Ensembl) and the run **fails loudly** if any panel gene
is missing.

Thresholds (`--min_cluster_cells`, `--rps4y1_min_fold`, `--min_vote_margin`, …) have
conservative starting defaults in `nextflow.config` and are **calibrated, then frozen**, during
validation (ROADMAP M3–M4).

## Outputs (per sample, under `--outdir`)

| Path | Description |
|---|---|
| `souporcell/<id>/` | Unmodified Souporcell output (clusters, genotypes, ambient, QC). |
| `<id>/cluster_marker_summary.tsv` | Per cluster × marker: detection, fraction, CP10K, enrichment. |
| `<id>/cluster_labels.tsv` | Per cluster: label, confidence, rule, supporting markers, reason, warnings. |
| `<id>/cell_assignments.tsv.gz` | Per barcode: Souporcell call + maternal/fetal label + reason. |
| `<id>/sample_qc_summary.tsv` | Verdict, cluster balance, doublet/unassigned rates, ambient, warnings. |
| `reports/<id>.html` | Human-readable report (verdict, evidence, demux QC, provenance). |

## Testing

```bash
python3 tests/test_decision_engine.py                 # decision-engine unit tests
python3 tests/make_fixture.py                          # build the synthetic fixture
nextflow run . -entry LABEL_ONLY -profile test         # end-to-end labeling smoke test
```

See [`docs/testing.md`](docs/testing.md) for details.

## Status

Early development (`v0.1.0-dev`). The Souporcell wrapper, labeling engine, tests, and reporting
are in place and validated on synthetic data; threshold calibration and biological validation
against real pilot data are the next milestones. See [`ROADMAP.md`](ROADMAP.md).
