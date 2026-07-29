# Roadmap — Maternal/Fetal Demultiplexing Pipeline (`nf-demux-preg`)

> **Hybrid roadmap**, drafted 2026-07-29. Merges the strongest parts of two prior drafts:
> - **`ROADMAP-opus48-v2.md`** → two-subworkflow architecture, upstream (nf-demultiplex) reuse
>   strategy, and the calibration harness.
> - **`ROADMAP-gpt56sol.md`** → minimal input contract, fixed RPS4Y1 precedence, scientific
>   safeguards, test plan, and definition of done.
>
> Working document — expect revision as pilot data arrives.

---

## 1. Goal

Demultiplex **maternal** and **fetal** cells from single-cell samples of pregnancies
(placenta, or fallopian tube for ectopic pregnancies). Each sample contains exactly two
individuals — mother and fetus — so the task is two clearly separated steps:

1. **Genotype demultiplexing** — separate cells into two genotype clusters with Souporcell (K=2).
2. **Cluster labeling** — decide which genotype cluster is fetal and which is maternal, using
   expression evidence from a pregnancy marker panel.

The pipeline must **preserve uncertainty**: it does not force every sample into a confident
maternal/fetal call when marker evidence is weak, conflicting, or biologically incompatible with
the K=2 assumption. It must also **retain the original Souporcell calls** rather than replacing
them with biological labels.

## 2. Key conceptual framing

Souporcell separates cells **by genotype** (SNPs in the reads), *not* by expression. Mother and
fetus are two individuals: the fetus carries paternal alleles the mother lacks, so the two are
cleanly separable even though they share ~50% of alleles.

**The marker panel does not do the separation — it only assigns a label** to each of the two
genotype clusters. Consequences:

- Separation quality depends on **genetic informativeness** (SNP coverage), not on whether the
  marker genes are expressed.
- **Female fetuses separate genetically just as well as male fetuses.** The only extra difficulty
  for a female fetus is *labeling*, because RPS4Y1 (Y chromosome) will not fire.
- **Maternal identity is inferred as the non-fetal cluster**, not from a dedicated maternal panel.

**Highest-value future enhancement:** obtain **maternal blood/DNA genotype** if the cohort has it.
Matching it against Souporcell's `cluster_genotypes.vcf` gives direct maternal-cluster
identification and true validation for female fetuses.

## 3. Locked design decisions (2026-07-29)

| Decision | Choice | Implication |
|---|---|---|
| Data modality | **10x scRNA-seq (3'/5')** | Best case for marker + RPS4Y1 capture. |
| Architecture | **End-to-end pipeline, two subworkflows** | `SOUPORCELL_DEMUX` (reused from nf-demultiplex) + `LABEL_MATERNAL_FETAL` (new). |
| Processing unit | **One library/sample, independently** | No pooled Souporcell, no cohort join in v1. |
| Genotypes per sample | **Always exactly K=2 (mother + fetus)** | Fixed K=2, not a routine user setting. Keep an "is this really 2 genotypes?" QC gate. |
| Souporcell vs. labeling | **Separate steps** | Souporcell calls retained; labeling never mutates them. |
| RPS4Y1 role | **Decisive when present; absence ≠ maternal** | Explicit, unit-tested precedence. |
| Validation ground truth | **RPS4Y1 (male fetuses) only; none external yet** | Male samples calibrate female thresholds; female calls internally consistent but **not independently validated**. |
| nf-demultiplex ownership | **We own/control it** | Pin as submodule; wrap the Souporcell path. |

## 4. Architecture  *(from Opus-v2)*

**End-to-end pipeline** (`BAM + barcodes + expression → labeled cells`, one command), composed of
two DSL2 subworkflows:

1. **`SOUPORCELL_DEMUX`** — runs Souporcell at K=2, publishes unmodified outputs.
2. **`LABEL_MATERNAL_FETAL`** — the new labeling + QC + report layer.

**Key design rule:** `LABEL_MATERNAL_FETAL`'s input interface is the **Souporcell output files**
(`clusters.tsv`, `ambient_rna.txt`, `cluster_genotypes.vcf`) + the expression matrix — *not*
nf-demultiplex internals. This gives all of: end-to-end UX; the labeling step is testable
standalone against fixture `clusters.tsv` (no Souporcell run needed in dev); anyone with existing
Souporcell output can run just the labeling step; and the coupling point is the **stable Souporcell
output format**.

The decision axis is not "one pipeline vs. two" — it is **where the single source of truth for the
Souporcell invocation lives**. The thing to avoid is copy-pasting Souporcell logic into two
independently-maintained places. ("Re-label without re-running Souporcell" is provided for free by
Nextflow `-resume`, so it was never a reason to split.)

### Sourcing the Souporcell step (we own nf-demultiplex)  *(from Opus-v2)*

nf-demultiplex is DSL2 but *not* modularized for import: all processes live in `main.nf` (no
`include`s), and Souporcell is `process run_souporcell` calling `bin/souporcell.sh` +
`bin/soup_qc.sh`. Verified upstream interface:

```
process run_souporcell {                       # tuple(name, barcodes, bam, index, K)
  bin/souporcell.sh  name barcodes bam known_genotypes soc_vcf soc_fasta K skip_remap no_umi
  bin/soup_qc.sh     name
}
```
Outputs consumed downstream: `clusters.tsv`, `cluster_genotypes.vcf`, `ambient_rna.txt`,
`doublets.err`, `qc.txt`.

- **Now (MVP):** pin nf-demultiplex as a **git submodule**; define a thin `SOUPORCELL` process in
  `nf-demux-preg/modules/` that calls the submodule's `bin/souporcell.sh` + `bin/soup_qc.sh`. Heavy
  logic stays single-source; only the ~10-line wrapper is local. K is hard-set to 2.
- **Later (optional):** promote `run_souporcell` into a proper `modules/souporcell.nf` inside
  nf-demultiplex (backward-compatible), then both pipelines `include` the same module. Do not block
  MVP on this.

## 5. Upstream foundation — nf-demultiplex

Repo: <https://github.com/cellgeni/nf-demultiplex>

- Sample sheet columns: `sample_id, BAM, barcodes.gz, K` → **K=2 per line, no upstream code change.**
- Reference defaults (from upstream `nextflow.config`, reused verbatim):
  - `soc_fasta = /nfs/cellgeni/STAR/human/2020A/GRCh38_v32_modified.fa`
  - `soc_vcf   = /nfs/cellgeni/pipeline-files/demultiplex/reference_vcfs/chr_filtered_2p_1kgenomes_GRCh38.vcf`
  - `skip_remap = true`, `no_umi = false`, `known_genotypes = false`
  - container: `/nfs/cellgeni/singularity/images/souporcell_v2.5-python_3.10-vireo_v0.5.8-cellsnp-lite_1.2.3-nextflow-compatible.sif`
- Already has a `--check_sex` step using **XIST/RPS4Y1** expression — reusable prior art for the
  RPS4Y1 logic.
- Souporcell outputs consumed by `LABEL_MATERNAL_FETAL`:
  - `clusters.tsv` — `barcode, status{singlet/doublet/unassigned}, assignment, log_prob_singleton, log_prob_doublet, cluster0, cluster1`
  - `cluster_genotypes.vcf` — per-cluster genotypes (for the genotype-matching enhancement)
  - `ambient_rna.txt` — ambient RNA fraction estimate (feeds confidence)

## 6. The marker gene panel

| Gene | Note |
|---|---|
| **RPS4Y1** | **Y chromosome — highest priority.** Cell-type-independent; fires in all fetal cells of a *male* fetus. Absent (beyond ambient) in maternal cells. |
| HAPLN1 | Fetal/placental. |
| PITX2 | Fetal/placental. |
| GAPLINC | Fetal/placental (lncRNA). |
| HPGDS | Fetal/placental. |
| XAGE2 | Fetal/placental (cancer-testis / X-antigen family). |
| XAGE3 | Fetal/placental (cancer-testis / X-antigen family). |
| HLA-G | Extravillous trophoblast — fetal. |

- **RPS4Y1 is decisive for a *male* fetal cluster** when cleanly cluster-specific, but
  **absence of RPS4Y1 must NOT be read as maternal evidence** — the fetus may be female, or the
  gene simply undetected. This is enforced structurally in the decision engine (§8).
- **Caveat:** all genes except RPS4Y1 are largely **trophoblast/placental-lineage,
  cell-type-specific** markers. For a female fetus whose captured fetal cells are not
  trophoblast/stromal, several may legitimately stay silent — so aggregate scoring beats a
  single-gene rule.
- **Symbol/alias resolution required:** CellRanger matrices key on Ensembl IDs + symbols;
  `HLA-G` and `GAPLINC` in particular have alias/formatting pitfalls. The pipeline resolves each
  panel gene (symbol → Ensembl ID, with alias fallback) at preflight and **fails loudly** if a gene
  cannot be found in the matrix. The marker table is a small, version-controlled asset.

## 7. Inputs — minimal contract  *(from GPT-5.6)*

Each library/sample is processed independently. The v1 samplesheet contains **only the fields
needed to run one sample** — no pregnancy-level grouping field is required.

```text
sample_id    bam    barcodes    expression
```

Required fields:

| Field | Description |
|---|---|
| `sample_id` | Unique library/sample identifier. |
| `bam` | CellRanger (or equivalent) aligned BAM with cell barcodes (index alongside, or created). |
| `barcodes` | Filtered cell-barcode file matching the expression matrix. |
| `expression` | Gene-expression matrix — CellRanger `filtered_feature_bc_matrix` is the v1 primary format; `.h5ad` is a possible secondary format decided during implementation. |

K is **not** a samplesheet column in this pipeline — it is fixed at 2 internally.

No optional/metadata fields are required in v1. `pregnancy_id`, `tissue`, and `known_fetal_sex` are
**deferred**: they add cross-library and validation bookkeeping but are not needed to make a
per-library call, and keeping them out keeps the v1 contract minimal. They return as opt-in columns
once the core is validated (see §13).

`LABEL_MATERNAL_FETAL` subworkflow interface (its stable contract, for standalone runs):
`souporcell_dir` (`clusters.tsv`, `ambient_rna.txt`, `cluster_genotypes.vcf`) + `expression`.

Potential future fields: `pregnancy_id`, `tissue`, `known_fetal_sex`, `known_genotypes_vcf`,
`maternal_genotype_vcf`, `expected_fetal_fraction`, `reference_genome`, `barcode_prefix`.

## 8. The labeling logic — explicit two-tier decision system  *(GPT-5.6 rule, Opus-v2 ambient robustness)*

The naive rule — "expresses one or more markers → fetal" — is prone to **ambient-RNA false
positives** (placenta is high-ambient; fetal transcripts leak into the maternal cluster). The rule
is therefore a **relative comparison between the two clusters** with RPS4Y1 as an override, kept as
an **auditable scoring table, not a black-box classifier**. Only **singlets** feed the decision;
doublets and unassigned cells inherit no maternal/fetal label.

**Tier 1 — Highest-priority fetal evidence (RPS4Y1, cell-type-independent).**
If one cluster has credible, *specific* RPS4Y1 expression — fraction-positive ≥ *floor* (start
~10%) **AND** ≥ *N*-fold enriched over the other cluster — label that cluster `fetal` and the other
`maternal`. Confidence = **high**, inferred sex = **male**. Overrides all other logic.
*(The enrichment requirement stops ambient RPS4Y1 in the maternal cluster from tripping the call.)*

**Tier 2 — Other fetal evidence (female fetuses / weak-Y).**
If RPS4Y1 is not informative, each of the other 7 genes votes for whichever cluster it is
sufficiently enriched in (≥ ~2-fold, with a fraction-positive floor). Higher tally → `fetal`, other
→ `maternal`; confidence = **provisional** (higher if multiple markers agree). Require a clear vote
margin, else → unresolved.

**Tier 3 — Insufficient or conflicting evidence → `unresolved`.**
Never force a label when neither cluster passes, both pass similarly, secondary markers disagree
strongly, or a QC gate has failed.

**Structural safeguard (from §6):** absence of RPS4Y1 can only *move the decision to Tier 2* — it
can never contribute maternal evidence. This is a code invariant, covered by a unit test.

Every label carries: **confidence tier** (`high` / `provisional` / `unresolved`); **the rule that
fired**; **supporting marker(s)**; **key cluster expression statistics**; and **any warning** that
reduced confidence.

**Thresholds above are starting points — calibrated in Phase 4 (§11), frozen only after validation.**

### Evidence → interpretation (first-pass)

| Evidence | Interpretation |
|---|---|
| Strong RPS4Y1 in one cluster only | Fetal, high confidence (male). |
| One non-Y marker strongly enriched in one cluster | Fetal, provisional. |
| Multiple non-Y markers enriched in one cluster | Fetal, provisional→high. |
| Markers detected in both clusters at similar strength | Unresolved (ambient). |
| No fetal markers detected in either cluster | Unresolved / maternal-only / low-signal warning. |
| RPS4Y1 detected in **both** clusters | Warning: ambient RNA, mixed cluster, or demux issue. |

### Tunable parameters (conservative defaults; calibrated on pilots)

| Parameter | Role |
|---|---|
| `min_cluster_cells` | Minimum singlet cells in a cluster to make a confident call. |
| `min_marker_detected_fraction` | Minimum fraction of cells expressing a marker in the fetal candidate cluster. |
| `min_marker_fraction_delta` | Minimum difference in expressing-fraction between clusters. |
| `min_marker_logfc` | Minimum expression enrichment between clusters. |
| `rps4y1_fraction_floor` | RPS4Y1 fraction-positive floor for the Tier-1 override. |
| `rps4y1_min_fold` | RPS4Y1 enrichment fold required for the Tier-1 override. |
| `min_vote_margin` | Minimum Tier-2 vote margin before a female call is made. |

## 9. Scientific & technical safeguards  *(from GPT-5.6 §5)*

- **RPS4Y1 is decisive when present, but absence is not maternal evidence.** A female fetus supplies
  no Y expression; a missing RPS4Y1 moves the decision to Tier 2, never to a maternal label.
- **Expression must be cluster-specific.** Combine detection with contrast between the two clusters.
  Weak expression in both clusters likely reflects ambient RNA, not two fetal clusters.
- **Marker coverage depends on cell type and tissue.** Placental and ectopic samples contain
  different fetal cell types; a marker may not behave equally across tissues. Validation must
  include both placenta and fallopian-tube material and stratify by major cell type where reference
  annotation is available.
- **K=2 is a biological assumption, not just a parameter.** Low fetal fraction, severe genotype
  imbalance, multiple gestation, contamination, or lack of genetic variation can all make a
  two-cluster result unreliable. The report must **expose cluster sizes, assignment rate, doublet
  rate and other Souporcell QC** rather than presenting two biological labels regardless of
  clustering quality.

### Edge cases (must be handled explicitly)

- **Effectively one genotype** (near-pure sample): K=2 splits noise → flag, don't label.
- **Extremely imbalanced cluster sizes** → flag; may indicate low fetal fraction or one genotype.
- **Both clusters marker-positive** → ambient contamination; lean on RPS4Y1 / relative enrichment.
- **Neither cluster positive** → fetal cells absent or panel off → flag.
- **Doublets / unassigned** → kept as their own labels; never coerced into maternal/fetal.

## 10. Outputs

Stable per-sample result layout.

### Demultiplexing outputs
- `souporcell/<sample_id>/` — original Souporcell assignments, `cluster_genotypes.vcf`, logs,
  `ambient_rna.txt`, `qc.txt`; the exact effective parameters and provenance.

### Pregnancy-interpretation outputs

| Output | Description |
|---|---|
| `cluster_marker_summary.tsv` | One row per cluster × marker: detection count/fraction, pseudobulk & normalised expression, enrichment vs. the other cluster, cells used, pass/fail flag. |
| `cluster_labels.tsv` | Souporcell cluster → maternal/fetal/unresolved label + confidence + rule that fired + supporting markers + warnings. |
| `cell_assignments.tsv.gz` | Per barcode: sample_id, barcode, original Souporcell cluster & call, maternal/fetal/unresolved label, confidence, decision reason. |
| `sample_qc_summary.tsv` | Per-sample QC, warnings, final verdict. |
| `reports/<sample_id>.html` | Human-readable report (see below). |

### Label vocabulary (separates biological call from confidence)

| Label | Meaning |
|---|---|
| `fetal_high_confidence_male` | RPS4Y1 evidence clearly identifies the fetal cluster. |
| `fetal_marker_supported` | Non-Y fetal markers support the fetal cluster (provisional). |
| `maternal_likely` | The opposite cluster after the fetal cluster is identified. |
| `unresolved` | Evidence insufficient or conflicting. |
| `unassigned` | Souporcell unassigned barcode. |
| `doublet` | Souporcell doublet barcode. |

### Human-readable report
Leads with the verdict (fetal/maternal cluster IDs, confidence), then: cell counts per Souporcell
call and cluster; cluster balance, singlet-assignment and doublet rates; marker expression &
detection by cluster; a clearly highlighted **RPS4Y1** result; final interpretation and evidence;
warnings and unresolved cases; and software / reference / parameter provenance. **Denominators are
explicit** (total barcodes, singlets, cells per cluster, marker-positive cells per cluster) and the
report must be interpretable from tables/text alone, not only color-coded plots.

An annotated expression object (`.h5ad`) can be considered after the tabular contract is stable; the
TSV outputs remain the primary machine-readable interface.

## 11. Milestones / build phases  *(Opus-v2 milestones with GPT-5.6 exit criteria)*

### M0 — Repository skeleton
`main.nf`, `nextflow.config`, `conf/`, `bin/`, `subworkflows/`, `modules/`,
`assets/schema_input.json`, `assets/marker_panel.tsv`, `docs/`, `test_data/`, `README.md`;
nf-demultiplex pinned as a submodule. Keep the first implementation deliberately narrow.
- **Exit:** a short architecture decision, fixed input contract, and named validation data are agreed.

### M1 — Pilot data (before writing labeling code)
Run Souporcell at K=2 on a handful of pilot samples, **including some expected to be male**.
Produces real `clusters.tsv` / `cluster_genotypes.vcf` / `ambient_rna.txt` to develop and calibrate
against.
- **Exit:** one representative sample completes reproducibly; Souporcell outputs match the chosen
  upstream baseline.

### M2 — Minimal viable pipeline (`LABEL_MATERNAL_FETAL` core)
Developed/tested standalone against fixture Souporcell outputs, then wired behind `SOUPORCELL_DEMUX`
end-to-end. Samplesheet schema; preflight validation; marker-aggregation script; two-tier decision
engine; cluster- and barcode-level TSVs; basic HTML report.
- **Exit:** one pilot sample runs sample-sheet → maternal/fetal output; the final call is traceable
  from the report back to raw cluster and marker evidence; all synthetic decision tests pass.

### M3 — Confidence & QC hardening + calibration harness  *(from Opus-v2)*
Configurable thresholds; explicit confidence tiers; warnings for weak evidence / invalid K=2;
per-sample QC summary; regression tests on tiny synthetic inputs. **Calibration harness:** on male
samples (RPS4Y1 = truth), sweep thresholds and ask — *would the marker-only logic reproduce the
Y-based label if RPS4Y1 were hidden?* Freeze thresholds for female samples; document that female
calls carry this as an assumption, not a validation.
- **Exit:** ambiguous cases are flagged, not forced; reports explain why.

### M4 — Pilot validation
Run pilot placenta and ectopic/fallopian samples; compare with known fetal sex / genotype where
available; review marker behavior by tissue; tune thresholds conservatively.
- **Exit:** performance is quantified; thresholds are justified; method limits documented. Male
  samples with clear RPS4Y1 correctly identified; female samples called when non-Y support is
  strong; weak/contradictory cases remain unresolved.

### M5 — Production readiness / v1 release
Full docs + examples; CI or local smoke tests; regression fixtures; pin containers, upstream
revision, and reference resources; test restart/resume and failure modes on FARM; release checklist.
Keep `LABEL_MATERNAL_FETAL` keyed on Souporcell outputs so it stays independently runnable.
- **Exit:** a new user runs a sample from a standard samplesheet, reproduces the test result, and
  interprets every published field from documentation alone.

## 12. Validation & test plan  *(from GPT-5.6 §7)*

### Unit tests — the decision engine
Small synthetic tables covering:
- RPS4Y1 in cluster 0, and in cluster 1;
- **RPS4Y1 taking precedence over conflicting secondary markers** (the precedence invariant);
- **absent RPS4Y1 never producing a maternal label** (moves to Tier 2);
- a female-fetus pattern supported only by secondary markers;
- one, several, or no secondary markers passing;
- weak or equal marker expression in both clusters → unresolved;
- missing marker genes (fail loudly at preflight);
- an empty or very small cluster;
- non-singlet cells excluded from evidence;
- deterministic propagation of cluster labels to barcodes.

### Integration tests
- Run a minimal fixture through input validation → Souporcell output parsing → marker aggregation →
  reporting.
- Confirm barcode identity and counts at every hand-off.
- Confirm reruns produce **identical** tabular labels from identical inputs.
- Confirm a failed/ambiguous classification is published as `unresolved` and does **not** fail the
  completed Souporcell analysis.

### Biological validation cohort
Include, where available: male pregnancies (RPS4Y1 positive control); female pregnancies (secondary
panel); both placental and fallopian-tube ectopic samples; a range of fetal fractions and cluster
imbalances; samples with independent genotype/metadata truth; and known difficult/failed demux
examples. Evaluate: agreement with independent truth; rate of correct high-confidence labels; rate
of unresolved; false fetal labeling from ambient; marker sensitivity by tissue/cell type; call
stability across plausible thresholds; Souporcell assignment/doublet/unassigned rates. Thresholds
and confidence policy are frozen only after this analysis.

## 13. Definition of done for version 1  *(from GPT-5.6 §10)*

Version 1 is complete when:
- every sample is independently demultiplexed with Souporcell K=2;
- upstream-equivalent Souporcell results are retained and reproducible;
- the eight markers are summarised by cluster from correctly matched **singlets**;
- **RPS4Y1 has explicit, tested precedence**;
- female-fetus and ambiguous cases are handled **without treating absent RPS4Y1 as proof of maternal
  origin**;
- maternal/fetal/unresolved decisions are **deterministic and auditable**;
- per-cluster and per-cell machine-readable outputs are stable;
- the report exposes **both** biological evidence **and** demultiplexing quality;
- performance has been evaluated against independent truth across relevant tissue/pregnancy contexts;
- the limitations of the K=2 assumption and marker panel are documented.

## 14. Open questions

- **Fetch maternal genotypes now vs. later?** Drives whether v1 includes a genotype-matching path
  (direct maternal-cluster ID + real female-fetus validation).
- **Expression input format:** `filtered_feature_bc_matrix` only in v1, or `.h5ad` too?
- **Credible-expression definition per marker:** minimum cells, proportion detected, pseudobulk
  count, enrichment, or a combination?
- **Secondary markers equal in v1**, or prior validation supporting different weights?
- **Unresolved output granularity:** strictly `unresolved`, or also show a best-effort candidate
  label kept clearly distinct from a call?
- **Which upstream nf-demultiplex revision** is the pinned reproducible baseline?
- **When do `pregnancy_id` / `tissue` / `known_fetal_sex` re-enter** as opt-in columns?

## 15. Possible later extensions

- Known-genotype matching where maternal or fetal genotypes are available.
- A Vireo comparison / consensus demultiplexing mode.
- Evidence from additional sex-linked genes while preserving RPS4Y1's special role.
- Tissue- or cell-type-aware fetal marker panels.
- Support for multiple gestations or additional contaminating genotypes.
- Cohort-level summaries after the independent per-sample workflow is stable.
