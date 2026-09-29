# Handoff — implement hybrid group-aware demultiplexing in Nextflow

> **You're picking this up cold.** Items marked 🟡 were done but never checked.
> Dead ends is not optional.

**Anchor** · as of 2026-09-29 · validation finished and the design decided;
no production code written; repository checked out on
`validation/group-matrix-core` at `0833dd0` plus this brief.
**Confidence** · written from one long conversation, not from a transcript.
Numbers below were re-read from result files on 2026-09-29 where marked; the
early part of the conversation (2026-09-18/19) is recalled, and the detailed
record for it is `docs/group_aware_redesign.md`. Unclear whether the user has
decided what, if anything, from the validation branches should reach `main`.

## Orientation

`nf-demux-preg` demultiplexes maternal and fetal cells in pregnancy scRNA-seq.
Today it runs Souporcell at K=2 **per library**, which on REQ-68805 gave one
usable library in nine; the others were fake 50:50 splits. A validation
programme (all on branches, none in the pipeline) showed that inferring
genotypes **once per pregnancy group** fixes this, and chose a **hybrid**:
Souporcell discovers the two genotypes on the group; Vireo assigns every
library's cells against those fixed genotypes; a composition gate decides
whether a second genotype is real. Validation is closed by agreement.

**Next action:** implement the hybrid in the Nextflow workflow on a new
integration branch — *based on: "make sure all of these results are noted down
so a new Claude Code session can start implementation in the Nextflow"*. The
first step (inferred, not requested) should be to read the current workflow,
propose the new process DAG and samplesheet contract to the user, and get it
agreed before writing code: this replaces the core of the pipeline.

Branching for that work: create an integration worktree under
`/nfs/cellgeni/batu/_repos/worktrees` from `validation/group-matrix-core`. Note
that the two adapters are **not** on that base: the Souporcell adapter is only
on `validation/group-souporcell` and the hybrid-capable Vireo adapter only on
`validation/group-vireo`. Bring in what you need from both (the original plan
was to build the integration branch from the selected backend rather than
merging both experimental arms wholesale).

## Where each thread got to

### The design — settled

1. **Per-library allele counting**, unchanged in spirit: remap + Vartrix per
   library as the current `SOUPORCELL` process already does (it keeps
   `alt.mtx`, `ref.mtx`, `barcodes.tsv`, `common_variants_covered.vcf`).
2. **Group matrix**: merge a group's per-library matrices on the union of
   variants, barcodes prefixed `<sample_id>::<raw_barcode>`. Code exists:
   `validation/group_matrix/build_group_matrix.py` (+ `check_group_matrix.py`).
3. **Discovery with Souporcell** on the group matrix: clustering binary →
   `troublet` → `consensus.py`, K=2. Produces `cluster_genotypes.vcf` (the
   profiles). Adapter: `backends/souporcell_adapter.py` on branch
   `validation/group-souporcell`.
4. **Projection with Vireo**: `vireo --vartrixData=alt,ref,barcodes,variants.vcf.gz
   --donorFile cluster_genotypes.vcf --genoTag GT --nDonor 2`, **never**
   `--forceLearnGT`. Adapter: `backends/vireo_adapter.py` on branch
   `validation/group-vireo` (commit `e6a9c34` made it accept Souporcell profiles).
5. **Composition gate** on the projection output (see below), then per-library
   results joined back through `barcode_provenance.tsv`.
6. **Maternal/fetal orientation** afterwards, from marker evidence — the
   existing `LABEL_MATERNAL_FETAL` subworkflow, which will need adapting (it
   currently reads a per-library Souporcell `clusters.tsv`).

### Why this design — settled, with evidence

All from `<ticket>/group_matrix/` result files; ticket =
`/lustre/scratch127/cellgen/cellgeni/tickets/REQ-68805`.

| Finding | Evidence |
|---|---|
| Group inference fixes fake splits | CW5 joint: 363 → 7,383 / 9 cells, 364 → 3,754 / 3 (per-library baseline 3,277 / 3,124 and 1,583 / 1,683); fitted ambient 1.8% (was 84–96%); 362 matches the trusted baseline at 100% |
| Backends equally accurate where signal is good | CW5 Souporcell vs Vireo: ARI 1.0000, genotype concordance 0.985 (`cross_arm` result, recorded in doc) |
| **Souporcell cannot project** | `--known_genotypes` relearns: returned profiles matched supplied at 0.726/0.747 (CW5) and 0.0 (CW4); assignments unrelated to joint (ARI ≈ 0) |
| Vireo projects well **only with a well-trained minority profile** | 362 retention 99.6% with profile learned on all of 362 vs 31.8% on half, lost cells → doublets |
| **Souporcell discovers far better** | spike-in: Souporcell detects fetal at 0.5% (35 in 7,000); Vireo only ≥5% |
| **Neither backend declines** when there is no second genotype | 0% spike-in: both split ~50:50, evenly across all three libraries |
| **Hybrid keeps sensitivity and does not invent a fetus** ✅ | re-read 2026-09-29 from `spikein/report/detection_curve.tsv`: recall/precision 1% 0.94–0.98 / 0.93–1.00, 2% 0.98–0.99 / 0.98–1.00, 10% 0.97 / 1.00; 0% → only 0.7% of cells confident (no false pair); 0.5% → one profile (misses, safely) |

### The composition gate — leaning, thresholds provisional

For the hybrid, gate on the **projection's own signals**:

- **confident-cell fraction** (singlet, assignment probability ≥ 0.90): 0.007
  at 0% versus 0.94–0.99 when a real fetal genotype exists;
- **minority-profile size**: at 0.5% the minority collapses to ≤1 cell;
- **doublet inflation** relative to discovery: flags an under-trained fetal
  profile (the 5% dip, below).

Below the ~1% floor, report one profile or unresolved, never a guessed second
genotype. Thresholds come from 12 synthetic 7,000-cell groups of **one**
pregnancy (CW5) — confirm on real production groups before treating as fixed.
**Do not** use split-half locus replication on Souporcell profiles (see Dead
ends); it works only on Vireo's likelihood-based profiles.

### Acceptance tests for the implementation — agreed in outline

- **REQ-68805 CW5** must reproduce: 362 ≈ 915 fetal / 5,427 maternal;
  363 and 364 single-profile.
- **CW4** (5 libraries) must come out unresolved or one profile. Neither
  backend found a reproducible second genotype there in 5 discovery subsets.
- **A larger group (10–15 libraries)** is the scale test, folded into
  acceptance rather than a separate experiment.
- If an acceptance result contradicts validation, treat it as a bug to fix,
  not a reason to reopen experiments (stopping rule agreed 2026-09-28).

## State of play

| Item | Status | Evidence / location |
|---|---|---|
| Group-matrix build + 20 integrity checks | ✅ verified | CW4/CW5/D6 20/20, byte-identical rebuilds; CW5→362 subset reproduced 12,617,422 entries exactly |
| Common assignment table, comparison metrics, retention, locus replication | ✅ verified | `validation/group_matrix/common/`; `tests/test_common.py`, `tests/test_spikein.py` pass |
| Souporcell adapter (discovery) | ✅ verified | parser reproduces 362 baseline exactly (911/5,465/81/790) |
| Vireo adapter (discovery + projection) | ✅ verified | used in every gate run |
| Hybrid path (Souporcell profiles → Vireo `GT`) | ✅ verified on synthetic groups | 12 spike-in groups; Vireo matched 164,587 of 597,962 variants to Souporcell's VCF |
| Hybrid on a real, un-synthesised group (CW5, CW4) | 🟡 never run | only spike-in groups built from CW5 cells |
| Hybrid on > 7,800 cells or > 3 libraries | 🟡 never run | scale unknown for the combined path |
| Gate thresholds | 🟡 provisional | one pregnancy, synthetic groups |
| Production Nextflow implementation | ⬜ planned | nothing written |
| Samplesheet with `group_id`, `sample_role`, `expected_genotypes` | ⬜ planned | contract in doc *Maternal reference libraries* and `validation/group_matrix/groups_validation.csv` |

Measured resources (use for process labels): Souporcell group discovery CW5
(19,336 cells) 7.6 GB / 16 min, CW4 (36,800) 18.4 GB / 68 min; half-group
spike-ins 2–4 GB / 4–20 min; **hybrid projection 0.9–1.3 GB / 1–2 min** per
7,000-cell group; Vireo discovery CW4 trio 25 GB / 86 min (not used in the
hybrid). Existing per-library Souporcell with remapping: 3.5–8 h, 15–16 GB.

## Decisions and positions

- **Hybrid over either backend alone** — firm. Neither alone meets both needs:
  Souporcell discovers but relearns fixed profiles; Vireo holds profiles but
  discovers only at ≥5%.
- **Discovery must include the fetal-rich library in full** — firm. Thinning
  it (e.g. a balanced per-library cell cap) collapses retention. This argues
  against the "balanced cells from every library" option in the doc's scaling
  section.
- **`--genoTag GT`** for Souporcell-derived profiles — firm. Neither backend
  writes `GP`; Souporcell's VCF has `GT:AO:RO:T:E:GO:GN`, samples `0`/`1`,
  `./.` for missing, no `##FORMAT` headers; Vireo parses it anyway.
- **Keep remapping** — default, not re-tested. Gate 6 (remap vs cellsnp-lite
  on original BAMs) deferred as a cost optimisation only.
- **Acceptance criteria amendments** — agreed: add minority retention; do not
  accept genotype concordance alone (it passed spuriously, 0.99, on CW4's unreal
  splits).

## ⚠️ Dead ends

- **Souporcell projection via `--known_genotypes`** · relearns the supplied
  profiles, re-creating fake splits · revisit only with a different Souporcell
  release that demonstrably holds profiles fixed.
- **Vireo as the discovery engine** · floor ≥5% fetal cells; below that it
  splits along an arbitrary line (precision 0.01–0.06) · revisit if typical
  fetal fractions prove to be ≥5%.
- **Split-half locus replication as the hybrid's gate** · on Souporcell profiles
  it reads 0.75–0.93 at every fraction including 0% · works on Vireo profiles
  (0.08 arbitrary vs 0.86 real); revisit only if discovery moves to Vireo.
- **Genotype concordance across runs as a stability test** · two near-copies of
  one genotype always "agree" (D6 0.947, CW4 0.99) · keep only alongside locus
  replication.
- **Cell agreement alone for projection quality** · read 1.0000 while 72% of
  fetal cells had become doublets (demoted cells leave the comparison) · use
  minority retention.
- **K sweep / K=1 vs K=2 model selection** · the user ruled it out: *"Why would
  we need to run K=1 at all, there is no such case for maternal/fetal??"*; Vireo
  0.5.8 also crashes on `--nDonor 1` · revisit for K=4 (two pairs) or K=3
  (contamination) checks only.
- **Passing a shorter barcode list as a "subset"** · both tools need one barcode
  per matrix column; misaligns cells · materialise subsets with
  `common/subset_group_matrix.py`.
- **`cellgen/apptainer` as container runtime** · non-suid, fails "Failed to
  create user namespace" on exec hosts · use `cellgen/singularity`.
- **LSF dependencies on job IDs for resumption** · finished jobs leave the queue
  and dependencies become unsatisfiable · resume on output files.

## Preferences and constraints

- D6 premise: *"you did many things while considering D6 as a single donor,
  that's not true at all!!"* — do not describe FCA7167231 as a single-donor
  control; its composition is unresolved.
- Branching: *"commit until what we have done in a branch, then this would be
  the shared base for the new branches. We can later decide what to push to
  main."* — base new work on `validation/group-matrix-core`; don't push or merge
  to `main` unasked.
- Stopping: *"we will need to stop at some point"* → agreed that validation is
  closed; further runs are acceptance tests.
- Submission: *"Consider that I'll submit them as cellgeni-su"* — the user
  submits LSF jobs as the shared `cellgeni-su` account; give them the command.
  Check jobs with `bjobs -u cellgeni-su -w`.
- Thresholds: never tune a threshold to make a run pass.
- Summaries: the user asked for shorter write-ups (*"This is too detailed, give
  me a single paragraph"*) when reporting.
- The user's global rules in `~/.claude/CLAUDE.md` apply (ask before `bsub`,
  bootstrap LSF commands, heavy work on LSF not the head node, worktrees under
  `/nfs/cellgeni/batu/_repos/worktrees`, single-line commit messages).

## Open questions

- **What reaches `main`** — undecided; ask. The pre-existing pipeline fixes
  (container enforcement, `skip_remap=false`) are commits `608563a`, `36ff355`
  on the base branch.
- **Where group-matrix building lives in production** — the harness is stdlib
  Python run on LSF; whether it becomes a Nextflow process as-is is a design
  choice for the implementation plan.
- **Marker labelling on group output** — `LABEL_MATERNAL_FETAL` expects
  per-library Souporcell `clusters.tsv`; decide whether to adapt it to the
  common assignment table or emit a compatible file.
- **D6's true composition** — needs the user; affects only the write-up.
- **Souporcell 5% recall dip** — seed-independent and inherited by the hybrid;
  fetal cells called doublets, precision still 1.00. Cause unknown; deferred.
- **Maternal reference libraries and K=4 pools** — in the original contract
  (doc *Maternal reference libraries*), never tested.
- **nf-demultiplex submodule** is dirty (`bin/souporcell.sh`, corrected boolean
  handling) and uncommitted in its own repo; preserve it.

## Pointers

| Where | What for |
|---|---|
| `docs/group_aware_redesign.md` | full evidence; read *Follow-up result: the hybrid works*, *Detection-limit result*, *Fetal retention…*, *Gate 4 result*, *Operational notes*. Sections marked *Superseded* are pre-validation |
| `validation/group_matrix/README.md` | harness usage |
| `validation/group_matrix/common/` | build/subset/select, assignment table, genotype comparison, metrics, spike-in design and scoring |
| branch `validation/group-souporcell` (worktree `/nfs/cellgeni/batu/_repos/worktrees/group-souporcell`) | Souporcell adapter; real binary paths `/opt/souporcell/souporcell/target/release/souporcell`, `/opt/souporcell/troublet/target/release/troublet`, `/opt/souporcell/consensus.py` |
| branch `validation/group-vireo` (worktree `…/group-vireo`) | Vireo adapter incl. hybrid support (`e6a9c34`) |
| `main.nf`, `subworkflows/local/souporcell_demux.nf`, `modules/local/souporcell.nf`, `subworkflows/local/label_maternal_fetal.nf` | the current per-library pipeline to be replaced |
| `<ticket>/group_matrix/` | all validation outputs; `spikein/report/detection_curve.tsv` is the decisive table |
| image `/nfs/cellgeni/singularity/images/souporcell_v2.5-python_3.10-vireo_v0.5.8-cellsnp-lite_1.2.3-nextflow-compatible.sif` | pinned tools; clustering binary reports **2.4** despite the v2.5 name |
| doc *Operational notes for running this harness on the farm* | container module and binds, LSF pitfalls, runner rules — read before writing any job |
