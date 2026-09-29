# Group-aware maternal/fetal demultiplexing redesign

Status: Gate 1 complete for CW4, CW5 and D6; pinned-image CLI captured;
Gate 2 adapters not yet written  
Last updated: 2026-09-19

This is the living handoff document for the redesign prompted by the first
REQ-68805 runs. It records decisions, hypotheses, validation plans, and open
questions. A proposal in this document is not a locked implementation decision
unless it is explicitly labelled as such.

## Authority and scope

For the group-aware redesign, this document supersedes the per-library-only,
fixed-K assumptions in `ROADMAP.md`. The existing roadmap remains useful for the
marker panel, uncertainty-preserving label logic, and original pipeline
architecture, but its statements that every library is processed independently
and K=2 must always be recovered per library are no longer current.

No production redesign has been implemented: the pipeline is still the
per-library Souporcell K=2 baseline with downstream marker labelling. The
group-aware work so far is validation-only code on branches, outside the public
workflow.

## Why the current design needs to change

The current pipeline runs Souporcell independently at K=2 for every library and
then uses expression markers to orient the two clusters as maternal and fetal.
The REQ-68805 pilot showed that successful task completion is not sufficient:
libraries with an under-represented fetal population can produce unstable or
artificial K=2 clusters. Enabling Souporcell remapping improved the genotype
results, but independent per-library inference and the downstream biological
labelling remain unreliable.

Libraries from the same pregnancy share the same maternal and fetal genotypes.
The redesign should therefore infer those profiles at group level, then assign
cells in each library against the shared profiles. It must permit a library to
contain only one detectable profile rather than forcing two clusters.

## Current working decisions

> **Superseded 2026-09-29.** Written before validation. The Souporcell-versus-Vireo question is settled (hybrid), and the five-library cap no longer applies to matrix construction. The decision and its evidence are in *Follow-up result: the hybrid works* below; the implementation brief is `docs/HANDOFF_hybrid_implementation.md`.

- Add a biological `group_id` (pregnancy/donor group) to the input contract.
- Keep `expected_genotypes` group-specific: normally 2, with later K=4 support
  for pools containing two maternal/fetal pairs.
- Establish shared genotype profiles once per group and retain per-library cell
  assignments and QC.
- Keep maternal/fetal orientation separate from genotype demultiplexing.
- Use known maternal reference libraries, when supplied, as genotype anchors.
- Preserve unresolved results when evidence is inadequate; do not force a
  fetal cluster or a maternal/fetal label.
- Treat five libraries as a provisional maximum for a single pooled
  **Souporcell discovery** run, not as a limit on group size.
- Do not choose Souporcell versus Vireo until the validation data have been
  compared directly.

## Handoff: resume here

Status on 2026-09-29. **Validation is complete and the design is decided**: the
hybrid -- Souporcell group discovery, Vireo fixed-profile projection
(`--genoTag GT`), and a composition gate on the projection's confident fraction,
minority-profile size and doublet inflation. It detects fetal genotypes down to
~1% of cells and, uniquely among the configurations tested, does not invent a
second genotype when none exists. See the follow-up result section.

Next, in order:

1. Implement the hybrid in the production Nextflow workflow on an integration
   branch from `validation/group-matrix-core`.
2. Acceptance-test it on REQ-68805 (CW5 must reproduce; CW4 must report
   unresolved or one profile), then on a larger group, which doubles as the scale
   test.
3. Confirm the gate thresholds on those real groups.
4. Resolve D6's composition for the write-up.
5. Deferred: Gate 6 (remapping A/B) as a cost optimisation; the cause of the
   Souporcell 5% dip.

Runners, all resumable -- re-running submits only work whose output is missing:

| Experiment | Command |
|---|---|
| Gates 2, 3, 5 | `bash <worktree>/validation/group_matrix/lsf/run_gates_2_3_5.sh` |
| Gate 4 | `bash <worktree>/validation/group_matrix/lsf/run_gate4_cw4.sh` |
| Detection limit | `bash /nfs/cellgeni/batu/_repos/nf-demux-preg/validation/group_matrix/lsf/run_spikein.sh` |

`<worktree>` is `/nfs/cellgeni/batu/_repos/worktrees/group-vireo` or
`.../group-souporcell`.

Current locations:

- repository: `/nfs/cellgeni/batu/_repos/nf-demux-preg`;
- ticket: `/lustre/scratch127/cellgen/cellgeni/tickets/REQ-68805`;
- validation samplesheet: `actions/samplesheet_validation.csv` under the ticket;
- current outputs: `results/` under the ticket;
- previous outputs: `results_old/` under the ticket;
- Nextflow log: `.nextflow.log` under the ticket;
- per-library allele inputs: `results/souporcell/<sample_id>/{alt.mtx,ref.mtx,barcodes.tsv,common_variants_covered.vcf}`.

The completed baseline was launched from the ticket directory with:

```bash
nextflow run /nfs/cellgeni/batu/_repos/nf-demux-preg/main.nf \
  --input actions/samplesheet_validation.csv \
  -profile sanger \
  -resume
```

Run provenance: Nextflow 25.04.4, session
`a32f6868-c359-4760-8d6e-e03aa0e4eb72`, Nextflow script revision
`c9b188d7f1` (a Nextflow content revision, not a Git commit). The run used
Souporcell K=2, a common 1000 Genomes variant panel, `skip_remap=false`, UMIs
enabled, and the pinned container listed below.

Do not:

- treat task success as biological validation;
- lower marker thresholds to force labels;
- interpret fitted ambient fractions literally when genotype profiles are
  unstable;
- use D6/FCA7167231 as a maternal anchor for CW4 or CW5;
- repeat BAM remapping before testing the already-generated allele matrices;
- modify or discard the current dirty worktree without first reviewing its
  diffs.

## Scaling to groups of 10-15 libraries

> **Superseded 2026-09-29.** Its balanced-cell-cap option is risky: thinning the fetal-rich library's cells collapses fixed-profile retention (99.6% to 31.8%). Discovery must include that library in full. The decision and its evidence are in *Follow-up result: the hybrid works* below; the implementation brief is `docs/HANDOFF_hybrid_implementation.md`.

Do not merge every BAM into one large Souporcell job. Remapping, variant
calling, allele counting, temporary storage, and failure cost all grow with the
combined reads and cells.

Proposed architecture:

1. Generate allele-count matrices independently for each library, in parallel,
   using one common SNP panel and compatible reference coordinates.
2. Select discovery evidence from the group. This can be either:
   - all cells from at most five suitably sized libraries; or
   - preferably for a Vireo design, a balanced cap of informative cells from
     every library, so rare fetal evidence is not lost just because its library
     was outside an arbitrary five-library subset.
3. Infer `expected_genotypes` shared profiles once from the discovery evidence.
4. Assign every cell from every library against the fixed profiles in parallel.
5. Orient the profiles as maternal/fetal using known maternal controls first,
   then group-aggregated expression evidence.

The five-library rule should therefore become a configurable resource safeguard
for the Souporcell discovery implementation. A total discovery-cell/read cap is
also needed because five very large libraries may be more expensive than ten
small libraries.

### Groups with five or fewer libraries

If the group is within the cell/read resource cap, use all its libraries in one
shared discovery analysis. This is one genotype-inference run, followed by
splitting assignments and QC back into libraries. A single low-information
library can still fail to support the second genotype; the shared design does
not manufacture evidence that is absent.

### Selecting discovery evidence without enrichment metadata

Run the cheap, reusable per-library SNP pileup first and rank libraries using
objective QC rather than filenames or acquisition order:

- number and fraction of cells with usable SNP observations;
- median covered SNPs per cell and total allele observations;
- number of informative polymorphic loci and balance of allele evidence;
- cell count, while preventing the largest library from dominating;
- optional fetal-marker expression as a weak tie-breaker, never as genotype
  truth, because these markers are tissue/cell-type dependent.

If quality is similar, select a diverse set across any available tissue, sort,
batch, or collection metadata. A preliminary independent K=2 result may be
reported diagnostically but must not be trusted as the selector when the fetal
fraction is low.

A Vireo-based design may avoid selecting five whole libraries altogether: use a
balanced informative-cell sample from all 10-15 per-library pileups for shared
profile discovery, then project all cells. This must be benchmarked for runtime,
memory, and genotype stability.

## Souporcell versus Vireo

> **Superseded 2026-09-29.** This comparison predates any measurement. Measured: equal accuracy on CW5; Souporcell discovers to 0.5% but cannot project; Vireo projects but discovers only at >=5%. The decision and its evidence are in *Follow-up result: the hybrid works* below; the implementation brief is `docs/HANDOFF_hybrid_implementation.md`.

The existing nf-demultiplex submodule and pinned container already contain both
Vireo 0.5.8 and cellsnp-lite 1.2.3, so a Vireo experiment does not require
introducing an entirely new toolchain. The current upstream Vireo path is still
per-library and will need redesign for group-aware inference.

| Requirement | Souporcell | cellsnp-lite + Vireo |
|---|---|---|
| Genotype-free inference | Supported | Supported |
| Fixed known profiles | Pipeline exposes it, but the pinned implementation must be tested | Explicitly supported through a donor VCF |
| Only maternal genotype known | Poor fit: current interface expects K supplied profiles | Strong fit: Vireo supports genotype data for only a subset while learning the rest |
| Large group execution | Full pooled workflow is expensive | Per-library pileups can run in parallel; sparse inference can be performed on balanced cells |
| Remapping/reference-bias mitigation | Built-in minimap2 path; this improved the pilot | Not built into cellSNP/Vireo; original versus remapped BAM must be tested |
| Assignment uncertainty | Singlet/doublet likelihood outputs | Posterior singlet/doublet probabilities and unassigned calls |
| Ambient RNA | Souporcell estimates a group ambient fraction during consensus | Vireo has an experimental ambient option; do not rely on it without validation |

Vireo is currently the stronger architectural candidate because partial known
genotypes match the maternal-control use case and its pileup/inference stages
can be separated cleanly. It is not yet the selected replacement because:

- remapping materially improved the current Souporcell pilot;
- mother and fetus are related and may be highly imbalanced, so performance on
  conventional unrelated-donor pools cannot be assumed;
- high ambient RNA may blur allele evidence in either method;
- matrix harmonisation and group-level discovery/projection need implementation
  and testing even though upstream nf-demultiplex already has basic Vireo steps.

A plausible hybrid is to retain a reusable remapping stage, run cellsnp-lite on
the remapped BAMs, and use Vireo for shared inference and fixed-profile
assignment. This avoids treating the choice as an all-or-nothing replacement.

### What Vireo adds that Souporcell does not

The main advantage is not that Vireo is universally more accurate. It is that
its interfaces fit this project's desired information flow:

- Vireo can take a genotype for only the known mother while learning the second
  profile; Souporcell's pipeline-level known-genotype mode expects profiles for
  all K clusters.
- Vireo treats allele counting and genotype inference as separate stages, so
  per-library pileups can be reused for group discovery, held-out projection,
  alternative discovery subsets, and reruns of the statistical model.
- Vireo can write an inferred donor-genotype VCF and then consume donor VCFs for
  fixed-profile assignment.
- Balanced subsets of cells can be used for profile learning without discarding
  the full per-library allele matrices needed for final assignment.

Souporcell has different advantages: integrated remapping and variant calling,
an explicit mature ambient/genotype consensus model, and good robustness in
genotype-free discovery. The redesign should choose on validation evidence,
not on architecture alone.

### Why remapping is being considered

The Cell Ranger BAM is well aligned for expression quantification. The question
is narrower: whether those STAR alignments are sufficiently unbiased for SNP
allele counting. Souporcell's authors report splice-alignment, soft-clipping,
and indel-versus-mismatch behaviours that can create false variants or favour
the reference allele. Their remapping step uses parameters intended for more
reliable variant evidence.

Remapping is therefore not an unconditional requirement. Its importance is
lower when allele counting is restricted to a reliable common-variant panel or
known genotypes. Because it improved the REQ-68805 pilot, it remains an
experimental factor. The redesign should test it once on a representative
group and then choose one default, rather than maintaining permanent remapped
and non-remapped production paths.

## Maternal reference libraries

A known maternal-only library such as peripheral blood should be represented
explicitly, not mixed into the input as an ordinary tissue library.

Proposed samplesheet roles:

| Field | Example | Meaning |
|---|---|---|
| `group_id` | `CW4` | Libraries and references sharing one pregnancy |
| `sample_id` | `CW4_maternal_blood` | Unique identifier; hypothetical matched reference |
| `sample_role` | `mixture` or `maternal_reference` | Whether cells require demultiplexing or anchor a known profile |
| `expected_genotypes` | `2` | Number of shared profiles for the group |

For each `maternal_reference` library:

1. Perform SNP pileup in pseudo-bulk or aggregate high-quality cells.
2. Create a maternal genotype-likelihood VCF with depth, quality, missingness,
   and contamination QC.
3. Combine multiple compatible maternal references for the group if available.
4. With Vireo, supply the maternal profile as the one known donor and request
   K=2 so the fetal profile is learned from the mixture libraries.
5. Use the reference match to orient genotype labels directly; use expression
   markers as confirmation and to investigate conflicts.

Reference libraries should normally contribute genotype evidence but should not
be included in tissue-level fetal-fraction reporting. A reference/profile
mismatch must fail group labelling rather than silently choosing the nearest
cluster. For K=4, controls need an explicit relationship to the relevant mother
or pair; a generic `maternal_reference` label is not sufficient.

## Validation plan

> **Superseded 2026-09-29.** The plan was run (Gates 1-5 plus the detection-limit experiments); Gates 6 and 7 are deferred. The description of D6 as a single-donor negative control was questioned by the user and is unresolved. The decision and its evidence are in *Follow-up result: the hybrid works* below; the implementation brief is `docs/HANDOFF_hybrid_implementation.md`.

Validation set available under REQ-68805:

- CW4: Embryo-Imp15289123 through Embryo-Imp15289127 (five libraries, including
  different enrichments).
- CW5: Embryo-Imp15350362 through Embryo-Imp15350364 (independent pregnancy).
- D6/FCA7167231: unrelated maternal peripheral-blood, used only as a
  single-donor negative control. There is no matched maternal reference for CW4
  or CW5 in the current validation set.

Current staged comparisons, deliberately avoiding four separate production
implementations:

1. **Gate 1 -- matrix integrity:** harmonise the existing per-library Vartrix
   matrices for CW5 and prove that row/column identities and allele counts are
   preserved.
2. **Gate 2 -- inference comparison with identical evidence:** run Souporcell
   clustering and Vireo on the same combined CW5 matrices. Stop a backend if it
   is unstable, cannot consume the contract safely, or fails the criteria below.
3. **Gate 3 -- shared-profile behaviour:** learn profiles using the informative
   CW5/Embryo-Imp15350362 library, project 363 and 364, and compare this with
   joint three-library inference. No library is required to contain both
   profiles.
4. **Gate 4 -- difficult-group replication:** apply the surviving approach to
   CW4 and test profile stability across alternative discovery subsets. CW4 has
   no trustworthy independent anchor.
5. **Gate 5 -- one-donor guard:** run D6 as a separate K=2 stress test and verify
   that QC rejects a credible second-profile interpretation. D6 is not projected
   against the unrelated CW4/CW5 profiles.
6. **Gate 6 -- remapping decision:** only if Vireo remains promising, run one
   representative original-BAM cellsnp-lite/Vireo comparison against Vireo on
   the remapped Vartrix evidence. Select one preprocessing default.
7. **Gate 7 -- scale:** simulate 10-15-library discovery using balanced cell
   caps and measure wall time, peak memory, temporary disk, and profile
   stability.
8. **Implementation:** build only the winning path into the production Nextflow
   workflow. Keep validation adapters outside the public pipeline.

Primary evaluation criteria:

- shared-genotype concordance across discovery subsets;
- held-out cell-assignment concordance and posterior/likelihood separation;
- informative SNP coverage and unassigned/doublet rates;
- ability to leave a genuinely single-profile library as single-profile;
- maternal-reference match strength and absence of false fetal calls in the
  reference;
- consistency with group-level fetal marker evidence, without treating those
  markers as genotype ground truth;
- runtime, memory, disk, and retry cost.

## Parameters proposed for the redesign

> **Superseded 2026-09-29.** `demux_backend` as a switch is no longer needed: the design fixes Souporcell discovery plus Vireo projection. The decision and its evidence are in *Follow-up result: the hybrid works* below; the implementation brief is `docs/HANDOFF_hybrid_implementation.md`.

Names are provisional:

- `group_id` and `sample_role` in the samplesheet.
- `expected_genotypes` per group.
- `max_discovery_libraries = 5` for a Souporcell discovery backend.
- `max_discovery_cells` and/or `max_discovery_cells_per_library`.
- `discovery_samples`: optional manual override; otherwise QC-based selection.
- `demux_backend`: initially an experimental Souporcell/Vireo comparison, not a
  permanent user-facing switch until validation identifies a default.
- maternal-reference genotype quality and profile-match thresholds.
- fixed-profile assignment confidence and minimum informative-SNP thresholds.

## Open decisions

> **Superseded 2026-09-29.** Most of these are now decided; the remaining open items are listed in the handoff brief. The decision and its evidence are in *Follow-up result: the hybrid works* below; the implementation brief is `docs/HANDOFF_hybrid_implementation.md`.

- Replace Souporcell with Vireo, retain Souporcell, or use a remap + Vireo
  hybrid.
- Whether discovery should select up to five complete libraries or balanced
  cells across every library.
- Safe discovery cell caps and genotype/profile QC thresholds.
- Whether production matrix generation should use Vartrix after remapping or
  cellsnp-lite on the original Cell Ranger BAM.
- How a future matched maternal reference will be supplied and validated; the
  current validation set contains no matched reference for CW4 or CW5.
- Ambient-RNA handling after observing the comparative validation results.
- K=4 grouping, maternal-reference identities, and mother/fetus pairing logic.

## Repository and branch state

As recorded on 2026-09-21. `main` is still at `2fc50e4` and nothing has been
pushed; every decision about what reaches `main` is still open.

Branches, all in this repository, with worktrees below
`/nfs/cellgeni/batu/_repos/worktrees`:

- `validation/group-matrix-core` at `9f004c2` -- the reviewed shared base. It
  carries the pre-existing pipeline work that was uncommitted on `main` plus the
  matrix contract and the `common/` layer. Commits, oldest first:
  `608563a` enforce container execution and keep remapping enabled,
  `36ff355` Souporcell wrapper argument regression test,
  `ac13a0e` rendered pipeline DAG documentation,
  `803bcf6` group-matrix harness and this handoff document,
  `18bd31a` verified genotype tag choice,
  `fe93051` shared assignment/genotype/comparison contract,
  `ef81c27` deterministic discovery-cell selection,
  `47c4e11` subsets as column views,
  `2234029` per-library baseline conversion,
  `9f004c2` wave B state and the preliminary one-donor signal.
- `validation/group-vireo` -- Vireo discovery/projection adapter, LSF wrappers,
  Gate 2/3/5 runner, merged up to date with the base.
- `validation/group-souporcell` -- Souporcell clustering/troublet/consensus
  adapter with the projection-contract check, same gate structure.

A final integration branch should be created from the selected backend rather
than merging both experimental arms into production.

### Still uncommitted, deliberately

- The nf-demultiplex submodule is dirty: `bin/souporcell.sh` differs from
  submodule commit `36f2d7f` by +36/-20 lines and contains the corrected
  Souporcell boolean handling. It is a separate repository, so it needs its own
  commit there before this repository can record a new pointer. Preserve it; do
  not clean, reset, or reimplement it.
- The repository root holds untracked home-style dotfiles that appear to have
  landed there by accident: `.bashrc`, `.zshrc`, `.profile`, `.bash_profile`,
  `.gitconfig`, `.ripgreprc`, plus `.idea/`, `.vscode/`, `.claude/` and
  `.mcp.json`. They are not part of this work and were left alone; they are
  candidates for `.gitignore` or removal.

## First expanded baseline run: 2026-09-18

The nine-library individual K=2 run completed successfully: 24 tasks completed,
12 tasks were resumed from cache, and there were no failures or retries. This
was an expanded **per-library baseline**, not yet a group-aware Souporcell/Vireo
comparison.

| Group | Library | Singlet split | Fitted ambient | Marker result |
|---|---|---:|---:|---|
| CW4 | Embryo-Imp15289123 | 2,047 / 2,122 | 95.938% | unresolved |
| CW4 | Embryo-Imp15289124 | 6,993 / 6,554 | 83.921% | unresolved |
| CW4 | Embryo-Imp15289125 | 1,605 / 1,578 | 89.255% | unresolved |
| CW4 | Embryo-Imp15289126 | 3,443 / 3,276 | 90.098% | unresolved |
| CW4 | Embryo-Imp15289127 | 3 / 4,151 | 99.815% | unresolved |
| CW5 | Embryo-Imp15350362 | 911 / 5,465 | 2.138% | fetal cluster 0, provisional female |
| CW5 | Embryo-Imp15350363 | 3,277 / 3,124 | 91.578% | unresolved |
| CW5 | Embryo-Imp15350364 | 1,583 / 1,683 | 90.315% | unresolved |
| D6 | FCA7167231 | 47 / 1,193 | 16.272% | unresolved, consistent with one dominant profile |

Interpretation:

- CW5/Embryo-Imp15350362 is the only strong independent two-profile result. Its
  911-cell candidate fetal cluster is supported by five markers (HAPLN1, PITX2,
  XAGE2, XAGE3, and HLA-G), while the 5,465-cell cluster is maternal-like.
- CW5/Embryo-Imp15350363 and 364 produce nearly balanced K=2 splits with about
  90% fitted ambient and no marker separation. Because they share a pregnancy
  with 362, they are ideal held-out libraries for fixed-profile projection.
- Four CW4 libraries produce near-50:50 splits with 84-96% fitted ambient and
  no marker separation. The HLA-G-positive 125 library has HLA-G-positive cells
  in both inferred clusters (22.4% versus 18.9%), supporting an artificial split
  rather than two biologically oriented profiles.
- CW4/Embryo-Imp15289127 and D6 show the opposite failure mode: one dominant
  cluster plus a tiny forced K=2 cluster. D6 is a different donor group and is a
  one-donor negative control, not a matched maternal anchor for CW4 or CW5.
- The fitted ambient estimates must not be interpreted literally when the K=2
  profiles are unstable; here they are primarily a diagnostic of model/profile
  failure.

The CW4 group contains 36,800 barcodes and CW5 contains 19,336, both below
Souporcell's 200,000-barcode sanity limit. The expensive part is nevertheless
BAM processing: individual Souporcell jobs took approximately 3.5-8 hours,
used about 15-16 GB peak RSS each, and read/wrote hundreds of GB. This supports
reusing the completed remapping and allele-count work instead of immediately
running another monolithic pooled-BAM workflow.

## Gate 1 result: 2026-09-19

Group matrices were built from the existing per-library Vartrix outputs under
`<ticket>/group_matrix/`, with no BAM remapping or allele counting repeated.
Every group passed all 20 integrity tests, including the byte-identical rebuild.

| Group | Libraries | Union variants | Cells | Entries per matrix | Checks |
|---|---:|---:|---:|---:|---|
| CW5 | 3 | 597,962 | 19,336 | 45,506,737 | 20/20 |
| CW4 | 5 | 863,707 | 36,800 | 130,623,167 | 20/20 |
| D6 | 1 | 211,925 | 1,501 | 2,534,274 | 20/20 |

Per-library columns extract back out of each group matrix with identical
dimensions, coordinate/value triples after row remapping, and ALT/REF sums.
`variants.vcf.gz`/`.tbi` exist for both replicate builds of every group. Cost is
minutes and a few hundred MB of RSS per group, far below the requested limits,
so the five-library discovery cap is not a matrix-construction constraint.

### Verified environment facts

- The container runtime is the `cellgen/singularity` module (SingularityCE
  4.1.0). `cellgen/apptainer` is non-suid and fails on the execution hosts with
  `Failed to create user namespace`; this is an environment property, not a
  broken image, and no user-namespace workaround was applied.
- Container calls must bind `/lustre` and `/nfs` explicitly. Without the bind,
  tools report a missing file for paths that plainly exist.
- LSF job scripts must source `/etc/profile` and load modules **before**
  `set -euo pipefail`. Sourcing it under `set -u` aborts the job instantly, and
  a redirect on the source hides the reason entirely.
- The queue is `normal`; there is no `small` queue on this farm.

### Pinned-image CLI, captured under LSF

Help output is stored under `<ticket>/group_matrix/tool_help/`. There is **no
`souporcell` executable on `PATH`** -- that name is a directory. The real
binaries are:

- `/opt/souporcell/souporcell/target/release/souporcell` -- reports version
  **2.4**, although the image name says v2.5. Accepts `--alt_matrix`,
  `--ref_matrix`, `--barcodes`, `--num_clusters`, `--seed`, `--restarts`,
  `--threads`, `--initialization_strategy`, `--min_alt/--min_ref`, and for fixed
  profiles `--known_genotypes` plus `--known_genotypes_sample_names`, or
  `--known_cell_assignments`.
- `/opt/souporcell/troublet/target/release/troublet` -- `--alts`, `--refs`,
  `--clusters`, `--doublet_prior`, `--singlet_threshold`, `--doublet_threshold`.
- `/opt/souporcell/consensus.py` -- `-c/--clusters`, `-a/--alt_matrix`,
  `-r/--ref_matrix`, `-v/--vcf`, `--soup_out`, `--vcf_out`, `--output_dir`,
  `-p/--ploidy`.
- `vireo` -- `--vartrixData`, `--donorFile`, `--genoTag` (GT/GP/PL, default PL),
  `--forceLearnGT` (donor GT as prior only), `--nInit`, `--randSeed`, `--nproc`,
  `--noDoublet`, `--extraDonor`, `--ASEmode`, `--cellRange`.
- `cellsnp-lite` 1.2.3, `vartrix`, and `minimap2` are also present.

Consequences for Gate 2 and Gate 3:

- Both backends can consume the harmonised group matrices directly, so the
  head-to-head comparison on identical evidence is feasible as designed.
- `consensus.py` takes the plain `variants.vcf`, which the contract already
  emits, and returns cluster genotypes and an ambient estimate.
- Souporcell does expose a fixed-genotype interface, but `--known_genotypes`
  takes a population VCF with named samples and the doc-level expectation of all
  K profiles still stands. Whether it truly projects without relearning must be
  demonstrated empirically in Gate 3, not assumed from the flag's existence.
- The clustering binary reporting 2.4 inside a v2.5-named image should be
  recorded in the validation report rather than silently accepted.

### Resolved: which `--genoTag` to use

Verified on 2026-09-19 by producing a real donor VCF, not from the manual.

- `vireoSNP.vcf.parse_donor_GPb` accepts exactly `GT`, `GP`, or `PL`.
- Vireo's own discovery output `GT_donors.vireo.vcf.gz` writes `GT:AD:DP:PL`.
- Souporcell's `consensus.py` output `cluster_genotypes.vcf` writes
  `GT:AO:RO:T:E:GO:GN`, with samples named `0`/`1` and `./.` for missing calls.
- **Neither backend writes `GP`.** The earlier preference for "GP when present,
  otherwise validated GT" is therefore not achievable and is superseded:
  - Vireo discovery -> Vireo projection: use `--genoTag PL`, which is both
    Vireo's default and the most informative field available, since it carries
    genotype likelihoods rather than hard calls.
  - Souporcell discovery -> Vireo projection: `--genoTag GT` is the only option.
    `GO`/`GN` are Souporcell-specific and are not parsed by Vireo. Missing `./.`
    calls are mapped to a uniform (1/3, 1/3, 1/3) prior, which is correct
    behaviour but means the missing-call fraction must be reported as a QC
    statistic for any Souporcell-derived profile.
- Neither VCF declares `##FORMAT` header lines. Vireo parses them regardless,
  but any third-party tool added later must be checked against this.

The same job confirmed that Vireo consumes the harmonised group matrices
directly via `--vartrixData=alt.mtx,ref.mtx,barcodes.tsv,variants.vcf.gz`, so
the Gate 2 evidence contract holds end to end for the Vireo arm.

### Revised next validation gate

Build harmonised group matrices from the existing per-library Vartrix outputs:

1. Prefix cell barcodes with the library ID.
2. Align variant rows by chromosome, position, reference, and alternate allele
   across each group's `common_variants_covered.vcf` files.
3. Combine `alt.mtx` and `ref.mtx` columns without repeating remapping or allele
   counting.
4. Run Souporcell clustering and Vireo inference on the same CW5 matrices first.
5. Test whether the profiles anchored by Embryo-Imp15350362 recover coherent
   assignments in 363 and 364, including allowing one-profile libraries.
6. If CW5 passes, repeat group discovery on CW4. CW4 has no trustworthy
   independent anchor, so stability across subsets is mandatory.
7. Use D6 to verify that the selected method does not manufacture a credible
   second genotype in a one-donor library.

This matrix-level comparison is now preferable to developing four complete
pipeline branches or rerunning all BAM preprocessing. It tests Souporcell versus
Vireo on identical remapped allele evidence. Only if Vireo is promising should
one representative original-BAM cellsnp-lite run be used to decide whether the
production workflow can omit remapping.

### Group-matrix technical contract

The validation harness must implement this contract before either backend is
compared.

Input invariants for every library:

- `alt.mtx` and `ref.mtx` are Matrix Market coordinate matrices with variants
  as rows and cells as columns. For example, CW5/Embryo-Imp15350362 is
  288,421 variants by 7,247 cells.
- Matrix row count equals the number of non-header records in
  `common_variants_covered.vcf`; column count equals lines in `barcodes.tsv`.
- ALT and REF dimensions and coordinate pairs must be identical. Preserve the
  Vartrix convention of explicit zero values at observed coordinates; do not
  independently sparsify one matrix.
- Allele counts must be finite, non-negative integers even though the Matrix
  Market header declares `real`.

Merge rules:

1. Define a variant by `(CHROM, POS, REF, ALT)`. Reject multiallelic records,
   duplicate keys within a library, or conflicting alleles at one genomic
   position.
2. Use the **union** of covered variants across the group, ordered according to
   the common source panel
   `/nfs/cellgeni/pipeline-files/demultiplex/reference_vcfs/chr_filtered_2p_1kgenomes_GRCh38.vcf`.
   A locus absent from one library contributes no matrix entries for that
   library; zero depth means missing evidence, not a homozygous-reference call.
3. Concatenate cells in validation-samplesheet order and retain each library's
   original barcode order. Use an internal globally unique ID of
   `<sample_id>::<raw_barcode>` and write an explicit mapping table; downstream
   code must not recover fields by parsing that delimiter.
4. Remap each library's row indices to the group union and offset its column
   indices. Preserve identical ALT/REF coordinate patterns and values.
5. Keep the complete all-cell matrices immutable. Balanced discovery-cell
   subsets are separate column-index views; final projection always uses all
   cells.

Required group output:

```text
<group_id>/
  alt.mtx
  ref.mtx
  variants.vcf
  variants.vcf.gz           # bgzip derivative for tools that require it
  variants.vcf.gz.tbi
  barcodes.tsv
  barcode_provenance.tsv   # group barcode, sample_id, raw barcode, source column
  library_offsets.tsv      # sample_id, first/last group column, n_cells
  matrix_qc.tsv            # dimensions, entries, allele sums, validation status
  manifest.json            # inputs, sizes, mtimes, hashes, tool version, parameters
```

Integrity tests are pass/fail, not warnings:

- output rows equal VCF records and output columns equal barcode/provenance rows;
- barcodes and variant keys are unique;
- ALT and REF coordinate pairs are identical;
- extracting each library's output columns exactly reproduces its input
  dimensions, coordinate/value triples after row remapping, and ALT/REF sums;
- group columns equal the sum of library columns and no source library is lost;
- two runs over the same inputs produce byte-identical matrices and metadata.

The first implementation is `validation/group_matrix/build_group_matrix.py`,
a validation utility rather than a Nextflow process, with
`check_group_matrix.py` as the independent pass/fail checker. Large real-data tests must run through LSF. The harness may use tiny
synthetic fixtures locally to exercise union, missing-locus, duplicate,
conflicting-allele, explicit-zero, and barcode-collision cases.

Backend command shapes, to be finalised from the pinned versions' help output:

```bash
# Genotype-free Vireo on the common matrices
vireo \
  --vartrixData=alt.mtx,ref.mtx,barcodes.tsv,variants.vcf.gz \
  --nDonor 2 \
  --nInit 50 \
  --randSeed <seed> \
  --nproc <cpus> \
  --outDir <output>

# Souporcell clustering on the same matrices
souporcell \
  --alt_matrix alt.mtx \
  --ref_matrix ref.mtx \
  --barcodes barcodes.tsv \
  --num_clusters 2 \
  --seed <seed> \
  --threads <cpus> \
  > clusters_tmp.tsv
```

The Souporcell adapter must then run the pinned Troublet and consensus stages so
the common output includes singlet/doublet status and cluster genotypes. Do not
guess their CLI from a different release. Before the first real-data job,
capture `-h` output for `vireo`, `souporcell`, `troublet`, and `consensus.py`
from the pinned image in an LSF probe job and store it with the validation
report.

An interactive image probe on the current head node failed before tool launch
because the loaded Singularity installation required `FUSE_3.12` and its
non-setuid extraction fallback was unavailable. The completed pipeline proves
the image works on LSF execution hosts. Treat this as a head-node runtime issue,
not a broken container, and do not apply user-namespace workarounds. Preparing
the probe script is safe; obtain user approval immediately before `bsub`.

### Provisional validation acceptance criteria

These thresholds select a backend for further development; they are not yet
production biological thresholds. Any change must be justified in the
validation report rather than adjusted until a tool passes.

- Every backend adapter must emit a common table containing group barcode,
  sample ID, singlet/doublet/unassigned status, profile, assignment probability,
  doublet probability, and informative-locus count. A high-confidence cell has
  normalized assignment probability at least 0.90.
- On CW5/Embryo-Imp15350362, the group result must agree with the strong baseline
  assignment for at least 95% of cells called high-confidence singlets by both
  methods, after permutation-invariant profile matching. The marker-supported
  911-cell profile must remain the fetal-oriented profile; markers are an
  orientation check, not genotype ground truth.
- Across three random initialisations, CW5 cell assignments must have pairwise
  adjusted Rand index at least 0.95 on common high-confidence singlets, and
  inferred profiles must have at least 95% genotype concordance over loci called
  confidently in both runs.
- Projection and joint inference for CW5/363 and 364 must agree for at least 90%
  of common high-confidence singlets after profile matching. A library may
  contain only one profile; no minimum minority fraction is imposed.
- Each accepted profile pair must have at least 100 well-supported loci with
  discordant genotype calls. Fewer loci is insufficient evidence for a locked
  two-profile interpretation.
- D6 must end as `one_profile` or `unresolved`, never a confident two-profile
  biological call. Treat the stress test as failed only if forced K=2 produces
  both a minor high-confidence singlet fraction of at least 5% and at least 100
  well-supported discordant loci. This is a QC stress test, not projection
  against CW4/CW5.
- CW4 must reproduce two profiles across at least three discovery subsets with
  at least 95% genotype concordance and adjusted Rand index at least 0.90 on
  common high-confidence cells. If it cannot, CW4 remains unresolved and must
  not be used to tune the method into passing.
- Report runtime, peak RSS, temporary disk, unassigned fraction, doublet
  fraction, and informative loci per cell. Scientific stability and correct
  one-profile behaviour take precedence over speed.

Stop conditions:

- If neither backend passes CW5, do not start production implementation; audit
  matrix construction, relatedness, coverage, and ambient contamination.
- If only one backend passes, advance only that backend.
- If both pass, prefer the simpler fixed-profile projection path, then use the
  remapping A/B and resource evidence as tie-breakers.

### Discovery and projection contract

Discovery outputs must use neutral profile names (`profile_1`, `profile_2`, ...)
until a maternal reference or group-level marker evidence orients them. The
discovery adapter must emit a donor/profile genotype VCF containing genotype
probabilities or calls at the harmonised variant keys and record which cells
contributed to learning it.

Projection is performed separately for each library against the fixed discovery
profiles:

- use the library's columns from the harmonised group matrices, so variant rows
  and alleles match the discovery VCF exactly;
- do not update or relearn the supplied profiles during the primary projection
  result;
- retain every cell, including unassigned and doublet calls;
- join results back to libraries through `barcode_provenance.tsv`, never by
  splitting the synthetic barcode;
- permit zero high-confidence cells for a profile in a library.

For Vireo, use the discovery `GT_donors.vireo.vcf.gz` as `--donorFile` with the
appropriate existing genotype tag (`GP` preferred when present, otherwise
validated `GT`) and do not enable `--forceLearnGT` for fixed projection. Confirm
the exact pinned CLI in the LSF help probe.

For Souporcell, the pipeline-level `--known_genotypes` mode expects all K
profiles, while documentation for the standalone clustering binary is
inconsistent about fixed-genotype support. The validation adapter must prove
that the pinned 2.5 implementation can project the harmonised matrices without
relearning profiles. If it cannot, Souporcell fails the required projection
contract unless a small, independently validated likelihood-assignment adapter
is justified; do not hide this limitation by rerunning de novo K=2 per library.

## Wave B implementation state: 2026-09-21

Branches, all from the reviewed shared base `validation/group-matrix-core`:

- `validation/group-matrix-core` holds the matrix contract plus the shared
  `common/` layer: the assignment-table schema and its validator, informative
  locus counting from the matrices (so the metric is backend-independent),
  profile-genotype comparison, permutation-invariant run comparison with the
  acceptance thresholds, deterministic discovery-cell selection, column-view
  subsetting, and conversion of the per-library baseline into the common table.
- `validation/group-vireo` adds `backends/vireo_adapter.py` (discovery and
  fixed-profile projection) with LSF wrappers and a Gate 2/3/5 DAG runner.
- `validation/group-souporcell` adds `backends/souporcell_adapter.py`
  (clustering, troublet, consensus) with the same gate structure.

Worktrees live under `/nfs/cellgeni/batu/_repos/worktrees`. Nothing is pushed
and `main` is untouched.

Design points worth keeping:

- A discovery subset is materialised as a real column view, not passed as a
  shorter barcode list. Both backends require one barcode per matrix column, so
  a bare subset list would silently misalign cells. The subset keeps the parent
  variant axis, including uncovered loci, so a profile VCF learned on one subset
  applies to any other without remapping. Verified: subsetting the CW5 group
  matrix back to Embryo-Imp15350362's cells reproduces exactly its original
  7,247 columns and 12,617,422 entries.
- Informative loci per cell are counted from the group matrices, treating
  Vartrix explicit zeros as absence of coverage, rather than taken from each
  backend's own field. On D6 this independently reproduced Vireo's own median
  n_vars of 1,776.
- Souporcell's projection mode tests rather than assumes the fixed-profile
  question: it supplies the discovery VCF through `--known_genotypes` and then
  compares the returned consensus genotypes with the profiles supplied. If they
  moved, it exits 3 and records a projection-contract failure. The gate runner
  waits on `ended(...)` for projections so that finding is still compared and
  reported instead of dropping out of the DAG.
- The souporcell adapter parser reproduces the documented baseline exactly for
  CW5/Embryo-Imp15350362: 911 and 5,465 singlets, 81 doublets, 790 unassigned.

### Where the artefacts are

Under the ticket, `/lustre/scratch127/cellgen/cellgeni/tickets/REQ-68805`:

```text
group_matrix/
  run_wave_a.sh              # submits the Gate 1 builds and the CLI probe
  logs/                      # LSF stdout/stderr for every job
  tool_help/                 # captured -h output, plus the image's binary list
  genotag_probe/vireo_out/   # the donor VCF that resolved --genoTag
  CW5/ CW4/ D6/              # each with build_a/ and build_b/ (determinism)
  gate2_vireo/               # created by the Vireo gate runner
  gate2_souporcell/          # created by the Souporcell gate runner
```

Each `build_*/<group_id>/` holds `alt.mtx`, `ref.mtx`, `variants.vcf`(`.gz`,
`.tbi`), `barcodes.tsv`, `barcode_provenance.tsv`, `library_offsets.tsv`,
`matrix_qc.tsv`, `manifest.json`, and `loci_per_cell.tsv` once an adapter or the
cell selector has computed it.

Observed Gate 1 cost, for planning: D6 (1 library) 28 s and 200 MB peak RSS;
CW5 (3 libraries) and CW4 (5 libraries) a few minutes each, well inside the
16 GB and 32 GB requested. Output size is roughly 650 MB per matrix for CW5 and
1.9 GB for CW4, doubled because two replicate builds are written.

### Preliminary one-donor signal, not yet a Gate 5 verdict

**Premise under revision (2026-09-22).** The description of D6/FCA7167231 as a
single-donor control came from this document and was questioned by the user, who
states it is not a single donor. What the measurements show is that the K=2
split Vireo found in D6 is not reproducible across disjoint halves and that its
profiles explain held-out cells as blends; that stands either way. The
*interpretation* does not: if D6 contains more than one donor, then Souporcell's
51-versus-1,256 result is a missed genotype rather than correct restraint, and
Gate 5 was testing sensitivity, not specificity. Resolve D6's true composition
and tissue before treating any Gate 5 verdict as settled.

The `--genoTag` probe forced Vireo to K=2 on D6, the single-donor control, with
default settings, one seed, and no locus filtering. Applying the Gate 5
criteria to that output gives: minor high-confidence fraction 0.462 (threshold
0.05) and 1,637 well-supported discordant loci (threshold 100). Both conditions
are met, which is the documented definition of a failed one-donor guard.

This must be reproduced through the adapter before it is treated as a result,
because the probe was an interface test: no minimum-coverage or minimum-locus
filtering was applied, only one initialisation was run, and no ambient handling
was configured. It is recorded here because it is the first evidence bearing
directly on the backend choice, and because the contrast is notable: the
per-library Souporcell baseline on D6 produced one dominant cluster plus a tiny
forced second cluster (47 versus 1,193), which is closer to correct
single-profile behaviour than an even split.

## Gates 2, 3 and 5 result: 2026-09-21

Both arms ran to completion on the CW5 group matrices and the D6 guard, 16 jobs
in total, roughly 30-40 minutes wall clock per arm.

### The redesign's core hypothesis is confirmed

Joint group-level inference resolves the artificial per-library splits, in both
backends independently, and correctly leaves a library with one profile:

| Library | Per-library baseline | Souporcell joint | Vireo joint |
|---|---|---|---|
| Embryo-Imp15350362 | 911 / 5,465 | 915 / 5,427 | 854 / 5,466 |
| Embryo-Imp15350363 | 3,277 / 3,124 | 7,383 / 9 | 7,567 / 9 |
| Embryo-Imp15350364 | 1,583 / 1,683 | 3,754 / 3 | 3,848 / 0 |

Souporcell's fitted ambient for the joint run is 1.80% (2.14% for the anchored
run), against 84-96% for the artificial per-library splits. Both backends
reproduce the trusted baseline on Embryo-Imp15350362 at 100% agreement over
common high-confidence singlets.

### Scoreboard against the acceptance criteria

| Criterion | Souporcell | Vireo |
|---|---|---|
| Seed stability (ARI >= 0.95) | pass, ARI 1.0000 | pass, ARI 1.0000 |
| Profile concordance across seeds (>= 0.95) | pass, 0.9999 | pass, 1.0000 |
| Agreement with the per-library baseline on 362 (>= 0.95) | pass, 1.0000 | pass, 1.0000 |
| Anchored discovery vs joint | pass, ARI 1.0000, genotypes 0.9850 | pass, ARI 1.0000, genotypes 0.9978 |
| Projection vs joint (>= 0.90) | **fail, 0.52, ARI ~0** | pass, 1.0000, ARI 1.0000 |
| Profiles held fixed during projection | **fail, 0.726 / 0.747** | pass |
| D6 one-donor guard | pass, minor 0.034 | **fail, minor 0.462, 1,637 discordant loci** |

### Souporcell cannot do infer-once/project-many

Supplying both K=2 profiles through `--known_genotypes` did not hold them: the
returned profiles matched the supplied ones at only 0.726 and 0.747, and the
resulting assignments are unrelated to joint inference (ARI ~0). The
distributions show what happened -- projection re-created the artificial split
(3,332 / 3,033 in 363, 1,620 / 1,667 in 364) exactly where joint inference found
a single profile. The profile-concordance figure alone would be confounded,
because `consensus.py` recomputes genotypes from whichever cells it sees, but
the assignment-level evidence is not. Souporcell's viable architecture here is
joint inference over the whole group, which makes Gate 7 (scale to 10-15
libraries) decisive for it rather than optional.

### What Vireo's D6 failure does and does not mean

Vireo split the single-donor control 694 / 804 with 1,637 well-supported
discordant loci, meeting both conditions for a failed guard, and it reproduced
under the adapter with `nInit 50`, so it is not a probe artefact.

It is a failure of **K selection**, not of assignment quality: in the same run
set, Vireo left CW5/363 and 364 as single-profile libraries when they were part
of a group whose profiles were learned jointly. The D6 scenario corresponds to a
pregnancy in which *no* library carries fetal cells, so the group itself has one
genotype. That is a real risk for this project, but a narrower one than "Vireo
splits single-profile libraries".

Consequence for either backend: a two-profile call must be gated by explicit QC
rather than by trusting the tool to decline. The obvious test, not yet run, is
to compare K=1 against K=2 on D6 and check whether the model actually prefers
the simpler explanation.

## Is Vireo's D6 split real? 2026-09-22

Two tests, with CW5 as a positive control so the numbers are interpretable. The
K sweep originally planned alongside them was dropped as decision-relevant: a
pregnancy group always contains a mother, so K=1 is not a production scenario,
and the question "is the second profile justified?" is answered directly by the
tests below without any penalised-likelihood assumptions.

### Split-half locus replication: the split does not replicate

Learn K=2 independently on two disjoint halves of the cells, then ask whether
the *same loci* separate the profiles in both halves. Genotype concordance
between halves is **not** the right statistic -- with one donor both profiles are
near-copies, so any pairing looks concordant (D6 scores 0.947 there).

| Group | Separating loci, half A | half B | in both | Jaccard | Fraction of A replicated |
|---|---:|---:|---:|---:|---:|
| CW5 (positive control) | 11,549 | 11,581 | 10,774 | **0.872** | **0.933** |
| D6 (negative control) | 453 | 479 | 114 | **0.139** | **0.252** |

A real second genotype replicates on 93% of its separating loci; D6's split
replicates on 25%. The statistic separates the two cases cleanly and is the
basis for a usable Gate 5 procedure.

Do not use the "enrichment over chance" column as the criterion: chance overlap
scales with set size, so D6 scores 8.6x and CW5 only 4.7x, inverting the
ranking. Jaccard and the replicated fraction are the metrics to keep.

### Held-out projection: consistent with the same conclusion

Projecting each group's half B against half A's fixed profiles:

| Group | Confident singlets | Doublet fraction | Unassigned |
|---|---|---:|---:|
| CW5 | 8,832 (117 / 8,715) | 0.038 | 0.049 |
| D6 | 12 (4 / 8) | **0.936** | 0.048 |

D6's profiles explain almost no held-out cell as a singlet: 93.6% look like
blends, which is what two near-copies of one genotype produce. Vireo's D6 split
is therefore an artefact by two independent measures, which is the better
outcome for Vireo -- the one-donor guard becomes a QC step rather than a
disqualification.

Also recorded: Vireo 0.5.8 **cannot fit K=1 at all**. `--nDonor 1` crashes in
`vireo_doublet.add_doublet_GT` (indexing a 1-D array as 2-D). Protection against
a manufactured second profile can never come from the tool declining, so it must
be external QC regardless of which backend is chosen.

### The positive control exposed a projection defect

Among the half-B cells that joint inference calls fetal, fixed-profile
projection recovers only a minority as singlets:

| Fate under projection of the 456 joint-called fetal cells | Count |
|---|---:|
| singlet on the fetal profile | 114 |
| **doublet** | **328** |
| unassigned | 14 |

Projection converts about 72% of the fetal population into doublets. For
CW5/Embryo-Imp15350362 it calls 116 fetal cells where joint inference calls 447.

Critically, the Gate 3 criterion did not notice: agreement on common
high-confidence singlets is **1.0000**, because cells that projection demotes to
doublet stop being "common high-confidence singlets" and silently leave the
comparison. The earlier Gate 3 pass for Vireo on 363 and 364 is therefore weaker
evidence than it appeared -- those libraries have no minority population to lose,
so the metric could not detect this failure mode there.

Consequences:

- The acceptance criteria need a companion metric to cell agreement: retention
  of the minority profile, defined as the fraction of joint-called minority
  cells that projection keeps as minority singlets, reported alongside the
  doublet-rate inflation. Agreement alone can be 1.0 while the fetal population
  is discarded, which is precisely the biology this pipeline exists to measure.
- A fair test of fetal retention needs the profile learned from a full minority
  population, not half of it: here half A held only ~430 fetal cells. The clean
  follow-up is a split-half restricted to Embryo-Imp15350362, the one library
  with a real minority population, learning on one half and projecting the other.

## Fetal retention under fixed-profile projection: 2026-09-22

Restricted to Embryo-Imp15350362, the only library with a genuine minority
population, varying how much of that population the profile is learned from.
Half B (3,623 cells) is projected in both arms; the reference is joint CW5
inference restricted to the same cells.

| Profiles learned on | Fetal cells in training | Retention | Fetal cells lost to doublet | Doublet fraction |
|---|---:|---:|---:|---|
| all of 362 (7,247 cells) | ~915 | **0.9955** (445/447) | 2 | 0.011 -> 0.008 |
| half A only (3,624 cells) | ~458 | **0.3177** (142/447) | 289 | 0.011 -> 0.085 |

Both arms score cell agreement **1.0000**, which is exactly why the retention
metric was needed.

The conclusion is not that fixed-profile projection is broken. With a
well-learned minority profile it is close to perfect. What breaks it is thin
training support for the minority: halving the fetal cells available to
discovery (about 915 to about 458) collapses retention from 99.6% to 31.8%, and
the lost cells become doublets rather than unassigned, because a noisy minority
profile makes a genuine fetal cell look like a blend of the two profiles.

Caveat on the first arm: half B's cells are inside the full-362 training set, so
0.9955 is an upper bound rather than a held-out result. The clean design for
separating "more training cells" from "training overlap" is a training-size
curve -- learn on 25%, 50% and 75% of the library and project the same held-out
25% each time. The 31.8% figure is the honest held-out number at half strength,
and it agrees with the 25% (114/456) seen in the earlier CW5 split-half, which
also trained at half strength.

Consequences for the design:

- Discovery must include enough of the minority population. The group-aware
  architecture already does this when the fetal-rich library is part of
  discovery, which is the intended behaviour; it argues against discovery
  subsets that thin out the fetal evidence, and so against a balanced per-library
  cell cap that would sample down the one library carrying fetal cells.
- **Doublet inflation is a usable production alarm.** A projection that raises
  the doublet fraction (here 0.011 to 0.085) is signalling an under-trained
  profile, not finding real doublets. Report it alongside retention.
- Vireo's Gate 3 pass stands, with the qualification that it holds only when the
  minority profile has adequate support.

## Backend head-to-head on CW5: 2026-09-23

The Gate 2 comparison neither arm runner could do, since neither can depend on
the other's jobs. Souporcell joint seed 1 against Vireo joint seed 1, same
matrices, same cells.

| Measure | Value |
|---|---|
| Common high-confidence singlets | 17,209 |
| Cell agreement / ARI | 1.0000 / 1.0000 |
| Genotype concordance | 0.9854 over 239,228 loci |
| Minority retention (Souporcell as reference) | 0.8727 (809 of 927) |
| Doublet fraction | 0.0069 -> 0.0041 |

**On CW5 the two backends give the same answer.** Of the 118 fetal cells they do
not share, 112 are unassigned by Vireo rather than assigned differently, and six
are called doublet. No cell is assigned to the opposite profile.

This shapes how the backend gets chosen: accuracy on CW5 does not separate them,
so the decision rests on the architectural differences already measured --
Souporcell cannot hold profiles fixed and must re-infer jointly over the whole
group, while Vireo projects well provided the minority profile has adequate
training support.

## Gate 4 result: CW4 stays unresolved, 2026-09-23

CW4 discovery on all five libraries plus four subsets chosen so no library
appears in every trio (`trio_early` 123/124/125, `trio_late` 125/126/127,
`trio_mixed` 123/125/127, `balanced2k` = 2,000 most informative cells per
library), in both backends, with projection of all CW4 cells against the
`trio_early` profiles. Reports: `<ticket>/group_matrix/gate4_{vireo,souporcell}/report/`.

Acceptance criterion: two profiles reproduced over at least three subsets with
genotype concordance >= 0.95 and ARI >= 0.90 on common high-confidence cells.
**Neither backend meets it, so CW4 remains unresolved, as the criteria require.**

### Souporcell: one profile in four of five subsets

| Subset | Result | Discordant loci between profiles |
|---|---|---:|
| `all5` | one profile in every library (1-22 minority cells each) | 0 |
| `trio_late`, `trio_mixed` | one profile | 0 |
| `balanced2k` | one profile, every cell | - |
| `trio_early` | ~50:50 in each library, 90.9% fitted ambient | 16,660 |

Four subsets collapse to a single genotype, and the second cluster is a copy of
the first. Fitted ambient of ~99.9% in those runs is the signature of that
collapse, not a measurement. The fifth subset reproduces the per-library
baseline's artificial split. Fixed-profile projection failed its contract again
(profile concordance 0.0), so the CW5 finding that Souporcell cannot hold
profiles fixed reproduces on the hard group.

### Vireo: even splits that do not replicate

Every Vireo subset splits every library roughly evenly, but:

| Measure | CW4 | CW5 (real genotype) | D6 |
|---|---|---|---|
| High-confidence fraction, full group | **0.20** (7,294 / 36,800) | 0.92 | - |
| Unassigned, full group | 0.67 | - | - |
| Separating-locus replication (Jaccard) | **0.19-0.28** | 0.872 | 0.139 |
| Cross-subset ARI | **0.008-0.124** | 1.000 | - |
| Cross-subset cell agreement | 0.55-0.68 (chance ~0.5) | 1.000 | - |

Genotype concordance across subsets is 0.988-0.990, above the 0.95 threshold,
**and must not be read as a pass**. This is the same trap found on D6: when both
profiles are near-copies of one dominant genotype, any pairing across runs looks
concordant. The discriminating statistics -- locus replication, ARI and the
confidence fraction -- all fail. Projection against `trio_early` profiles calls
42% of cells doublets and assigns almost nothing in libraries 126 and 127.

### What this does and does not show

Two independent methods, across five discovery subsets, find no reproducible
second genotype in CW4's evidence. That is a finding about what is *separable*
from these allele matrices, not proof that CW4 contains no fetal cells: a low
fetal fraction, high ambient contamination or insufficient per-cell coverage
would produce the same result. The per-library baseline already showed
HLA-G-positive cells distributed across both artificial clusters of
Embryo-Imp15289125 (22.4% versus 18.9%), which is consistent with either reading.

Consequences:

- CW4 does not discriminate between the backends: both fail it, for their
  different characteristic reasons.
- The acceptance criterion needs amending: genotype concordance across runs
  should not count toward Gate 4 on its own. Locus replication and the
  high-confidence fraction should be required alongside it.
- A composition QC gate now has three calibration points: CW5 accepts (Jaccard
  0.87, 92% confident), and D6 and CW4 do not (Jaccard 0.14 and 0.19-0.28, CW4
  20% confident). D6's own interpretation stays open until its composition is
  confirmed.
- The open scientific question is sensitivity: what fetal fraction and coverage
  would be needed for a second genotype to become separable. A detection-limit
  experiment is the way to answer it, since neither CW4 nor D6 can.

## Detection-limit experiment: design, 2026-09-23

Neither CW4 nor D6 can say whether "no reproducible second genotype" means no
fetal cells or a fetal population below the detection floor. This experiment
measures the floor directly, using only CW5 data both backends agree on.

- **Truth** is consensus, not ground truth: a cell is fetal or maternal only
  when Souporcell and Vireo joint CW5 runs both call it so at high confidence.
  Pool: 809 fetal (804 from Embryo-Imp15350362, 5 from 363) and 16,400
  maternal, with **zero** confident disagreements between the backends.
- **Background**: 7,000 consensus-maternal cells drawn from all three CW5
  libraries. Because the fetal cells come almost entirely from 362, drawing the
  background from 362 as well stops library identity from predicting genotype.
- **Spike**: consensus-fetal cells at 0, 0.5, 1, 2, 5 and 10% of the group (0 to
  778 cells), two replicates each. Mother and fetus are genuinely related, which
  is the realistic hard case. At 10% the replicates share 778 of 809 fetal cells,
  so replicate variance there is understated.
- **Per group**: Vireo discovery on the whole group and on two disjoint halves
  (composition gate), and Souporcell discovery on the whole group; 60 jobs plus
  a scoring job.
- **Scored separately**: detection (recall and precision of the spiked cells,
  and false positives in the 0% control) and the composition gate
  (separating-locus replication between halves, confident-cell fraction). A
  backend can detect the fetus while the gate rejects it, and the experiment is
  built to show where that happens.

Runner: `validation/group_matrix/lsf/run_spikein.sh` on the shared base, driving
both adapters through their worktrees. It was dry-run end to end with a stub
`bsub` before submission, and nested subsetting was verified exactly on D6.
Output: `<ticket>/group_matrix/spikein/report/detection_curve.tsv`.

## Detection-limit result, 2026-09-24

All 60 jobs completed; `<ticket>/group_matrix/spikein/report/detection_curve.tsv`.
Recall and precision are of the spiked consensus-fetal cells, from K=2 discovery
on the whole group. The gate is Vireo split-half locus replication (Jaccard).

| Fetal fraction (cells) | Souporcell recall / precision | Vireo recall / precision | Vireo gate Jaccard |
|---|---|---|---|
| 0% (0) | false split ~50:50 | false split ~50:50 | 0.08-0.09 |
| 0.5% (35) | **1.00 / 0.83-0.85** | 0.51-0.57 / 0.01 | 0.08 |
| 1% (71) | 0.99 / 0.89-0.99 | 0.46-0.69 / 0.01-0.02 | 0.09 |
| 2% (143) | 0.99-1.00 / 0.96-1.00 | 0.73-0.76 / 0.06 | 0.09-0.10 |
| 5% (368) | 0.64-0.69 / 1.00 | **1.00 / 1.00** | 0.38-0.43 |
| 10% (778) | 0.97 / 1.00 | 1.00 / 1.00 | **0.87** |

### Findings

1. **Souporcell is roughly ten times more sensitive in discovery.** It isolates
   the fetal genotype from 35 cells in 7,000 (0.5%) as a profile of 41-42 cells.
   Vireo's floor lies between 2% and 5%. Below it, Vireo still splits the
   group, but along an arbitrary line: precision 0.01-0.06 means the fetal cells
   are scattered across both halves, not found.
2. **Neither backend declines when there is nothing to find.** In the 0%
   control both split ~50:50, and both profiles draw evenly from all three
   libraries, so the split is not a batch or ambient artefact but an arbitrary
   halving of one genotype. This corrects the earlier reading of Souporcell's
   D6 and CW4 behaviour as restraint: on a clean single-genotype background it
   manufactures a split too. The composition gate is required for **both**
   backends.
3. **The Vireo gate tracks whether the found split is real.** Jaccard sits at
   0.08-0.10 wherever Vireo's split is not the fetal genotype (0-2%), rises to
   0.38-0.43 at 5%, and reaches 0.87 at 10%, matching CW5's 0.872. The gate's
   discriminating power depends on how many minority cells each half holds
   (~184 at 5%, ~389 at 10%), consistent with the retention finding.
4. **A Souporcell anomaly at 5%.** Recall drops to 0.64-0.69 there, with the
   lost fetal cells called doublets (110-129), while 2% and 10% are near
   perfect. Unexplained; worth a targeted look before relying on Souporcell
   around that fraction.

### What this means for earlier results

- **CW4.** Souporcell detects a related fetal genotype at 0.5% in CW5 data, yet
  on CW4 it found no genotype difference at all (0 discordant loci between its
  profiles, minority clusters of 1-22 cells). If CW4 carries fetal cells, they
  are most likely below ~0.5% of cells or undetectable for reasons this
  experiment does not model: a different pregnancy, coverage or ambient level.
  The floor was measured on CW5 cells and may not transfer exactly.
- **Gate thresholds.** For Vireo halves of ~3,500 cells, a Jaccard threshold
  between 0.10 and 0.38 separates real from arbitrary splits. CW4's trio
  replication (0.19-0.28) falls inside that gap, but it came from a different
  design (overlapping trios of 14,000-25,000 cells), so it cannot be classified
  with this calibration.
- **Backend choice.** The picture is now split by stage: Souporcell is the
  better discoverer and cannot project; Vireo projects well but discovers only
  at >=5%. That points to the hybrid already named in this document: Souporcell
  discovery feeding Vireo fixed-profile projection with `--genoTag GT`. It is
  untested, and it is the next experiment.

### Follow-up submitted, 2026-09-28

`validation/group_matrix/lsf/run_spikein_followup.sh`, on the same 12 groups
(40 jobs plus scoring): Vireo projection of each group against its Souporcell
discovery profiles with `--genoTag GT` (the hybrid); Souporcell discovery on
each half (the Souporcell-side gate); and Souporcell seeds 2 and 3 on the 5%
groups, to test whether the recall dip is a local optimum. The Vireo adapter
was changed to accept Souporcell's `0`/`1` profile names and its
`cluster_to_profile` map, which the hybrid needs; this path has not run before.

### Gaps

- No Souporcell-side gate was run: split halves were Vireo only. A Souporcell
  split-half at the same fractions is needed before a gate can cover the hybrid.
- Two replicates per fraction; at 10% they share 778 of 809 fetal cells.

## Follow-up result: the hybrid works, 2026-09-29

All 40 jobs completed. Hybrid = Souporcell discovery profiles, projected by
Vireo with `--genoTag GT`. Recall / precision of the spiked fetal cells, two
replicates each:

| Fetal fraction | Souporcell alone | **Hybrid** | Vireo alone | Hybrid confident fraction |
|---|---|---|---|---|
| 0% | false ~50:50 split | **no confident profile pair** (19-23 cells) | false ~50:50 split | **0.007** |
| 0.5% (35 cells) | 1.00 / 0.83-0.85 | 0.00-0.02 (one profile) | 0.51-0.57 / 0.01 | 0.94-0.96 |
| 1% (71) | 0.98 / 0.88-0.98 | **0.94-0.98 / 0.93-1.00** | 0.46-0.69 / 0.01-0.02 | 0.99 |
| 2% (143) | 0.99-1.00 / 0.95-0.99 | **0.98-0.99 / 0.98-1.00** | 0.73-0.76 / 0.05 | 0.99 |
| 5% (368) | 0.64-0.69 / 1.00 | 0.64-0.68 / 1.00 | 1.00 / 1.00 | 0.97 |
| 10% (778) | 0.97 / 1.00 | **0.97 / 1.00** | 1.00 / 1.00 | 0.98-0.99 |

### Findings

1. **The hybrid detects down to 1% of cells** (~70 fetal cells in 7,000) with
   precision 0.93-1.00, against >=5% for Vireo alone. It keeps almost all of
   Souporcell's sensitivity while using a projection step that holds profiles
   fixed.
2. **It is the only configuration that does not invent a second genotype.**
   With no fetal cells, both backends alone split the group ~50:50, but the
   hybrid assigns almost nobody confidently (0.7% of cells), because Souporcell's
   two arbitrary halves of one genotype give Vireo nothing to separate. When a
   real fetal genotype exists the confident fraction is 0.94-0.99. That contrast
   is a built-in composition check.
3. **At 0.5% it fails conservatively.** It reports one profile rather than a
   wrong one: the fetal profile learned from 35 cells is too weak, and the fetal
   cells are called doublets. A missed fetus, not a false one.
4. **The split-half gate does not work on Souporcell profiles.** Its Jaccard is
   0.75-0.93 at every fraction, **including 0%**, so it cannot tell a real
   genotype from an arbitrary split there. The statistic works on Vireo's
   likelihood-based profiles (0.08 versus 0.86) but not on Souporcell's GT
   calls. The hybrid's gate must therefore use its own projection signal: the
   confident fraction, the minority-profile size, and doublet inflation.
5. **The 5% dip is real, not a local optimum.** Seeds 2 and 3 reproduce it
   exactly (recall 0.63-0.64 and 0.69), and the hybrid inherits it with a
   different assignment engine, with the lost fetal cells again called doublets
   (111-130). The defect is therefore in Souporcell's discovered fetal profile
   for those groups, not in clustering randomness or troublet. Precision stays
   1.00, so it loses fetal cells rather than mislabelling maternal ones. Doublet
   inflation flags it. Cause unexplained.

### Decision

**Adopt the hybrid:** Souporcell group discovery, Vireo fixed-profile
projection with `--genoTag GT`, and a composition gate on the projection's
confident fraction, minority-profile size and doublet inflation. Floor ~1% fetal
cells; below it, report one profile or unresolved. This meets the stopping rule
agreed on 2026-09-28: the backend is chosen and the gate is set from data, so
validation ends here and further runs are acceptance tests of the production
implementation.

Thresholds were measured on 7,000-cell groups from one pregnancy and should be
confirmed on the first real production groups rather than treated as final.

## Operational notes for running this harness on the farm

Lessons that each cost a failed or misread run. They apply to any new runner.

**Submitting and watching**

- Runs are submitted as `cellgeni-su`, a shared service account. `bjobs` as
  another user shows nothing; use `bjobs -u cellgeni-su -w`.
- `bjobs` occasionally returns empty output. A watcher must treat an empty reply
  as unknown, and stop only when the `JOBID` header is printed with no matching
  jobs, or on `No unfinished job found`.
- `bhist -l` ends every unfinished job with a "Summary of time ... by <date>"
  line stamped with the time of the query. It is not a termination record; the
  job may still be running.
- A running job's log is often silent for its whole runtime: the adapters print
  nothing until the backend finishes. A stale log is not evidence of death.
- `TERM_OWNER` with exit 130 is a `bkill` from the owning account (it sends
  SIGINT first), not a resource or queue limit. Job 798449 (Vireo on full CW4)
  was killed this way at 38 GB of 64 GB and 71 minutes in; who issued it could
  not be recovered. On a shared account, tell co-users before long runs.
- Resume on output, never on job IDs: a finished job leaves the queue, and a
  dependency on it can become unsatisfiable. All runners now check for their
  result files and submit only what is missing.

**Writing runners**

- `bash -n` parses but never executes, so it misses runtime word-expansion bugs.
  Dry-run the submission path with a stub `bsub` that logs each job's
  arguments and environment (`BSUB=<stub> bash run_*.sh`), as done for the
  detection-limit runner.
- `VAR=x "$@" bsub ...` does not work: bash parses assignments before expansion,
  so an expanded `MODE=libraries` is treated as a command name. Use
  `env VAR=x "$@" bsub ...`.
- Pass a value that starts with `--` using `=`, as in
  `--extra-args=--noDoublet`. Otherwise argparse reads it as a new flag.
- Heavy data work (subsetting CW4 streams ~3.8 GB per matrix) goes in LSF jobs,
  not on the submitting host.
- Caches written by concurrent jobs must be written atomically (temporary file
  plus rename); `loci_per_cell.tsv` is.

**Permissions and paths**

- Output directories are group-writable with ACLs; mixed `bc8` / `cellgeni-su`
  ownership under the ticket is expected and harmless.
- Container calls need `cellgen/singularity` and `-B /lustre,/nfs`; see the
  verified environment facts under Gate 1.

## Evidence and implementation notes

- Souporcell performs remapping, variant calling, allele counting, clustering,
  doublet calling, and consensus/ambient estimation in one pipeline:
  <https://github.com/wheaton5/souporcell>.
- Vireo supports no genotype reference, complete donor genotypes, and genotype
  data for only a subset of donors:
  <https://github.com/single-cell-genetics/vireo/blob/master/doc/manual.rst>.
- cellsnp-lite supports targeted pileup from a BAM and barcode list using a
  common variant VCF:
  <https://cellsnp-lite.readthedocs.io/en/latest/main/manual.html>.
- Current pinned image:
  `souporcell_v2.5-python_3.10-vireo_v0.5.8-cellsnp-lite_1.2.3-nextflow-compatible.sif`.

## Session log

- 2026-09-18: agreed that five is a provisional Souporcell discovery batch
  limit rather than a biological/group-size limit.
- 2026-09-18: proposed infer-once/project-many architecture for 10-15 libraries.
- 2026-09-18: identified Vireo partial-genotype support as a strong fit for
  maternal blood/reference libraries.
- 2026-09-18: kept the backend choice open pending CW4/CW5/D6 validation.
- 2026-09-18: clarified that remapping targets SNP-alignment bias, not general
  Cell Ranger mapping quality, and is a one-time validation decision.
- 2026-09-18: replaced four parallel implementation branches with staged gates
  that reuse Souporcell allele matrices and implement only the selected design.
- 2026-09-18: initially paused architecture-specific development pending the
  expanded baseline; that hold was superseded after the run completed. Future
  experimental branches will share a common contract, and only the
  evidence-supported winner will become the maintained workflow.
- 2026-09-18: reviewed the completed nine-library individual baseline. Only
  CW5/Embryo-Imp15350362 produced a convincing two-profile result; the remaining
  mixtures confirmed the need for group-aware inference.
- 2026-09-18: changed the next validation step to harmonise and combine existing
  Vartrix matrices, compare Souporcell and Vireo on identical evidence, and
  avoid repeating expensive BAM remapping.
- 2026-09-18: audited this document for cold-start handoff, removed superseded
  validation state, recorded exact run/repository provenance, specified the
  group-matrix contract and experimental branches, and added provisional
  pass/fail criteria.
- 2026-09-19: implemented the group-matrix harness under
  `validation/group_matrix/` (builder, independent integrity/determinism
  checker, synthetic fixtures, validation groups samplesheet, LSF wrappers for
  the CW5 build and the pinned-image tool help probe). Implemented in the
  working tree on `main`; no experimental branch or worktree was created,
  because capturing the existing dirty changes is still a user decision.
- 2026-09-19: ran Gate 1 on LSF for CW5, CW4 and D6. All three passed every
  integrity and determinism test; results and the captured pinned-image CLI are
  recorded in the Gate 1 section above.
- 2026-09-29: follow-up complete. The hybrid detects down to ~1% with precision
  0.93-1.00 and assigns almost no cell confidently when there is no fetal
  genotype; the split-half gate fails on Souporcell profiles; the 5% dip is
  seed-independent. Adopted the hybrid and closed validation.
- 2026-09-24: detection-limit results. Souporcell discovers at 0.5%, Vireo only
  at >=5%; both split a fetal-free control ~50:50 across all libraries, so the
  gate is needed for both. Leaning moved to Souporcell discovery plus Vireo
  projection.
- 2026-09-23: designed and submitted the fetal-fraction detection-limit
  experiment; refreshed the handoff and recorded the farm operational lessons.
- 2026-09-23: ran Gate 4 on CW4 in both backends. Neither reproduces two
  profiles across subsets (Souporcell collapses to one profile in four of five;
  Vireo's splits replicate at Jaccard 0.19-0.28 with ARI near zero), so CW4
  remains unresolved. Souporcell's projection failure reproduced. Genotype
  concordance passed spuriously and should not count toward Gate 4 alone.
- 2026-09-23: ran the cross-arm head-to-head on CW5. The backends agree exactly
  (ARI 1.0000, genotype concordance 0.9854), so CW5 accuracy does not separate
  them and the choice rests on architecture.
- 2026-09-22: established that fixed-profile projection retains 99.6% of fetal
  cells when the profile is learned from the full population but only 31.8% at
  half strength, with the lost cells becoming doublets. Added a
  minority-retention criterion, since cell agreement scored 1.0000 in both arms.
- 2026-09-22: showed Vireo's D6 split is not reproducible (locus replication 0.139
  against 0.872 for CW5; 93.6% of held-out cells called doublet), dropped the
  K sweep as decision-relevant, and found that fixed-profile projection converts
  about 72% of joint-called fetal cells into doublets while still scoring 1.0000
  on the Gate 3 criterion.
- 2026-09-21: ran Gates 2, 3 and 5 for both backends on CW5 and D6. Group-level
  inference resolved the artificial per-library splits in both backends.
  Souporcell failed fixed-profile projection; Vireo failed the one-donor guard.
  Neither backend passes every criterion, so the stop conditions do not yet
  select a winner.
- 2026-09-21: implemented the shared contract layer and both backend adapters
  on their branches, with self-tests passing and every parser verified against
  real backend output. Gate 2/3/5 DAGs are prepared but not submitted.
- 2026-09-19: resolved the `--genoTag` question empirically (PL within Vireo,
  GT across backends; GP is written by neither) and committed the harness on
  branch `validation/group-matrix-core` as the shared base for the backend
  branches.
- 2026-09-19: dry-ran the variant-key contract over the real CW5 metadata
  without touching the matrices. All three libraries are internally consistent
  (VCF records equal matrix rows, barcodes equal columns, ALT and REF share
  dimensions) and no duplicate, multiallelic, or cross-library conflicting
  alleles were found. The CW5 union is 597,962 variants x 19,336 cells with
  45,506,737 entries per matrix; the VCF stage peaked at about 720 MB RSS.
