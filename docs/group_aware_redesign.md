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

No production redesign has been implemented. The present pipeline is still the
per-library Souporcell K=2 baseline with downstream marker labelling.

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

The validation-only group-matrix harness now exists at
`validation/group_matrix/` (see its README). It implements the group-matrix
technical contract below, passes its synthetic self-test, and has LSF wrappers
prepared but **not submitted**. The immediate task is therefore to run the tool
help probe and the CW5 build under LSF, then use that evidence to decide whether
Souporcell or Vireo merits a production implementation. The public Nextflow
workflow is still untouched, as intended.

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

As recorded on 2026-09-18, the repository is on `main` at `2fc50e4` with
uncommitted work. Modified paths are `ROADMAP.md`, `docs/testing.md`, `main.nf`,
`modules/local/souporcell.nf`, `nextflow.config`, and the nf-demultiplex
submodule. Untracked paths include `docs/dag.md`, `docs/dag.png`, `docs/dag.svg`,
this document, and `tests/test_souporcell_wrapper.sh`. Inside nf-demultiplex,
`bin/souporcell.sh` is modified relative to submodule commit `36f2d7f`.

These changes include the corrected Souporcell boolean handling and belong to
the current work. Preserve and inspect them before creating branches; do not
clean, reset, or silently reimplement them.

No experimental redesign branch has been created. When implementation begins,
use worktrees below `/nfs/cellgeni/batu/_repos/worktrees` and share one reviewed
matrix-contract base. Proposed short-lived branch roles are:

- `validation/group-matrix-core`: harmonisation, integrity checks, common QC,
  and fixtures only;
- `validation/group-souporcell`: Souporcell clustering/consensus adapter on top
  of the shared contract;
- `validation/group-vireo`: Vireo discovery/projection adapter on top of the
  same contract;
- a final integration branch created from the selected backend, rather than
  merging both experimental workflows into production.

Before creating worktrees, decide with the user how the current dirty changes
will be captured. Do not create branches from an assumed clean `main`.

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
