# Group-matrix validation harness (Gate 1)

Validation-only utilities for the group-aware redesign described in
[docs/group_aware_redesign.md](../../docs/group_aware_redesign.md). Nothing here
is part of the public Nextflow workflow, and nothing here re-runs BAM remapping
or allele counting: it reuses the per-library Souporcell/Vartrix outputs that
the completed baseline already produced.

## What it does

`build_group_matrix.py` harmonises one group's per-library
`alt.mtx`/`ref.mtx`/`barcodes.tsv`/`common_variants_covered.vcf` into the group
matrices required by the technical contract:

```text
<outdir>/<group_id>/
  alt.mtx  ref.mtx            # variants x cells, union rows, concatenated columns
  variants.vcf                # union of covered variants, reference-panel order
  variants.vcf.gz(.tbi)       # bgzip derivative, for tools that require it
  barcodes.tsv                # <sample_id>::<raw_barcode>
  barcode_provenance.tsv      # the mapping table; never parse the :: delimiter
  library_offsets.tsv         # per-library column block
  matrix_qc.tsv               # dimensions, entries, allele sums, status
  manifest.json               # inputs, sizes, mtimes, sha256, parameters
```

Contract points it enforces, all as hard failures (exit status 2):

- a variant is `(CHROM, POS, REF, ALT)`; multiallelic records, duplicate keys
  within a library, and conflicting alleles at one position are rejected;
- rows are the **union** of covered variants — a locus absent from a library
  contributes no entries there, because zero depth is missing evidence, not a
  homozygous-reference call;
- variant order follows the reference panel's contig order then coordinate, so
  the group rows are a subsequence of the panel;
- matrix dimensions must equal the VCF record count and barcode count;
- ALT and REF must share identical coordinate pairs, and Vartrix's explicit
  zeros are preserved rather than sparsified away;
- allele counts must be finite non-negative integers despite the `real` banner.

`check_group_matrix.py` re-derives every claim from the inputs and outputs
independently, including extracting each library's columns back out of the group
matrices and comparing dimensions, coordinate/value triples after row remapping,
and ALT/REF sums. It exits 0 on all-pass, 1 on any failure, 2 if it could not
run. `--compare-dir` adds the byte-identical-rebuild test.

## Run the self-test (safe on a head node)

```bash
bash tests/run_tests.sh
```

Tiny synthetic fixtures cover union, missing loci, explicit zeros, barcode
collisions between libraries, duplicate keys, conflicting alleles, multiallelic
records, dimension mismatch, ALT/REF coordinate mismatch, negative counts, and
an off-panel contig — plus a tampered output that the checker must catch.

## Run real data (LSF only)

The real matrices are hundreds of MB per library, so builds run under LSF:

```bash
source /etc/profile >/dev/null 2>&1 || true; source ~/.bashrc >/dev/null 2>&1 || true
GROUP_ID=CW5 OUTDIR=/lustre/scratch127/cellgen/cellgeni/tickets/REQ-68805/group_matrix \
  bsub -J group_matrix_CW5 < lsf/build_group.bsub
```

`lsf/probe_tool_help.bsub` captures `-h` for `vireo`, `souporcell`, `troublet`,
`consensus.py`, and `cellsnp-lite` from the pinned image. Run it before writing
any backend adapter; the head node cannot launch that image (its Singularity
requires `FUSE_3.12`), which is a head-node runtime issue, not a broken image.

`groups_validation.csv` defines CW4, CW5, and D6 over the REQ-68805 baseline
outputs. D6 is recorded with `expected_genotypes=1`: it is the one-donor
negative control, and its forced-K=2 stress test is a backend-level override,
not a claim about the group.

No maternal reference exists for CW4 or CW5 in this validation set, so every row
is `sample_role=mixture`. The `maternal_reference` role is accepted by the
harness and exercised by the fixtures, ready for a matched reference later.

## Order of work

1. self-test (above);
2. `probe_tool_help.bsub`;
3. CW5 build + integrity check — this is Gate 1;
4. only then Gate 2, the Souporcell/Vireo comparison on identical evidence.
