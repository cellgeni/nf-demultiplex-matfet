#!/usr/bin/env python3
"""Pass/fail integrity tests for a built group matrix.

Re-derives every claim from the inputs and the outputs, independently of the
build code paths. Any violation is a failure, never a warning. Exit status 0
means every test passed; 1 means at least one test failed; 2 means the check
could not run.
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from build_group_matrix import (  # noqa: E402
    ContractError,
    open_mtx_body,
    read_barcodes,
    read_contig_order,
    read_groups_csv,
    read_vcf,
    sha256,
)

DETERMINISM_FILES = (
    "alt.mtx",
    "ref.mtx",
    "variants.vcf",
    "barcodes.tsv",
    "barcode_provenance.tsv",
    "library_offsets.tsv",
    "matrix_qc.tsv",
)


class Report:
    def __init__(self) -> None:
        self.failures: list[str] = []
        self.passes: list[str] = []

    def check(self, name: str, condition: bool, detail: str = "") -> bool:
        if condition:
            self.passes.append(name)
            print(f"PASS  {name}")
        else:
            self.failures.append(f"{name}: {detail}")
            print(f"FAIL  {name}: {detail}")
        return condition


def read_group_vcf_keys(path: str) -> list[tuple[str, int, str, str]]:
    keys = []
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            keys.append((fields[0], int(fields[1]), fields[3].upper(), fields[4].upper()))
    return keys


def read_provenance(path: str) -> list[tuple[str, str, str, int, int]]:
    rows = []
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        expected = ["group_barcode", "sample_id", "raw_barcode", "group_column", "source_column"]
        if header != expected:
            raise ContractError(f"{path}: unexpected header {header}")
        for line in handle:
            group_barcode, sample_id, raw, group_col, source_col = line.rstrip("\n").split("\t")
            rows.append((group_barcode, sample_id, raw, int(group_col), int(source_col)))
    return rows


def read_offsets(path: str) -> list[tuple[str, str, int, int, int]]:
    rows = []
    with open(path) as handle:
        handle.readline()
        for line in handle:
            sample_id, role, first, last, n_cells = line.rstrip("\n").split("\t")
            rows.append((sample_id, role, int(first), int(last), int(n_cells)))
    return rows


def run_checks(groups_csv: str, group_id: str, outdir: str, panel_vcf: str | None,
               compare_dir: str | None) -> Report:
    report = Report()
    group_dir = os.path.join(outdir, group_id)

    libraries = read_groups_csv(groups_csv, group_id)
    for library in libraries:
        read_vcf(library)
        read_barcodes(library)

    contigs = read_contig_order(panel_vcf)
    rank = {name: index for index, name in enumerate(contigs)}

    # --- axis identities -------------------------------------------------
    group_keys = read_group_vcf_keys(os.path.join(group_dir, "variants.vcf"))
    provenance = read_provenance(os.path.join(group_dir, "barcode_provenance.tsv"))
    with open(os.path.join(group_dir, "barcodes.tsv")) as handle:
        group_barcodes = [line.strip() for line in handle if line.strip()]
    offsets = read_offsets(os.path.join(group_dir, "library_offsets.tsv"))

    alt_in, alt_dims = open_mtx_body(os.path.join(group_dir, "alt.mtx"))
    ref_in, ref_dims = open_mtx_body(os.path.join(group_dir, "ref.mtx"))
    n_rows, n_cols, n_entries = alt_dims

    report.check("alt/ref declared dimensions identical", alt_dims == ref_dims,
                 f"{alt_dims} vs {ref_dims}")
    report.check("rows equal VCF records", n_rows == len(group_keys),
                 f"{n_rows} rows vs {len(group_keys)} records")
    report.check("columns equal barcodes", n_cols == len(group_barcodes),
                 f"{n_cols} columns vs {len(group_barcodes)} barcodes")
    report.check("columns equal provenance rows", n_cols == len(provenance),
                 f"{n_cols} columns vs {len(provenance)} provenance rows")
    report.check("barcodes unique", len(set(group_barcodes)) == len(group_barcodes),
                 "duplicate group barcodes")
    report.check("variant keys unique", len(set(group_keys)) == len(group_keys),
                 "duplicate variant keys")
    report.check(
        "variants in reference-panel order",
        all(
            (rank.get(a[0], -1), a[1], a[2], a[3]) < (rank.get(b[0], -1), b[1], b[2], b[3])
            for a, b in zip(group_keys, group_keys[1:])
        ),
        "variant order is not strictly increasing in panel order",
    )
    report.check(
        "barcodes match provenance column order",
        group_barcodes == [row[0] for row in provenance]
        and [row[3] for row in provenance] == list(range(1, len(provenance) + 1)),
        "barcodes.tsv and barcode_provenance.tsv disagree",
    )

    # --- union correctness ----------------------------------------------
    union_from_inputs = set()
    for library in libraries:
        union_from_inputs.update(library.keys)
    report.check("group variants equal union of library variants",
                 set(group_keys) == union_from_inputs,
                 f"{len(set(group_keys))} group keys vs {len(union_from_inputs)} input keys")

    # --- library coverage -------------------------------------------------
    report.check("every source library present",
                 [row[0] for row in offsets] == [lib.sample_id for lib in libraries],
                 "library_offsets.tsv does not list the samplesheet libraries in order")
    report.check("group columns equal sum of library columns",
                 n_cols == sum(lib.n_cells for lib in libraries),
                 f"{n_cols} vs {sum(lib.n_cells for lib in libraries)}")
    contiguous = True
    cursor = 1
    for (sample_id, _role, first, last, n_cells), library in zip(offsets, libraries):
        if (first, last, n_cells) != (cursor, cursor + library.n_cells - 1, library.n_cells):
            contiguous = False
        cursor += library.n_cells
    report.check("library column blocks contiguous and correctly sized", contiguous,
                 "library_offsets.tsv blocks do not tile the group columns")

    prov_ok = True
    index = 0
    for library in libraries:
        for source_column, barcode in enumerate(library.barcodes, 1):
            row = provenance[index]
            if row[1:3] != (library.sample_id, barcode) or row[4] != source_column:
                prov_ok = False
            index += 1
    report.check("provenance reproduces per-library barcode order", prov_ok,
                 "provenance rows disagree with the library barcodes.tsv files")

    # --- round trip through the group matrices ----------------------------
    union_index = {key: i + 1 for i, key in enumerate(group_keys)}
    column_owner: list[int] = []
    for library_index, library in enumerate(libraries):
        column_owner.extend([library_index] * library.n_cells)

    rowmaps = [[union_index.get(key, 0) for key in lib.keys] for lib in libraries]
    col_offsets = []
    cursor = 0
    for library in libraries:
        col_offsets.append(cursor)
        cursor += library.n_cells

    streams = []
    for library in libraries:
        a_in, a_dims = open_mtx_body(library.path("alt.mtx"))
        r_in, r_dims = open_mtx_body(library.path("ref.mtx"))
        streams.append([a_in, r_in, a_dims, r_dims, 0, 0, 0])  # +alt_sum, ref_sum, seen

    coordinate_identity = True
    roundtrip = True
    detail = ""
    group_alt_sum = 0
    group_ref_sum = 0
    seen_entries = 0

    try:
        for alt_line, ref_line in zip(alt_in, ref_in):
            a = alt_line.split()
            r = ref_line.split()
            seen_entries += 1
            if a[0] != r[0] or a[1] != r[1]:
                coordinate_identity = False
                detail = detail or f"group ALT/REF coordinates differ at entry {seen_entries}"
                break
            row = int(a[0])
            col = int(a[1])
            if not 1 <= row <= n_rows or not 1 <= col <= n_cols:
                roundtrip = False
                detail = detail or f"entry {seen_entries} outside declared dimensions"
                break
            group_alt_sum += int(float(a[2]))
            group_ref_sum += int(float(r[2]))

            library_index = column_owner[col - 1]
            state = streams[library_index]
            src_alt = state[0].readline().split()
            src_ref = state[1].readline().split()
            if not src_alt or not src_ref:
                roundtrip = False
                detail = detail or (
                    f"{libraries[library_index].sample_id}: group holds more entries than the input"
                )
                break
            expected_row = rowmaps[library_index][int(src_alt[0]) - 1]
            expected_col = int(src_alt[1]) + col_offsets[library_index]
            if (row, col) != (expected_row, expected_col) or a[2] != src_alt[2] or r[2] != src_ref[2]:
                roundtrip = False
                detail = detail or (
                    f"{libraries[library_index].sample_id}: entry mismatch, group "
                    f"({row},{col},{a[2]}/{r[2]}) vs input "
                    f"({expected_row},{expected_col},{src_alt[2]}/{src_ref[2]})"
                )
                break
            if src_alt[0] != src_ref[0] or src_alt[1] != src_ref[1]:
                coordinate_identity = False
                detail = detail or f"{libraries[library_index].sample_id}: input ALT/REF coordinates differ"
                break
            state[4] += int(float(src_alt[2]))
            state[5] += int(float(src_ref[2]))
            state[6] += 1

        leftovers = [
            libraries[i].sample_id
            for i, state in enumerate(streams)
            if state[0].readline()
        ]
        if leftovers and roundtrip:
            roundtrip = False
            detail = f"library entries missing from the group matrices: {leftovers}"
    finally:
        alt_in.close()
        ref_in.close()
        for state in streams:
            state[0].close()
            state[1].close()

    report.check("alt/ref coordinate pairs identical", coordinate_identity, detail)
    report.check("per-library entries reproduced exactly", roundtrip, detail)
    report.check("entry count matches header", seen_entries == n_entries,
                 f"header declares {n_entries}, read {seen_entries}")

    per_library_sums_ok = True
    for library, state in zip(libraries, streams):
        _a, _r, a_dims, r_dims, alt_sum, ref_sum, seen = state
        if a_dims != r_dims or a_dims[0] != library.n_variants or a_dims[1] != library.n_cells:
            per_library_sums_ok = False
        if seen != a_dims[2]:
            per_library_sums_ok = False
    report.check("library dimensions and entry counts reproduced", per_library_sums_ok,
                 "a library's dimensions or entry count was not reproduced")
    report.check("group allele sums equal library allele sums",
                 group_alt_sum == sum(state[4] for state in streams)
                 and group_ref_sum == sum(state[5] for state in streams),
                 "ALT/REF sums differ between group and libraries")

    # --- reported QC consistency ------------------------------------------
    qc_path = os.path.join(group_dir, "matrix_qc.tsv")
    qc_ok = False
    with open(qc_path) as handle:
        handle.readline()
        for line in handle:
            scope, sample, n_var, n_cell, n_ent, alt_sum, ref_sum, _status = line.rstrip("\n").split("\t")
            if scope == "group" and sample == group_id:
                qc_ok = (
                    int(n_var) == n_rows
                    and int(n_cell) == n_cols
                    and int(n_ent) == seen_entries
                    and int(alt_sum) == group_alt_sum
                    and int(ref_sum) == group_ref_sum
                )
    report.check("matrix_qc.tsv group row agrees with the matrices", qc_ok,
                 "matrix_qc.tsv group row does not match the observed matrices")

    # --- determinism -------------------------------------------------------
    if compare_dir:
        other = os.path.join(compare_dir, group_id)
        identical = []
        for name in DETERMINISM_FILES:
            a_path = os.path.join(group_dir, name)
            b_path = os.path.join(other, name)
            identical.append(
                os.path.isfile(b_path) and sha256(a_path) == sha256(b_path)
            )
        report.check("two builds are byte-identical", all(identical),
                     f"differing files: "
                     f"{[n for n, ok in zip(DETERMINISM_FILES, identical) if not ok]}")

    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--groups", required=True)
    parser.add_argument("--group-id", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--panel-vcf", default=None)
    parser.add_argument("--compare-dir", default=None,
                        help="a second independent build of the same group, for the determinism test")
    args = parser.parse_args(argv)

    try:
        report = run_checks(args.groups, args.group_id, args.outdir, args.panel_vcf,
                            args.compare_dir)
    except (ContractError, OSError) as exc:
        print(f"CHECK COULD NOT RUN: {exc}", file=sys.stderr)
        return 2

    print(f"\n{len(report.passes)} passed, {len(report.failures)} failed")
    for failure in report.failures:
        print(f"  FAILED: {failure}")
    return 1 if report.failures else 0


if __name__ == "__main__":
    sys.exit(main())
