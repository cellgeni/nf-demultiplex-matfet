#!/usr/bin/env python3
"""Materialise a column subset of a built group matrix.

A discovery subset is conceptually a column-index view (group-matrix contract,
merge rule 5), but both backends expect `barcodes.tsv` to have exactly one entry
per matrix column: handing them a shorter barcode list silently misaligns cells
or fails. So a subset is written out as a real, smaller build.

The variant axis is deliberately left untouched, including loci with no coverage
in the subset. That keeps row indices identical to the parent and to every other
subset, so a profile VCF learned on one subset can be applied to another without
any remapping.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import json
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from build_group_matrix import MM_HEADER, open_mtx_body  # noqa: E402


class SubsetError(Exception):
    pass


def link_or_copy(source: str, target: str) -> str:
    if os.path.exists(target):
        os.remove(target)
    try:
        os.link(source, target)
        return "hardlink"
    except OSError:
        shutil.copy2(source, target)
        return "copy"


def write_subset_matrix(group_dir: str, name: str, keep: list[bool],
                        new_index: list[int], out_dir: str) -> tuple[int, int]:
    """Filter one matrix to the kept columns. Returns (n_entries, allele_sum)."""
    source = os.path.join(group_dir, name)
    target = os.path.join(out_dir, name)
    body_path = target + ".body"

    handle, dims = open_mtx_body(source)
    n_rows, n_cols, _declared = dims
    if len(keep) != n_cols:
        raise SubsetError(
            f"{source}: {n_cols} columns but the parent barcode list has {len(keep)}"
        )

    entries = 0
    allele_sum = 0
    try:
        with open(body_path, "w") as body:
            for line in handle:
                fields = line.split()
                column = int(fields[1])
                if not keep[column - 1]:
                    continue
                body.write(f"{fields[0]} {new_index[column - 1]} {fields[2]}\n")
                entries += 1
                allele_sum += int(float(fields[2]))
    finally:
        handle.close()

    n_kept = sum(keep)
    with open(target, "w") as out:
        out.write(MM_HEADER + "\n")
        out.write(f"% column subset of {os.path.abspath(source)}\n")
        out.write(f"{n_rows} {n_kept} {entries}\n")
        with open(body_path) as body:
            shutil.copyfileobj(body, out)
    os.remove(body_path)
    return entries, allele_sum


def subset(group_dir: str, cells_file: str, outdir: str, group_id: str | None = None) -> dict:
    group_id = group_id or os.path.basename(os.path.normpath(group_dir))
    out_dir = os.path.join(outdir, group_id)
    os.makedirs(out_dir, exist_ok=True)

    with open(os.path.join(group_dir, "barcodes.tsv")) as handle:
        parent_barcodes = [line.strip() for line in handle if line.strip()]
    with open(cells_file) as handle:
        wanted = {line.strip() for line in handle if line.strip()}

    unknown = wanted - set(parent_barcodes)
    if unknown:
        raise SubsetError(
            f"{cells_file}: {len(unknown)} barcodes are not in the parent build, "
            f"e.g. {sorted(unknown)[:3]}"
        )
    if not wanted:
        raise SubsetError(f"{cells_file}: no barcodes selected")

    keep = [barcode in wanted for barcode in parent_barcodes]
    new_index: list[int] = []
    cursor = 0
    for flag in keep:
        if flag:
            cursor += 1
            new_index.append(cursor)
        else:
            new_index.append(0)

    kept_barcodes = [barcode for barcode, flag in zip(parent_barcodes, keep) if flag]
    with open(os.path.join(out_dir, "barcodes.tsv"), "w") as out:
        for barcode in kept_barcodes:
            out.write(barcode + "\n")

    # Provenance is re-derived from the parent so sample_id never comes from
    # parsing the synthetic barcode.
    parent_provenance: dict[str, dict] = {}
    with open(os.path.join(group_dir, "barcode_provenance.tsv"), newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            parent_provenance[row["group_barcode"]] = row

    per_library: dict[str, int] = {}
    with open(os.path.join(out_dir, "barcode_provenance.tsv"), "w") as out:
        out.write("group_barcode\tsample_id\traw_barcode\tgroup_column\tsource_column\n")
        for column, barcode in enumerate(kept_barcodes, 1):
            row = parent_provenance[barcode]
            per_library[row["sample_id"]] = per_library.get(row["sample_id"], 0) + 1
            out.write(f"{barcode}\t{row['sample_id']}\t{row['raw_barcode']}\t"
                      f"{column}\t{row['source_column']}\n")

    with open(os.path.join(out_dir, "library_offsets.tsv"), "w") as out:
        out.write("sample_id\tfirst_group_column\tlast_group_column\tn_cells\n")
        cursor = 0
        for sample_id, count in per_library.items():
            out.write(f"{sample_id}\t{cursor + 1}\t{cursor + count}\t{count}\n")
            cursor += count

    alt_entries, alt_sum = write_subset_matrix(group_dir, "alt.mtx", keep, new_index, out_dir)
    ref_entries, ref_sum = write_subset_matrix(group_dir, "ref.mtx", keep, new_index, out_dir)
    if alt_entries != ref_entries:
        raise SubsetError(
            f"{out_dir}: ALT kept {alt_entries} entries but REF kept {ref_entries}"
        )

    # The variant axis is shared with the parent, so the VCF is reused as is.
    variant_files = {}
    for name in ("variants.vcf", "variants.vcf.gz", "variants.vcf.gz.tbi"):
        source = os.path.join(group_dir, name)
        if os.path.isfile(source):
            variant_files[name] = link_or_copy(source, os.path.join(out_dir, name))

    # Informative-locus counts are inherited where the parent has them.
    parent_loci = os.path.join(group_dir, "loci_per_cell.tsv")
    if os.path.isfile(parent_loci):
        with open(parent_loci) as handle, open(os.path.join(out_dir, "loci_per_cell.tsv"), "w") as out:
            out.write(handle.readline())
            for line in handle:
                if line.split("\t", 1)[0] in wanted:
                    out.write(line)

    with open(os.path.join(out_dir, "matrix_qc.tsv"), "w") as out:
        out.write("scope\tsample_id\tn_variants\tn_cells\tn_entries\talt_sum\tref_sum\tstatus\n")
        n_rows = open_mtx_body(os.path.join(out_dir, "alt.mtx"))[1][0]
        out.write(f"subset\t{group_id}\t{n_rows}\t{len(kept_barcodes)}\t{alt_entries}\t"
                  f"{alt_sum}\t{ref_sum}\tpass\n")

    manifest = {
        "tool": "subset_group_matrix.py",
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "parent_group_dir": os.path.abspath(group_dir),
        "cells_file": os.path.abspath(cells_file),
        "group_id": group_id,
        "n_cells": len(kept_barcodes),
        "n_entries": alt_entries,
        "alt_sum": alt_sum,
        "ref_sum": ref_sum,
        "per_library_cells": per_library,
        "variant_axis": "identical to parent; row indices and variants.vcf reused unchanged",
        "variant_files": variant_files,
    }
    with open(os.path.join(out_dir, "manifest.json"), "w") as out:
        json.dump(manifest, out, indent=2, sort_keys=True)
        out.write("\n")
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--group-dir", required=True, help="parent built group directory")
    parser.add_argument("--cells", required=True, help="barcode subset file (cells.tsv)")
    parser.add_argument("--outdir", required=True, help="<outdir>/<group_id> is written")
    parser.add_argument("--group-id", default=None)
    args = parser.parse_args(argv)

    try:
        manifest = subset(args.group_dir, args.cells, args.outdir, args.group_id)
    except (SubsetError, OSError) as exc:
        print(f"SUBSET FAILURE: {exc}", file=sys.stderr)
        return 2
    print(f"subset written: {manifest['n_cells']} cells, {manifest['n_entries']} entries")
    print(f"  per library: {manifest['per_library_cells']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
