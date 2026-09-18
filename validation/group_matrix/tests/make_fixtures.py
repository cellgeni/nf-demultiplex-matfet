#!/usr/bin/env python3
"""Write tiny synthetic Vartrix-shaped fixtures for the group-matrix harness.

Each fixture mimics a Souporcell output directory closely enough to exercise
the contract: variant union, missing loci, duplicate keys, conflicting alleles,
explicit zeros, barcode collisions across libraries, and malformed matrices.
"""

from __future__ import annotations

import os
import shutil
import sys

VCF_HEADER = [
    "##fileformat=VCFv4.2",
    "##source=make_fixtures.py",
] + [f"##contig=<ID=chr{c}>" for c in list(range(1, 23)) + ["X"]] + [
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
]


def write_library(root: str, sample_id: str, variants, barcodes, entries,
                  alt_dims=None, ref_coords=None):
    """entries: list of (row, col, alt_value, ref_value), 1-based."""
    path = os.path.join(root, sample_id)
    os.makedirs(path, exist_ok=True)

    with open(os.path.join(path, "common_variants_covered.vcf"), "w") as out:
        out.write("\n".join(VCF_HEADER) + "\n")
        for chrom, pos, ref, alt in variants:
            out.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\tAF=0.5\n")

    with open(os.path.join(path, "barcodes.tsv"), "w") as out:
        for barcode in barcodes:
            out.write(barcode + "\n")

    dims = alt_dims or (len(variants), len(barcodes), len(entries))
    for name, value_index in (("alt.mtx", 2), ("ref.mtx", 3)):
        with open(os.path.join(path, name), "w") as out:
            out.write("%%MatrixMarket matrix coordinate real general\n")
            out.write("% written by make_fixtures.py\n")
            out.write(f"{dims[0]} {dims[1]} {dims[2]}\n")
            for index, entry in enumerate(entries):
                row, col = entry[0], entry[1]
                if name == "ref.mtx" and ref_coords is not None:
                    row, col = ref_coords[index]
                out.write(f"{row} {col} {entry[value_index]}\n")
    return path


def build(root: str) -> str:
    if os.path.isdir(root):
        shutil.rmtree(root)
    os.makedirs(root)

    rows = []

    def add(group_id, sample_id, role="mixture", expected=2):
        rows.append(f"{group_id},{sample_id},{role},{expected},{os.path.join(root, sample_id)}")

    # --- GOOD: union, missing loci, explicit zeros, shared barcode strings ---
    # libA covers chr1:100, chr1:300, chr2:50 ; libB covers chr1:100, chr1:200, chr2:50
    write_library(
        root, "libA",
        variants=[("chr1", 100, "A", "G"), ("chr1", 300, "T", "C"), ("chr2", 50, "G", "A")],
        barcodes=["AAA-1", "CCC-1"],
        entries=[(1, 1, 3, 0), (1, 2, 0, 4), (2, 2, 1, 1), (3, 1, 0, 0)],
    )
    write_library(
        root, "libB",
        variants=[("chr1", 100, "A", "G"), ("chr1", 200, "C", "T"), ("chr2", 50, "G", "A")],
        # deliberately reuses libA barcode strings: the group must disambiguate them
        barcodes=["AAA-1", "GGG-1", "CCC-1"],
        entries=[(1, 3, 2, 2), (2, 1, 5, 0), (3, 2, 0, 7)],
    )
    add("GOOD", "libA")
    add("GOOD", "libB", role="maternal_reference")

    # --- DUP: duplicate variant key inside one library ---
    write_library(
        root, "libDup",
        variants=[("chr1", 100, "A", "G"), ("chr1", 100, "A", "G")],
        barcodes=["AAA-1"],
        entries=[(1, 1, 1, 1), (2, 1, 1, 1)],
    )
    add("DUP", "libDup")

    # --- CONFLICT: same position, different alleles across libraries ---
    write_library(
        root, "libC1",
        variants=[("chr1", 100, "A", "G")], barcodes=["AAA-1"], entries=[(1, 1, 1, 1)],
    )
    write_library(
        root, "libC2",
        variants=[("chr1", 100, "A", "T")], barcodes=["TTT-1"], entries=[(1, 1, 1, 1)],
    )
    add("CONFLICT", "libC1")
    add("CONFLICT", "libC2")

    # --- MULTI: multiallelic record ---
    write_library(
        root, "libMulti",
        variants=[("chr1", 100, "A", "G,T")], barcodes=["AAA-1"], entries=[(1, 1, 1, 1)],
    )
    add("MULTI", "libMulti")

    # --- DIMS: matrix rows disagree with the VCF ---
    write_library(
        root, "libDims",
        variants=[("chr1", 100, "A", "G")], barcodes=["AAA-1"], entries=[(1, 1, 1, 1)],
        alt_dims=(5, 1, 1),
    )
    add("DIMS", "libDims")

    # --- COORD: ALT and REF disagree on a coordinate ---
    write_library(
        root, "libCoord",
        variants=[("chr1", 100, "A", "G"), ("chr1", 200, "C", "T")],
        barcodes=["AAA-1"],
        entries=[(1, 1, 1, 1)],
        ref_coords=[(2, 1)],
    )
    add("COORD", "libCoord")

    # --- NEG: negative allele count ---
    write_library(
        root, "libNeg",
        variants=[("chr1", 100, "A", "G")], barcodes=["AAA-1"], entries=[(1, 1, -1, 1)],
    )
    add("NEG", "libNeg")

    # --- UNKNOWN_CONTIG: contig outside the panel order ---
    write_library(
        root, "libContig",
        variants=[("chrZ", 100, "A", "G")], barcodes=["AAA-1"], entries=[(1, 1, 1, 1)],
    )
    add("UNKNOWN_CONTIG", "libContig")

    groups_csv = os.path.join(root, "groups.csv")
    with open(groups_csv, "w") as out:
        out.write("group_id,sample_id,sample_role,expected_genotypes,souporcell_dir\n")
        out.write("\n".join(rows) + "\n")
    return groups_csv


if __name__ == "__main__":
    target = sys.argv[1] if len(sys.argv) > 1 else "fixtures"
    print(build(os.path.abspath(target)))
