#!/usr/bin/env python3
"""Build barcode subset files for discovery evidence selection.

The full group matrices stay immutable: a discovery subset is a column view,
expressed as a barcode list, never a new matrix (group-matrix technical
contract, merge rule 5).

Modes:
  libraries   all cells from the named libraries (Gate 3's single-library
              anchor, and Gate 4's alternative discovery subsets)
  balanced    up to N cells per library across every library, ranked by
              informative loci, so rare evidence in a small library is not lost
              because its library was outside an arbitrary subset
  cap         a total cell cap, allocated evenly across libraries

Selection is deterministic: cells are ranked by informative loci, then by
barcode, so the same inputs always yield the same subset.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from assignment_table import read_locus_qc, write_locus_qc  # noqa: E402


def read_library_cells(group_dir: str) -> dict[str, list[str]]:
    path = os.path.join(group_dir, "barcode_provenance.tsv")
    cells: dict[str, list[str]] = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            cells.setdefault(row["sample_id"], []).append(row["group_barcode"])
    return cells


def locus_counts(group_dir: str, outdir: str) -> dict[str, int]:
    for candidate in (os.path.join(group_dir, "loci_per_cell.tsv"),
                      os.path.join(outdir, "loci_per_cell.tsv")):
        if os.path.isfile(candidate):
            return read_locus_qc(candidate)
    try:
        return write_locus_qc(group_dir, os.path.join(group_dir, "loci_per_cell.tsv"))
    except OSError:
        return write_locus_qc(group_dir, os.path.join(outdir, "loci_per_cell.tsv"))


def rank(cells: list[str], loci: dict[str, int]) -> list[str]:
    """Most informative first; barcode breaks ties so the order is stable."""
    return sorted(cells, key=lambda barcode: (-loci.get(barcode, 0), barcode))


def select(group_dir: str, outdir: str, mode: str, libraries: list[str] | None,
           per_library: int | None, total: int | None, min_loci: int) -> dict:
    os.makedirs(outdir, exist_ok=True)
    by_library = read_library_cells(group_dir)
    loci = locus_counts(group_dir, outdir)

    if mode == "libraries":
        if not libraries:
            raise SystemExit("--libraries is required for mode 'libraries'")
        unknown = [name for name in libraries if name not in by_library]
        if unknown:
            raise SystemExit(f"libraries absent from this group: {unknown}")
        chosen = {name: [c for c in by_library[name] if loci.get(c, 0) >= min_loci]
                  for name in libraries}
    elif mode == "balanced":
        if not per_library:
            raise SystemExit("--per-library is required for mode 'balanced'")
        chosen = {
            name: rank([c for c in cells if loci.get(c, 0) >= min_loci], loci)[:per_library]
            for name, cells in by_library.items()
        }
    elif mode == "cap":
        if not total:
            raise SystemExit("--total is required for mode 'cap'")
        share = max(1, total // max(1, len(by_library)))
        chosen = {
            name: rank([c for c in cells if loci.get(c, 0) >= min_loci], loci)[:share]
            for name, cells in by_library.items()
        }
    else:
        raise SystemExit(f"unknown mode {mode!r}")

    # Keep the group's column order so downstream tools see a consistent view.
    order = {}
    with open(os.path.join(group_dir, "barcodes.tsv")) as handle:
        for index, line in enumerate(handle):
            order[line.strip()] = index
    selected = sorted({c for cells in chosen.values() for c in cells},
                      key=lambda barcode: order[barcode])

    cells_path = os.path.join(outdir, "cells.tsv")
    with open(cells_path, "w") as out:
        for barcode in selected:
            out.write(barcode + "\n")

    report = {
        "group_dir": os.path.abspath(group_dir),
        "mode": mode,
        "libraries_requested": libraries,
        "per_library": per_library,
        "total": total,
        "min_informative_loci": min_loci,
        "n_selected": len(selected),
        "cells_file": os.path.abspath(cells_path),
        "per_library_selected": {name: len(cells) for name, cells in sorted(chosen.items())},
        "per_library_available": {name: len(cells) for name, cells in sorted(by_library.items())},
    }
    with open(os.path.join(outdir, "selection.json"), "w") as out:
        json.dump(report, out, indent=2, sort_keys=True)
        out.write("\n")
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--group-dir", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--mode", required=True, choices=["libraries", "balanced", "cap"])
    parser.add_argument("--libraries", default=None,
                        help="comma-separated sample_ids, for mode 'libraries'")
    parser.add_argument("--per-library", type=int, default=None)
    parser.add_argument("--total", type=int, default=None)
    parser.add_argument("--min-loci", type=int, default=1,
                        help="drop cells with fewer informative loci than this")
    args = parser.parse_args(argv)

    report = select(
        args.group_dir, args.outdir, args.mode,
        args.libraries.split(",") if args.libraries else None,
        args.per_library, args.total, args.min_loci,
    )
    print(f"selected {report['n_selected']} cells -> {report['cells_file']}")
    print(f"  per library: {report['per_library_selected']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
