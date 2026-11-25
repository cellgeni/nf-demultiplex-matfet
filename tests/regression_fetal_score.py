#!/usr/bin/env python3
"""Regression check: bin/fetal_score.py against the frozen prototype it was ported from.

The prototype (score_v3.py) produced the published embryo results
and the D8 benchmark. The port must reproduce its outputs for the same inputs:
  - cells:   the same cells in the same order; library, barcode, n_sites, n_umis, rho_hat,
             status and y_umis equal; every score column equal within --tol (the prototype ran
             in a conda environment with scipy 1.16, the container has scipy 1.15);
  - summary: the same rows; counts equal, the other columns within --tol;
  - removed sites: identical.
The prototype's atlas-only `supergroup` column is not compared.

Usage: regression_fetal_score.py <reference prefix> <port prefix> [--tol 1e-6]
Exits 0 when everything matches, 1 otherwise, printing the first differences.
"""
import argparse
import sys

import numpy as np
import pandas as pd

EXACT_CELLS = ["library", "barcode", "n_sites", "n_umis", "rho_hat", "status", "y_umis"]
EXACT_SUMMARY = ["scope", "cells", "deep_cells", "flagged", "insufficient"]
SKIP = {"supergroup"}


def compare(name, ref, new, exact, tol):
    problems = []
    if len(ref) != len(new):
        return [f"{name}: {len(ref)} rows in the reference, {len(new)} in the port"]
    cols = [c for c in ref.columns if c not in SKIP]
    missing = [c for c in cols if c not in new.columns]
    if missing:
        problems.append(f"{name}: columns missing from the port: {missing}")
    for col in cols:
        if col in missing:
            continue
        r, n = ref[col].reset_index(drop=True), new[col].reset_index(drop=True)
        if col in exact or not pd.api.types.is_float_dtype(r):
            bad = ~((r == n) | (r.isna() & n.isna()))
            if bad.any():
                i = int(np.flatnonzero(bad.values)[0])
                problems.append(f"{name}.{col}: {int(bad.sum()):,} rows differ, first at row {i}: {r[i]!r} vs {n[i]!r}")
        else:
            diff = (r - n).abs()
            scale = np.maximum(r.abs(), 1.0)
            worst = float((diff / scale).max(skipna=True)) if len(r) else 0.0
            nan_mismatch = (r.isna() != n.isna()).sum()
            if worst > tol or nan_mismatch:
                problems.append(f"{name}.{col}: largest relative difference {worst:.2e} (tolerance {tol:g}), "
                                f"{int(nan_mismatch)} NaN mismatches")
            else:
                print(f"  {name}.{col}: largest relative difference {worst:.1e}")
    return problems


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reference", help="prototype output prefix (…/<tag>, without .cells.tsv.gz)")
    ap.add_argument("port", help="bin/fetal_score.py output prefix")
    ap.add_argument("--tol", type=float, default=1e-6, help="relative tolerance for score columns")
    args = ap.parse_args(argv)

    problems = []
    rc = pd.read_csv(f"{args.reference}.cells.tsv.gz", sep="\t")
    nc = pd.read_csv(f"{args.port}.cells.tsv.gz", sep="\t")
    problems += compare("cells", rc, nc, EXACT_CELLS, args.tol)
    rs = pd.read_csv(f"{args.reference}.summary.tsv", sep="\t")
    ns = pd.read_csv(f"{args.port}.summary.tsv", sep="\t")
    problems += compare("summary", rs, ns, EXACT_SUMMARY, args.tol)
    with open(f"{args.reference}.removed_sites.tsv") as a, open(f"{args.port}.removed_sites.tsv") as b:
        if a.read() != b.read():
            problems.append("removed_sites: files differ")

    status = rc.status.value_counts().to_dict()
    print(f"{args.port}: {len(rc):,} cells; reference status counts {status}")
    if problems:
        print("FAIL")
        for p in problems[:20]:
            print("  " + p)
        return 1
    print("PASS: the port reproduces the prototype")
    return 0


if __name__ == "__main__":
    sys.exit(main())
