#!/usr/bin/env python3
"""Regression check: bin/souporcell_distance.py against the paper's range check C'.

The reference is the loss-heatmap pairs.tsv (compute.py, job 287926),
which gave the published heatmaps. Running souporcell_distance.py on the same embryo Souporcell runs
must give, for every pair of clusters of those runs: the same matched_sites and msd_sites, and msd
equal within --tol. The paper found no mother in any embryo, so no cluster may be marked for review.

Usage: regression_souporcell_distance.py <reference pairs.tsv> <port outdir> [--tol 1e-9]
The port's library names must be the reference's run names. Exits 0 when everything matches.
"""
import argparse
import sys

import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reference")
    ap.add_argument("port")
    ap.add_argument("--tol", type=float, default=1e-9)
    args = ap.parse_args()

    new = pd.read_csv(f"{args.port}/pairs.tsv", sep="\t")
    ref = pd.read_csv(args.reference, sep="\t").rename(columns={"run_a": "library_a", "run_b": "library_b"})
    key = ["library_a", "cluster_a", "library_b", "cluster_b"]
    libs = set(new.library_a) | set(new.library_b)
    ref = ref[ref.library_a.isin(libs) & ref.library_b.isin(libs)]
    m = ref.merge(new, on=key, how="outer", suffixes=("_ref", "_new"), indicator=True)
    problems = []
    if (m._merge != "both").any():
        problems.append(f"{(m._merge != 'both').sum()} cluster pairs in only one of the two:\n"
                        f"{m[m._merge != 'both'][key + ['_merge']].head().to_string()}")
    b = m[m._merge == "both"]
    for col in ("matched_sites", "msd_sites"):
        bad = b[b[f"{col}_ref"] != b[f"{col}_new"]]
        if len(bad):
            problems.append(f"{col}: {len(bad)} pairs differ\n{bad[key + [col + '_ref', col + '_new']].head().to_string()}")
    diff = (b.msd_ref - b.msd_new).abs()
    if (diff > args.tol).any():
        problems.append(f"msd: {(diff > args.tol).sum()} pairs differ by more than {args.tol}; max {diff.max():.3g}")
    lib = pd.read_csv(f"{args.port}/libraries.tsv", sep="\t")
    print(lib.groupby("group").status.value_counts().to_string())
    if (lib.status == "review").any():
        problems.append(f"{(lib.status == 'review').sum()} clusters marked for review; the paper found none")
    print(f"{len(b)} cluster pairs compared over {len(libs)} runs; max |msd difference| {diff.max():.3g}")
    print("\n".join(problems) if problems else "OK: the port reproduces the reference")
    sys.exit(1 if problems else 0)


if __name__ == "__main__":
    main()
