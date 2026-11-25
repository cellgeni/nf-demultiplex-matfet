#!/usr/bin/env python3
"""Per-cell total, XIST and Y-gene UMIs from one library's filtered expression matrix.

Fetal mode uses them for a genotype-free check on male fetuses: a maternal cell expresses
XIST and no Y gene (fetal_summary.py). Same counting as the prototype check that bounded
F147 (xist_check.py).

Writes <prefix>.sex_counts.tsv.gz: library, barcode, total_umis, xist, y_umis.
"""
import argparse
import gzip
import os

import numpy as np
import pandas as pd

Y_GENES = ["RPS4Y1", "DDX3Y", "UTY", "KDM5D", "EIF1AY", "USP9Y", "ZFY"]


def first_existing(d, name):
    for f in (f"{d}/{name}.gz", f"{d}/{name}"):
        if os.path.exists(f):
            return f
    raise SystemExit(f"expression directory {d}: no {name}(.gz)")


def per_cell(lib, d):
    genes = pd.read_csv(first_existing(d, "features.tsv"), sep="\t", header=None)[1].tolist()
    idx = {}
    for i, g in enumerate(genes):
        idx.setdefault(g, i + 1)
    bcs = pd.read_csv(first_existing(d, "barcodes.tsv"), header=None)[0].values
    mf = first_existing(d, "matrix.mtx")
    with (gzip.open(mf, "rt") if mf.endswith(".gz") else open(mf)) as fh:
        line = fh.readline()
        while line.startswith("%"):
            line = fh.readline()
        _, n_cells, _ = map(int, line.split())
        m = pd.read_csv(fh, sep=" ", header=None, names=["g", "c", "v"],
                        dtype={"g": np.int32, "c": np.int32, "v": np.int32}, engine="c")
    total = np.bincount(m.c.values, weights=m.v.values, minlength=n_cells + 1)[1:]
    xm = m.g.values == idx.get("XIST", -1)
    xist = np.bincount(m.c.values[xm], weights=m.v.values[xm], minlength=n_cells + 1)[1:]
    ym = m.g.isin([idx[g] for g in Y_GENES if g in idx]).values
    y = np.bincount(m.c.values[ym], weights=m.v.values[ym], minlength=n_cells + 1)[1:]
    return pd.DataFrame({"library": lib, "barcode": bcs, "total_umis": total.astype(int),
                         "xist": xist.astype(int), "y_umis": y.astype(int)})


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--library", required=True)
    ap.add_argument("--expression", required=True, help="filtered matrix directory (CellRanger or STARsolo)")
    ap.add_argument("--prefix", required=True)
    args = ap.parse_args(argv)
    t = per_cell(args.library, args.expression)
    t.to_csv(f"{args.prefix}.sex_counts.tsv.gz", sep="\t", index=False)
    print(f"{args.library}: {len(t):,} cells; Y-positive {np.mean(t.y_umis > 0):.1%}; "
          f"XIST >= 2 and no Y: {int(((t.xist >= 2) & (t.y_umis == 0)).sum())}")


if __name__ == "__main__":
    main()
