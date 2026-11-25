#!/usr/bin/env python3
"""Fetal mode, per fetus: combine the batch scores into the result, its bound and its warnings.

Inputs are the fetus's batch outputs (fetal_score.py) and per-library sex counts
(sex_counts.py). The settings and their evidence are in docs/fetal_mode_design.md.

Per cell, fetal_mode is one of
  fetal             fetal_consistent: >= --min-sites informative sites, not flagged
  too_little_data   insufficient_evidence: fewer sites; not confirmed fetal, kept
  maternal          flagged as a maternal cell
  maternal_in_pair  flagged as a fetal-maternal pair (one droplet)
Flagged cells (maternal, maternal_in_pair) are the removal candidates.

Upper bound: the one-sided 95% bound on the share of maternal cells among the cells cleared as
fetal is 3.0 / cleared cells (the exact zero-event bound; no cleared cell is known to be
maternal), pooled over the fetus's batches. It holds only while the fetus is in fetal mode's
valid range: the pooled genotype is the fetus's only if maternal cells are a small share of each
batch (D8 benchmark: valid up to ~2%; above ~5-10% nothing is flagged and the bound is wrong).
So the result comes with in-range evidence:
  - sex: male if >= 50% of cells carry Y-gene UMIs, female if <= 5%, otherwise unclear (a
    possible mix of individuals: warned);
  - male fetus: cells with XIST >= 2 UMIs and no Y-gene UMI are maternal candidates found
    without genotypes; with P(XIST >= 2 | female cell) by depth from the female embryo F158
    (5' v2, prototype phase0/xist_check.py), (candidates + 3) / (P * cells) bounds the
    maternal share at any share. Warned if that bound is >= --max-maternal-share;
  - female fetus: no genotype-free check exists; stated in the outputs;
  - any batch whose flagged share of cells is >= --warn-flagged-fraction is warned: fetal mode
    is then near or past the end of its range, where it catches only part of the maternal cells.

Writes, in --outdir:
  fetal_summary.tsv          one row for the fetus and one per batch
  sex_check.tsv              sex call and the XIST/Y check
  warnings.txt               one warning per line (empty file when none)
  cells.tsv.gz               every cell: library, barcode, batch, fetal_mode and the evidence
  flagged_cells.tsv          the removal candidates
  per_library/<library>/fetal_mode_cells.tsv.gz   barcode, fetal_mode, n_sites, llr, llr_doublet
"""
import argparse
import os
import re

import numpy as np
import pandas as pd

FETAL_MODE = {"fetal_consistent": "fetal", "insufficient_evidence": "too_little_data",
              "maternal": "maternal", "fetal_maternal_doublet": "maternal_in_pair"}
FLAGGED = ("maternal", "maternal_in_pair")
MIN_XIST = 2
# P(XIST >= 2 UMIs | female cell) by total UMIs per cell, measured on F158 (female, 5' v2)
XIST_DETECT = [(0, 2000, 0.070), (2000, 5000, 0.181), (5000, 10000, 0.321), (10000, 20000, 0.577),
               (20000, 10**9, 0.871)]
MALE_MIN_Y_FRACTION = 0.5
FEMALE_MAX_Y_FRACTION = 0.05


def exact_upper95(cleared):
    """One-sided 95% upper bound on a share with zero events among `cleared` cells."""
    return 3.0 / cleared if cleared else float("nan")


def tally(c):
    n = len(c)
    counts = c.fetal_mode.value_counts()
    k = {m: int(counts.get(m, 0)) for m in FETAL_MODE.values()}
    flagged = k["maternal"] + k["maternal_in_pair"]
    return {"cells": n, "fetal": k["fetal"], "too_little_data": k["too_little_data"],
            "maternal": k["maternal"], "maternal_in_pair": k["maternal_in_pair"], "flagged": flagged,
            "flagged_fraction": flagged / n if n else float("nan"),
            "fetal_fraction_of_cells": k["fetal"] / n if n else float("nan"),
            "upper95_maternal_among_fetal": exact_upper95(k["fetal"])}


def xist_detection(total_umis):
    p = np.zeros(len(total_umis))
    for lo, hi, v in XIST_DETECT:
        p[(total_umis > lo) & (total_umis <= hi)] = v
    return float(p.mean()) if len(p) else float("nan")


def sex_check(sex, max_share):
    n = len(sex)
    y_frac = float((sex.y_umis > 0).mean()) if n else float("nan")
    call = "male" if y_frac >= MALE_MIN_Y_FRACTION else "female" if y_frac <= FEMALE_MAX_Y_FRACTION else "unclear"
    row = {"cells": n, "y_positive_fraction": y_frac, "sex_call": call,
           "xist_candidates": np.nan, "p_xist_detect": np.nan, "xist_upper95_maternal_share": np.nan}
    warnings = []
    if call == "male":
        k = int(((sex.xist >= MIN_XIST) & (sex.y_umis == 0)).sum())
        p = xist_detection(sex.total_umis.values)
        bound = (k + 3) / (p * n) if p > 0 else float("nan")
        row.update(xist_candidates=k, p_xist_detect=p, xist_upper95_maternal_share=bound)
        if not bound < max_share:
            warnings.append(f"XIST/Y check: maternal share up to {bound:.2%} ({k} cells with XIST >= {MIN_XIST} and no Y "
                            f"reads); fetal mode is valid only up to ~{max_share:.0%} maternal cells")
    elif call == "female":
        row["note"] = "female fetus: no genotype-free check of the maternal share exists"
    else:
        warnings.append(f"sex unclear: {y_frac:.1%} of cells carry Y-gene reads (male fetuses ~90%, female ~0%); "
                        "the group may mix individuals, which fetal mode assumes it does not")
    return row, warnings


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--group-id", required=True)
    ap.add_argument("--cells", nargs="+", required=True, help="<batch>.cells.tsv.gz from fetal_score.py")
    ap.add_argument("--sex-counts", nargs="+", required=True, help="<library>.sex_counts.tsv.gz from sex_counts.py")
    ap.add_argument("--warn-flagged-fraction", type=float, default=0.02)
    ap.add_argument("--max-maternal-share", type=float, default=0.02)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)

    parts = []
    for f in args.cells:
        batch = re.sub(r"\.cells\.tsv\.gz$", "", os.path.basename(f))
        parts.append(pd.read_csv(f, sep="\t").assign(batch=batch))
    cells = pd.concat(parts, ignore_index=True)
    unknown = set(cells.status) - set(FETAL_MODE)
    if unknown:
        raise SystemExit(f"unknown status values in the batch outputs: {sorted(unknown)}")
    cells["fetal_mode"] = cells.status.map(FETAL_MODE)

    sex = pd.concat([pd.read_csv(f, sep="\t") for f in args.sex_counts], ignore_index=True)
    sex["library"] = sex.library.astype(str)
    cells = cells.drop(columns=["y_umis"], errors="ignore").merge(sex, on=["library", "barcode"], how="left")

    rows = [{"scope": "fetus", "id": args.group_id, **tally(cells)}]
    warnings = []
    for batch, sub in cells.groupby("batch", sort=True):
        t = tally(sub)
        rows.append({"scope": "batch", "id": batch, **t})
        if t["flagged_fraction"] >= args.warn_flagged_fraction:
            warnings.append(f"batch {batch}: {t['flagged_fraction']:.2%} of cells flagged; at maternal shares near or above "
                            f"~2% fetal mode catches only part of them; check the batch with --mode pregnancy (Souporcell)")
    summary = pd.DataFrame(rows)
    summary.to_csv(f"{args.outdir}/fetal_summary.tsv", sep="\t", index=False, float_format="%.6g")

    srow, swarn = sex_check(sex, args.max_maternal_share)
    warnings = swarn + warnings
    pd.DataFrame([{"id": args.group_id, **srow}]).to_csv(f"{args.outdir}/sex_check.tsv", sep="\t", index=False,
                                                          float_format="%.6g")
    with open(f"{args.outdir}/warnings.txt", "w") as fh:
        fh.write("".join(w + "\n" for w in warnings))

    cols = ["library", "barcode", "batch", "fetal_mode", "status", "n_sites", "n_umis", "llr", "rho_hat",
            "llr_doublet", "total_umis", "xist", "y_umis"]
    cells[cols].to_csv(f"{args.outdir}/cells.tsv.gz", sep="\t", index=False, float_format="%.6g")
    cells[cells.fetal_mode.isin(FLAGGED)].sort_values("llr", ascending=False)[cols].to_csv(
        f"{args.outdir}/flagged_cells.tsv", sep="\t", index=False, float_format="%.6g")
    for lib, sub in cells.groupby("library", sort=False):
        d = f"{args.outdir}/per_library/{lib}"
        os.makedirs(d, exist_ok=True)
        sub[["barcode", "fetal_mode", "n_sites", "llr", "llr_doublet"]].to_csv(
            f"{d}/fetal_mode_cells.tsv.gz", sep="\t", index=False, float_format="%.6g")

    f = summary.iloc[0]
    print(f"{args.group_id}: {int(f.cells):,} cells; fetal {int(f.fetal):,}, too little data {int(f.too_little_data):,}, "
          f"maternal {int(f.maternal)}, in pair {int(f.maternal_in_pair)}; "
          f"upper 95% bound on maternal cells among fetal: {f.upper95_maternal_among_fetal:.2e}; "
          f"sex {srow['sex_call']}; {len(warnings)} warning(s)")
    for w in warnings:
        print("  WARNING " + w)


if __name__ == "__main__":
    main()
