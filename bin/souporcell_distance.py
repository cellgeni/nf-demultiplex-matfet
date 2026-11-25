#!/usr/bin/env python3
"""Fetal mode, Souporcell cross-library check: genotype distance between every pair of Souporcell
clusters, per fetus and across fetuses, with SVG heatmaps for the report.

Each selected library was run through Souporcell at k = 2. With no mother, k = 2 still splits the
library in two, but both halves are the fetus. A mother cluster is the odd one out: it sits far from
the clusters of the fetus's other libraries. The distance is the per-site genotype distance used for
the paper's range check C' (the loss-heatmap compute.py, ported unchanged):

  msd = sum over sites of (ref-allele fraction in cluster A - in cluster B)^2 / number of sites
  counts from alt.mtx/ref.mtx, singlet cells only; sites matched by chr, pos, alt and kept where every
  cluster of both runs has >= 8 reads (shared_samples.py's rule).

Reference levels (embryos F137/F147/F158/F218 and the D8-D10 benchmark, docs/fetal_mode_design.md):
  same person, benchmark 0.0000-0.0065; any two clusters within one embryo 0.0048-0.0169;
  mother vs her own fetus 0.041-0.047; different pregnancies 0.068-0.089.
A cluster whose median distance to the fetus's other libraries reaches --warn-msd (default 0.03,
between 0.0169 and 0.041) is reported for review. Nothing is decided automatically: a warning means
"look at the heatmap", not "mother found".

Input: a manifest TSV with columns group, library, dir (a Souporcell output directory with
alt.mtx, ref.mtx, clusters.tsv, common_variants_covered.vcf).
Writes to --outdir: clusters.tsv, pairs.tsv, libraries.tsv, groups.tsv, warnings.tsv,
<group>.heatmap.svg per fetus, and groups.heatmap.svg when there are two fetuses or more.
"""
import argparse
import csv
import html
import os
import re

import numpy as np
import pandas as pd

MIN_READS = 8          # shared_samples.py / compute.py
VMAX = 0.09            # colour scale top, as the paper heatmaps (different pregnancies reach 0.089)
REF_LEVELS = [(0.0169, "same fetus", "max 0.017"), (0.041, "mother-fetus", "min 0.041"), (0.068, "unrelated", "min 0.068")]
BLUES = ["#f7fbff", "#deebf7", "#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5", "#08519c", "#08306b"]


def load_run(d):
    """Per-cluster ref and alt read sums at every site of one Souporcell run (singlets only)."""
    keys = []
    with open(os.path.join(d, "common_variants_covered.vcf")) as fh:
        for line in fh:
            if line[0] != "#":
                f = line.split("\t", 5)
                keys.append((f[0], f[1], f[4]))
    cl = pd.read_csv(os.path.join(d, "clusters.tsv"), sep="\t", dtype=str)
    k = sum(1 for c in cl.columns if re.fullmatch(r"cluster\d+", c)) or 2
    assign = np.full(len(cl) + 1, -1)                 # mtx cell index is the 1-based row of clusters.tsv
    sing = (cl["status"] == "singlet").to_numpy()
    assign[1:][sing] = cl["assignment"][sing].astype(int)
    alt = pd.read_csv(os.path.join(d, "alt.mtx"), sep=" ", skiprows=3, header=None).to_numpy()
    ref = pd.read_csv(os.path.join(d, "ref.mtx"), sep=" ", skiprows=3, header=None).to_numpy()
    if not (alt[:, :2] == ref[:, :2]).all():
        raise SystemExit(f"{d}: alt.mtx and ref.mtx list different entries")
    a = assign[alt[:, 1].astype(int)]
    loc = alt[:, 0].astype(int)
    R = np.zeros((k, len(keys) + 1))
    A = np.zeros((k, len(keys) + 1))
    for c in range(k):
        m = a == c
        R[c] = np.bincount(loc[m], weights=ref[m, 2], minlength=len(keys) + 1)
        A[c] = np.bincount(loc[m], weights=alt[m, 2], minlength=len(keys) + 1)
    sizes = [int((assign == c).sum()) for c in range(k)]
    return keys, R[:, 1:], A[:, 1:], sizes, len(keys)


def msd_pair(x, y):
    """Loss and site count for every cluster pair of runs x, y (compute.py's msd_pair)."""
    _, ix, iy = np.intersect1d(x["gids"], y["gids"], return_indices=True)
    Rx, Ax, Ry, Ay = x["R"][:, ix], x["A"][:, ix], y["R"][:, iy], y["A"][:, iy]
    keep = ((Rx + Ax) >= MIN_READS).all(0) & ((Ry + Ay) >= MIN_READS).all(0)
    fx = Rx[:, keep] / (Rx + Ax)[:, keep]
    fy = Ry[:, keep] / (Ry + Ay)[:, keep]
    loss = {(a, b): float(((fx[a] - fy[b]) ** 2).sum()) for a in range(len(fx)) for b in range(len(fy))}
    return loss, int(keep.sum()), len(ix)


def compute(manifest):
    gid_of, runs, clusters = {}, [], []
    for row in manifest:
        keys, R, A, sizes, n_sites = load_run(row["dir"])
        gids = np.array([gid_of.setdefault(k, len(gid_of)) for k in keys])
        runs.append({**row, "gids": gids, "R": R, "A": A})
        for c, n in enumerate(sizes):
            clusters.append({"group": row["group"], "library": row["library"], "cluster": c,
                             "label": f"{row['library']} c{c}", "singlets": n, "sites_vcf": n_sites})
        print(f"loaded {row['library']} ({row['group']}): {n_sites} sites, singlets {sizes}", flush=True)
    pairs = []
    for i, x in enumerate(runs):
        for y in runs[i:]:
            loss, n_used, n_match = msd_pair(x, y)
            for (a, b), v in loss.items():
                if x is y and b <= a:
                    continue
                pairs.append({"group_a": x["group"], "library_a": x["library"], "cluster_a": a,
                              "group_b": y["group"], "library_b": y["library"], "cluster_b": b,
                              "matched_sites": n_match, "msd_sites": n_used, "msd_loss": v,
                              "msd": v / n_used if n_used else np.nan})
    return pd.DataFrame(clusters), pd.DataFrame(pairs)


def per_library(cl, p, min_singlets, warn_msd):
    """One row per cluster: distance to its sister cluster and median distance to the clusters of the fetus's
    other libraries; plus the warnings for review."""
    rows, warns = [], []
    usable = {(r.library, r.cluster) for r in cl.itertuples() if r.singlets >= min_singlets}
    # both directions of every within-fetus pair: (library, cluster) -> (other library, other cluster)
    within = p[p.group_a == p.group_b]
    swap = {"library_a": "library_b", "cluster_a": "cluster_b", "library_b": "library_a", "cluster_b": "cluster_a"}
    d = pd.concat([within, within.rename(columns=swap)], ignore_index=True)
    d = d[[(lib, c) in usable for lib, c in zip(d.library_b, d.cluster_b)]]
    for g, gc in cl.groupby("group", sort=False):
        if gc.library.nunique() < 2:
            warns.append((g, "", "Only one library was run through Souporcell, so the cross-library check is not "
                          "possible. Run with --fetal_souporcell all, or give souporcell_dir for more libraries."))
        for r in gc.itertuples():
            mine = d[(d.library_a == r.library) & (d.cluster_a == r.cluster)]
            sister = mine[mine.library_b == r.library].msd
            others = mine[mine.library_b != r.library].msd.dropna()
            status = ("too_few_cells" if (r.library, r.cluster) not in usable
                      else "no_other_library" if others.empty else "ok")
            med = float(others.median()) if len(others) else np.nan
            if status == "ok" and med >= warn_msd:
                status = "review"
                warns.append((g, r.library, f"{r.library} cluster {r.cluster} ({r.singlets} singlets) has median "
                              f"distance {med:.4f} to the clusters of {g}'s other libraries (warn at {warn_msd}; "
                              f"same fetus <= 0.017, mother-fetus 0.041-0.047). Review the heatmap: this is how "
                              f"a maternal cluster would look, but it is not a decision."))
            rows.append({"group": g, "library": r.library, "cluster": r.cluster, "singlets": r.singlets,
                         "msd_to_sister": float(sister.iloc[0]) if len(sister) else np.nan,
                         "msd_to_other_libraries_median": med,
                         "msd_to_other_libraries_max": float(others.max()) if len(others) else np.nan,
                         "status": status})
    return pd.DataFrame(rows), pd.DataFrame(warns, columns=["group", "library", "warning"])


def group_matrix(cl, p, min_singlets):
    """Fetus x fetus median msd; diagonal over pairs of different clusters within the fetus (plot_condensed.py)."""
    usable = {(r.library, r.cluster) for r in cl.itertuples() if r.singlets >= min_singlets}
    q = p[[(a, b) in usable and (c, d) in usable
           for a, b, c, d in zip(p.library_a, p.cluster_a, p.library_b, p.cluster_b)]]
    groups = list(dict.fromkeys(cl.group))
    mat = pd.DataFrame(np.nan, index=groups, columns=groups)
    for a in groups:
        for b in groups:
            m = q[((q.group_a == a) & (q.group_b == b)) | ((q.group_a == b) & (q.group_b == a))].msd.dropna()
            if len(m):
                mat.loc[a, b] = m.median()
    return mat


# --- SVG ---------------------------------------------------------------------------------------
def colour(v):
    if v != v:
        return "#e6e6e6"
    t = min(max(v / VMAX, 0.0), 1.0) * (len(BLUES) - 1)
    i = min(int(t), len(BLUES) - 2)
    f = t - i
    c0, c1 = BLUES[i], BLUES[i + 1]
    rgb = [round(int(c0[j:j + 2], 16) * (1 - f) + int(c1[j:j + 2], 16) * f) for j in (1, 3, 5)]
    return "#%02x%02x%02x" % tuple(rgb)


def heatmap_svg(labels, mat, title, blocks=None, lower_only=False, bold=()):
    """labels: n row/column names; mat: n x n array (nan = no value); blocks: row indices where a new block starts."""
    n = len(labels)
    cell = 34 if n <= 8 else (22 if n <= 20 else 14)
    show_vals = n <= 8                                 # beyond that, values are in the tooltips and tables
    pad_l = 8 + 7 * max(len(s) for s in labels)
    pad_t = 26
    by = pad_t + n * cell + pad_l + 22                 # colour bar top, below the rotated column labels
    bw = max(n * cell, 300)
    w = max(pad_l + bw + 40, 16 + 7.2 * len(title))
    h = by + 44
    e = html.escape
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w:.0f}" height="{h:.0f}" font-family="sans-serif" '
           f'font-size="11" role="img" aria-label="{e(title)}">',
           f'<text x="8" y="16" font-size="12" font-weight="bold">{e(title)}</text>']
    for i in range(n):
        for j in range(n):
            if lower_only and j > i:
                continue
            v = mat[i][j]
            x, y = pad_l + j * cell, pad_t + i * cell
            out.append(f'<rect x="{x}" y="{y}" width="{cell}" height="{cell}" fill="{colour(v)}">'
                       f'<title>{e(labels[i])} vs {e(labels[j])}: {"n/a" if v != v else f"{v:.4f}"}</title></rect>')
            if show_vals and v == v:
                out.append(f'<text x="{x + cell / 2}" y="{y + cell / 2 + 3}" text-anchor="middle" font-size="8" '
                           f'fill="{"#ffffff" if v > 0.5 * VMAX else "#1b2126"}">{v:.3f}</text>')
    for b in blocks or []:
        z = b * cell
        out.append(f'<line x1="{pad_l}" y1="{pad_t + z}" x2="{pad_l + n * cell}" y2="{pad_t + z}" stroke="#fff" stroke-width="2"/>')
        out.append(f'<line x1="{pad_l + z}" y1="{pad_t}" x2="{pad_l + z}" y2="{pad_t + n * cell}" stroke="#fff" stroke-width="2"/>')
    for i, s in enumerate(labels):
        fw = ' font-weight="bold"' if s in bold else ""
        out.append(f'<text x="{pad_l - 4}" y="{pad_t + i * cell + cell / 2 + 4}" text-anchor="end"{fw}>{e(s)}</text>')
        cx, cy = pad_l + i * cell + cell / 2 + 4, pad_t + n * cell + 4
        out.append(f'<text x="{cx}" y="{cy}" text-anchor="end" transform="rotate(-90 {cx} {cy})"{fw}>{e(s)}</text>')
    # colour bar with the reference levels
    steps = 60
    for k in range(steps):
        out.append(f'<rect x="{pad_l + k * bw / steps:.1f}" y="{by}" width="{bw / steps + 0.5:.1f}" height="10" '
                   f'fill="{colour(VMAX * (k + 0.5) / steps)}"/>')
    for v, lab, val in REF_LEVELS:
        x = pad_l + bw * v / VMAX
        out.append(f'<line x1="{x:.1f}" y1="{by - 3}" x2="{x:.1f}" y2="{by + 13}" stroke="#1b2126"/>')
        out.append(f'<text x="{x:.1f}" y="{by + 24}" text-anchor="middle" font-size="9">{e(lab)}</text>')
        out.append(f'<text x="{x:.1f}" y="{by + 35}" text-anchor="middle" font-size="9">{val}</text>')
    out.append(f'<text x="{pad_l}" y="{by - 6}" font-size="9">genotype distance per SNP (0 to {VMAX})</text>')
    out.append("</svg>")
    return "\n".join(out)


def cluster_heatmap(g, cl, p):
    gc = cl[cl.group == g]
    keys = list(zip(gc.library, gc.cluster))
    labels = [f"{lib} c{c} ({n})" for lib, c, n in zip(gc.library, gc.cluster, gc.singlets)]
    idx = {k: i for i, k in enumerate(keys)}
    mat = [[np.nan] * len(keys) for _ in keys]
    pg = p[(p.group_a == g) & (p.group_b == g)]
    for r in pg.itertuples():
        i, j = idx[(r.library_a, r.cluster_a)], idx[(r.library_b, r.cluster_b)]
        mat[i][j] = mat[j][i] = r.msd
    blocks = [i for i in range(1, len(keys)) if keys[i][0] != keys[i - 1][0]]
    return heatmap_svg(labels, mat, f"{g}: Souporcell clusters of each checked library (singlets)", blocks)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True, help="TSV: group, library, dir")
    ap.add_argument("--min-singlets", type=int, default=20, help="clusters with fewer singlets are shown, not judged")
    ap.add_argument("--warn-msd", type=float, default=0.03)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args(argv)

    with open(args.manifest) as fh:
        manifest = list(csv.DictReader(fh, delimiter="\t"))
    os.makedirs(args.outdir, exist_ok=True)
    cl, p = compute(manifest)
    lib, warns = per_library(cl, p, args.min_singlets, args.warn_msd)
    mat = group_matrix(cl, p, args.min_singlets)

    o = lambda name: os.path.join(args.outdir, name)  # noqa: E731
    cl.to_csv(o("clusters.tsv"), sep="\t", index=False)
    p.to_csv(o("pairs.tsv"), sep="\t", index=False, float_format="%.6g")
    lib.to_csv(o("libraries.tsv"), sep="\t", index=False, float_format="%.5f")
    mat.to_csv(o("groups.tsv"), sep="\t", float_format="%.5f")
    warns.to_csv(o("warnings.tsv"), sep="\t", index=False)
    for g in mat.index:
        with open(o(f"{g}.heatmap.svg"), "w") as fh:
            fh.write(cluster_heatmap(g, cl, p))
    if len(mat) > 1:
        with open(o("groups.heatmap.svg"), "w") as fh:
            fh.write(heatmap_svg(list(mat.index), mat.to_numpy().tolist(),
                                 "Fetuses: median distance between clusters (diagonal: within the fetus)",
                                 lower_only=True))
    print(lib.to_string(index=False))
    print(mat.round(4).to_string())
    print(f"{len(warns)} warning(s)")


if __name__ == "__main__":
    main()
