#!/usr/bin/env python3
"""Aggregate pregnancy-panel marker expression per Souporcell cluster (singlets only).

Reads a CellRanger filtered_feature_bc_matrix directory directly with scipy (no scanpy
dependency) plus a Souporcell clusters.tsv, resolves each panel gene (symbol -> Ensembl,
with alias fallback; FAIL LOUDLY if unresolved), and emits a long-format
cluster_marker_summary.tsv consumed by label_clusters.py.

Per (cluster, marker):
  n_cluster_cells, n_pos, frac_pos, cp10k (library-size-normalised pseudobulk, CP10K),
  mean_counts.  Mirrors the nf-demultiplex check_sex.R normalisation (CP10K + fraction
  of cells).
"""
import argparse
import csv
import gzip
import os
import re
import sys

import numpy as np
from scipy import io as sio
from scipy import sparse

STATUS_SINGLET = "singlet"


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def _find(matrix_dir, names):
    for n in names:
        for cand in (n, n + ".gz"):
            p = os.path.join(matrix_dir, cand)
            if os.path.exists(p):
                return p
    sys.exit("ERROR: none of %s found in %s" % (names, matrix_dir))


def load_matrix(matrix_dir):
    """Return (matrix CSR genes x cells, barcodes list, features list of (id, symbol))."""
    feat_p = _find(matrix_dir, ["features.tsv", "genes.tsv"])
    bc_p = _find(matrix_dir, ["barcodes.tsv"])
    mtx_p = _find(matrix_dir, ["matrix.mtx"])

    features = []
    with _open(feat_p) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            gid = row[0]
            sym = row[1] if len(row) > 1 else row[0]
            features.append((gid, sym))
    with _open(bc_p) as fh:
        barcodes = [line.strip() for line in fh if line.strip()]

    mat = sio.mmread(mtx_p)
    mat = sparse.csr_matrix(mat)
    if mat.shape[0] != len(features) or mat.shape[1] != len(barcodes):
        # Some matrices are cells x genes; transpose if that orientation matches.
        if mat.shape[1] == len(features) and mat.shape[0] == len(barcodes):
            mat = sparse.csr_matrix(mat.T)
        else:
            sys.exit(
                "ERROR: matrix shape %s does not match %d features x %d barcodes"
                % (mat.shape, len(features), len(barcodes))
            )
    return mat, barcodes, features


def read_panel(path):
    """Return list of dicts: {symbol, tier, aliases:[...]}."""
    panel = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if parts[0] == "symbol":
                continue
            symbol = parts[0].strip()
            tier = parts[1].strip() if len(parts) > 1 and parts[1].strip() else "secondary"
            aliases = []
            if len(parts) > 2 and parts[2].strip():
                aliases = [a.strip() for a in parts[2].split(",") if a.strip()]
            panel.append({"symbol": symbol, "tier": tier, "aliases": aliases})
    return panel


def resolve_genes(panel, features):
    """Resolve each panel gene to a matrix row index. FAIL LOUDLY if any unresolved."""
    sym_to_row = {}
    id_to_row = {}
    for i, (gid, sym) in enumerate(features):
        sym_to_row.setdefault(sym.upper(), i)
        id_to_row.setdefault(gid.upper(), i)

    resolved = []
    missing = []
    for g in panel:
        candidates = [g["symbol"]] + g["aliases"]
        hit = None
        via = None
        for name in candidates:
            key = name.upper()
            if key in sym_to_row:
                hit, via = sym_to_row[key], ("symbol" if name == g["symbol"] else "alias:%s" % name)
                break
            if key in id_to_row:
                hit, via = id_to_row[key], "ensembl"
                break
        if hit is None:
            missing.append(g["symbol"])
        else:
            gid, sym = features[hit]
            resolved.append({
                "symbol": g["symbol"], "tier": g["tier"], "row": hit,
                "matched_symbol": sym, "ensembl_id": gid, "via": via,
            })
    if missing:
        sys.exit(
            "ERROR: %d panel gene(s) could not be resolved in the expression matrix: %s\n"
            "Add an alias to the marker panel or check the reference."
            % (len(missing), ", ".join(missing))
        )
    return resolved


def _normalise_barcodes(bcs):
    return [re.sub(r"-\d+$", "", b) for b in bcs]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--clusters", required=True)
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--resolution-out", required=True)
    args = ap.parse_args(argv)

    mat, barcodes, features = load_matrix(args.matrix)
    panel = read_panel(args.panel)
    resolved = resolve_genes(panel, features)

    with open(args.resolution_out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["symbol", "tier", "resolved_via", "matched_symbol", "ensembl_id"])
        for r in resolved:
            w.writerow([r["symbol"], r["tier"], r["via"], r["matched_symbol"], r["ensembl_id"]])

    # Map singlet barcode -> cluster from Souporcell.
    bc_cluster = {}
    with open(args.clusters) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        idx = {n: i for i, n in enumerate(header)}
        bc_i, st_i, as_i = idx.get("barcode", 0), idx.get("status", 1), idx.get("assignment", 2)
        for row in reader:
            if row and row[st_i] == STATUS_SINGLET:
                bc_cluster[row[bc_i]] = row[as_i]

    # Align expression barcodes to Souporcell barcodes.
    bc_index = {b: i for i, b in enumerate(barcodes)}
    overlap = set(bc_index) & set(bc_cluster)
    used_map = bc_cluster
    if not overlap:
        # Retry after stripping the "-1"-style suffix on both sides.
        norm_expr = _normalise_barcodes(barcodes)
        norm_expr_index = {}
        for i, b in enumerate(norm_expr):
            norm_expr_index.setdefault(b, i)
        norm_clusters = {re.sub(r"-\d+$", "", b): c for b, c in bc_cluster.items()}
        overlap2 = set(norm_expr_index) & set(norm_clusters)
        if not overlap2:
            sys.exit(
                "ERROR: no overlap between expression-matrix barcodes and Souporcell "
                "singlet barcodes (checked raw and suffix-stripped). Check inputs."
            )
        bc_index = norm_expr_index
        used_map = norm_clusters
        sys.stderr.write(
            "[aggregate_markers] matched barcodes after suffix-stripping (%d overlap)\n"
            % len(overlap2)
        )
    else:
        sys.stderr.write("[aggregate_markers] %d barcodes overlap\n" % len(overlap))

    # Column indices per cluster.
    cluster_cols = {}
    for bc, cl in used_map.items():
        j = bc_index.get(bc)
        if j is not None:
            cluster_cols.setdefault(cl, []).append(j)

    # Precompute per-cell total counts (library size) and per-marker dense rows.
    col_sums = np.asarray(mat.sum(axis=0)).ravel()
    marker_rows = {r["symbol"]: np.asarray(mat.getrow(r["row"]).todense()).ravel() for r in resolved}

    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "sample_id", "cluster", "n_cluster_cells", "marker", "tier",
            "ensembl_id", "n_pos", "frac_pos", "cp10k", "mean_counts",
        ])
        for cl in sorted(cluster_cols):
            cols = np.array(cluster_cols[cl], dtype=int)
            n = len(cols)
            total = float(col_sums[cols].sum())
            for r in resolved:
                vals = marker_rows[r["symbol"]][cols]
                n_pos = int((vals > 0).sum())
                frac_pos = n_pos / n if n else 0.0
                gene_sum = float(vals.sum())
                cp10k = (gene_sum / total * 1e4) if total > 0 else 0.0
                mean_counts = gene_sum / n if n else 0.0
                w.writerow([
                    args.sample_id, cl, n, r["symbol"], r["tier"], r["ensembl_id"],
                    n_pos, "%.6f" % frac_pos, "%.4f" % cp10k, "%.4f" % mean_counts,
                ])

    print("[aggregate_markers] %s -> %d clusters, %d markers"
          % (args.sample_id, len(cluster_cols), len(resolved)))


if __name__ == "__main__":
    main()
