#!/usr/bin/env python3
"""Aggregate pregnancy-panel marker expression per group genotype cluster (singlets only).

The genotype clusters come from one Souporcell run over a whole pregnancy group, so
marker evidence is pooled across every library in the group: the fetal cells may sit
almost entirely in one library. Reads each library's CellRanger
filtered_feature_bc_matrix directly with scipy (no scanpy dependency), the group
clusters.tsv, and barcode_provenance.tsv to map each group barcode back to its library
and raw barcode. Resolves each panel gene (symbol -> Ensembl, with alias fallback; FAIL
LOUDLY if unresolved or if libraries resolve it differently), and emits a long-format
cluster_marker_summary.tsv consumed by label_clusters.py.

Per (cluster, marker), summed over the group's libraries:
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


def read_provenance(path):
    """group_barcode -> (sample_id, raw_barcode)."""
    mapping = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            mapping[row["group_barcode"]] = (row["sample_id"], row["raw_barcode"])
    if not mapping:
        sys.exit("ERROR: %s has no provenance rows" % path)
    return mapping


def read_matrices(path):
    """Two-column TSV: sample_id, expression matrix directory (in group order)."""
    matrices = []
    with open(path) as fh:
        for line in fh:
            if line.strip():
                sample_id, matrix_dir = line.rstrip("\n").split("\t")
                matrices.append((sample_id, matrix_dir))
    return matrices


def align(barcodes, wanted):
    """Map wanted raw barcodes -> expression column, retrying suffix-stripped."""
    index = {b: i for i, b in enumerate(barcodes)}
    hits = {b: index[b] for b in wanted if b in index}
    if hits:
        return hits, "raw"
    norm_index = {}
    for i, b in enumerate(_normalise_barcodes(barcodes)):
        norm_index.setdefault(b, i)
    hits = {}
    for b in wanted:
        j = norm_index.get(re.sub(r"-\d+$", "", b))
        if j is not None:
            hits[b] = j
    return hits, "suffix-stripped"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--group-id", required=True)
    ap.add_argument("--clusters", required=True, help="group Souporcell clusters.tsv")
    ap.add_argument("--provenance", required=True, help="barcode_provenance.tsv")
    ap.add_argument("--matrices", required=True,
                    help="TSV of sample_id and expression matrix directory")
    ap.add_argument("--panel", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--resolution-out", required=True)
    args = ap.parse_args(argv)

    panel = read_panel(args.panel)
    provenance = read_provenance(args.provenance)
    matrices = read_matrices(args.matrices)

    # Singlet raw barcode -> cluster, per library.
    singlets = {}
    with open(args.clusters) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["status"] != STATUS_SINGLET:
                continue
            if row["barcode"] not in provenance:
                sys.exit("ERROR: barcode %r is not in %s" % (row["barcode"], args.provenance))
            sample_id, raw = provenance[row["barcode"]]
            singlets.setdefault(sample_id, {})[raw] = row["assignment"]

    missing = sorted(set(sample for sample, _raw in provenance.values())
                     - set(sample for sample, _dir in matrices))
    if missing:
        sys.exit("ERROR: no expression matrix supplied for %s" % ", ".join(missing))

    resolution = None
    # cluster -> {"n": cells, "total": library-size sum, "pos": {m: n}, "sum": {m: counts}}
    acc = {}
    for sample_id, matrix_dir in matrices:
        mat, barcodes, features = load_matrix(matrix_dir)
        resolved = resolve_genes(panel, features)
        signature = [(r["symbol"], r["ensembl_id"]) for r in resolved]
        if resolution is None:
            resolution = resolved
        elif signature != [(r["symbol"], r["ensembl_id"]) for r in resolution]:
            sys.exit("ERROR: %s resolves the marker panel to different genes than the other "
                     "libraries in group %s; check the references." % (sample_id, args.group_id))

        wanted = singlets.get(sample_id, {})
        if not wanted:
            sys.stderr.write("[aggregate_markers] %s: no singlets\n" % sample_id)
            continue
        hits, how = align(barcodes, wanted)
        if not hits:
            sys.exit("ERROR: %s: no overlap between expression-matrix barcodes and Souporcell "
                     "singlet barcodes (checked raw and suffix-stripped)." % sample_id)
        sys.stderr.write("[aggregate_markers] %s: %d of %d singlets matched (%s)\n"
                         % (sample_id, len(hits), len(wanted), how))

        col_sums = np.asarray(mat.sum(axis=0)).ravel()
        marker_rows = {r["symbol"]: np.asarray(mat.getrow(r["row"]).todense()).ravel()
                       for r in resolved}
        cluster_cols = {}
        for raw, j in hits.items():
            cluster_cols.setdefault(wanted[raw], []).append(j)
        for cl, cols_list in cluster_cols.items():
            cols = np.array(cols_list, dtype=int)
            bucket = acc.setdefault(cl, {"n": 0, "total": 0.0, "pos": {}, "sum": {}})
            bucket["n"] += len(cols)
            bucket["total"] += float(col_sums[cols].sum())
            for r in resolved:
                vals = marker_rows[r["symbol"]][cols]
                bucket["pos"][r["symbol"]] = bucket["pos"].get(r["symbol"], 0) + int((vals > 0).sum())
                bucket["sum"][r["symbol"]] = bucket["sum"].get(r["symbol"], 0.0) + float(vals.sum())

    if not acc:
        sys.exit("ERROR: group %s has no singlets matched to expression" % args.group_id)

    with open(args.resolution_out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["symbol", "tier", "resolved_via", "matched_symbol", "ensembl_id"])
        for r in resolution:
            w.writerow([r["symbol"], r["tier"], r["via"], r["matched_symbol"], r["ensembl_id"]])

    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "group_id", "cluster", "n_cluster_cells", "marker", "tier",
            "ensembl_id", "n_pos", "frac_pos", "cp10k", "mean_counts",
        ])
        for cl in sorted(acc, key=lambda c: (len(c), c)):
            b = acc[cl]
            n = b["n"]
            for r in resolution:
                n_pos = b["pos"][r["symbol"]]
                gene_sum = b["sum"][r["symbol"]]
                w.writerow([
                    args.group_id, cl, n, r["symbol"], r["tier"], r["ensembl_id"], n_pos,
                    "%.6f" % (n_pos / n if n else 0.0),
                    "%.4f" % (gene_sum / b["total"] * 1e4 if b["total"] > 0 else 0.0),
                    "%.4f" % (gene_sum / n if n else 0.0),
                ])

    print("[aggregate_markers] %s -> %d clusters, %d markers, %d libraries"
          % (args.group_id, len(acc), len(resolution), len(matrices)))


if __name__ == "__main__":
    main()
