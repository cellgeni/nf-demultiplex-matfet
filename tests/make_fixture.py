#!/usr/bin/env python3
"""Generate a tiny synthetic 10x fixture + Souporcell outputs for a male-fetus sample.

Creates under test_data/<sample>/:
  filtered_feature_bc_matrix/{features,barcodes}.tsv.gz + matrix.mtx.gz
  souporcell/{clusters.tsv, ambient_rna.txt, qc.txt}

Cluster 0 = maternal, cluster 1 = fetal (RPS4Y1 + secondary markers enriched).
Deterministic (fixed RNG seed) so integration output is reproducible.
"""
import argparse
import gzip
import os

import numpy as np
from scipy import io as sio
from scipy import sparse

PANEL = ["RPS4Y1", "HAPLN1", "PITX2", "GAPLINC", "HPGDS", "XAGE2", "XAGE3", "HLA-G"]
BACKGROUND = ["ACTB", "GAPDH", "MALAT1", "XIST", "RPL13", "TMSB4X"]


def build(sample_dir, n_per_cluster=120, seed=0):
    rng = np.random.default_rng(seed)
    genes = BACKGROUND + PANEL
    # Fake but plausible Ensembl-style IDs.
    ids = ["ENSGFAKE%08d" % i for i in range(len(genes))]
    gi = {g: i for i, g in enumerate(genes)}

    n0 = n_per_cluster  # maternal
    n1 = n_per_cluster  # fetal
    ncells = n0 + n1
    ngenes = len(genes)
    M = np.zeros((ngenes, ncells), dtype=int)

    # Background housekeeping in all cells.
    for g in BACKGROUND:
        M[gi[g], :] = rng.poisson(20, ncells)

    fetal_cols = np.arange(n0, n0 + n1)
    mat_cols = np.arange(0, n0)

    # Fetal cluster: RPS4Y1 strong + secondary markers enriched.
    M[gi["RPS4Y1"], fetal_cols] = rng.poisson(6, n1) * (rng.random(n1) < 0.55)
    for g in ["HAPLN1", "PITX2", "GAPLINC", "HPGDS", "XAGE2", "XAGE3", "HLA-G"]:
        M[gi[g], fetal_cols] = rng.poisson(4, n1) * (rng.random(n1) < 0.40)

    # Maternal cluster: only faint ambient leakage of fetal markers.
    M[gi["RPS4Y1"], mat_cols] = rng.poisson(1, n0) * (rng.random(n0) < 0.02)
    for g in ["HAPLN1", "PITX2", "GAPLINC", "HPGDS", "XAGE2", "XAGE3", "HLA-G"]:
        M[gi[g], mat_cols] = rng.poisson(1, n0) * (rng.random(n0) < 0.02)

    barcodes = ["CELL%05d-1" % i for i in range(ncells)]

    mdir = os.path.join(sample_dir, "filtered_feature_bc_matrix")
    os.makedirs(mdir, exist_ok=True)
    with gzip.open(os.path.join(mdir, "features.tsv.gz"), "wt") as fh:
        for i, g in enumerate(genes):
            fh.write("%s\t%s\tGene Expression\n" % (ids[i], g))
    with gzip.open(os.path.join(mdir, "barcodes.tsv.gz"), "wt") as fh:
        fh.write("\n".join(barcodes) + "\n")
    with gzip.open(os.path.join(mdir, "matrix.mtx.gz"), "wb") as fh:
        sio.mmwrite(fh, sparse.csr_matrix(M))

    # Souporcell outputs.
    sdir = os.path.join(sample_dir, "souporcell")
    os.makedirs(sdir, exist_ok=True)
    with open(os.path.join(sdir, "clusters.tsv"), "w") as fh:
        fh.write("barcode\tstatus\tassignment\tlog_prob_singleton\tlog_prob_doublet\tcluster0\tcluster1\n")
        for i, bc in enumerate(barcodes):
            cl = "0" if i < n0 else "1"
            # Sprinkle a few doublets/unassigned so propagation is exercised.
            if i % 60 == 0 and i > 0:
                fh.write("%s\tdoublet\t0/1\t-10\t-1\t-5\t-5\n" % bc)
            elif i % 61 == 0:
                fh.write("%s\tunassigned\t%s\t-2\t-2\t-2\t-2\n" % (bc, cl))
            else:
                fh.write("%s\tsinglet\t%s\t-1\t-10\t-1\t-20\n" % (bc, cl))
    with open(os.path.join(sdir, "ambient_rna.txt"), "w") as fh:
        fh.write("ambient RNA estimated as 4.20%\n")
    with open(os.path.join(sdir, "qc.txt"), "w") as fh:
        fh.write("4.20% ambient RNA\n0 doublets removed\n(synthetic fixture)\n")

    print("Fixture written to", sample_dir)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "..", "test_data", "sampleM"))
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    build(os.path.abspath(args.out), seed=args.seed)
