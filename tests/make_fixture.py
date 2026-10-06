#!/usr/bin/env python3
"""Generate a tiny synthetic two-group fixture for the end-to-end test profile.

Creates under <out>/ (default test_data/fixture/):
  samplesheet.csv            relative paths, so the fixture can live anywhere
  panel.vcf                  reference-panel header (##contig lines only)
  <sample>/barcodes.tsv.gz
  <sample>/truth.tsv                           barcode -> maternal | fetal
  <sample>/possorted_genome_bam.bam            empty placeholder; never read, because
                                               every row supplies souporcell_dir
  <sample>/filtered_feature_bc_matrix/         features, barcodes, matrix
  <sample>/souporcell/{alt.mtx,ref.mtx,barcodes.tsv,common_variants_covered.vcf}

Groups:
  G1  libA (180 maternal + 70 fetal) and libB (200 maternal + 8 fetal): one
      mother and one male fetus, so the gate should accept two profiles and
      RPS4Y1 should orient them.
  G2  libC (240 cells of one genotype): Souporcell will still split it at K=2,
      and the gate should report it unresolved.

Allele counts mimic Vartrix: an entry for every covered (locus, cell) pair, with
explicit zeros, identical coordinates in ALT and REF. Deterministic (fixed seed).
"""
import argparse
import gzip
import os

import numpy as np
from scipy import io as sio
from scipy import sparse

PANEL = ["RPS4Y1", "HAPLN1", "PITX2", "GAPLINC", "HPGDS", "XAGE2", "XAGE3", "HLA-G"]
BACKGROUND = ["ACTB", "GAPDH", "MALAT1", "XIST", "RPL13", "TMSB4X"]
N_LOCI = 3000
COVER_FRACTION = 0.15     # loci observed per cell
ERROR = 0.01              # per-read allele error


def hwe_genotypes(rng, freqs):
    return rng.binomial(2, freqs)


def child_genotypes(rng, mother, father):
    from_mother = rng.binomial(1, mother / 2.0)
    from_father = rng.binomial(1, father / 2.0)
    return from_mother + from_father


def write_expression(sample_dir, rng, barcodes, is_fetal):
    genes = BACKGROUND + PANEL
    ids = ["ENSGFAKE%08d" % i for i in range(len(genes))]
    gi = {g: i for i, g in enumerate(genes)}
    n = len(barcodes)
    M = np.zeros((len(genes), n), dtype=int)
    for g in BACKGROUND:
        M[gi[g], :] = rng.poisson(20, n)
    fetal = np.flatnonzero(is_fetal)
    maternal = np.flatnonzero(~is_fetal)
    M[gi["RPS4Y1"], fetal] = rng.poisson(6, len(fetal)) * (rng.random(len(fetal)) < 0.55)
    M[gi["RPS4Y1"], maternal] = rng.poisson(1, len(maternal)) * (rng.random(len(maternal)) < 0.02)
    for g in PANEL[1:]:
        M[gi[g], fetal] = rng.poisson(4, len(fetal)) * (rng.random(len(fetal)) < 0.40)
        M[gi[g], maternal] = rng.poisson(1, len(maternal)) * (rng.random(len(maternal)) < 0.02)

    mdir = os.path.join(sample_dir, "filtered_feature_bc_matrix")
    os.makedirs(mdir, exist_ok=True)
    with gzip.open(os.path.join(mdir, "features.tsv.gz"), "wt") as fh:
        for i, g in enumerate(genes):
            fh.write("%s\t%s\tGene Expression\n" % (ids[i], g))
    with gzip.open(os.path.join(mdir, "barcodes.tsv.gz"), "wt") as fh:
        fh.write("\n".join(barcodes) + "\n")
    with gzip.open(os.path.join(mdir, "matrix.mtx.gz"), "wb") as fh:
        sio.mmwrite(fh, sparse.csr_matrix(M))


def write_alleles(sample_dir, rng, barcodes, genotypes, loci, alleles):
    """genotypes: cells x N_LOCI alt-allele copies. Only a library-specific subset
    of loci is 'covered', as in a real Souporcell run."""
    covered = np.flatnonzero(rng.random(N_LOCI) < 0.85)
    rows, cols, alt_vals, ref_vals = [], [], [], []
    for j in range(len(barcodes)):
        seen = covered[rng.random(len(covered)) < COVER_FRACTION]
        depth = rng.poisson(1.0, len(seen)) + 1
        p_alt = genotypes[j, seen] / 2.0 * (1 - 2 * ERROR) + ERROR
        alt = rng.binomial(depth, p_alt)
        row_index = np.searchsorted(covered, seen) + 1
        rows.extend(row_index.tolist())
        cols.extend([j + 1] * len(seen))
        alt_vals.extend(alt.tolist())
        ref_vals.extend((depth - alt).tolist())

    order = np.lexsort((np.array(rows), np.array(cols)))
    sdir = os.path.join(sample_dir, "souporcell")
    os.makedirs(sdir, exist_ok=True)
    for name, vals in (("alt.mtx", alt_vals), ("ref.mtx", ref_vals)):
        with open(os.path.join(sdir, name), "w") as fh:
            fh.write("%%MatrixMarket matrix coordinate real general\n% written by make_fixture.py\n")
            fh.write("%d %d %d\n" % (len(covered), len(barcodes), len(order)))
            for k in order:
                fh.write("%d %d %d\n" % (rows[k], cols[k], vals[k]))
    with open(os.path.join(sdir, "barcodes.tsv"), "w") as fh:
        fh.write("\n".join(barcodes) + "\n")
    with open(os.path.join(sdir, "common_variants_covered.vcf"), "w") as fh:
        fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for i in covered:
            ref, alt = alleles[i]
            fh.write("chr1\t%d\t.\t%s\t%s\t.\tPASS\t.\n" % (loci[i], ref, alt))


def build(out, seed=0):
    rng = np.random.default_rng(seed)
    os.makedirs(out, exist_ok=True)
    loci = [10000 + 500 * i for i in range(N_LOCI)]
    bases = "ACGT"
    alleles = []
    for _ in range(N_LOCI):
        ref = bases[rng.integers(4)]
        alt = rng.choice([b for b in bases if b != ref])
        alleles.append((ref, alt))
    freqs = rng.uniform(0.1, 0.5, N_LOCI)

    with open(os.path.join(out, "panel.vcf"), "w") as fh:
        fh.write("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=248956422>\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    mother1 = hwe_genotypes(rng, freqs)
    fetus1 = child_genotypes(rng, mother1, hwe_genotypes(rng, freqs))
    mother2 = hwe_genotypes(rng, freqs)

    libraries = [
        ("G1", "libA", [(mother1, 180, False), (fetus1, 70, True)]),
        ("G1", "libB", [(mother1, 200, False), (fetus1, 8, True)]),
        ("G2", "libC", [(mother2, 240, False)]),
    ]

    sheet = ["group_id,sample_id,bam,barcodes,expression,souporcell_dir"]
    for group_id, sample_id, parts in libraries:
        sample_dir = os.path.join(out, sample_id)
        os.makedirs(sample_dir, exist_ok=True)
        genotypes = np.vstack([np.tile(g, (n, 1)) for g, n, _f in parts])
        is_fetal = np.concatenate([np.full(n, f) for _g, n, f in parts])
        shuffle = rng.permutation(len(is_fetal))
        genotypes, is_fetal = genotypes[shuffle], is_fetal[shuffle]
        barcodes = ["".join(rng.choice(list(bases), 16)) + "-1" for _ in range(len(is_fetal))]
        assert len(set(barcodes)) == len(barcodes)

        write_expression(sample_dir, rng, barcodes, is_fetal)
        write_alleles(sample_dir, rng, barcodes, genotypes, loci, alleles)
        with gzip.open(os.path.join(sample_dir, "barcodes.tsv.gz"), "wt") as fh:
            fh.write("\n".join(barcodes) + "\n")
        with open(os.path.join(sample_dir, "truth.tsv"), "w") as fh:
            fh.write("barcode\tgenotype\n")
            for barcode, fetal in zip(barcodes, is_fetal):
                fh.write("%s\t%s\n" % (barcode, "fetal" if fetal else "maternal"))
        open(os.path.join(sample_dir, "possorted_genome_bam.bam"), "w").close()
        sheet.append(",".join([
            group_id, sample_id,
            "%s/possorted_genome_bam.bam" % sample_id,
            "%s/barcodes.tsv.gz" % sample_id,
            "%s/filtered_feature_bc_matrix" % sample_id,
            "%s/souporcell" % sample_id,
        ]))

    with open(os.path.join(out, "samplesheet.csv"), "w") as fh:
        fh.write("\n".join(sheet) + "\n")
    print("Fixture written to", out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "..", "test_data", "fixture"))
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    build(os.path.abspath(args.out), seed=args.seed)
