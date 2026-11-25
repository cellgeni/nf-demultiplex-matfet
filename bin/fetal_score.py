#!/usr/bin/env python3
"""Fetal mode: flag maternal cells in libraries from one fetus, with no maternal sample.

Port of the frozen prototype v3 (score_v3.py) that scored the
F147, F158 and F218 embryos and the D8 public benchmark. The scoring is unchanged; only the
interface is: libraries and their cellsnp-lite directories, the variant panel and the output
prefix are arguments instead of hard-coded paths and environment variables, Y-gene counts are
read from each library's expression matrix (optional) instead of a fixed STARsolo folder, and
the atlas cell-group annotation is gone. tests/regression_fetal_score.py checks the port
against the prototype's outputs.

Method, per batch of libraries from one fetus (design: docs/fetal_mode_design.md):
- Fetal genotype: every site's genotype is called from all cells pooled, which assumes the
  batch is mostly fetal (valid up to about 2% maternal cells; blind above about 5-10%).
- Read model: which alleles appear among a cell's UMIs at a site (ref only, alt only, both),
  with the per-cell allele fraction drawn from a Beta around the expected value; the Beta
  concentration phi is fitted on the data for fetal-het and homozygous sites.
- Maternal hypothesis: the mother's genotype is unknown but shares one allele with the fetus;
  it is marginalised over her possible genotypes given the panel frequency, with a share rho
  of fetal background RNA (best of --rhos kept); a fetal-maternal pair is rho = 0.5.
- Cut-offs: the larger of 10 and the quantile giving at most 1 false flag per 10,000 cells in
  simulated fetal cells at each cell's own coverage.
- Status: maternal / fetal_maternal_doublet when above the cut-off; otherwise
  fetal_consistent with >= --min-sites informative sites, else insufficient_evidence.

Writes <prefix>.cells.tsv.gz, <prefix>.removed_sites.tsv and <prefix>.summary.tsv.
"""
import argparse
import gzip
import os

import numpy as np
import pandas as pd
from scipy.special import betaln

Y_GENES = ["RPS4Y1", "DDX3Y", "UTY", "KDM5D", "EIF1AY", "USP9Y", "ZFY"]
MIN_BULK_DP = 20      # pooled UMIs needed to use a site
TARGET_FPR = 1e-4     # at most 1 false flag per 10,000 fetal cells
LLR_FLOOR = 10.0      # never flag below a likelihood ratio of e^10
N_NULL = 20           # simulated fetal replicates per cell
N_MOTHERS = 10        # simulated mothers (maternal replicates per cell)
MERGE_BP = 10         # sites closer than this are one event; keep the deepest
DISP_MAX = 3.0        # drop fetal-het sites with dispersion (chi2 / cells) above this
DISP_MIN_CELLS = 20   # ...measured over at least this many cells with >= 2 UMIs
RHO_DOUBLET = 0.5
RHO_SIM_MAX = 0.1     # simulated mothers carry U(0, this) fetal background
MU_MIN = 1e-4         # expected allele fractions are kept within [MU_MIN, 1 - MU_MIN]
PHI_GRID = np.array([0.3, 0.5, 1, 2, 3, 5, 10, 20, 50, 100, 300, 1000], dtype=float)


def read_mtx(path):
    with open(path) as fh:
        line = fh.readline()
        while line.startswith("%"):
            line = fh.readline()
        ns, nc, _ = map(int, line.split())
        m = pd.read_csv(fh, sep="\t", header=None, names=["s", "c", "v"], dtype=np.int32)
    return ns, nc, m


def first_existing(d, name):
    for f in (f"{d}/{name}.gz", f"{d}/{name}"):
        if os.path.exists(f):
            return f
    raise SystemExit(f"expression directory {d}: no {name}(.gz)")


def y_counts(lib, d):
    """Y-gene UMIs per barcode from a CellRanger/STARsolo filtered matrix directory."""
    genes = pd.read_csv(first_existing(d, "features.tsv"), sep="\t", header=None)[1].tolist()
    bcs = pd.read_csv(first_existing(d, "barcodes.tsv"), header=None)[0].tolist()
    yidx = {i + 1 for i, x in enumerate(genes) if x in Y_GENES}
    yc = np.zeros(len(bcs) + 1)
    mf = first_existing(d, "matrix.mtx")
    with (gzip.open(mf, "rt") if mf.endswith(".gz") else open(mf)) as fh:
        line = fh.readline()
        while line.startswith("%"):
            line = fh.readline()
        for ch in pd.read_csv(fh, sep=" ", header=None, names=["g", "c", "v"], chunksize=5_000_000):
            ch = ch[ch.g.isin(yidx)]
            np.add.at(yc, ch.c.values, ch.v.values)
    return {(lib, b): v for b, v in zip(bcs, yc[1:])}


def pairs(values, what):
    out = []
    for v in values or []:
        if "=" not in v:
            raise SystemExit(f"--{what} expects LIBRARY=PATH, got {v!r}")
        out.append(tuple(v.split("=", 1)))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--library", action="append", required=True, metavar="LIB=CELLSNP_DIR",
                    help="a library and its cellsnp-lite output directory; repeat, in batch order")
    ap.add_argument("--expression", action="append", metavar="LIB=MATRIX_DIR",
                    help="optional filtered expression matrix per library, for the y_umis column")
    ap.add_argument("--panel", required=True, help="common-variant panel VCF with AF= in INFO")
    ap.add_argument("--min-sites", type=int, default=1000,
                    help="informative sites a cell needs to be fetal_consistent (default 1000)")
    ap.add_argument("--rhos", default="0,0.1,0.25",
                    help="fetal background shares tried for a maternal cell (default 0,0.1,0.25)")
    ap.add_argument("--prefix", required=True, help="output file prefix")
    args = ap.parse_args(argv)

    libraries = pairs(args.library, "library")
    expression = dict(pairs(args.expression, "expression"))
    LIBS = [lib for lib, _ in libraries]
    if len(set(LIBS)) != len(LIBS):
        raise SystemExit(f"duplicate library in --library: {LIBS}")
    RHOS = [float(x) for x in args.rhos.split(",")]
    MIN_SITES = args.min_sites
    rng = np.random.default_rng(1)

    # --- load and align all libraries on a shared site index ------------------------------
    site_id, cell_tab, ent = {}, [], []
    c0 = 0
    for lib, cs in libraries:
        vcf = pd.read_csv(f"{cs}/cellSNP.base.vcf.gz", sep="\t", comment="#", header=None,
                          usecols=[0, 1, 3, 4], names=["chrom", "pos", "ref", "alt"])
        keys = list(zip(vcf.chrom, vcf.pos, vcf.ref, vcf.alt))
        gidx = np.array([site_id.setdefault(k, len(site_id)) for k in keys], dtype=np.int64)
        ns, nc, dp = read_mtx(f"{cs}/cellSNP.tag.DP.mtx")
        _, _, ad = read_mtx(f"{cs}/cellSNP.tag.AD.mtx")
        kd = dp.s.values.astype(np.int64) * (nc + 1) + dp.c.values
        ka = ad.s.values.astype(np.int64) * (nc + 1) + ad.c.values
        o = np.argsort(kd); kd, dpv, ds, dc = kd[o], dp.v.values[o], dp.s.values[o], dp.c.values[o]
        pos = np.searchsorted(kd, ka)
        adv = np.zeros(len(kd), dtype=np.int32); adv[pos] = ad.v.values
        bcs = pd.read_csv(f"{cs}/cellSNP.samples.tsv", header=None)[0].values
        cell_tab.append(pd.DataFrame({"library": lib, "barcode": bcs}))
        ent.append(pd.DataFrame({"g": gidx[ds - 1], "c": dc - 1 + c0, "a": adv, "n": dpv}))
        c0 += nc
        print(f"{lib}: {nc:,} cells, {len(kd):,} cell-site readings", flush=True)
    cells = pd.concat(cell_tab, ignore_index=True)
    E = pd.concat(ent, ignore_index=True)
    del ent
    NC, NS = len(cells), len(site_id)
    keys = pd.DataFrame(list(site_id.keys()), columns=["chrom", "pos", "ref", "alt"])

    # --- pooled fetal genotype ---------------------------------------------------------------
    bdp = np.bincount(E.g, weights=E.n, minlength=NS)
    bad = np.bincount(E.g, weights=E.a, minlength=NS)
    af = (bad + 0.5) / (bdp + 1.0)
    gf = np.full(NS, -1)
    gf[af <= 0.03] = 0; gf[(af >= 0.3) & (af <= 0.7)] = 1; gf[af >= 0.97] = 2
    use = (bdp >= MIN_BULK_DP) & (gf >= 0) & keys.chrom.str.fullmatch(r"chr\d+").values

    need = set(zip(keys.chrom[use], keys.pos[use]))
    qd = {}
    for ch in pd.read_csv(args.panel, sep="\t", comment="#", header=None, usecols=[0, 1, 7],
                          names=["chrom", "pos", "info"], chunksize=1_000_000):
        ch = ch[[k in need for k in zip(ch.chrom, ch.pos)]]
        for c, p, info in ch.itertuples(index=False):
            qd[(c, p)] = float(info.split("AF=", 1)[1].split(";", 1)[0])
    q = np.array([qd.get(k, np.nan) for k in zip(keys.chrom, keys.pos)])
    use &= ~np.isnan(q)
    q = np.clip(np.nan_to_num(q, nan=0.5), 0.01, 0.99)
    print(f"sites: {NS:,} covered in any library, {use.sum():,} used "
          f"(hom-ref {((gf == 0) & use).sum():,}, het {((gf == 1) & use).sum():,}, hom-alt {((gf == 2) & use).sum():,})", flush=True)

    # --- site filters ---------------------------------------------------------------------------
    # 1. Over-dispersed fetal-het sites. Under the fetal model every cell's alt count at a site is
    #    binomial around the pooled fraction; paralog pseudo-SNPs, imprinted and monoallelically
    #    expressed genes spread far wider. Only het sites are tested: at hom sites a rare maternal
    #    cell is exactly the outlier we want to keep.
    Eu = E[use[E.g.values]]
    m2 = Eu[(Eu.n.values >= 2) & (gf[Eu.g.values] == 1)]
    p2 = af[m2.g.values]
    chi = (m2.a.values - m2.n.values * p2) ** 2 / (m2.n.values * p2 * (1 - p2))
    ncell = np.bincount(m2.g.values, minlength=NS)
    disp = np.bincount(m2.g.values, weights=chi, minlength=NS) / np.maximum(ncell, 1)
    overdisp = use & (gf == 1) & (ncell >= DISP_MIN_CELLS) & (disp > DISP_MAX)
    use &= ~overdisp
    del Eu, m2, p2, chi

    # 2. Adjacent sites: variants within MERGE_BP are one event (e.g. the four HBG1 positions);
    #    keep only the deepest of each run.
    ui = np.flatnonzero(use)
    kk = keys.iloc[ui].assign(d=bdp[ui]).sort_values(["chrom", "pos"])
    run = ((kk.chrom != kk.chrom.shift()) | (kk.pos - kk.pos.shift() > MERGE_BP)).cumsum().values
    keep = kk.groupby(run).d.idxmax().values
    adjacent = use.copy(); adjacent[keep] = False
    use &= ~adjacent
    removed = pd.concat([keys[overdisp].assign(reason="overdispersed_het", dispersion=disp[overdisp], cells=ncell[overdisp]),
                         keys[adjacent].assign(reason="adjacent_merged", dispersion=np.nan, cells=np.nan)])
    removed["pooled_af"], removed["bulk_dp"] = af[removed.index], bdp[removed.index]
    removed.sort_values(["reason", "dispersion"], ascending=[True, False]).to_csv(
        f"{args.prefix}.removed_sites.tsv", sep="\t", index=False)
    print(f"site filters: {overdisp.sum():,} over-dispersed het sites dropped "
          f"(of {((gf == 1) & (ncell >= DISP_MIN_CELLS)).sum():,} testable), {adjacent.sum():,} adjacent sites merged away; "
          f"{use.sum():,} sites left", flush=True)

    # P(G_m = k | G_f, q): mother shares one allele with the fetus
    prior = np.stack([(1 - q) ** 2, 2 * q * (1 - q), q ** 2], 1)
    lik = np.zeros_like(prior)
    gsafe = np.clip(gf, 0, 2)
    for k, mk in enumerate((0.0, 0.5, 1.0)):
        pfk = np.stack([(1 - mk) * (1 - q), mk * (1 - q) + (1 - mk) * q, mk * q], 1)
        lik[:, k] = pfk[np.arange(NS), gsafe]
    wm = prior * lik; wm /= wm.sum(1, keepdims=True)
    logw = np.log(np.clip(wm, 1e-300, None))

    # mother's background at her homozygous sites: what the fetal pool shows at its own
    err_hom = float(np.clip(np.median(np.minimum(af, 1 - af)[use & (gf != 1)]), 1e-3, 0.05))
    afm = np.tile([err_hom, 0.5, 1 - err_hom], (NS, 1))
    same = np.arange(3)[None, :] == gsafe[:, None]
    afm[same] = np.repeat(af[:, None], 3, 1)[same]          # mother's genotype = fetus's

    E = E[use[E.g.values]].reset_index(drop=True)
    g, c, a, n = E.g.values, E.c.values, E.a.values, E.n.values.astype(np.int64)
    del E
    pf = np.clip(af[g], MU_MIN, 1 - MU_MIN)
    het_f = gf[g] == 1

    def pattern(a_, n_):
        """0 = ref only, 1 = alt only, 2 = both alleles among the cell's UMIs at the site."""
        return np.where(a_ == 0, 0, np.where(a_ == n_, 1, 2)).astype(np.int8)

    def logp(pat, mu, phi):
        """log P(pattern | n UMIs, per-cell allele fraction ~ Beta(mean mu, concentration phi))."""
        mu = np.clip(mu, MU_MIN, 1 - MU_MIN)
        phi = np.broadcast_to(np.asarray(phi, dtype=float), mu.shape)
        al, be = mu * phi, (1 - mu) * phi
        one = n == 1
        lr = np.empty(len(n)); la = np.empty(len(n))
        lr[one] = np.log1p(-mu[one]); la[one] = np.log(mu[one])
        m = ~one
        b0 = betaln(al[m], be[m])
        lr[m] = betaln(al[m], be[m] + n[m]) - b0
        la[m] = betaln(al[m] + n[m], be[m]) - b0
        out = np.where(pat == 0, lr, la)
        both = pat == 2
        out[both] = np.log(np.clip(-np.expm1(np.logaddexp(lr[both], la[both])), 1e-300, None))
        return out

    def sim_pattern(mu, phi):
        mu = np.clip(mu, MU_MIN, 1 - MU_MIN)
        p = rng.beta(mu * phi, (1 - mu) * phi)
        return pattern(rng.binomial(n, p), n)

    pat_obs = pattern(a, n)

    def fit_phi(sel):
        idx = np.flatnonzero(sel & (n >= 2))
        if len(idx) > 2_000_000:
            idx = np.sort(rng.choice(idx, 2_000_000, replace=False))
        if len(idx) == 0:
            return 10.0
        sub_n = n[idx]; sub_p = pat_obs[idx]; sub_mu = pf[idx]
        ll = []
        for ph in PHI_GRID:
            mu = np.clip(sub_mu, MU_MIN, 1 - MU_MIN)
            al, be = mu * ph, (1 - mu) * ph
            b0 = betaln(al, be)
            lr = betaln(al, be + sub_n) - b0; la = betaln(al + sub_n, be) - b0
            lb = np.log(np.clip(-np.expm1(np.logaddexp(lr, la)), 1e-300, None))
            ll.append(np.where(sub_p == 0, lr, np.where(sub_p == 1, la, lb)).sum())
        return float(PHI_GRID[int(np.argmax(ll))])

    phi_het, phi_hom = fit_phi(het_f), fit_phi(~het_f)
    PHIK = np.array([phi_hom, phi_het, phi_hom])
    phi_f = np.where(het_f, phi_het, phi_hom)
    print(f"read model: phi het {phi_het:g}, phi hom {phi_hom:g}; mother's hom-site background {err_hom:.4f}", flush=True)

    def llr_parts(pat):
        """Per cell: maternal LLR (best rho in RHOS), that rho, doublet LLR, and the LLR at rho = 0.1."""
        lf = logp(pat, pf, phi_f)
        res = {}
        for rho in RHOS + [RHO_DOUBLET]:
            lm = np.full(len(n), -np.inf)
            for k in range(3):
                lm = np.logaddexp(lm, logw[g, k] + logp(pat, (1 - rho) * afm[g, k] + rho * pf, PHIK[k]))
            res[rho] = np.bincount(c, weights=lm - lf, minlength=NC)
        lr = np.vstack([res[r] for r in RHOS])
        return lr.max(0), np.array(RHOS)[lr.argmax(0)], res[RHO_DOUBLET], res[0.1]

    llr, rho_hat, llr_d, llr_10 = llr_parts(pat_obs)
    cells["n_sites"] = np.bincount(c, minlength=NC)
    cells["n_umis"] = np.bincount(c, weights=n, minlength=NC).astype(int)
    cells["llr"], cells["rho_hat"], cells["llr_doublet"], cells["llr_rho10"] = llr, rho_hat, llr_d, llr_10

    # --- null: simulated fetal cells at each cell's own coverage ---------------------------------
    null = np.empty((N_NULL, NC)); null_d = np.empty((N_NULL, NC))
    for r in range(N_NULL):
        null[r], _, null_d[r], _ = llr_parts(sim_pattern(pf, phi_f))
    # Separate cut-offs for the two hypotheses, each at half the target so that together a
    # fetal cell is flagged (as either) at most TARGET_FPR of the time.
    thr = max(LLR_FLOOR, float(np.quantile(null.ravel(), 1 - TARGET_FPR / 2)))
    thr_d = max(LLR_FLOOR, float(np.quantile(null_d.ravel(), 1 - TARGET_FPR / 2)))
    null_exceed = ((null > thr) | (null_d > thr_d)).mean()
    cells["null_q99"] = np.quantile(null, 0.99, axis=0)
    del null, null_d

    # --- positive control: simulated mothers, same read model, some fetal background -------------
    mat = np.empty((N_MOTHERS, NC)); hit = np.zeros(NC)
    for r in range(N_MOTHERS):
        u = rng.random(NS)
        gm = (u > wm[:, 0]).astype(int) + (u > wm[:, 0] + wm[:, 1]).astype(int)
        k = gm[g]
        rho_e = rng.uniform(0, RHO_SIM_MAX, NC)[c]
        m_, _, d_, _ = llr_parts(sim_pattern((1 - rho_e) * afm[g, k] + rho_e * pf, PHIK[k]))
        mat[r] = m_
        hit += (m_ > thr) | (d_ > thr_d)
    cells["maternal_sim_llr"] = np.median(mat, 0)
    cells["power"] = hit / N_MOTHERS
    del mat

    # --- status ---------------------------------------------------------------------------------
    st = np.where(cells.n_sites < MIN_SITES, "insufficient_evidence", "fetal_consistent").astype(object)
    st[(cells.llr_doublet > thr_d) & (cells.llr_doublet > cells.llr)] = "fetal_maternal_doublet"
    st[(cells.llr > thr) & (cells.llr >= cells.llr_doublet)] = "maternal"
    cells["status"] = st

    # --- batch-level share of maternal cells ------------------------------------------------------
    # Kept from the prototype for comparison only. The reported bound is the exact zero-event
    # bound pooled per fetus (fetal_summary.py); this likelihood-ratio bound is too optimistic at
    # zero events.
    def maternal_share(l):
        """MLE and one-sided 95% upper bound of the maternal share pi in a mixture of cells,
        each with likelihood ratio exp(l) for maternal (rho = 0.1) against fetal."""
        grid = np.concatenate([[0.0], np.logspace(-6, np.log10(0.5), 400)])
        ll = np.array([0.0] + [np.logaddexp(np.log1p(-p), np.log(p) + l).sum() for p in grid[1:]])
        i = int(ll.argmax())
        up = grid[(grid >= grid[i]) & (ll >= ll[i] - 1.353)].max()   # chi2_1(0.90) / 2
        return float(grid[i]), float(up)

    # --- Y-gene UMIs (male fetuses: a cell with Y reads is fetal) ---------------------------------
    ycount = {}
    for lib in LIBS:
        if lib in expression:
            ycount.update(y_counts(lib, expression[lib]))
    cells["y_umis"] = [ycount.get(k, 0) for k in zip(cells.library, cells.barcode)]
    cells.to_csv(f"{args.prefix}.cells.tsv.gz", sep="\t", index=False)

    srows = []
    for scope, sub in [("all", cells)] + list(cells.groupby("library")):
        deep = sub.n_sites >= MIN_SITES
        est, up = maternal_share(sub.llr_rho10.values[deep])
        est_all, up_all = maternal_share(sub.llr_rho10.values)
        srows.append({"scope": scope, "cells": len(sub), "deep_cells": int(deep.sum()),
                      "flagged": int(sub.status.isin(["maternal", "fetal_maternal_doublet"]).sum()),
                      "insufficient": int((sub.status == "insufficient_evidence").sum()),
                      "maternal_share_est": est, "maternal_share_upper95": up,
                      "all_cells_share_est": est_all, "all_cells_share_upper95": up_all})
    summ = pd.DataFrame(srows).assign(phi_het=phi_het, phi_hom=phi_hom, err_hom=err_hom, thr=thr, thr_d=thr_d,
                                      null_exceed=null_exceed)
    summ.to_csv(f"{args.prefix}.summary.tsv", sep="\t", index=False)

    # --- log ---------------------------------------------------------------------------------------
    ypos = cells.y_umis > 0
    print(f"\ncells {NC:,}; Y-positive {ypos.mean():.1%}")
    print(f"cut-offs: maternal LLR > {thr:.1f}, doublet LLR > {thr_d:.1f} (floor {LLR_FLOOR}); "
          f"simulated fetal cells above either: {null_exceed:.2e} ({N_NULL} x {NC:,} draws)")
    print("status counts:", cells.status.value_counts().to_dict())
    flag = cells.status.isin(["maternal", "fetal_maternal_doublet"])
    if ypos.any():
        print(f"flagged {flag.sum()} | Y-positive among flagged {int((flag & ypos).sum())} "
              f"| false-flag rate among Y+ cells {(flag & ypos).sum() / ypos.sum():.2e}")
    print("calibration, real cells vs their own simulated fetal twins: share above the twins' 99th percentile "
          f"{(cells.llr > cells.null_q99).mean():.4f} (expected about 0.01 if the read model fits)")
    if flag.any():
        print(cells[flag].sort_values("llr", ascending=False)
              [["library", "barcode", "n_sites", "llr", "rho_hat", "llr_doublet", "y_umis", "status"]].head(40).to_string(index=False))


if __name__ == "__main__":
    main()
