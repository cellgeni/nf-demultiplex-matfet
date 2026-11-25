#!/usr/bin/env python3
"""Compact, self-contained per-fetus HTML report for fetal mode.

Leads with the result and its bound, states the range in which it is valid next to it, then
the warnings, the in-range evidence (sex and the XIST/Y check, and the Souporcell cross-library check
with its heatmaps when it ran), the per-batch table and the flagged cells. Reads fetal_summary.py's
outputs and, optionally, souporcell_distance.py's output directory.
"""
import argparse
import csv
import os

from make_report import STYLE, esc, table


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def num(x, fmt="{:,}"):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return ""
    if v != v:
        return ""
    return fmt.format(int(v)) if fmt == "{:,}" else fmt.format(v)


def souporcell_section(soc, gid):
    libs = [r for r in read_tsv(os.path.join(soc, "libraries.tsv")) if r["group"] == gid]
    out = ['<p>Souporcell (k = 2) on %d libraries of this fetus. Every library splits in two even with no mother; '
           'what matters is whether one cluster stands apart from the clusters of the <em>other</em> libraries. '
           'Reference distances per SNP: one person 0.000-0.017; mother vs her fetus 0.041-0.047; unrelated '
           '0.068-0.089. A cluster at about 0.04 from the rest would be a candidate mother and needs review '
           '(barcode, marker and XIST/Y checks); the check never changes the calls above.</p>'
           % len({r["library"] for r in libs})]
    svg = os.path.join(soc, "%s.heatmap.svg" % gid)
    if os.path.exists(svg):
        with open(svg) as fh:
            out.append('<div class="fig">%s</div>' % fh.read())
    rows = [{"library": r["library"], "cluster": r["cluster"], "singlets": num(r["singlets"]),
             "to sister cluster": num(r["msd_to_sister"], "{:.4f}"),
             "to other libraries (median)": num(r["msd_to_other_libraries_median"], "{:.4f}"),
             "to other libraries (max)": num(r["msd_to_other_libraries_max"], "{:.4f}"), "status": r["status"]}
            for r in libs]
    out.append(table(rows, ["library", "cluster", "singlets", "to sister cluster", "to other libraries (median)",
                            "to other libraries (max)", "status"], highlight_col="status", highlight_val="review"))
    svg = os.path.join(soc, "groups.heatmap.svg")
    if os.path.exists(svg):
        with open(svg) as fh:
            out.append('<p>All fetuses of this run. Each fetus should be light on the diagonal and dark against the '
                       'others; a fetus near 0.04 to another would point to a relative or a sample mix-up.</p>'
                       '<div class="fig">%s</div>' % fh.read())
    out.append('<p class="muted">Tables: fetal/souporcell_check/ (libraries.tsv, pairs.tsv, groups.tsv). '
               'Souporcell outputs: souporcell/&lt;sample_id&gt;/.</p>')
    return "".join(out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--group-id", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--sex-check", required=True)
    ap.add_argument("--warnings", required=True)
    ap.add_argument("--flagged", required=True)
    ap.add_argument("--souporcell-check", help="souporcell_distance.py output directory (optional)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    summary = read_tsv(args.summary)
    fetus = next(r for r in summary if r["scope"] == "fetus")
    batches = [r for r in summary if r["scope"] == "batch"]
    sex = read_tsv(args.sex_check)[0]
    with open(args.warnings) as fh:
        warnings = [w.strip() for w in fh if w.strip()]
    flagged = read_tsv(args.flagged)
    soc = args.souporcell_check
    if soc:
        warnings += [w["warning"] for w in read_tsv(os.path.join(soc, "warnings.tsv")) if w["group"] == args.group_id]

    conf_cls = "conf-prov" if warnings else "conf-high"
    parts = ["<h1>Fetal Mode Report</h1>",
             '<p class="muted">Fetus: <strong>%s</strong> &nbsp; Batches: %d</p>' % (esc(args.group_id), len(batches))]
    parts.append('<div class="verdict %s"><p><strong>%s maternal cells and %s fetal-maternal pairs flagged</strong> '
                 'among %s cells (removal candidates).</p>'
                 '<p>Cells cleared as fetal: %s (%s of cells); too little data to tell: %s, kept and not confirmed fetal.</p>'
                 '<p>One-sided 95%% upper bound on maternal cells among the cleared cells: <strong>%s</strong> '
                 '(3 &divide; cleared cells).</p>'
                 '<p class="muted">Valid only if maternal cells are a small share of each batch (benchmark: up to ~2%%). '
                 'Above ~5-10%% fetal mode flags nothing and this bound is wrong; see the in-range evidence below.</p></div>'
                 % (conf_cls, num(fetus["maternal"]), num(fetus["maternal_in_pair"]), num(fetus["cells"]), num(fetus["fetal"]),
                    num(fetus["fetal_fraction_of_cells"], "{:.1%}"), num(fetus["too_little_data"]),
                    num(fetus["upper95_maternal_among_fetal"], "{:.2e}")))

    parts.append("<h2>Warnings</h2>")
    parts.append("".join('<div class="warn">%s</div>' % esc(w) for w in warnings) if warnings else "<p>None.</p>")

    parts.append("<h2>In-range evidence</h2>")
    parts.append("<p>Sex call: <strong>%s</strong> (%s of cells carry Y-gene reads).</p>"
                 % (esc(sex["sex_call"]), num(sex["y_positive_fraction"], "{:.1%}")))
    if sex["sex_call"] == "male":
        parts.append("<p>XIST/Y check, independent of genotypes: %s cells with XIST &ge; 2 UMIs and no Y-gene reads; "
                     "a female cell at these depths shows XIST &ge; 2 with probability %s, so the maternal share is at most "
                     "<strong>%s</strong> (one-sided 95%%, all candidates counted as maternal).</p>"
                     % (num(sex["xist_candidates"]), num(sex["p_xist_detect"], "{:.3f}"),
                        num(sex["xist_upper95_maternal_share"], "{:.2e}")))
    elif sex["sex_call"] == "female":
        parts.append("<p>Female fetus: no genotype-free check of the maternal share exists, so the valid range is assumed "
                     "from the tissue, not measured.</p>")

    parts.append("<h2>Souporcell cross-library check</h2>")
    if not soc:
        parts.append("<p>Not run (<code>--fetal_souporcell off</code>).</p>")
    else:
        parts.append(souporcell_section(soc, args.group_id))

    parts.append("<h2>Batches</h2>")
    rows = [{"batch": b["id"], "cells": num(b["cells"]), "fetal": num(b["fetal"]),
             "too little data": num(b["too_little_data"]), "maternal": num(b["maternal"]),
             "in pair": num(b["maternal_in_pair"]), "flagged share": num(b["flagged_fraction"], "{:.3%}"),
             "upper 95% bound": num(b["upper95_maternal_among_fetal"], "{:.2e}")} for b in batches]
    parts.append(table(rows, ["batch", "cells", "fetal", "too little data", "maternal", "in pair", "flagged share",
                              "upper 95% bound"]))

    parts.append("<h2>Flagged cells</h2>")
    parts.append(table(flagged, ["library", "barcode", "fetal_mode", "n_sites", "llr", "llr_doublet", "rho_hat",
                                 "total_umis", "xist", "y_umis"]))
    parts.append('<p class="muted">Method: docs/fetal_mode_design.md. Per-cell calls: fetal/%s/cells.tsv.gz and '
                 '&lt;sample_id&gt;/fetal_mode_cells.tsv.gz.</p>' % esc(args.group_id))

    with open(args.out, "w") as fh:
        fh.write('<!doctype html><html><head><meta charset="utf-8"><title>Fetal mode: %s</title><style>%s%s</style>'
                 '</head><body>%s</body></html>' % (esc(args.group_id), STYLE, ".fig{overflow-x:auto;margin:.5rem 0}", "".join(parts)))


if __name__ == "__main__":
    main()
