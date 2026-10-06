#!/usr/bin/env python3
"""Compact, self-contained per-group HTML report.

Leads with the maternal/fetal verdict, then the composition gate that allowed or
withheld it, per-library cell counts, cluster labels, marker evidence (RPS4Y1
highlighted) and provenance. Interpretable from tables/text alone — no reliance on
colour-coded plots.
"""
import argparse
import csv
import html


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def esc(x):
    return html.escape(str(x))


def table(rows, cols, highlight_col=None, highlight_val=None):
    if not rows:
        return "<p><em>No rows.</em></p>"
    head = "".join("<th>%s</th>" % esc(c) for c in cols)
    body = []
    for r in rows:
        cls = ""
        if highlight_col and r.get(highlight_col) == highlight_val:
            cls = ' class="hl"'
        cells = "".join("<td>%s</td>" % esc(r.get(c, "")) for c in cols)
        body.append("<tr%s>%s</tr>" % (cls, cells))
    return "<table><thead><tr>%s</tr></thead><tbody>%s</tbody></table>" % (
        head, "".join(body))


CONF_CLASS = {"high": "conf-high", "provisional": "conf-prov", "unresolved": "conf-unres"}

STYLE = """
body{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;margin:2rem;color:#1a1a1a;max-width:60rem}
h1{margin-bottom:0}h2{border-bottom:2px solid #eee;padding-bottom:.3rem;margin-top:2rem}
.verdict{padding:1rem 1.25rem;border-radius:8px;margin:1rem 0;border:1px solid #ddd}
.conf-high{background:#e6f4ea;border-color:#34a853}.conf-prov{background:#fef7e0;border-color:#f9ab00}
.conf-unres{background:#fce8e6;border-color:#ea4335}
table{border-collapse:collapse;width:100%;margin:.5rem 0;font-size:.9rem}
th,td{border:1px solid #ddd;padding:.35rem .6rem;text-align:left}th{background:#f5f5f5}
tr.hl{background:#fff3cd;font-weight:600}
.warn{background:#fce8e6;border-left:4px solid #ea4335;padding:.5rem .8rem;margin:.4rem 0}
.muted{color:#666;font-size:.85rem}
"""


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--group-id", required=True)
    ap.add_argument("--marker-summary", required=True)
    ap.add_argument("--cluster-labels", required=True)
    ap.add_argument("--qc-summary", required=True)
    ap.add_argument("--gate", required=True)
    ap.add_argument("--gate-libraries", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    markers = read_tsv(args.marker_summary)
    labels = read_tsv(args.cluster_labels)
    qc = read_tsv(args.qc_summary)
    qc0 = qc[0] if qc else {}
    gate = read_tsv(args.gate)[0]
    libraries = read_tsv(args.gate_libraries)

    conf = qc0.get("confidence", "unresolved")
    verdict = qc0.get("verdict", "unresolved")
    conf_cls = CONF_CLASS.get(conf, "conf-unres")
    warnings = [w for w in (qc0.get("warnings", "") or "").split(" | ") if w]
    label_of = {l["cluster"]: l["label"] for l in labels}

    parts = []
    parts.append("<h1>Maternal/Fetal Demultiplexing Report</h1>")
    parts.append('<p class="muted">Pregnancy group: <strong>%s</strong> &nbsp; Libraries: %d</p>'
                 % (esc(args.group_id), len(libraries)))

    # Verdict box
    parts.append('<div class="verdict %s">' % conf_cls)
    parts.append("<h2 style='margin-top:0;border:none'>Verdict: %s (%s confidence)</h2>"
                 % (esc(verdict), esc(conf)))
    parts.append("<p><strong>Composition gate:</strong> %s &nbsp; "
                 "<strong>Fetal cluster:</strong> %s &nbsp; "
                 "<strong>Maternal cluster:</strong> %s &nbsp; "
                 "<strong>Inferred fetal sex:</strong> %s</p>"
                 % (esc(gate["verdict"]), esc(qc0.get("fetal_cluster", "NA")),
                    esc(qc0.get("maternal_cluster", "NA")),
                    esc(qc0.get("inferred_sex", "NA"))))
    reason = next((l.get("reason") for l in labels if l.get("reason")), qc0.get("reason", ""))
    parts.append("<p><strong>Reason:</strong> %s</p>" % esc(reason))
    parts.append("</div>")

    if warnings:
        parts.append("<h2>Warnings</h2>")
        for w in warnings:
            parts.append('<div class="warn">%s</div>' % esc(w))

    # Composition gate
    parts.append("<h2>Composition gate</h2>")
    parts.append("<p>%s</p>" % esc(gate["reason"]))
    parts.append(table([gate], ["discordant_fraction", "separating_loci", "loci_called_in_both",
                                "cluster_singlets", "minority_cluster", "minority_singlets",
                                "n_cells", "n_doublet", "n_unassigned", "doublet_fraction"]))
    parts.append('<p class="muted">Decides on the discordant fraction (share of loci called '
                 'in both profiles whose genotypes differ; accept &ge; %s) and minority-cluster '
                 'singlets (&ge; %s). Reported as QC only: fitted ambient RNA %s%%, missing '
                 'genotype calls %s. See docs/pipeline_design.md.</p>'
                 % (esc(gate["min_discordant_fraction"]), esc(gate["min_minority_singlets"]),
                    esc(gate["ambient_rna_pct"]), esc(gate["missing_gt_fraction"])))

    # Per-library counts, cluster columns renamed to their labels
    parts.append("<h2>Cells per library</h2>")
    cluster_cols = [c for c in (libraries[0].keys() if libraries else [])
                    if c.startswith("singlets_cluster_")]
    lib_rows = []
    for row in libraries:
        out = dict(row)
        for col in cluster_cols:
            cluster = col[len("singlets_cluster_"):]
            out["cluster %s (%s)" % (cluster, label_of.get(cluster, "unresolved"))] = row[col]
        lib_rows.append(out)
    parts.append(table(lib_rows, ["sample_id", "n_cells", "n_singlet", "n_doublet", "n_unassigned"]
                       + ["cluster %s (%s)" % (c[len("singlets_cluster_"):],
                                               label_of.get(c[len("singlets_cluster_"):],
                                                            "unresolved"))
                          for c in cluster_cols]))
    parts.append('<p class="muted">Singlet counts per genotype cluster. A library can hold '
                 'zero cells of a cluster.</p>')

    parts.append("<h2>Cluster labels</h2>")
    parts.append(table(labels, ["cluster", "label", "confidence", "inferred_sex",
                                "rule", "supporting_markers", "n_singlet_cells"]))

    marker_cols = ["cluster", "marker", "tier", "n_cluster_cells", "n_pos",
                   "frac_pos", "cp10k", "mean_counts"]
    parts.append("<h2>Marker evidence by cluster (all libraries pooled)</h2>")
    parts.append('<p class="muted">Denominators explicit: frac_pos = positive singlets / '
                 'singlets in cluster; cp10k = library-size-normalised pseudobulk. '
                 'RPS4Y1 rows highlighted.</p>')
    parts.append(table(markers, marker_cols, highlight_col="marker", highlight_val="RPS4Y1"))

    parts.append('<h2>Provenance</h2>')
    parts.append('<p class="muted">Per-library allele counts (Souporcell remapping + Vartrix) '
                 'merged into one group matrix; Souporcell K=2 discovery on the group; '
                 'composition gate; labels from the version-controlled marker panel via the '
                 'two-tier decision engine (RPS4Y1 override → secondary vote → unresolved). '
                 'Thresholds are recorded in the run parameters (pipeline_info/). '
                 'Female-fetus calls are internally consistent but not independently '
                 'validated. See docs/pipeline_design.md.</p>')

    doc = "<!doctype html><html><head><meta charset='utf-8'><title>%s — maternal/fetal</title>" \
          "<style>%s</style></head><body>%s</body></html>" % (
              esc(args.group_id), STYLE, "".join(parts))
    with open(args.out, "w") as fh:
        fh.write(doc)
    print("[make_report] wrote %s" % args.out)


if __name__ == "__main__":
    main()
