#!/usr/bin/env python3
"""Compact, self-contained per-sample HTML report.

Leads with the verdict, then marker evidence (RPS4Y1 highlighted), Souporcell demux QC,
and provenance. Interpretable from tables/text alone — no reliance on colour-coded plots.
"""
import argparse
import csv
import html


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_text(path):
    try:
        with open(path) as fh:
            return fh.read()
    except OSError:
        return "(not available)"


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
pre{background:#f5f5f5;padding:.8rem;border-radius:6px;overflow-x:auto;font-size:.85rem}
.muted{color:#666;font-size:.85rem}
"""


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--marker-summary", required=True)
    ap.add_argument("--cluster-labels", required=True)
    ap.add_argument("--qc-summary", required=True)
    ap.add_argument("--souporcell-qc", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    markers = read_tsv(args.marker_summary)
    labels = read_tsv(args.cluster_labels)
    qc = read_tsv(args.qc_summary)
    qc0 = qc[0] if qc else {}
    soup_qc = read_text(args.souporcell_qc)

    conf = qc0.get("confidence", "unresolved")
    verdict = qc0.get("verdict", "unresolved")
    conf_cls = CONF_CLASS.get(conf, "conf-unres")

    warnings = [w for w in (qc0.get("warnings", "") or "").split(" | ") if w]

    # Order marker columns; keep RPS4Y1 first and highlight it.
    marker_cols = ["cluster", "marker", "tier", "n_cluster_cells", "n_pos",
                   "frac_pos", "cp10k", "mean_counts"]

    parts = []
    parts.append("<h1>Maternal/Fetal Demultiplexing Report</h1>")
    parts.append('<p class="muted">Sample: <strong>%s</strong></p>' % esc(args.sample_id))

    # Verdict box
    parts.append('<div class="verdict %s">' % conf_cls)
    parts.append("<h2 style='margin-top:0;border:none'>Verdict: %s (%s confidence)</h2>"
                 % (esc(verdict), esc(conf)))
    parts.append("<p><strong>Fetal cluster:</strong> %s &nbsp; "
                 "<strong>Maternal cluster:</strong> %s &nbsp; "
                 "<strong>Inferred fetal sex:</strong> %s</p>"
                 % (esc(qc0.get("fetal_cluster", "NA")),
                    esc(qc0.get("maternal_cluster", "NA")),
                    esc(qc0.get("inferred_sex", "NA"))))
    reason = next((l.get("reason") for l in labels if l.get("reason")), qc0.get("reason", ""))
    parts.append("<p><strong>Reason:</strong> %s</p>" % esc(reason))
    parts.append("</div>")

    # Warnings
    if warnings:
        parts.append("<h2>Warnings</h2>")
        for w in warnings:
            parts.append('<div class="warn">%s</div>' % esc(w))

    # Cluster labels
    parts.append("<h2>Cluster labels</h2>")
    parts.append(table(labels, ["cluster", "label", "confidence", "inferred_sex",
                                 "rule", "supporting_markers", "n_singlet_cells"]))

    # Marker evidence (RPS4Y1 highlighted)
    parts.append("<h2>Marker evidence by cluster</h2>")
    parts.append('<p class="muted">Denominators explicit: frac_pos = positive singlets / '
                 'singlets in cluster; cp10k = library-size-normalised pseudobulk. '
                 'RPS4Y1 rows highlighted.</p>')
    parts.append(table(markers, marker_cols, highlight_col="marker", highlight_val="RPS4Y1"))

    # Demux QC
    parts.append("<h2>Souporcell demultiplexing QC</h2>")
    parts.append("<p>Singlets: %s &nbsp; Doublets: %s &nbsp; Unassigned: %s &nbsp; "
                 "Cluster imbalance: %s &nbsp; Ambient RNA: %s</p>"
                 % (esc(qc0.get("n_singlets", "?")), esc(qc0.get("n_doublets", "?")),
                    esc(qc0.get("n_unassigned", "?")), esc(qc0.get("cluster_imbalance_ratio", "?")),
                    esc(qc0.get("ambient_rna", "?"))))
    parts.append("<pre>%s</pre>" % esc(soup_qc))

    parts.append('<h2>Provenance</h2>')
    parts.append('<p class="muted">Souporcell K=2. Labels derived from the version-controlled '
                 'marker panel via the two-tier decision engine (RPS4Y1 override → secondary '
                 'vote → unresolved). Thresholds are recorded in the run parameters '
                 '(pipeline_info/). Female-fetus calls are internally consistent but not '
                 'independently validated. See ROADMAP.md.</p>')

    doc = "<!doctype html><html><head><meta charset='utf-8'><title>%s — maternal/fetal</title>" \
          "<style>%s</style></head><body>%s</body></html>" % (
              esc(args.sample_id), STYLE, "".join(parts))
    with open(args.out, "w") as fh:
        fh.write(doc)
    print("[make_report] wrote %s" % args.out)


if __name__ == "__main__":
    main()
