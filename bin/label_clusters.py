#!/usr/bin/env python3
"""Two-tier maternal/fetal cluster decision engine + barcode propagation.

Runs once per pregnancy group, on the group's Souporcell clusters, and only orients
clusters the composition gate accepted (verdict two_profiles); for any other verdict
every cluster stays unresolved. Design (see docs/pipeline_design.md, "Orientation"):
  Tier 1  RPS4Y1 override  -> fetal cluster, HIGH confidence, inferred male.
  Tier 2  secondary-marker vote -> fetal cluster, PROVISIONAL, likely female.
  Tier 3  unresolved       -> never forces a call; never labels from RPS4Y1 ABSENCE.

The decision functions (`build_cluster_stats`, `decide`) are pure and importable so
they can be unit-tested without Nextflow or a matrix (see tests/test_decision_engine.py).
Only singlets contribute evidence; doublets/unassigned inherit no maternal/fetal label.
"""
import argparse
import csv
import gzip
import math
import sys

EPS = 1e-9
RPS4Y1 = "RPS4Y1"

# Label vocabulary (docs/pipeline_design.md, "Outputs").
LABEL_FETAL_MALE = "fetal_high_confidence_male"
LABEL_FETAL_MARKER = "fetal_marker_supported"
LABEL_MATERNAL = "maternal_likely"
LABEL_UNRESOLVED = "unresolved"

# The only composition-gate verdict under which clusters are oriented.
GATE_ACCEPTED = "two_profiles"

DEFAULT_PARAMS = {
    "min_cluster_cells": 50,
    "min_marker_detected_fraction": 0.10,
    "min_marker_fraction_delta": 0.05,
    "min_marker_logfc": 1.0,
    "rps4y1_fraction_floor": 0.10,
    "rps4y1_min_fold": 5.0,
    "min_vote_margin": 2,
    "imbalance_warn_ratio": 10.0,
}


def _fold(a, b):
    return (a + EPS) / (b + EPS)


def _log2fc(a, b):
    return math.log2((a + EPS) / (b + EPS))


def build_cluster_stats(rows):
    """Build the structured stats dict the engine consumes from long-format rows.

    Each row is a dict with at least:
      cluster, n_cluster_cells, marker, tier, frac_pos, cp10k, n_pos
    Returns: {clusters, n_cells, markers, tiers}.
    """
    clusters = []
    n_cells = {}
    markers = {}
    tiers = {}
    for r in rows:
        c = str(r["cluster"])
        m = r["marker"]
        if c not in markers:
            markers[c] = {}
            clusters.append(c)
        n_cells[c] = int(float(r["n_cluster_cells"]))
        markers[c][m] = {
            "frac_pos": float(r["frac_pos"]),
            "cp10k": float(r["cp10k"]),
            "n_pos": int(float(r.get("n_pos", 0))),
        }
        tiers[m] = r.get("tier", "secondary")
    clusters = sorted(clusters, key=lambda x: (len(x), x))
    return {"clusters": clusters, "n_cells": n_cells, "markers": markers, "tiers": tiers}


def _fmt_pct(x):
    return f"{x * 100:.1f}%"


def decide(stats, params=None):
    """Apply the two-tier rule. Returns a result dict (see module docstring).

    Invariant: absence of RPS4Y1 can only move the decision to Tier 2 — it can never
    assign a maternal label. Maternal is only ever the *other* cluster once a fetal
    cluster is positively identified.
    """
    p = dict(DEFAULT_PARAMS)
    if params:
        p.update({k: v for k, v in params.items() if v is not None})

    clusters = stats["clusters"]
    n_cells = stats["n_cells"]
    markers = stats["markers"]
    tiers = stats["tiers"]
    warnings = []

    result = {
        "fetal_cluster": None,
        "maternal_cluster": None,
        "tier": "unresolved",
        "confidence": "unresolved",
        "inferred_sex": "undetermined",
        "rule": "no_rule_fired",
        "supporting_markers": [],
        "votes": {c: 0 for c in clusters},
        "warnings": warnings,
        "reason": "",
        "labels": {c: LABEL_UNRESOLVED for c in clusters},
    }

    # --- K=2 sanity gate ---
    if len(clusters) < 2:
        warnings.append(
            "Fewer than two singlet clusters present; the K=2 maternal/fetal "
            "assumption is not supported. No label made."
        )
        result["reason"] = "Only %d singlet cluster(s) present." % len(clusters)
        return result
    if len(clusters) > 2:
        warnings.append(
            "More than two singlet clusters present (%d); expected exactly 2 for K=2."
            % len(clusters)
        )

    # Work with the two largest clusters if more than two somehow appear.
    ordered = sorted(clusters, key=lambda c: n_cells.get(c, 0), reverse=True)
    a, b = ordered[0], ordered[1]

    # --- QC gates (warnings; downgrade HIGH -> PROVISIONAL, never fabricate a call) ---
    small = [c for c in (a, b) if n_cells.get(c, 0) < p["min_cluster_cells"]]
    if small:
        warnings.append(
            "Cluster(s) %s have < %d singlet cells; confidence downgraded."
            % (", ".join(small), p["min_cluster_cells"])
        )
    sizes = sorted([n_cells.get(a, 0), n_cells.get(b, 0)])
    imbalance = (sizes[1] + EPS) / (sizes[0] + EPS)
    if imbalance > p["imbalance_warn_ratio"]:
        warnings.append(
            "Cluster sizes are highly imbalanced (%.0fx); may indicate low fetal "
            "fraction or a single genotype." % imbalance
        )
    qc_fail_high = bool(small)

    def finalize(fetal, maternal, tier, confidence, sex, rule, supporting, reason):
        result.update(
            fetal_cluster=fetal,
            maternal_cluster=maternal,
            tier=tier,
            confidence=confidence,
            inferred_sex=sex,
            rule=rule,
            supporting_markers=list(supporting),
            reason=reason,
        )
        label_fetal = LABEL_FETAL_MALE if tier == "rps4y1" else LABEL_FETAL_MARKER
        result["labels"] = {c: LABEL_UNRESOLVED for c in clusters}
        result["labels"][fetal] = label_fetal
        result["labels"][maternal] = LABEL_MATERNAL
        return result

    # --- Tier 1: RPS4Y1 override (decisive when present & cluster-specific) ---
    def marker_of(cluster, m):
        return markers.get(cluster, {}).get(m, {"frac_pos": 0.0, "cp10k": 0.0, "n_pos": 0})

    ya, yb = marker_of(a, RPS4Y1), marker_of(b, RPS4Y1)
    if ya["frac_pos"] >= yb["frac_pos"]:
        cand, other, y_cand, y_other = a, b, ya, yb
    else:
        cand, other, y_cand, y_other = b, a, yb, ya

    both_y_positive = (
        y_cand["frac_pos"] >= p["rps4y1_fraction_floor"]
        and y_other["frac_pos"] >= p["rps4y1_fraction_floor"]
    )
    if both_y_positive:
        warnings.append(
            "RPS4Y1 detected in BOTH clusters (%s vs %s); possible ambient RNA, mixed "
            "cluster, or demux issue."
            % (_fmt_pct(y_cand["frac_pos"]), _fmt_pct(y_other["frac_pos"]))
        )

    y_fold = _fold(y_cand["frac_pos"], y_other["frac_pos"])
    if y_cand["frac_pos"] >= p["rps4y1_fraction_floor"] and y_fold >= p["rps4y1_min_fold"]:
        confidence = "provisional" if qc_fail_high else "high"
        reason = (
            "RPS4Y1 detected in %s of cluster %s and %s of cluster %s "
            "(%.1fx enrichment; CP10K %.1f vs %.1f)."
            % (
                _fmt_pct(y_cand["frac_pos"]), cand,
                _fmt_pct(y_other["frac_pos"]), other,
                y_fold, y_cand["cp10k"], y_other["cp10k"],
            )
        )
        return finalize(
            cand, other, "rps4y1", confidence, "male",
            "tier1_rps4y1_override", [RPS4Y1], reason,
        )

    # --- Tier 2: secondary-marker vote (female fetus / weak-Y) ---
    secondary = [m for m, t in tiers.items() if t == "secondary"]
    votes = {a: 0, b: 0}
    supporters = {a: [], b: []}
    for m in secondary:
        ma, mb = marker_of(a, m), marker_of(b, m)
        c_cand, c_other, m_cand, m_other = (
            (a, b, ma, mb) if ma["frac_pos"] >= mb["frac_pos"] else (b, a, mb, ma)
        )
        passes = (
            m_cand["frac_pos"] >= p["min_marker_detected_fraction"]
            and (m_cand["frac_pos"] - m_other["frac_pos"]) >= p["min_marker_fraction_delta"]
            and _log2fc(m_cand["cp10k"], m_other["cp10k"]) >= p["min_marker_logfc"]
        )
        if passes:
            votes[c_cand] += 1
            supporters[c_cand].append(m)
    result["votes"] = dict(votes)

    winner = a if votes[a] >= votes[b] else b
    loser = b if winner == a else a
    margin = votes[winner] - votes[loser]

    if votes[winner] >= 1 and margin >= p["min_vote_margin"]:
        # v1 caps secondary-marker calls at PROVISIONAL: female calls are internally
        # consistent but not independently validated.
        sup = supporters[winner]
        reason = (
            "%d fetal marker(s) enriched in cluster %s (%s); secondary-marker vote "
            "%d–%d." % (len(sup), winner, ", ".join(sup), votes[winner], votes[loser])
        )
        return finalize(
            winner, loser, "markers", "provisional", "likely_female",
            "tier2_secondary_marker_vote", sup, reason,
        )

    # --- Tier 3: unresolved ---
    if votes[a] == 0 and votes[b] == 0:
        result["reason"] = "No fetal marker passed the configured thresholds in either cluster."
    elif margin < p["min_vote_margin"]:
        result["reason"] = (
            "Secondary-marker votes did not reach the required margin of %d "
            "(cluster votes: %s)." % (p["min_vote_margin"], dict(votes))
        )
    return result


# ---------------------------------------------------------------------------
# CLI (I/O layer around the pure engine)
# ---------------------------------------------------------------------------
def _read_marker_summary(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _read_clusters(path):
    """Yield dicts with barcode/status/assignment from a Souporcell clusters.tsv."""
    with open(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        idx = {name: i for i, name in enumerate(header)}
        bc_i = idx.get("barcode", 0)
        st_i = idx.get("status", 1)
        as_i = idx.get("assignment", 2)
        for row in reader:
            if not row:
                continue
            yield {"barcode": row[bc_i], "status": row[st_i], "assignment": row[as_i]}


def _parse_ambient(path):
    if not path:
        return "NA"
    try:
        with open(path) as fh:
            text = fh.read()
        for tok in text.replace("%", " ").split():
            try:
                return "%.3f" % float(tok)
            except ValueError:
                continue
    except OSError:
        pass
    return "NA"


def _read_gate(path):
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if len(rows) != 1:
        sys.exit("ERROR: %s must hold exactly one gate row" % path)
    return rows[0]


def apply_gate(res, clusters, verdict, gate_reason):
    """Withhold every label unless the composition gate accepted two profiles."""
    if verdict == GATE_ACCEPTED:
        return res
    res = dict(res)
    res.update(
        fetal_cluster=None,
        maternal_cluster=None,
        tier="unresolved",
        confidence="unresolved",
        inferred_sex="undetermined",
        rule="composition_gate",
        supporting_markers=[],
        labels={c: LABEL_UNRESOLVED for c in clusters},
        reason="Composition gate verdict %s: %s" % (verdict, gate_reason),
    )
    return res


def _open_out(path):
    return gzip.open(path, "wt", newline="") if path.endswith(".gz") else open(path, "w", newline="")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--group-id", required=True)
    ap.add_argument("--marker-summary", required=True)
    ap.add_argument("--clusters", required=True)
    ap.add_argument("--ambient", default=None)
    ap.add_argument("--gate", required=True, help="gate.tsv from composition_gate.py")
    ap.add_argument("--cluster-labels-out", required=True)
    ap.add_argument("--cell-assignments-out", required=True)
    ap.add_argument("--qc-out", required=True)
    for k in DEFAULT_PARAMS:
        ap.add_argument("--" + k.replace("_", "-"), type=float, default=None)
    args = ap.parse_args(argv)

    params = {k: getattr(args, k) for k in DEFAULT_PARAMS}
    rows = _read_marker_summary(args.marker_summary)
    if not rows:
        sys.exit("ERROR: empty marker summary: %s" % args.marker_summary)

    stats = build_cluster_stats(rows)
    gate = _read_gate(args.gate)
    res = apply_gate(decide(stats, params), stats["clusters"], gate["verdict"], gate["reason"])

    warnings_str = " | ".join(res["warnings"]) if res["warnings"] else ""

    # cluster_labels.tsv
    with open(args.cluster_labels_out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "group_id", "cluster", "label", "confidence", "inferred_sex", "rule",
            "supporting_markers", "n_singlet_cells", "reason", "warnings",
        ])
        for c in stats["clusters"]:
            is_fetal = c == res["fetal_cluster"]
            w.writerow([
                args.group_id, c, res["labels"][c],
                res["confidence"] if res["labels"][c] != LABEL_UNRESOLVED else "unresolved",
                res["inferred_sex"] if is_fetal else ("undetermined" if res["fetal_cluster"] is None else "NA"),
                res["rule"],
                ";".join(res["supporting_markers"]) if is_fetal else "",
                stats["n_cells"].get(c, 0),
                res["reason"], warnings_str,
            ])

    # group_qc_summary.tsv
    counts = {"singlet": 0, "doublet": 0, "unassigned": 0}
    per_cluster = {}
    cluster_rows = list(_read_clusters(args.clusters))
    for r in cluster_rows:
        counts[r["status"]] = counts.get(r["status"], 0) + 1
        if r["status"] == "singlet":
            per_cluster[r["assignment"]] = per_cluster.get(r["assignment"], 0) + 1
    sizes = sorted(per_cluster.values())
    imbalance = "%.1f" % ((sizes[-1] + EPS) / (sizes[0] + EPS)) if len(sizes) >= 2 else "NA"

    with open(args.qc_out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "group_id", "gate_verdict", "verdict", "confidence", "inferred_sex",
            "fetal_cluster", "maternal_cluster", "n_singlets", "n_doublets", "n_unassigned",
            "n_singlet_clusters", "cluster_imbalance_ratio", "ambient_rna",
            "n_warnings", "reason", "warnings",
        ])
        w.writerow([
            args.group_id, gate["verdict"], res["tier"], res["confidence"], res["inferred_sex"],
            res["fetal_cluster"] if res["fetal_cluster"] is not None else "NA",
            res["maternal_cluster"] if res["maternal_cluster"] is not None else "NA",
            counts.get("singlet", 0), counts.get("doublet", 0), counts.get("unassigned", 0),
            len(per_cluster), imbalance, _parse_ambient(args.ambient),
            len(res["warnings"]), res["reason"], warnings_str,
        ])

    # cell_assignments.tsv.gz, keyed by group barcode; split_by_library.py resolves
    # each row to its library and raw barcode through barcode_provenance.tsv.
    with _open_out(args.cell_assignments_out) as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "group_barcode", "souporcell_status", "souporcell_cluster",
            "label", "confidence", "reason",
        ])
        for r in cluster_rows:
            if r["status"] == "singlet":
                label = res["labels"].get(r["assignment"], LABEL_UNRESOLVED)
                conf = res["confidence"] if label != LABEL_UNRESOLVED else "unresolved"
                reason = res["reason"]
            else:
                label = r["status"]  # 'doublet' or 'unassigned' — never coerced
                conf = "NA"
                reason = "Souporcell %s; no maternal/fetal label assigned." % r["status"]
            w.writerow([
                r["barcode"], r["status"], r["assignment"],
                label, conf, reason,
            ])

    print("[label_clusters] %s -> gate=%s verdict=%s confidence=%s fetal=%s" % (
        args.group_id, gate["verdict"], res["tier"], res["confidence"], res["fetal_cluster"]))


if __name__ == "__main__":
    main()
