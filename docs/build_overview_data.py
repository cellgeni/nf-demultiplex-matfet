#!/usr/bin/env python3
"""Regenerate the data block of docs/pipeline_overview.html.

Reads a finished nf-demultiplex-matfet run (gate, labels and per-library summaries), the
validation runs under the REQ-68805 group_matrix directory (re-scored with the
pipeline's own composition_gate functions), the spike-in detection curve, and the
library metadata table, and rewrites the page's <script id="pipeline-data"> block in
place. Paths are for REQ-68805; edit them for another ticket.

    python3 docs/build_overview_data.py
"""
import csv
import glob
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "bin"))
from composition_gate import parse_ambient, profile_discordance  # noqa: E402

TICKET = "/lustre/scratch127/cellgen/cellgeni/tickets/REQ-68805"
RESULTS = f"{TICKET}/results_all"
VALIDATION = f"{TICKET}/group_matrix"
PERLIB = f"{TICKET}/results/souporcell"
ACTIONS = "/nfs/cellgeni/tickets/REQ-68805/actions"
SAMPLESHEET = f"{ACTIONS}/samplesheet_all.csv"
METADATA = f"{ACTIONS}/Souporcell_ECTatlas_scRNAseq_libraries - scRNAseq_by_library.csv"


def tsv_row(path):
    return next(csv.DictReader(open(path), delimiter="\t"))


def singlet_minority(clusters_tsv):
    counts = {}
    for r in csv.DictReader(open(clusters_tsv), delimiter="\t"):
        if r["status"] == "singlet":
            counts[r["assignment"]] = counts.get(r["assignment"], 0) + 1
    return sorted(counts.values())[0] if len(counts) >= 2 else 0


def validation_runs():
    runs = []
    for d in sorted(glob.glob(f"{VALIDATION}/spikein/groups/*/souporcell_*")):
        grp, arm = d.split("/")[-2], d.split("/")[-1]
        f = float(re.match(r"f([0-9p]+)_r", grp).group(1).replace("p", "."))
        truth = "fake" if f == 0 else ("weak" if f == 0.5 else "real")
        runs.append((f"spike {grp} {arm.replace('souporcell_', '')}", f"{d}/consensus", f"{d}/clusters.tsv", truth))
    for n, t in [("joint_seed1", "real"), ("joint_seed2", "real"), ("joint_seed3", "real"),
                 ("anchor_discovery", "real"), ("d6_forced_k2", "unknown")]:
        d = f"{VALIDATION}/gate2_souporcell/{n}"
        runs.append((f"CW5/D6 {n}", f"{d}/consensus", f"{d}/clusters.tsv", t))
    for n in ["discovery_all5", "discovery_balanced2k", "discovery_trio_early",
              "discovery_trio_late", "discovery_trio_mixed"]:
        d = f"{VALIDATION}/gate4_souporcell/{n}"
        runs.append((f"CW4 {n}", f"{d}/consensus", f"{d}/clusters.tsv", "fake"))
    for d in sorted(glob.glob(f"{PERLIB}/*")):
        s = os.path.basename(d)
        t = "real" if s.endswith("362") else ("unknown" if s.startswith("FCA") else "fake")
        runs.append((f"perlib {s}", d, f"{d}/clusters.tsv", t))
    out = []
    for name, cdir, clusters, truth in runs:
        vcf = f"{cdir}/cluster_genotypes.vcf"
        if not os.path.exists(vcf):
            continue
        dis = profile_discordance(vcf)
        amb = parse_ambient(f"{cdir}/ambient_rna.txt")
        called = dis["loci_called_in_both"]
        out.append({"run": name, "truth": truth, "ambient": None if amb is None else round(amb, 1),
                    "minority": singlet_minority(clusters),
                    "frac": round(dis["separating_loci"] / called, 3) if called else None})
    return out


def production():
    md = {r["Library_id"]: r for r in csv.DictReader(open(METADATA, encoding="utf-8-sig"))}
    sheet = list(csv.DictReader(open(SAMPLESHEET)))
    groups, libs = [], []
    for g in dict.fromkeys(r["group_id"] for r in sheet):
        gd = f"{RESULTS}/groups/{g}"
        members = [r["sample_id"] for r in sheet if r["group_id"] == g]
        if not os.path.exists(f"{gd}/gate.tsv"):
            groups.append({"group": g, "pending": True, "n": len(members)})
            continue
        gate = tsv_row(f"{gd}/gate.tsv")
        labels = {l["cluster"]: l for l in csv.DictReader(open(f"{gd}/cluster_labels.tsv"), delimiter="\t")}
        fetal = [c for c, l in labels.items() if l["label"].startswith("fetal")]
        first = labels[fetal[0]] if fetal else next(iter(labels.values()))
        m = md[members[0]]
        frac = gate["discordant_fraction"]
        groups.append({
            "group": g, "n": len(members), "frac": None if frac == "NA" else round(float(frac), 3),
            "ambient": round(float(gate["ambient_rna_pct"]), 1), "minority": int(gate["minority_singlets"]),
            "cells": int(gate["n_cells"]), "verdict": gate["verdict"],
            "label": first["label"] if fetal else "unresolved", "conf": first["confidence"], "rule": first["rule"],
            "dataset": m["Dataset_id"], "chem": m["10X_chemistry"], "tissue": m["Tissue_general"]})
        per_lib = {r["sample_id"]: r for r in csv.DictReader(open(f"{gd}/gate_libraries.tsv"), delimiter="\t")}
        for s in members:
            x, mm = per_lib[s], md[s]
            sing = int(x["n_singlet"])
            fc = sum(int(x.get(f"singlets_cluster_{c}", 0)) for c in fetal)
            libs.append({"lib": s, "group": g, "region": mm["Tissue_region"], "sort": mm["Cell_enrichment"],
                         "tissue": mm["Tissue_general"], "singlets": sing, "fetal": fc,
                         "pct": round(100 * fc / sing, 1) if sing else 0.0, "oriented": bool(fetal),
                         "verdict": gate["verdict"]})
    return groups, libs


def spikein(validation):
    full = {v["run"].split()[1]: v for v in validation if v["run"].startswith("spike") and v["run"].endswith("full")}
    out = []
    for r in csv.DictReader(open(f"{VALIDATION}/spikein/report/detection_curve.tsv"), delimiter="\t"):
        if r["backend"] != "souporcell":
            continue
        f, rep = float(r["fraction"]), int(r["replicate"])
        key = ("f0" if f == 0 else "f" + ("%g" % (f * 100)).replace(".", "p")) + f"_r{rep}"
        v = full.get(key, {})
        out.append({"fraction": f, "rep": rep, "n_fetal": int(r["n_spiked_fetal"]),
                    "recall": float(r["recall"]) if r["recall"] else None,
                    "precision": float(r["precision"]) if r["precision"] else None,
                    "ambient": v.get("ambient"), "frac": v.get("frac")})
    return out


if __name__ == "__main__":
    val = validation_runs()
    groups, libs = production()
    data = {"validation": val, "groups": groups, "libraries": libs, "spikein": spikein(val)}
    block = ('<script id="pipeline-data">\n'
             "// Generated by docs/build_overview_data.py; do not edit by hand.\n"
             "window.PIPELINE_DATA = " + json.dumps(data, separators=(",", ":")) + ";\n"
             "</script>")
    page = os.path.join(HERE, "pipeline_overview.html")
    html = open(page).read()
    html, n = re.subn(r'<script id="pipeline-data">.*?</script>', lambda _: block, html, flags=re.DOTALL)
    if n != 1:
        sys.exit(f"ERROR: expected one pipeline-data block in {page}, found {n}")
    open(page, "w").write(html)
    print(f"updated {page}", file=sys.stderr)
