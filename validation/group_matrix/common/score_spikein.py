#!/usr/bin/env python3
"""Score the fetal-fraction detection-limit experiment.

For each synthetic group and backend, two separate questions:

1. **Detection** -- did K=2 discovery recover the spiked fetal cells?
   recall     spiked fetal cells called as confident singlets on the fetal
              profile, over all spiked fetal cells
   precision  truly fetal cells among those confidently called fetal
   The fetal profile is whichever profile best matches the spiked cells, since
   profile names carry no meaning. In the 0% control there is nothing to
   recall, so the number of cells called to the smaller profile is reported
   as the false-positive count instead.

2. **Composition gate** (where split halves were run) -- would the QC accept a
   second genotype here? Separating-locus replication between the two halves,
   the confident-cell fraction, and the loci separating the two profiles.

A backend can detect the fetus while the gate still rejects it, or the reverse.
Reporting both keeps that disagreement visible, and it is what sets the gate's
thresholds on evidence rather than on three calibration points.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from assignment_table import is_high_confidence, read_table  # noqa: E402
from genotype_vcf import discordant_locus_replication, within_run_discordance  # noqa: E402

PROFILE_VCF = {
    "vireo": os.path.join("vireo_out", "GT_donors.vireo.vcf.gz"),
    "souporcell": os.path.join("consensus", "cluster_genotypes.vcf"),
}


def read_truth(path: str) -> dict[str, str]:
    with open(path, newline="") as handle:
        return {row["group_barcode"]: row["truth"]
                for row in csv.DictReader(handle, delimiter="\t")}


def detection(rows: list[dict], truth: dict[str, str]) -> dict:
    confident = [row for row in rows if is_high_confidence(row)]
    profiles: dict[str, int] = {}
    for row in confident:
        profiles[row["profile"]] = profiles.get(row["profile"], 0) + 1

    n_fetal = sum(1 for v in truth.values() if v == "fetal")
    n_cells = len(rows)
    result = {
        "n_cells": n_cells,
        "n_spiked_fetal": n_fetal,
        "n_high_confidence": len(confident),
        "high_confidence_fraction": len(confident) / n_cells if n_cells else 0.0,
        "profiles": profiles,
        "n_doublet": sum(1 for r in rows if r["status"] == "doublet"),
        "n_unassigned": sum(1 for r in rows if r["status"] == "unassigned"),
    }

    if len(profiles) < 2:
        result.update({"fetal_profile": None, "recall": 0.0 if n_fetal else None,
                       "precision": None, "false_positive_cells": 0,
                       "note": "fewer than two confident profiles"})
        return result

    if n_fetal == 0:
        smaller = min(profiles, key=profiles.get)
        result.update({"fetal_profile": None, "recall": None, "precision": None,
                       "false_positive_cells": profiles[smaller],
                       "note": "0% control: cells on the smaller profile are false positives"})
        return result

    best = None
    for profile in profiles:
        called = [r for r in confident if r["profile"] == profile]
        tp = sum(1 for r in called if truth.get(r["group_barcode"]) == "fetal")
        recall = tp / n_fetal
        precision = tp / len(called) if called else 0.0
        f1 = (2 * recall * precision / (recall + precision)) if recall + precision else 0.0
        candidate = {"fetal_profile": profile, "recall": recall, "precision": precision,
                     "f1": f1, "true_positive": tp, "called_fetal": len(called)}
        if best is None or f1 > best["f1"]:
            best = candidate
    result.update(best)

    fates = {"fetal_singlet": 0, "maternal_singlet": 0, "doublet": 0,
             "unassigned": 0, "low_confidence": 0}
    for row in rows:
        if truth.get(row["group_barcode"]) != "fetal":
            continue
        if row["status"] == "doublet":
            fates["doublet"] += 1
        elif row["status"] == "unassigned":
            fates["unassigned"] += 1
        elif not is_high_confidence(row):
            fates["low_confidence"] += 1
        elif row["profile"] == best["fetal_profile"]:
            fates["fetal_singlet"] += 1
        else:
            fates["maternal_singlet"] += 1
    result["fetal_cell_fates"] = fates
    return result


def gate(base_dir: str, backend: str) -> dict:
    """Composition-gate statistics, if split halves were run for this backend."""
    vcf_a = os.path.join(base_dir, f"{backend}_A", PROFILE_VCF[backend])
    vcf_b = os.path.join(base_dir, f"{backend}_B", PROFILE_VCF[backend])
    if not (os.path.isfile(vcf_a) and os.path.isfile(vcf_b)):
        return {"evaluated": False}
    rep = discordant_locus_replication(vcf_a, vcf_b)
    return {
        "evaluated": True,
        "jaccard": rep.get("jaccard"),
        "fraction_replicated": rep.get("fraction_of_a_replicated"),
        "n_separating_a": rep.get("n_discordant_a"),
        "n_separating_b": rep.get("n_discordant_b"),
    }


def score(base: str) -> list[dict]:
    rows_out = []
    with open(os.path.join(base, "spikein_design.tsv"), newline="") as handle:
        design = list(csv.DictReader(handle, delimiter="\t"))

    for d in design:
        gdir = os.path.join(base, "groups", d["label"])
        truth = read_truth(os.path.join(gdir, "truth.tsv"))
        for backend in ("vireo", "souporcell"):
            run_dir = os.path.join(gdir, f"{backend}_full")
            table = os.path.join(run_dir, "assignments.tsv")
            entry = {"label": d["label"], "backend": backend,
                     "fraction": float(d["fraction_target"]),
                     "replicate": int(d["replicate"])}
            if not os.path.isfile(table):
                entry["error"] = "no assignments produced"
                rows_out.append(entry)
                continue
            entry["detection"] = detection(read_table(table), truth)
            vcf = os.path.join(run_dir, PROFILE_VCF[backend])
            if os.path.isfile(vcf):
                entry["separating_loci"] = within_run_discordance(vcf).get(
                    "min_discordant_loci_over_pairs")
            entry["gate"] = gate(gdir, backend)
            rows_out.append(entry)
    return rows_out


def write_curve(rows: list[dict], path: str) -> None:
    header = ["backend", "fraction", "replicate", "n_spiked_fetal", "recall",
              "precision", "false_positive_cells", "high_confidence_fraction",
              "separating_loci", "gate_jaccard", "gate_fraction_replicated", "error"]
    with open(path, "w") as out:
        out.write("\t".join(header) + "\n")
        for r in sorted(rows, key=lambda x: (x["backend"], x["fraction"], x["replicate"])):
            det = r.get("detection", {})
            g = r.get("gate", {})

            def fmt(v):
                return "" if v is None else (f"{v:.4f}" if isinstance(v, float) else str(v))
            values = [r["backend"], f"{r['fraction']:g}", str(r["replicate"]),
                      fmt(det.get("n_spiked_fetal")), fmt(det.get("recall")),
                      fmt(det.get("precision")), fmt(det.get("false_positive_cells")),
                      fmt(det.get("high_confidence_fraction")), fmt(r.get("separating_loci")),
                      fmt(g.get("jaccard")), fmt(g.get("fraction_replicated")),
                      r.get("error", "")]
            out.write("\t".join(values) + "\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--base", required=True, help="spike-in experiment directory")
    args = parser.parse_args(argv)

    rows = score(args.base)
    out_dir = os.path.join(args.base, "report")
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "spikein_scores.json"), "w") as handle:
        json.dump(rows, handle, indent=2, sort_keys=True)
        handle.write("\n")
    curve = os.path.join(out_dir, "detection_curve.tsv")
    write_curve(rows, curve)
    with open(curve) as handle:
        print(handle.read())
    print(f"reports: {out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
