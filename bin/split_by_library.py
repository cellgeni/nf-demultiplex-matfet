#!/usr/bin/env python3
"""Split a group's labelled cell assignments back into its libraries.

Each group barcode is resolved to its library and raw barcode through
barcode_provenance.tsv, never by splitting the synthetic barcode. Every cell of
every library is written, including doublets and unassigned cells, and a library
may hold zero cells of a profile.

Writes, per library:
  <outdir>/<sample_id>/cell_assignments.tsv.gz
  <outdir>/<sample_id>/library_summary.tsv
"""
import argparse
import csv
import gzip
import os
import sys

LABELS = ["fetal_high_confidence_male", "fetal_marker_supported", "maternal_likely",
          "unresolved", "doublet", "unassigned"]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--group-id", required=True)
    ap.add_argument("--cell-assignments", required=True, help="group cell_assignments.tsv.gz")
    ap.add_argument("--provenance", required=True, help="barcode_provenance.tsv")
    ap.add_argument("--gate", required=True, help="gate.tsv")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args(argv)

    with open(args.gate) as fh:
        verdict = next(csv.DictReader(fh, delimiter="\t"))["verdict"]

    provenance = {}
    order = []
    with open(args.provenance) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            provenance[row["group_barcode"]] = (row["sample_id"], row["raw_barcode"])
            if row["sample_id"] not in order:
                order.append(row["sample_id"])

    rows = {sample_id: [] for sample_id in order}
    seen = set()
    with gzip.open(args.cell_assignments, "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            barcode = row["group_barcode"]
            if barcode not in provenance:
                sys.exit("ERROR: %s is not in %s" % (barcode, args.provenance))
            seen.add(barcode)
            sample_id, raw = provenance[barcode]
            rows[sample_id].append((raw, row))
    missing = len(provenance) - len(seen)
    if missing:
        sys.exit("ERROR: %d cells of group %s have no assignment row" % (missing, args.group_id))

    for sample_id in order:
        sample_dir = os.path.join(args.outdir, sample_id)
        os.makedirs(sample_dir, exist_ok=True)
        counts = {"singlet": 0, "doublet": 0, "unassigned": 0}
        labels = {label: 0 for label in LABELS}
        with gzip.open(os.path.join(sample_dir, "cell_assignments.tsv.gz"), "wt",
                       newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["group_id", "sample_id", "barcode", "souporcell_status",
                        "souporcell_cluster", "label", "confidence", "reason"])
            for raw, row in rows[sample_id]:
                counts[row["souporcell_status"]] = counts.get(row["souporcell_status"], 0) + 1
                labels[row["label"]] = labels.get(row["label"], 0) + 1
                w.writerow([args.group_id, sample_id, raw, row["souporcell_status"],
                            row["souporcell_cluster"], row["label"], row["confidence"],
                            row["reason"]])

        with open(os.path.join(sample_dir, "library_summary.tsv"), "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["group_id", "sample_id", "gate_verdict", "n_cells", "n_singlet",
                        "n_doublet", "n_unassigned"] + ["label_" + label for label in LABELS])
            w.writerow([args.group_id, sample_id, verdict, len(rows[sample_id]),
                        counts["singlet"], counts["doublet"], counts["unassigned"]]
                       + [labels[label] for label in LABELS])

    print("[split_by_library] %s -> %d libraries" % (args.group_id, len(order)))


if __name__ == "__main__":
    main()
