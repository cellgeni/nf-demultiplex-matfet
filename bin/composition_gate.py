#!/usr/bin/env python3
"""Composition gate: is the second genotype Souporcell found in a group real?

Souporcell at K=2 always returns two clusters, including when the group holds a
single genotype; it then halves that genotype arbitrarily, and the two consensus
profiles are near-copies of each other. The gate measures how different the two
profiles are: the discordant fraction, the share of loci called in both profiles
whose genotype calls differ. Arbitrary halves of one genotype differ at 0.17-0.30
of loci in validation (0.27-0.36 on the first production groups); a real
mother/fetus pair at 0.53-0.66 (0.44-0.72 in production).

Verdicts (the evidence and the reasoning behind each threshold are in
docs/pipeline_design.md, "Composition gate"):

  unresolved     the profiles share no called loci (Souporcell collapsed to one
                 genotype), or their discordant fraction is below the threshold
                 (near-copies of one genotype)
  one_profile    the profiles differ, but the minority cluster has too few
                 singlets to be a usable genotype profile
  two_profiles   the profiles differ and the minority cluster is large enough;
                 maternal/fetal orientation then runs

consensus.py's fitted ambient RNA, the doublet fraction and the missing-call
fraction are reported as QC, not used. Ambient was the first gate; it did not
transfer to 5' and 3' v2/v3 data (real fetuses read 34-72%).
"""
import argparse
import csv
import re
import sys

DEFAULTS = {
    "min_discordant_fraction": 0.40,
    "min_minority_singlets": 20,
}

VERDICT_TWO = "two_profiles"
VERDICT_ONE = "one_profile"
VERDICT_UNRESOLVED = "unresolved"

GENOTYPE_CODES = {"0/0": 0, "0|0": 0, "0/1": 1, "1/0": 1, "0|1": 1, "1|0": 1,
                  "1/1": 2, "1|1": 2}


def decide_gate(separating_loci, loci_called_in_both, minority_singlets, params=None):
    """Return (verdict, reason). Pure, so the tests can call it directly."""
    p = dict(DEFAULTS)
    if params:
        p.update({k: v for k, v in params.items() if v is not None})

    if loci_called_in_both == 0:
        return VERDICT_UNRESOLVED, (
            "The two profiles share no called loci: Souporcell collapsed the group to one "
            "genotype. No second genotype is reported.")
    fraction = separating_loci / loci_called_in_both
    if fraction < p["min_discordant_fraction"]:
        return VERDICT_UNRESOLVED, (
            "The profiles differ at %.2f of the %d loci called in both (< %.2f): they are "
            "near-copies of one genotype, the fake-split signature. No second genotype is "
            "reported." % (fraction, loci_called_in_both, p["min_discordant_fraction"]))
    if minority_singlets < p["min_minority_singlets"]:
        return VERDICT_ONE, (
            "The profiles differ at %.2f of called loci, but the minority cluster has %d "
            "singlets (< %d), too few for a usable second profile. One profile is reported."
            % (fraction, minority_singlets, p["min_minority_singlets"]))
    return VERDICT_TWO, (
        "The profiles differ at %.2f of the %d loci called in both (>= %.2f) and the "
        "minority cluster has %d singlets (>= %d): two genotype profiles accepted."
        % (fraction, loci_called_in_both, p["min_discordant_fraction"], minority_singlets,
           p["min_minority_singlets"]))


def parse_ambient(path):
    """'ambient RNA estimated as 1.80%' -> 1.80; None if absent or unparseable."""
    try:
        with open(path) as fh:
            text = fh.read()
    except OSError:
        return None
    match = re.search(r"([0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\s*%", text)
    return float(match.group(1)) if match else None


def read_provenance(path):
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        mapping = {row["group_barcode"]: row["sample_id"] for row in reader}
    if not mapping:
        sys.exit("ERROR: %s has no provenance rows" % path)
    return mapping


def read_clusters(path, provenance):
    """Troublet clusters.tsv -> list of (sample_id, status, assignment)."""
    cells = []
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sample_id = provenance.get(row["barcode"])
            if sample_id is None:
                sys.exit("ERROR: barcode %r in %s is not in barcode_provenance.tsv"
                         % (row["barcode"], path))
            cells.append((sample_id, row["status"], row["assignment"]))
    if len(cells) != len(provenance):
        sys.exit("ERROR: %s has %d cells but the group matrix has %d"
                 % (path, len(cells), len(provenance)))
    return cells


def profile_discordance(path):
    """Loci where both profiles have a call and the calls differ (QC only)."""
    samples = []
    called_both = 0
    discordant = 0
    missing = 0
    records = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = fields[9:]
                continue
            records += 1
            keys = fields[8].split(":")
            calls = []
            for sample_field in fields[9:]:
                gt = dict(zip(keys, sample_field.split(":"))).get("GT", "./.")
                code = GENOTYPE_CODES.get(gt)
                if code is None:
                    missing += 1
                calls.append(code)
            if len(calls) >= 2 and calls[0] is not None and calls[1] is not None:
                called_both += 1
                if calls[0] != calls[1]:
                    discordant += 1
    n_calls = records * len(samples)
    return {
        "n_profiles": len(samples),
        "separating_loci": discordant,
        "loci_called_in_both": called_both,
        "missing_gt_fraction": missing / n_calls if n_calls else 0.0,
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--group-id", required=True)
    ap.add_argument("--clusters", required=True, help="troublet clusters.tsv")
    ap.add_argument("--ambient", required=True, help="consensus.py ambient_rna.txt (QC only)")
    ap.add_argument("--profiles", required=True, help="consensus.py cluster_genotypes.vcf")
    ap.add_argument("--provenance", required=True, help="barcode_provenance.tsv")
    ap.add_argument("--min-discordant-fraction", type=float, default=None)
    ap.add_argument("--min-minority-singlets", type=int, default=None)
    ap.add_argument("--out", required=True)
    ap.add_argument("--libraries-out", required=True)
    args = ap.parse_args(argv)

    params = {
        "min_discordant_fraction": args.min_discordant_fraction,
        "min_minority_singlets": args.min_minority_singlets,
    }
    effective = dict(DEFAULTS)
    effective.update({k: v for k, v in params.items() if v is not None})

    provenance = read_provenance(args.provenance)
    cells = read_clusters(args.clusters, provenance)
    ambient = parse_ambient(args.ambient)
    discordance = profile_discordance(args.profiles)

    singlets = {}
    status_counts = {"singlet": 0, "doublet": 0, "unassigned": 0}
    for _sample, status, assignment in cells:
        status_counts[status] = status_counts.get(status, 0) + 1
        if status == "singlet":
            singlets[assignment] = singlets.get(assignment, 0) + 1
    clusters = sorted(singlets, key=lambda c: (len(c), c))
    ordered = sorted(clusters, key=lambda c: singlets[c], reverse=True)
    minority_cluster = ordered[1] if len(ordered) >= 2 else "NA"
    minority_singlets = singlets[ordered[1]] if len(ordered) >= 2 else 0

    separating = discordance["separating_loci"]
    called = discordance["loci_called_in_both"]
    verdict, reason = decide_gate(separating, called, minority_singlets, params)
    n_cells = len(cells)

    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow([
            "group_id", "verdict", "reason", "discordant_fraction", "separating_loci",
            "loci_called_in_both", "n_cells", "n_singlet", "n_doublet", "n_unassigned",
            "doublet_fraction", "cluster_singlets", "minority_cluster", "minority_singlets",
            "ambient_rna_pct", "missing_gt_fraction", "min_discordant_fraction",
            "min_minority_singlets",
        ])
        w.writerow([
            args.group_id, verdict, reason,
            "%.4f" % (separating / called) if called else "NA", separating, called,
            n_cells, status_counts.get("singlet", 0), status_counts.get("doublet", 0),
            status_counts.get("unassigned", 0),
            "%.4f" % (status_counts.get("doublet", 0) / n_cells if n_cells else 0.0),
            ";".join("%s:%d" % (c, singlets[c]) for c in clusters),
            minority_cluster, minority_singlets,
            "NA" if ambient is None else "%.3f" % ambient,
            "%.4f" % discordance["missing_gt_fraction"],
            effective["min_discordant_fraction"], effective["min_minority_singlets"],
        ])

    per_library = {}
    order = []
    for sample_id, status, assignment in cells:
        if sample_id not in per_library:
            order.append(sample_id)
            per_library[sample_id] = {"n_cells": 0, "singlet": 0, "doublet": 0,
                                      "unassigned": 0, "clusters": {}}
        bucket = per_library[sample_id]
        bucket["n_cells"] += 1
        bucket[status] = bucket.get(status, 0) + 1
        if status == "singlet":
            bucket["clusters"][assignment] = bucket["clusters"].get(assignment, 0) + 1

    with open(args.libraries_out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["group_id", "sample_id", "n_cells", "n_singlet", "n_doublet",
                    "n_unassigned"] + ["singlets_cluster_%s" % c for c in clusters])
        for sample_id in order:
            b = per_library[sample_id]
            w.writerow([args.group_id, sample_id, b["n_cells"], b["singlet"], b["doublet"],
                        b["unassigned"]] + [b["clusters"].get(c, 0) for c in clusters])

    print("[composition_gate] %s -> %s (discordant fraction %s, minority %d singlets)"
          % (args.group_id, verdict, "%.2f" % (separating / called) if called else "NA",
             minority_singlets))


if __name__ == "__main__":
    main()
