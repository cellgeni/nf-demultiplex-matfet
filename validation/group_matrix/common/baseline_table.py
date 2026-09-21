#!/usr/bin/env python3
"""Convert a per-library Souporcell baseline result into the common table.

The acceptance criteria compare group-level results against the per-library
baseline on CW5/Embryo-Imp15350362, the one library that produced a convincing
two-profile result. That baseline predates the harness, so its barcodes are raw
and its profiles are its own cluster ids.

This lives in the shared base because both backends are measured against the
same baseline, and the profile ids must be mapped consistently.

Baseline clusters.tsv columns, as produced by the completed pilot:
  barcode, status, assignment, singlet_posterior, doublet_posterior,
  log_prob_singleton, log_prob_doublet, cluster0..clusterN
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from assignment_table import (  # noqa: E402
    AdapterError,
    profile_name,
    read_locus_qc,
    summarise,
    write_table,
)


def raw_to_group(group_dir: str, sample_id: str) -> dict[str, str]:
    """raw barcode -> group barcode, for one library, via the mapping table."""
    path = os.path.join(group_dir, "barcode_provenance.tsv")
    mapping: dict[str, str] = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["sample_id"] == sample_id:
                mapping[row["raw_barcode"]] = row["group_barcode"]
    if not mapping:
        raise AdapterError(f"{path}: no rows for sample_id {sample_id!r}")
    return mapping


def convert(clusters_tsv: str, group_dir: str, sample_id: str, outdir: str,
            profile_map_path: str | None) -> dict:
    os.makedirs(outdir, exist_ok=True)
    mapping = raw_to_group(group_dir, sample_id)

    loci: dict[str, int] = {}
    cached = os.path.join(group_dir, "loci_per_cell.tsv")
    if os.path.isfile(cached):
        loci = read_locus_qc(cached)

    profile_map: dict[str, str] | None = None
    if profile_map_path:
        with open(profile_map_path) as handle:
            profile_map = json.load(handle).get("cluster_to_profile")

    with open(clusters_tsv, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"barcode", "status", "assignment", "singlet_posterior", "doublet_posterior"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise AdapterError(f"{clusters_tsv}: missing columns {sorted(missing)}")
        records = list(reader)

    if profile_map is None:
        clusters = sorted(
            {r["assignment"] for r in records
             if r["status"] == "singlet" and r["assignment"].isdigit()},
            key=int,
        )
        profile_map = {cluster: profile_name(int(cluster)) for cluster in clusters}

    rows = []
    unmapped = 0
    for record in records:
        group_barcode = mapping.get(record["barcode"])
        if group_barcode is None:
            unmapped += 1
            continue
        status = record["status"]
        if status == "singlet":
            profile = profile_map.get(record["assignment"])
            if profile is None:
                raise AdapterError(
                    f"{clusters_tsv}: cluster {record['assignment']!r} absent from the profile map"
                )
        else:
            profile = ""
            status = "doublet" if status == "doublet" else "unassigned"

        rows.append({
            "group_barcode": group_barcode,
            "sample_id": sample_id,
            "status": status,
            "profile": profile,
            "assign_prob": min(1.0, max(0.0, float(record["singlet_posterior"]))),
            "doublet_prob": min(1.0, max(0.0, float(record["doublet_posterior"]))),
            "n_informative_loci": loci.get(group_barcode, 0),
            "backend": "souporcell",
            "mode": "per_library_baseline",
            "native_status": record["status"],
            "native_profile": record["assignment"] if record["status"] == "singlet" else "",
        })

    if not rows:
        raise AdapterError(
            f"{clusters_tsv}: no barcodes matched library {sample_id!r} in {group_dir}"
        )

    table_path = os.path.join(outdir, "assignments.tsv")
    write_table(table_path, rows)
    manifest = {
        "source_clusters_tsv": os.path.abspath(clusters_tsv),
        "group_dir": os.path.abspath(group_dir),
        "sample_id": sample_id,
        "mode": "per_library_baseline",
        "n_cells": len(rows),
        "n_unmapped_barcodes": unmapped,
        "cluster_to_profile": profile_map,
        "note": "profile ids are this baseline's own clusters; comparisons are "
                "permutation-invariant so the labels need not match another run's",
        "summary": summarise(rows),
    }
    with open(os.path.join(outdir, "adapter_manifest.json"), "w") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--clusters", required=True, help="per-library clusters.tsv")
    parser.add_argument("--group-dir", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--profile-map", default=None,
                        help="optional profile_map.json to reuse cluster->profile names")
    args = parser.parse_args(argv)

    try:
        manifest = convert(args.clusters, args.group_dir, args.sample_id, args.outdir,
                           args.profile_map)
    except (AdapterError, OSError) as exc:
        print(f"BASELINE CONVERSION FAILURE: {exc}", file=sys.stderr)
        return 2
    summary = manifest["summary"]["group"]
    print(f"baseline table: {manifest['n_cells']} cells "
          f"({summary['n_singlet']} singlet, {summary['n_doublet']} doublet, "
          f"{summary['n_unassigned']} unassigned), {manifest['n_unmapped_barcodes']} unmapped")
    print(f"  profiles: {summary['profiles']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
