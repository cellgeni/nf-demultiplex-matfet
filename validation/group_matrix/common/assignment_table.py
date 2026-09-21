#!/usr/bin/env python3
"""The common assignment table that every backend adapter must emit.

Defined once here so Souporcell and Vireo results are comparable on identical
terms (docs/group_aware_redesign.md, "Provisional validation acceptance
criteria"). Backends differ in what they report natively, so native fields are
preserved in extra columns while the common columns carry a single meaning.

Columns:
  group_barcode         <sample_id>::<raw_barcode>, as built by the harness
  sample_id             resolved through barcode_provenance.tsv, never by
                        splitting the synthetic barcode
  status                singlet | doublet | unassigned
  profile               profile_1..profile_K for a singlet, else empty.
                        Neutral names: orientation to maternal/fetal happens
                        later, from a reference or marker evidence.
  assign_prob           normalized assignment probability in [0, 1]
  doublet_prob          doublet probability in [0, 1]
  n_informative_loci    loci with any coverage for this cell, counted from the
                        group matrices so the definition is backend-independent
  backend               vireo | souporcell
  mode                  discovery | projection
  native_status         the backend's own label, unmodified
  native_profile        the backend's own cluster/donor id, unmodified
"""

from __future__ import annotations

import csv
import os

COLUMNS = [
    "group_barcode",
    "sample_id",
    "status",
    "profile",
    "assign_prob",
    "doublet_prob",
    "n_informative_loci",
    "backend",
    "mode",
    "native_status",
    "native_profile",
]

STATUSES = {"singlet", "doublet", "unassigned"}

# A high-confidence cell, per the acceptance criteria: a singlet whose
# normalized assignment probability is at least 0.90.
HIGH_CONFIDENCE_MIN_PROB = 0.90


class AdapterError(Exception):
    """A backend adapter could not produce a contract-conforming result."""


def profile_name(index: int) -> str:
    """0-based backend cluster/donor index -> neutral profile name."""
    return f"profile_{index + 1}"


def read_provenance(group_dir: str) -> dict[str, str]:
    """group_barcode -> sample_id, from the harness mapping table."""
    path = os.path.join(group_dir, "barcode_provenance.tsv")
    mapping: dict[str, str] = {}
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            mapping[row["group_barcode"]] = row["sample_id"]
    if not mapping:
        raise AdapterError(f"{path}: no barcode provenance rows")
    return mapping


def count_informative_loci(group_dir: str) -> dict[str, int]:
    """Loci with any allele coverage per cell, from the group matrices.

    Streams ALT and REF together and counts, per column, the entries where
    alt + ref > 0. Vartrix writes explicit zeros, so a stored entry is not by
    itself evidence of coverage.
    """
    import sys

    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from build_group_matrix import open_mtx_body  # noqa: E402

    with open(os.path.join(group_dir, "barcodes.tsv")) as handle:
        barcodes = [line.strip() for line in handle if line.strip()]

    counts = [0] * len(barcodes)
    alt_in, alt_dims = open_mtx_body(os.path.join(group_dir, "alt.mtx"))
    ref_in, ref_dims = open_mtx_body(os.path.join(group_dir, "ref.mtx"))
    try:
        if alt_dims != ref_dims:
            raise AdapterError(f"{group_dir}: ALT/REF dimensions differ")
        if alt_dims[1] != len(barcodes):
            raise AdapterError(f"{group_dir}: matrix columns do not match barcodes.tsv")
        for alt_line, ref_line in zip(alt_in, ref_in):
            alt_fields = alt_line.split()
            ref_fields = ref_line.split()
            if alt_fields[1] != ref_fields[1]:
                raise AdapterError(f"{group_dir}: ALT/REF coordinate mismatch while counting loci")
            if float(alt_fields[2]) + float(ref_fields[2]) > 0:
                counts[int(alt_fields[1]) - 1] += 1
    finally:
        alt_in.close()
        ref_in.close()

    return dict(zip(barcodes, counts))


def write_locus_qc(group_dir: str, path: str) -> dict[str, int]:
    counts = count_informative_loci(group_dir)
    with open(path, "w") as out:
        out.write("group_barcode\tn_informative_loci\n")
        for barcode, count in counts.items():
            out.write(f"{barcode}\t{count}\n")
    return counts


def read_locus_qc(path: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    with open(path) as handle:
        handle.readline()
        for line in handle:
            barcode, count = line.rstrip("\n").split("\t")
            counts[barcode] = int(count)
    return counts


def write_table(path: str, rows: list[dict]) -> None:
    """Write the common table, validating each row before it is committed."""
    for row in rows:
        missing = [column for column in COLUMNS if column not in row]
        if missing:
            raise AdapterError(f"assignment row missing columns {missing}: {row}")
        if row["status"] not in STATUSES:
            raise AdapterError(f"unknown status {row['status']!r}")
        if row["status"] == "singlet" and not row["profile"]:
            raise AdapterError(f"singlet without a profile: {row['group_barcode']}")
        if row["status"] != "singlet" and row["profile"]:
            raise AdapterError(
                f"non-singlet {row['group_barcode']} carries profile {row['profile']!r}"
            )
        for field in ("assign_prob", "doublet_prob"):
            value = float(row[field])
            if not 0.0 <= value <= 1.0:
                raise AdapterError(f"{field}={value} outside [0, 1] for {row['group_barcode']}")

    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t",
                                extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_table(path: str) -> list[dict]:
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = [column for column in COLUMNS if column not in (reader.fieldnames or [])]
        if missing:
            raise AdapterError(f"{path}: missing columns {missing}")
        rows = list(reader)
    for row in rows:
        row["assign_prob"] = float(row["assign_prob"])
        row["doublet_prob"] = float(row["doublet_prob"])
        row["n_informative_loci"] = int(row["n_informative_loci"])
    return rows


def is_high_confidence(row: dict) -> bool:
    return row["status"] == "singlet" and float(row["assign_prob"]) >= HIGH_CONFIDENCE_MIN_PROB


def summarise(rows: list[dict]) -> dict:
    """Per-group and per-library counts the validation report needs."""
    def blank() -> dict:
        return {
            "n_cells": 0,
            "n_singlet": 0,
            "n_doublet": 0,
            "n_unassigned": 0,
            "n_high_confidence": 0,
            "profiles": {},
            "loci_sum": 0,
        }

    overall = blank()
    per_library: dict[str, dict] = {}

    for row in rows:
        for bucket in (overall, per_library.setdefault(row["sample_id"], blank())):
            bucket["n_cells"] += 1
            bucket[f"n_{row['status']}"] += 1
            bucket["loci_sum"] += int(row["n_informative_loci"])
            if is_high_confidence(row):
                bucket["n_high_confidence"] += 1
                bucket["profiles"][row["profile"]] = bucket["profiles"].get(row["profile"], 0) + 1

    for bucket in [overall, *per_library.values()]:
        n = bucket["n_cells"] or 1
        bucket["doublet_fraction"] = bucket["n_doublet"] / n
        bucket["unassigned_fraction"] = bucket["n_unassigned"] / n
        bucket["median_informative_loci"] = None
        bucket["mean_informative_loci"] = bucket["loci_sum"] / n

    loci = sorted(int(row["n_informative_loci"]) for row in rows)
    if loci:
        overall["median_informative_loci"] = loci[len(loci) // 2]
    for sample_id, bucket in per_library.items():
        sample_loci = sorted(
            int(row["n_informative_loci"]) for row in rows if row["sample_id"] == sample_id
        )
        if sample_loci:
            bucket["median_informative_loci"] = sample_loci[len(sample_loci) // 2]

    return {"group": overall, "libraries": per_library}
