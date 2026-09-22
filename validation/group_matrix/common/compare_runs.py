#!/usr/bin/env python3
"""Compare two backend runs against the provisional acceptance criteria.

Used for every pairwise comparison the gates require: seed against seed (Gate 2
stability), projection against joint inference (Gate 3), backend against
backend or against the per-library baseline (Gate 2), and discovery subset
against discovery subset (Gate 4).

The comparison is permutation-invariant: profile names are arbitrary labels
until a maternal reference or marker evidence orients them, so the best
label matching is chosen before any agreement is computed.

Reports numbers and each criterion's verdict. It does not adjust a threshold to
make a run pass; thresholds live in THRESHOLDS and any change must be justified
in the validation report.
"""

from __future__ import annotations

import argparse
import itertools
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from assignment_table import is_high_confidence, read_table, summarise  # noqa: E402
from genotype_vcf import between_run_concordance, within_run_discordance  # noqa: E402

THRESHOLDS = {
    "cell_concordance_min": 0.95,      # Gate 2, against the strong baseline
    "projection_concordance_min": 0.90,  # Gate 3, projection vs joint
    "ari_min": 0.95,                   # Gate 2, across initialisations
    "ari_min_subsets": 0.90,           # Gate 4, across discovery subsets
    "genotype_concordance_min": 0.95,
    "discordant_loci_min": 100,        # per accepted profile pair
    # Provisional. Cell agreement can be 1.0 while a projection discards the
    # minority population, because demoted cells stop being "common
    # high-confidence singlets" and leave the comparison. Retention is measured
    # separately for that reason. Any change must be justified in the report.
    "minority_retention_min": 0.80,
}


def adjusted_rand_index(labels_a: list[str], labels_b: list[str]) -> float:
    """ARI over paired labels. Pure stdlib; the pair counts are exact integers."""
    if len(labels_a) != len(labels_b):
        raise ValueError("label vectors differ in length")
    if not labels_a:
        return float("nan")

    contingency: dict[tuple[str, str], int] = {}
    rows: dict[str, int] = {}
    columns: dict[str, int] = {}
    for a, b in zip(labels_a, labels_b):
        contingency[(a, b)] = contingency.get((a, b), 0) + 1
        rows[a] = rows.get(a, 0) + 1
        columns[b] = columns.get(b, 0) + 1

    def choose2(n: int) -> int:
        return n * (n - 1) // 2

    sum_ij = sum(choose2(n) for n in contingency.values())
    sum_i = sum(choose2(n) for n in rows.values())
    sum_j = sum(choose2(n) for n in columns.values())
    total = choose2(len(labels_a))
    if total == 0:
        return float("nan")
    expected = sum_i * sum_j / total
    maximum = (sum_i + sum_j) / 2
    if maximum == expected:
        return 1.0 if sum_ij == expected else 0.0
    return (sum_ij - expected) / (maximum - expected)


def best_profile_matching(rows_a: list[dict], rows_b: list[dict]) -> dict:
    """Choose the profile label mapping that maximises agreement."""
    by_barcode_a = {row["group_barcode"]: row for row in rows_a}
    by_barcode_b = {row["group_barcode"]: row for row in rows_b}
    common = [
        barcode for barcode in by_barcode_a
        if barcode in by_barcode_b
        and is_high_confidence(by_barcode_a[barcode])
        and is_high_confidence(by_barcode_b[barcode])
    ]

    profiles_a = sorted({by_barcode_a[b]["profile"] for b in common})
    profiles_b = sorted({by_barcode_b[b]["profile"] for b in common})

    if not common:
        return {"mapping": {}, "n_common_high_confidence": 0, "concordance": float("nan"),
                "n_agree": 0, "profiles_a": profiles_a, "profiles_b": profiles_b}

    best = None
    # Pad so the permutation search works when one run found fewer profiles.
    padded_b = profiles_b + [None] * max(0, len(profiles_a) - len(profiles_b))
    for permutation in itertools.permutations(padded_b, len(profiles_a)):
        mapping = dict(zip(profiles_a, permutation))
        agree = sum(
            1 for barcode in common
            if mapping.get(by_barcode_a[barcode]["profile"]) == by_barcode_b[barcode]["profile"]
        )
        candidate = {
            "mapping": {key: value for key, value in mapping.items() if value is not None},
            "n_common_high_confidence": len(common),
            "n_agree": agree,
            "concordance": agree / len(common),
            "profiles_a": profiles_a,
            "profiles_b": profiles_b,
        }
        if best is None or candidate["concordance"] > best["concordance"]:
            best = candidate

    labels_a = [by_barcode_a[barcode]["profile"] for barcode in common]
    labels_b = [by_barcode_b[barcode]["profile"] for barcode in common]
    best["adjusted_rand_index"] = adjusted_rand_index(labels_a, labels_b)
    best["common_barcodes"] = common
    return best


def minority_retention(rows_ref: list[dict], rows_test: list[dict],
                       mapping: dict) -> dict:
    """Does the test run keep the reference run's minority population?

    The minority profile is the biologically interesting one here: it is the
    candidate fetal population. Cell agreement does not measure this, so the
    fate of every reference minority cell is tracked explicitly.
    """
    by_barcode_test = {row["group_barcode"]: row for row in rows_test}

    counts: dict[str, int] = {}
    for row in rows_ref:
        if is_high_confidence(row):
            counts[row["profile"]] = counts.get(row["profile"], 0) + 1
    if len(counts) < 2:
        return {"evaluable": False,
                "note": "the reference run has fewer than two profiles, so there "
                        "is no minority population to retain"}

    minority = min(counts, key=counts.get)
    majority = max(counts, key=counts.get)
    expected = mapping.get(minority)

    fates = {"singlet_minority": 0, "singlet_other_profile": 0, "doublet": 0,
             "unassigned": 0, "low_confidence_singlet": 0, "absent": 0}
    for row in rows_ref:
        if not (is_high_confidence(row) and row["profile"] == minority):
            continue
        other = by_barcode_test.get(row["group_barcode"])
        if other is None:
            fates["absent"] += 1
        elif other["status"] == "doublet":
            fates["doublet"] += 1
        elif other["status"] == "unassigned":
            fates["unassigned"] += 1
        elif not is_high_confidence(other):
            fates["low_confidence_singlet"] += 1
        elif expected is not None and other["profile"] == expected:
            fates["singlet_minority"] += 1
        else:
            fates["singlet_other_profile"] += 1

    n_minority = counts[minority]
    common = [row for row in rows_ref if row["group_barcode"] in by_barcode_test]
    ref_doublets = sum(1 for row in common if row["status"] == "doublet")
    test_doublets = sum(1 for row in common
                        if by_barcode_test[row["group_barcode"]]["status"] == "doublet")

    return {
        "evaluable": True,
        "minority_profile_in_reference": minority,
        "majority_profile_in_reference": majority,
        "n_reference_minority_high_confidence": n_minority,
        "matched_profile_in_test": expected,
        "fates": fates,
        "retention": fates["singlet_minority"] / n_minority if n_minority else 0.0,
        "n_common_cells": len(common),
        "reference_doublet_fraction": ref_doublets / len(common) if common else 0.0,
        "test_doublet_fraction": test_doublets / len(common) if common else 0.0,
        "doublet_inflation": (
            (test_doublets / len(common)) - (ref_doublets / len(common))
            if common else 0.0
        ),
    }


def per_library_agreement(rows_a: list[dict], rows_b: list[dict], mapping: dict) -> dict:
    by_barcode_b = {row["group_barcode"]: row for row in rows_b}
    result: dict[str, dict] = {}
    for row in rows_a:
        other = by_barcode_b.get(row["group_barcode"])
        if other is None or not (is_high_confidence(row) and is_high_confidence(other)):
            continue
        bucket = result.setdefault(row["sample_id"], {"n_common": 0, "n_agree": 0})
        bucket["n_common"] += 1
        if mapping.get(row["profile"]) == other["profile"]:
            bucket["n_agree"] += 1
    for bucket in result.values():
        bucket["concordance"] = bucket["n_agree"] / bucket["n_common"] if bucket["n_common"] else float("nan")
    return result


def compare(table_a: str, table_b: str, vcf_a: str | None, vcf_b: str | None,
            comparison: str) -> dict:
    rows_a = read_table(table_a)
    rows_b = read_table(table_b)

    matching = best_profile_matching(rows_a, rows_b)
    report = {
        "comparison": comparison,
        "inputs": {"table_a": os.path.abspath(table_a), "table_b": os.path.abspath(table_b),
                   "vcf_a": os.path.abspath(vcf_a) if vcf_a else None,
                   "vcf_b": os.path.abspath(vcf_b) if vcf_b else None},
        "thresholds": THRESHOLDS,
        "summary_a": summarise(rows_a),
        "summary_b": summarise(rows_b),
        "profile_matching": {key: value for key, value in matching.items()
                             if key != "common_barcodes"},
        "per_library_agreement": per_library_agreement(rows_a, rows_b, matching["mapping"]),
        # A is the reference (e.g. joint inference), B the run under test.
        "minority_retention": minority_retention(rows_a, rows_b, matching["mapping"]),
    }

    if vcf_a:
        report["within_run_a"] = within_run_discordance(vcf_a)
    if vcf_b:
        report["within_run_b"] = within_run_discordance(vcf_b)
    if vcf_a and vcf_b:
        report["genotype_concordance"] = between_run_concordance(vcf_a, vcf_b)

    # Verdicts. A criterion that cannot be evaluated is reported as such rather
    # than silently passing.
    verdicts: dict[str, dict] = {}

    cell_min = (THRESHOLDS["projection_concordance_min"]
                if comparison.startswith("projection") else THRESHOLDS["cell_concordance_min"])
    verdicts["cell_concordance"] = _verdict(
        matching["concordance"], cell_min,
        f"high-confidence cell agreement >= {cell_min}",
        matching["n_common_high_confidence"] > 0,
    )

    ari_min = (THRESHOLDS["ari_min_subsets"]
               if comparison.startswith("subset") else THRESHOLDS["ari_min"])
    verdicts["adjusted_rand_index"] = _verdict(
        matching.get("adjusted_rand_index"), ari_min,
        f"adjusted Rand index >= {ari_min}",
        matching["n_common_high_confidence"] > 0,
    )

    retention = report["minority_retention"]
    verdicts["minority_retention"] = _verdict(
        retention.get("retention"), THRESHOLDS["minority_retention_min"],
        f"minority-profile retention >= {THRESHOLDS['minority_retention_min']}",
        retention.get("evaluable", False),
    )

    if "genotype_concordance" in report:
        verdicts["genotype_concordance"] = _verdict(
            report["genotype_concordance"]["concordance"],
            THRESHOLDS["genotype_concordance_min"],
            f"genotype concordance >= {THRESHOLDS['genotype_concordance_min']}",
            report["genotype_concordance"]["n_loci"] > 0,
        )

    for label in ("within_run_a", "within_run_b"):
        if label in report and report[label].get("n_profiles", 0) >= 2:
            observed = report[label]["min_discordant_loci_over_pairs"]
            verdicts[f"discordant_loci_{label[-1]}"] = _verdict(
                observed, THRESHOLDS["discordant_loci_min"],
                f"well-supported discordant loci >= {THRESHOLDS['discordant_loci_min']}",
                True,
            )

    report["verdicts"] = verdicts
    report["all_evaluated_criteria_passed"] = all(
        entry["passed"] for entry in verdicts.values() if entry["evaluated"]
    )
    report["unevaluated_criteria"] = [
        name for name, entry in verdicts.items() if not entry["evaluated"]
    ]
    return report


def _verdict(observed, threshold, description: str, evaluated: bool) -> dict:
    if not evaluated or observed is None:
        return {"description": description, "observed": observed, "threshold": threshold,
                "evaluated": False, "passed": None,
                "note": "not evaluable from these inputs"}
    return {"description": description, "observed": observed, "threshold": threshold,
            "evaluated": True, "passed": bool(observed >= threshold)}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--table-a", required=True, help="common assignment table, run A")
    parser.add_argument("--table-b", required=True, help="common assignment table, run B")
    parser.add_argument("--vcf-a", default=None, help="profile genotype VCF, run A")
    parser.add_argument("--vcf-b", default=None, help="profile genotype VCF, run B")
    parser.add_argument("--comparison", default="seed",
                        help="label: seed, projection-vs-joint, subset, backend, baseline")
    parser.add_argument("--out", required=True, help="JSON report to write")
    args = parser.parse_args(argv)

    report = compare(args.table_a, args.table_b, args.vcf_a, args.vcf_b, args.comparison)
    with open(args.out, "w") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")

    matching = report["profile_matching"]
    print(f"comparison: {report['comparison']}")
    print(f"  common high-confidence singlets: {matching['n_common_high_confidence']}")
    print(f"  cell agreement: {matching['concordance']:.4f}")
    print(f"  adjusted Rand index: {matching.get('adjusted_rand_index', float('nan')):.4f}")
    if "genotype_concordance" in report:
        genotype = report["genotype_concordance"]
        print(f"  genotype concordance: {genotype['concordance']:.4f} over {genotype['n_loci']} loci")
    retention = report["minority_retention"]
    if retention.get("evaluable"):
        print(f"  minority retention: {retention['retention']:.4f} of "
              f"{retention['n_reference_minority_high_confidence']} reference minority cells "
              f"(fates: {retention['fates']})")
        print(f"  doublet fraction: {retention['reference_doublet_fraction']:.4f} -> "
              f"{retention['test_doublet_fraction']:.4f}")
    for name, entry in report["verdicts"].items():
        state = "PASS" if entry["passed"] else ("n/a " if not entry["evaluated"] else "FAIL")
        print(f"  [{state}] {entry['description']} (observed {entry['observed']})")
    if report["unevaluated_criteria"]:
        print(f"  not evaluable: {report['unevaluated_criteria']}")
    print(f"report: {args.out}")
    return 0 if report["all_evaluated_criteria_passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
