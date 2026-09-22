#!/usr/bin/env python3
"""Model-selection evidence for "is there really more than one genotype?".

Two independent lines of evidence, because neither alone is sufficient:

1. **Penalised likelihood over a K sweep.** Raw log-likelihood always improves
   with more donors, so it cannot select K on its own. AIC and BIC add a
   parameter penalty. Both are heuristics here: Vireo fits a variational Bayes
   model, its reported `logLik` is not an exact marginal likelihood, and the
   parameter count below is an approximation. Treat this as supporting evidence,
   never as the verdict.

2. **Split-half profile reproducibility.** Split a group's cells into two
   disjoint halves, learn K=2 independently on each, and compare the inferred
   genotypes. A real second donor is present in both halves, so the profiles
   agree. A split driven by noise, coverage structure or ambient RNA has no
   stable genotype behind it, so the two halves disagree. This is the primary
   evidence, it needs no penalty assumptions, and it is interpreted against a
   positive control (a group known to contain two genotypes).

Parameter count for the penalty, per K:
    2 * K * n_variants   genotype probabilities (3 categories, 2 free) per donor
    + (K - 1)            donor proportions
    + 3                  the shared beta parameters over genotype classes
Observations are the (variant, cell) pairs with any allele coverage, since a
stored Vartrix zero is not evidence.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from build_group_matrix import open_mtx_body  # noqa: E402
from genotype_vcf import between_run_concordance, within_run_discordance  # noqa: E402

LOGLIK_RE = re.compile(r"logLik[:=]\s*(-?[\d.eE+]+)")


def parse_vireo_log(path: str) -> dict:
    """Pull logLik out of Vireo's _log.txt. Returns {} when absent."""
    if not os.path.isfile(path):
        return {}
    with open(path) as handle:
        text = handle.read()
    match = LOGLIK_RE.search(text)
    return {"logLik": float(match.group(1))} if match else {}


def count_observations(group_dir: str) -> dict:
    """(variant, cell) pairs with coverage, and variants with any coverage."""
    alt_in, alt_dims = open_mtx_body(os.path.join(group_dir, "alt.mtx"))
    ref_in, _ = open_mtx_body(os.path.join(group_dir, "ref.mtx"))
    covered_rows: set[int] = set()
    n_obs = 0
    try:
        for alt_line, ref_line in zip(alt_in, ref_in):
            a = alt_line.split()
            r = ref_line.split()
            if float(a[2]) + float(r[2]) > 0:
                n_obs += 1
                covered_rows.add(int(a[0]))
    finally:
        alt_in.close()
        ref_in.close()
    return {"n_observations": n_obs,
            "n_variants_with_coverage": len(covered_rows),
            "n_variants_total": alt_dims[0],
            "n_cells": alt_dims[1]}


def n_parameters(k: int, n_variants: int) -> int:
    return 2 * k * n_variants + (k - 1) + 3


def penalised(log_lik: float, k: int, n_variants: int, n_obs: int) -> dict:
    p = n_parameters(k, n_variants)
    return {
        "n_parameters": p,
        "aic": 2 * p - 2 * log_lik,
        "bic": p * math.log(n_obs) - 2 * log_lik if n_obs > 0 else None,
    }


def k_sweep(group_dir: str, runs: list[tuple[int, str]]) -> dict:
    counts = count_observations(group_dir)
    n_variants = counts["n_variants_with_coverage"]
    n_obs = counts["n_observations"]

    table = []
    for k, outdir in runs:
        entry = {"k": k, "vireo_out": os.path.abspath(outdir)}
        entry.update(parse_vireo_log(os.path.join(outdir, "_log.txt")))
        if "logLik" in entry:
            entry.update(penalised(entry["logLik"], k, n_variants, n_obs))
        else:
            entry["note"] = "no logLik found; did this run complete?"
        vcf = os.path.join(outdir, "GT_donors.vireo.vcf.gz")
        if os.path.isfile(vcf) and k >= 2:
            within = within_run_discordance(vcf)
            entry["min_discordant_loci_over_pairs"] = within.get(
                "min_discordant_loci_over_pairs"
            )
        table.append(entry)

    evaluable = [row for row in table if row.get("bic") is not None]
    best_bic = min(evaluable, key=lambda row: row["bic"])["k"] if evaluable else None
    best_aic = min(evaluable, key=lambda row: row["aic"])["k"] if evaluable else None
    monotone = all(
        evaluable[i]["logLik"] <= evaluable[i + 1]["logLik"] + 1e-6
        for i in range(len(evaluable) - 1)
    ) if len(evaluable) > 1 else None

    return {
        "counts": counts,
        "runs": table,
        "k_selected_by_bic": best_bic,
        "k_selected_by_aic": best_aic,
        "loglik_increases_with_k": monotone,
        "caveat": "variational logLik with an approximate parameter count; "
                  "supporting evidence only",
    }


def split_half(vcf_a: str, vcf_b: str, label: str) -> dict:
    """Profile agreement between two disjoint halves of the same group."""
    result = {"label": label, "vcf_a": os.path.abspath(vcf_a),
              "vcf_b": os.path.abspath(vcf_b)}
    try:
        concordance = between_run_concordance(vcf_a, vcf_b)
    except Exception as exc:  # noqa: BLE001 - record why, do not crash the report
        result["error"] = str(exc)
        return result
    result.update({
        "genotype_concordance": concordance["concordance"],
        "n_loci_compared": concordance["n_loci"],
        "profile_mapping": concordance["mapping"],
    })
    for key, path in (("within_a", vcf_a), ("within_b", vcf_b)):
        within = within_run_discordance(path)
        result[key] = {
            "min_discordant_loci_over_pairs": within.get("min_discordant_loci_over_pairs"),
            "missing_gt_fraction": within.get("missing_gt_fraction"),
        }
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--group-dir", default=None,
                        help="group build the K sweep ran on (for observation counts)")
    parser.add_argument("--run", action="append", default=[], metavar="K:DIR",
                        help="a K-sweep run, e.g. 2:/path/to/vireo_out; repeatable")
    parser.add_argument("--split-half", action="append", default=[], metavar="LABEL:VCF_A:VCF_B",
                        help="a split-half comparison; repeatable")
    parser.add_argument("--out", required=True)
    args = parser.parse_args(argv)

    report: dict = {}
    if args.run:
        if not args.group_dir:
            parser.error("--group-dir is required with --run")
        runs = []
        for item in args.run:
            k, _, path = item.partition(":")
            runs.append((int(k), path))
        report["k_sweep"] = k_sweep(args.group_dir, sorted(runs))

    if args.split_half:
        report["split_half"] = []
        for item in args.split_half:
            label, vcf_a, vcf_b = item.split(":", 2)
            report["split_half"].append(split_half(vcf_a, vcf_b, label))

    with open(args.out, "w") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")

    if "k_sweep" in report:
        sweep = report["k_sweep"]
        print("K sweep:")
        for row in sweep["runs"]:
            if "bic" in row and row["bic"] is not None:
                print(f"  K={row['k']}  logLik={row['logLik']:.4e}  "
                      f"AIC={row['aic']:.4e}  BIC={row['bic']:.4e}  "
                      f"discordant_loci={row.get('min_discordant_loci_over_pairs')}")
            else:
                print(f"  K={row['k']}  {row.get('note', 'incomplete')}")
        print(f"  BIC prefers K={sweep['k_selected_by_bic']}, "
              f"AIC prefers K={sweep['k_selected_by_aic']}")
    for entry in report.get("split_half", []):
        if "error" in entry:
            print(f"split-half {entry['label']}: ERROR {entry['error']}")
            continue
        print(f"split-half {entry['label']}: genotype concordance "
              f"{entry['genotype_concordance']:.4f} over {entry['n_loci_compared']} loci "
              f"(within-half discordant loci: "
              f"{entry['within_a']['min_discordant_loci_over_pairs']} / "
              f"{entry['within_b']['min_discordant_loci_over_pairs']})")
    print(f"report: {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
