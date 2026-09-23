#!/usr/bin/env python3
"""Design the fetal-fraction detection-limit experiment.

Question: at what fetal fraction does a real, related second genotype become
recoverable, and where does the composition gate switch from reject to accept?
CW4 and D6 cannot answer this because their composition is unknown. CW5 can,
because it holds a fetal population that both backends agree on.

Construction, per synthetic group:

- **Truth labels are consensus labels**, not ground truth: a cell is fetal or
  maternal only when both backends' joint CW5 runs call it the same way as a
  high-confidence singlet. Cells either backend calls doublet, unassigned or
  low-confidence are excluded.
- **Background**: a fixed number of consensus-maternal cells drawn from *all*
  CW5 libraries. The fetal cells come almost entirely from Embryo-Imp15350362,
  so drawing the background from every library, 362 included, stops library
  identity from predicting genotype. Otherwise a method could separate "362"
  from "not 362" and appear to detect the fetus.
- **Spike**: consensus-fetal cells added so they make up the target fraction of
  the group. Mother and fetus share half their alleles, which is the realistic
  and hard case.
- **Halves**: each group's cells split into two disjoint halves, alternating
  over sorted barcodes, for the split-half composition test.

Deterministic: the same seed always yields the same groups.

Outputs under <outdir>/groups/<label>/: cells.tsv, cells_A.tsv, cells_B.tsv,
truth.tsv; plus <outdir>/spikein_design.tsv listing every group.
"""

from __future__ import annotations

import argparse
import json
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from assignment_table import is_high_confidence, read_table  # noqa: E402


def minority_profile(rows: list[dict]) -> tuple[str, str]:
    counts: dict[str, int] = {}
    for row in rows:
        if is_high_confidence(row):
            counts[row["profile"]] = counts.get(row["profile"], 0) + 1
    if len(counts) != 2:
        raise SystemExit(f"expected two profiles in a truth run, found {counts}")
    minority = min(counts, key=counts.get)
    majority = max(counts, key=counts.get)
    return minority, majority


def consensus_labels(table_a: str, table_b: str) -> tuple[dict[str, str], dict]:
    """barcode -> 'fetal' | 'maternal' where both runs agree confidently."""
    rows_a = read_table(table_a)
    rows_b = read_table(table_b)
    minority_a, majority_a = minority_profile(rows_a)
    minority_b, majority_b = minority_profile(rows_b)

    def label(row, minority, majority):
        if not is_high_confidence(row):
            return None
        if row["profile"] == minority:
            return "fetal"
        if row["profile"] == majority:
            return "maternal"
        return None

    by_b = {row["group_barcode"]: row for row in rows_b}
    labels: dict[str, str] = {}
    disagree = 0
    for row in rows_a:
        other = by_b.get(row["group_barcode"])
        if other is None:
            continue
        la = label(row, minority_a, majority_a)
        lb = label(other, minority_b, majority_b)
        if la is None or lb is None:
            continue
        if la != lb:
            disagree += 1
            continue
        labels[row["group_barcode"]] = la

    sample_of = {row["group_barcode"]: row["sample_id"] for row in rows_a}
    per_library: dict[str, dict[str, int]] = {}
    for barcode, value in labels.items():
        bucket = per_library.setdefault(sample_of[barcode], {"fetal": 0, "maternal": 0})
        bucket[value] += 1
    stats = {
        "n_fetal": sum(1 for v in labels.values() if v == "fetal"),
        "n_maternal": sum(1 for v in labels.values() if v == "maternal"),
        "n_confident_disagreements": disagree,
        "per_library": per_library,
        "sample_of": sample_of,
    }
    return labels, stats


def fraction_label(fraction: float) -> str:
    return "f" + f"{fraction * 100:g}".replace(".", "p")


def build(table_a: str, table_b: str, outdir: str, fractions: list[float],
          replicates: int, background: int, seed: int) -> list[dict]:
    labels, stats = consensus_labels(table_a, table_b)
    fetal = sorted(b for b, v in labels.items() if v == "fetal")
    maternal = sorted(b for b, v in labels.items() if v == "maternal")
    if background > len(maternal):
        raise SystemExit(f"background {background} exceeds {len(maternal)} maternal cells")

    groups_dir = os.path.join(outdir, "groups")
    os.makedirs(groups_dir, exist_ok=True)
    design = []

    for fraction in fractions:
        n_fetal = round(fraction * background / (1 - fraction)) if fraction > 0 else 0
        capped = n_fetal > len(fetal)
        n_fetal = min(n_fetal, len(fetal))
        for rep in range(1, replicates + 1):
            rng = random.Random(f"{seed}:{fraction}:{rep}")
            chosen_bg = rng.sample(maternal, background)
            chosen_fetal = rng.sample(fetal, n_fetal)
            cells = sorted(chosen_bg + chosen_fetal)

            label = f"{fraction_label(fraction)}_r{rep}"
            gdir = os.path.join(groups_dir, label)
            os.makedirs(gdir, exist_ok=True)
            with open(os.path.join(gdir, "cells.tsv"), "w") as out:
                out.write("\n".join(cells) + "\n")
            for name, half in (("A", cells[0::2]), ("B", cells[1::2])):
                with open(os.path.join(gdir, f"cells_{name}.tsv"), "w") as out:
                    out.write("\n".join(half) + "\n")
            fetal_set = set(chosen_fetal)
            with open(os.path.join(gdir, "truth.tsv"), "w") as out:
                out.write("group_barcode\ttruth\tsample_id\n")
                for barcode in cells:
                    truth = "fetal" if barcode in fetal_set else "maternal"
                    out.write(f"{barcode}\t{truth}\t{stats['sample_of'][barcode]}\n")

            per_library: dict[str, dict[str, int]] = {}
            for barcode in cells:
                sample = stats["sample_of"][barcode]
                bucket = per_library.setdefault(sample, {"fetal": 0, "maternal": 0})
                bucket["fetal" if barcode in fetal_set else "maternal"] += 1
            entry = {
                "label": label,
                "fraction_target": fraction,
                "replicate": rep,
                "n_cells": len(cells),
                "n_fetal": n_fetal,
                "n_maternal": background,
                "fraction_actual": n_fetal / len(cells),
                "fetal_pool_capped": capped,
                "per_library": per_library,
            }
            with open(os.path.join(gdir, "design.json"), "w") as out:
                json.dump(entry, out, indent=2, sort_keys=True)
                out.write("\n")
            design.append(entry)

    with open(os.path.join(outdir, "spikein_design.tsv"), "w") as out:
        out.write("label\tfraction_target\treplicate\tn_cells\tn_fetal\tfraction_actual\tfetal_pool_capped\n")
        for d in design:
            out.write(f"{d['label']}\t{d['fraction_target']}\t{d['replicate']}\t{d['n_cells']}\t"
                      f"{d['n_fetal']}\t{d['fraction_actual']:.4f}\t{d['fetal_pool_capped']}\n")
    with open(os.path.join(outdir, "truth_pool.json"), "w") as out:
        json.dump({k: v for k, v in stats.items() if k != "sample_of"}, out, indent=2, sort_keys=True)
        out.write("\n")
    return design


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--truth-a", required=True, help="joint CW5 assignments, backend A")
    parser.add_argument("--truth-b", required=True, help="joint CW5 assignments, backend B")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--fractions", default="0,0.005,0.01,0.02,0.05,0.10")
    parser.add_argument("--replicates", type=int, default=2)
    parser.add_argument("--background", type=int, default=7000)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args(argv)

    design = build(args.truth_a, args.truth_b, args.outdir,
                   [float(x) for x in args.fractions.split(",")],
                   args.replicates, args.background, args.seed)
    pool = json.load(open(os.path.join(args.outdir, "truth_pool.json")))
    print(f"consensus truth: {pool['n_fetal']} fetal, {pool['n_maternal']} maternal, "
          f"{pool['n_confident_disagreements']} confident disagreements excluded")
    print(f"  per library: {pool['per_library']}")
    for d in design:
        cap = "  (fetal pool exhausted: replicates overlap heavily)" if d["fetal_pool_capped"] else ""
        print(f"{d['label']:10} cells={d['n_cells']:5} fetal={d['n_fetal']:4} "
              f"actual={d['fraction_actual']:.4f}{cap}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
