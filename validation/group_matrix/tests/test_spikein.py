#!/usr/bin/env python3
"""Self-test for the spike-in scorer. Synthetic, stdlib only."""

from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "common"))

import score_spikein as ss  # noqa: E402

failures = []


def check(name, cond, detail=""):
    print(("PASS  " if cond else "FAIL  ") + name + ("" if cond else f": {detail}"))
    if not cond:
        failures.append(name)


def row(bc, status, profile, prob=0.99):
    return {"group_barcode": bc, "sample_id": "lib", "status": status,
            "profile": profile, "assign_prob": prob, "doublet_prob": 0.0}


# 10 spiked fetal cells, 90 maternal. Profile names deliberately swapped
# relative to size so the scorer must pick the fetal profile by overlap.
truth = {f"f{i}": "fetal" for i in range(10)} | {f"m{i}": "maternal" for i in range(90)}
rows = [row(f"f{i}", "singlet", "profile_1") for i in range(8)]     # 8 recovered
rows += [row("f8", "doublet", ""), row("f9", "singlet", "profile_2")]  # 1 doublet, 1 missed
rows += [row(f"m{i}", "singlet", "profile_2") for i in range(88)]
rows += [row("m88", "singlet", "profile_1"), row("m89", "unassigned", "")]  # 1 false positive

d = ss.detection(rows, truth)
check("fetal profile found by overlap", d["fetal_profile"] == "profile_1", d)
check("recall 8/10", abs(d["recall"] - 0.8) < 1e-9, d["recall"])
check("precision 8/9", abs(d["precision"] - 8 / 9) < 1e-9, d["precision"])
check("fetal fates tracked",
      d["fetal_cell_fates"] == {"fetal_singlet": 8, "maternal_singlet": 1, "doublet": 1,
                                "unassigned": 0, "low_confidence": 0}, d["fetal_cell_fates"])

# 0% control: the smaller profile is counted as false positives.
truth0 = {f"m{i}": "maternal" for i in range(100)}
rows0 = [row(f"m{i}", "singlet", "profile_1") for i in range(95)]
rows0 += [row(f"m{i}", "singlet", "profile_2") for i in range(95, 100)]
d0 = ss.detection(rows0, truth0)
check("0% control reports false positives", d0["false_positive_cells"] == 5 and d0["recall"] is None, d0)

# A run that collapses to one profile recovers nothing.
rows1 = [row(bc, "singlet", "profile_1") for bc in truth]
d1 = ss.detection(rows1, truth)
check("single-profile run has zero recall", d1["recall"] == 0.0, d1)

print()
if failures:
    print(f"{len(failures)} FAILED"); sys.exit(1)
print("ALL SPIKE-IN SCORER TESTS PASSED")
