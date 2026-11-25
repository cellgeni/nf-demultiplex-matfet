#!/usr/bin/env python3
"""Unit tests for the fetal-mode summary (bin/fetal_summary.py).

The bound and the XIST/Y numbers are the embryo results the method reported (F147, F158,
F218); the warnings are the in-range guards from the D8 benchmark.

Runnable two ways:
    pytest tests/test_fetal_summary.py
    python3 tests/test_fetal_summary.py      # built-in mini-runner, no pytest needed
"""
import math
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

import fetal_summary as fs  # noqa: E402


def cells(**counts):
    return pd.DataFrame({"fetal_mode": [m for m, k in counts.items() for _ in range(k)]})


# --- the pooled zero-event bound, as reported per embryo -----------------------------------
def test_bound_matches_reported_embryo_bounds():
    # cells cleared as fetal in F147, F158, F218; 3 / cleared = 8.27e-6, 9.91e-6, 7.66e-6
    for cleared, bound in ((362678, 8.27e-6), (302581, 9.91e-6), (391773, 7.66e-6)):
        assert abs(fs.exact_upper95(cleared) - bound) < 0.01e-6


def test_bound_without_cleared_cells_is_undefined():
    assert math.isnan(fs.exact_upper95(0))


def test_tally_counts_pairs_as_flagged_and_bounds_on_fetal_only():
    t = fs.tally(cells(fetal=1000, too_little_data=500, maternal=5, maternal_in_pair=3))
    assert t["flagged"] == 8 and t["cells"] == 1508
    assert t["upper95_maternal_among_fetal"] == 3.0 / 1000


def test_every_status_maps_to_one_atlas_value():
    assert set(fs.FETAL_MODE.values()) == {"fetal", "too_little_data", "maternal", "maternal_in_pair"}


# --- sex call and the XIST/Y check ------------------------------------------------------------
def sex(n_y, n_noy, xist_noy=0, total=8000):
    rows = [(total, 0, 1)] * n_y + [(total, 3, 0)] * xist_noy + [(total, 0, 0)] * (n_noy - xist_noy)
    return pd.DataFrame(rows, columns=["total_umis", "xist", "y_umis"])


def test_male_with_no_candidates_is_in_range():
    row, warnings = fs.sex_check(sex(9400, 600), 0.02)
    assert row["sex_call"] == "male" and row["xist_candidates"] == 0 and not warnings
    assert row["xist_upper95_maternal_share"] == 3 / (0.321 * 10000)


def test_male_with_many_xist_only_cells_is_warned():
    row, warnings = fs.sex_check(sex(9000, 1000, xist_noy=300), 0.02)
    assert row["xist_upper95_maternal_share"] > 0.02 and len(warnings) == 1


def test_female_has_no_genotype_free_check():
    row, warnings = fs.sex_check(sex(10, 9990), 0.02)
    assert row["sex_call"] == "female" and not warnings and "no genotype-free check" in row["note"]


def test_mixed_sex_signal_is_warned():
    row, warnings = fs.sex_check(sex(3000, 7000), 0.02)
    assert row["sex_call"] == "unclear" and "mix individuals" in warnings[0]


def test_xist_detection_follows_depth_bins():
    assert fs.xist_detection(pd.Series([1000, 30000]).values) == (0.070 + 0.871) / 2


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    failed = 0
    for t in tests:
        try:
            t()
            print(f"PASS {t.__name__}")
        except AssertionError as e:
            failed += 1
            print(f"FAIL {t.__name__}: {e}")
    print(f"{len(tests) - failed}/{len(tests)} passed")
    sys.exit(1 if failed else 0)
