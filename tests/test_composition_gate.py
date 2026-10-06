#!/usr/bin/env python3
"""Unit tests for the composition gate (bin/composition_gate.py).

The cases are real runs the threshold was set from or checked against
(docs/pipeline_design.md, "Composition gate"): separating loci, loci called in
both profiles and minority singlets as each run produced them. Each must land on
the verdict its known truth supports with the default thresholds.

Runnable two ways:
    pytest tests/test_composition_gate.py
    python3 tests/test_composition_gate.py      # built-in mini-runner, no pytest needed
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

import composition_gate as cg  # noqa: E402


def verdict(separating, called, minority):
    return cg.decide_gate(separating, called, minority)[0]


# --- known real splits (validation) ----------------------------------------
def test_real_cw5_group_is_accepted():
    assert verdict(41214, 78252, 927) == cg.VERDICT_TWO          # 0.527


def test_spikein_one_percent_is_accepted():
    assert verdict(17580, 27960, 71) == cg.VERDICT_TWO           # 0.629


def test_spikein_one_percent_half_group_is_accepted():
    assert verdict(8816, 13593, 30) == cg.VERDICT_TWO            # 0.649


# --- known fake splits (validation) ----------------------------------------
def test_no_fetus_spikein_is_unresolved():
    assert verdict(2756, 12278, 3057) == cg.VERDICT_UNRESOLVED   # 0.224


def test_highest_validation_fake_is_unresolved():
    # Embryo-Imp15289124 per-library fake split
    assert verdict(15926, 52597, 6554) == cg.VERDICT_UNRESOLVED  # 0.303


def test_cw4_trio_early_fake_split_is_unresolved():
    assert verdict(16660, 63597, 9115) == cg.VERDICT_UNRESOLVED  # 0.262


def test_collapse_to_one_genotype_is_unresolved():
    # CW4 all five libraries: no locus called in both profiles
    assert verdict(0, 0, 56) == cg.VERDICT_UNRESOLVED
    assert "collapsed" in cg.decide_gate(0, 0, 56)[1]


# --- first production groups (checked, not used to choose the threshold) ---
def test_production_fetuses_ambient_withheld_are_accepted():
    assert verdict(9521, 13860, 313) == cg.VERDICT_TWO           # CW15, 0.687
    assert verdict(15673, 21791, 619) == cg.VERDICT_TWO          # Hrv43, 0.719


def test_lowest_production_real_is_accepted():
    assert verdict(25893, 59304, 6566) == cg.VERDICT_TWO         # D8, 0.437


def test_highest_production_fake_is_unresolved():
    assert verdict(2534, 6972, 5232) == cg.VERDICT_UNRESOLVED    # Hrv46, 0.363


# --- rule mechanics ---------------------------------------------------------
def test_small_minority_is_one_profile():
    assert verdict(500, 1000, 19) == cg.VERDICT_ONE
    assert verdict(500, 1000, 0) == cg.VERDICT_ONE


def test_boundaries_are_inclusive():
    assert verdict(400, 1000, 20) == cg.VERDICT_TWO
    assert verdict(399, 1000, 500) == cg.VERDICT_UNRESOLVED


def test_thresholds_can_be_overridden():
    assert cg.decide_gate(350, 1000, 500, {"min_discordant_fraction": 0.30})[0] == cg.VERDICT_TWO
    assert cg.decide_gate(500, 1000, 30, {"min_minority_singlets": 50})[0] == cg.VERDICT_ONE


def test_parse_ambient_reads_consensus_output():
    import tempfile
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as fh:
        fh.write("ambient RNA estimated as 1.802275914240492%\n")
        path = fh.name
    try:
        assert abs(cg.parse_ambient(path) - 1.802275914240492) < 1e-12
    finally:
        os.remove(path)
    assert cg.parse_ambient("/nonexistent/ambient_rna.txt") is None


def _run():
    tests = [v for k, v in sorted(globals().items())
             if k.startswith("test_") and callable(v)]
    failures = 0
    for t in tests:
        try:
            t()
            print("  PASS %s" % t.__name__)
        except AssertionError as e:
            failures += 1
            print("  FAIL %s: %s" % (t.__name__, e))
        except Exception as e:  # noqa: BLE001
            failures += 1
            print("  ERROR %s: %r" % (t.__name__, e))
    print("\n%d/%d passed" % (len(tests) - failures, len(tests)))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(_run())
