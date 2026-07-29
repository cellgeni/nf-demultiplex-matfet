#!/usr/bin/env python3
"""Unit tests for the maternal/fetal decision engine (bin/label_clusters.py).

Covers ROADMAP.md §12: RPS4Y1 precedence, the absent-RPS4Y1 safeguard, the female
secondary-marker pattern, vote margins, QC gates, and determinism.

Runnable two ways:
    pytest tests/test_decision_engine.py
    python3 tests/test_decision_engine.py       # built-in mini-runner, no pytest needed
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

import label_clusters as lc  # noqa: E402

TIERS = {
    "RPS4Y1": "rps4y1",
    "HAPLN1": "secondary", "PITX2": "secondary", "GAPLINC": "secondary",
    "HPGDS": "secondary", "XAGE2": "secondary", "XAGE3": "secondary", "HLA-G": "secondary",
}


def make_stats(n_cells, data, tiers=TIERS):
    """data: {cluster: {marker: (frac_pos, cp10k)}}."""
    markers = {}
    for c, md in data.items():
        markers[c] = {
            m: {"frac_pos": fp, "cp10k": cp, "n_pos": int(round(fp * n_cells[c]))}
            for m, (fp, cp) in md.items()
        }
    return {
        "clusters": sorted(data),
        "n_cells": dict(n_cells),
        "markers": markers,
        "tiers": dict(tiers),
    }


def _all_markers(rps4y1=(0.0, 0.0), secondary=(0.0, 0.0)):
    d = {"RPS4Y1": rps4y1}
    for m in ("HAPLN1", "PITX2", "GAPLINC", "HPGDS", "XAGE2", "XAGE3", "HLA-G"):
        d[m] = secondary
    return d


# --------------------------------------------------------------------------- #
# Tier 1 — RPS4Y1 override
# --------------------------------------------------------------------------- #
def test_rps4y1_identifies_cluster0_fetal():
    stats = make_stats(
        {"0": 200, "1": 200},
        {"0": _all_markers(rps4y1=(0.30, 50.0)), "1": _all_markers(rps4y1=(0.01, 0.5))},
    )
    r = lc.decide(stats)
    assert r["fetal_cluster"] == "0"
    assert r["maternal_cluster"] == "1"
    assert r["tier"] == "rps4y1"
    assert r["confidence"] == "high"
    assert r["inferred_sex"] == "male"
    assert r["labels"]["0"] == lc.LABEL_FETAL_MALE
    assert r["labels"]["1"] == lc.LABEL_MATERNAL


def test_rps4y1_identifies_cluster1_fetal():
    stats = make_stats(
        {"0": 200, "1": 200},
        {"0": _all_markers(rps4y1=(0.005, 0.2)), "1": _all_markers(rps4y1=(0.40, 80.0))},
    )
    r = lc.decide(stats)
    assert r["fetal_cluster"] == "1"
    assert r["maternal_cluster"] == "0"
    assert r["tier"] == "rps4y1"


def test_rps4y1_takes_precedence_over_conflicting_secondary_markers():
    # RPS4Y1 in cluster 0; ALL secondary markers strong in cluster 1.
    # Precedence invariant: RPS4Y1 wins -> fetal is cluster 0.
    stats = make_stats(
        {"0": 200, "1": 200},
        {
            "0": _all_markers(rps4y1=(0.35, 60.0), secondary=(0.0, 0.0)),
            "1": _all_markers(rps4y1=(0.0, 0.0), secondary=(0.5, 80.0)),
        },
    )
    r = lc.decide(stats)
    assert r["fetal_cluster"] == "0"
    assert r["tier"] == "rps4y1"
    assert r["supporting_markers"] == ["RPS4Y1"]


def test_rps4y1_not_enriched_does_not_override():
    # Present in both near the floor, low fold -> Tier 1 must NOT fire.
    stats = make_stats(
        {"0": 200, "1": 200},
        {"0": _all_markers(rps4y1=(0.12, 3.0)), "1": _all_markers(rps4y1=(0.10, 2.5))},
    )
    r = lc.decide(stats)
    assert r["tier"] != "rps4y1"
    assert any("BOTH clusters" in w for w in r["warnings"])


def test_rps4y1_both_positive_but_enriched_fires_with_warning():
    stats = make_stats(
        {"0": 200, "1": 200},
        {"0": _all_markers(rps4y1=(0.60, 90.0)), "1": _all_markers(rps4y1=(0.11, 3.0))},
    )
    r = lc.decide(stats)
    assert r["tier"] == "rps4y1"
    assert r["fetal_cluster"] == "0"
    assert any("BOTH clusters" in w for w in r["warnings"])


# --------------------------------------------------------------------------- #
# The safeguard: absent RPS4Y1 is never maternal evidence
# --------------------------------------------------------------------------- #
def test_absent_rps4y1_with_no_other_evidence_is_unresolved_not_maternal():
    stats = make_stats(
        {"0": 200, "1": 200},
        {"0": _all_markers(), "1": _all_markers()},  # everything zero
    )
    r = lc.decide(stats)
    assert r["tier"] == "unresolved"
    assert r["fetal_cluster"] is None
    assert r["maternal_cluster"] is None
    # No cluster may be labeled maternal when nothing identified a fetal cluster.
    assert lc.LABEL_MATERNAL not in r["labels"].values()
    assert all(v == lc.LABEL_UNRESOLVED for v in r["labels"].values())


def test_female_pattern_supported_only_by_secondary_markers():
    stats = make_stats(
        {"0": 200, "1": 200},
        {
            "0": _all_markers(rps4y1=(0.0, 0.0), secondary=(0.02, 0.5)),
            "1": _all_markers(rps4y1=(0.0, 0.0), secondary=(0.35, 40.0)),
        },
    )
    r = lc.decide(stats)
    assert r["fetal_cluster"] == "1"
    assert r["tier"] == "markers"
    assert r["confidence"] == "provisional"
    assert r["inferred_sex"] == "likely_female"
    assert r["labels"]["1"] == lc.LABEL_FETAL_MARKER
    assert r["labels"]["0"] == lc.LABEL_MATERNAL
    assert len(r["supporting_markers"]) >= 2


# --------------------------------------------------------------------------- #
# Tier 2 — vote margins / ambiguity
# --------------------------------------------------------------------------- #
def test_single_secondary_marker_below_margin_is_unresolved():
    d1 = _all_markers(rps4y1=(0.0, 0.0), secondary=(0.0, 0.0))
    d1["HLA-G"] = (0.4, 50.0)  # only one marker enriched -> margin 1 < 2
    stats = make_stats({"0": 200, "1": 200}, {"0": _all_markers(), "1": d1})
    r = lc.decide(stats)
    assert r["tier"] == "unresolved"
    assert r["fetal_cluster"] is None


def test_two_secondary_markers_reach_margin_and_call():
    d1 = _all_markers(rps4y1=(0.0, 0.0), secondary=(0.0, 0.0))
    d1["HLA-G"] = (0.4, 50.0)
    d1["XAGE3"] = (0.3, 40.0)
    stats = make_stats({"0": 200, "1": 200}, {"0": _all_markers(), "1": d1})
    r = lc.decide(stats)
    assert r["fetal_cluster"] == "1"
    assert r["tier"] == "markers"
    assert set(r["supporting_markers"]) == {"HLA-G", "XAGE3"}


def test_markers_similar_in_both_clusters_is_unresolved():
    stats = make_stats(
        {"0": 200, "1": 200},
        {
            "0": _all_markers(secondary=(0.20, 20.0)),
            "1": _all_markers(secondary=(0.20, 20.0)),
        },
    )
    r = lc.decide(stats)
    assert r["tier"] == "unresolved"
    assert all(v == 0 for v in r["votes"].values())


# --------------------------------------------------------------------------- #
# QC gates
# --------------------------------------------------------------------------- #
def test_single_cluster_is_unresolved():
    stats = make_stats({"0": 300}, {"0": _all_markers(rps4y1=(0.5, 80.0))})
    r = lc.decide(stats)
    assert r["tier"] == "unresolved"
    assert any("Fewer than two" in w for w in r["warnings"])


def test_small_cluster_downgrades_high_to_provisional():
    stats = make_stats(
        {"0": 200, "1": 10},  # cluster 1 below default min_cluster_cells (50)
        {"0": _all_markers(rps4y1=(0.30, 50.0)), "1": _all_markers(rps4y1=(0.0, 0.0))},
    )
    r = lc.decide(stats)
    assert r["tier"] == "rps4y1"
    assert r["confidence"] == "provisional"
    assert any("singlet cells" in w for w in r["warnings"])


def test_cluster_size_imbalance_warns():
    stats = make_stats(
        {"0": 2000, "1": 50},
        {"0": _all_markers(rps4y1=(0.30, 50.0)), "1": _all_markers()},
    )
    r = lc.decide(stats)
    assert any("imbalanced" in w for w in r["warnings"])


# --------------------------------------------------------------------------- #
# Structure / determinism
# --------------------------------------------------------------------------- #
def test_build_cluster_stats_from_rows():
    rows = [
        {"cluster": "0", "n_cluster_cells": "100", "marker": "RPS4Y1", "tier": "rps4y1",
         "frac_pos": "0.3", "cp10k": "50", "n_pos": "30"},
        {"cluster": "1", "n_cluster_cells": "120", "marker": "RPS4Y1", "tier": "rps4y1",
         "frac_pos": "0.01", "cp10k": "0.5", "n_pos": "1"},
    ]
    stats = lc.build_cluster_stats(rows)
    assert stats["clusters"] == ["0", "1"]
    assert stats["n_cells"] == {"0": 100, "1": 120}
    assert stats["markers"]["0"]["RPS4Y1"]["frac_pos"] == 0.3
    assert stats["tiers"]["RPS4Y1"] == "rps4y1"


def test_decision_is_deterministic():
    stats = make_stats(
        {"0": 200, "1": 200},
        {"0": _all_markers(rps4y1=(0.30, 50.0)), "1": _all_markers(rps4y1=(0.01, 0.5))},
    )
    r1 = lc.decide(stats)
    r2 = lc.decide(stats)
    assert r1 == r2


# --------------------------------------------------------------------------- #
# Mini-runner so the file works without pytest installed
# --------------------------------------------------------------------------- #
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
