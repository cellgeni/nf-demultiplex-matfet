#!/usr/bin/env python3
"""Self-test for the shared assignment/genotype/comparison contract.

Synthetic only, stdlib only, safe on a head node. Exits 0 on success.
"""

from __future__ import annotations

import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "common"))

import assignment_table as at  # noqa: E402
import compare_runs as cr  # noqa: E402
import genotype_vcf as gv  # noqa: E402

failures: list[str] = []


def check(name: str, condition: bool, detail: str = "") -> None:
    if condition:
        print(f"PASS  {name}")
    else:
        failures.append(f"{name}: {detail}")
        print(f"FAIL  {name}: {detail}")


def row(barcode, sample, status, profile, prob, doublet=0.0, loci=100):
    return {
        "group_barcode": barcode, "sample_id": sample, "status": status,
        "profile": profile, "assign_prob": prob, "doublet_prob": doublet,
        "n_informative_loci": loci, "backend": "test", "mode": "discovery",
        "native_status": status, "native_profile": profile,
    }


def write_vcf(path, records, samples, fmt="GT"):
    with open(path, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
        for (chrom, pos, ref, alt), calls in records:
            out.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\t{fmt}\t" + "\t".join(calls) + "\n")


with tempfile.TemporaryDirectory() as tmp:
    # --- table validation rejects contract violations -----------------------
    good = [row("libA::AAA", "libA", "singlet", "profile_1", 0.99),
            row("libA::BBB", "libA", "doublet", "", 0.4, 0.95),
            row("libA::CCC", "libA", "unassigned", "", 0.2, 0.1)]
    path = os.path.join(tmp, "good.tsv")
    at.write_table(path, good)
    check("valid table round-trips", len(at.read_table(path)) == 3)

    for name, bad in (
        ("singlet without profile", [row("x::1", "x", "singlet", "", 0.9)]),
        ("non-singlet with profile", [row("x::1", "x", "doublet", "profile_1", 0.9, 0.9)]),
        ("bad status", [row("x::1", "x", "maybe", "", 0.9)]),
        ("probability out of range", [row("x::1", "x", "singlet", "profile_1", 1.5)]),
    ):
        try:
            at.write_table(os.path.join(tmp, "bad.tsv"), bad)
            check(f"rejects {name}", False, "no AdapterError raised")
        except at.AdapterError:
            check(f"rejects {name}", True)

    # --- high-confidence definition -----------------------------------------
    check("0.90 is high confidence", at.is_high_confidence(row("a", "a", "singlet", "profile_1", 0.90)))
    check("0.89 is not", not at.is_high_confidence(row("a", "a", "singlet", "profile_1", 0.89)))
    check("confident doublet is not high-confidence singlet",
          not at.is_high_confidence(row("a", "a", "doublet", "", 0.99, 0.99)))

    # --- summarise ----------------------------------------------------------
    summary = at.summarise(good)
    check("summary counts statuses",
          (summary["group"]["n_singlet"], summary["group"]["n_doublet"],
           summary["group"]["n_unassigned"]) == (1, 1, 1),
          str(summary["group"]))

    # --- ARI ----------------------------------------------------------------
    check("ARI identical labels is 1",
          abs(cr.adjusted_rand_index(["a", "a", "b", "b"], ["x", "x", "y", "y"]) - 1.0) < 1e-12)
    check("ARI relabelled is 1",
          abs(cr.adjusted_rand_index(["a", "b"], ["b", "a"]) - 1.0) < 1e-12)
    ari_mixed = cr.adjusted_rand_index(["a", "a", "b", "b"], ["x", "y", "x", "y"])
    check("ARI unrelated partition is near 0", abs(ari_mixed) < 0.5, str(ari_mixed))

    # --- permutation-invariant matching -------------------------------------
    run_a = [row(f"libA::{i}", "libA", "singlet", "profile_1" if i < 5 else "profile_2", 0.99)
             for i in range(10)]
    # same partition, labels swapped: agreement must be 1.0 after matching
    run_b = [row(f"libA::{i}", "libA", "singlet", "profile_2" if i < 5 else "profile_1", 0.99)
             for i in range(10)]
    matching = cr.best_profile_matching(run_a, run_b)
    check("swapped labels match perfectly", matching["concordance"] == 1.0, str(matching))
    check("matching reports the label map",
          matching["mapping"] == {"profile_1": "profile_2", "profile_2": "profile_1"},
          str(matching["mapping"]))

    # low-confidence cells are excluded from the comparison
    run_c = [dict(r) for r in run_b]
    run_c[0]["assign_prob"] = 0.5
    check("low-confidence cells excluded",
          cr.best_profile_matching(run_a, run_c)["n_common_high_confidence"] == 9)

    # a library present in only one run contributes nothing
    check("non-overlapping barcodes give no comparison",
          cr.best_profile_matching(run_a, [row("libZ::1", "libZ", "singlet", "profile_1", 0.99)])
          ["n_common_high_confidence"] == 0)

    # --- genotype VCF: GT only (souporcell shape) ---------------------------
    sc = os.path.join(tmp, "souporcell.vcf")
    write_vcf(sc, [
        (("chr1", 100, "A", "G"), ["0/0", "1/1"]),
        (("chr1", 200, "C", "T"), ["0/1", "0/1"]),
        (("chr1", 300, "G", "A"), ["./.", "1/1"]),
    ], ["0", "1"], fmt="GT")
    within = gv.within_run_discordance(sc)
    check("within-run discordance counts only loci confident in both",
          within["pairs"]["0|1"] == {"n_loci_confident_in_both": 2, "n_discordant_loci": 1,
                                     "discordant_fraction": 0.5},
          str(within["pairs"]))
    check("missing GT fraction reported",
          abs(within["missing_gt_fraction"] - 1 / 6) < 1e-9, str(within["missing_gt_fraction"]))

    # --- genotype VCF: PL-based confidence (vireo shape) --------------------
    vr = os.path.join(tmp, "vireo.vcf")
    write_vcf(vr, [
        # confident 0/0 for donor0 (PL 0 for first genotype), confident 1/1 for donor1
        (("chr1", 100, "A", "G"), ["0/0:0:10:0,60,120", "1/1:10:10:120,60,0"]),
        # donor0 ambiguous: two genotypes nearly equally likely -> not confident
        (("chr1", 200, "C", "T"), ["0/1:5:10:0,1,120", "0/0:0:10:0,60,120"]),
    ], ["donor0", "donor1"], fmt="GT:AD:DP:PL")
    parsed = gv.read_profile_vcf(vr)
    check("PL-confident call kept", parsed["calls"][("chr1", 100, "A", "G")][0][0] == 0)
    check("PL-ambiguous call dropped", parsed["calls"][("chr1", 200, "C", "T")][0] is None)

    # --- between-run concordance is permutation invariant -------------------
    vr_swapped = os.path.join(tmp, "vireo_swapped.vcf")
    write_vcf(vr_swapped, [
        (("chr1", 100, "A", "G"), ["1/1:10:10:120,60,0", "0/0:0:10:0,60,120"]),
        (("chr1", 200, "C", "T"), ["0/0:0:10:0,60,120", "0/1:5:10:0,1,120"]),
    ], ["donor0", "donor1"], fmt="GT:AD:DP:PL")
    concordance = gv.between_run_concordance(vr, vr_swapped)
    check("swapped profiles are perfectly concordant",
          concordance["concordance"] == 1.0 and concordance["n_loci"] > 0, str(concordance))
    check("concordance reports the profile map",
          concordance["mapping"] == {"donor0": "donor1", "donor1": "donor0"},
          str(concordance["mapping"]))

    # --- end-to-end compare, including verdicts -----------------------------
    table_a = os.path.join(tmp, "a.tsv")
    table_b = os.path.join(tmp, "b.tsv")
    at.write_table(table_a, run_a)
    at.write_table(table_b, run_b)
    report = cr.compare(table_a, table_b, sc, sc, "seed")
    check("identical partitions pass the cell criterion",
          report["verdicts"]["cell_concordance"]["passed"] is True)
    check("too few discordant loci fails the two-profile criterion",
          report["verdicts"]["discordant_loci_a"]["passed"] is False,
          str(report["verdicts"].get("discordant_loci_a")))
    check("overall verdict reflects a failed criterion",
          report["all_evaluated_criteria_passed"] is False)

    report_no_vcf = cr.compare(table_a, table_b, None, None, "seed")
    check("missing VCFs leave genotype criteria unevaluated rather than passed",
          "genotype_concordance" not in report_no_vcf["verdicts"])

    # --- minority retention: the blind spot that agreement alone misses ------
    # Reference: 20 majority + 5 minority cells, all high-confidence singlets.
    ref = [row(f"libA::{i}", "libA", "singlet", "profile_1", 0.99) for i in range(20)]
    ref += [row(f"libA::m{i}", "libA", "singlet", "profile_2", 0.99) for i in range(5)]
    # Test run keeps the majority but demotes every minority cell to doublet.
    test = [row(f"libA::{i}", "libA", "singlet", "profile_1", 0.99) for i in range(20)]
    test += [row(f"libA::m{i}", "libA", "doublet", "", 0.4, 0.99) for i in range(5)]

    matching = cr.best_profile_matching(ref, test)
    check("agreement is blind to a discarded minority",
          matching["concordance"] == 1.0, str(matching["concordance"]))
    retention = cr.minority_retention(ref, test, matching["mapping"])
    check("retention catches the discarded minority",
          retention["retention"] == 0.0 and retention["fates"]["doublet"] == 5,
          str(retention))
    check("retention reports doublet inflation",
          abs(retention["doublet_inflation"] - 0.2) < 1e-9,
          str(retention["doublet_inflation"]))

    # A test run that keeps the minority scores full retention.
    kept = [row(f"libA::{i}", "libA", "singlet", "profile_1", 0.99) for i in range(20)]
    kept += [row(f"libA::m{i}", "libA", "singlet", "profile_2", 0.99) for i in range(5)]
    retention_ok = cr.minority_retention(ref, kept, cr.best_profile_matching(ref, kept)["mapping"])
    check("retention is 1.0 when the minority survives", retention_ok["retention"] == 1.0,
          str(retention_ok["retention"]))

    # Swapped labels must not look like a lost minority.
    swapped = [row(f"libA::{i}", "libA", "singlet", "profile_2", 0.99) for i in range(20)]
    swapped += [row(f"libA::m{i}", "libA", "singlet", "profile_1", 0.99) for i in range(5)]
    retention_swap = cr.minority_retention(
        ref, swapped, cr.best_profile_matching(ref, swapped)["mapping"])
    check("retention is permutation invariant", retention_swap["retention"] == 1.0,
          str(retention_swap["retention"]))

    # A single-profile reference has no minority to retain.
    single = [row(f"libA::{i}", "libA", "singlet", "profile_1", 0.99) for i in range(10)]
    check("single-profile reference is not evaluable",
          cr.minority_retention(single, single, {"profile_1": "profile_1"})["evaluable"] is False)

    # End to end: the verdict must fail even though agreement is perfect.
    table_ref = os.path.join(tmp, "ref.tsv"); at.write_table(table_ref, ref)
    table_test = os.path.join(tmp, "test.tsv"); at.write_table(table_test, test)
    rep = cr.compare(table_ref, table_test, None, None, "projection-vs-joint")
    check("a discarded minority fails the criteria despite perfect agreement",
          rep["verdicts"]["cell_concordance"]["passed"] is True
          and rep["verdicts"]["minority_retention"]["passed"] is False,
          str(rep["verdicts"]))

    # --- informative loci counted from real matrix shapes -------------------
    group_dir = os.path.join(tmp, "GRP")
    os.makedirs(group_dir)
    with open(os.path.join(group_dir, "barcodes.tsv"), "w") as out:
        out.write("libA::AAA\nlibA::BBB\n")
    # cell 1 has one covered locus and one explicit zero; cell 2 has none
    for name, values in (("alt.mtx", [3, 0, 0]), ("ref.mtx", [1, 0, 0])):
        with open(os.path.join(group_dir, name), "w") as out:
            out.write("%%MatrixMarket matrix coordinate real general\n% test\n2 2 3\n")
            out.write(f"1 1 {values[0]}\n2 1 {values[1]}\n1 2 {values[2]}\n")
    counts = at.count_informative_loci(group_dir)
    check("explicit zeros are not counted as coverage",
          counts == {"libA::AAA": 1, "libA::BBB": 0}, str(counts))

print()
if failures:
    print(f"{len(failures)} FAILED")
    for failure in failures:
        print(f"  {failure}")
    sys.exit(1)
print("ALL COMMON-CONTRACT TESTS PASSED")
