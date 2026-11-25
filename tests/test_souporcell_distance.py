#!/usr/bin/env python3
"""Unit tests for the fetal-mode Souporcell cross-library check (bin/souporcell_distance.py).

Tiny Souporcell-shaped runs are written to a temporary folder: the fetus's libraries split into two
halves of the same genotype, and one library whose second cluster carries a different genotype
(a "mother"). The check must stay silent on the former and warn on the latter, and the report must
embed the heatmap. The real-data check against the paper's numbers is
tests/regression_souporcell_distance.py (LSF).

Runnable two ways:
    pytest tests/test_souporcell_distance.py
    python3 tests/test_souporcell_distance.py      # built-in mini-runner, no pytest needed
"""
import os
import sys
import tempfile

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

import make_fetal_report  # noqa: E402
import souporcell_distance as sd  # noqa: E402

N_SITES, CELLS, READS = 400, 30, 2
RNG = np.random.default_rng(1)
FETUS = RNG.choice([0.0, 0.5, 1.0], N_SITES)            # ref-allele fraction per site
MOTHER = np.where(RNG.random(N_SITES) < 0.5, FETUS, RNG.choice([0.0, 0.5, 1.0], N_SITES))


def write_run(d, genotypes):
    """One Souporcell-shaped run: genotypes[c] is cluster c's ref fraction per site; CELLS singlets each."""
    os.makedirs(d)
    with open(f"{d}/common_variants_covered.vcf", "w") as fh:
        fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\n")
        fh.writelines(f"chr1\t{1000 + 10 * i}\t.\tA\tG\n" for i in range(N_SITES))
    cells, alt, ref = [], [], []
    for c, g in enumerate(genotypes):
        for k in range(CELLS):
            cell = len(cells) + 1
            cells.append((f"BC{cell}", "singlet", str(c)))
            r = RNG.binomial(READS, g)
            for i in range(N_SITES):
                alt.append((i + 1, cell, READS - r[i]))
                ref.append((i + 1, cell, r[i]))
    cells.append((f"BC{len(cells) + 1}", "doublet", "0/1"))     # doublets are ignored
    pd.DataFrame(cells, columns=["barcode", "status", "assignment"]).assign(cluster0=0, cluster1=0).to_csv(
        f"{d}/clusters.tsv", sep="\t", index=False)
    for name, rows in (("alt.mtx", alt), ("ref.mtx", ref)):
        with open(f"{d}/{name}", "w") as fh:
            fh.write(f"%%MatrixMarket matrix coordinate integer general\n%\n{N_SITES} {len(cells)} {len(rows)}\n")
            fh.writelines(f"{i} {j} {v}\n" for i, j, v in rows)


def run_check(tmp, with_mother):
    rows = []
    for lib in ("L1", "L2", "L3"):
        g = [FETUS, MOTHER] if (with_mother and lib == "L2") else [FETUS, FETUS]
        write_run(f"{tmp}/{lib}", g)
        rows.append(f"F1\t{lib}\t{tmp}/{lib}\n")
    write_run(f"{tmp}/U1", [MOTHER, MOTHER])                    # a second, unrelated "fetus"
    rows.append(f"F2\tU1\t{tmp}/U1\n")
    with open(f"{tmp}/manifest.tsv", "w") as fh:
        fh.write("group\tlibrary\tdir\n" + "".join(rows))
    sd.main(["--manifest", f"{tmp}/manifest.tsv", "--outdir", f"{tmp}/out"])
    return f"{tmp}/out"


def test_one_fetus_split_in_two_is_silent():
    with tempfile.TemporaryDirectory() as tmp:
        out = run_check(tmp, with_mother=False)
        lib = pd.read_csv(f"{out}/libraries.tsv", sep="\t")
        f1 = lib[lib.group == "F1"]
        assert (f1.status == "ok").all(), f1
        assert f1.msd_to_other_libraries_median.max() < 0.03
        w = pd.read_csv(f"{out}/warnings.tsv", sep="\t")
        assert not (w.group == "F1").any()


def test_mother_cluster_is_warned_and_only_it():
    with tempfile.TemporaryDirectory() as tmp:
        out = run_check(tmp, with_mother=True)
        lib = pd.read_csv(f"{out}/libraries.tsv", sep="\t").set_index(["library", "cluster"])
        assert lib.loc[("L2", 1), "status"] == "review"
        assert (lib.drop(("L2", 1)).query("group == 'F1'").status == "ok").all()
        w = pd.read_csv(f"{out}/warnings.tsv", sep="\t")
        assert list(w[w.group == "F1"].library) == ["L2"]


def test_single_library_fetus_is_warned_as_not_checkable():
    with tempfile.TemporaryDirectory() as tmp:
        out = run_check(tmp, with_mother=False)
        w = pd.read_csv(f"{out}/warnings.tsv", sep="\t")
        assert "not possible" in w[w.group == "F2"].warning.iloc[0]
        assert (pd.read_csv(f"{out}/libraries.tsv", sep="\t").query("group == 'F2'").status == "no_other_library").all()


def test_fetus_matrix_and_heatmaps_are_written():
    with tempfile.TemporaryDirectory() as tmp:
        out = run_check(tmp, with_mother=False)
        m = pd.read_csv(f"{out}/groups.tsv", sep="\t", index_col=0)
        assert m.loc["F1", "F1"] < 0.03 < m.loc["F1", "F2"]
        for name in ("F1.heatmap.svg", "F2.heatmap.svg", "groups.heatmap.svg"):
            assert open(f"{out}/{name}").read().startswith("<svg")


def test_doublets_and_small_clusters_are_not_judged():
    cl = pd.DataFrame({"group": "F1", "library": ["A", "A", "B", "B"], "cluster": [0, 1, 0, 1],
                       "singlets": [500, 5, 400, 300]})
    p = pd.DataFrame([("F1", "A", 0, "F1", "A", 1, 0.09), ("F1", "A", 0, "F1", "B", 0, 0.005),
                      ("F1", "A", 0, "F1", "B", 1, 0.006), ("F1", "A", 1, "F1", "B", 0, 0.09),
                      ("F1", "A", 1, "F1", "B", 1, 0.09), ("F1", "B", 0, "F1", "B", 1, 0.008)],
                     columns=["group_a", "library_a", "cluster_a", "group_b", "library_b", "cluster_b", "msd"])
    lib, warns = sd.per_library(cl, p, min_singlets=20, warn_msd=0.03)
    assert lib.set_index(["library", "cluster"]).loc[("A", 1), "status"] == "too_few_cells"
    assert warns.empty                                        # the 5-cell cluster is far but not judged


def test_report_embeds_heatmap_and_check_warnings():
    with tempfile.TemporaryDirectory() as tmp:
        out = run_check(tmp, with_mother=True)
        for name, text in (("summary.tsv", "scope\tid\tcells\tfetal\ttoo_little_data\tmaternal\tmaternal_in_pair\t"
                                           "fetal_fraction_of_cells\tupper95_maternal_among_fetal\tflagged_fraction\n"
                                           "fetus\tF1\t100\t90\t10\t0\t0\t0.9\t0.033\t0\n"
                                           "batch\tF1_b1\t100\t90\t10\t0\t0\t0.9\t0.033\t0\n"),
                           ("sex.tsv", "sex_call\ty_positive_fraction\nfemale\t0.0\n"),
                           ("warnings.txt", ""),
                           ("flagged.tsv", "library\tbarcode\tfetal_mode\n")):
            with open(f"{tmp}/{name}", "w") as fh:
                fh.write(text)
        make_fetal_report.main(["--group-id", "F1", "--summary", f"{tmp}/summary.tsv", "--sex-check", f"{tmp}/sex.tsv",
                                "--warnings", f"{tmp}/warnings.txt", "--flagged", f"{tmp}/flagged.tsv",
                                "--souporcell-check", out, "--out", f"{tmp}/r.html"])
        html = open(f"{tmp}/r.html").read()
        assert "<svg" in html and "Souporcell cross-library check" in html
        assert "L2 cluster 1" in html and "conf-prov" in html


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
