#!/usr/bin/env python3
"""Profile genotype comparison for the validation gates.

Two different questions are asked of profile VCFs, and they must not be
confused:

1. Within one run: do the inferred profiles actually differ? The acceptance
   criteria require at least 100 well-supported loci with discordant genotype
   calls before a two-profile interpretation is accepted.
2. Between two runs (seeds, discovery subsets, or backends): do the profiles
   agree? The criteria require at least 95% genotype concordance over loci
   called confidently in both.

Handles both backends' formats, verified against real output:
  Vireo      GT:AD:DP:PL   samples donor0, donor1
  Souporcell GT:AO:RO:T:E:GO:GN   samples 0, 1, with ./. for missing
Neither writes GP, and neither declares ##FORMAT headers.
"""

from __future__ import annotations

import gzip
import math
import os

# A genotype call is confident when its posterior, derived from PL when
# available, reaches this. With GT only, a non-missing call counts as confident
# and the caller must report the missing-call fraction alongside.
CONFIDENT_MIN_POSTERIOR = 0.90

GENOTYPE_CODES = {"0/0": 0, "0|0": 0, "0/1": 1, "0|1": 1, "1/0": 1, "1|0": 1,
                  "1/1": 2, "1|1": 2}


class GenotypeError(Exception):
    pass


def _open(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def _posterior_from_pl(pl_field: str) -> tuple[int, float] | None:
    """PL is phred-scaled likelihood, lower is better. Returns (call, posterior)."""
    try:
        values = [float(x) for x in pl_field.split(",")]
    except ValueError:
        return None
    if len(values) != 3:
        return None
    likelihoods = [10 ** (-value / 10) for value in values]
    total = sum(likelihoods)
    if total <= 0 or math.isnan(total):
        return None
    posteriors = [value / total for value in likelihoods]
    best = max(range(3), key=lambda i: posteriors[i])
    return best, posteriors[best]


def read_profile_vcf(path: str) -> dict:
    """Parse a profile VCF into per-sample genotype calls keyed by variant.

    Returns {"samples": [...], "calls": {(chrom,pos,ref,alt): [(call, posterior) | None, ...]}}
    where call is 0/1/2 alt-allele copies and None means no confident call.
    """
    if not os.path.isfile(path):
        raise GenotypeError(f"{path}: not found")

    samples: list[str] = []
    calls: dict[tuple[str, int, str, str], list] = {}
    n_missing_gt = 0
    n_records = 0

    with _open(path) as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = fields[9:]
                continue
            if not samples:
                raise GenotypeError(f"{path}: records before the #CHROM header")
            n_records += 1
            key = (fields[0], int(fields[1]), fields[3].upper(), fields[4].upper())
            format_keys = fields[8].split(":")
            record: list = []
            for sample_field in fields[9:]:
                parts = sample_field.split(":")
                data = dict(zip(format_keys, parts))
                gt = data.get("GT", "./.")
                if gt in (".", "./.", ".|."):
                    n_missing_gt += 1
                    record.append(None)
                    continue
                if "PL" in data:
                    parsed = _posterior_from_pl(data["PL"])
                    if parsed is None or parsed[1] < CONFIDENT_MIN_POSTERIOR:
                        record.append(None)
                        continue
                    record.append(parsed)
                    continue
                code = GENOTYPE_CODES.get(gt)
                record.append(None if code is None else (code, 1.0))
            calls[key] = record

    if not samples:
        raise GenotypeError(f"{path}: no #CHROM header found")

    return {
        "samples": samples,
        "calls": calls,
        "n_records": n_records,
        "n_missing_gt": n_missing_gt,
        "missing_gt_fraction": n_missing_gt / (n_records * len(samples)) if n_records else 0.0,
    }


def within_run_discordance(vcf_path: str) -> dict:
    """Do this run's own profiles differ? Counts well-supported discordant loci."""
    parsed = read_profile_vcf(vcf_path)
    samples = parsed["samples"]
    if len(samples) < 2:
        return {
            "n_profiles": len(samples),
            "n_loci_confident_in_all": 0,
            "n_discordant_loci": 0,
            "pairs": {},
            "missing_gt_fraction": parsed["missing_gt_fraction"],
        }

    pairs: dict[str, dict] = {}
    for i in range(len(samples)):
        for j in range(i + 1, len(samples)):
            confident = 0
            discordant = 0
            for record in parsed["calls"].values():
                a, b = record[i], record[j]
                if a is None or b is None:
                    continue
                confident += 1
                if a[0] != b[0]:
                    discordant += 1
            pairs[f"{samples[i]}|{samples[j]}"] = {
                "n_loci_confident_in_both": confident,
                "n_discordant_loci": discordant,
                "discordant_fraction": discordant / confident if confident else 0.0,
            }

    worst = min(pairs.values(), key=lambda entry: entry["n_discordant_loci"])
    return {
        "n_profiles": len(samples),
        "pairs": pairs,
        "min_discordant_loci_over_pairs": worst["n_discordant_loci"],
        "missing_gt_fraction": parsed["missing_gt_fraction"],
    }


def between_run_concordance(path_a: str, path_b: str,
                            mapping: dict[str, str] | None = None) -> dict:
    """Genotype concordance between two runs' profiles.

    `mapping` maps sample names in A to sample names in B. Without it, the best
    permutation is chosen, which is what "permutation-invariant profile
    matching" means: profile identity is arbitrary, only the partition matters.
    """
    import itertools

    a = read_profile_vcf(path_a)
    b = read_profile_vcf(path_b)
    shared_keys = set(a["calls"]) & set(b["calls"])

    def concordance_for(pairs: list[tuple[int, int]]) -> tuple[int, int]:
        confident = 0
        agree = 0
        for key in shared_keys:
            record_a = a["calls"][key]
            record_b = b["calls"][key]
            for index_a, index_b in pairs:
                call_a = record_a[index_a]
                call_b = record_b[index_b]
                if call_a is None or call_b is None:
                    continue
                confident += 1
                if call_a[0] == call_b[0]:
                    agree += 1
        return agree, confident

    if mapping:
        pairs = [(a["samples"].index(key), b["samples"].index(value))
                 for key, value in mapping.items()]
        agree, confident = concordance_for(pairs)
        best = {"mapping": mapping, "n_agree": agree, "n_loci": confident,
                "concordance": agree / confident if confident else 0.0}
    else:
        if len(a["samples"]) != len(b["samples"]):
            raise GenotypeError(
                f"profile counts differ: {len(a['samples'])} vs {len(b['samples'])}"
            )
        best = None
        for permutation in itertools.permutations(range(len(b["samples"]))):
            pairs = list(enumerate(permutation))
            agree, confident = concordance_for(pairs)
            score = agree / confident if confident else 0.0
            candidate = {
                "mapping": {a["samples"][i]: b["samples"][permutation[i]]
                            for i in range(len(a["samples"]))},
                "n_agree": agree,
                "n_loci": confident,
                "concordance": score,
            }
            if best is None or score > best["concordance"]:
                best = candidate

    best["n_shared_variants"] = len(shared_keys)
    best["missing_gt_fraction_a"] = a["missing_gt_fraction"]
    best["missing_gt_fraction_b"] = b["missing_gt_fraction"]
    return best

def discordant_locus_set(vcf_path: str) -> dict:
    """Which loci actually separate this run's two profiles.

    Needed because between-run genotype concordance is *not* a valid test of
    whether a second profile is real. When a group holds one donor, both
    inferred profiles are near-copies of that donor's genotype, so any pairing
    across two runs looks highly concordant. The discriminating question is
    whether the *difference* between the profiles replicates: a real second
    donor is separated by the same loci in every subset of cells, whereas a
    noise-driven split is separated by different loci each time.
    """
    parsed = read_profile_vcf(vcf_path)
    samples = parsed["samples"]
    if len(samples) < 2:
        return {"n_profiles": len(samples), "loci": set(), "n_confident_loci": 0}

    loci = set()
    confident = 0
    for key, record in parsed["calls"].items():
        a, b = record[0], record[1]
        if a is None or b is None:
            continue
        confident += 1
        if a[0] != b[0]:
            loci.add(key)
    return {"n_profiles": len(samples), "loci": loci, "n_confident_loci": confident}


def discordant_locus_replication(vcf_a: str, vcf_b: str) -> dict:
    """Do the same loci separate the profiles in two disjoint cell subsets?

    Reports the overlap of the two discordant-locus sets, restricted to loci
    confidently called in both runs so that coverage differences do not masquer-
    ade as disagreement. `jaccard` near 1 means the separating genotype is
    reproducible; near 0 means each subset invented its own split.
    """
    a = discordant_locus_set(vcf_a)
    b = discordant_locus_set(vcf_b)
    if a["n_profiles"] < 2 or b["n_profiles"] < 2:
        return {"error": "both runs must contain at least two profiles"}

    parsed_a = read_profile_vcf(vcf_a)
    parsed_b = read_profile_vcf(vcf_b)

    def confident_keys(parsed):
        return {key for key, record in parsed["calls"].items()
                if record[0] is not None and record[1] is not None}

    shared = confident_keys(parsed_a) & confident_keys(parsed_b)
    set_a = a["loci"] & shared
    set_b = b["loci"] & shared
    union = set_a | set_b
    intersection = set_a & set_b

    # What overlap would chance alone give, if each set were drawn at random
    # from the shared loci? Reported so a modest Jaccard is interpretable.
    expected = (len(set_a) * len(set_b) / len(shared)) if shared else 0.0

    return {
        "n_shared_confident_loci": len(shared),
        "n_discordant_a": len(set_a),
        "n_discordant_b": len(set_b),
        "n_discordant_in_both": len(intersection),
        "jaccard": len(intersection) / len(union) if union else 0.0,
        "fraction_of_a_replicated": len(intersection) / len(set_a) if set_a else 0.0,
        "fraction_of_b_replicated": len(intersection) / len(set_b) if set_b else 0.0,
        "expected_overlap_by_chance": expected,
        "enrichment_over_chance": (len(intersection) / expected) if expected > 0 else None,
    }
