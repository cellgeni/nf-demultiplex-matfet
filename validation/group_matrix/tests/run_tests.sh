#!/usr/bin/env bash
# Local self-test for the group-matrix harness: tiny synthetic fixtures only.
# Safe on a head node. Real-data builds must run under LSF.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HARNESS="$(dirname "$HERE")"
WORK="${1:-${TMPDIR:-/tmp}/group_matrix_selftest}"
rm -rf "$WORK"
mkdir -p "$WORK"

GROUPS_CSV="$(python3 "$HERE/make_fixtures.py" "$WORK/fixtures")"
fails=0

expect_fail() {
  local group="$1"
  python3 "$HARNESS/build_group_matrix.py" --groups "$GROUPS_CSV" --group-id "$group" \
    --outdir "$WORK/out_$group" --no-hash >"$WORK/$group.log" 2>&1
  if [[ $? -eq 2 ]]; then
    echo "PASS  $group rejected: $(grep -m1 'CONTRACT FAILURE' "$WORK/$group.log")"
  else
    echo "FAIL  $group was not rejected"
    cat "$WORK/$group.log"
    fails=$((fails + 1))
  fi
}

echo "== negative cases =="
for group in DUP CONFLICT MULTI DIMS COORD NEG UNKNOWN_CONTIG; do
  expect_fail "$group"
done

echo
echo "== GOOD build =="
for run in a b; do
  python3 "$HARNESS/build_group_matrix.py" --groups "$GROUPS_CSV" --group-id GOOD \
    --outdir "$WORK/out_good_$run" || { echo "FAIL  GOOD build $run"; fails=$((fails + 1)); }
done

echo
echo "== GOOD expected values =="
check_line() {
  local description="$1" expected="$2" actual="$3"
  if [[ "$expected" == "$actual" ]]; then
    echo "PASS  $description"
  else
    echo "FAIL  $description: expected '$expected', got '$actual'"
    fails=$((fails + 1))
  fi
}
G="$WORK/out_good_a/GOOD"
check_line "group dimensions" "4 5 7" "$(sed -n '3p' "$G/alt.mtx")"
check_line "union variant order" "chr1 100|chr1 200|chr1 300|chr2 50" \
  "$(grep -v '^#' "$G/variants.vcf" | cut -f1,2 | tr '\t' ' ' | paste -sd'|')"
check_line "barcodes disambiguated" "libA::AAA-1|libA::CCC-1|libB::AAA-1|libB::GGG-1|libB::CCC-1" \
  "$(paste -sd'|' "$G/barcodes.tsv")"
check_line "libB missing chr1:300 contributes no entries" "0" \
  "$(awk 'NR>3 && $1==3 && $2>2' "$G/alt.mtx" | wc -l)"
check_line "explicit zeros preserved" "3" "$(awk 'NR>3 && $3==0' "$G/alt.mtx" | wc -l)"

echo
echo "== integrity check (includes determinism) =="
python3 "$HARNESS/check_group_matrix.py" --groups "$GROUPS_CSV" --group-id GOOD \
  --outdir "$WORK/out_good_a" --compare-dir "$WORK/out_good_b" || fails=$((fails + 1))

echo
echo "== corrupted output must be caught =="
CORRUPT="$WORK/out_corrupt"
cp -r "$WORK/out_good_a" "$CORRUPT"
awk 'NR==4 {$3=$3+1} {print}' "$CORRUPT/GOOD/alt.mtx" >"$CORRUPT/GOOD/alt.tmp" \
  && mv "$CORRUPT/GOOD/alt.tmp" "$CORRUPT/GOOD/alt.mtx"
python3 "$HARNESS/check_group_matrix.py" --groups "$GROUPS_CSV" --group-id GOOD \
  --outdir "$CORRUPT" >"$WORK/corrupt.log" 2>&1
if [[ $? -eq 1 ]]; then
  echo "PASS  tampered allele count detected"
else
  echo "FAIL  tampered allele count not detected"
  cat "$WORK/corrupt.log"
  fails=$((fails + 1))
fi

echo
if [[ $fails -eq 0 ]]; then
  echo "ALL TESTS PASSED (workdir: $WORK)"
else
  echo "$fails TEST(S) FAILED (workdir: $WORK)"
fi
exit $((fails > 0))
