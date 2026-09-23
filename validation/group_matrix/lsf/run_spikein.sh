#!/usr/bin/env bash
# Fetal-fraction detection-limit experiment, both backends, one dependency DAG.
#
# Neither CW4 nor D6 can say whether "no reproducible second genotype" means no
# fetal cells or a fetal population below the detection floor. This measures the
# floor directly: consensus-labelled fetal cells from CW5 are spiked into a
# consensus-maternal background drawn from all three CW5 libraries, at
# 0 / 0.5 / 1 / 2 / 5 / 10 %, two replicates each. See common/make_spikein.py.
#
# Per group:
#   prepare        group view + two disjoint halves            (1 job)
#   vireo_full     K=2 discovery on the whole group            (detection)
#   vireo_A/B      K=2 discovery on each half                  (composition gate)
#   souporcell_full K=2 discovery on the whole group           (detection)
# then one report job scores every group into report/detection_curve.tsv.
#
# This runner lives on the shared base but drives both backend adapters, which
# are on their own branches, through their worktrees. Re-running it submits only
# work whose output is missing.
#
#   bash run_spikein.sh

source /etc/profile >/dev/null 2>&1 || true
source ~/.bashrc >/dev/null 2>&1 || true
set -euo pipefail

CORE="${CORE:-/nfs/cellgeni/batu/_repos/nf-demux-preg/validation/group_matrix}"
VIREO_H="${VIREO_H:-/nfs/cellgeni/batu/_repos/worktrees/group-vireo/validation/group_matrix}"
SOC_H="${SOC_H:-/nfs/cellgeni/batu/_repos/worktrees/group-souporcell/validation/group_matrix}"
TICKET="${TICKET:-/lustre/scratch127/cellgen/cellgeni/tickets/REQ-68805}"
GM="$TICKET/group_matrix"
BASE="${BASE:-$GM/spikein}"
CW5="$GM/CW5/build_a/CW5"
QUEUE="${QUEUE:-normal}"
BSUB="${BSUB:-bsub}"          # overridable so the submission path can be dry-run

mkdir -p "$BASE/logs"
cd "$BASE/logs"

jobid() { sed -n 's/^Job <\([0-9]*\)>.*/\1/p'; }

echo "== design =="
if [[ ! -f "$BASE/spikein_design.tsv" ]]; then
  python3 "$CORE/common/make_spikein.py" \
    --truth-a "$GM/gate2_souporcell/joint_seed1/assignments.tsv" \
    --truth-b "$GM/gate2_vireo/joint_seed1/assignments.tsv" \
    --outdir "$BASE"
else
  echo "   reusing $BASE/spikein_design.tsv"
fi

declare -a ALL_IDS=()
while IFS=$'\t' read -r label _rest; do
  [[ "$label" == "label" ]] && continue
  GDIR="$BASE/groups/$label"
  echo "== $label =="

  PREP_ID=""
  if [[ ! -f "$GDIR/B/CW5/alt.mtx" ]]; then
    PREP_ID=$(HARNESS="$CORE" GROUP_DIR="$CW5" GDIR="$GDIR" \
      $BSUB -J "spk_prep_$label" -q "$QUEUE" < "$CORE/lsf/prepare_spikein.bsub" | tee /dev/stderr | jobid)
    ALL_IDS+=("$PREP_ID")
  fi
  WAIT=()
  [[ -n "$PREP_ID" ]] && WAIT=(-w "done($PREP_ID)")

  for part in full A B; do
    OUT="$GDIR/vireo_$part"
    [[ -f "$OUT/assignments.tsv" ]] && continue
    mem=16000; [[ "$part" != "full" ]] && mem=10000
    id=$(HARNESS="$VIREO_H" MODE=discovery GROUP_DIR="$GDIR/$part/CW5" OUTDIR="$OUT" \
      N_DONOR=2 SEED=1 N_INIT=50 THREADS=4 \
      $BSUB -J "spk_vireo_${label}_$part" -q "$QUEUE" -M "$mem" \
        -R "select[mem>=${mem}] rusage[mem=${mem}] span[hosts=1]" "${WAIT[@]}" \
        < "$VIREO_H/lsf/vireo_run.bsub" | tee /dev/stderr | jobid)
    ALL_IDS+=("$id")
  done

  OUT="$GDIR/souporcell_full"
  if [[ ! -f "$OUT/assignments.tsv" ]]; then
    id=$(HARNESS="$SOC_H" MODE=discovery GROUP_DIR="$GDIR/full/CW5" OUTDIR="$OUT" \
      NUM_CLUSTERS=2 SEED=1 THREADS=4 \
      $BSUB -J "spk_soc_${label}" -q "$QUEUE" -M 12000 \
        -R "select[mem>=12000] rusage[mem=12000] span[hosts=1]" "${WAIT[@]}" \
        < "$SOC_H/lsf/souporcell_run.bsub" | tee /dev/stderr | jobid)
    ALL_IDS+=("$id")
  fi
done < "$BASE/spikein_design.tsv"

echo "== report =="
cat >"$BASE/report.bsub" <<EOF
#BSUB -J spk_report
#BSUB -q normal
#BSUB -n 1
#BSUB -M 8000
#BSUB -R "select[mem>=8000] rusage[mem=8000]"
#BSUB -o %J.spk_report.out
#BSUB -e %J.spk_report.err
#BSUB -Ep /software/cellgen/cellgeni/etc/notify-slack.sh
source /etc/profile >/dev/null 2>&1 || true
python3 $CORE/common/score_spikein.py --base $BASE
EOF

REPORT_WAIT=()
if [[ ${#ALL_IDS[@]} -gt 0 ]]; then
  dep=""
  for id in "${ALL_IDS[@]}"; do dep="${dep:+$dep && }ended($id)"; done
  REPORT_WAIT=(-w "$dep")
fi
$BSUB "${REPORT_WAIT[@]}" < "$BASE/report.bsub"

echo
echo "submitted ${#ALL_IDS[@]} jobs plus the report; results under $BASE/report"
