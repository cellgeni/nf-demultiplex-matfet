#!/usr/bin/env bash
# Follow-up to the detection-limit experiment, on the same 12 spike-in groups.
# Reuses the existing group views, so nothing is re-subsetted.
#
#   hybrid        Vireo fixed-profile projection of the whole group against
#                 Souporcell's discovery profiles, --genoTag GT. Souporcell
#                 discovers down to 0.5% but cannot project; Vireo projects but
#                 only discovers at >=5%. This tests whether the combination
#                 keeps Souporcell's sensitivity. (12 jobs)
#   souporcell_A/B Souporcell discovery on each half, so the composition gate
#                 covers the discovery stage the hybrid would use. (24 jobs)
#   seeds         Souporcell on the 5% groups with seeds 2 and 3: its recall
#                 dipped to 0.64-0.69 there, with fetal cells called doublets,
#                 while 2% and 10% were near perfect. A seed-dependent dip points
#                 at a local optimum rather than a property of 5%. (4 jobs)
#
# Then one scoring job rewrites report/detection_curve.tsv with every run.
# Re-running submits only work whose output is missing.
#
#   bash run_spikein_followup.sh

source /etc/profile >/dev/null 2>&1 || true
source ~/.bashrc >/dev/null 2>&1 || true
set -euo pipefail

CORE="${CORE:-/nfs/cellgeni/batu/_repos/nf-demux-preg/validation/group_matrix}"
VIREO_H="${VIREO_H:-/nfs/cellgeni/batu/_repos/worktrees/group-vireo/validation/group_matrix}"
SOC_H="${SOC_H:-/nfs/cellgeni/batu/_repos/worktrees/group-souporcell/validation/group_matrix}"
TICKET="${TICKET:-/lustre/scratch127/cellgen/cellgeni/tickets/REQ-68805}"
BASE="${BASE:-$TICKET/group_matrix/spikein}"
QUEUE="${QUEUE:-normal}"
BSUB="${BSUB:-bsub}"

[[ -f "$BASE/spikein_design.tsv" ]] || { echo "ERROR: run run_spikein.sh first" >&2; exit 1; }
mkdir -p "$BASE/logs"
cd "$BASE/logs"

jobid() { sed -n 's/^Job <\([0-9]*\)>.*/\1/p'; }
declare -a ALL_IDS=()

while IFS=$'\t' read -r label _rest; do
  [[ "$label" == "label" ]] && continue
  GDIR="$BASE/groups/$label"
  echo "== $label =="
  for need in "$GDIR/full/CW5/alt.mtx" "$GDIR/souporcell_full/profile_map.json"; do
    [[ -f "$need" ]] || { echo "ERROR: $need missing; the first experiment must finish" >&2; exit 1; }
  done

  if [[ ! -f "$GDIR/hybrid_full/assignments.tsv" ]]; then
    id=$(HARNESS="$VIREO_H" MODE=projection GROUP_DIR="$GDIR/full/CW5" OUTDIR="$GDIR/hybrid_full" \
      DONOR_VCF="$GDIR/souporcell_full/consensus/cluster_genotypes.vcf" \
      PROFILE_MAP="$GDIR/souporcell_full/profile_map.json" GENO_TAG=GT THREADS=4 \
      $BSUB -J "spk_hybrid_$label" -q "$QUEUE" -M 10000 \
        -R "select[mem>=10000] rusage[mem=10000] span[hosts=1]" \
        < "$VIREO_H/lsf/vireo_run.bsub" | tee /dev/stderr | jobid)
    ALL_IDS+=("$id")
  fi

  for half in A B; do
    [[ -f "$GDIR/souporcell_$half/assignments.tsv" ]] && continue
    id=$(HARNESS="$SOC_H" MODE=discovery GROUP_DIR="$GDIR/$half/CW5" OUTDIR="$GDIR/souporcell_$half" \
      NUM_CLUSTERS=2 SEED=1 THREADS=4 \
      $BSUB -J "spk_soc_${label}_$half" -q "$QUEUE" -M 8000 \
        -R "select[mem>=8000] rusage[mem=8000] span[hosts=1]" \
        < "$SOC_H/lsf/souporcell_run.bsub" | tee /dev/stderr | jobid)
    ALL_IDS+=("$id")
  done

  if [[ "$label" == f5_* ]]; then
    for seed in 2 3; do
      [[ -f "$GDIR/souporcell_full_s$seed/assignments.tsv" ]] && continue
      id=$(HARNESS="$SOC_H" MODE=discovery GROUP_DIR="$GDIR/full/CW5" OUTDIR="$GDIR/souporcell_full_s$seed" \
        NUM_CLUSTERS=2 SEED=$seed THREADS=4 \
        $BSUB -J "spk_soc_${label}_s$seed" -q "$QUEUE" -M 12000 \
          -R "select[mem>=12000] rusage[mem=12000] span[hosts=1]" \
          < "$SOC_H/lsf/souporcell_run.bsub" | tee /dev/stderr | jobid)
      ALL_IDS+=("$id")
    done
  fi
done < "$BASE/spikein_design.tsv"

echo "== report =="
cat >"$BASE/report_followup.bsub" <<EOF
#BSUB -J spk_report_followup
#BSUB -q normal
#BSUB -n 1
#BSUB -M 8000
#BSUB -R "select[mem>=8000] rusage[mem=8000]"
#BSUB -o %J.spk_report_followup.out
#BSUB -e %J.spk_report_followup.err
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
$BSUB "${REPORT_WAIT[@]}" < "$BASE/report_followup.bsub"

echo
echo "submitted ${#ALL_IDS[@]} jobs plus the report; results under $BASE/report"
