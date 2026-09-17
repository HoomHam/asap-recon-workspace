#!/bin/zsh
# Re-recon the four sessions that have a GOOD second free-breathing dynamic, with both
# dynamics merged into one input (Hooman's decision 2026-09-16, relayed by the XeCS
# session; list from XeCS outputs/calspec/k0_sessions.csv `extra_free`).
# Methods s + d only (p is unsupported on merged input, see merge_dyn.py).
# Detach:  nohup zsh workspace/pipeline/merged_batch.sh > /dev/null 2>&1 &
# Logs:    /Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/pipeline_runs/logs/<date>_<id>_merged.log
set -u
cd "$(dirname "$0")/.."                                    # workspace/
PY=/opt/homebrew/Caskroom/miniforge/base/bin/python3       # has mrd + mapvbvd (F15)
SPEC=${SPEC:-pipeline/recon_codespec_40d23a4.yml}
LOGS=/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/pipeline_runs/logs
mkdir -p "$LOGS"

# date id mids(acquisition order) ref(breath-hold MID or -)
WORK=(
  "2024-03-06 030DN MID00543,MID00545 MID00515"
  "2024-08-12 013VM MID00408,MID00410 -"
  "2023-11-02 025VP MID02741,MID02743 -"
  "2024-01-31 008CR MID00194,MID00199 MID00157"
)

for line in "${WORK[@]}"; do
  set -- ${=line}
  date=$1; id=$2; mids=$3; ref=$4
  log="$LOGS/${date}_${id}_merged.log"
  args=(--date "$date" --id "$id" --merge "$mids" --methods s,d --codespec "$SPEC")
  [[ "$ref" != "-" ]] && args+=(--ref "$ref")
  echo "=== $(date '+%F %T') START ${date}_${id} merge=$mids ref=$ref" >> "$log"
  "$PY" pipeline/dyn_recon.py "${args[@]}" >> "$log" 2>&1
  rc=$?                       # capture BEFORE the $(date) substitution below resets $?
  echo "=== $(date '+%F %T') EXIT $rc ${date}_${id}" >> "$log"
done
echo "=== $(date '+%F %T') ALL DONE" >> "$LOGS/merged_batch_SUMMARY.txt"
