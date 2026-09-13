#!/bin/zsh
# Reconstruct every v2/v3 dynamic session that had no output in AIkill_Dynamic
# (inventory 2026-09-12). One subject at a time, all methods (s,p,d); a failure
# is logged and the loop moves on. Detach with nohup — each subject can take
# an hour+ on the slow link.
#
#   nohup zsh todo_batch.sh <codespec.yml> [worklist-filter-regex] > <log> 2>&1 &
#
# Per-subject logs + SUMMARY.txt -> $LOGS. Re-running is safe: dyn_recon skips
# methods whose output.mrd already exists. FORCE=1 in the environment re-runs
# them anyway (e.g. after an image change):
#   FORCE=1 nohup zsh todo_batch.sh <spec> '011CN|025VP|007IT' > <log> 2>&1 &
set -u
SPEC=${1:?codespec yml}
FILTER=${2:-.}
HERE=${0:A:h}
PY=/opt/homebrew/Caskroom/miniforge/base/bin/python3   # has mrd + mapvbvd (F15)
R=/Volumes/HoomHamExt/_5t_images_roundtrip/Images
S=/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/staged_src   # symlinked loose-layout sessions
OUT=/Volumes/HoomHamExt/AIkill_Dynamic
LOGS=/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/pipeline_runs/logs
mkdir -p $LOGS

# date id data-root   (small / single-coil first, 8-coil 1-2 GB last)
WORKLIST=(
  "2024-03-06 030DN $R"
  "2024-04-29 011CN $R"
  "2023-11-02 025VP $R"
  "2024-08-28 007IT $R"
  "2024-01-17 042DR $R"
  "2023-07-24 019WR $S"
  "2023-01-31 000LL $R"
  "2023-10-27 000LL $R"
  "2023-04-11 000KR $R"
  "2023-04-14 000KR $R"
  "2023-05-15 000SK $R"
  "2023-05-18 000SV $R"
  "2023-11-03 000LL $R"
  # "2023-11-02 000LL $S"   # bigmac copy = byte-identical to 2023-11-03/000LL; output quarantined 2026-09-13
  "2024-06-18 000LL $R"
  "2024-09-09 000LL $R"
  "2024-09-20 000LL $R"
  "2024-02-02 005DS $R"
  "2024-02-02 000LL $R"
  "2023-03-08 017AK $R"
  "2023-03-09 023DB $R"
  "2023-03-10 004LR $R"
  "2023-03-24 029CK $R"
  "2023-03-27 02BB $R"
  "2023-03-31 032WS $R"
)

for line in $WORKLIST; do
  parts=(${=line}); d=$parts[1]; id=$parts[2]; root=$parts[3]
  [[ "${d}_${id}" =~ $FILTER ]] || continue
  log=$LOGS/${d}_${id}.log
  echo "[batch $(date +%H:%M:%S)] START ${d}_${id}"
  $PY $HERE/dyn_recon.py --date $d --id $id --data-root $root --out-root $OUT \
      --methods s,p,d --codespec $SPEC ${FORCE:+--force} >> $log 2>&1
  rc=$?
  done_m=""
  for m in s p d; do [[ -s $OUT/${d}_${id}/$m/output.mrd ]] && done_m+=$m || done_m+=-; done
  echo "$(date +%F_%T) ${d}_${id} rc=$rc outputs=$done_m" >> $LOGS/SUMMARY.txt
  echo "[batch $(date +%H:%M:%S)] END   ${d}_${id} rc=$rc outputs=$done_m"
done
echo "[batch $(date +%H:%M:%S)] ALL DONE"
