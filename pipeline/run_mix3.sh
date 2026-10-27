#!/bin/zsh
# 2026-10-27 — Hooman's vote: apply MIX binning ONCE to three subjects and study it against s/d/p ("mix-binning").
# Image = fork 36f67fb (508b329 + int-params fix) (MIX mode + binsigma/binpou knobs + bin export). Runs m (MIX, sigma 0.5, partition of unity)
# and p (PNEUMOTACH, same image, default kernel) for each session; s and d exist on the same image lineage in
# AIkill_Dynamic_b44d_sig/<key>/s and AIkill_Dynamic_b44d/<key>/d (262b16f == 508b329 for those modes: defaults
# reproduce the old kernel bit for bit).
# Output (mirror law): /Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/AIkill_Dynamic_mix/<key>/{m,p}/
#   zsh run_mix3.sh                      (3 sessions, sequential, ~3 min each per mode on Tyger)
set -u
HERE=${0:A:h}
SPEC=${SPEC:-$HERE/recon_codespec_40c248b.yml}
PY=/opt/homebrew/Caskroom/miniforge/base/bin/python3
R=/Volumes/HoomHamExt/_5t_images_roundtrip/Images
OUT=/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/AIkill_Dynamic_mix
LOGS=$OUT/_logs
KEYS=("$@")
MODES=(${=MODES:-m p})          # MODES="m" or "p" to run one mode only
[[ ${#KEYS[@]} -eq 0 ]] && KEYS=(2024-08-19_015HK 2024-08-23_016PG 2024-08-30_014KS)
[[ -f $SPEC ]] || { echo "codespec $SPEC missing (image not built yet?)"; exit 2; }
mount | grep -q " /Volumes/HoomHamExt " || { echo "EXT NOT MOUNTED"; exit 3; }
free_gib=$(df -g /Volumes/HoomHamExt | tail -1 | awk '{print $4}')
[[ $free_gib -gt 40 ]] || { echo "Ext free ${free_gib} GiB < 40 GiB floor"; exit 4; }
mode=$(python3 -c "import json;print(json.load(open('$HOME/Hooman/System/state/compute_mode.json'))['mode'])")
mkdir -p $LOGS
for key in $KEYS; do
  d=${key%%_*}; id=${key#*_}
  for m in $MODES; do
    log=$LOGS/${key}_$m.log
    echo "[$(date +%H:%M:%S)] $key/$m: START (compute_mode=$mode)"
    cmd=($PY $HERE/dyn_recon.py --date $d --id $id --data-root $R --out-root $OUT --methods $m --codespec $SPEC)
    [[ $m == m ]] && cmd+=(--binsigma 0.5 --binpou 1)
    [[ $m == p ]] && cmd+=(--pneumoclock)
    if [[ "$mode" == "comfort" ]]; then taskpolicy -b "${cmd[@]}" >> $log 2>&1; else "${cmd[@]}" >> $log 2>&1; fi
    echo "[$(date +%H:%M:%S)] $key/$m: exit $? -> $OUT/$key/$m"
  done
done
