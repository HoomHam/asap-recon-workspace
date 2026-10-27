#!/bin/zsh
# PREPARED 2026-10-26 for the PCA Step-1 work — NOT RUN. Needs Hooman's go (Tyger compute + ~1 GB/session on Ext).
#
# SIGNAL-binned (self-gated, time-phase bins) recon with the PRODUCTION image, so the result is the b44d image
# (calcb phase low-pass 4.4 vox + deapod + calcb guard + header exclude ranges = fork 262b16f) binned by Steve's
# k0-signal EE minima instead of the diaphragm navigator. Output tree (mirror law, same project folder on Ext):
#     /Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/AIkill_Dynamic_b44d_sig/<date>_<id>/s/{input.mrd,output.mrd,tyger.log,codespec.yml,recon.mat,...}
# dyn_recon.py does convert (raw .dat -> input.mrd with binning=SIGNAL) + tyger submit + post (recon.mat).
#
# If the bin-export fork patch (notes/ForkPatch_BinExport_2026-10-26.md) is applied and a new image built, point
# SPEC at that image's codespec instead; the output then also carries ilvbin / ee_times / per-bin time + nav z.
#
#   GO=1 zsh run_sig_b44d.sh 2024-08-19_015HK 2024-08-30_014KS ...        (one session at a time, ~10-20 min each
#                                                                          wall on the slow upload link)
set -u
HERE=${0:A:h}
SPEC=${SPEC:-$HERE/recon_codespec_262b16f.yml}      # NOT the default recon_codespec.yml (that pins the old d136eb1 image)
PY=/opt/homebrew/Caskroom/miniforge/base/bin/python3   # has mrd + mapvbvd (F15)
R=/Volumes/HoomHamExt/_5t_images_roundtrip/Images
OUT=/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/AIkill_Dynamic_b44d_sig
LOGS=$OUT/_logs
[[ "${GO:-0}" == "1" ]] || { echo "prepared only; set GO=1 after Hooman's go"; exit 2; }
mount | grep -q " /Volumes/HoomHamExt " || { echo "EXT NOT MOUNTED"; exit 3; }
free_gib=$(df -g /Volumes/HoomHamExt | tail -1 | awk '{print $4}')
[[ $free_gib -gt 40 ]] || { echo "Ext free ${free_gib} GiB < 40 GiB floor (System rule: stop before 30)"; exit 4; }
mode=$(python3 -c "import json;print(json.load(open('$HOME/Hooman/System/state/compute_mode.json'))['mode'])")
mkdir -p $LOGS
for key in "$@"; do
  d=${key%%_*}; id=${key#*_}
  log=$LOGS/$key.log
  echo "[$(date +%H:%M:%S)] $key: START (compute_mode=$mode)"
  cmd=($PY $HERE/dyn_recon.py --date $d --id $id --data-root $R --out-root $OUT --methods s --codespec $SPEC)
  if [[ "$mode" == "comfort" ]]; then taskpolicy -b "${cmd[@]}" >> $log 2>&1; else "${cmd[@]}" >> $log 2>&1; fi
  echo "[$(date +%H:%M:%S)] $key: exit $? -> $OUT/$key/s"
done
