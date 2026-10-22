#!/bin/zsh
# 2022 saturation family (fa_spiral_dyn_fancy_20220602): gas imaging every line + DP_Gauss saturation pulses, no DP
# readout -> gas-only ASAP (arm = counter over all lines mod 640, F117). Trajectory = entry 4 = v2 gp npy
# (seqnames.txt + gp symlink 2026-10-21; no dp symlink on purpose). Spec 6b01071, method s,d, all coils, root reader.
# Waits for run_failed4 (one Tyger job at a time).
cd /Users/hoomham/Hooman/Work/Codes/2026_ASAP_Recon/workspace/pipeline
OUT=/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/AIkill_Dynamic_2022
IMG=/Volumes/HoomHamExt/_5t_images_roundtrip/Images
LOG=$OUT/_logs
unset ASAP_DROP_CH ASAP_CONVERTER
export MPLBACKEND=Agg
ts() { date +%H:%M:%S; }
until grep -q FAILED4DONE /Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/AIkill_Dynamic_b44d/_logs/_failed4.log 2>/dev/null; do sleep 30; done

run() {  # date id
  echo "[$(ts)] START $1_$2 s,d" >> $LOG/_sat2.log
  python3 dyn_recon.py --data-dir $IMG/$1/$2 --date $1 --id $2 --out-root $OUT \
      --codespec recon_codespec_6b01071.yml --methods s,d > $LOG/$1_$2.log 2>&1
  echo "[$(ts)] END $1_$2 exit $?" >> $LOG/_sat2.log
}

run 2022-06-28 019WR
run 2022-07-05 021JM
echo "[$(ts)] SAT2DONE" >> $LOG/_sat2.log
