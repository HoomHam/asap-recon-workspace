#!/bin/zsh
# Retry of the sessions that failed calcb before the fork 6b01071 guard (flat-phase fallback, empty-mask safe).
# 2023-07-24 019WR (v2, 1 channel, 2026-09-12 run 573: calcb `ydata` must not be empty). Methods s,d (no pneumotach).
# Production tree (b44 + deapod; 6b01071 = ee3c91f + guard), root converter, ASAP_DROP_CH not set.
cd /Users/hoomham/Hooman/Work/Codes/2026_ASAP_Recon/workspace/pipeline
OUT=/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/AIkill_Dynamic_b44d
IMG=/Volumes/HoomHamExt/_5t_images_roundtrip/Images
LOG=$OUT/_logs
unset ASAP_DROP_CH ASAP_CONVERTER
export MPLBACKEND=Agg
ts() { date +%H:%M:%S; }

run() {  # date id datadir methods
  echo "[$(ts)] START $1_$2 $4" >> $LOG/_failed4.log
  python3 dyn_recon.py --data-dir $3 --date $1 --id $2 --out-root $OUT \
      --codespec recon_codespec_6b01071.yml --methods $4 > $LOG/$1_$2.log 2>&1
  echo "[$(ts)] END $1_$2 exit $?" >> $LOG/_failed4.log
}

run 2023-07-24 019WR $IMG/2023-07-24 s,d
echo "[$(ts)] FAILED4DONE" >> $LOG/_failed4.log
