#!/bin/zsh
# PNEUMOTACH (p) binning for the one diffusion session that has a pneumotach file (2023-11-09 000HH; run_adc4_2026-10-21
# ran s,d only). Same set-up as run_adc4: convert_mixed_adc.py ASAP_BSET=b0|b, root reader, v2 gp traj (gas-only),
# spec 6b01071. b-set time axis = index * 2 TR (12-block rule keeps every other scan, jitter <= 1 TR) -> pneumotach
# alignment holds. Output <key>_b0/p, <key>_b/p next to the existing s/ d/. Ext floor: >= 30 GiB free before each run.
cd /Users/hoomham/Hooman/Work/Codes/2026_ASAP_Recon/workspace/pipeline
OUT=/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/AIkill_Dynamic_2022
IMG=/Volumes/HoomHamExt/_5t_images_roundtrip/Images
GP=$PWD/../data/xe/human/traj/fa_spiral_dyn_fancy_v2_20230131_gp.npy
LOG=$OUT/_logs
unset ASAP_DROP_CH
export MPLBACKEND=Agg ASAP_CONVERTER=$PWD/convert_mixed_adc.py ASAP_ROOT_READER=1
ts() { date +%H:%M:%S; }

for b in b0 b; do
  free=$(df -k /Volumes/HoomHamExt | tail -1 | awk '{print int($4/1048576)}')
  if (( free < 32 )); then echo "[$(ts)] STOP: Ext free ${free} GiB (floor 30)" >> $LOG/_adc_p.log; break; fi
  echo "[$(ts)] START 2023-11-09_000HH_$b p (Ext free ${free} GiB)" >> $LOG/_adc_p.log
  ASAP_BSET=$b python3 dyn_recon.py --data-dir $IMG/2023-11-09/000HH --date 2023-11-09 --id 000HH_$b --out-root $OUT \
      --gp-traj $GP --codespec recon_codespec_6b01071.yml --methods p > $LOG/2023-11-09_000HH_${b}_p.log 2>&1
  echo "[$(ts)] END 2023-11-09_000HH_$b exit $?" >> $LOG/_adc_p.log
done
echo "[$(ts)] ADCPDONE" >> $LOG/_adc_p.log
