#!/bin/zsh
# 2022 'ADC' (diffusion) dynamics: one Tyger run per b-set (convert_mixed_adc.py ASAP_BSET=b0|b, 12-scan block rule
# F117, TR doubled to the in-set spacing; root reader). Trajectory = entry 4 = v2 gp npy, gas-only (seqnames +
# 0421/0422 gp symlinks 2026-10-21; 2023-11-09 runs on seq 20221103 -> --gp-traj forces gas-only there too).
# Spec 6b01071, methods s,d, all coils. Output <key>_b0 / <key>_b. Waits for run_sat2 (one Tyger job at a time).
cd /Users/hoomham/Hooman/Work/Codes/2026_ASAP_Recon/workspace/pipeline
OUT=/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/AIkill_Dynamic_2022
IMG=/Volumes/HoomHamExt/_5t_images_roundtrip/Images
GP=$PWD/../data/xe/human/traj/fa_spiral_dyn_fancy_v2_20230131_gp.npy
LOG=$OUT/_logs
unset ASAP_DROP_CH
export MPLBACKEND=Agg ASAP_CONVERTER=$PWD/convert_mixed_adc.py ASAP_ROOT_READER=1
ts() { date +%H:%M:%S; }
until grep -q SAT2DONE $LOG/_sat2.log 2>/dev/null; do sleep 30; done

run() {  # date id bset
  echo "[$(ts)] START $1_$2_$3 s,d" >> $LOG/_adc4.log
  ASAP_BSET=$3 python3 dyn_recon.py --data-dir $IMG/$1/$2 --date $1 --id $2_$3 --out-root $OUT --gp-traj $GP \
      --codespec recon_codespec_6b01071.yml --methods s,d > $LOG/$1_$2_$3.log 2>&1
  echo "[$(ts)] END $1_$2_$3 exit $?" >> $LOG/_adc4.log
}

for b in b0 b; do
  run 2023-11-09 000HH $b
  run 2022-04-22 000HH $b
  run 2022-04-21 000HH $b
  run 2022-04-29 018BI $b
done
echo "[$(ts)] ADC4DONE" >> $LOG/_adc4.log
