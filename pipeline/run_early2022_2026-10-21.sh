#!/bin/zsh
# Early-2022 (pre-entry-4) BigMac sessions, read in place from WorkVault (read-only). Calibrations from Faraz's library
# calibrations_3D_20220308 (entry 1 2022/02/02 axial, entry 2 2022/02/15 coronal; 32 leaves x 15 reps = 480 arms,
# FOV 350, dwell 5 us = the scans) and the local axial repeat traj_2022-02-24.mat (054SK folder; = entry 1 within 0.9 %).
# Test 2026-10-21: 027CB 02-04 with entry 1 matches Faraz's dynvent anatomy. Gas-only (no DP in these set-ups);
# 2022-02-24 000HH = ADC with plain b0/b alternation -> ASAP_BSET. IDs unknown on 02-02_2 / 02-09 / 02-15 -> 'unk-MID<first>'.
# Skipped: 2022-01-28 000HH (1H), 2022-02-01 (FOV 415, no calib), 2022-02-08 000HH + 02-09 dyn_20220208 (11 leaves, no calib).
cd /Users/hoomham/Hooman/Work/Codes/2026_ASAP_Recon/workspace/pipeline
W=/Volumes/WorkVault/Work/_inbox/_from_BigMac_HD4_2026_09_14/Images/MRI/Human
O=/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/AIkill_Dynamic_2022
T=$PWD/../data/xe/human/traj
E1=$T/calib_faraz_entry1_20220202_axial_gp.npy; E2=$T/calib_faraz_entry2_20220215_coronal_gp.npy; L24=$T/calib_local_20220224_axial_gp.npy
CONV=$PWD/convert_mixed_adc.py
LOG=$O/_logs
unset ASAP_DROP_CH ASAP_CONVERTER ASAP_ROOT_READER ASAP_BSET ASAP_GASPHASE
export MPLBACKEND=Agg
ts() { date +%H:%M:%S; }

one() {  # date key dir traj [env...]
  local d=$1 key=$2 dir=$3 tr=$4; shift 4
  if [[ ! -d $dir ]]; then echo "[$(ts)] SKIP ${d}_$key: dir missing" >> $LOG/_early.log; return; fi
  echo "[$(ts)] START ${d}_$key ($(basename $tr) $*)" >> $LOG/_early.log
  env "$@" python3 dyn_recon.py --data-dir $dir --date $d --id $key --out-root $O --gp-traj $tr \
      --codespec recon_codespec_6b01071.yml --methods s,d > $LOG/${d}_$key.log 2>&1
  echo "[$(ts)] END ${d}_$key exit $?" >> $LOG/_early.log
}

one 2022-02-04 027CB $W/2022-02-04_027CB $E1
one 2022-02-02 AW $W/2022-02-02 $E1
one 2022-02-03 000ID $W/2022-02-03_000ID $E1
one 2022-02-25 054SK $W/2022-02-25_054SK $L24
one 2022-03-10 006RL $W/2022-03-10_006RL $L24
one 2022-03-09 015MG $W/2022-03-09_015MG $E2
one 2022-02-24 000HH_b0 $W/2022-02-24_000HH $L24 ASAP_CONVERTER=$CONV ASAP_ROOT_READER=1 ASAP_BSET=b0
one 2022-02-24 000HH_b $W/2022-02-24_000HH $L24 ASAP_CONVERTER=$CONV ASAP_ROOT_READER=1 ASAP_BSET=b
one 2022-02-02 unk-MID14 $W/2022-02-02_2 $E1
one 2022-02-09 unk-MID38 $W/2022-02-09 $E1
one 2022-02-09 unk-MID74 $W/2022-02-09_2 $E1
one 2022-02-15 unk-MID159 $W/2022-02-15_2 $E2
echo "[$(ts)] EARLYDONE" >> $LOG/_early.log
