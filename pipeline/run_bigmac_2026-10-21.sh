#!/bin/zsh
# BigMac-only 2022/v2 sessions, read IN PLACE from WorkVault (read-only; Estate custody 2026-10-21: no writes, no copy;
# run dirs + outputs on HoomHamExt). Trajectory entry 4/10 = v2 npy via seqnames + symlinks (2026-10-21).
# Family per session from the WorkVault header/k0 probe (scratchpad wv_probe.txt): std = root converter (pattern
# auto-detected by raw.py); bset = 'ADC' b0/b split (convert_mixed_adc ASAP_BSET, rule auto) -> <ID>_b0 / <ID>_b;
# gas3 = GP-DP-DP with gas-laden DP scans, gas phase only (ASAP_GASPHASE=3); acq7 = mixed-ADC reader (1024 spectra).
# Spec 6b01071 (= b44 + deapod + calcb guard), all coils, methods s,p,d (p only with a pneumotach file).
cd /Users/hoomham/Hooman/Work/Codes/2026_ASAP_Recon/workspace/pipeline
O22=/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/AIkill_Dynamic_2022
Ob44d=/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/AIkill_Dynamic_b44d
GP=$PWD/../data/xe/human/traj/fa_spiral_dyn_fancy_v2_20230131_gp.npy
CONV=$PWD/convert_mixed_adc.py
LOG=$O22/_logs
unset ASAP_DROP_CH ASAP_CONVERTER ASAP_ROOT_READER ASAP_BSET ASAP_GASPHASE
export MPLBACKEND=Agg
ts() { date +%H:%M:%S; }

one() {  # date id tree dir key [env...] -- [extra args]
  local d=$1 i=$2 tree=$3 dir=$4 key=$5; shift 5
  local out=$O22; [[ $tree == b44d ]] && out=$Ob44d
  local envs=() extra=()
  while [[ $# -gt 0 && $1 != -- ]]; do envs+=$1; shift; done; [[ $1 == -- ]] && shift; extra=("$@")
  if [[ ! -d $dir ]]; then echo "[$(ts)] SKIP ${d}_$key: WorkVault dir missing (unplugged?)" >> $LOG/_bigmac.log; return; fi
  echo "[$(ts)] START ${d}_$key ($envs)" >> $LOG/_bigmac.log
  env $envs python3 dyn_recon.py --data-dir $dir --date $d --id $key --out-root $out \
      --codespec recon_codespec_6b01071.yml $extra > $LOG/${d}_$key.log 2>&1
  echo "[$(ts)] END ${d}_$key exit $?" >> $LOG/_bigmac.log
}

run() {  # date id mode tree dir
  case $3 in
    std)  one $1 $2 $4 $5 $2 ;;
    acq7) one $1 $2 $4 $5 $2 ASAP_CONVERTER=$CONV ;;
    gas3) one $1 $2 $4 $5 $2 ASAP_CONVERTER=$CONV ASAP_ROOT_READER=1 ASAP_GASPHASE=3 -- --gp-traj $GP ;;
    bset) for b in b0 b; do one $1 $2 $4 $5 ${2}_$b ASAP_CONVERTER=$CONV ASAP_ROOT_READER=1 ASAP_BSET=$b -- --gp-traj $GP; done ;;
  esac
}

run 2022-03-25 016ID std 22 "/Volumes/WorkVault/Work/_inbox/_from_BigMac_HD4_2026_09_14/Images/MRI/Human/2022-03-25_016ID"
run 2022-03-30 020WC std 22 "/Volumes/WorkVault/Work/_inbox/_from_BigMac_HD4_2026_09_14/Images/MRI/Human/2022-03-30_020WC"
run 2022-04-13 000LL std 22 "/Volumes/WorkVault/Work/_inbox/_from_BigMac_HD4_2026_09_14/Images/MRI/Human/2022-04-13_000LL"
run 2022-04-19 000HH std 22 "/Volumes/WorkVault/Work/_inbox/_from_BigMac_HD4_2026_09_14/Images/MRI/Human/2022-04-19_000HH"
run 2022-04-20 017AK bset 22 "/Volumes/WorkVault/Work/_inbox/_from_BigMac_HD4_2026_09_14/Images/MRI/Human/2022-04-20_017AK"
run 2022-04-25 000LL bset 22 "/Volumes/WorkVault/Work/_inbox/_from_BigMac_HD4_2026_09_14/Images/MRI/Human/2022-04-25_000LL"
run 2022-04-26 000ID bset 22 "/Volumes/WorkVault/Work/_inbox/_from_BigMac_HD4_2026_09_14/Images/MRI/Human/2022-04-26_000ID"
run 2022-05-02 020ST bset 22 "/Volumes/WorkVault/Work/_inbox/_from_BigMac_HD4_2026_09_14/Images/MRI/Human/2022-05-02_020ST"
run 2022-06-02 008LL std 22 "/Volumes/WorkVault/Work/_inbox/_from_BigMac_HD4_2026_09_14/Images/MRI/Human/2022-06-02-008LL"
run 2022-06-27 000LL std 22 "/Volumes/WorkVault/Work/_inbox/_from_BigMac_HD4_2026_09_14/Images/MRI/Human/2022-06-27_000LL"
run 2022-07-06 000HH std 22 "/Volumes/WorkVault/Work/_inbox/_from_BigMac_HD4_2026_09_14/Images/MRI/Human/2022-07-06-000HH"
run 2022-07-07 001SW std 22 "/Volumes/WorkVault/Work/_inbox/_from_BigMac_HD4_2026_09_14/Images/MRI/Human/2022-07-07-001SW"
run 2022-07-08 000HH std 22 "/Volumes/WorkVault/Work/_inbox/_from_BigMac_HD4_2026_09_14/Images/MRI/Human/2022-07-08-000HH"
run 2022-07-18 025VP std 22 "/Volumes/WorkVault/Work/_inbox/_from_BigMac_HD4_2026_09_14/Images/MRI/Human/2022-07-18_025VP"
run 2022-07-20 023DB std 22 "/Volumes/WorkVault/Work/_inbox/_from_BigMac_HD4_2026_09_14/Images/MRI/Human/2022-07-20-023DB"
run 2022-07-21 014JW std 22 "/Volumes/WorkVault/Work/_inbox/_from_BigMac_HD4_2026_09_14/Images/MRI/Human/2022-07-21-014JW"
run 2022-08-12 000LL acq7 22 "/Volumes/WorkVault/Work/_inbox/_from_BigMac_HD4_2026_09_14/Images/MRI/Human/2022-08-12-000LL"
run 2022-10-31 000LL std 22 "/Volumes/WorkVault/Work/_inbox/_from_BigMac_HD4_2026_09_14/Images/MRI/Human/2022-10-31_000LL"
run 2022-11-04 000HH gas3 22 "/Volumes/WorkVault/Work/_inbox/_from_BigMac_HD4_2026_09_14/Images/MRI/Human/2022-11-04-000HH"
run 2022-11-07 000HH std 22 "/Volumes/WorkVault/Work/_inbox/_from_BigMac_HD4_2026_09_14/Images/MRI/Human/2022-11-07_000HH"
run 2023-01-24 000LL std 22 "/Volumes/WorkVault/Work/_inbox/_from_BigMac_HD4_2026_09_14/Images/MRI/Human/2023-01-24_000LL"
run 2023-02-17 000LL std b44d "/Volumes/WorkVault/Work/_inbox/_from_BigMac_HD4_2026_09_14/Images/MRI/Human/2023-02-17_000LL"
run 2023-03-02 000LL std b44d "/Volumes/WorkVault/Work/_inbox/_from_BigMac_HD4_2026_09_14/Images/MRI/Human/2023-03-02_000LL"
echo "[$(ts)] BIGMACDONE" >> $LOG/_bigmac.log
