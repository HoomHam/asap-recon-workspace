#!/bin/zsh
# Early-2022 dynamics the scanner split into 100-rep chunks (MDH gaps 0.1-5 s, every chunk restarts at Lin 0 / Rep 0):
# one merged input per series via dyn_recon --merge (merge_dyn trims every non-last chunk to whole 480-arm cycles,
# fix 2026-10-21). Replaces the single-chunk recons of run_early2022 for these keys. Same calibrations as there.
# Plus 2022-02-03 000ID MID00142 (the other of two separate 98 s scans; MID00143 already done) via a staged symlink.
cd /Users/hoomham/Hooman/Work/Codes/2026_ASAP_Recon/workspace/pipeline
W=/Volumes/WorkVault/Work/_inbox/_from_BigMac_HD4_2026_09_14/Images/MRI/Human
O=/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/AIkill_Dynamic_2022
T=$PWD/../data/xe/human/traj
E1=$T/calib_faraz_entry1_20220202_axial_gp.npy; E2=$T/calib_faraz_entry2_20220215_coronal_gp.npy; L24=$T/calib_local_20220224_axial_gp.npy
LOG=$O/_logs
unset ASAP_DROP_CH ASAP_CONVERTER ASAP_ROOT_READER ASAP_BSET ASAP_GASPHASE
export MPLBACKEND=Agg
ts() { date +%H:%M:%S; }

mrg() {  # date key dir traj mids
  echo "[$(ts)] START ${1}_$2 merge $5 ($(basename $4))" >> $LOG/_early_merge.log
  python3 dyn_recon.py --data-dir $3 --date $1 --id $2 --out-root $O --gp-traj $4 --merge $5 \
      --codespec recon_codespec_6b01071.yml --methods s,d > $LOG/${1}_${2}_merged.log 2>&1
  echo "[$(ts)] END ${1}_$2 exit $?" >> $LOG/_early_merge.log
}
mids() { local p=$1 a=$2 b=$3 out=() i; for i in $(seq $a $b); do out+=$(printf 'MID%05d' $i); done; echo ${(j:,:)out}; }

mrg 2022-03-10 006RL_A $W/2022-03-10_006RL $L24 $(mids x 1264 1268)
mrg 2022-03-10 006RL_B $W/2022-03-10_006RL $L24 $(mids x 1289 1293)
mrg 2022-02-25 054SK $W/2022-02-25_054SK $L24 $(mids x 410 414)
mrg 2022-03-09 015MG $W/2022-03-09_015MG $E2 $(mids x 1072 1075)
mrg 2022-02-09 unk-MID38 $W/2022-02-09 $E1 $(mids x 38 42)
mrg 2022-02-09 unk-MID74 $W/2022-02-09_2 $E1 $(mids x 74 78)
mrg 2022-02-15 unk-MID159 $W/2022-02-15_2 $E2 $(mids x 159 163)
mrg 2022-02-02 unk-MID16 $W/2022-02-02_2 $E1 $(mids x 16 24)
echo "[$(ts)] START 2022-02-03_000ID_MID142" >> $LOG/_early_merge.log
python3 dyn_recon.py --data-dir /Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/staged_src/2022-02-03/000ID_MID142 --date 2022-02-03 \
    --id 000ID_MID142 --out-root $O --gp-traj $E1 --codespec recon_codespec_6b01071.yml --methods s,d > $LOG/2022-02-03_000ID_MID142.log 2>&1
echo "[$(ts)] END 2022-02-03_000ID_MID142 exit $?" >> $LOG/_early_merge.log
echo "[$(ts)] EARLYMERGEDONE" >> $LOG/_early_merge.log
