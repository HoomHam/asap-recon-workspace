#!/bin/zsh
# Two-subject Tyger check of fork ee3c91f (deapod) from stored d/input.mrd; compared against the offline deapod tree.
SPEC=/Users/hoomham/Hooman/Work/Codes/2026_ASAP_Recon/workspace/pipeline/recon_codespec_ee3c91f.yml
SRC=/Volumes/HoomHamExt/AIkill_Dynamic
DST=/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/tyger_deapod_check_2026-09-28
# wait for the GitHub image build of ee3c91f
while true; do
  st=$(gh run list --repo HoomHam/asap_recon --limit 5 --json headSha,status,conclusion --jq '.[] | select(.headSha | startswith("ee3c91f")) | .status + ":" + (.conclusion // "")')
  [[ $st == completed:success ]] && { echo "[$(date +%H:%M:%S)] image build OK"; break; }
  [[ $st == completed:* ]] && { echo "[$(date +%H:%M:%S)] image build FAILED ($st)"; exit 2; }
  sleep 20
done
for key in "$@"; do
  mount | grep -q " /Volumes/HoomHamExt " || { echo "[$(date +%H:%M:%S)] $key: EXT NOT MOUNTED, abort"; exit 3; }
  out=$DST/$key/d; mkdir -p $out
  echo "[$(date +%H:%M:%S)] $key: submitting"
  cat $SRC/$key/d/input.mrd | ~/bin/tyger run exec -f $SPEC --logs > $out/output.mrd 2> $out/tyger.log
  echo "[$(date +%H:%M:%S)] $key: exit $? size $(du -h $out/output.mrd | cut -f1)"
  cp $SPEC $out/codespec.yml
  /opt/homebrew/Caskroom/miniforge/base/bin/python mrd_to_mat.py $DST/$key
done
echo "[$(date +%H:%M:%S)] CHECK RUNS DONE"
