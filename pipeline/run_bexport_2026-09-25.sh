#!/bin/zsh
# Re-run DIAPHRAGM recons with the fork image f5ac7c5 (complex gas export) from the stored d/input.mrd.
# Usage: run_gascplx_2026-09-25.sh <key> [<key> ...]   e.g. 2024-09-10_045VS
SPEC=/Users/hoomham/Hooman/Work/Codes/2026_ASAP_Recon/workspace/pipeline/recon_codespec_c8366c3.yml
SRC=/Volumes/HoomHamExt/AIkill_Dynamic
DST=/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/tyger_bexport_2026-09-25
for key in "$@"; do
  # MOUNT GUARD: never write into an unmounted /Volumes path
  mount | grep -q " /Volumes/HoomHamExt " || { echo "[$(date +%H:%M:%S)] $key: EXT NOT MOUNTED, abort"; exit 3; }
  out=$DST/$key/d; mkdir -p $out
  if [ -f $out/output.mrd ] && [ $(stat -f%z $out/output.mrd) -gt 300000000 ]; then echo "[$(date +%H:%M:%S)] $key: exists, skip"; continue; fi
  echo "[$(date +%H:%M:%S)] $key: submitting $SRC/$key/d/input.mrd"
  cat $SRC/$key/d/input.mrd | ~/bin/tyger run exec -f $SPEC --logs > $out/output.mrd 2> $out/tyger.log
  echo "[$(date +%H:%M:%S)] $key: exit $? size $(du -h $out/output.mrd | cut -f1)"
  cp $SPEC $out/codespec.yml
done
