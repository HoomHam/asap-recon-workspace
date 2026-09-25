#!/bin/zsh
# PRODUCTION rerun with calcb phase low-pass 4.4 vox (fork 04f445b) from stored d/input.mrd.
# Output tree: Ext AIkill_Dynamic_b44/<key>/d/{output.mrd, tyger.log, codespec.yml}; recon.mat via mrd_to_mat.py.
SPEC=/Users/hoomham/Hooman/Work/Codes/2026_ASAP_Recon/workspace/pipeline/recon_codespec_04f445b.yml
SRC=/Volumes/HoomHamExt/AIkill_Dynamic
DST=/Volumes/HoomHamExt/AIkill_Dynamic_b44
for key in "$@"; do
  mount | grep -q " /Volumes/HoomHamExt " || { echo "[$(date +%H:%M:%S)] $key: EXT NOT MOUNTED, abort"; exit 3; }
  out=$DST/$key/d; mkdir -p $out
  if [ -f $out/output.mrd ] && [ $(stat -f%z $out/output.mrd) -gt 300000000 ] && python3 -c "import mrd,sys
with mrd.BinaryMrdReader('$out/output.mrd') as r:
    r.read_header(); n=sum(1 for _ in r.read_data())" 2>/dev/null; then echo "[$(date +%H:%M:%S)] $key: exists, skip"; continue; fi
  echo "[$(date +%H:%M:%S)] $key: submitting"
  cat $SRC/$key/d/input.mrd | ~/bin/tyger run exec -f $SPEC --logs > $out/output.mrd 2> $out/tyger.log
  echo "[$(date +%H:%M:%S)] $key: exit $? size $(du -h $out/output.mrd | cut -f1)"
  cp $SPEC $out/codespec.yml
done
