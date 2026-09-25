#!/bin/zsh
# Re-run the atlas sessions with the raw-b export image (c8366c3), 3 in flight, from stored d/input.mrd.
# Outputs: Ext Work/Codes/2026_ASAP_Recon/tyger_bexport_2026-09-25/<key>/d/{output.mrd,tyger.log,codespec.yml}
cd /Users/hoomham/Hooman/Work/Codes/2026_ASAP_Recon/workspace
LOG=/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/tyger_bexport_2026-09-25_batch.log
echo "[$(date +%H:%M:%S)] batch start $(wc -l < pipeline/bexport_todo.txt) sessions" >> $LOG
cat pipeline/bexport_todo.txt | xargs -P 3 -I{} zsh -c './pipeline/run_bexport_2026-09-25.sh {} >> '$LOG' 2>&1'
echo "[$(date +%H:%M:%S)] batch end" >> $LOG
