#!/bin/sh
# waits for prep merge, then runs refine+overlay+video for every manifest key, 3 in parallel
cd /Users/hoomham/Hooman/Work/Codes/2026_ASAP_Recon/workspace
until [ -s outputs/ct_overlay/logs/prep_merge.log ]; do sleep 20; done
helpers/.venv/bin/python helpers/ct/ct_cohort_batch.py keys > outputs/ct_overlay/logs/cohort_keys.txt
xargs -P 3 -I{} sh -c 'helpers/.venv/bin/python -u helpers/ct/ct_cohort_batch.py run {} > outputs/ct_overlay/logs/run_{}.log 2>&1; echo "DONE {} $?" >> outputs/ct_overlay/logs/run_summary.log' < outputs/ct_overlay/logs/cohort_keys.txt
echo ALL_DONE >> outputs/ct_overlay/logs/run_summary.log
