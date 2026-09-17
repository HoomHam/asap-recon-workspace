# Recon input rule — check every ASAP twix of a session before reconstructing (Hooman, 2026-09-16)

**Status: STANDING RULE for every future recon (Hooman 2026-09-16: "make this also clear for all the next recons").**

1. Before reconning a session, look it up in `outputs/calspec/k0_sessions.csv` (built by
   `helpers/asap_k0_classify.py` + `asap_k0_summary.py`; contact sheet `outputs/calspec/k0_contact_sheet.png`):
   - `recon_twix` = the file the AIkill recon used (== the largest dynamic .dat; audited 104/104 by the ASAP
     session 2026-09-16, table `2026_ASAP_Recon/workspace/outputs/twix_audit_2026-09-16/twix_audit.csv`).
   - `hold_twix` / `hold_grade` / `hold_T1eff_s` = the single-breath breath-hold (same sequence) and its T1,eff.
   - `extra_free` = other free-breathing dynamics (sequence stopped and restarted).
2. **If a session has a GOOD second free-breathing dynamic, COMBINE both sets into one recon input** (one
   k-space set on one clock; the cal block sits only at the head of the first file). Sessions decided by Hooman:
   **025VP 2023-11-02, 008CR 2024-01-31, 030DN 2024-03-06, 013VM 2024-08-12.** Short aborted starts (10–26 s) are
   not worth merging.
   **Outcome (ASAP session, 2026-09-17, image 40d23a4):** merged re-recons adopted for **030DN, 013VM, 008CR**
   (gas SNR at end-inspiration 24.9→27.7, 27.4→30.1, 29.4→31.9; DP 15.8→18.2, 11.9→14.8); **025VP REJECTED by
   Hooman** (19.1→14.7, DP 10.9→8.4). **Condition added: merge only pays when the two doses have comparable k0
   SNR** — Steve's bins average interleaves with equal weight, so a 2.4× weaker dose (025VP MID02741, k0 SNR 53
   vs 126) dilutes the mean while adding its noise. Merged outputs: `AIkill_Dynamic/<date>_<id>_merged/{s,d}`
   (originals untouched; 025VP's kept as `..._merged_REJECTED`). Tool: ASAP `workspace/pipeline/merge_dyn.py` +
   `dyn_recon.py --merge/--ref`; file A trimmed to whole arm cycles (832 v3 / 640 v2 gas interleaves); p-binning
   unsupported on merged input. Comparison: `2026_ASAP_Recon/workspace/outputs/merged_dyn_2026-09-16/`.
   ⚠ downstream tables keyed on recon.mat (registration, SNR table, calspec session_table `recon_twix`) now have
   two versions for those three sessions.
3. Breathing quality inside the recon's own dynamic (irregular / collapsing / dose-out) is visible on the k0 sheet;
   dropping bad-breathing reps before recon is on the table (Hooman) — not yet implemented.
4. Breath-hold grades are rule-based and fuzzy: some HOLD_BAD are salvageable; T1,eff from `k0_classify.csv`
   supersedes `session_db.stat_T1_s` where the DB fitted a free-breathing file (009JT 02-27, 023DB 05-02,
   012EB, 044DY 03-12, 013VM 08-12, 004DS 01-16, 043AS 12-05).

The same rule was sent to the ASAP recon session (2026-asap-recon-31) to record on its side.
