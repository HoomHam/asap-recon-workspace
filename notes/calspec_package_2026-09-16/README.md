# CALSPEC data package — 2026-09-16 (XeCS session)

Built by `2026_XeCS_Recon/workspace/helpers/calspec_package.py`. Source of truth = the XeCS workspace
`outputs/calspec/`; this folder is a COPY. Scripts: `asap_k0_classify.py` (k0), `calspec_extract.py` (stage 1),
`calspec_cohort.py` (stage 2). Raw caches (complex cal FIDs, k0 per rep) on Ext:
`/Volumes/HoomHamExt/Work/Codes/2026_XeCS_Recon/calspec/{cache,k0}/`.

## session_table.csv — one row per dynamic session (104)
| column | meaning |
|---|---|
| session, sid, cohort, seq, ncha | cohort from `cohort_map.csv` (HC young, HC_OLD >45, EBV incl. pre-valve, LTX, RT); seq = spiral sequence version; ncha = coil channels |
| recon_twix | MID of the twix the AIkill recon used (= largest dynamic .dat; audited 104/104 by the ASAP session) |
| hold_twix, hold_grade, hold_T1eff_s, hold_snr | the single-breath breath-hold (same sequence): HOLD_GOOD / HOLD_BAD (jump/drop or poor log-linear fit; fuzzy, some salvageable) / HOLD_LOW_SNR; T1,eff of the gas from the k0 decay |
| extra_free_dynamics | other free-breathing dynamics in the session (MID and duration). Rule: combine the GOOD ones with the recon twix before the next recon (025VP 11-02, 008CR 01-31, 030DN 03-06, 013VM 08-12) |
| ambiguous_twix | k0 trace neither clean hold nor regular breathing — look at k0_contact_sheet.png |
| db_stat_T1_s | session DB value, for comparison; garbage where no breath-hold exists |
| spec_block | dissolved-phase spectroscopy (cal) block on the recon twix: NONE / XE_IN_BLOCK (xenon arrived inside it, usable) / XE_AFTER / NOISE |
| spec_note, spec_n_cal_reps | arrival rep / block length (v3: 30 reps = 780 FIDs = 26.7 s; v2: 25 reps = 500 FIDs = 16.5 s) |
| spec_carrier, spec_carrier_ppm_above_gas, spec_gas_snr | where the spectroscopy demod sat relative to the gas line (ON_RBC = 218 in almost all); gas-line SNR |
| pooled_* | BREATHING-IGNORED spectrum: every FID from xenon arrival to the end of the block, time-domain M2 fit (gas L, RBC L, membrane Voigt). Ratios are UNCORRECTED for the spectrally selective pulse (RBC/gas, mem/gas are sequence numbers; RBC/mem is comparable within v3). Shifts in ppm from the fitted gas line; FWHM in ppm; T2* = 1/(pi FWHM) in ms. `_te` = amplitude extrapolated to the pulse (TE 0.62 ms + 60 us). M3 = two membrane lines (Robertson 2017) |
| n_breath_cycles, phase16_status, phase8_status | complete trough-to-trough breathing cycles inside the block; OK = phase-resolved curves exist in cohort/subject_bins_b16.csv / _b8.csv |
| pooled_rbc_snr | RBC-line SNR proxy = pooled dissolved SNR x RBC/mem |
| phase8_valid_membrane | yes when the per-subject 8-bin MEMBRANE curves (shift, width, mem/gas) are trustworthy: pooled dissolved SNR >= 30.0 (chosen from the data: membrane shift sd across bins <= 0.3 ppm for every such subject) |
| phase8_valid_rbc | yes when the per-subject 8-bin RBC curves (RBC shift/width, RBC/mem, RBC/gas) are trustworthy: RBC-line SNR >= 15.0 (below it the RBC shift scatters > 1 ppm across bins) |

## other files
- `k0_classify.csv` (201 twix), `k0_sessions.csv` (104), `k0_contact_sheet.png` — every ASAP twix by its gas k0.
- `extract_summary.csv` — stage-1 cal-block extraction per twix; `contact_sheet_dissolved.png`, `failed_blocks_raw.png`.
- `cohort/subjects_pooled.csv` — full pooled fit per usable block (M2 + M3, all parameters).
- `cohort/subjects_b16.csv`, `subject_bins_b16.csv`, `cohort_bins_b16.csv` (and `_b8`) — phase-resolved per subject and cohort-average fits.
- `cohort/cohort_pooled.csv`, `spectra_pooled.png` — breathing-ignored cohort-average spectra (M2 + M3).
- `cohort/spectra_<cohort>_<seq>_b16.png`, `params_<cohort>_<seq>_b16.png` (and `_b8`) — phase-resolved figures.
- `CalSpec_Suggestion_2026-09-16.md` (method), `CalSpec_Cohort_Report_2026-09-16.md` (results), `Recon_Input_Rule_2026-09-16.md`.
