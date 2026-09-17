# CalSpec cohort report — phase-resolved dissolved-phase spectroscopy from the cal block

**STATUS: FIRST FULL RUN, 2026-09-16 03:40. Numbers are UNCORRECTED for the Spec pulse excitation profile
(RBC/gas, mem/gas are sequence-weighted; RBC/mem is also profile-weighted, same weight for every v3 session).
v2 is never pooled with v3 for ratios.** Method: `CalSpec_Suggestion_2026-09-16.md` §2. Scripts:
`helpers/calspec_extract.py` → `calspec_cohort.py` (uses `calspec_fit.py`). Outputs: `outputs/calspec/cohort/`.

## 1. Who is in

| | usable block (stage 1) | ≥ 1 complete breathing cycle + SNR ≥ 8 (stage 2) |
|---|---|---|
| v3 | 57 | **43 sessions / 36 subjects** (HC 5, EBV 8, LTX 23, RT 7) |
| v2 | 18 | 4 (HC 1, LTX 1, RT 2) |

28 blocks dropped as NO_CYCLE: xenon arrived too late in the block for one trough-to-trough cycle
(all but 4 of the v2 blocks — 16.5 s long — and 15 v3). A relaxed rule (phase from one trough + the
imaging period) would recover some; not done.

Phase surrogate: the block's own gas band (v3, 41/43) or dissolved band; period prior from the
imaging k0 autocorrelation; trough (end-expiration) → trough = one cycle; 16 bins → ~16 FIDs per bin per
session, ~700 per bin in the v3 cohort average.

## 2. Per-subject numbers (pooled over the cycle, M2 fit, median [IQR]) — v3

| cohort | n | RBC/mem | RBC ppm | mem ppm | RBC FWHM ppm | mem FWHM ppm | T2* RBC ms | T2* mem ms | gas FWHM ppm (lower bd) |
|---|---|---|---|---|---|---|---|---|---|
| **HC** | 5 | **0.45 [0.43–0.69]** | 217.7 | 197.6 | 10.7 | 8.8 | 1.68 | 2.05 | 1.33 |
| EBV | 8 | 0.27 [0.20–0.33] | 217.8 | 197.3 | 11.0 | 9.1 | 1.65 | 2.00 | 1.44 |
| LTX | 23 | 0.30 [0.22–0.35] | 218.0 | 197.7 | 10.8 | 8.6 | 1.68 | 2.09 | 1.57 |
| RT | 7 | 0.33 [0.26–0.33] | 217.6 | 197.5 | 10.1 | 8.3 | 1.79 | 2.18 | 1.49 |

Cohort-average-FID fits agree: RBC/mem HC 0.57 (M2) / 0.43 (M3); EBV 0.25 / 0.26; LTX 0.28 / 0.23; RT 0.29 / 0.26.
Chemical shifts and linewidths are the same in every cohort to ±0.3 ppm: **RBC 217.9, membrane 197.6 ppm;
RBC FWHM 11.0, membrane 8.7 ppm; T2* RBC 1.64 ms, membrane 2.09 ms (1.5 T)** — vs ~1.0 ms at 3 T in the
literature, i.e. the expected ≈2× for susceptibility-dominated lines. TE-extrapolated RBC/mem is ~10 % higher
(RBC decays faster than membrane over the 0.68 ms to the first fitted sample).

M3 (two membrane lines) on the v3 cohort average: **mem1 200.8 ppm / 12.4 ppm FWHM, mem2 196.0 / 8.1**,
BIC −715 vs M2 (strongly preferred). Robertson 2017 (3 T): 201.3 / 14.6 and 196.3 / 9.3. The split reproduces in
every cohort (EBV 201.4/195.7, LTX 200.5/196.0, RT 201.1/196.1, HC 200.6/196.2).

## 3. Across the breathing cycle (v3 cohort average, phase 0 = end-expiration trough)

| quantity | trough (0) | mid-cycle (~0.4, end-inspiration) | trough (1) |
|---|---|---|---|
| membrane amplitude (a.u.) | 18 | **37** | 25 |
| RBC amplitude | 4 | 7.5 | 5 |
| RBC/mem | 0.32 | 0.28 | 0.32 |
| membrane shift, ppm | 197.7 | **197.4** | 197.7 |
| RBC shift, ppm | 218.1 | 217.8 | 217.8 |
| gas frequency, ppm rel. mean | −0.07 | **+0.05** | −0.07 |
| membrane T2*, ms | 2.15 | **2.03** | 2.15 |
| RBC T2*, ms | 1.65 | 1.65 | 1.65 |
| mem/gas, RBC/gas | 13.5, 4.5 | **10, 3** | 14, 4.7 |

So: on inspiration the dissolved signal doubles (uptake follows the fresh gas), the membrane line shifts
−0.3 ppm and broadens (T2* 2.15 → 2.03 ms) at full inflation, the gas line moves +0.12 ppm peak-to-trough,
RBC shift and RBC T2* do not move, RBC/mem dips ~10 % at end-inspiration. HC shows the same shapes with
RBC/mem ≈ 0.55 flat. The wash-in staircase is removed by the per-cycle detrend only in the SURROGATE; the
amplitudes themselves are per-bin means over 2–4 breaths of a rising signal — a slow drift across the
cycle (bin 15 < bin 0 is possible) is that, not physiology.

Figures: `spectra_<cohort>_v3.png` (16 bin-averaged spectra), `params_<cohort>_v3.png` (cohort M2 thick,
per-subject median ± IQR band, subject spaghetti). Tables: `subjects.csv` (pooled per twix, M2 + M3),
`subject_bins.csv` (per twix × bin, M2), `cohort_bins.csv` (per cohort × seq × bin, M2 + M3).

## 4. Caveats (all stated in the doc, repeated here)
- Excitation profile unknown (no pulse shape/duration in the header) → ratios uncorrected.
- 16 FIDs/bin per subject: per-subject per-bin curves are noisy (see spaghetti); the cohort curves are not.
- Wash-in, not steady state: the first 2–4 breaths after arrival.
- Gas linewidth is a lower bound on T2* (30.7 ms readout); gas FREQUENCY is well determined.
- v2: 4 sessions only; carrier between/on RBC, gas 10× weaker.
