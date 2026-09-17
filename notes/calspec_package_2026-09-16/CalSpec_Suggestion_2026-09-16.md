# CalSpec — cohort-wide dissolved-phase spectroscopy from the cal block

**STATUS: SUGGESTION (stage-2 method), 2026-09-16. Stage 1 (extraction) BUILT and run — see §1.**
Hooman's order: *"done for all the 104 subjects who have it … an average for all … phase-resolved
(16 bins if possible, else 8) … RBC/gas, tissue/gas, RBC/tissue, T2* (or linewidth), chemical shift
for RBC and tissue, for one breathing cycle … gas may shift a little … one param each for RBC and
tissue misses the shoulders … look at the last 7 years of literature … if noisy ignore … almost all
had the pulse on RBC (or between the two peaks), a few the opposite."*

Precedent: `CalBlock_Recipe.md` (n=1, 025JC, F100–F102). Everything there still holds: cal block =
gradient-free dissolved-phase spectroscopy, **60 µs dwell** (NOT the 5 µs spiral ADC), 512 pts,
carrier on/near RBC, gas −218 ppm from carrier (narrow), membrane −22 ppm from RBC (broad).

---

## 1. Stage 1 — what the cohort actually contains (`helpers/calspec_extract.py`)

Every `meas_*fa_spiral_dyn_fancy*.dat` in the 104 AIkill session dirs (drive folder is now
`/Volumes/HoomHamExt/AIdiff_FMIG_Drive_Data/dynamic/human/` — the recipe's `FMIG_Drive_Data` path is
stale). Per twix: complex read → cal block **by cadence** → 60 µs spectra → gas line auto-located
(narrowest strong line that has the dissolved pair 185–235 ppm above it) → carrier position from the
MEASURED RBC/membrane peaks → cache on Ext (`/Volumes/HoomHamExt/Work/Codes/2026_XeCS_Recon/calspec/`,
mirror-structure law), QC png per twix, one row per twix in `outputs/calspec/extract_summary.csv`.

Cross-check that costs nothing: `len(fid_signal) − 2·len(nav_ilvtime) == n_cal·nlin` (recon.mat vs
twix). Holds on 025JC (780) and 000LL v2 (500).

Sequence variants seen in the headers:

| protocol | coil | lines/rep | cal reps | Spec pulse | carrier (measured) | gas SNR (late reps) |
|---|---|---|---|---|---|---|
| `v3_20240130` | 1 | 26 | 30 (25/22 rare) | `Spec_ex` | **ON RBC** (gas −218) | ~600 |
| `v2_20230131` | 1 | 20 | 25 (30 rare) | `Spec_ex` | **BETWEEN** (RBC +5.5, mem −17, gas −214) | ~50 — gas 10× less excited than in v3 |
| `v3_20230821` | 8 | 26 | **0** | none (`GP_ex` only) | — | no spectroscopy block at all |
| `v2_20230131` | 8 (early 2023) | 20 | **0–1** | | — | none |

⇒ **v3 single-coil (~55 sessions) is the main cohort; v2 single-coil (~25) is a second group with a
different carrier and a weaker gas line; the 8-channel sessions have no cal block.** Counts and the
full per-twix table: §1b (filled from the run).

Header items worth knowing: `sWipMemBlock.adFree[2] = 218` (dissolved offset, ppm) in every variant;
`adFree[0]` = 3.2 (v3) / 3.6 (v2) / 9–10 (8-ch); `adFree[4]` = 1.28 / 1.024; TE 620/600 µs
(`alTE[0..1]`); Larmor 17.612 MHz → **17.61 Hz/ppm**; Spec pulse amplitude 133 V (v3) / 131 V (v2).
The Spec pulse DURATION is not in MeasYaps — Hooman knows it (needed for the profile caveat, §2.4).

### 1b. Cohort table — the run (2026-09-16 01:00–01:30)

201 dynamic twix in 104 sessions, all read. `outputs/calspec/extract_summary.csv` = one row per twix.

| status | twix | what |
|---|---|---|
| NO_CAL | 106 | no cadence-detectable cal block: every 8-ch file (v2 early-2023 + `v3_20230821`), the SECOND dynamic of most v3 sessions (repeat acquisitions carry NO cal block — so "two cal blocks per session" does not happen), 011CN 04-29 |
| OK, xenon arrived in-block | **75 (74 sessions)** | **the usable set: 57 `v3_20240130` + 18 `v2_20230131`, all single-coil** |
| OK, no arrival | 2 | 002ZS 11-06 (v2), 005JJ 05-13 (v3): block is noise, xenon came after it |
| DWELL_OR_LINE_FAIL | 17 | 12 v2 (see the lone-line group below), 2 v3 (030DN 03-06, 017ED: SNR < 40, noise), 2 8-ch v3 (000LL/005DS 02-02: flat after a 300× sample-0 spike — **parked**), 1 v2 |

Carrier position (ppm above the detected gas line), usable set: **ON RBC (215–222) in 70/75**;
BETWEEN in 2 (000LL 10-27 at 214, 041WF 12-19 at 206); ON MEM 1 (004DS 01-16 at 197);
BELOW MEM 1 (038RL 11-16 at 175); ABOVE RBC 1 (006MM 08-15 at 231). **No carrier-on-gas session
among the usable ones.**

★ **The lone-line group (v2, 2023-05 → 2023-12: 000SK, 002BB 10-10, 002JM, 039CP 11-21, 003PM 12-06,
+ ~5 noise-level):** ONE broad line at −20…−31 ppm from carrier (FWHM 8.5 ppm, T2* 2.1 ms —
dissolved-like), no partner at 0 or −218, and the imaging k0 says the lung had no gas until the
block ended (002BB: 31 s). Either (a) membrane with no RBC and no gas (unphysical), or (b) the
"opposite" carrier (on gas) and this is xenon in the BAG beside the coil, fringe-field broadened.
**Hooman's call — the v2 block is 16.5 s, xenon usually arrives at 14–30 s, so most v2 blocks miss
the lung either way.**

★★ **THE CONTACT SHEET (`outputs/calspec/contact_sheet_dissolved.png`): in almost every session the
MEMBRANE line dominates and RBC is a shoulder** — RBC/membrane ≈ 0.1–0.3 by eye, vs ≈ 1 on 025JC
(the recipe's n=1) and on 000LL (both sessions), 016PG, 021CH. The pulse is centred ON RBC, so this
is not the profile suppressing RBC. Patients (RT/EBV/LTX) with low RBC transfer is the obvious
reading (Robertson: 0.18 in IPF vs 0.58 healthy); the cross-ref tags will say. 025JC was an
unrepresentative n=1.

**Retraction (same night):** the "lone-line" group is NOT bag gas. `outputs/calspec/failed_blocks_raw.png`
(full spectrum + band time courses + imaging k0 for all 19 failed/no-arrival blocks) shows it is the
MEMBRANE line alone, RBC too weak to see, the gas line weak (v2 pulse), xenon arriving in the last
3–5 s of the 16.5 s v2 block. The rest of the failures are pure noise (xenon after the block).

**How the carrier position is read (Hooman's Q, 2026-09-16):** the GAS line is the only narrow line
(T2* ≥ 10 ms → ~1 bin) — its offset from 0 ppm says where the demod sits relative to gas; RBC and
membrane are then expected at gas+218 and gas+197. Demod = F0 + `sWipMemBlock.adFree[2]` ppm. Measured
carrier ≠ 218 on five usable twix, each explained by the header:

| session | measured carrier (ppm above gas) | `adFree[2]` | F0 vs neighbouring days | reading |
|---|---|---|---|---|
| 038RL 2023-11-16 | **175** | 218 | **−800 Hz (−45 ppm)** | frequency adjust missed the gas that day |
| 006MM 2024-08-15 | **231** | 218 | **+236 Hz (+13 ppm)** | same, other sign |
| 041WF 2023-12-19 | 206 | **208** | ok | "between the peaks" was SET (208) on this day |
| 004DS 2024-01-16 | **197 (on membrane)** | **208** | −11 ppm | 208 set + F0 low → landed on membrane |
| 000LL 2023-10-27 | 214 | 218 | −4 ppm | small F0 miss |

⇒ a 45 / 13 ppm F0 error also sits on the IMAGING of 038RL 11-16 and 006MM (spiral off-resonance
blur, gas k-space centre-frequency offset) — **note for the recon side; not pursued here.** `adFree[2]`
= 208 exists only on 2023-12-19 and 2024-01-16.

**The 8-ch v3 blocks (000LL / 005DS 2024-02-02):** cadence finds 30 cal reps, but the combined FIDs are
a 300× spike at sample 0 followed by noise on every coil (per-coil check done); no xenon line anywhere.
Parked: 2 sessions, and the 8-ch coil's spec block may simply not have been received. Not a combine bug.

Also measured on 025JC (`helpers/calspec_fit.py`, n=1 draft `outputs/calspec/n1_*_bins16.png`):
sample 0 of every cal FID is an ADC transient (|FID| 2.6e-5 vs 7.4e-5 at sample 1) → fit from
sample 1 (30× lower residual). M3 reproduces Robertson's two barriers almost exactly: mem1 200.1 ppm
/ 12.5 ppm FWHM, mem2 196.7 / 8.9, 99° apart (Robertson 201.3/14.6, 196.3/9.3, ~100°). RBC 216.8–217.5
ppm, 9–10 ppm FWHM; **dissolved T2* 1.8–2.6 ms at 1.5 T** (≈ 2× the 3 T literature, as expected);
gas T2* ≥ 13 ms (lower bound). Two complete breathing cycles after arrival → 16 FIDs/bin/session at
16 bins. RBC/mem flat across the cycle (≈1.0), RBC/gas and mem/gas swing 8→13 (the gas breathes),
gas frequency swings +0.15 ppm across the cycle, linewidths flat.

---

## 2. Stage 2 — the fit (suggestion)

### 2.1 Literature (2017–2026), what people actually fit

| paper | model | notes |
|---|---|---|
| **Robertson 2017 MRM** (PMC5422139) "third dissolved-phase resonance" | time-domain, complex, sum of `a·e^{iφ}·e^{2πift}·e^{−πwt}` (Lorentzian); **gas ×2** (airway 0, alveolar −2.3 ppm), **RBC ×1, barrier ×2**; trust-region LSQ | healthy: RBC 216.5 / B1 201.3 / B2 196.3 ppm, FWHM 11.5 / 14.6 / 9.3 ppm, RBC:barriers 0.44. B1 = tissue near the alveolar wall (broad), B2 = plasma/interstitial (narrow, ~7 ppm like in-vitro plasma); ~100° phase between B1 and B2. 2-peak dissolved fits leave **patterned residuals**. |
| **Bier 2019 NMR Biomed** (PMC6447038) dynamic spectroscopy, cardiogenic oscillations | time-domain complex; **gas L, RBC L, barrier Voigt** (L⊗G, one extra parameter instead of a whole extra line) | chosen for STABILITY at low per-FID SNR in dynamic data — exactly our regime. Preproc: SIFT filter along time + 5-FID boxcar. RBC 218 ppm, RBC:barrier 0.58 healthy. |
| **Niedbalski 2021 MRM** consortium position paper | gas 0, RBC 218, "tissue/plasma" 197 ppm; recommends dissolved calibration spectroscopy | terminology + reference shifts. |
| **Kammerman 2020 MRM** | dissolved R2* as IPF biomarker | linewidth is a biomarker in its own right. |
| **Leewiwatwong 2025 MRM** (multisite 3 T reference) | global T2* membrane 0.99 ± 0.04 ms, RBC 1.04 ± 0.03 ms at 3 T | at our **1.5 T** expect ≈ 2× longer if susceptibility-dominated — a number this cohort can measure. |
| **Costelle 2026 NMR Biomed** | CSI at 3 T, T2* membrane 1.07 / RBC 1.13 ms | same. |
| 2023 (PMC10380088) glycation third RBC line | extra RBC resonance at high glucose | a 4th dissolved line is on the menu if the residual demands it; not by default. |

Consensus: **fit the complex FID in the time domain**, never a phased magnitude spectrum; Lorentzian
for gas and RBC; membrane needs more than one Lorentzian — either Voigt (Bier, cheaper) or two
Lorentzians (Robertson, interpretable). Chemical shifts referenced to the fitted gas line = 0 ppm.

### 2.2 Recommended model ladder (BIC-selected on the cohort average, fixed per bin)

Model `S(t) = Σ_n a_n e^{iφ_n} e^{2πi f_n t} L_n(t)`, `t` from the first ADC sample:

| | gas | RBC | membrane | params |
|---|---|---|---|---|
| **M1** | L | L | L | 12 |
| **M2** (Bier) | L | L | **Voigt** (L⊗G) | 13 |
| **M3** (Robertson) | L | L | **L + L** | 16 |
| M4 | L | **L + L** | L + L | 20 — only if M3's residual still shows structure at RBC |

Fit M1–M3 on the cohort-average complex FID (all bins pooled, then per bin). Report BIC and the
residual spectrum for each. **Use M2 for the per-subject / per-bin work** (13 params, stable), and
**M3 on the cohort average** where SNR affords it — that is where the "shoulders" Hooman mentions get
their own numbers. Per-line outputs: amplitude (at t=0 of the FID, and TE-corrected — §2.4), f (ppm
from gas), FWHM (Hz and ppm), T2* = 1/(π·FWHM_L), phase.

Gas caveat: gas T2* ~10 ms ≥ readout window (30.7 ms) only partially — the gas linewidth from a 30 ms
FID is a **lower bound on T2*** (resolution 32.6 Hz = 1.85 ppm/bin; the time-domain fit does better
than a bin but not by much). Report it as such. Gas frequency, on the other hand, is well determined
(sub-Hz from a 30 ms coherent line) — that is the "gas may shift a little" track, per FID.

### 2.3 Phase-resolved binning

- **Usable window** = reps after xenon arrival (arrival rep 15–20 of 30 on the two tested; ~10–15 s,
  **2–4 breaths** per twix). Sessions with two dynamics contribute two cal blocks; both are kept,
  tagged primary / repeat (Hooman's split-acquisition rule).
- **Surrogate**: the FID's own band integrals — no cross-clock step at all. v3: gas band (SNR ~600);
  v2: dissolved band (gas too weak). Detrend the wash-in staircase with a per-breath running median,
  troughs → within-breath fraction `(t − t_trough)/(t_next − t_trough)` exactly as `piston_phase` does
  for the imaging block. Pneumotach (clock settled, b=0 on most, `clk_b_s` in the session DB) is the
  QC, not the driver.
- **Budget**: per-FID dissolved SNR ≈ 1.4 (025JC: peak/noise 20 on a 208-FID mean). 16 bins × ~30
  FIDs/bin per v3 twix → per-subject-per-bin SNR ≈ 8 (marginal for linewidths, OK for the RBC/mem
  ratio); cohort of ~55 v3 twix → ~1600 FIDs/bin → **SNR ≈ 60 per bin: 16 bins is comfortable for the
  cohort average**. Per subject: 8 bins (SNR ≈ 11) with an SNR gate; noisy sessions dropped and listed.
- **Averaging across subjects** — do it on the complex FIDs, not on magnitude spectra:
  1. demodulate each FID so its OWN fitted gas frequency sits at 0 Hz (removes per-subject shim /
     carrier offsets; keeps RBC/mem-from-gas exact),
  2. remove the per-subject zero-order phase (phase of the gas line at t=0),
  3. normalise by the subject's mean dissolved amplitude over the window (equal weight per subject),
  4. average per bin, then fit. Also fit per subject per bin and report median ± IQR of the parameters
     — the two estimators must agree or the average is hiding a mixture.
- **Carrier groups are NOT pooled for ratios** (§2.4): v3 (ON RBC) and v2 (BETWEEN) averaged
  separately; shifts and linewidths may be pooled (profile-independent), ratios may not.

### 2.4 Caveats that change numbers — Hooman's calls

1. **Excitation profile.** `Spec_ex` is a long, spectrally selective pulse. RBC/mem/gas amplitudes are
   weighted by its profile at 0 / −22 / −218 ppm (v3) or +5.5 / −17 / −214 (v2). **RBC/gas and mem/gas
   are sequence numbers, not physiology** (the recipe already says so: dissolved/gas 8.6 on 025JC);
   **RBC/mem is biased too** — the two lines sit at different points of the profile, and differently in
   v2 vs v3. Correction needs the pulse shape + duration (not in MeasYaps). Options: (a) report
   uncorrected per group and say so; (b) Hooman supplies the pulse (duration, shape) → simulate the
   profile → divide; (c) calibrate the v2↔v3 profile ratio empirically from subjects scanned under both
   (several were: 000LL, 023DB, 030DN, 032WS, 036RL, 037GD, 038RL, 039CP, 003KH, 003PM, 004DS, 005DS,
   002ZS, 041WF, 025VP, 029CK…). (c) is free and is my recommendation as a first pass.
2. **TE.** First sample at TE ≈ 0.6 ms; dissolved T2* ~1–2 ms ⇒ 25–45 % of the dissolved amplitude is
   already gone at t=0 of the FID, and RBC vs membrane decay differently. Report amplitudes both as
   fitted (t=0 of FID, what everyone publishes) and extrapolated to the pulse (`a·e^{TE/T2*}`).
3. **Flip-angle depletion inside the cal block.** 26 (v3) / 20 (v2) FIDs per rep at 34 ms spacing, with
   the Spec pulse: the dissolved pool is replenished by exchange (fast, ~ms) but the gas pool is
   not — the gas line is depleted along a rep if the Spec pulse tips it at all. Check: gas amplitude vs
   line index within a rep; if there is a slope, the gas reference is line-dependent and must be
   modelled (or only line 0 used as the gas reference).
4. **Wash-in is not steady state.** RBC/mem during the first three breaths of xenon may not equal the
   breath-hold value (dissolved uptake lags). The phase-resolved curve is for "the first breathing
   cycles after arrival" — say so in every figure.
5. **Gas line: two components?** Robertson resolves airway (0) and alveolar (−2.3 ppm) gas at 3 T with
   a 1.7 ms T2*-limited line. At 1.85 ppm/bin we cannot; single gas Lorentzian.

### 2.5 Deliverables (per carrier group; v3 = main)

- `outputs/calspec/cohort_bins_<group>.csv`: per bin × model — a, f (ppm from gas), FWHM (Hz, ppm),
  T2*, φ for gas / RBC / mem(1,2); RBC/gas, mem/gas, RBC/mem (fitted and TE-corrected); n FIDs, n
  subjects, SNR.
- `cohort_average_spectrum_<group>.png`: waterfall of the 16 bin-averaged spectra with the fit and the
  residual; `cohort_params_vs_phase_<group>.png`: each parameter vs breathing phase, cohort curve ±
  IQR with per-subject spaghetti underneath (Hooman's "plot the inputs" rule).
- `per_subject_bins.csv` (8 bins) + SNR-gate list of excluded sessions with the reason.
- The stage-1 QC contact sheet (all cal spectra on one page) — inputs first.

### 2.5b Cohort membership (`outputs/calspec/cohort_map.csv`, built 2026-09-16 from the clinical
subject sheets: RT = the RT overview list; LTX = sheet says "Transplant"; HC = "healthy control" +
Hooman's named young healthy 021CH / 000LL; EBV = sheet says "Valve"; 000SK/000SV/000KR/017AK/02BB unknown).
Usable cal-block sessions: **LTX 42 · EBV 15 · RT 9 · HC 7 · unknown 1**. Hooman 2026-09-16: 16 bins
for the whole cohort, **separately for HC / EBV / LTX / RT**; v2 never pooled with v3 for ratios.

---

## 3. Side project — every ASAP twix by its gas k0 (`helpers/asap_k0_classify.py`, `asap_k0_summary.py`)

Hooman's "small project": breath-hold (decay → T1,eff of the gas) vs free-breathing (fluctuation);
grade the holds; find second real free-breathing dynamics. k0 = |sample 1| gas parity per rep, cal
block removed. `outputs/calspec/k0_classify.csv` (201 twix), `k0_sessions.csv` (104 sessions),
`k0_contact_sheet.png` (all 201 traces, colour = class, * = the recon's twix).

| | n |
|---|---|
| sessions with a breath-hold | 83 / 104 |
| HOLD_GOOD (r² ≥ 0.95 on log k0 from peak to 5 %, no jump > 20 % up / 30 % down above 15 % of peak, SNR ≥ 20) | 63 |
| HOLD_BAD (jump/drop or poor fit — list in `k0_sessions.csv`) | 18 |
| HOLD_LOW_SNR | 2 (000LL 11-03, 043AS 12-05) |
| no breath-hold at all | 21 |
| sessions with an EXTRA free-breathing dynamic | 9 |

**T1,eff validation:** where the session DB has `stat_T1_s`, the k0 fit reproduces it to 0.01 s (46/46).
**Six DB values are garbage** because those sessions have no breath-hold and the DB fitted a free-breathing
file: 009JT 02-27 (376 s), 023DB 05-02 (209), 012EB (131), 044DY 03-12 (90), 013VM 08-12 (89), 004DS
01-16 (14). 004DS 07-29: DB 7.63 vs k0 4.79 (DB used a different file). 043AS 12-05: DB −238.

**Extra free-breathing dynamics** (= the "sequence stopped, restarted" cases): short aborted starts
(10–26 s): 02BB 03-27, 025VP 11-02, 004DS 01-16, 023DB 05-02, 007IT 08-28, 008CR 01-31 ×2. Long
enough to be a REAL second dynamic: **030DN 2024-03-06 MID00545 (101 s), 008CR 01-31 MID00194 (58 s),
013VM 08-12 MID00408 (45 s), 009JT 02-27 MID00333 (43 s)**. None of the extras carries a cal block
(NO_CAL in stage 1) — Hooman's prediction holds. Ambiguous traces (eye): 02BB 03-27 MID00429, 002BB
10-10, 003KH 12-15, 002ZS 05-13, 014KS 08-30 — the recon's own twix in four of those.

### 2.6 What I need from Hooman before stage 2 runs
- Pool v2 with v3 for ratios: no (my default) / yes with empirical profile factor / yes after pulse
  simulation (then: pulse shape + duration).
- TE-corrected amplitudes: report both (my default) or one.
- 16 bins cohort + 8 bins per subject (my default) or 16 everywhere.
- Repeat cal blocks (second dynamic of the same session): include as extra data (my default, tagged) or
  primary only.
