# Bin timing, session pick, header facts — for PCA Step 1 (global A→B→C) — 2026-10-26

Status: VALID (code-read + measured). Requested by the PCA session (2026-pca-registration-04) on Hooman's order.
Script: `helpers/_bin_timing.py <study> …` (DIAPHRAGM from b44d recon.mat nav arrays; SIGNAL from the old-tree
`signal_pneumo.npz` k0 signal with Steve's EE finder re-implemented). Headers: Siemens ASCCONV of the raw .dat
(`/Volumes/HoomHamExt/_5t_images_roundtrip/Images/<date>/<id>/`).

## 1. Binning mechanics (both recon types share the kernel)

| item | value | code |
|---|---|---|
| soft kernel | `binwt = exp(-d²/bindist0sq)`, `bindist0sq = 2` → σ = 1 bin, FWHM 2.35 bins; w(1) 0.61, w(2) 0.135, w(3) 0.011 | `recon.py:13, 104-111` |
| distance | cyclic: `min(|ib-c|, |ib+16-c|, |ib-16-c|)`; c = `bins*16` continuous in [0,16) | `recon.py:109-110`, `results.py:291` |
| bin centre | integer bin ib centred at fraction ib/16 (NOT (ib+½)/16) → bin 0 straddles the wrap | — |
| per-bin image | per k-space cell Σw·raw / Σw (weighted MEAN; amplitude independent of count) | `recon.py:56-63` (cudarenorm) |
| excluded ilv | mask 0 → raw 0 → skipped in sum and norm; bin −1 → skipped | `results.py:292-296`, `recon.py:101,107` |
| effective count | w_eff = σ√2π = 2.51 × hard count per bin (measured 2.51 on every session) | — |
| no DCF, no time term, no per-breath weighting | — | — |

## 2. DIAPHRAGM (b44d production) — rank bins of the navigator

| step | what | code |
|---|---|---|
| navigator frame | 26 consecutive gas interleaves (one Thomson set) = 26 × 17.2 ms = 0.447 s, low-res image MS = IS+4 | `tyger_recon.py:47-63` |
| position | coronal |projection| (middle half in y) → z-profile → half-max crossing by parabola on the 25–75 % transition, walking from the high-z (diaphragm) end | `tyger_recon.py:72-98` |
| smoothing | Savgol(5,2) over frames; (x−min)/(max−min) over the WHOLE scan | `tyger_recon.py:100-103` |
| to interleaves | `np.interp` to each gas-interleave time (frame time = mean of its 26 ilv), nan outside the first/last valid frame | `tyger_recon.py:102-103` |
| `nav_volume` | exactly that array (ntotalilvs,), nan = no navigator (before first frame centre; after the last valid frame, i.e. inside the auto low-signal exclude range) | `tyger_recon.py:113` |
| bins | `raw.bin()` (verbatim Steve): direction = sign of v[i]−v[i−1] on consecutive valid entries (equality → inhale); inhale sorted ascending, exhale reflected `maxiv·1.00001 + maxev − v`; bin = rank/len → equal counts per bin, pooled over ALL breaths | `raw.py:157-181` |

Consequences: bin = absolute diaphragm position (direction-split), not phase within own breath. EI boundary sits at
rank = fraction of interleaves classed inhale (015HK/039CP → bin 7, 014KS → bin 8). Bin 0 = EE pause from both
directions. Direction is piecewise-constant over 0.447 s frames (linear interp). Last inhale bin pools the EI
plateau (volume max, bolus already in for a while) → part of the "gas peaks one bin before volume" in this binning.

## 3. SIGNAL (self-gated) — time-phase bins

| step | what | code |
|---|---|---|
| signal | Σ|first 8 samples| of each gas line, summed over coils, min-max normalised | `raw.py:461-463` |
| EE | strict local minimum within ±N = 50 gas interleaves (±0.86 s), two passes (`sigminN` MRD param, fork 262b16f) | `raw.py:470-490` |
| bin | `(t − t_EE,k) / (t_EE,k+1 − t_EE,k)` → pure time phase of own breath; before first / after last EE → −1 | `raw.py:500-510` |

Bin b is centred at b/16 of that breath's period: seconds per bin = T_k/16, exact. Bin 0 = k0 minimum = instant
before the bolus signal rises = phase A by definition (not a volume minimum). Count per bin ≈ equal.

## 4. Measured (helpers/_bin_timing.py)

| study | recon | used / gas ilv | breaths | period median (IQR) s | hard/bin | w_eff/bin | s per bin | phase IQR in bin |
|---|---|---|---|---|---|---|---|---|
| 015HK | DIAPHRAGM | 5772 / 7017 | 23 | 4.02 (4.02–4.45) | 361 | 904 | ~0.25 | 0.05–0.07 |
| 015HK | SIGNAL | 6655 / 7017 | 28 | 4.08 (3.89–4.26) | 416 | 1043 | 0.255 | 0.031 (ideal) |
| 014KS | DIAPHRAGM | 8268 / 8945 | 25 | 4.94 (4.49–5.37) | 517 | 1295 | ~0.30 | 0.06–0.10; bins 3–5 0.10–0.39 |
| 014KS | SIGNAL | 8400 / 8946 | 26 | 5.05 (4.80–5.31) | 525 | 1316 | 0.316 | 0.031 |
| 039CP | DIAPHRAGM | 8268 / 8876 | 29 | 4.90 (4.47–4.94) | 517 | 1295 | ~0.30 | 0.06–0.09 |
| 039CP | SIGNAL | 7669 / 8877 | 28 | 4.64 (4.33–5.13) | 479 | 1201 | 0.29 | 0.031 |

015HK DIAPHRAGM median seconds after EE per bin: b1 0.26, b2 0.52, b3 0.77, b4 1.00, b5 1.26, b6 1.50, b7 1.75,
b8 2.08, b9 2.27, b10 2.53, b11 2.79, b12 3.03, b13 3.29, b14 3.54, b15 3.85.
015HK: 7017 − 5772 = 1245 ilv = 21.4 s = 99.5 → 120.7 s (auto exclude range 100–129 s → nav frames ≥ 223 skipped).

CAVEAT: the s/ (SIGNAL) recons on Ext are the FIRST batch (`AIkill_Dynamic/<study>/s/recon.mat`, fork
472fbc9/40d23a4): no deapod, calcb phase low-pass not 4.4 vox — NOT the b44d image. Same-session ratios are
comparable, not identical. A b44d-image SIGNAL rerun is ~5 min/session on Tyger (needs Hooman's go).

## 5. Pneumotach recons — located, NOT validated

93 in `AIkill_Dynamic/<study>/p/recon.mat` (all logged `binned dynamic recon with PNEUMOTACH binning`), none in
b44/b44d. Why not to trust: `raw.py:184-202` zeroes the pneumo clock at its FIRST PACKET and crops to [0, scan end]
— alignment to the scanner is an assumption (no trigger). Pneumo span 3.8–162 s (015HK), 4.6–197 s (014KS),
4.8–433 s (039CP) vs 121/154/153 s gas clock. Volume = cumsum(pressure) with a deg-8 polynomial through the EE
minima removed (`raw.py:491-493`).

## 6. Header facts (Siemens ASCCONV)

v3 protocol `fa_spiral_dyn_fancy_v3_20240130` (015HK, 014KS, 016PG, 025JC, 039CP — identical except amplitudes):

| item | value | field |
|---|---|---|
| TR | 17.2 ms = one gas + one dissolved excitation per TR (dual-phase kernel) → gas-to-gas 17.2 ms, each pool 58.1 pulses/s | `alTR[0] = 17200` |
| TE | 0.62 ms | `alTE[0] = 620` |
| gas flip (nominal) | 3.5° | `adFlipAngleDegree[0] = 3.5` |
| WIP adFree[0] | 3.4 (015HK), 3.3 (014KS), 3.9 (016PG), 3.2 (025JC), 3.8 (039CP) — per-session, meaning UNKNOWN (calibrated gas flip? ask Steve) | `sWipMemBlock.adFree[0]` |
| dissolved flip (nominal) | 30° | `sWipMemBlock.alFree[3] = 30` (also alFree[13] = 43 ?) |
| dissolved offset | +218.0 ppm from gas = RBC line (ON_RBC, spec carrier 218.14 ppm in the calspec table) | `sWipMemBlock.adFree[2] = 218.0` |
| pulse amplitudes | GP_ex 14.43 V, DP_ex 98.98 V, Spec_ex 141.87 V (015HK); DP/GP = 6.857 on every v3 session | `sTXSPEC.aRFPULSE[*]` |
| DP pulse duration | 1.28 (v3) vs 1.024 (v2); ratio 1.25; 6.857 × 1.25 = 8.571 = 30/3.5 EXACTLY → amplitude × duration reproduces 30° for a 3.5° gas pulse. INFERRED: adFree[4] = DP pulse length (ms), GP pulse 1.024 ms, same shape. Pulse SHAPE not in the header (Steve's source needed). If hard/rectangular: first null ±781 Hz, gas is 3839 Hz away (218 ppm × 17.61 Hz/ppm) → sinc sidelobe few % → ≤ ~1° effective on the gas per DP pulse; if shaped, negligible | `sWipMemBlock.adFree[4]` |
| spectral block | 600 lines (alFree[11]) × 20 → 780 lines after the ×26/20 kluge, dwell 60 µs (alFree[12]) | — |
| readout | dwell 5 µs, 512 samples (510 after killpts 2) = 2.56 ms; FOV 350 mm; Larmor 17.612 MHz (1.5 T) | — |
| interleaves | 832 unique (26 × 32 Thomson sets); acquired gas ilv: 015HK 7017 (120.7 s), 014KS 8945 (153.9 s), 039CP 8876 (152.7 s), 016PG 11664 (200.6 s), 025JC 6980 (120.1 s) → 8.4–14 passes; one 4 s breath = 232 gas ilv = 0.28 of a full set | `lTotalScanTimeSec = 928` is the nominal max; stopped by hand |

v2 protocol `fa_spiral_dyn_fancy_v2_20230131` (2023 HC sessions): TR 14.9 ms (14.86 on 000LL 01-31), gas flip
3.0° (000SV, 000SK) / 4.0° (000LL), dissolved 30° (40° on 000LL 2023-01-31), DP/GP amplitude = 10.0 exactly, DP
duration 1.024, +218 ppm. 000KR (2023-04) = gas-only, TR 7.43 ms, 2.5°. 000LL 2024-02-02 = v3 8-ch coil,
TR 25.9 ms, 4°/25°, 2 DP per gas (alFree[9] = 2) — different protocol.

## 7. Healthy sessions — the pool is 4, not 10

cohort_map (notes/calspec_package_2026-09-16/cohort_map.csv) + inventory F120: HC on the v3 17.2 ms protocol =
015HK, 014KS, 016PG, 025JC only. **039CP is LTX (lung transplant), not healthy** — the PCA doc calls it healthy.

| session | cohort | protocol | gas SNR insp/exp | DP SNR insp/exp | period s | b44d recon | s recon (old) | p recon | breath-hold |
|---|---|---|---|---|---|---|---|---|---|
| 2024-08-19_015HK | HC | v3 17.2 ms 3.5°/30° | 38.0 / 22.6 | 26.7 / 20.6 | 4.02 | b44d/…/d/recon.mat | yes | yes | MID00111 GOOD, T1eff 5.77 s |
| 2024-08-30_014KS | HC | v3 | 38.2 / 20.7 | 46.8 / 30.1 | 4.94 | yes | yes | yes | MID00118 GOOD, 5.32 s |
| 2024-08-23_016PG | HC | v3 | 21.5 / 9.2 (max bin 4) | 16.5 / 7.7 | 6.71 (IQR 4.8–9.4, irregular; 200 s scan, lungvox 154k) | yes | yes | yes | MID00273 GOOD, 4.52 s |
| 2024-11-13_025JC | HC | v3 | 31.8 / 21.5 | 26.2 / 25.3 | 3.58 (3.58–4.01) | yes | yes | yes | MID00063 GOOD, 4.80 s |
| 2024-09-24_039CP | **LTX** | v3 | 38.5 / 21.0 | 36.4 / 17.2 | 4.90 | yes | yes | yes | MID00085 GOOD, 4.51 s |
| 2023-05-18_000SV | HC_OLD | v2 14.9 ms 3.0°/30° | 28.3 / 13.9 | 16.7 / 11.3 | 5.05 | yes | yes | yes | none |
| 2023-10-27_000LL | HC | v2 14.9 ms 4.0°/30° | 31.4 / 24.2 | 23.0 / 23.7 | 3.26 | old tree only (AIkill_Dynamic/…/d) | yes | yes | none |
| 2023-01-31_000LL | HC | v2 14.86 ms 4.0°/40° | 29.7 / 11.1 | 21.1 / 8.7 | 4.75 (IQR 2.0–5.6) | yes | yes | — | none |
| 2023-05-15_000SK | HC_OLD | v2 3.0°/30° | 14.9 / 8.1 (weak) | 6.1 / 5.2 | 3.59 (1.8–4.8) | old tree only | yes | yes | MID00103 GOOD, 5.53 s |
| 2023-03-02_000LL | HC | v2 | — | — | 3.88 (1.5–4.5, irregular) | yes | — | — | — |
| 2023-02-17_000LL | HC | v2 | — | — | 1.33 (fast-breathing session) | yes | — | — | — |

## 8. Recon-free per-breath k0 (helpers/_bin_timing_export.py → outputs/bin_timing_2026-10-26/)

k0 (Σ|first 8 samples| per gas line, from the old-tree `signal_pneumo.npz`, affine scan-normalised) = coil-weighted
total magnetisation = `M_L` without mask / Rician / registration. Per-breath peak/trough, median (IQR):
015HK 4.41 (4.09–4.81) · 014KS 4.59 (4.41–4.92) · 039CP 4.84 (4.26–5.53) · 025JC 2.80 (2.62–2.92) ·
016PG 5.9 (4.5–12) · 000SV 5.9 (5.4–6.5) · 000LL 10-27 2.75 · 000LL 01-31 10.0 (v2 4°/40°).
Kernel cost: 015HK kernel-weighted k0 in SIGNAL bins = 0.749/0.195 = 3.84 vs per-breath 4.41 (σ = 1 bin smooths
~13 % of the swing); PCA's b44d mask sum 3.88. k0 peak at 1.3–1.4 s after EE vs navigator plateau 1.75–2.1 s →
Δt (B before C) survives in raw time. Files: `<study>_ilv.csv` (per gas interleave: t, k0 gas/diss, SIGNAL and
DIAPHRAGM bin fraction, nav_volume, sub-voxel nav z), `<study>_bins.csv` (per bin: n_hard, w_eff, t after EE,
kernel-weighted nav z / nav_volume / k0), `<study>_breaths.csv`.

## 8b. RAW k0 (helpers/_bin_timing_rawk0.py, conda base python, `<study>_{ilv,bins,breaths}_raw.csv`)

Straight from d/input.mrd, samples 2:10, sum8 / coh8 / s2. Raw B/A medians: 015HK 4.47, 014KS 4.54, 039CP 4.94,
016PG 5.57, 025JC 2.84. Per breath also C = navigator max (k0_C, t_C), nav z at A/B/C, dissolved k0 at A/B/C.
- **D/G at k0 level at A**: 0.116 / 0.151 / 0.091 / 0.098 / 0.125 (015HK/014KS/039CP/016PG/025JC). Image-mask D/G
  (PCA, b44d EI mask) = 0.50–0.64 of that, every session: dissolved k-space apodised with dplb 40 vs gplb 300
  (results.py:336 → recon.py:98-99) → wide dissolved PSF → mask collects less. Magnetisation ratio = k0 ratio ÷
  sin30°/sin3.5° = 8.19 → 0.011–0.018.
- **Gas excited by the DP pulse** (calspec table pooled_rbc_gas / pooled_mem_gas): gas/(RBC+mem) = 0.067 / 0.037 /
  0.053 / 0.039 / 0.047 → 4–10 % of the dissolved k0 over the cycle; θ_gas ≈ 0.03° → no XTC on the gas pool from
  the DP pulses; the DP pulse is shaped (hard-pulse sidelobe estimate in §6 retracted). alFree[13] = 43 = Spec_ex
  flip (141.9/99.0 V = sin43/sin30), alFree[3] = 30 = DP flip, adFree[4] = DP pulse ms: all consistent.
- **Drift**: pneumo tidal amplitude trend 015HK −13.7, 014KS +17.1, 025JC +7.6 %/min (own clock, crude) vs k0 trend
  −36 / +24 / −4 %/min → 014KS's rising k0 = deeper breathing, not dose. 039CP/016PG pneumo detrend failed.

## 8c. Audit of the per-breath picks (PCA budget §4d)

- C = argmax of the piecewise-linear nav z (one frame centre, 0.447 s grid). 95 % plateau columns added
  (t_C95_start/mid/end); C/B moves ≤ 1 % except 016PG (+5 %, 1.7 s plateaus).
- Rician negligible: trough = 42–64× the sum8 noise floor; sum8/coh8 = 8.00 at A and B (fully coherent 8 samples).
- Real bias = extreme-value picks (A = min, B = max of a 1.4–6 % scatter series). Savgol(11,2)-smoothed picks
  (`*_smsm` columns): B/A 4.34 / 4.41 / 4.78 / 5.17 / 2.76 (015HK/014KS/039CP/016PG/025JC), C/B 0.955 / 0.910 /
  0.937 / 0.789 / 0.965, nextEE/C 0.236 / 0.251 / 0.232 / 0.224 / 0.375 — ρ up ~2.5 %, m_L(A) up ~3 % vs raw picks.
  Drop breaths with period > 1.5× median (missed/merged EE).
- Outside-5σ-mask gas on b44d (Σ real(gas) whole volume, unbiased): 1.5–2.5 % at B, 5–13 % at bin 0 (mask shrinkage).
  FOV 350 mm isotropic, trachea inside, mouthpiece/bag outside the coil's sensitive region.

## 9. Old s/ tree dissolved ≠ b44d dissolved (PCA found D/G 1.3–2.9× higher on old s/)

Old image (472fbc9/40d23a4) split with Steve's old basis angle 2π·df·TEeff (dphi_old −17° 015HK / −36° 014KS /
−8° 039CP → gain 1/|sin| 3.4 / 1.7 / 7.5), stored real/imag = aRBC/aTP → `dissolved_phase_magnitude` =
|aRBC + i aTP| ≠ |z|. Inversion (linear): |z| = |real·exp(i·dphi_old) + imag|, dphi_old per session in
`outputs/resplit_2026-09-24/resplit_summary.csv`. Deapod and calcb low-pass are ratio-neutral.

## 10. Prepared, NOT run (Hooman's go): SIGNAL rerun on the b44d image + bin-export fork patch

`pipeline/run_sig_b44d.sh` (GO=1 gate, Ext ≥ 40 GiB, compute_mode) → `AIkill_Dynamic_b44d_sig/<key>/s/`, image
262b16f (= b44d; the default `recon_codespec.yml` still pins d136eb1 — never use it bare). 10–20 min wall and
~1 GB per session. `notes/ForkPatch_BinExport_2026-10-26.md`: export_nav flag + ilvbin / ee_times / per-bin
kernel-weighted time and nav z in every binning mode (needs fork commit + image build).

Paths: `/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/AIkill_Dynamic_b44d/<session>/d/recon.mat` (production),
`…/AIkill_Dynamic/<session>/{d,s,p}/recon.mat` (first batch). Breath-hold MIDs = same dual-phase sequence, graded in
`notes/calspec_package_2026-09-16/session_table.csv` (hold_grade, hold_T1eff_s = global T1* of the hold).

## 11. Image phase convention (answer to PCA, 2026-10-26 evening) — same for b44d d/ and 262b16f s/

- `gas_phase[bin] = Σ_ch Re(F_bin,ch · b_ch)` (`results.py` dyn_recon). `b` = conj(whole-scan pooled gas image)/noise per
  channel, phase inside bmask (pooled > 10× corner noise) Gaussian low-passed σ 4.4 vox (~5° rms removed), outside bmask
  phase = quadratic polynomial extrapolation, then `b /= max_ch|b|` per voxel → 1-ch sessions: |b| ≡ 1, pure phase
  reference. gas_phase = real projection onto the pooled-scan phase; imag discarded. Noise = zero-mean Gaussian
  (σ_real = σ_cplx/√2). Σ_mask gas_phase unbiased, no Rician correction. `gas_phase_magnitude` = |F·b| (1-ch, QC) is Rician.
  Caveat: pooled phase reference → per-bin cos(Δφ) loss where bin phase ≠ pooled phase (rim outside bmask, 2steve/06).
  Diagnostic: Σ gas_phase / Σ gas_phase_magnitude over |.|>10σ voxels per bin ≈ 1.
- Dissolved: `thisrspace = F_dp,bin · b`; 1-point Dixon per bin: global ph_bin s.t. masked RBC/TP = spectroscopic target;
  stored real = aRBC, imag = aTP with `thisrspace·e^{i ph} = aRBC·e^{i dphi} + aTP` (oblique basis). z = real·e^{i dphi} + imag.
  dphi (b44 d/tyger.log, 16/16 bins split): 015HK 53.5°, 014KS 38.5°, 039CP 53.3°, 016PG 53.2°, 025JC 51.6°.
- Global total dissolved = Σ_mask(aRBC + aTP) = Σ_mask(real + imag): linear, unbiased, phase-spread independent.
  |Σ z| under-counts (in-lung phase spread up to dphi from RBC-fraction variation); Σ|z| over-counts (Rician).
- SIGNAL rerun launched 19:15 (`run_sig_b44d.sh`, 5 sessions, Ext `AIkill_Dynamic_b44d_sig/<key>/s/`), Tyger run 1360 = 015HK.

## 12. b44d_sig SIGNAL rerun DONE + Tyger's exact membership + auto exclusions (2026-10-26, 19:15–19:29 run)

- `run_sig_b44d.sh` ran 5 sessions (015HK, 014KS, 039CP 09-24, 016PG, 025JC): 2–3.5 min each (upload 6 s, Tyger run
  ~2 min). Tree `Ext/AIkill_Dynamic_b44d_sig/<key>/s/`. recon.mat re-exported with `pipeline/mrd_to_mat.py` (now takes a
  dir holding output.mrd and merges existing keys) → carries calcb_b, gas_phase_magnitude, rbc_tp_* + post_process keys.
- Scale: `calcb_b` s/ vs d/ max|diff| 3e-7 (|b| ≤ 1) → identical pooled phase reference, absolute scale equal.
- **Exact membership**: `helpers/_tyger_bins_local.py` runs the fork's `raw.load_from_arr` on `<s>/input.mrd` offline
  (CPU) → `<key>_ilvbin_tyger.csv` (ilvbin_frac, excluded, k0_tyger_norm), `<key>_ee_tyger.csv`. Matches my raw-k0
  derivation (`bins_raw.csv`) exactly on all 5; the fid-affine derivation (`bins.csv`, from post_process fid_signal =
  samples 0–7 incl. ADC transients vs Tyger killpts=2) is off on 014KS (33 % of ilv), 025JC (13 %), 016PG (8 %).
- **Auto exclude ranges** (raw.py:432-445: whole 832-ilv sets with Σk0 < mean/2; kernel skips zeroed samples →
  dropped from every image incl. d/): 015HK 100–129 s; 014KS 14–29 (set 1 = 9 % of mean, no xenon) + 143–158;
  039CP 143–158 (none binned); 016PG 186–215; 025JC 100–129. Binned ilv lost 1055/1207/0/705/995. Not printed in
  tyger.log unless manual ranges are given.
- Flatness proof Σ_FOV gas_phase(bin)/k0_wmean(bin) with exact membership: ±1.3 / ±2 / ±2 / ±1 / ±0.6 % (PCA's 014KS
  ramp 0.81→1.14 was the mis-split breath in bins.csv). Residual ±2 % = |k0| (magnitude) vs coherent DC at troughs.
- New per-bin table `helpers/_bins_from_tyger.py` → `<key>_bins_tyger.csv` incl. `t_after_anchorEE_wmean_s` (kernel-
  weighted real time relative to the EE preceding the bin's phase centre, wrap members folded; ≈ b/16 × median T).
- RBC/TP per bin: `rbc_tp_R` = target ± 0.003 in all bins/sessions → ratio fixed by construction; A→C changes seen by
  PCA with a fixed EI mask are mask effects. Total dissolved = Σ(aRBC + aTP).
- Bin 0 vs bin 15 (k0 bin15/bin0 = 0.90 on 014KS with exact membership): soft kernel + steep post-EE rise vs flat
  pre-EE plateau; A on signal bins = min(bin 0, bin 15) or per-breath k0_A.
- Proposal discussed (not run): absolute-time-since-EE binning in raw.py (T_ref = median period) + converter flag.

## 13. DIAPHRAGM-bin exact tables (PCA ask, 2026-10-26 late)

- `helpers/_bins_diaphragm_tyger.py` → `<key>_bins_diaphragm_tyger.csv`, `<key>_ilvbin_diaphragm_tyger.csv`. Uses the
  ROOT fork's `raw.bin` on exported `nav_volume` (= post-rescale `ilvvol[DIAPHRAGM]`, exactly what Tyger binned).
  My `_bin_timing.steve_bin` copy agrees to 0.0000 bins on all 5 sessions.
- **Trap (014KS)**: mid-scan auto-excluded window 14–29 s is linearly interpolated in `nav_volume` (np.interp inside
  the measured range), so its 873 interleaves get diaphragm bins 2–6 while the gridding kernel skips them (zeroed).
  Per-bin kernel averages must zero those weights (flag from `<key>_ilvbin_tyger.csv`). With the flag: Σ_FOV image /
  k0_wmean flat 0.977–1.019 (was 0.87–1.43). End-of-scan exclusions are NaN in nav (outside measured range) → fine.
- Diaphragm bins 0–2 and 14–15 pool the exhale plateau across breaths: t_after_EE wsd 1.2–1.8 s vs 0.3 s mid-bins;
  bin 0 mean time ≈ 2.6 s (015HK), not EE. Physical axis for those bins = nav_z / nav_volume, not time.
