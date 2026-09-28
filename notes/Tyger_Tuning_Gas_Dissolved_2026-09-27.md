# Tyger tuning, step 2 — grid/kernel/LB knobs against a known object

**Status:** RESULT (Opus 5.5, session 14, 2026-09-27; revised 2026-09-28 after the SNR / soft-bin / registration checks, F96–F98) · recommendation = deapod ON, everything else unchanged · facts F90–F98 · card C54
Step 2 of `notes/Two_Targets_Registration_vs_Truth.md`, for the T (truth) image. Subjects: 2024-11-13_025JC
(healthy, trust A) and 2024-10-02_011CN (EBV, trust B); 015KH was never reconstructed, 002BB has no cal block.

## Method (C54 `helpers/tune_gas.py`)
Known lung object (all-data Steve image, zeroed outside a 12-vox lung dilation) sampled with finufft at the REAL
unique trajectory; real per-bin soft weights, exclude mask and EPS keep rule; noise per sample calibrated to the real
far background. Steve's normalized-convolution gridding collapsed exactly onto the unique samples (production
rebuild relerr 4e-7 gas, corr 1.000000 dissolved). Per setting: blur (all data), per-bin hole aliasing, noise, total
in-lung rms error, PSF from delta probes. Gas bracketed with two objects (built without / with apod).

## Result
| knob | derivation / measurement | verdict |
|---|---|---|
| IS | FOV 350 mm hardcoded → pixel 3.5 mm; kmax 40.4 Δk → native ≈ 4.3 mm; IS 100 = ×1.24 oversampled | keep 100 |
| MS | 240 → 200: holes 0.56 → 0.45, gas total 0.87× prod with deapod on these two SNR-20 sessions; dissolved neutral. BUT SNR-conditional (F97, 2026-09-28): crosses 1.0 at SNR ≈ 9–12, 1.11–1.17× at SNR 5; 23/61 sessions have bins < 10 | **keep 240** (revised) |
| deapod | divide by the Gaussian kernel rolloff (0.92 at the crop face): 0.88–0.95× on every setting, gas and dissolved; ratio maps unchanged | **on** |
| kdist0sq | wider than σ ≈ 0.2 Δk always worse (k-space smoothing bias > hole savings) | keep 0.2 |
| gplb | optimum bracketed 300–450 by the two objects; flat | keep 300 |
| dplb | F92: e⁻¹ at n = 40 = 24° RBC−TP drift; ≈ 27–30 mm PSF. It is the one-point-Dixon limit, not an SNR knob | keep 40 (see below) |

Gain from knobs: ≈ 10–18 % lower total in-lung error. Eye (real bin 0, `outputs/tune_gas/candidates_bin0.png`):
subtle; MS200+deapod ≈ prod with a slightly brighter periphery (lung mean +4–5 %), MS160 visibly noisier.

## The bigger lever (F91)
On the same samples, noise-free, Steve's gridding misses the object by 8–15 % while weighted-LS CG gets 2.2 % (all
data) / 3.8–3.9 % (bin 0). With noise, plain CG is no better than Steve (13.5 % / 15.2 % best). My first draft called
a DCF-preconditioned regularised solve "the predicted ~2×": RETRACTED 2026-09-28 (F95). XeCS's record (F195/F201/
F202) shows no gain over the best windowed linear arm under any tested configuration (CG+DCF, wavelet/TV CS, 4D CS);
Steve's gridder is the SNR-best linear arm at equal FWHM (1.7×, ACR). The noise-free floor is real; converting it is
an open question, not a pending gain.
For dissolved resolution, the equivalent lever is putting the RBC−TP drift into the forward model (two species with
known Δf) instead of raising dplb.

## Not covered
Motion (static object → bindist0sq × nbins untouched, step 4); calcb b per setting; the R (registration) image;
more than two subjects; a Tyger run with the new settings (step 3 would add MS/deapod to the fork's spec).
