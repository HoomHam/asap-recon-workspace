---
tags: [learn]
status: VALID (2026-09-13)
---
# RBC/TP split — how we knew it was broken (noise diagnostics)

**Context.** While adding RBC and TP SNR columns (`helpers/snr_calc.py`), many "split" sessions
gave a negative mean RBC or TP inside the lung. A negative amplitude can't be physical. The four
checks below turned that suspicion into proof, and led to a dissolved SNR that doesn't depend on
Steve's split. Steve's side of the story: `../../2steve/04_RBCTP_Split_PhaseStop_Conditioning.md`.

## The setup in one line
Steve stores the dissolved image as `aRBC + 1j·aTP` (real = RBC, imag = TP). He gets them by
inverting a 2×2 system per voxel, after rotating the image by a phase `ph` from a sweep:

    Re(z) = aRBC·cos φR + aTP·cos φT
    Im(z) = aRBC·sin φR + aTP·sin φT        with Δφ = φR − φT = 2π(fRBC − fTP)·TEeff

## Test 1 — negative amplitudes (the first flag)
- 045VS 2024-09-10, end-inspiration bin 6: lung mean aRBC **+6202**, aTP **−4226**; aTP < 0 in 87 % of lung voxels.
- Across sessions the negative part varies: in some it's TP, in others RBC. The aRBC and aTP maps look
  like mirror images of each other.

## Test 2 — did the "solve" actually hit its target?
- `tyger.log` prints `RBC/TP phase solved at ph=… (R=… crossed target RBCTPratio=…)`.
- 045VS logs **R = 37 against a target of 0.293**; 034LR logs R = 227 against 0.661.
  "Solved" was printed, but R never came near the target.
- Over all 86 split sessions: 19 on target in every bin, **60 in no bin**, 7 mixed.
- Why: the stop test `(R−target)(lastR−target) < 0` also passes when ΣaTP crosses zero
  (R jumps −∞ → +∞). And R uses whole-volume sums, not the lung.

## Test 3 — the noise that didn't add up (std of the parts vs std of the sum)
In the background corners, where there is only noise:

| | std |
|---|---|
| aRBC | 341 |
| aTP | 346 |
| **aRBC + aTP** | **37** |

If the two noises were independent, the sum would be **√2 ≈ 1.4× noisier** than each part
(about 485). Instead it is **9× quieter**. That can only happen if the two noises are nearly
exact opposites of each other.

## Test 4 — the correlation itself
- corr(aRBC, aTP) in the background: **−0.994** (045VS), −0.981 (030DN), −0.778 (034LR).
  Panel (b) of the figure is a straight line.
- **Theory.** Complex noise in z is isotropic, σ per component. Push it through the 2×2 inversion:
  - corr(aRBC, aTP) = −cos Δφ
  - std(aRBC) = std(aTP) = σ / |sin Δφ|
  - std(aRBC + aTP) ≈ σ
- **Check.** For 045VS, corr −0.994 gives Δφ ≈ 6°, so 1/sin 6° ≈ 9.1. The measured
  341 / 37 = 9.2. The theory and the measurement agree.
- **Meaning.** RBC and TP were nearly in phase at this TEeff. A one-point Dixon split needs about
  90°; at 6° the inversion is almost singular. Each "component" is 9× amplified noise that cancels
  in the sum, which is why the maps mirror each other.

## Test 5 — the "2": proving the stored pair is just a linear map of the complex image
- **Idea.** The stored v = (aRBC, aTP) = M·(Re z, Im z) for some 2×2 matrix M (M = identity if
  unsplit). Background noise of v has covariance C = σ²·M·Mᵀ. So per voxel:

      q = vᵀ C⁻¹ v = |z|² / σ²          (the rotation part of M cancels; only |z| survives)

- **Why 2.** In pure noise, whitening gives two independent unit-variance Gaussians (real and
  imaginary). q is then a sum of two squared standard normals, χ² with 2 degrees of freedom, so
  **E[q] = 2**.
- **Measured** background mean q: **2.04** (045VS), 2.07 (030DN), 2.03 (034LR), 2.07 (007IT,
  unsplit control). Per-bin range 2.00–2.22.
  → Confirms (a) the stored pair really is a linear transform of isotropic complex noise, and
  (b) the covariance estimate is right. So the whitened magnitude is trustworthy.
- **Rician correction.** |z|²/σ² has a noise bias of exactly that 2, so the corrected magnitude
  SNR per voxel is √max(q − 2, 0). Mean over the lung mask gives the dissolved SNR.

## Test 6 — control on an unsplit session
- 007IT (unsplit, no RBC/TP fit) had its dissolved SNR computed two ways:
  - plain Rician magnitude: EI 15.3
  - whitened magnitude: EI 14.7
- Same answer, as expected when M = identity, so the new method doesn't invent signal.
- On split sessions it rose above the old aRBC + aTP value (045VS 52.9 → 57.7; 034LR 24.4 → 27.8).
  The sum is a projection at a wrong angle and loses signal.

## Takeaways
1. A **negative mean amplitude** is a flag. **std(sum) ≪ std(parts)** plus **corr ≈ −1** is the proof:
   the split is ill-conditioned noise amplification.
2. **corr = −cos Δφ** and **gain = 1/|sin Δφ|** let you read the RBC–TP phase difference straight
   off the background noise, with no spectra needed.
3. **E[vᵀC⁻¹v] = 2** is a clean self-test that any stored 2-channel image is an honest linear map
   of complex Gaussian noise.
4. For most cohort sessions, `dissolved_phase_real/imag` in `recon.mat` are **not** usable RBC/TP.
   Use the whitened magnitude for total dissolved. Check `rbc_tp_solve_check.csv` before any RBC:TP work.

## Where things are
- Figure: `outputs/snr_2026-09-13/04_rbctp_split_evidence.png` (script `helpers/_rbctp_fig.py`)
- Per-session solve validity: `outputs/snr_2026-09-13/rbc_tp_solve_check.csv`; Excel col S
- SNR code (whitened magnitude): `helpers/snr_calc.py`
