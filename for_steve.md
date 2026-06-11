# ASAP Reconstruction — Report for Steve

**From:** Hooman
**Date:** 2026-06-11
**Purpose:** I spent the last days reading your reconstruction code line by line, reimplementing it, and validating it against an independent reference on real phantom data. This report lays out what I believe your implementation does and why, what we measured, and a claim I'd like you to attack: **the GPU is not necessary — the same recon, same math, runs in comparable time on a laptop CPU.** Please play devil's advocate; I'll respond to every objection.

---

## 1. What I understand your implementation to be

Files read: `recon.py`, `raw.py`, `results.py`, `gtypes.py`, `main.py` (plus `asap/asap.c` for orientation). The recon core, as I read it:

**Gridding (`cudarecon`):** each sample is distributed to a 4³ box of grid cells around its trajectory point with Gaussian weights `exp(-d²/0.2)`. With σ ≈ 0.32 cells, the weight one cell away is already e⁻⁵ ≈ 0.7%, so effectively each sample votes into its nearest cell with slight sub-cell smoothing. The Kaiser–Bessel machinery (`bessi0`) is present but commented out.

**Normalization (`cudarenorm`):** every cell is divided by its accumulated weight sum (`knorm`). This is the step I think is under-appreciated: it turns the gridder into a **kernel-regression estimator of k-space** —

```
k̂(cell) = Σᵢ C(kᵢ − cell)·sᵢ / Σᵢ C(kᵢ − cell)
```

— a per-cell *average* rather than a sum. That single design choice is why you need no explicit density compensation: averaging is self-normalizing at cell resolution, the spiral center's thousandfold oversampling simply averages down, and any zeroed sample (spike, exclusion) drops out with the bookkeeping done automatically. It also explains why the narrow Gaussian is not a defect but a requirement: widen the kernel under this normalization and neighboring samples bleed into each cell's average, low-passing k-space. Narrow-kernel + normalize is one self-consistent architecture; wide-KB + iterative DCF + plain sum (Faraz's) is the other. The kernel width and the density handling are a single coupled choice — they can't be mixed and matched.

**Everything else I traced:** killpts=2 dropping the ramp samples; the per-channel noise normalization via Gaussian fit to the Re(raw) histogram; FID rephasing to zero phase; the 10× spike filter against cross-interleave statistics; the readout filter `exp(-(t/gplb)²)` inside the kernel; soft binning `exp(-Δbin²/2)` letting every interleave contribute to every bin; the b-map coil combine with the polynomial phase fill in low-signal voxels; forward-FFT convention; MS=240 → crop IS=100 (oversampling 2.4 baked into the trajectory rescale). I wrote all of this up in detail — happy to share the full comparison document.

If I've misread any of this, that's the first thing I want to hear.

## 2. What we built to test it

1. **A faithful CPU re-implementation** of `cudarecon` + `cudarenorm` in vectorized numpy — same filter, same asymmetric 4-cell box (`range(cx-2, cx+2)`), same Gaussian, same knorm divide including the 1e-5 floor, same F-order reshape, FFT, and crop. Not bit-identical (float64 vs float32 atomics, summation order), but algorithmically line-for-line.
2. **An independent reference reconstruction** with a different architecture entirely: FINUFFT (Flatiron's non-uniform FFT library, kernel accuracy ~1e-9) as the forward/adjoint operator pair, and a conjugate-gradient least-squares inverse on top. The operator pair passes the adjointness dot-product test at machine precision (6×10⁻¹⁵). This recon shares *nothing* with your code path — different kernel, different density handling (none — the solver handles it), different library.
3. **Test data:** the v3 FOV250 ACR phantom acquisition (512×832, gas-only pattern, single channel) with the measured trajectory calibration. Inputs were produced by *your own* `raw.py`/`traj` loaders running headless on CPU — your noise normalization, rephasing, spike filtering, and npts inference all ran unmodified (your periodicity-based npts inference correctly found 510, and the gas-only pattern detection worked first try).

## 3. Results

**Your recon agrees with the rigorous inverse at correlation 0.984.** The CG least-squares reconstruction and the faithful re-implementation of your kernel produce nearly identical volumes on the phantom. All ACR structures (resolution inserts, grid patterns, rings) resolve equally in both.

**Your geometry is physically correct.** Measured object extent from your recon: 178×180×150 mm against the ACR phantom's true 190 mm diameter × 148 mm length (the in-plane numbers read slightly low due to thresholding; the length is within 2 mm). I note this specifically because the same dataset reconstructed through another path in our group came out magnified ×1.205 due to a k-rescaling step in the analysis script — your pipeline was the one that matched physical truth, and our independent recon confirmed it. The arbiter role works.

**Interpretation I'd defend:** for fully-sampled static data, your normalized-average gridder is a very good approximation of the true least-squares inverse — which is exactly what the kernel-regression view predicts. The textbook objections to Gaussian-kernel gridding (aliasing sidelobes, rolloff, no de-apodization) are largely neutralized by the per-cell normalization plus the 2.4× oversample-and-crop. Where I'd expect the approximation to genuinely degrade is heavily *undersampled* bins, where cell-resolution density compensation and the narrow kernel's lack of interpolation both start to cost. That's a testable prediction and on our list.

## 4. The claim to attack: the GPU is not necessary

### 4a. The cost is architectural, not physical

Counting operations in `results.dyn_recon`: for each of the 16 bins × each channel, the kernel traverses the **entire acquired dataset** — every repetition of the trajectory over the whole scan — because the soft-bin weight is evaluated per sample inside the kernel. For a scan where the trajectory repeats R times, that's `16 × nch × R` passes over one trajectory's worth of samples, with 64 cells × 3 atomic adds each. The atomics also serialize hardest exactly where a spiral has the most samples — k-center. None of this cost buys image quality; it buys implementation uniformity (one kernel for static, binned, sliding-window, T1RF — which I agree is genuinely valuable).

The same output is computable in **one pass**: accumulate 16 per-bin (k, knorm) grid pairs simultaneously (16 × two 240³ arrays ≈ 3.5 GB, fits in RAM), adding each sample's contribution `binwt(i,b)·wt·sᵢ` to each bin's grids. Mathematically identical result — the sums are just reordered. Additionally, `binwt = exp(-Δ²/2)` is below 10⁻³ for |Δ| > 3.7 bins, so ~60% of the bin-sample pairs can be skipped with no visible change.

### 4b. Measured numbers (Apple M4 Pro laptop, CPU only, no CUDA anywhere)

The vectorized re-implementation processes your per-sample loop as 64 array-wide scatter-add operations — Python never touches individual samples.

| Operation (10M samples ≈ 3-min dynamic scan, 70 MB .dat) | Measured |
|---|---|
| One (bin × channel) pass of your kernel, 240³ grid | **22.4 s** |
| FINUFFT adjoint of the same 10M samples → 100³ | **0.3 s** |

Projected full recons from those measurements:

| 16-bin dynamic recon | 1 ch | 8 ch |
|---|---|---|
| Your architecture, verbatim (16 full re-passes) | ~6 min | ~48 min |
| Same math, one-pass restructure | ~35 s | ~4–5 min |
| FINUFFT adjoint per bin | ~10 s | ~1 min |
| FINUFFT CG inverse per bin (15 iterations) | ~2.5 min | ~20 min |

The phantom recon end-to-end (your loaders → your kernel → comparison figures) runs in about a minute on the laptop.

### 4c. What I am *not* claiming

- Not that the GPU version was a mistake — in pure-Python-loop form this recon would take days, and CUDA was the right 2-line-decorator escape hatch. The observation is that **vectorized CPU is the same escape hatch**, and on modern laptop silicon it lands within the same order of magnitude as the GPU while removing the CUDA dependency entirely (relevant for Kento's deployment question too: the cloud GPU node may be solving a problem that no longer exists).
- Not bit-exactness. Float32 atomics vs float64 ordered sums will differ at the ~1e-6 level. A one-time certification run of your GPU output against the CPU re-implementation on identical input would close that — I'd like to do it.
- Not that this holds at any scale. At 100× the data, or for the iterative CS reconstruction we're building next, GPU acceleration becomes attractive again — but as an optimization, not a requirement.

## 5. Questions where I'd value your advocacy

1. **The kernel-regression reading of knorm** (§1) — is that how you think about it, or am I retrofitting theory onto a pragmatic choice? Specifically: was `kdist0sq = 0.2` tuned (against what criterion?), or chosen as "small enough to be nearest-cell"?
2. **gplb = 300** — derived from Xe T2* and the dwell time, or empirical? It dominates the effective PSF, and for our cross-implementation comparisons we need to either zero it or replicate it exactly.
3. **killpts = 2** — ADC settling, gradient ramp distrust, or both? Faraz keeps all samples; the difference lands at k-center. We plan to test both ways.
4. **The disabled Hermitian symmetrization** in `cudarenorm` (`oldway = 1`) — abandoned, or waiting for something?
5. **The knorm < 1e-5 → 0 kluge** — empty cells become k-space zeros rather than interpolated values. At 2.4× oversampling with full sampling this looks harmless; do you see a regime where it bites?
6. **`dyn_usimg_recon` omits the inner `ifftshift`** that `calcb`/`dyn_recon` have, and flips axis 0. Intentional (cancels elsewhere) or vestigial?
7. **Where would you bet your recon beats the CG inverse?** Undersampled bins is where I'd bet against it — but you know failure modes of this data I don't.

## 6. Where this is going

Next step is a compressed-sensing reconstruction: `min ‖A·ρ − s‖² + λ·R(ρ)` with the FINUFFT operator as A, your soft-bin weights surviving as a diagonal weighting in the data term, and temporal regularization across bins replacing data-domain bin sharing. Your implementation stays load-bearing in that plan in three ways: as the validated loader/preprocessing front end, as the fast preview recon, and as the baseline any CS result must beat.

All code, the full Steve-vs-Faraz comparison document, and the measurement scripts are in my workspace repo — happy to walk through any of it.

— Hooman
