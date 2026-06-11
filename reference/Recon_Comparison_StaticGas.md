# Steve vs Faraz — Static Gas-Phase Reconstruction Comparison

**Date:** 2026-06-10
**Scope:** Single bin, static phantom, gas phase only. Excludes: dissolved phase / Dixon, temporal binning strategies, dynamic recon.
**Method:** Fresh line-by-line read of both codebases (not derived from previous comparison markdowns).

**Sources read:**

| Side | Files |
|------|-------|
| Steve (Python/CUDA) | `recon.py`, `raw.py`, `results.py`, `gtypes.py` |
| Faraz (MATLAB) | `fa_spiral_dyn_recon.m`, `gridrecon_fa_20230113.m`, `grid_lookup_20230113.m`, `createKBkernel.m`, `iterative_dcf_fa_20190910.m`, `combinecoils_fa.m`, `loadtrajectory3D.m`, `spiral_human_20240227.m` |

---

## Pipeline at a glance

| Step | Steve | Faraz |
|------|-------|-------|
| Trajectory source | `.npy` file × FOV → delta-k units (`traj.load`, raw.py:30) | Calibration `.mat`, matched by freq/FOV/nsamples/nleaves/nreps (`loadtrajectory3D.m`) |
| Trajectory units at gridding | Grid-index units: `k·MS/IS + MS/2` (raw.py:68) | Physical 1/mm; grid spacing `1/(fov·os)` |
| Initial-point kill | First 2 points per interleave dropped (`killpts=2`) | None |
| Channel noise normalization | Raw domain: Gaussian fit to histogram of Re(data) per channel (raw.py:211-215) | Image domain: corner noise region in `combinecoils_fa` |
| Noise spike removal | Zero points > 10× mean of that readout index across interleaves (raw.py:416-427) | 3σ over moving average (window 3) across interleaves; replaced by neighbor mean, not zeroed (spiral_human:136-151) |
| Global phasing | FID start rephased to zero phase per channel (mean of first 4 pts; raw.py:257-260) | None at raw stage (phase handled in coil combine); optional FOV-shift phase ramp |
| Effective oversampling | MS/IS = 240/100 = **2.4** | `os` = **3** |
| Kernel | **Gaussian** `exp(-d²/0.2)`, d in grid cells (recon.py:124-129) | **Kaiser–Bessel**, width 5 cells, Beatty β ≈ 12.78, 10000-pt lookup table |
| Kernel support | Box `cx-2 … cx+1` per axis (4 cells, asymmetric) | 25 nearest grid points (radius ≈ 1.8 cells) |
| Density compensation | **None explicit.** Per-cell divide by accumulated kernel weight (`knorm`) | **Iterative Pipe–Menon DCF** (5 iter, `sqrt(KB)` kernel) pre-applied, *plus* per-cell divide by Σ(wi·kernel) |
| De-apodization | None | None (the "partial de-apodization" divide is weight normalization, not kernel-FT correction) |
| FFT | `fftshift(fftn(ifftshift(k)))`, forward FFT | Per-dim `fftshift(fft(fftshift(X,i),[],i),i)`, forward FFT |
| Crop | Central IS=100 from MS=240 | Central `NumK·zfill` from oversampled grid |
| Coil combine | `Σ_ch real(img_ch · b_ch)`, b from full-data recon with polynomial phase fill (results.py:59-153) | Bydder MRM 2002: `real(b'·img)`, b = Σimg/Σ\|img\| (combinecoils_fa) |
| Output | Real image, zero-mean noise (negative background) | Real image, noise-normalized |
| Compute | numba.cuda kernels, atomic adds | MATLAB sample loop, vectorized channels, parfor knnsearch, grid lookup cache |

---

## Step-by-step detail

### 1. Trajectory

Both use **measured (calibrated) trajectories**, not analytic ones — same philosophy, different storage.

- **Steve** (`raw.py traj.load`): loads `.npy`, multiplies by FOV (350 mm) so coordinates are in delta-k units. Infers `npts` per interleave from the periodicity of |k|² via FFT (fragile if trajectory changes — known pitfall). Drops first `killpts=2` samples of every interleave (gradient ramp-up garbage). Rescales to grid coordinates `k·MS/IS + MS/2`, i.e. one delta-k maps to 2.4 grid cells.
- **Faraz** (`loadtrajectory3D.m`): looks up a calibration entry in `calibrations_3D_20220308.mat` matched on frequency (±5%), FOV (±5%), nsamples, nleaves, nreps, imgsize, orientation, nucleus. Returns `[kx ky kz]` in physical units. No initial-point kill.

**Difference that matters for static phantom:** Steve discards 2 points/interleave near k-space center (highest signal); Faraz keeps them. Steve's center-of-k-space data is otherwise re-phased to zero (below), so the dropped points are partially compensated by design.

### 2. Raw preprocessing

- **Steve:**
  1. Per-channel noise normalization: histogram of Re(raw), Gaussian fit, divide by fitted σ (raw.py:211-215).
  2. Global rephasing: each channel multiplied by `conj(mean(first 4 points of all gas FIDs))/|·|` so FID starts at zero phase (raw.py:257-260).
  3. Noise spikes: any sample > 10× the mean magnitude of that readout index across all interleaves → zeroed; the CUDA kernel later skips zeros (`abs(rawval) < eps → continue`).
  4. Readout apodization inside the gridder: `fwt = exp(-((idx % npts)/smoothing)²)` — a Gaussian filter along the readout (line-broadening analog). For the full gas recon `smoothing = gplb = 300` points.
- **Faraz:**
  1. Noise spikes: per readout-position moving average (window 3 across interleaves), threshold 1 + 3σ; bad points replaced by mean of neighbors (interpolated, not zeroed).
  2. No raw-domain channel normalization, no global rephasing, no readout apodization. Channel normalization deferred to `combinecoils_fa` (image-domain corner noise).
  3. Optional FOV-shift via linear phase ramp on raw data (`shiftk`).

**Difference:** Steve filters (apodizes) k-space along the readout — mild resolution loss for SNR; Faraz grids unfiltered data. Spike handling: zero-and-skip (Steve) vs interpolate (Faraz) — for Steve a zeroed sample simply contributes nothing (normalization absorbs it), so the two behaviors are closer than they look.

### 3. Grid geometry

- **Steve:** fixed MS = 240 oversampled grid, IS = 100 final image → oversampling 2.4. Grid index space 0…239, center 120.
- **Faraz:** grid spans `±ceil(zfill·imgsize)/2/fov` in steps of `1/(fov·os)` with os = 3, zfill = 1, imgsize = matsize from protocol header. Oversampling 3.

Both crop the central target-matrix region after FFT; oversampling pushes gridding-aliasing sidelobes outside the crop. Comparable approach, slightly different factor.

### 4. Gridding kernel — the largest conceptual difference

- **Faraz:** true Kaiser–Bessel. `createKBkernel.m`: `besseli(0, β·sqrt(1−(2u/w)²))`, normalized to peak 1. β from Beatty et al. Eq (5): `π·sqrt(k²/os²·(os−0.5)² − 0.8)`; with k = 5, os = 3 → **β ≈ 12.78**. Width = 5 grid cells. Kernel evaluated by linear interpolation of a 10000-entry table at the **radial Euclidean distance** to each neighbor (radially symmetric application, not separable per-axis).
- **Steve:** Gaussian `exp(-dsq/0.2)` on radial squared distance in grid cells, i.e. σ² = 0.1, **σ ≈ 0.32 cells**. At d = 1 cell the weight is e⁻⁵ ≈ 0.7%; at d = 2 it's e⁻²⁰ ≈ 2×10⁻⁹. Effectively this is **nearest-cell binning with slight smoothing**, despite the 4-cell-wide loop. The Kaiser–Bessel code (`bessi0`) exists in `recon.py` but is commented out.

In common delta-k units (Steve cell = 0.417 Δk, Faraz cell = 0.333 Δk): Steve's kernel σ ≈ 0.13 Δk; Faraz's KB half-width = 0.83 Δk. **Faraz's kernel is roughly 5–6× wider** in effective support.

### 5. Neighborhood / support

- **Steve:** `range(max(0,cx-2), min(MS,cx+2))` → cells `cx−2 … cx+1` per axis. **Asymmetric box** (excludes `cx+2`). Harmless given the kernel is ≈0 beyond 1 cell, but it is a quirk.
- **Faraz:** `knnsearch(..., 'K', kernelsize^2)` → **25 nearest grid points in 3D** (k², not k³ = 125). A 25-point ball has radius ≈ 1.8 cells, so the nominal 2.5-cell KB support is truncated. Second quirk: the kernel tail beyond ~1.8 cells never contributes.

### 6. Density compensation & normalization — the most consequential split

- **Faraz:** classic two-part scheme.
  1. **Iterative DCF** (`iterative_dcf_fa_20190910.m`, adapted from Zwart/Pipe MRM 2011): 5 iterations of grid→degrid with `sqrt(KB)` as the DCF kernel; `wi ← wi / (C·wi)` per sample. Samples are multiplied by `wi` **before** gridding.
  2. **Per-cell normalization:** gridded k-space divided by `Auxiliary = Σ_samples wi·kernel` per cell. Net effect: each grid cell holds a **wi-weighted average** of nearby samples.
- **Steve:** single mechanism. Grid accumulates `Σ kernel·data` and `Σ kernel` (`knorm`); `cudarenorm` divides: each cell holds an **unweighted kernel-average** of nearby samples. Cells with `knorm < 1e-5` are zeroed (labeled `TEMPORARY KLUGE`). There is a disabled Hermitian-symmetrization path (`oldway = 0`) in `cudarenorm`.

**Interpretation:** both are "normalized gridding" (cell = average of contributing samples), which inherently compensates sampling density *at cell resolution*. Because Steve's kernel barely reaches past one cell, his scheme reduces to bin-averaging — robust, no DCF needed, but the k-space weighting is effectively boxcar-per-cell. Faraz's wider KB kernel *would* smear density errors across cells, which is exactly why he needs the explicit iterative DCF inside the average. For a fully-sampled static phantom the two should converge to similar k-space estimates; differences show up where sampling density varies fast (spiral center, edge of k-space) and in the noise-vs-resolution tradeoff.

### 7. De-apodization — both omit it

Gridding convolves k-space with the kernel, multiplying the image by the kernel's Fourier transform (rolloff). Standard recon divides the image by this FT.

- **Faraz:** the comment "Partial de-apodization" on `KSpace_out./Auxiliary` is misleading — that divide is sample-weight normalization, not image-domain rolloff correction. No `FT(KB)` divide anywhere.
- **Steve:** nothing either; the Gaussian kernel's image-domain rolloff is uncorrected.

**Practical impact:** Steve's σ ≈ 0.32-cell Gaussian → image rolloff is a very wide Gaussian, nearly flat over the cropped FOV (per-cell normalization also re-flattens the kernel response locally). Faraz's 5-cell KB → FT rolloff is significant toward FOV edges, but the cell-normalization (divide by Auxiliary) largely cancels the kernel weighting in the same way. Neither image carries textbook apodization correction; both rely on normalization instead. For phantom comparison, expect mild shading differences near FOV edges.

### 8. FFT and crop

- **Steve:** `np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(k)))` then crop `[ISLL:ISUL]³` (central 100³ from 240³). Note: forward `fftn` from k-space to image (sign convention → spatial flip/conjugate relative to `ifftn`); `dyn_usimg_recon` additionally omits the inner `ifftshift` and flips axis 0 — but for the static path (`calcb`/`dyn_recon`) shifts are consistent.
- **Faraz:** per-dimension `fftshift(fft(fftshift(X,i),[],i),i)`, then crop central `NumK·zfill` cube. Also forward FFT.

Both use forward FFT k→image, so both carry the same parity/conjugation convention. Equivalent up to global flips. (Known pitfall: Steve's background is complex/zero-mean → looks negative; take `abs()` or compare real parts carefully.)

### 9. Coil combination

Strikingly similar in spirit — both are phased real-part combinations:

- **Steve** (`results.calcb`): reconstructs each channel from *all* data, normalizes by image-corner noise, `b = conj(img)`; masks low-signal voxels and fills their phase with a 10-term 3D quadratic polynomial fit (`curve_fit` on `angle(b)`); normalizes `b` by max magnitude across channels. Final: `Σ_ch real(img_ch · b_ch)`.
- **Faraz** (`combinecoils_fa.m`): normalizes each channel by corner-noise mean, `b = Σ_rep img / Σ_rep |img|` (coherence-weighted phase map), combination `real(b'·img)` per voxel, citing Bydder et al., MRM 47:539 (2002). Final image renormalized by corner noise.

Differences: Steve's `b` carries fitted smooth phase in low-SNR regions (less noise-bias in background phase); Faraz's `b` magnitude encodes inter-repetition phase coherence. **For a single-channel static phantom both collapse to: take the real part of the phase-referenced image.** Steve guarantees the zero-phase reference via raw-domain rephasing; Faraz gets it through `b`.

### 10. Averaging repeated interleaves

- **Faraz** (fancy_v3 path, `fa_spiral_dyn_recon.m:76-91`): explicitly averages repeated acquisitions of the same interleave/rotation (`rawdata2./weights`) **before** gridding, then grids one averaged set.
- **Steve:** grids every acquired interleave directly; the trajectory index wraps modulo `nuniquesmp` (`kidx = (idx+idxoff) % nuniquesmp`), and the `knorm` divide averages repeats implicitly at the cell level.

Mathematically near-equivalent for uniform repeat counts; differs if some interleaves are spike-zeroed (Steve's average weight then adapts per cell automatically).

---

## Quantitative kernel summary

| Quantity | Steve | Faraz |
|----------|-------|-------|
| Oversampled cell size | 0.417 Δk | 0.333 Δk |
| Kernel type | Gaussian, σ² = 0.1 cell² | KB, β ≈ 12.78, width 5 cells |
| Kernel weight at 1 cell | e⁻⁵ ≈ 0.007 | ≈ 0.5 (KB at u = 0.2·width/2... order 0.5) |
| Effective support radius | < 1 cell | ≈ 1.8 cells (knnsearch-truncated from 2.5) |
| Explicit DCF | none | Pipe–Menon, 5 iterations |
| Per-cell normalization | yes (`knorm`) | yes (`Auxiliary`) |
| De-apodization | no | no |

---

## What is genuinely different vs cosmetically different

**Genuinely different (expect measurable image differences):**
1. Kernel: near-nearest-neighbor Gaussian binning vs proper truncated Kaiser–Bessel interpolation. Faraz should show smoother k-space interpolation / fewer gridding artifacts at equal sampling; Steve trades that for simplicity and inherent density handling.
2. Density compensation: implicit bin-average only vs iterative DCF + weighted average. Largest impact where spiral sampling density changes rapidly.
3. Readout apodization: Steve Gaussian-filters the FID (gplb = 300); Faraz does not. Steve's PSF is slightly broadened by design.
4. killpts: Steve discards the first 2 samples of each interleave; Faraz keeps them.

**Cosmetically different (same idea, different mechanics):**
5. Coil combination — both phased real-part combinations with noise normalization; details of the `b` map differ.
6. Spike rejection — zero-and-skip vs interpolate; both threshold against cross-interleave statistics.
7. Repeat averaging — pre-gridding (Faraz) vs in-normalization (Steve).
8. FFT conventions — both forward FFT k→image with shift sandwiches; equivalent up to flips.
9. Oversampling 2.4 vs 3 with center crop.

**Shared properties (identical on both sides, relevant when comparing against anything else):**
10. De-apodization — neither side does it (detailed above).
11. Orientation/layout — Faraz permutes/flips to a display convention in his driver scripts (e.g., permute [X Z Y], L–R flip); Steve returns grid order (with an axis-0 flip in the dynamic path only). Side-by-side comparisons must match orientation first or images look mirrored/transposed.
12. B0 / off-resonance correction — **neither** pipeline corrects for off-resonance during the long spiral readout. Benign for a phantom on-resonance; would cause blurring/warping in vivo. Both stacks could add time-segmented or CG-based correction later; not active in either.

**Quirks/bugs noted in passing:**
- Steve: asymmetric gridding box (`cx−2…cx+1`); `knorm < 1e-5 → 0` kluge; disabled Hermitian symmetrization; `dyn_usimg_recon` missing inner `ifftshift` (dynamic path only, out of scope here).
- Faraz: 3D neighbor count is k² = 25 (not k³), silently truncating KB support; "partial de-apodization" comment mislabels weight normalization; `iterative_dcf` header admits "3D was modified and needs to be double-checked".

---

## Scoring

Criteria: theoretical correctness, SNR behavior, artifact risk, robustness, compute cost. Scores are per-approach out of 10; "Impact" rates how much the difference matters for the static-phantom comparison.

### Core four (genuinely different)

#### 1. Kernel — KB (Faraz) vs narrow Gaussian (Steve)

**Faraz: 8/10 · Steve: 5/10 · Winner: Faraz · Impact: HIGH**

| | Faraz KB | Steve Gaussian |
|---|---|---|
| Pros | Near-optimal interpolation (Beatty-optimized β); minimal aliasing energy folded into FOV; smooth k-space estimate between samples | Dead simple; no kernel table; no de-apodization debt (rolloff ≈ flat); inherently density-safe; fast on GPU |
| Cons | knnsearch truncation to 25 neighbors cuts the KB tail, so the β optimality is partly wasted; needs DCF to work; uncorrected KB rolloff shades FOV edges | σ ≈ 0.32 cells means bin-average gridding; noisier k-space estimate between samples; sample position quantized to cell → ~half-cell blur; resolution loss baked in |

Context: hyperpolarized gas is low-SNR with smooth objects, so Steve's effective box-binning costs high-frequency fidelity that is mostly absent anyway. But for sharp phantom edges — exactly the comparison case — KB wins visibly: edge sharpness and ringing will differ. For the CS pipeline downstream, a KB/NUFFT-style operator is the right base; adopt Faraz's geometry but fix the k³ neighbor truncation.

**Theory note — kernels exist only as a computational compromise.** The exact interpolation kernel for band-limited resampling is the sinc, with infinite support. Every practical gridding kernel is a truncation of it, and truncation buys speed at the cost of aliasing sidelobes. Kaiser–Bessel is the provably near-optimal truncation — Beatty's result is precisely that KB minimizes aliasing energy inside the FOV for a given support width. The quality ordering is therefore: sinc (exact, infinite support) → KB width 5 (good) → KB width 2 (worse) → Steve's σ ≈ 0.32-cell Gaussian (close to nearest-neighbor binning). Despite the GPU and the 4³ loop, Steve sits at the *most* truncated end of this chain, not the exact end: his per-sample loop visits 64 cells, but the Gaussian weight at 1 cell is already e⁻⁵ ≈ 0.7%, so only the nearest ~1–8 cells contribute meaningfully. His *effective* support is smaller than Faraz's truncated KB.

The genuinely kernel-free recon — "every sample contributes to every voxel" — is the **direct adjoint DFT**: `img(r) = Σᵢ wᵢ·dᵢ·exp(2πi kᵢ·r)`, evaluated in the image domain with no intermediate grid. For ~10⁶ samples × 100³ voxels this is ~10¹² operations — minutes on a modern GPU. Neither implementation does this. Two caveats apply even there: (a) it still needs density compensation, because the adjoint is not the inverse — the uncompensated adjoint gives a k-center-weighted blur; (b) the actually-exact reconstruction is the iterative inverse (CG-NUFFT / compressed sensing), where the kernel survives only inside the fast forward operator as a speed trick — or disappears entirely if an exact GPU DFT per iteration is affordable. **This is the natural baseline for the CS pipeline: exact DFT forward/adjoint on GPU, no kernel, no DCF, solver handles density and inversion.** It removes both implementations' approximations from the comparison at once.

One coupling to respect: Steve cannot simply widen his Gaussian to climb this quality ladder. His `knorm` normalization makes each cell a *weighted average* of nearby samples (kernel regression in k-space), not a convolution sum. Widening the kernel under that scheme bleeds neighboring samples into each cell's average → low-pass-filtered k-space → blurred image. Narrow-kernel + normalize is one coherent design; wide-KB + DCF + sum is the other. The kernel width and the density-handling scheme (finding 2) must be switched together.

**Why Steve's recon is the computationally heavy one — and it is not the kernel.** Per sample-visit the two are comparable (Steve: 64 cells × exp() × 3 atomic adds = 192 atomic memory ops; Faraz: 25 table-lookup multiply-adds, ≈ 2.5–7× lighter). The order-of-magnitude difference comes from *how much data gets gridded how many times*:

1. **No pre-averaging.** Steve grids the entire acquired dataset — all `ntotalilvs` interleaves, i.e. every repetition of the trajectory over the whole scan — every time (`results.py` reshapes the full `[npts, nch, ntotalilvs]` block per pass). Faraz collapses all repetitions into one averaged k-space of fixed size `nsamples × nleaves × nreps` *before* gridding (cheap accumulation), then grids only one trajectory's worth. For a long dynamic acquisition where the trajectory repeats R times, Steve's gridding workload per pass is R× Faraz's (R ~ 10–50).
2. **Soft binning inside the gridder.** Steve's bin weight `binwt = exp(-bindist²/2)` is evaluated per sample per bin pass, so each of the `nbins = 16` passes traverses the *full* dataset again: 16 × R × one-trajectory samples. Faraz pays for soft binning once, in the cheap accumulate step (weights `e^-d` applied while summing into `rawdata2`), then runs 16 gridding passes each on the small averaged set.
3. **No precomputation.** Faraz's neighbor indices (`Ind`), distances (`Dist`), kernel values, and DCF (`wi`) are computed once and cached to a `.mat` (`grid_lookup_20230113`) — across sessions this cost amortizes to zero. Steve recomputes nearest-cell indices and Gaussian weights for every sample on every pass (each individually cheap, but multiplied by points 1–2).
4. **GPU frictions.** Atomic adds serialize where many samples hit the same cells — which for a spiral is exactly the oversampled k-center, the worst case; plus a `cudarezero` + device-to-host copy of the 240³ complex grid per bin per channel.

Net: total cell-update count is roughly `R × nbins × 2.5` times Faraz's — easily 10²–10³× — which is why it needs a GPU *and can still be slower* than MATLAB looping over a precomputed, pre-averaged, cached problem. The irony: this cost does not buy interpolation quality (the kernel is effectively single-cell). It buys flexibility — soft bins, per-sample exclusion, no precompute step, identical code path for every recon variant. Restructuring Steve's pipeline to pre-average per bin would remove most of the cost without touching image quality for the static case.

#### 2. Density compensation — iterative DCF (Faraz) vs knorm-only (Steve)

**Faraz: 7/10 · Steve: 6/10 · Winner: Faraz by a slim margin · Impact: MEDIUM**

| | Faraz Pipe–Menon + normalization | Steve knorm normalization |
|---|---|---|
| Pros | Correct approach for wide kernels; handles spiral-center oversampling properly; 5 iterations converge well for smooth spirals; DCF cached in lookup .mat | Zero tuning; cannot diverge; adapts automatically when samples are zeroed (spikes, bins); exact density compensation at cell resolution |
| Cons | His own header: "3D was modified and needs to be double-checked"; sqrt(KB) DCF kernel choice is ad hoc; double-counting risk — applies wi *and* divides by Σ(wi·kernel), so the DCF partly cancels itself (net effect: wi acts as reliability weights inside a cell-average, not classic DCF) | Boxcar density weighting per cell → noise amplification at sparsely sampled k-edge; `knorm < 1e-5 → 0` kluge leaves k-space holes instead of interpolating |

Subtle point: because both divide by accumulated weights, both are "averaging gridders," and Faraz's DCF inside an average mostly redistributes within-cell weights. The practical difference is smaller than "DCF vs no DCF" sounds. Faraz's combo is also non-standard — convention is either DCF + plain sum, or no DCF + normalize, not both. For CS this finding evaporates: an iterative solver with an explicit NUFFT forward model replaces DCF entirely.

#### 3. Readout apodization — gplb=300 Gaussian (Steve) vs none (Faraz)

**Steve: 6/10 · Faraz: 7/10 · Winner: situational · Impact: HIGH for comparison validity**

| | Steve filter | Faraz unfiltered |
|---|---|---|
| Pros | Suppresses late-FID noise (low-SNR tail of each spiral arm); kills T2*-decayed garbage; better apparent SNR | True acquired resolution; no hidden PSF broadening; what is measured is what is gridded |
| Cons | PSF broadened and undocumented — protocol-stated resolution ≠ delivered; gplb = 300 is a hardcoded magic number, not derived from T2* or trajectory; conflates filtering with reconstruction | Late-readout noise enters the image; for short-T2* gas the k-edge is mostly noise anyway |

Steve's filter is defensible physics (the signal genuinely decays along the readout) but should be explicit, parameterized, and reported — not buried in the CUDA kernel. **For the phantom comparison: set Steve's smoothing to 0 or apply the identical filter to Faraz's data first; otherwise the comparison measures filters, not recons.**

#### 4. killpts — drop 2 samples/interleave (Steve) vs keep all (Faraz)

**Steve: 7/10 · Faraz: 5/10 · Winner: Steve · Impact: LOW–MEDIUM**

| | Steve drop | Faraz keep |
|---|---|---|
| Pros | First samples sit in ADC settling / gradient ramp distortion; trajectory calibration is least accurate there; cheap insurance | No data loss; k-center samples have the highest SNR; trusts the trajectory calibration |
| Cons | Discards the highest-signal k-center points; hardcoded 2, not derived from hardware timing; only partially compensated by rephasing | If the first samples are corrupted, the error lands at k-center → global scale/shading error, the worst possible place |

Small in count, big in leverage: 2 points × all interleaves at k ≈ 0 dominate the DC term. If Faraz's calibration accurately measures the actual ramp trajectory, keeping them is strictly better; if not, Steve's caution wins. Testable cheaply: recon Faraz's pipeline both ways on the phantom and compare DC/shading.

### Cosmetic five (items 5–9)

#### 5. Coil combination — polynomial-phase b map (Steve) vs Bydder Σimg/Σ|img| (Faraz)

**Steve: 8/10 · Faraz: 7/10 · Winner: Steve, marginal · Impact: LOW (zero for single channel)**

| | Steve | Faraz |
|---|---|---|
| Pros | b from full-data (highest SNR) recon; polynomial phase fill prevents noise-driven phase in low-signal voxels from corrupting the combine; channel sensitivity magnitude retained | Standard, citable method (Bydder MRM 2002); self-contained — no separate calibration recon; coherence weighting across repetitions is a free reliability metric |
| Cons | Complex pipeline (mask, 10-term 3D quadratic curve_fit per channel); fit can misbehave with odd coil geometry; bespoke, harder to validate | b magnitude from Σ/Σ\|·\| degrades when per-repetition phase drifts; noisy phase in low-signal regions enters the combine unfiltered |

Both are phased real-part combinations — same family. For a single-channel phantom both collapse to "take the real part," so this difference is invisible in the planned comparison.

#### 6. Spike rejection — zero-and-skip (Steve) vs neighbor-interpolate (Faraz)

**Steve: 7/10 · Faraz: 6/10 · Winner: Steve, marginal · Impact: LOW**

| | Steve | Faraz |
|---|---|---|
| Pros | A zeroed sample contributes nothing and knorm renormalizes around it — statistically clean (missing data stays missing); threshold (10× cross-interleave mean) is simple and conservative | Moving-average threshold (1 + 3σ) is more sensitive — catches smaller spikes; interpolation preserves sample count |
| Cons | Conservative threshold misses moderate spikes; zeroing is only clean *because* of the knorm averaging — would bias a plain-sum gridder | Interpolated values are fabricated data; neighbor mean along readout assumes smooth k-space signal, which fails at spiral center; more aggressive threshold risks clipping genuine signal transients |

Philosophically: discard vs impute. With an averaging gridder, discard is the safer default. Difference only matters on noisy acquisitions; negligible for a clean phantom scan.

#### 7. Repeat averaging — pre-gridding (Faraz) vs in-normalization (Steve)

**Faraz: 7/10 · Steve: 7/10 · Tie · Impact: NEGLIGIBLE**

| | Faraz pre-average | Steve implicit |
|---|---|---|
| Pros | Smaller gridding workload (one averaged set); explicit weights array makes the averaging auditable | No bookkeeping; weight adapts per cell automatically when samples are dropped (spikes/exclusions) |
| Cons | NaN risk where weights = 0 (handled later by NaN→0); fixed equal weighting across repeats even if SNR varies | Averaging hidden inside knorm — harder to inspect intermediate |

Mathematically near-equivalent for uniform repeat counts. Only diverges when samples are unevenly dropped, where Steve's per-cell adaptation is slightly more graceful.

#### 8. FFT conventions — both forward FFT k→image

**Both: 6/10 · Tie · Impact: NEGLIGIBLE (if handled)**

Both use forward FFT instead of inverse for k→image, carrying the same conjugation/parity convention (global spatial flip vs ifft-based recon). Since both sides share the quirk, it cancels in comparison — but anyone comparing either against a third implementation (e.g., a CS solver using standard ifft) must account for the flip and conjugation. Steve additionally flips axis 0 in the dynamic path and omits an inner ifftshift there (out of scope for static).

#### 9. Oversampling — 2.4× (Steve, MS/IS) vs 3× (Faraz, os)

**Faraz: 7/10 · Steve: 6/10 · Winner: Faraz, marginal · Impact: LOW**

| | Faraz os = 3 | Steve MS/IS = 2.4 |
|---|---|---|
| Pros | More guard band — gridding aliasing replicas pushed further outside the crop; os explicit and tunable | Smaller grid (240³ fixed) → less memory/compute; with a sub-cell kernel, aliasing from kernel sidelobes is minimal anyway |
| Cons | 1.95× more grid voxels than Steve's (720³-equivalent spacing on same FOV basis) — memory and FFT cost | Oversampling factor implicit in MS/IS ratio — changing IS silently changes gridding behavior; 2.4 is below the common ≥2.5 comfort zone for wide kernels (fine for his narrow kernel, marginal if KB were re-enabled) |

Literature (Beatty) says os ≥ 2 suffices with a properly sized KB kernel; both are adequate for their respective kernels. Note the coupling: if Steve's commented-out KB were re-enabled, the 2.4 factor plus the asymmetric 4-cell box would need rechecking together.

#### 10. (Shared) De-apodization — both omit it

**Both: 4/10 · Impact on cross-comparison: LOW · Impact vs ground truth: MEDIUM (Faraz side)**

Neither divides the image by the kernel's Fourier transform. Steve gets away with it — a σ ≈ 0.32-cell Gaussian has a nearly flat rolloff across the cropped FOV. Faraz's 5-cell KB has a real rolloff toward FOV edges that his Auxiliary normalization only partially flattens. Consequence: Faraz's phantom may show edge-of-FOV shading that is recon-induced, not physical. When comparing both against a known phantom geometry (rather than against each other), add KB de-apodization to Faraz's image or restrict ROIs to the central FOV.

### Scoreboard

| # | Finding | Steve | Faraz | Winner | Impact on phantom comparison |
|---|---------|-------|-------|--------|------------------------------|
| 1 | Kernel | 5 | 8 | Faraz | **High** — edge sharpness, ringing |
| 2 | Density compensation | 6 | 7 | Faraz (slim) | Medium — k-edge noise behavior |
| 3 | Readout apodization | 6 | 7 | situational | **High** — confounder, equalize first |
| 4 | killpts | 7 | 5 | Steve | Low–medium — DC/shading |
| 5 | Coil combination | 8 | 7 | Steve (marginal) | Low; zero for single channel |
| 6 | Spike rejection | 7 | 6 | Steve (marginal) | Low |
| 7 | Repeat averaging | 7 | 7 | tie | Negligible |
| 8 | FFT convention | 6 | 6 | tie | Negligible if flips handled |
| 9 | Oversampling | 6 | 7 | Faraz (marginal) | Low |
| 10 | De-apodization (shared omission) | 4 | 4 | — | Low between them; medium vs ground truth |

**Comparison protocol implications:** neutralize 3 (set gplb = 0 or filter both) and 4 (test killpts both ways) before attributing image differences to 1 and 2, which are the genuine algorithmic comparison. Items 5–9 are implementation style, not algorithm. Item 10 only matters when judging against known phantom geometry.

**CS pipeline implications:** take Faraz's KB geometry but with full k³ kernel support; drop DCF (the iterative solver's forward model handles density); make any readout filtering an explicit, reported parameter; settle killpts empirically on phantom data.

---

## Why each beats the textbook — estimator objectives

*(Added 2026-06-11 after the ACR phantom measurements: Steve's recon measured
SNR 28.7 vs the unbiased CG inverse's 19.6; Faraz's measured the flattest
interiors, lowfreq-CV 0.093 vs CG's 0.110. Both "beat" the textbook-optimal
least-squares reconstruction. Neither result is paradoxical.)*

The resolution of the paradox is that "textbook optimal" means optimal **for
the criterion the textbook chose** — minimum variance *among unbiased
estimators*. Total error is MSE = bias² + variance, and a biased estimator
beats the unbiased one whenever its bias is cheap for the object at hand. Both
implementations quietly optimize different criteria, and each wins exactly the
metric its design implicitly targets.

### Steve's SNR: minimum variance, bias accepted

1. **The object's k-space is mostly empty.** Phantom and lung spectra
   concentrate at low |k|; the outer half of k-space is essentially pure
   noise. Steve's cell-averaging and gplb filter suppress fluctuations there
   at near-zero signal cost. The unbiased inverse must faithfully reconstruct
   that noise — unbiasedness *requires* it.
2. **Kernel regression averages within cells**: variance ÷ effective samples
   per cell, with bias (k-space assumed constant over 0.4 Δk) negligible for
   any object smaller than the FOV.
3. **The SNR metric rewards low-pass** (object mean / background std improves
   under smoothing until visible blur). Steve sits just below that threshold.

Measured decomposition: his gplb filter applied to our CG data buys +2 SNR
(19.5 → 21.5); the remaining ~7 points are the cell-averaging bias–variance
trade. Confirmed negatively by the λ sweep (`cg_tune.py`): Tikhonov over four
decades moves CG's SNR only 19.6 → 20.9 — λI penalizes amplitude, not
roughness, so it cannot buy what Steve's smoothing buys.

### Faraz's homogeneity: flatness is literally his objective function

The Pipe–Menon fixed point is `w ⊛ C = 1` across sampled k-space — iterate
until the **net spectral weighting of the whole acquisition+recon chain is
flat**. A flat transfer function is precisely what the eye reads as a uniform
image; he runs five iterations of an algorithm whose convergence criterion
*is* flatness. Secondary contributors: smooth KB interpolation (no
cell-quantization ripple — Steve's nearest-cell binning staircases sample
positions, producing the faint low-frequency swirl), and 3.125 mm voxels
averaging more per voxel.

The unbiased LS solution optimizes per-sample data fidelity instead; its
spectral response is not constrained flat, and finite-iteration CG makes it
worse — CG converges eigencomponents unevenly (extremes first), leaving
clustered low-frequency modes partially equalized → the measured lowfreq-CV
~0.12 that persists even at 30 iterations.

### The unified picture

| Recon | Implicit objective | Wins | Pays |
|---|---|---|---|
| Steve | minimum variance (smoothing bias accepted) | SNR | resolution tax, k-edge fidelity |
| Faraz | flat spectral transfer (DCF fixed point) | homogeneity | DCF heuristics; grid-overflow fragility (the `resizing` zoom bug) |
| Textbook CG | unbiased data fidelity | resolution, faithfulness, geometry | keeps all the noise; finite-iteration low-freq shading |

Nobody beat the textbook — each beat it on the axis the textbook never
optimized, by paying on an axis it did. Both veterans encoded application
knowledge as estimator bias without writing the objective down.

**Consequence for the CS step:** compressed sensing is the framework where the
objective *is* written down — `min ‖Aρ − s‖² + λ·R(ρ)` — so the bias is chosen
deliberately (sparsity/smoothness matched to lungs) rather than inherited from
the gridder's plumbing. Steve's SNR advantage and Faraz's flatness both become
tunable terms of one explicit functional; the expectation is to match or beat
both on their home metrics simultaneously, at equal or better resolution.
