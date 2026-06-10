# Non-Cartesian Reconstruction — From FID to Image

**Date:** 2026-06-10
**Audience:** knows MRI basics, the FID, and the 2D FFT/DFT. Everything else built from there.
**Goal:** understand *why* every step in Steve's (`recon.py`/`results.py`) and Faraz's (`gridrecon_fa_20230113.m`) code exists.
**Companions:** `Recon_Overview_Steve.md`, `Recon_Overview_Faraz.md`, `Recon_Comparison_StaticGas.md`

---

## 1. The bridge: an FID sample is a k-space sample

Start from what you know. After excitation, every spin at position **r** precesses and contributes to the FID. With no gradients, all spins precess at the same frequency and the FID is just a decaying oscillation. Turn on a gradient **G**(t) and the precession frequency becomes position-dependent: spins at **r** accumulate extra phase

```
φ(r, t) = -2π k(t)·r        where        k(t) = (γ/2π) ∫₀ᵗ G(τ) dτ
```

The receiver coil sums all spins, so the demodulated signal at time t is

```
s(t) = ∫ ρ(r) · e^(-i2π k(t)·r) dr
```

Read that carefully: **the signal at time t is the Fourier transform of the object, evaluated at the single frequency-space point k(t).** The FID is not "a signal that gets Fourier transformed" — each individual ADC sample of the FID *is already one point of the object's Fourier transform*. The gradient waveform steers a pen through Fourier space (k-space), and the ADC reads off the transform value wherever the pen happens to be.

This is the entire conceptual content of MRI spatial encoding. Everything downstream — Cartesian, spiral, radial — differs only in *where the pen goes*.

Two immediate consequences:

1. **The trajectory is part of the measurement.** You cannot reconstruct without knowing k(t) — and k(t) is whatever the gradients *actually did*, not what you programmed. (This is why both codebases load *measured* trajectories — section 9.)
2. **The data is complex.** s(t) comes from quadrature demodulation; both magnitude and phase of each k-space sample carry information.

## 2. Why Cartesian recon is just an FFT (recap, with the parts that will break)

In a spin-warp (Cartesian) sequence, phase-encode steps and a constant readout gradient place the samples on a *uniform rectangular grid* in k-space: spacing Δk, extent ±k_max. Then the imaging equations are exactly the DFT's assumptions, and two famous relations fall out:

```
FOV = 1/Δk                 (sample spacing sets unaliased field of view)
δx ≈ 1/(2·k_max)           (k-space extent sets resolution)
```

The 2D FFT works because of three properties that are *easy to take for granted*:

- **Uniform spacing** — the FFT's input slots correspond to equally spaced frequencies. Period.
- **Uniform density** — every region of k-space is sampled exactly once; no part of the transform is over-counted.
- **Grid alignment** — samples sit exactly where the DFT expects them, no interpolation needed.

Non-Cartesian acquisition violates all three at once. The whole subject of non-Cartesian reconstruction is repairing these three violations efficiently.

## 3. Why bother: what spirals buy (especially for hyperpolarized Xe)

A spiral interleave starts at the k-space center and winds outward. For ordinary ¹H imaging spirals buy speed (long efficient readouts, few excitations). For **hyperpolarized ¹²⁹Xe** the case is much stronger, and it is worth understanding because it explains design choices in both codebases:

1. **Non-renewable magnetization.** Hyperpolarized magnetization does not recover by T1 — every RF pulse spends an irreplaceable resource. You want maximal k-space coverage per excitation. A spiral collects a 2D/3D swath per shot; Cartesian collects one line.
2. **Center-out = signal-first.** The k-space center (overall brightness, contrast) is sampled at the very start of the FID, before T2* decay eats it. Edges (fine detail) get the decayed tail. Graceful degradation.
3. **Center oversampling = self-navigation.** Every interleave re-measures k ≈ 0. The first samples of each FID track total signal ∝ gas in lungs → free respiratory signal. Both pipelines exploit this (Steve: `volmeasvol[SIGNAL]` from first 8 points, `raw.py:429`; Faraz: `findpeaks` on first sample, `spiral_human:184`). Out of scope for static recon but it explains why the data layout is what it is.
4. **Motion robustness & undersampling tolerance.** Spiral aliasing smears into diffuse swirls instead of coherent ghosts — and that incoherence is exactly what compressed sensing wants (section 12).

Cost: samples land *nowhere near* a uniform grid, density varies enormously (center heavily oversampled, edge sparse), and the long readout makes the image sensitive to off-resonance and trajectory errors. Hence this document.

## 4. The reconstruction problem, stated honestly

Discretize the object into N voxels ρ ∈ ℂᴺ. The acquisition is a linear operator:

```
s = A ρ        where        A[i, j] = e^(-i2π k_i · r_j)
```

A is an M×N matrix of complex exponentials — a **Fourier transform evaluated at arbitrary frequencies** (a "non-uniform DFT", NUDFT). Reconstruction means inverting this. Three escalating answers:

**(a) The adjoint ("conjugate phase" recon).** Apply Aᴴ:

```
ρ̂_adj(r_j) = Σ_i  s_i · e^(+i2π k_i · r_j)
```

Every sample contributes to every voxel. Cost O(M·N) — for 10⁶ samples × 100³ voxels ≈ 10¹² operations; minutes on a GPU, weeks in 1990. This is the "use all the data, no kernel" recon. **But Aᴴ ≠ A⁻¹.** For Cartesian sampling AᴴA = identity (orthogonality of DFT exponentials — this is why the FFT *is* the inverse there). For non-uniform sampling AᴴA is not the identity: the result is the true image *convolved with the point spread function of the sampling pattern*, and over-sampled regions (spiral center) are counted many times → a heavy low-frequency bias, i.e. a blurry, shaded image. Fixing this weighting is the **density compensation** problem (section 5).

**(b) Gridding.** Approximate the same adjoint in O(M·w³ + N log N) by interpolating samples onto a Cartesian grid and using the FFT (sections 6–8). This is what both Steve and Faraz do. Gridding is not a different *idea* from (a) — it is a fast approximation of the density-compensated adjoint.

**(c) The actual inverse.** Solve `min_ρ ‖Aρ − s‖² (+ regularization)` iteratively (CG, FISTA…). Each iteration applies A and Aᴴ (via gridding/NUFFT or exact DFT). This is where compressed sensing lives, and it makes the DCF question disappear because the solver itself accounts for AᴴA (section 12). Neither current pipeline does this — it is the planned destination.

## 5. Density compensation: the quadrature-weights problem

Look again at the adjoint sum and compare it with the continuous inverse transform it tries to approximate:

```
ρ(r) = ∫ ŝ(k) e^(+i2π k·r) dk     ≈     Σ_i  w_i · s_i · e^(+i2π k_i·r)
```

Approximating an integral by a sum over scattered points needs **quadrature weights**: each sample must be weighted by the amount of k-space it "represents". Samples crowded together near the spiral center each represent a tiny patch (small w_i); lonely samples at the k-space edge represent a big one (large w_i). Without w_i, the center is massively over-counted → image dominated by a blurry low-frequency blob.

This is the **density compensation function (DCF)**. Ways to get it:

- **Analytic / geometric.** Voronoi tessellation: w_i = volume of sample i's Voronoi cell — the literal "area it represents". Exact in principle, expensive in 3D, fragile at duplicated points (spiral center!). For analytically designed spirals there are closed-form approximations (Meyer / Hoge). Faraz's generic loader exposes these as options (`recon_20210622.m` `dcf` option: 0 = gridding-based, 1 = Meyer, 2/3 = Voronoi, 4 = none) — historical record of him trying them all.
- **Iterative (Pipe & Menon 1999).** Insight: the correct w makes the gridded sampling function flat. Fixed-point iteration:

  ```
  w ← w / (w ⊛ C)|at sample locations
  ```

  Grid the current weights with kernel C, read the gridded density back at each sample location (degridding), divide. Where samples crowd, (w ⊛ C) > 1 → weight shrinks. Converges in a handful of iterations. **This is exactly `iterative_dcf_fa_20190910.m`**: grid → degrid → `wi = wi ./ w1i`, 5 iterations, with `sqrt(KB)` as the convolving kernel C. (His header admits the 3D version "needs to be double-checked"; the sqrt choice is folklore, not theory.)
- **Implicit, by normalization.** Don't weight samples at all; instead make each *grid cell* the **average** (not sum) of the samples that land in it. Averaging is self-normalizing: a cell hit by 100 center samples and a cell hit by 1 edge sample both end up with one representative value. **This is exactly Steve's `knorm`** (`recon.py`: accumulate `Σ wt·s` and `Σ wt`; `cudarenorm`: divide). It is density compensation *at the resolution of one grid cell* — perfectly adequate when the kernel is narrower than a cell (Steve's is), increasingly wrong as the kernel widens (which is why Faraz, with a wide kernel, needs the explicit DCF).

The deep point: **DCF is only needed because we stopped at the adjoint.** The true inverse (option c) has no DCF — the solver figures out the weighting implicitly by minimizing the data misfit. DCF is a patch that makes a one-shot adjoint look like an inverse.

## 6. Gridding: the fast adjoint, step by step

We want the adjoint but with an FFT. The FFT demands samples *on* a uniform grid. So: move the data onto a grid first. Nearest-neighbor assignment would be crude (you'd shift each sample by up to half a cell — phase errors, blur, aliasing). The principled version is **convolution interpolation**:

### 6.1 The idea as a convolution

The measured data is the continuous transform ρ̂(k) multiplied by the sampling pattern (a sum of delta functions at the k_i, with DCF weights folded in). Convolve this with a smooth, compact kernel C(k):

```
(measured, weighted samples) ⊛ C   →   smooth function defined everywhere   →   evaluate on grid
```

Concretely, for each sample i: visit grid cells g near k_i, add `w_i · s_i · C(k_i − g)` into cell g. Then FFT.

The convolution theorem tells you the price: convolving by C in k-space **multiplies the image by c(r) = FT⁻¹{C}**. A smooth bump in k-space ↔ a broad rolloff (apodization) across the image. So textbook gridding has four moves:

```
1. weight       s_i → w_i·s_i                  (density compensation)
2. spread       onto grid via kernel C          (convolution interpolation)
3. FFT          grid → image
4. de-apodize   image ÷ c(r)                    (undo the kernel's rolloff)
```

### 6.2 Kernel choice: why Kaiser–Bessel, and what "optimal" means

What would the *perfect* C be? One whose image-domain footprint c(r) is a perfect box covering exactly the FOV: then the convolution does ideal band-limited interpolation and de-apodization is trivial. The k-space function with a box transform is the **sinc — with infinite support.** Every sample would touch every grid cell: we're back to O(M·N), the thing gridding exists to avoid.

So truncate. Truncation in k-space = ringing/sidelobes in c(r) = energy from *outside* the FOV folding back in = **aliasing artifact specific to the kernel**, on top of any undersampling aliasing. The kernel design problem is: *for a support of w cells, minimize the aliasing energy that lands inside the FOV*. The solution family is the prolate spheroidal wave functions; the **Kaiser–Bessel** function is their closed-form near-twin:

```
C(κ) = I₀( β·√(1 − (2κ/w)²) ) ,   |κ| ≤ w/2
```

with the shape parameter β tuned to the support w and the grid oversampling factor α. Beatty et al. 2005 give the standard recipe:

```
β = π·√( w²/α² · (α − ½)² − 0.8 )
```

**This exact formula is in Faraz's code** (`grid_lookup_20230113.m:28`), with w = 5, α = 3 → β ≈ 12.78. `createKBkernel.m` tabulates C at 10000 points; samples look up kernel values by their distance to each neighbor cell (`interp1` on `Dist`). The analytic de-apodization function for KB is also closed-form (a sinh-like expression) — *neither codebase applies it* (see 6.4).

A Gaussian kernel (Steve's `exp(−d²/0.2)`, `recon.py:129`) is the "cheap and cheerful" alternative: trivially evaluated on the fly (no table), smooth, but its transform is another Gaussian — never flat in the passband, never compactly supported — so for equal support it leaks more aliasing and shades more. Steve's, however, is so narrow (σ ≈ 0.32 cells) that it is better understood not as an interpolation kernel at all but as **bin-averaging**: each sample effectively votes into its nearest cell. With his knorm averaging, the whole scheme is the statistician's *kernel regression* (Nadaraya–Watson) estimate of k-space on the grid:

```
k̂(g) = Σ_i C(k_i − g)·s_i  /  Σ_i C(k_i − g)
```

This single formula **is** `cudarecon` + `cudarenorm`. (Faraz's `./Auxiliary` divide makes his pipeline a wi-weighted version of the same estimator — both codebases are normalized-average gridders, not textbook sum-gridders. See comparison doc §6 for why this means Faraz's DCF partly cancels itself.)

### 6.3 Oversampling the grid: why MS > IS, why os = 3

Aliasing from kernel truncation, and replicas from gridding itself, fold into the image with periodicity = grid FOV. Trick: grid onto a *finer* k-space grid (factor α) → image FOV is α× larger → the junk lands in the outer ring → **crop the center** and throw the ring away.

- Faraz: explicit `os = 3`, grid spacing `1/(fov·os)`, crops `NumK·zfill` after FFT (`gridrecon_fa_20230113.m:268`).
- Steve: implicit α = MS/IS = 240/100 = 2.4; trajectory rescaled by `MS/IS` (`raw.py:68`), crops `ISLL:ISUL` after FFT (`results.py:103`).

Beatty's analysis says α ≥ ~2 with a well-chosen KB kernel is ample. Larger α costs memory and FFT time cubically in 3D — the 240³ complex grid is most of Steve's GPU memory traffic.

`zfill` (Faraz) is a different, purely cosmetic knob: zero-padding to interpolate the *display* grid. Zero-filling adds no information — it sinc-interpolates the image.

### 6.4 De-apodization — the step both codebases skip

After the FFT the image is multiplied by c(r), the kernel's transform: bright center, dim edges. Textbook fix: divide by the known c(r). Neither pipeline does. Why they (mostly) get away with it:

- Steve: a σ ≈ 0.32-cell Gaussian in k-space → c(r) is a *very wide* Gaussian in image space → essentially flat over the cropped 100³ region.
- Faraz: w = 5 KB → real rolloff toward the crop edge, *partially* flattened by his per-cell normalization (the `./Auxiliary` divide removes the kernel's local *weighting*, though not its *smoothing*). Residual edge shading is recon-induced — relevant when validating against known phantom geometry (comparison doc §10).

### 6.5 The FFT itself: shifts, directions, and flips

Three conventions to keep straight whenever comparing recons:

1. **fftshift sandwiches.** The DFT indexes frequencies 0…N−1 with the "center" at index 0; images want k = 0 in the middle. `ifftshift → fft → fftshift` moves between conventions. Get one shift wrong → checkerboard phase (±1 alternation) in the other domain. Steve: `fftshift(fftn(ifftshift(k)))` (`results.py:99`). Faraz: per-dimension `fftshift(fft(fftshift(X,i),[],i),i)`.
2. **Forward vs inverse.** Going k → image "should" be the inverse FFT. Both codebases use the **forward** FFT. Harmless for one pipeline alone — it conjugates/flips the image (and both deal in real parts/magnitudes at the end) — but it matters the moment you compare against a third implementation that uses `ifft`: expect mirrored geometry and conjugated phase.
3. **Layout.** Reshape order (`order='F'` in Steve matching MATLAB column-major), axis permutes and flips for display (Faraz `fliplr(permute(img,[1 3 2 4]))`). Pure bookkeeping; first thing to rule out when two recons of the same data "look different".

## 7. Multi-coil combination: why it's not just a sum

Each receive channel sees the object through its own complex sensitivity b_ch(r) (magnitude = proximity, phase = geometry):

```
img_ch(r) = b_ch(r) · ρ(r) + noise_ch
```

Summing channels naively lets their differing phases cancel. The classic fixes:

- **Sum-of-squares (SoS):** `√Σ|img_ch|²`. No phase needed, magnitude-only, background noise becomes strictly positive (Rician floor).
- **Phased / matched-filter combine (Roemer; Bydder 2002):** estimate b_ch, compute `Σ Re(b*_ch · img_ch)`. SNR-optimal-ish, *keeps signed/complex information*, background noise stays zero-mean (negative pixels are normal!).

Both pipelines use the phased flavor, differing in how they estimate b:

- Steve `results.calcb`: reconstruct each channel from *all* data (highest SNR), b = conj of that; where signal is too weak to trust the phase, replace it with a smooth 3D quadratic polynomial fitted to the well-lit voxels (`curve_fit` over 10 terms, `results.py:145-149`). The final dyn/static images are `Σ Re(img_ch · b_ch)`.
- Faraz `combinecoils_fa`: normalize channels by corner noise, b = `Σ_rep img / Σ_rep |img|` (phase-coherence-weighted average over repetitions), combine `real(b'·img)`.

Practical consequence that bit us already: a phased-combined image has **zero-mean noise** — the background flickers negative. That is correct behavior, not a bug. Take `abs()` only for display. (The old archived comparison doc got this wrong for Faraz — it claimed magnitude/Rician.)

## 8. Hyperpolarized-specific physics in the pipeline

These show up directly as code:

1. **Signal varies across interleaves.** RF consumption + T1 decay + gas redistribution make later interleaves dimmer. For static recon this acts as a smooth weighting across k-space coverage (mild blur); both pipelines monitor it (Steve's exclusion of low-SNR fully-sampled blocks, `raw.py:397-412`; Faraz's manual `range` pick).
2. **T2* decay along each readout.** The FID decays *during* the spiral. k-space edge (late samples) is weighted down by physics no matter what you do. Steve leans into it: explicit Gaussian readout filter (`fwt = exp(−(t/gplb)²)`, `recon.py:98-99`) — accepting deliberate extra apodization for noise suppression. Faraz grids the raw decay-weighted data as-is. (Comparison confounder #3.)
3. **Noise spikes.** Hardware bursts stand out brutally in low-signal hyperpolarized data. Steve: zero them (the averaging gridder auto-renormalizes around holes). Faraz: interpolate from readout neighbors. Both threshold against cross-interleave statistics of the same readout index — the right reference population, since each readout index has its own decay-set signal level.
4. **Channel noise scaling.** Gas-phase signal levels vary scan to scan; meaningful channel combination needs channels on a common noise scale. Steve: fit a Gaussian to the histogram of Re(raw) per channel and divide by σ (`raw.py:211-215`). Faraz: divide by image-corner noise per channel inside the combine.

## 9. Trajectory calibration: why nobody trusts the programmed spiral

The reconstruction needs k_i = where the samples *actually* were. Real gradient chains low-pass the demanded waveform: amplifier bandwidth, eddy currents, group delays (μs-scale timing offsets between gradient axes and ADC). For Cartesian imaging these mostly shift the image; for spirals they *warp the trajectory* — and a wrong trajectory in gridding is a wrong **position** for every sample: blur, swirl artifacts, ghosting.

Hence both pipelines use **measured trajectories** (the `build_calibration_from_xyz` work — see Obsidian notes): play the spiral on a phantom, measure the actual k(t) (e.g. Duyn-style thin-slice phase method), store it.

- Faraz: `calibrations_3D_20220308.mat`, entries matched on (frequency ±5%, FOV ±5%, nsamples, nleaves, nreps, imgsize, orientation, nucleus) — `loadtrajectory3D.m`.
- Steve: `.npy` per protocol, × FOV → delta-k units; infers samples-per-interleave from the periodicity of |k(t)|² (`raw.py:41-47` — clever, fragile); drops the first `killpts = 2` samples where the measured ramp is least trustworthy.

Related: the first ADC samples sit on the gradient ramp-up where trajectory error is worst — but they are also the k ≈ 0 samples with the most signal leverage. Steve discards (insurance), Faraz keeps (trusts calibration). Comparison confounder #4.

## 10. One sample's journey through each codebase

The same physics, two dialects. Take ADC sample number i, value s_i, of channel ch.

**Steve (`cudarecon`, one GPU thread per sample):**
```
1. fwt   = exp(-((i mod npts)/gplb)²)        readout T2*-matched filter      [§8.2]
2. s     = raw[i]·fwt; skip if |s|≈0          spike holes skipped             [§8.3]
3. kidx  = (i + idxoff) mod nuniquesmp        trajectory repeats wrap         [§9]
4. (kx,ky,kz) = traj[kidx]                    already in grid units (·MS/IS + MS/2)
5. binwt: nbins=1 → 1.0                       soft binning inert for static
6. visit 4³ cells around nearest cell:
      wt = exp(-d²/0.2)                       narrow Gaussian kernel          [§6.2]
      atomic: k[cell]     += wt·s
      atomic: knorm[cell] += wt
   …later, cudarenorm:  k /= knorm            kernel-regression normalize     [§5, §6.2]
7. fftshift(fftn(ifftshift(k))) → crop IS³    oversample-and-crop             [§6.3, §6.5]
8. × b_ch, take real, Σ over channels         phased combine                  [§7]
```

**Faraz (`gridrecon_fa_20230113`, MATLAB loop over samples):**
```
0. (once, cached) Ind/Dist via knnsearch(K=25); KB table; wi via 5 Pipe–Menon iters   [§5, §6.2]
1. s averaged into rawdata2 with repeat-mates   pre-averaging                  [comparison §10]
2. M_i = s_i · wi_i                             explicit DCF                   [§5]
3. for the 25 neighbor cells:
      k[Ind(i,:)]   += M_i · KB(Dist(i,:))      wide KB spread                 [§6.2]
      Aux[Ind(i,:)] += wi_i · KB(Dist(i,:))
   …later:  k ./= Aux                           weighted-average normalize     [§5]
4. per-dim fftshift(fft(fftshift(·))) → crop    oversample-and-crop            [§6.3, §6.5]
5. combinecoils_fa: real(b'·img)                phased combine                 [§7]
6. fliplr(permute(...))                         display orientation            [§6.5]
```

Same skeleton: **(weight) → spread → normalize → FFT → crop → combine.** Every difference between the two is a choice *within* a step, never a different skeleton. That is the most useful single sentence for reading either codebase.

## 11. Concept → code lookup table

| Concept | Section | Steve | Faraz |
|---------|---------|-------|-------|
| Signal eq. / k-space sample | §1 | (implicit in data) | (implicit in data) |
| Measured trajectory | §9 | `raw.py traj.load` | `loadtrajectory3D.m` |
| Trajectory → grid units | §6.3 | `traj.rescale_to_MS` | grid built in 1/mm, `grid_lookup:116-118` |
| DCF — iterative | §5 | — | `iterative_dcf_fa_20190910.m` |
| DCF — implicit (averaging) | §5 | `cudarecon`/`cudarenorm` knorm | `./Auxiliary` divide |
| Kernel — KB | §6.2 | `bessi0` (commented out) | `createKBkernel.m` + Beatty β |
| Kernel — Gaussian | §6.2 | `exp(-d²/0.2)`, `recon.py:129` | — |
| Neighbor search | §6.2 | nearest cell + 4³ box | `knnsearch`, K = 25, cached |
| Oversample + crop | §6.3 | MS=240 → IS=100 | os=3 → NumK·zfill |
| De-apodization | §6.4 | absent | absent |
| FFT shifts / direction | §6.5 | `results.py:99` | `gridrecon:262-264` |
| Readout filter (T2*) | §8 | `gplb=300`, `recon.py:98` | — |
| Spike rejection | §8 | `raw.py:416-427` (zero) | `spiral_human:136-151` (interp) |
| Channel noise norm | §8 | `raw.py:211-215` (raw domain) | `combinecoils_fa` (image domain) |
| Coil combine | §7 | `calcb` + `Σ Re(img·b)` | `combinecoils_fa` |
| Orientation | §6.5 | grid order (+flip in dyn) | `fliplr(permute(...))` |

## 12. Beyond gridding: where compressed sensing picks up

Gridding answers "give me *an* image fast." CS answers "give me the *most consistent* image with what I know about images." The reframing:

```
ρ̂ = argmin_ρ  ‖A ρ − s‖²  +  λ·R(ρ)
```

- **A** = the NUDFT of section 4 — *the same operator gridding approximates.* In practice implemented as a **NUFFT**: apodize → oversampled FFT → kernel interpolation, i.e. gridding's steps reused as a fast, *accurate-by-construction* forward model (kernel error controlled to ~10⁻⁶ rather than accepted as image artifact). Or, at our scales, the exact DFT on GPU (§4a: ~10¹² ops) — no kernel approximation at all.
- **Aᴴ** = the adjoint = gridding without DCF. Inside an iterative solver, *DCF disappears as a concept* — AᴴA's density weighting is handled by the optimization. (At most a DCF-like preconditioner accelerates convergence.)
- **R(ρ)** = prior: sparsity (wavelets, TV) for static; temporal sparsity / low-rank for dynamic ventilation series.
- **Incoherence**: CS wants undersampling artifacts that look like noise in the sparse domain. Spiral/ASAP sampling is naturally good at this — that is what the PSF/coherence measurements in the Obsidian notes (`3D PSF Incoherence Measurement`, `Measuring Coherence in K-space`) are quantifying.

What carries over from the current codebases: trajectory calibration (§9), preprocessing (§8), coil maps (§7 — b becomes part of A in a SENSE-style model), and gridding itself (as Aᴴ, as the solver's initializer x₀, and as the sanity-check baseline). What gets *discarded*: DCF heuristics, kernel-quality compromises, and the one-shot mindset.

## 13. Reading list

| Paper | What it gives |
|-------|---------------|
| O'Sullivan, IEEE TMI 1985 | Gridding (convolution-interpolation) original |
| Jackson et al., IEEE TMI 1991 | Kernel selection analysis, KB as practical optimum |
| Rasche et al. / Meyer et al. | Analytic & Voronoi DCF approaches |
| Pipe & Menon, MRM 1999 | Iterative DCF (Faraz's `iterative_dcf`) |
| Bydder et al., MRM 47:539 2002 | Coil combination (Faraz's `combinecoils_fa`) |
| Fessler & Sutton, IEEE TSP 2003 | NUFFT: gridding as controlled-accuracy operator |
| Beatty et al., IEEE TMI 2005 | Minimal-oversampling gridding; the β formula in Faraz's code |
| Pruessmann et al., MRM 2001 | CG-SENSE: iterative recon with the A/Aᴴ machinery |
| Zwart et al., MRM 2012 | Modern iterative DCF refinement (the paper `iterative_dcf` cites) |
| Lustig et al., MRM 2007 | Sparse MRI — the CS endpoint |
