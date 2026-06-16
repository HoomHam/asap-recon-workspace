# BART `pics` vs our CS — independence test (roadmap step 2)

**Date:** 2026-06-15  **Status:** first sweep landed, ONE open artifact (grid texture).
**Data:** `data/v3_fov250/recon_io` — ACR phantom, static gas, single bin, 1 coil.
**Code (all new, nothing in the past touched):**
`helpers/recon/bart_compare.py`, `bartio.py`, `metrics_v2.py`.

## Why BART
Our CS (`cs_recon.py`) is built on sigpy = Lustig-lab (Frank Ong); the MATLAB
Lustig comparison shares that lineage (personal baseline, already done — ours
won). BART (Uecker/Lustig, C) is the **independent** codebase: different
implementation, native 3D non-Cartesian `pics`. This is the real independence
test.

## Run table

| # | Pipeline | Regularizer | Swept param | Values |
|---|----------|-------------|-------------|--------|
| 1 | ours (sigpy FISTA) | ℓ1-wavelet db4 | t (finalists) | 0.003, 0.01 |
| 2 | ours CG-20 | none (ℓ2) | — | anchor |
| 3 | BART `pics` `-R W:7:0` | ℓ1-wavelet | **λ** | 1e-4, 3e-4, 1e-3, 3e-3, 1e-2 |
| 4 | BART `pics` `-R T:7:0` | TV | **λ** | same 5 |

Primary swept parameter = **regularization strength λ** × 2 regularizer families.
All scored with `metrics_v2` (the fixed metrics — see below). 100 iters each.

## BART setup (the conventions that matter)
- **Trajectory units:** ours is grid-index; BART wants the N-grid to span
  [-N/2, N/2]. Conversion: `bart_traj = traj_rad * N/(2π)` (N=IS=100), dims
  `(3, M, 1)`. ksp `(1, M, 1)`. Single coil → `sens = ones(N,N,N,1)`.
- **cfl/hdr I/O:** `bartio.py` (complex64, column-major) — no dependency on
  BART's toolbox python.
- **Wavelet MUST use ADMM (`-m`).** pics' default FISTA path **diverged to
  ~1e32** on the wavelet prox (bad Lipschitz step from the unnormalized
  ones-sensitivity). `-m` (ADMM) is stable and matches the TV solver. `-e`
  (max-eig scaling) also fixes it; `-S` does not. Fix is in `run_pics`.
- **Orientation:** BART vs finufft differ by axis/flip/conjugate convention.
  `orient_to_ours` picks the flip/transpose of |BART| best-correlating with our
  CG-20 (magnitude, so conjugation drops out). Best cc ≈ 0.92–0.95.

## Results (metrics_v2, lower lfCV = smoother at equal structure)

| run | SNR | lowfreq_CV | edge_sharp | extent_mm | orient_cc |
|-----|-----|-----------|-----------|-----------|-----------|
| ours_cg20 | 20.1 | 0.194 | 6.5 | 172/168/128 | — |
| **ours_wav_t0.003** | 37.2 | **0.157** | 10.3 | 170/168/132 | — |
| **ours_wav_t0.01** | 64.9 | **0.157** | 10.2 | 170/168/132 | — |
| bart_wav_l1e-4 | 26.3 | 0.177 | 8.8 | 172/168/130 | 0.92 |
| bart_wav_l1e-3 | 36.1 | 0.178 | 9.1 | 172/168/130 | 0.92 |
| bart_wav_l1e-2 | 70.5 | 0.175 | 8.6 | 168/168/132 | 0.93 |
| bart_tv_l1e-3 | 68.3 | 0.181 | 9.2 | 172/165/130 | 0.93 |
| bart_tv_l1e-2 | 95.7 | 0.192 | 7.7 | 172/165/122 | 0.95 |

**Reads:**
1. BART **independently reproduces** the phantom — same extent (≈170/168/130
   vs our 170/168/132), orientation cc ≈ 0.92. The geometry agrees across two
   unrelated codebases. That is the independence test passing at the structural
   level.
2. At equal structure, **our wavelet lfCV (0.157) sits below BART's floor
   (~0.175–0.18)** across the whole λ sweep — ours is smoother. BUT see the
   open artifact: BART's CV is inflated by a texture ours lacks, so do NOT yet
   bank this as "ours wins." Confounded.
3. BART SNR climbs with λ the same way ours does; high-λ SNR (95.7 at TV 1e-2)
   is the over-smoothing inflation the new metric still reports honestly via
   the falling edge_sharp (7.7 < 9.2) — the guard working.

## BART grid texture — it is the REAL ACR resolution insert (see VERDICT below)
**Superseded read (kept for the trail):** I first called the cross-hatch a
BART-specific texture / artifact and chased a cause (DCF, Toeplitz, kernel,
fftmod — all ruled out). WRONG framing. It is the phantom's resolution/grid
insert at z≈50, which BART RESOLVES and our wavelet BLURS. The "ours clean, BART
gridded" reads in `wavelet_twoway_montage.png` / `texture_compare.png` were
ours OVER-SMOOTHING real structure, not BART adding it. See VERDICT CORRECTED.

**Scalars failed to capture it — twice — don't trust them here:**
1. First HF probe (0.69) was a code BUG (FFT over the singleton trailing dims of
   an unsqueezed 16-D cfl).
2. Fixed HF band (r 6-22) and spectral peakiness (max/median) both read ~equal
   for BART vs ours (HF 0.05-0.07, peakiness ~6.7) — because the phantom's own
   mid-freq spectrum overlaps the grid band. These scalars CANNOT separate a
   regular grid from random mottle on this object. The eye can.

**Cause: BART's nufft gridding internals — UNRESOLVED.** Ruled out by test:
traj over-range (|k|max 47 < 50), Toeplitz embedding (`-U` no change), KB kernel
width/oversampling (`-o3 -w8` no change), fftmod (`bart fftmod` no-op), and DCF /
unweighting (ours-WITHOUT-DCF is still grid-free, so it is NOT the unweighted-LS
density signature — an earlier guess, disproven). `-p` does not wire per-sample
noncart DCF in this build. Next ideas: BART rolloff/apodization correction; a
point-source PSF dump (`nufft` of a delta traj) to see the kernel's spatial
replica.

## FINAL VERDICT (2026-06-15, z-matched) — all three resolve the grid; ours is sharp
The cross-hatch is the ACR **resolution/grid insert** (z≈50-54), NOT an artifact.

**The methodological key (Hooman): the phantoms are slightly different SIZES, so
a flip-only orientation aligns the bulk but NOT the thin insert slice.** Fitted
z-affine to ours: BART scale **1.02** (+1.5 slice), Lustig scale **0.93** (Lustig
is 7% smaller in z). A fixed-z comparison is therefore INVALID — at z=50 BART sat
on its grid peak while ours/Lustig were on uniform planes, which faked a "BART
resolves, ours blurs" result.

After z-affine registration (`z_register_compare.py`) the grid lines up in the
SAME column across all rows (`zreg_finez_montage.png`), and at the matched insert
slice zoomed (`grid_slice_zoom_z52.png`):
- **ours wav t0.003: grid sharply resolved.**
- **Lustig TV w0.001: grid sharply resolved.**
- **BART wav l3e-4: grid present but softer here.**

**So all three resolve the real ACR grid; ours is as sharp or sharper than BART.**
Every prior "winner" call (ours-wins-lfCV; then BART-resolves-best) was an artifact
of mismatched slices / smoothness-rewarding metrics — RETRACTED. The pipelines are
genuinely comparable at matched resolution.

**Six-way (Faraz/Steve/our CG/our wav/BART/Lustig), z-registered to our CG**
(`zreg_sixway_montage.py` → `zreg_sixway_montage.png`): grid insert aligns at
z≈52 across all six; all resolve it. z-scales rel. our CG: Faraz **1.16**,
Steve 1.01, our wav 0.98, BART 1.04, Lustig **0.94** — phantom size spans **~23%**
across pipelines (real, pipeline-dependent FOV/scaling difference). Faraz slice
corr is only 0.75 (in-plane recon difference, not misregistration).

**Mandatory protocol going forward:** z-affine register (scale+shift) before ANY
slice comparison — sizes differ up to ~23% across pipelines. Judge resolution at the insert slice by eye;
metrics_v2 is blind to it. Lesson (eye beats single scalar, 3rd+ time):
`[[eye_vs_metric]]`. Figures: `zreg_finez_montage.png`, `grid_slice_zoom_z52.png`,
`slice_matched_montage.png`, `resolution_sweep_z50.png`.

## Lustig TV sweep + three-way TV (added 2026-06-15)

Earlier the Lustig baseline was ONE point (TVWeight=xfmWeight=0.01, both terms on,
identity XFM — neither a clean TV nor a wavelet run) scored with the OLD metric.
Fixed both: swept pure TV (`run_lustig_sweep.py` → `run_cs_sweep.m`, xfmWeight=0)
and scored with `metrics_v2`, matching BART's λ sweep.

**Wavelet stays ours-vs-BART only:** Lustig's `Wavelet` operator is `FWT2_PO`
(2D). On the 3D volume it's the documented 2D-on-3D bug. A real 3D-wavelet Lustig
run needs a new operator — not built.

Three-way TV on ONE ruler (`metrics_v2`, matched plane via orientation search),
`tv_threeway.py` → `recon_io/tv_threeway_{montage.png,metrics.json}`:

| pipeline (TV) | lfCV range | note |
|---|---|---|
| ours (PDHG) | 0.185–0.194 | clean, no texture |
| BART (`-R T`) | 0.180–0.192 | + cross-hatch moiré (the artifact) |
| Lustig (fnlCg) | 0.178–0.206 | softens/dims as λ↑; NCG `obj` flat across 15 iters (final≈init) |

**Verdict:** on TV the three are a **wash** (~0.18 lfCV) — TV is commoditized
across implementations. Differentiation lives in the **wavelet** prior (ours
0.157 < BART 0.177), where Lustig can't compete (2D operator). Lustig's NCG
barely moves from its DCF-seeded init (objective flat ~1094, RMS 0.0335 all
iters) — confirms the stalled-solver / "final ≈ init" story.

**Metric caveat (new):** `edge_sharp` is inflated by texture/noise — BART's moiré
and Lustig's grain both raise it. It guards against blur, NOT texture; don't read
their higher edge_sharp as finer detail.

New code: `lustig_oneshot/run_cs_sweep.m`, `run_lustig_sweep.py`,
`recon/tv_threeway.py` (full 48-orientation search for Lustig). Originals
untouched.

## Files written
- `recon_io/bart/bart_vs_ours_montage.png` — central-slice montage, shared axis
- `recon_io/bart/bart_compare_metrics.json` — all metrics_v2 numbers
- `recon_io/bart/{traj,ksp,sens,bart_*}.{cfl,hdr}` — BART workspace

## Reproduce
```
cd helpers/recon
../.venv/bin/python bart_compare.py ../../data/v3_fov250/recon_io \
    --bart /Users/hoomham/bin/bart-src/bart --iters 100
```
BART binary built from source at `/Users/hoomham/bin/bart-src/bart`
(gcc-15, OPENBLAS=1, PNG=0; macOS needs `gmake`, a `cblas_openblas.h`→`cblas.h`
symlink in the openblas include dir, and `FFTW_BASE`/`CPATH` set to brew prefixes).
```
