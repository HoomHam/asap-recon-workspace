# Faraz's Recon — Code Map

**Date:** 2026-06-10 (from full code read, same session as `Recon_Overview_Steve.md`)
**Repo:** `/Users/hoomham/Hooman/Work/Codes/2023_Faraz_Recon/`
**Hooman's fork:** `workspace/codes/2023_Faraz_Recon_HH/`
**Function map (deep):** Obsidian `Action/MRI/ASAP Recon/ASAP Recon Faraz Approach.md`
**Steve vs Faraz diff:** `reference/Recon_Comparison_StaticGas.md`

## File Roles

| File | Role |
|------|------|
| `spiral_human_20240227.m` | Top-level driver script (cell-mode, per-study): load → spike filter → peak/bin → combine → grid recon → coil combine → orient → save. The de facto entry point for human studies |
| `recon_20210622.m` | Generic loader/recon wrapper: parses TWIX (`mapVBVD.m` + `parseheader.m`), exposes all recon Options (os=3, k=5, zfill=1, dcf mode, calib files); `'recon',0` = load-only (how spiral_human uses it) |
| `fa_spiral_dyn_recon.m` | Protocol-dispatching recon: fancy_v2/v3 path (current) → reshape, average repeats, grid, combine; legacy paths for older protocols (GA rotation, radial-interleaved DP, sequential) |
| `gridrecon_fa_20230113.m` | Core gridder: neighbor lookup, KB kernel table, DCF, accumulate, normalize by Auxiliary, FFT, crop. 2D + 3D paths |
| `grid_lookup_20230113.m` | Grid cache: computes or loads `Ind`/`Dist`/`wi` keyed on (trajectory, imgsize, fov, os, kernel, β, zfill) → `grid_lookup_20230113.mat` |
| `createKBkernel.m` | Kaiser–Bessel lookup table: `besseli(0, β·sqrt(1−(2u/w)²))`, 10000 entries, peak-normalized |
| `iterative_dcf_fa_20190910.m` | Pipe–Menon iterative DCF (grid→degrid, `wi ← wi/(C·wi)`), 5 iterations, `sqrt(KB)` kernel. Header: "3D was modified and needs to be double-checked" |
| `combinecoils_fa.m` | Coil combine: per-channel corner-noise normalize, `b = Σ_rep img / Σ_rep \|img\|`, `real(b'·img)` per voxel (Bydder MRM 2002) |
| `loadtrajectory3D.m` | Measured-trajectory lookup from `calibrations_3D_20220308.mat`, matched on freq/FOV (±5%), nsamples, nleaves, nreps, imgsize, orientation, nucleus |
| `calcspiralgrad_20190918.m`, `calcradialtraj_20220713.m` | Analytic gradient/trajectory design (acquisition side) |
| `mapVBVD.m`, `twix_map_obj.m`, `read_twix_hdr.m`, `parseheader.m` | Siemens TWIX parsing stack |
| `resize_fa.m`, `wrapImage_fa.m`, `threshold_faraz.m`, `printParam.m`, `ticker.m` | Utilities (resize, display wrap, threshold, param print, progress) |
| `Fit4LorentzianCplxPh_Dixon_*.m`, `fitVoigt.m`, `evalVoigt.m` | Spectra fitting / Dixon (dissolved phase — out of static-gas scope) |

## Data Flow (static gas recon path)

```
Siemens .dat ──► recon_20210622(files,'recon',0)    mapVBVD parse + header → data struct, twix raw
                  │
                 reshape → [nsamples, nch, nacq]
                  │  split off leading spectra acqs if bSpectra
                  │
                 spike filter      movmean(3) ratio > 1+3σ → replace by neighbor mean
                  │
                 findpeaks(|first sample|)      breath detection (out of scope for static)
                  │
                 accumulate into rawdata2       average repeats of same interleave/rotation
                  │  [nsamples, nleaves, nreps, nch(, nphases)] ./ weights
                  │  optional FOV shift: phase ramp exp(i·2π·k·Δx)
                  │
calib .mat ────► loadtrajectory3D()             measured k(x,y,z), physical 1/mm units
                  │
                 grid_lookup_20230113()         Ind/Dist via knnsearch (K = k² = 25),
                  │                             KB table, iterative DCF wi — all cached to .mat
                  │
                 gridrecon_fa_20230113()        M = raw·wi; per sample:
                  │                             k(Ind) += M·KB(Dist); Aux(Ind) += wi·KB(Dist)
                  │                             k ./= Aux (NaN→0)
                  │
                 fftshift(fft(fftshift(k,i),[],i),i) per dim     forward FFT, os·matsize grid
                  │
                 crop center matsize·zfill
                  │
                 combinecoils_fa()              noise-normalize channels, b map, real(b'·img)
                  │
                 fliplr(permute(img,[1 3 2 4]))  orientation to display convention
                  │
                 save img_dyn_*.mat (+ _ch per-channel version)
```

## Entry Points

| Task | Start at |
|------|----------|
| Human study recon (current) | `spiral_human_20240227.m` — edit file list + flags, run cells top to bottom |
| Generic/older protocol recon | `recon_20210622(files, Options)` with `'recon',1` → dispatches `fa_spiral_dyn_recon` |
| Static single-bin gas | spiral_human path with `nphases=1`, `img_flag=0` (binning collapses to all-data average) |
| Grid/DCF inspection | `grid_lookup_20230113.mat` entries (trajectory-keyed struct array) |

## Compute Model

CPU MATLAB. Heavy steps amortized: `knnsearch` neighbor lookup (parfor for grids >180–200/dim) and iterative DCF computed **once per trajectory** and cached to `grid_lookup_20230113.mat` (441 MB — grows with every new trajectory/matrix combo). Gridding itself: serial loop over samples, vectorized over channels/phases; single-precision accumulators (`k_real`/`k_imag` separate, combined after). Per-phase gridding operates on the pre-averaged `nsamples×nleaves×nreps` set, not the full acquisition — the key reason it stays tractable without GPU (see comparison doc, compute-cost note).

## What's NOT in Faraz's Code (vs Steve)

- No GPU path — all CPU MATLAB
- No raw-domain channel noise normalization (image-domain corner-noise instead)
- No global FID rephasing (phase handled entirely in coil-combine b map)
- No readout apodization / line-broadening filter (Steve's gplb has no equivalent)
- No initial-point kill (keeps all samples incl. ramp; trusts trajectory calibration)
- No low-SNR full-image exclusion ranges (manual `range` selection in driver instead)
- Full k³ kernel support missing — knnsearch K = k² = 25 truncates KB ball to ≈1.8 cells
- No de-apodization (same omission as Steve; "partial de-apodization" comment is mislabeled normalization)

## Key Magic Numbers

| Value | Where | Meaning |
|-------|-------|---------|
| `os=3, k=5, zfill=1` | spiral_human / recon_20210622 defaults | Oversampling, kernel size, zero-fill |
| `β = π·√(k²/os²·(os−0.5)²−0.8)` ≈ 12.78 | grid_lookup / gridrecon | Beatty Eq 5 KB shape |
| `klength=10000` | createKBkernel call sites | Kernel lookup-table resolution |
| `iter=5` | gridrecon/grid_lookup | DCF iterations |
| `dcftable = kernel.^(1/2)` | gridrecon/grid_lookup | sqrt(KB) as DCF kernel (ad hoc) |
| `K = kernelsize²` (=25) | knnsearch calls | 3D neighbor count — **k², not k³** (quirk) |
| `1 + 3σ` over `movmean(…,3)` | spiral_human spike filter | Spike threshold |
| `NoiseSize=6` | combinecoils_fa | Corner noise-region size |
| grid >180 (lookup) / >200 (gridrecon) per dim | grid_lookup vs gridrecon | Threshold switching knnsearch → per-axis + parfor (inconsistent between the two files) |
| `tolerance=0.05` | loadtrajectory3D | Freq/FOV match tolerance for calibration entries |
