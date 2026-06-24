# CS Adaptation Plan — ASAP → 3D Static Recon

**For:** Claude Code CLI sessions working in 2026_Xe129_CS  
**Source:** spiral3d_cs_Dt.m, spiral3d_cs_DzDt.m (2025_Xe129_CS), ASAP Recon Pipeline.md

---

## Goal (Phase 1)

Adapt the Lustig CS reconstruction (`2025_Xe129_CS`) to reconstruct a **static 3D ACR water phantom** from ASAP-acquired data. No temporal dimension yet — that is Phase 2.

---

## The CS Code (2025_Xe129_CS)

Two main scripts:
- `spiral3d_cs_Dt.m` — CS with temporal finite-difference (Dt) TV regularization
- `spiral3d_cs_DzDt.m` — CS with both spatial-z (Dz) and temporal (Dt) TV

Supporting classes:
- `@TVOPDt/` — temporal TV operator
- `@TVOPDzDt/` — spatial+temporal TV operator

**External dependencies (must be on MATLAB path):**
- `irt/` (Jeff Fessler's IRT toolbox) — provides NUFFT. Run `setup.m` first.
- `mapVBVD-main/` — Siemens raw data reader (add with subfolders)
- `sparseMRI_v0.2/utils/` — provides `voronoidens` for DCF

---

## How the Original Code Works

```matlab
% 1. Load pre-packaged data
data_frames = load('xe129_frames.mat');
ktrajs = data_frames.ktrajs;   % k-space coordinates
kdatas = data_frames.kdatas;   % raw k-space signal  
kcomps = data_frames.kcomps;   % density compensation

% 2. Build NUFFT operator
FT = NUFFT4D(ktrajs, 1, ph, 0, imSize, 2);

% 3. DC image (density-compensated gridding)
im_dc = FT' * (kdatas .* kcomps);

% 4. CS iterations
param.FT = FT;
param.TV = TVOPDt;        % temporal TV
param.data = kdatas;
param.TVWeight = 0.01;
res = fnlCg(res, param);  % nonlinear conjugate gradient
```

---

## What Changes for 3D Static ASAP ACR Phantom

### 1. Data loading
Replace `load('xe129_frames.mat')` with a script that:
- Reads raw Siemens data using mapVBVD
- Loads the ASAP k-space trajectory (from kasap.c output or MATLAB equivalent)
- Computes DCF for the ASAP trajectory
- Packages into `ktrajs`, `kdatas`, `kcomps`

### 2. Image dimensions
Original: `imSize = [256, 128, 256, 16]` (4D: x, z, y, frames)  
New: `imSize = [N, N, N]` (3D: confirm matrix size from scan params)  
Use `NUFFT3D` (or `NUFFT4D` with 1 frame) accordingly.

### 3. TV regularization
Original uses temporal TV (`TVOPDt`) — meaningless for static phantom.  
New: use spatial TV only. Options:
- Set `TVWeight = 0` and use wavelet (`xfmWeight > 0`) only
- Or implement/use a 3D spatial TV operator

### 4. k-trajectory normalization
NUFFT expects trajectories normalized to [-0.5, 0.5].  
kasap.c outputs in physical units (cycles/m). Must normalize:
```matlab
ktrajs = ktrajs / (2 * kmax);  % where kmax = ms/(2*FOV)
```
Verify exact normalization convention used by the IRT NUFFT.

### 5. DCF computation
For ASAP (non-Cartesian), compute DCF iteratively:
```matlab
% Using voronoidens or iterative method
% Warm-start from analytical approximation, then 3-5 CG iterations
```

### 6. Fix hardcoded paths
All save paths in original code are Windows (`C:\Users\P53-LOCAL\...`).  
Replace with relative paths or Mac paths under `~/Hooman/Work/Analysis/`.

---

## Raw Data Format

> **TODO — fill in before starting CLI session:**
> - Scanner: Siemens [model?]
> - Raw data file: `.dat` file read by mapVBVD?
> - Key variables after mapVBVD load: [twix.image.data dimensions?]
> - Data ordering: how do reps/interleaves/samples map to array indices?
> - Number of coils: [?]
> - Acquisition parameters for ACR scan: NI=?, NPTS=?, NREPS=?, FOV=?, matrix=?

---

## File / Folder Conventions

```
Work/Codes/2026_Xe129_CS/         ← git repo, new work lives here
  CLAUDE.md                        ← project context for CLI sessions
  wiki_ASAP_trajectory.md          ← ASAP reference (this + companion)
  wiki_CS_adaptation.md            ← this file
  recon_ASAP_3D.m                  ← main adapted reconstruction script
  load_asap_data.m                 ← raw data → ktrajs/kdatas/kcomps
  compute_dcf_asap.m               ← DCF for ASAP trajectory
  @TVOPspatial/                    ← 3D spatial TV operator (if needed)

Work/Analysis/[date]_ACR/          ← data lives here, NOT in Codes/
  raw/                             ← original .dat files
  results/                         ← output images
```

---

## Phase 2 (Future)

Once 3D static works: add temporal dimension for dynamic Xe129.
- `kdatas` becomes `[Nframes × Ntotal]`
- `ktrajs` may vary per frame or be shared
- Reintroduce `TVOPDt` or `TVOPDzDt`
- Image size: `[N, N, N, Nframes]`
