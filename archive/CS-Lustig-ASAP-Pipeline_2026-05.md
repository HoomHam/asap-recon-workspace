---
tags: [action, mri, compressed-sensing, asap, pipeline, code, matlab, python, work, claude-supported]
updated: 2026-05-22
---

# CS–ASAP Pipeline (Code Chain)

## Overview

Full end-to-end code chain for running CS reconstruction on ASAP data, as confirmed from the working ACR phantom implementation in `Work/Analysis/2025-09-24_ACR/`.

- Adaptation context: [[CS-ASAP-Adaptation]]
- Trajectory + data format: [[ASAP-Reconstruction]]
- Current status: [[CS-ASAP-Status]]

---

## MATLAB Path Setup (required before running any stage)

```matlab
run('/Users/hoomham/Hooman/Work/Codes/2025_CS/irt/setup.m')
addpath(genpath('/Users/hoomham/Hooman/Work/Codes/2025_CS/sparseMRI_v0.2'))
% mapVBVD.m lives in 2025-09-24_ACR/ — already on path when running from there
```

`2025_CS/` structure:
```
2025_CS/
  irt/              ← Fessler IRT: setup.m, NUFFT, etc.
  sparseMRI_v0.2/   ← Lustig CS: fnlCg.m, @NUFFT3D, @TVOP3D, utils/voronoidens.m
```

---

## Pipeline at a Glance

```
Siemens .dat file
      ↓  load_rawdata_20250816  (MATLAB)
twixs + allimages (metadata)
      ↓  loadtrajectory3D / build_calibration_from_xyz  (MATLAB)
KSpaceCoor [Nseg×3] in 1/mm  +  rawdata2 [Ntotal×Ncoils]
      ↓  scale + aggregate  (MATLAB)
ACR_data.mat  [425984×5]  — Re/Im/kx/ky/kz
      ↓  spiral3d_frames_mat_hoom.ipynb  (Python)
ACR_test.mat  — ktrajs[1×3×N], kdatas[1×N], kcomps[1×N]
      ↓  spiral3d_cs_3D_hoom.m  (MATLAB)
CS reconstructed 3D volume
```

---

## Stage 1 — Load Raw Data (MATLAB)

**Script:** `cs_spiral_gpt_ACR_hoom_20250924.m` (top section)  
**Function:** `load_rawdata_20250816(files, 'nStudies', 1:1)`

**What it does:**
- Opens `.dat` via `mapVBVD`
- Extracts header: `nCol`(=NPTS), `nLin`(=NI), `nRep`(=NREPS), `fovPE`, `MatSize`, `dwelltime`, `nCha`
- Detects trajectory type from filename (`spiral`, `fancy` → 3D)
- Returns `allimages` (metadata struct) and `twixs` (raw TWIX handles)

**Raw data access:**
```matlab
raw_flat = squeeze(twixs.file1.image(:,:,:,:,:,:,:,:,:,:,:));
rawdata  = reshape(raw_flat, nsamples, nchannels, nblocks);
% shape: [512 × 1 × 832]   for ACR (1 coil, 832 interleaf blocks)
```

---

## Stage 2 — Trajectory Calibration (MATLAB)

**Function:** `loadtrajectory3D('BuildFromXYZ', struct(...))`  
→ calls `build_calibration_from_xyz(ROFile, PEFile, SSFile, outdir, version)`

**What it does:**
1. Loads 3 calibration `.dat` files (X/Y/Z gradient axes)
2. Reads 2 slices × 3 gradient states (on/inv/off)
3. Unwraps phases, subtracts baseline, applies two-slice ±D difference-of-differences formula:
   `k = [(φ1_on − φ1_inv) − (φ2_on − φ2_inv)] / (4·D·2π)`
4. Normalizes interleaves, combines coils (mean), applies loess smoothing
5. Saves `cal` struct to `.mat` (filename encodes version, FOV, NS, NL, etc.)
6. Returns `KSpaceCoor [Nseg×3]` in **1/mm**

**Key calibration struct fields:**
```
cal.kx, cal.ky, cal.kz   — [Nseg×1] in 1/mm
cal.kmax_measured         — scalar in 1/mm
cal.nsamples, nleaves     — 512, 26 (per-rep)
cal.FOV_mm                — 250
cal.dwell_us              — 5
```

**Trajectory scaling** (in main script):
```matlab
kmax_meas  = max(sqrt(sum(KSpaceCoor.^2, 2)));
kmax_nyq   = matsize / (2 * FOV_mm);          % 80/(2×250) = 0.16 1/mm
alpha      = kmax_nyq / kmax_meas;
KSpaceCoor = alpha * KSpaceCoor;              % now fills Nyquist grid
```

---

## Stage 3 — Aggregate & Save Intermediate (MATLAB)

**Script:** `cs_spiral_gpt_ACR_hoom_20250924.m` (aggregation section)

**Aggregation** (time-averaging across repeated blocks):
```matlab
raw_gp = zeros(nsamples, nleaves, nreps, nchannels, 'single');
% fill loop: each block → (leaf index, rep index) via modular indexing
rawdata2 = reshape(raw_gp, [], nchannels);   % [425984 × 1]
```

**Save:**
```matlab
a = [real(rawdata2), imag(rawdata2), KSpaceCoor];  % [425984 × 5]
save ACR_data a
```

`ACR_data.mat` is the handoff between MATLAB and Python.

---

## Stage 4 — Repackage for CS (Python)

**Notebook:** `spiral3d_frames_mat_hoom.ipynb`  
**Libraries:** numpy, scipy.io, torchkbnufft

**What it does:**
```python
data = sio.loadmat('ACR_data.mat')['a']   # [425984 × 5]

ktrajx = data[:,2]; ktrajy = data[:,3]; ktrajz = data[:,4]
rawdata = data[:,0] + 1j*data[:,1]

# Normalize trajectory to [-π, π]
dk_mag = np.sqrt(ktrajx**2 + ktrajy**2 + ktrajz**2)
ktraj = np.zeros((3, 425984))
ktraj[0,:] = ktrajx / dk_mag.max() * np.pi
ktraj[1,:] = ktrajy / dk_mag.max() * np.pi
ktraj[2,:] = ktrajz / dk_mag.max() * np.pi

# DCF (pipe method)
dcomp = tkbn.calc_density_compensation_function(ktraj=ktraj, im_size=(256,64,256))

# Save
sio.savemat('ACR_test.mat', {'ktrajs':[ktraj], 'kdatas':[rawdata], 'kcomps':[dcomp]})
```

**Output shapes in `ACR_test.mat`:**
| Variable | Shape | Units |
|----------|-------|-------|
| `ktrajs` | `[1, 3, 425984]` | radians [−π, π] |
| `kdatas` | `[1, 425984]` | complex a.u. |
| `kcomps` | `[1, 425984]` | pipe DCF weights |

⚠️ **BUG: `im_size=(256,64,256)` is wrong.** This is leftover from a cardiac dataset. The MATLAB CS reconstruction runs at `[90,90,90]` (or `MS_recon_desired × 3`). The DCF is being computed for the wrong grid. **Must fix before any serious evaluation of CS output quality.** Decision on how to restructure the notebook is pending.

---

## Stage 5 — CS Reconstruction (MATLAB)

**Script:** `spiral3d_cs_3D_hoom.m`  
**Dependencies:** IRT toolbox (`setup.m`), `sparseMRI_v0.2`

```matlab
% Load
ktrajs = data_frames.ktrajs;   % [1×3×425984]
kdatas = data_frames.kdatas;   % [1×425984]
kcomps = data_frames.kcomps;   % [1×425984]

% Normalize trajectory to [-0.5, 0.5] for IRT NUFFT
dkx = ktrajs(1,1,:)/(2*pi);
dky = ktrajs(1,2,:)/(2*pi);
dkz = ktrajs(1,3,:)/(2*pi);
k   = [dkx(:), dky(:), dkz(:)];

% Setup NUFFT and density-compensated init
imSize = [90, 90, 90];   % MS_recon_desired
FT = NUFFT3D(k, 1, 1, 0, imSize, 2);
w  = kcomps(1,:)' / max(kcomps(:));
data = kdatas(1,:)' / max(abs(kdatas(:)));
im_dc = FT' * (data .* w);   % gridded init

% CS parameters
param.FT        = FT;
param.TV        = TVOP3D;
param.data      = data;
param.TVWeight  = 0.01;
param.xfmWeight = 0.01;
param.Itnlim    = 15;

% Iterate
res = im_dc;
for n = 1:15
    res = fnlCg(res, param);
end
```

---

## Key Open Questions

- 🔲 What `im_size` should the Python notebook use? Currently `(256,64,256)` which doesn't match the MATLAB reconstruction grid.
- 🔲 Compare KB-iterative DCF (MATLAB) vs pipe DCF (Python) — does it affect CS convergence?
- 🔲 Tune CS parameters: `TVWeight`, `xfmWeight`, `Itnlim` for ACR phantom.
- 🔲 Wavelet variant (`spiral3d_cs_3D_hoom_Wavelet.m`) — when to use vs pure TV?

---

## Log

### 2026-05-22 — Page created
- Source: Cowork session — full code trace of `2025-09-24_ACR/` scripts + mat file inspection
- Confirmed entire pipeline is functional end-to-end
