---
tags: [action, mri, compressed-sensing, asap, reconstruction, matlab, work, claude-supported]
updated: 2026-05-22
---

# CS–ASAP Adaptation

## Overview

Adapting Lustig's 3D spiral CS reconstruction (`Work/Codes/2025_Xe129_CS`) to work with ASAP-acquired data. Phase 1: static 3D ACR water phantom — **done and running**. Phase 2: dynamic 4D Xe129.

- ASAP trajectory reference: [[ASAP-Reconstruction]]
- Full code pipeline: [[CS-ASAP-Pipeline]]
- Current status: [[CS-ASAP-Status]]
- Theory background: [[MRI-Learning]]
- Working code lives in: `Work/Analysis/2025-09-24_ACR/`

---

## The CS Code (2025_Xe129_CS)

Two main scripts:
- `spiral3d_cs_Dt.m` — temporal TV regularization
- `spiral3d_cs_DzDt.m` — spatial-z + temporal TV

**External dependencies (must be on MATLAB path):**

| Toolbox | Purpose | Setup |
|---------|---------|-------|
| `irt/` (Fessler IRT) | NUFFT operator | run `setup.m` |
| `mapVBVD-main/` | Siemens raw data reader | add with subfolders |
| `sparseMRI_v0.2/utils/` | `voronoidens` for DCF | add to path |

### Core interface
```matlab
FT = NUFFT4D(ktrajs, 1, ph, 0, imSize, 2);
im_dc = FT' * (kdatas .* kcomps);   % density-compensated gridding
res = fnlCg(res, param);             % CS iterations
```

---

## Phase 1: 3D Static ACR Phantom

### Changes from original code

| What | Original | New |
|------|----------|-----|
| Data source | `xe129_frames.mat` | load from raw ASAP .dat |
| Image size | [256, 128, 256, 16] (4D) | [N, N, N] (3D) |
| NUFFT | NUFFT4D | NUFFT3D (or 4D with 1 frame) |
| TV regularization | temporal (Dt) | spatial only, or TV=0 + wavelet |
| Trajectory | pre-packed in .mat | generated from kasap.c / MATLAB equivalent |
| Save paths | `C:\Users\P53-LOCAL\...` | `~/Hooman/Work/Analysis/[date]_ACR/` |
| gamma | 42.577 MHz/T | same (proton for ACR) |

### Actual data packaging (confirmed, `cs_spiral_gpt_ACR_hoom_20250924.m`)

```matlab
% 1. Load raw data
[data, twixs] = load_rawdata_20250816(files, 'nStudies', 1:1);

% 2. Extract scan params
nsamples=512; nleaves=26; nreps=32; FOV_mm=250; nchannels=1;

% 3. Flatten TWIX → [nsamples × nchannels × nblocks]
raw_flat = squeeze(twixs.file1.image(:,:,:,:,:,:,:,:,:,:,:));
rawdata  = reshape(raw_flat, nsamples, nchannels, []);  % [512×1×832]

% 4. Load trajectory from XYZ calibration scans
[KSpaceCoor, ~, ~, ~] = loadtrajectory3D('BuildFromXYZ', struct(...
    'RO', ROFileName, 'PE', PEFileName, 'SS', SSFileName, ...
    'OutDir', outdir, 'Version', version));
% KSpaceCoor: [Nseg × 3] in 1/mm; tiled to match rawdata2 if needed

% 5. Scale trajectory to grid
kmax_meas  = max(sqrt(sum(KSpaceCoor.^2, 2)));
kmax_nyq   = matsize / (2*FOV_mm);
KSpaceCoor = (kmax_nyq / kmax_meas) * KSpaceCoor;

% 6. Aggregate into [nsamples × nleaves × nreps × ncoils], flatten
raw_gp   = zeros(nsamples, nleaves, nreps, nchannels, 'single');
% ... (loop accumulates per-leaf, per-rep averages)
rawdata2 = reshape(raw_gp, [], nchannels);   % [425984 × 1]

% 7. Save intermediate
a = [real(rawdata2), imag(rawdata2), KSpaceCoor];   % [425984 × 5]
save ACR_data a
```

### k-trajectory normalization (confirmed two-step)

**Step 1 — MATLAB** (physical units → grid-normalized 1/mm):
```matlab
alpha      = (matsize/2) / (FOV_mm * kmax_meas);   % 1-voxel margin variant
KSpaceCoor = alpha * KSpaceCoor;   % now fills grid to ~Nyquist
```

**Step 2 — Python notebook** (1/mm → radians for torchkbnufft):
```python
dk_mag = np.sqrt(ktrajx**2 + ktrajy**2 + ktrajz**2)
ktraj[0,:] = ktrajx / dk_mag.max() * np.pi   # → [-π, π]
ktraj[1,:] = ktrajy / dk_mag.max() * np.pi
ktraj[2,:] = ktrajz / dk_mag.max() * np.pi
```

**Step 3 — MATLAB CS script** (radians → IRT NUFFT range):
```matlab
dkx = ktrajs(frame,1,:);
dkx = dkx(:) / (2*pi);   % → [-0.5, 0.5] for NUFFT3D
```

### DCF

Currently: **torchkbnufft pipe method** (Python notebook, `calc_density_compensation_function`).  
Available but unused in CS path: iterative KB DCF from `grid_lookup_20230113` (MATLAB, `iterative_dcf_fa_20190910`).  
🔲 Open: compare the two — iterative KB may give better CS init.

---

## Raw Data — Confirmed

- **Scanner:** Siemens (3T proton, ACR phantom)
- **Raw file:** `meas_MID00123_FID12098_fa_spiral_dyn_fancy_v2_20230131.dat`
- **mapVBVD access:** `squeeze(twixs.file1.image(:,:,:,...))` → reshape to `[nsamples × nchannels × nblocks]`
- **Dimension ordering:** samples fast, channels middle, interleaf-blocks slow
- **Ncoils:** 1 (single channel)
- **ACR scan params:** NI=26, NPTS=512, NREPS=32, FOV=250mm, matrix=80, Ntotal=425,984
- **Trajectory:** built fresh from XYZ calibration `.dat` files via `build_calibration_from_xyz`

---

## File Structure (actual, as of 2026-05-22)

```
Work/Analysis/2025-09-24_ACR/        ← working folder (data + code together)
  cs_spiral_gpt_ACR_hoom_20250924.m  ← MATLAB: .dat → ACR_data.mat
  spiral3d_frames_mat_hoom.ipynb     ← Python: ACR_data.mat → ACR_test.mat
  spiral3d_cs_3D_hoom.m              ← MATLAB CS: ACR_test.mat → CS recon
  spiral3d_cs_3D_hoom_Wavelet.m      ← Wavelet variant
  ACR_data.mat                       ← intermediate [425984×5]
  ACR_test.mat                       ← CS input (ktrajs/kdatas/kcomps)
  calibrations_3D_20220308.mat       ← legacy calibration (reference)
```

No separate `2026_Xe129_CS` folder yet — all ACR work lives in `2025-09-24_ACR`.

---

## Phase 2: Dynamic 4D Xe129 (Future)

- `kdatas` → `[Nframes × Ntotal]`
- `ktrajs` — shared or per-frame
- Reintroduce `TVOPDt` / `TVOPDzDt`
- `imSize` → `[N, N, N, Nframes]`
- gamma → 11.777 MHz/T (Xe129)

---

## Log

### 2026-05-22 — Wiki page created
- Cowork session planning the 2026_Xe129_CS project
- Source: spiral3d_cs_Dt.m, ASAP Recon Pipeline.md, kasap.c

### 2026-05-22 — Updated with confirmed pipeline
- Source: Cowork session — read all scripts + mat files in `2025-09-24_ACR/`
- ✅ Raw data format confirmed (see above)
- ✅ Normalization chain confirmed (MATLAB → Python → MATLAB CS)
- ✅ CS pipeline already running end-to-end on ACR phantom
- ✅ File structure updated to reflect reality
- 🔲 Open: compare KB-iterative vs pipe DCF quality for CS init
- 🔲 Open: tune CS parameters (TVWeight, xfmWeight, Itnlim)
- 🔲 Open: Phase 2 — adapt for Xe129 (4D, gamma correction)
