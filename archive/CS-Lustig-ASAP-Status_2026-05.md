---
tags: [action, mri, compressed-sensing, asap, status, work, claude-supported]
updated: 2026-05-22
---

# CS–ASAP — Current Status

> **This is the orientating page.** Upload this + the relevant detail page(s) at the start of a Claude session to get full context fast.

---

## Where Things Stand (2026-05-22)

**Phase 1 (ACR phantom, static 3D): ✅ End-to-end pipeline working.**

The full chain from raw `.dat` → CS-reconstructed 3D volume is running. Code lives in `Work/Analysis/2025-09-24_ACR/`. No new infrastructure needed for basic operation.

**Phase 2 (Xe129, dynamic 4D): 🔲 Not started.**

---

## What's Working

- `cs_spiral_gpt_ACR_hoom_20250924.m` — loads raw ASAP data, builds trajectory from XYZ calibration, aggregates, saves `ACR_data.mat`
- `spiral3d_frames_mat_hoom.ipynb` — converts `ACR_data.mat` → `ACR_test.mat` (CS format, torchkbnufft DCF)
- `spiral3d_cs_3D_hoom.m` — runs NUFFT3D + `fnlCg` CS iterations on `ACR_test.mat`
- `spiral3d_cs_3D_hoom_Wavelet.m` — wavelet variant (exists, not fully evaluated)

Scan parameters confirmed: NI=26, NPTS=512, NREPS=32, FOV=250mm, matrix=80, 1 coil, Ntotal=425,984.

---

## Open Issues (in priority order)

1. **⚠️ `im_size` mismatch — must fix before trusting DCF quality**  
   Python notebook (`spiral3d_frames_mat_hoom.ipynb`) hardcodes `im_size=(256,64,256)` — leftover from a cardiac dataset. MATLAB CS script runs at `[90,90,90]` (or whatever `MS_recon_desired` is). The torchkbnufft DCF is computed for the wrong grid. This is not "close enough" — fix `im_size` in the notebook to match the actual recon grid before any serious evaluation.  
   **Status:** flagged, not yet fixed. Decision pending on how to restructure the notebook.

2. **DCF method comparison**  
   Pipeline uses torchkbnufft pipe DCF (Python). MATLAB has iterative KB DCF (`iterative_dcf_fa_20190910`) from the standard gridded recon path. Neither has been compared for CS initialization quality.

3. **CS parameter tuning**  
   `TVWeight=0.01`, `xfmWeight=0.01`, `Itnlim=15` are Lustig demo defaults — not tuned for ASAP spiral on ACR phantom.

4. **Phase 2 (Xe129)**  
   Same MATLAB pipeline; changes: `nchannels` >1, `gamma=11.777 MHz/T`, `kdatas → [Nframes × Ntotal]`, temporal TV (`TVOPDt`).

---

## MATLAB Setup (run at start of every session)

```matlab
% 1. Fessler IRT — NUFFT operators
run('/Users/hoomham/Hooman/Work/Codes/2025_CS/irt/setup.m')

% 2. Lustig sparseMRI — fnlCg, NUFFT3D, TVOP3D, voronoidens
addpath(genpath('/Users/hoomham/Hooman/Work/Codes/2025_CS/sparseMRI_v0.2'))

% 3. mapVBVD is already in Analysis/2025-09-24_ACR/ — no extra step
```

`2025_CS/` contains: `irt/` (Fessler IRT with `setup.m`) and `sparseMRI_v0.2/` (Lustig CS classes + `fnlCg.m` + `utils/voronoidens.m`).  
There is **no mapVBVD in `2025_CS/`** — it lives alongside the analysis scripts.

---

## Session Guide

### "I want to work on the CS code"
Upload: **this page** + [[CS-ASAP-Pipeline]]

### "I want to understand the trajectory / raw data"
Upload: **this page** + [[ASAP-Reconstruction]]

### "I want to plan the Xe129 adaptation"
Upload: **this page** + [[CS-ASAP-Adaptation]] + [[ASAP-Reconstruction]]

### "Full context, fresh start"
Upload: **this page** + [[CS-ASAP-Pipeline]] + [[ASAP-Reconstruction]]

---

## Key Files

| File | Location | Purpose |
|------|----------|---------|
| `cs_spiral_gpt_ACR_hoom_20250924.m` | `2025-09-24_ACR/` | MATLAB: raw → ACR_data.mat |
| `spiral3d_frames_mat_hoom.ipynb` | `2025-09-24_ACR/` | Python: ACR_data → ACR_test.mat |
| `spiral3d_cs_3D_hoom.m` | `2025-09-24_ACR/` | MATLAB CS recon |
| `ACR_data.mat` | `2025-09-24_ACR/` | Intermediate [425984×5] |
| `ACR_test.mat` | `2025-09-24_ACR/` | CS input (ktrajs/kdatas/kcomps) |
| calibration `.mat` | saved to `outdir` by script | Built from XYZ .dat files |

**Raw data location:** `Work/Images/2025/2025-08-16_ACR/RAW/`  
**Calibration scans:** `Work/Images/2025/2025-08-18_ACR/` (X/Y/Z, v2)

---

## Log

### 2026-05-22 — Page created
- Source: Cowork session — first full code audit of `2025-09-24_ACR/`
- Discovered the pipeline was already running (not just planned)
- Traced full chain from .dat → CS output
- Created [[CS-ASAP-Pipeline]] to document the code chain
- Updated [[ASAP-Reconstruction]] and [[CS-ASAP-Adaptation]] with confirmed details
