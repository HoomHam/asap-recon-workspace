# Steve's ASAP Recon — Code Map

**Date:** 2026-06-10 (rewritten fresh from code read; original 2026-06-07 version deleted, partial copy in `archive/screenshots/`)
**Repo:** `/Users/hoomham/Hooman/Work/Codes/2026_ASAP_Recon/`
**Full pipeline notes:** Obsidian `Action/MRI/ASAP Recon/ASAP Recon Pipeline.md`
**Faraz function map:** Obsidian `Action/MRI/ASAP Recon/ASAP Recon Faraz Approach.md`
**Steve vs Faraz diff:** `reference/Recon_Comparison_StaticGas.md`

## File Roles

| File | Role |
|------|------|
| `main.py` | tkinter GUI, global state (`gvar`, `raw`, `traj`, `results` instances), button callbacks. Three recon flows: `calcLVcb` (lung-volume tracking via undersampled recons), `dyn_recon_function` (binned recon), `calc_T1RF_function` (T1/RF decay map from breath-hold ref) |
| `raw.py` | `raw` class: loads Siemens TWIX via `mapvbvd`, gas/dissolved split, noise normalization, spike filtering, binning signals. `traj` class: loads `.npy` trajectory, infers `npts` from \|k\|² periodicity, rescales to grid units |
| `recon.py` | All CUDA kernels via `numba.cuda`: `cudarecon` (gridding), `cudarenorm` (per-cell normalize), `cudarezero`, `bessi0` (unused — KB path commented out) |
| `results.py` | `results` class: `calcb` (coil sensitivity/phase map), `dyn_recon` (binned), `dyn_usimg_recon` (undersampled dynamic), `T1RF_recon`; FFT, crop, coil combine |
| `gtypes.py` | Enums + containers: `gvar` (MS=240, IS=100, nbins=16, gplb=300, dplb=40), `gdir`, `imgtype` (GPDYN/DPDYN/GPREF/DPREF), `bintype` (PNEUMOTACH/DIAPHRAGM/SIGNAL), `graddir`, `species` |
| `asap/asap.c` | C-level ASAP kernel — Kento's cloud GPU deployment target (differs from `workspace/codes/kasap.c`, compare before assuming equivalence) |

## Data Flow (static gas recon path)

```
Siemens .dat ──► raw.load()        mapvbvd parse → [npts, nch, nilvs] complex64
                  │                 kill first 2 pts/ilv (killpts, via trajectory)
                  │                 per-channel noise norm (Gaussian fit to Re histogram)
                  │                 split gas/dissolved by FFT of first-sample pattern
                  │                 global rephase: FID start → zero phase
                  │                 spike filter: >10× cross-ilv mean → zeroed
                  │                 exclude low-SNR fully-sampled ranges
                  │
.npy traj ─────► traj.load()       × FOV → delta-k units; npts from |k|² periodicity
                  │
                 traj.rescale_to_MS()   k·MS/IS + MS/2 → grid-index units (OS = 2.4)
                  │
                 results.calcb()    per channel: full-data gridding → b map
                  │                 (phase-fit polynomial fill in low-signal voxels)
                  │
                 cudarezero ──► cudarecon ──► cudarenorm     (per bin, per channel)
                  │    k=0, knorm=eps   Σ wt·data, Σ wt        k /= knorm
                  │                 kernel: Gaussian exp(-d²/0.2), 4³ cell box
                  │                 optional readout filter exp(-(t/gplb)²)
                  │                 soft bin weight exp(-Δbin²/2) when nbins>1
                  │
                 np.fft.fftshift(fftn(ifftshift(k)))     240³ forward FFT
                  │
                 crop center IS³ (100³) ──► × b_ch, take real ──► Σ channels
                  │
                 result: real image, zero-mean noise (background looks negative;
                 abs() for MATLAB-like visuals)
```

## Recon Entry Points (main.py)

| Button/flow | Calls | Purpose |
|-------------|-------|---------|
| `calcLVcb` | `calcb` → loop `dyn_usimg_recon(iusimg)` | Per-undersampled-image recon → diaphragm/lung-volume tracking |
| `dyn_recon_cb` | `calcb` → `dyn_recon(bins)` | Binned dynamic recon (16 soft bins), gas + dissolved |
| `calc_T1RF_cb` | `calcb` → `T1RF_recon` ×2 halves of breath-hold ref | T1/RF decay map, saved to `T1RF.mat` |

Static single-bin gas recon = `dyn_recon` path with `nbins=1` (binwt=1, soft binning inert).

## GPU Dispatch

`@conditional_decorator` wraps CUDA kernels — `usegpu=False` runs same code on CPU (slow fallback; some paths just `barf`). Kernel launch: `threadsperblock=256` (calcb/dyn_recon) or 16 (usimg/T1RF). Grid accumulation via `cuda.atomic.add` on `k.real`, `k.imag`, `knorm` — contention worst at spiral center. Per-pass: `cudarezero` + full 240³ device-to-host copy. Cloud deployment: `asap/asap.c` is Kento's target — C kernel replacing the Python CUDA path.

## What's NOT in Steve's Code (vs Faraz)

- No iterative DCF (per-cell knorm normalization instead)
- No explicit de-apodization after gridding (Faraz lacks it too — see comparison doc)
- No adaptive (Roemer/Bydder-style published) coil combination — uses own b-map real-part combine (similar in spirit; see comparison doc §9 — old claim "returns per-coil only" was wrong)
- No orientation flip/permute to display convention (dynamic path flips axis 0 only)
- KB kernel commented out (`bessi0` present, unused)
- No grid/DCF caching — everything recomputed per pass
- No pre-averaging of repeated interleaves — full dataset gridded every pass (main compute cost; see comparison doc "Why Steve's recon is computationally heavy")

## Key Magic Numbers

| Value | Where | Meaning |
|-------|-------|---------|
| `MS=240, IS=100` | gtypes.py | Grid / image size → OS = 2.4 |
| `kdist0sq=0.2` | recon.py | Gaussian kernel width (σ ≈ 0.32 cells) |
| `bxsz=2` | recon.py | Gridding box half-size (asymmetric: cx−2…cx+1) |
| `bindist0sq=2` | recon.py | Soft-bin Gaussian width |
| `gplb=300, dplb=40` | gtypes.py | Readout filter lengths (gas / dissolved) |
| `killpts=2` | raw.py | Dropped samples per interleave start |
| `noisespikethresh=10` | raw.py | Spike rejection threshold |
| `knorm<1e-5→0` | recon.py | Empty-cell kluge (k-space holes, not interpolated) |
| TR/extratime kluges | raw.py:350-358 | Hardcoded 22.26 ms TR, 2.191 ms extratime — header rounding workarounds |
