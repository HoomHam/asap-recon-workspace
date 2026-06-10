# 2023_Faraz_Recon_HH — Hooman's MATLAB Recon Fork

## Purpose
Owns: Hooman's active MATLAB reconstruction — spiral k-space gridding, DCF, coil combine for Xe-129 lung MRI
Does not own: Steve's Python/CUDA pipeline (`../../../../`), compressed sensing (→ `../../helpers/`)

## Code Map

### Find It Fast
| Looking for... | Go to |
|----------------|-------|
| Entry point / study script | `spiral_human_20240227.m` |
| Core gridding (KB kernel + neighbor lookup) | `gridrecon_fa_20230113.m` |
| Dynamic recon orchestrator | `fa_spiral_dyn_recon.m` |
| Neighbor lookup table builder | `grid_lookup_20230113.m` |
| KB kernel definition | `createKBkernel.m` |
| Iterative density compensation | `iterative_dcf_fa_20190910.m` |
| Adaptive coil combination | `combinecoils_fa.m` |
| Raw data loader (Siemens TWIX) | `recon_20210622.m` → calls `mapVBVD.m` |
| 3D trajectory loader | `loadtrajectory3D.m` |
| Colormap (scientific colormaps) | `crameri.m` + `CrameriColourMaps6.0.mat` |
| Lorentzian spectral fit | `Fit4LorentzianCplxPh_Dixon_Con_20230126.m` |
| Spiral gradient design | `calcspiralgrad_20190918.m` |
| Radial trajectory | `calcradialtraj_20220713.m` |

### Key Relationships
```
spiral_human_20240227.m (entry)
    └─► recon_20210622.m       [raw data load + initial parse]
            └─► mapVBVD.m      [Siemens TWIX parser]
    └─► loadtrajectory3D.m     [load calibrated traj from .mat]
    └─► grid_lookup_20230113.m [build/load neighbor lookup]
    └─► gridrecon_fa_20230113.m [gridding: KB kernel + DCF]
            └─► createKBkernel.m
            └─► iterative_dcf_fa_20190910.m
    └─► fa_spiral_dyn_recon.m  [dynamic recon for fancy_v3 protocol]
    └─► combinecoils_fa.m      [adaptive Roemer coil combine → magnitude]
```

## Design Rationale

- **Problem**: Steve's Python recon uses Gaussian kernel + per-cell normalization → intensity shading, no iterative DCF
- **Core insight**: KB kernel + iterative DCF gives textbook-quality PSF and uniform intensity; Faraz's stack implements this in MATLAB
- **Constraints**: MATLAB-only; neighbor lookup cached to `.mat` files for speed (expensive to recompute)

## Key Parameters

| Parameter | Location | Default | Effect |
|-----------|----------|---------|--------|
| `os` (oversampling) | `grid_lookup_20230113.m` | 3 | Higher → fewer aliases, more memory |
| `k` (kernel size) | `grid_lookup_20230113.m` | 5 | KB support width; k=5 → k²=25 neighbors per sample |
| `zfill` | `gridrecon_fa_20230113.m` | 1 | Zero-fill factor for display resolution |
| `nreps` | `fa_spiral_dyn_recon.m` | from acq_params | Number of dynamic frames |
| `spokes` | `spiral_human_20240227.m` | 26 | Interleaves per frame |

## Contracts

- MATLAB filename = function name — never rename a file without updating all callers
- `grid_lookup_20230113.m` saves lookup `.mat` keyed by (traj, imgsize, FOV) — cached file must be deleted if any of these change
- `combinecoils_fa.m` returns **magnitude** output (Rician noise, positive floor) — not complex
- `loadtrajectory3D.m` returns KSpaceCoor in **1/mm** — all downstream expects this unit

## Pitfalls

- **Hardcoded paths**: `spiral_human_20240227.m` has paths to Faraz's Windows machine (`C:\Users\faraz\...`) and Hooman's `/Users/hoomham/Hooman/Images/...` — update before running on new data
- **Stale lookup table**: If trajectory or FOV changes, delete `grid_lookup_20220418.mat` and re-run `grid_lookup_20230113.m` — stale cache causes silent wrong gridding
- **Rename = break**: MATLAB resolves functions by filename; renaming `combinecoils_fa.m` without updating `fa_spiral_dyn_recon.m` → silent call to the old version on path
- **fancy_v3 protocol branch**: `fa_spiral_dyn_recon.m` has separate logic for `fancy_v2`/`fancy_v3` protocols — check `contains(protocol,'fancy_v3')` branch if data looks wrong
- **GA (golden angle) flag**: `acq_params.GA` triggers different interleaf ordering — missing field handled by fallback, but verify `GAperiod` is set for GA acquisitions

## Pre-flight Checks

Before running `spiral_human_20240227.m` on new data:
1. Update `files` array with correct `.dat` path
2. Confirm `spokes` matches the acquisition parameter
3. If trajectory changed: delete `grid_lookup_20220418.mat` (will be regenerated)
4. Confirm MATLAB path includes this folder (all functions must be findable)

## Entry Points

| Task | Start Here |
|------|------------|
| Run reconstruction on new data | Edit `files` in `spiral_human_20240227.m`, then run |
| Debug gridding artifacts | `gridrecon_fa_20230113.m` → check `os`, `k`, DCF convergence |
| Change dynamic frame grouping | `fa_spiral_dyn_recon.m` → `nreps`, `intind` logic |
| Understand KB kernel math | `createKBkernel.m` (6 lines) + Eq(5) in Beatty 2005 |
| Port iterative DCF to Python | `iterative_dcf_fa_20190910.m` — 5 CG iterations |
