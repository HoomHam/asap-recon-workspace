# Steve vs Faraz — Reconstruction Comparison

> **Revival note (2026-06-10):** Original file (created 2026-06-07 10:40) was accidentally deleted. Reconstructed from ScreenPipe screenshots (deleted after revival). Sections the screenshots cut off are marked. **Superseded** by `reference/Recon_Comparison_StaticGas.md` (2026-06-10), which also corrects errors in this doc (MATLAB combine is real-part, not magnitude; Steve does combine coils via b-map; Faraz does not de-apodize).

**Full deep-dive:** `~/Library/Mobile Documents/iCloud~md~obsidian/Documents/Action/MRI/ASAP Recon/ASAP Faraz vs Steve.md`

## Stack Identity

| Dimension | Steve (ASAP / Python) | Faraz (MATLAB) |
|-----------|----------------------|----------------|
| Language | Python + `numba.cuda` | MATLAB |
| Entry point | `main.py` (tkinter GUI) | `spiral_human_20240227.m` → `recon_20210622.m` |
| Raw data loader | `raw.py` via `mapvbvd` | `mapVBVD.m` |
| Gridding kernel | **Gaussian** (KB path present but commented out) | **Kaiser–Bessel** via `createKBkernel.m` |
| Density compensation | Per-cell k/knorm normalization | **Iterative DCF** (5 iters, `iterative_dcf_fa_20190910.m`) |
| Trajectory source | `.npy` files (normalized, kmax inferred from k[cut off]) | [cut off — calibrations `.mat`, 1/mm] |
| Oversampling | `MS` parameter | `os` + `zfill` in `grid_lookup_20230113.m` |

## 10 Key Differences (summary — see Obsidian for full)

1. **Kernel:** KB (MATLAB) vs Gaussian (Python) → different PSF/apodization/alias behavior
2. **DCF:** Iterative (MATLAB) vs knorm division (Python) → intensity shading, less center shading
3. **Background stats:** Magnitude-Rician (MATLAB) vs complex/zero-mean (Python) — this explains visual difference
4. **De-apodization:** Explicit in MATLAB; knorm-only in Python
5. **Trajectory units:** MATLAB in 1/mm from `.mat`; Python normalizes `.npy`, rescales to MS
6. **OS/ZF:** Both support, but defaults may differ → replica levels and apparent sharpness differ
7. **Neighborhood:** MATLAB k² lookup table (KB support); Python bxsz cube (Gaussian)
8. **Coil combine:** MATLAB adaptive → SNR-optimal magnitude; Python per-coil complex
9. **Orientation:** MATLAB explicit permute+flip; Python grid order (may be mirrored/transposed)
10. **B0 correction:** Neither applies off-resonance correction currently

## Why Images Look Different

*[section body not captured in screenshots]*

## Quality Priority

*[section header not captured; bullets visible:]*

- KB kernel + matched de-apodization
- Iterative DCF (even 3–5 iters)
- Adaptive coil combine
- Matched OS/ZF + orientation

## Speed Priority (quick preview / on-scanner)

Current Python path fine — Gaussian GPU gridding is fast. Take SoS magnitude for MATLAB-like visuals.

## To Align Them

1. Enable KB path in `recon.py` (`bessi0` is already there, commented)
2. Add iterative DCF pass in Python (port from MATLAB or use sigpy)
3. Match OS, ZF, neighborhood size
4. Apply same coil combine + orientation flip
