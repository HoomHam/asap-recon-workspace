# Lustig MATLAB CS — Hooman's Pre-Project Baseline

**Status:** active reference (2026-06-12). Replaces four 2026-05 wiki pages, now in
`archive/` (`CS-Lustig-ASAP-{Pipeline,Status,Adaptation}_2026-05.md`,
`ASAP-Lustig-Reconstruction_2026-05.md`) — full pipeline detail lives there.
Code now mirrored in `workspace/codes/2025_CS/` (Fessler IRT + sparseMRI_v0.2)
and `workspace/codes/2025_Xe129_CS/` (temporal-TV variants). Working scripts +
data remain in `~/Hooman/Work/Analysis/2025-09-24_ACR/`.

**Role in roadmap:** Step-2 personal-baseline comparison (handoff D-step 2).
This is the CS Hooman ran *before* this project — benchmarks where he was
standing, not an independence test (that's BART).

---

## Pipeline in one screen

```
Siemens .dat ──load_rawdata_20250816──▶ twix raw [512×1×832]
   ──loadtrajectory3D('BuildFromXYZ')──▶ KSpaceCoor [Nseg×3] 1/mm
   ──scale α = kmax_nyq/kmax_meas, aggregate──▶ ACR_data.mat [425984×5]
   ──Python: spiral3d_frames_mat_hoom.ipynb (torchkbnufft)──▶ ACR_test.mat
        ktrajs [1×3×N] rad [−π,π],  kdatas [1×N],  kcomps [1×N] (pipe DCF)
   ──MATLAB: spiral3d_cs_3D_hoom.m──▶ CS volume [90,90,90]
        k = ktrajs/(2π) → [−0.5,0.5];  FT = NUFFT3D(k,1,1,0,imSize,2)
        im_dc = FT'*(data.*w);  res = fnlCg ×(15 outer × Itnlim 15)
        TVOP3D + L1, TVWeight = xfmWeight = 0.01, XFM = 1 (identity, NOT wavelet)
```

Three trajectory unit conversions in a row: 1/mm → [−π,π] (Python, ÷max|k|·π)
→ [−0.5,0.5] (MATLAB, ÷2π) → ×2π inside NUFFT3D for IRT `om`. Net: max|k| → π.

## How densities were generated (and the framework mismatch)

`kcomps` come from a **third framework**: torchkbnufft's Pipe-Menon
`calc_density_compensation_function`, computed in the Python notebook on the
[−π,π] trajectory — not from IRT (which does the recon) and not from the
MATLAB gridded path (Faraz's `iterative_dcf_fa_20190910`, never compared).
Then `w = w/max(w)`.

Two structural facts:
1. **Known bug:** pipe DCF computed with `im_size=(256,64,256)` — cardiac
   leftover; recon grid is [90,90,90]. DCF is for the wrong (anisotropic!)
   grid. Flagged 2026-05-22, never fixed.
2. **DCF enters ONLY the init.** `NUFFT3D(k, 1, …)` is built with w=1; fnlCg
   minimizes unweighted ‖Fx−y‖² + λ·smoothed-L1. `kcomps` touch nothing but
   `im_dc = FT'*(data.*w)`.

Consequence (matches our measured conditioning trap, `CS_Implementation.md`
§4): on the spiral operator with ~1.8e6 density spread, an unpreconditioned
first-order/NCG solver moves the solution very little per iteration. The
final image ≈ init + small correction → **output quality is dominated by the
DCF used in the init**, even though the objective never sees it.

## Why "Steve's densities" gave no good result

Steve's framework has no per-sample DCF to lend. His "density" is `knorm`
(`recon.py: cudarecon`/`cudarenorm`): per-**grid-cell** accumulated Gaussian
kernel weights on the MS³ Cartesian grid, used as post-gridding division
`k[idx] /= knorm[idx]` (Jackson-style weight normalization). Plugging it (or
anything derived from it) in as `kcomps` fails structurally:

| Mismatch | Lustig pipeline expects | Steve's knorm is |
|---|---|---|
| Domain | per-k-sample weight, length 425 984 | per-grid-cell array, MS³ |
| Operation | pre-weight data before adjoint | divide grid after accumulation — not equivalent |
| Kernel/units | DCF for IRT KB NUFFT, Kd=2·Nd, traj [−0.5,0.5] | Gaussian kdist0sq=0.2, traj in grid-index units, no oversampling |
| Edge behavior | smooth pipe weights | knorm→eps floor at empty cells → 1/knorm explodes at k-space edge |
| Sample ordering | rawdata2 order [npts×leaves×reps] | modular `kidx=(idx+idxoff)%nuniquesmp`, one-rep trajectory tiled |

Any one of these breaks the init image; combined with the "final ≈ init"
dynamic above, garbage densities → garbage CS output. So the bad result was
expected, not a tuning failure. (Same reason the wrong-grid pipe DCF needs
fixing before this baseline is trusted at all.)

## vs our CS (`helpers/recon/cs_recon.py`)

| Aspect | Lustig baseline | Ours |
|---|---|---|
| Operator | Fessler IRT table NUFFT (Jd=6, Kd=2Nd) | FINUFFT (validated wrappers) |
| Solver | fnlCg: nonlinear CG, backtracking, smoothed L1 (`l1Smooth`), 15×15 restarts | FISTA (exact L1 prox, wavelet) / PDHG (TV); CG baseline |
| DCF role | init only; data term unweighted | preconditioner inside data term, w = 1/\|AAᴴ1\| — D6 contract |
| λ semantics | objective units on max-normalized data (TVWeight=xfmWeight=0.01, Lustig demo defaults, untuned) | soft-threshold in coefficient units rel. p99 of \|W·x_cg\| — D5 contract |
| Sparsity | TV3D + L1 on **image** (XFM=1; wavelet variant exists, unevaluated) | db4 wavelet L1, or TV |
| DCF source | torchkbnufft pipe, wrong grid (256,64,256) | finufft on correct grid |
| Convergence safeguard | none (plausible-looking ≠ converged) | conditioning trap documented; CG cross-check |

Note the Lustig λ scheme is exactly our measured failure mode #1 (λ in
objective units) — Lustig makes it survivable by normalizing `data` to
max 1, but it stays grid- and data-dependent.

## Before using as comparison baseline (step 2)

1. Fix `im_size` in the notebook → recon grid (or regenerate kcomps with
   finufft pipe on our v3_fov250 trajectory — cleaner: same acq inputs).
2. Feed it the same `recon_io/` acq + trajectory (convert grid-index →
   [−0.5,0.5]) so it sees identical data, then `orient_to_match`, fixed
   metrics, montage + eye.
3. Decide wavelet variant (`spiral3d_cs_3D_hoom_Wavelet.m`) vs identity-L1 —
   identity-L1 is not comparable to our wavelet prior.
