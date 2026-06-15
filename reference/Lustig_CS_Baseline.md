# Lustig MATLAB CS — Hooman's Pre-Project Baseline

**Status:** active reference (updated 2026-06-15). Replaces four 2026-05 wiki
pages, now in `archive/` (`CS-Lustig-ASAP-{Pipeline,Status,Adaptation}_2026-05.md`,
`ASAP-Lustig-Reconstruction_2026-05.md`) — full pipeline detail lives there.
Code mirrored in `workspace/codes/2025_CS/` (Fessler IRT + sparseMRI_v0.2) and
`workspace/codes/2025_Xe129_CS/` (temporal-TV variants); a manually-run copy
with three saved runs is in `workspace/codes/2025-09-24_ACR/`. Working scripts +
data also in `~/Hooman/Work/Analysis/2025-09-24_ACR/`.

**Now runnable in one command** — `workspace/helpers/lustig_oneshot/` collapses
the old MATLAB→Colab→Drive→MATLAB dance into `run_lustig.py <recon_io>`. See
the one-shot section below.

**Role in roadmap:** Step-2 personal-baseline comparison (handoff D-step 2).
This is the CS Hooman ran *before* this project — benchmarks where he was
standing, not an independence test (that's BART). **Verdict (2026-06-15): does
not beat our CS** — see results below.

---

## Pipeline in one screen

```
Siemens .dat ──load_rawdata_20250816──▶ twix raw [512×1×832]
   ──loadtrajectory3D('BuildFromXYZ')──▶ KSpaceCoor [Nseg×3] 1/mm
   ──scale α = kmax_nyq/kmax_meas, aggregate──▶ ACR_data.mat [425984×5]
   ──Python: spiral3d_frames_mat_hoom.ipynb (torchkbnufft)──▶ ACR_test.mat
        ktrajs [1×3×N] rad [−π,π],  kdatas [1×N],  kcomps [1×N] (pipe DCF)
   ──MATLAB: spiral3d_cs_3D_hoom.m──▶ CS volume [100,100,100]  (was 90 in 2026-05)
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
1. **`im_size` — bug was in an OLDER notebook only.** The 2026-05 docs flagged
   `im_size=(256,64,256)` (cardiac leftover, anisotropic). The notebook version
   actually used for the saved runs / one-shot sets `im_size=(100,100,100)` —
   correct, isotropic, matches the recon grid. So there is **nothing to fix**
   here; the earlier warning applied to a superseded notebook.
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
expected, not a tuning failure.

## vs our CS (`helpers/recon/cs_recon.py`)

| Aspect | Lustig baseline | Ours |
|---|---|---|
| Operator | Fessler IRT table NUFFT (Jd=6, Kd=2Nd) | FINUFFT (validated wrappers) |
| Solver | fnlCg: nonlinear CG, backtracking, smoothed L1 (`l1Smooth`), 15×15 restarts | FISTA (exact L1 prox, wavelet) / PDHG (TV); CG baseline |
| DCF role | init only; data term unweighted | preconditioner inside data term, w = 1/\|AAᴴ1\| — D6 contract |
| λ semantics | objective units on max-normalized data (TVWeight=xfmWeight=0.01, Lustig demo defaults, untuned) | soft-threshold in coefficient units rel. p99 of \|W·x_cg\| — D5 contract |
| Sparsity | TV3D + L1 on **image** (XFM=1; wavelet variant exists, unevaluated) | db4 wavelet L1, or TV |
| DCF source | torchkbnufft pipe, `im_size=(100,100,100)` (correct) | finufft `w=1/\|AAᴴ1\|` |
| Convergence safeguard | none (plausible-looking ≠ converged) | conditioning trap documented; CG cross-check |

Note the Lustig λ scheme is exactly our measured failure mode #1 (λ in
objective units) — Lustig makes it survivable by normalizing `data` to
max 1, but it stays grid- and data-dependent.

## One-shot tool (`workspace/helpers/lustig_oneshot/`)

Single command, reads a `recon_io` folder directly, reproduces the old analysis:

```bash
.venv_lustig/bin/python run_lustig.py ../../data/v3_fov250/recon_io
# -> <recon_io>/lustig/ : ACR_test.mat, lustig_cs.mat  gas(15,100,100,100)
```

Three pieces (README has full detail):
- `build_acrtest.py` — npy → ACR_test.mat. Exact replica of
  `spiral3d_frames_mat_hoom.ipynb`: recenter grid-index traj (−MS/2),
  max-radius normalize to [−π,π], torchkbnufft pipe DCF at `(100,100,100)`.
- `run_cs.m` — headless byte-for-byte `spiral3d_cs_3D_hoom.m` (NUFFT3D + fnlCg,
  TV+L1 0.01/0.01, 15 iters, per-slice rot90), auto-adds IRT + sparseMRI paths.
- `run_lustig.py` — driver (build → `matlab -batch` → metrics).

Exactness hinges on two things: trajectory **center offset** (recon_io's
max-radius normalization is scale-invariant, only −MS/2 matters), and the DCF
staying **torchkbnufft pipe** (MATLAB `voronoidens` is a different algorithm) —
hence the dedicated `.venv_lustig` (torch+torchkbnufft, arm64 CPU), separate
from the recon `.venv`. Runtime ≈ 32 s DCF + ~4 min CS.

## Results — same-data comparison on v3_fov250 (2026-06-15)

Ran the one-shot on v3_fov250 (424 320 samples), compared to our CS on the
**same** data. Montage: `data/v3_fov250/recon_io/montage_lustig_v3_vs_ours.png`
(rebuild with `helpers/recon/lustig_compare.py`).

| | SNR | lowfreq-CV | extent mm |
|---|---|---|---|
| Lustig tv01 final iter (15) | 99.5 | **0.115** | 162/160/128 |
| ours wavelet t0.003 | 32.8 | **0.084** | 212/175/250 |
| ours wavelet t0.01 | 47.0 | **0.085** | 180/172/185 |

**Our CS wins.** Lustig is visibly softer (bars/edges smeared); its high SNR is
the background-zeroing artifact (sparsity collapses σ_bg — same blind spot as
`helpers/recon/AGENTS.md`), while the honest metric (lowfreq-CV) is worse than
ours *and* worse than the Faraz bar (0.093). Across all 15 fnlCg iters the trend
is monotonic — no early sweet spot, it just smooths toward the init. Matches the
"final ≈ init" prediction: the stalled NCG never engages real sparsity.

Earlier caveat resolved: the first compare used `tv01g01.mat`, which was the
**old 2025-08-16 ACR scan** (425 984 samples), not v3. The one-shot fixes that —
both pipelines now on identical v3 data, conclusion unchanged.

### If pushing the Lustig side before closing it
Lever is the **prior**, not iterations or data: its wavelet variant
(`spiral3d_cs_3D_hoom_Wavelet.m` / `_Adop`) is closer to our db4. The
identity-L1 TV run here is not a prior-controlled match to ours. The `wlet_ph1`
run was abandoned (diverging, slow); `wlet_ph1_adopt.mat` is the candidate to
re-check via `lustig_compare.py`.
