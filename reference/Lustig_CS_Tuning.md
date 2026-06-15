# Lustig MATLAB CS — Parameters, Conditioning, What to Sweep

Companion to `Lustig_CS_Baseline.md` (that = pipeline + why-it-loses; this =
how it's wired and what knobs exist). Grounded in the actual code under
`workspace/codes/2025_CS/sparseMRI_v0.2/` and the three run scripts in
`workspace/codes/2025-09-24_ACR/`. Solver objective (from `fnlCg.m` header):

```
Phi(x) = || FT * XFM' * x  −  data ||²  +  xfmWeight·|x|₁  +  TVWeight·TV(XFM'·x)
```

`x` lives in the **transform domain**: with `XFM=1`, `x` is the image; with
`XFM=Wavelet`, `x` is the wavelet coefficient vector and `XFM'` is the inverse
transform. fnlCg is nonlinear conjugate gradient + backtracking line search.

---

## The three runs at a glance

| Script | XFM (sparsity) | TVWeight | xfmWeight | imSize | Itnlim × outer | init | extras | saved |
|---|---|---|---|---|---|---|---|---|
| `spiral3d_cs_3D_hoom.m` | `1` (identity → L1 on image) | **0.01** | **0.01** | 100³ | 15 × 15 | `FT'*(data.*w)` | data `/max|data|` | `tv01g01.mat` |
| `spiral3d_cs_3D_hoom_Wavelet.m` | `Wavelet('Daubechies',4,4)` (db4, 4 lvl) | 0.00 | **0.01** | **128³** (2ⁿ) | 200 × 10 | `FT'*(data.*w)` | drops traj pts 3:end−5; `data` NOT normalized; fov=350 | `wlet_ph1.mat` (abandoned — slow/diverging) |
| `spiral3d_cs_3D_hoom_Wavelet_Adop.m` | `1` (identity; wavelet line commented) | 0.00→ | **0.00** | 100³ | 5 × 15 | `FT'*(data.*w)` | drops traj pts 3:end; raised-cosine FOV taper (n=1); TV "continuation" (no-op, starts 0) | `wlet_ph1_adopt.mat` |

Read carefully: **`_Wavelet_Adop` is misnamed** — its `XFM = 1` and both
weights are 0, so fnlCg runs **pure least-squares data consistency** (no
sparsity at all) + a one-time FOV edge window. It's an unregularized iterative
gridded recon, not a wavelet CS. The only script that actually uses the wavelet
transform is `_Wavelet.m`, and that one was abandoned.

---

## How the CS is initialized

Yes — all three start from
```matlab
im_dc = FT' * (data .* w);          % density-compensated adjoint ("zero-fill w/ DCF")
res   = XFM * im_dc;                 % move init into the transform domain
```
`im_dc` is the gridded recon: adjoint NUFFT of the DCF-weighted k-space data.
`res` is the optimization variable (coeffs if XFM=Wavelet, else the image).
**The DCF `w` is used ONLY here** — see conditioning §below.

---

## Operator / parameter reference

### `FT = NUFFT3D(k, 1, ph, 0, imSize, 2)`
Fessler-IRT table NUFFT (`Jd=[6 6 6]`, `Kd=2·Nd` oversampling, `nufft_init`).
Arg-by-arg:

| Arg | Value in scripts | Meaning |
|---|---|---|
| `k` | `[dkx dky dkz]` in **[−0.5,0.5]** (`ktrajs/(2π)`) | nonuniform k-coords; ×2π → `om` inside |
| `w` | **1** | internal DCF — set to 1, i.e. NUFFT itself is **unweighted**. (The DCF used for init is applied separately as `data.*w`.) |
| `phase` (`ph`) | **1** | image phase-correction / spatial mask. `1` = none. Commented code shows it can be a circular-FOV mask `sqrt(x²+y²)<1` — a lever to kill out-of-FOV energy. |
| `shift` | 0 | image-space shift (FOV recentering) |
| `imSize` | `[MS MS MS]` | recon grid |
| `mode` | **2** | 1 = constrain image real, **2 = complex image** |

Adjoint (`FT'`) divides by `sqrt(prod(imSize))` and multiplies by `conj(phase)`.

### `w` (DCF, from `kcomps`)
`w = kcomps(frame,:)'; w = w/max(w(:))`. This is the **torchkbnufft pipe DCF**
imported via `ACR_test.mat`. Only enters `im_dc`. Commented alternative:
`w = voronoidens(k)` (MATLAB Voronoi DCF) — a drop-in swap if you want to test
DCF sensitivity without the Python step.

### `data` normalization
`data = data / max(abs(data))` (present in `tv01g01` & `_Adop`, **absent** in
`_Wavelet`). Scales the data term so objective-unit weights (0.01) land at a
predictable magnitude. Omitting it (as `_Wavelet` does) changes the effective
regularization strength — a hidden confound in that run.

### `XFM` — sparsity transform
- `1`: identity. `|x|₁` penalizes image-pixel magnitude (soft pull to 0).
- `Wavelet('Daubechies',4,4)`: `MakeONFilter('Daubechies',4)` filter,
  coarsest level `L=4`. **Caveat (real bug):** the transform is `FWT2_PO` —
  a strictly **2D, n×n dyadic** Wavelab routine — but it's handed the **3D**
  volume. It does **not** correctly transform the through-slice axis, and it
  **requires power-of-two** size (hence 128³, not 100³). This is almost
  certainly why `_Wavelet` was slow/unstable and got abandoned. A correct 3D
  wavelet CS needs a 3D transform (or an explicit slice loop), not this object.

### `TV = TVOP3D` — total variation
Genuine **3D** forward finite differences along x,y,z (`private/D.m`:
`cat(4,Dx,Dy,Dz)`), with matched adjoint. `TV(·)` in the objective is the
isotropic-ish ℓ1 of these gradients (smoothed by `l1Smooth`). Higher
`TVWeight` → flatter, cartoon/stair-stepped image; edges preserved but texture
and fine bars wash out.

### Weights
- `xfmWeight` — strength of `|x|₁`. With identity XFM, shrinks pixels toward 0
  (sparser foreground, darker background). With wavelet XFM, shrinks
  coefficients (denoise/sharpen-ish if the transform were correct).
- `TVWeight` — strength of TV. Smoothing/piecewise-constant prior.
- Both are in **objective units** on the (max-normalized) data — see §λ note.

### Iteration controls
- `Itnlim` — **inner** NLCG iterations inside one `fnlCg` call.
- outer loop `for n=1:N` — **warm restarts**: each pass re-enters fnlCg from the
  previous `res`. Effective work ≈ `Itnlim × N` (e.g. tv01g01 = 15×15 = 225;
  `_Wavelet` = 200×10 = 2000 → why it was slow).
- `gas(n,...)` stores the per-outer-iteration image (the 15-frame stacks you
  already have).

### `init.m` defaults (rarely touched, but they exist)
| Field | Default | Role |
|---|---|---|
| `l1Smooth` | 1e-15 | ε in `(|w|²+ε)^(p/2)` — smooths the L1 corner so it's differentiable for NLCG |
| `pNorm` | 1 | norm type (1 = L1) |
| `gradToll` | 1e-30 | step-size stop tol (effectively unused) |
| `lineSearchItnlim` | 150 | max backtracking steps |
| `lineSearchAlpha` | 0.01 | Armijo sufficient-decrease constant |
| `lineSearchBeta` | 0.6 | backtracking shrink factor |
| `lineSearchT0` | 1 | initial step (auto-adapts: shrinks if >2 LS iters, grows if <1) |

### `_Adop`-only extras
- **FOV edge taper**: raised-cosine window `edgeWin` (flat to `r0=0.85` of
  radius, cosine roll-off over `bw=0.5` to 0), multiplied into the image at
  `n==1` to suppress boundary brightening / out-of-FOV ringing. A blunt but
  effective spatial mask — the principled version is `ph` (the NUFFT phase/mask
  arg) so it's enforced every iteration, not once.
- **TV continuation**: `param.TVWeight *= 0.25` every 5 outer loops — but
  `TVWeight0=0`, so it multiplies zero → **no-op** as written. The scaffolding
  is there if you set a nonzero start.
- `dropoint` trims leading (and `_Wavelet` also trailing) trajectory samples —
  drops the very center / very edge of each readout (DC glitch, gradient
  ramp/decay points).

---

## Conditioning — why this solver behaves as it does

1. **DCF lives only in the init (`im_dc`), never in the objective.** The NUFFT
   is built with `w=1`, so fnlCg minimizes an **unpreconditioned** data term.
   On this spiral the sample-density spread is ~1.8e6 (see
   `CS_Implementation.md` §4), so NLCG moves the solution very little per
   iteration → **final image ≈ DCF-gridded init + small correction**. The DCF
   choice dominates output quality; the regularizer barely engages. This is the
   central limitation and the reason all three runs look like polished gridded
   recons rather than aggressive CS.
2. **λ in objective units.** `TVWeight`/`xfmWeight` are absolute penalty
   weights, not threshold-relative. They only behave predictably because `data`
   is max-normalized; change the data scale (or drop the normalization, as
   `_Wavelet` does) and the same 0.01 means something different. Our pipeline
   avoids this by parameterizing λ as a coefficient-space threshold (D5).
3. **Wavelet transform is 2D applied to a 3D volume** → through-plane
   incoherence not exploited, and power-of-two forced. Treat `_Wavelet` results
   as unsound until the transform is fixed.
4. **Unregularized LS can pass for CS.** `_Adop` (both weights 0) still yields a
   clean-looking image — plausible ≠ regularized. Always check whether a
   penalty is actually active before crediting it.

---

## What to sweep, in priority order

Assuming the goal is the best Lustig-side image to put next to ours:

1. **Sparsity weight on a *correct* transform.** First fix the wavelet (use a
   true 3D transform or slice-loop FWT2_PO), then sweep `xfmWeight` ∈
   {0.003, 0.01, 0.03, 0.1} on max-normalized data. This is the single biggest
   image-quality lever and the closest match to our db4 wavelet CS.
2. **`TVWeight`** ∈ {0, 0.001, 0.01, 0.05}. With identity XFM this is the only
   real regularizer; tv01g01 used 0.01. Higher = smoother/cartoon, watch the
   ACR resolution bars.
3. **Precondition the data term** (biggest *structural* win): fold `w` into the
   NUFFT (or use ADMM-with-inner-CG) so the objective — not just the init —
   sees the DCF. Without this, knobs 1–2 are fighting a stalled solver. This is
   the same lesson as our CG-vs-FISTA conditioning trap.
4. **`imSize`** — 100³ for TV/identity; **128³ (or any 2ⁿ) is mandatory** for
   the wavelet path.
5. **`Itnlim × outer`** — raise only after preconditioning; otherwise more
   iterations just buy more smoothing toward the init (confirmed: tv01g01
   metrics are monotonic across all 15, no sweet spot).
6. **DCF source** — swap torchkbnufft pipe ↔ MATLAB `voronoidens(k)` to test
   init sensitivity (cheap, one-line).
7. **FOV mask via `ph`** — pass a circular/spherical mask as the NUFFT `phase`
   arg to enforce FOV support every iteration (cleaner than `_Adop`'s one-shot
   `edgeWin`).
8. **`dropoint`** — trim leading/trailing readout samples (DC spike, ramp
   points); small effect, but cheap to test 0 / 3 / 5.

Quickest meaningful experiment: fix the wavelet to a real 3D transform + add
DCF preconditioning, then sweep `xfmWeight`. Without 3 first, expect the same
"looks like a gridded recon" outcome regardless of the weights.
