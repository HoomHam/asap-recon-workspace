# Compressed Sensing on the ASAP Operator — Implementation Notes

*2026-06-12. Companion to `Physics_Notes.md` (picks up where its §12 leaves
off) and `Recon_Comparison_StaticGas.md` (whose estimator-objectives capstone
motivates everything here). Code: `helpers/recon/cs_recon.py`,
`helpers/recon/cs_montage.py`. Data: static ACR phantom,
`data/v3_fov250/recon_io/` (M = 424,320 samples, 100³ image grid).*

This document is educational in the same spirit as `Physics_Notes.md`: every
concept is mapped to a line of code, and every claim that came from a measured
failure says so. Two of the three bugs hit during implementation are more
instructive than the working code — they are kept in full.

---

## 1. Why CS, restated in one paragraph

The comparison phase ended with a diagnosis: Steve's gridder wins SNR (28.7 vs
our CG 19.6) because cell-averaging is a *smoothing estimator* — an implicit
regularizer he never wrote down. Faraz wins flatness (lowfreq-CV 0.093)
because Pipe–Menon's fixed point `w ⊛ C = 1` makes flatness *literally his
objective*. Each pipeline optimizes one criterion by accident of its
algorithm. Compressed sensing writes both down explicitly:

```
min_ρ   ‖W^(1/2) (A ρ − s)‖²   +   λ R(ρ)
        └── data fidelity ──┘       └ prior ┘
```

The data term keeps the recon honest against the measurements (resolution,
geometry); the prior term buys noise suppression the way Steve's smoothing
does — but tunable, and chosen rather than inherited. The bar set by the
handoff: beat Steve's SNR **and** Faraz's flatness *simultaneously*, at equal
or better resolution. First sweep result: L1-wavelet at t = 0.003–0.01 does
(SNR 32.8–47.2, lowfreq-CV 0.084–0.085). §8 qualifies how seriously to take
those numbers.

## 2. The pieces and where they live

| Piece | Symbol | Code |
|---|---|---|
| NUFFT forward (image → samples) | A | `asap_recon.forward` → `cs_recon.FinufftForward` |
| NUFFT adjoint (samples → image) | Aᴴ | `asap_recon.adjoint` → `cs_recon.FinufftAdjoint` |
| DCF preconditioning weights | W | computed inline in `cs_recon.main` (`dens`, `w`) |
| Wavelet transform | Ψ | `sigpy.linop.Wavelet(..., wave_name="db4")` |
| Finite differences | G | `sigpy.linop.FiniteDifference` |
| L1 proximal operator | prox | `sigpy.prox.L1Reg` (complex soft-threshold) |
| FISTA / PDHG | — | `sigpy.app.LinearLeastSquares` (solver auto-selected) |
| Quality metrics | — | `cg_tune.metrics` (SNR, CV, lowfreq-CV, extents) |

The operator pair is *ours* — the same finufft functions validated to 6e-15
adjointness, with the Steve-units (grid-index) trajectory convention and
isign = −1 baked in. Wrapping them as sigpy `Linop`s (15 lines each) was
chosen over sigpy's own NUFFT precisely to keep that single validated
convention; sigpy only contributes the solvers, transforms, and proxes.

A `Linop` is just a shape-annotated function with a `.H`. The only subtlety:
`Linop(oshape, ishape)` takes **output shape first** — `FinufftForward` is
`super().__init__((M,), (n, n, n))`.

## 3. Lesson one: λ is meaningless without the solver's step size

FISTA solves `min_x ½‖Ax − y‖² + λ‖x‖₁` by iterating

```
x ← soft_threshold( x − α Aᴴ(Ax − y),  α·λ )
```

with step size α = 1/maxeig(AᴴA). The quantity that actually touches the
image is **α·λ — the soft-threshold in coefficient units** — not λ itself.
λ only means something relative to the curvature of the data term.

This was learned the measured way, twice in one afternoon:

- **Attempt 1:** λ scaled to `max|Aᴴy|` (= 3.8e5, feels like "the natural
  data scale"). maxeig(AᴴA) = 6.6e9, so the effective threshold was
  α·λ ≈ 1.7e-8 … 5.7e-7 — four to five orders below the wavelet coefficient
  scale (~3e-3). Result: TV output **byte-identical across four decades of
  λ**; wavelet barely moved. A perfect no-op that still produced
  plausible-looking images and metrics.
- **Attempt 2:** threshold pinned to coefficient units, swept at
  {0.1, 0.3, 1, 3}× the p99 coefficient magnitude. Now over-regularized:
  t = 0.1 already shrank the object, t = 1 emptied it entirely (crashing the
  extent metric — empty mask guard added to `cg_tune.metrics` as a result).

The final parameterization (`cs_recon.py`, `T_RELS`): pick the threshold *t*
directly, relative to the p99 coefficient magnitude of the CG-20 solution —
a quantity with physical meaning ("kill coefficients smaller than t× the
typical solution coefficient") — and derive λ from it. After the operator
normalization of §5, λ = t exactly. Working range found: t ∈ [0.003, 0.1],
i.e. the interesting regime sits *well below* the p99 coefficient level,
down near the noise floor of the coefficients. That is Donoho-style
wavelet denoising intuition: threshold at the noise level, not at the
signal level.

> **Rule worth keeping:** never sweep λ in objective-function units. Sweep
> the threshold (or equivalent prox strength) in image/coefficient units and
> back out λ. Objective units change every time the operator, FOV, or sample
> count changes; coefficient units don't.

## 4. Lesson two: CG forgives ill-conditioning, gradient methods don't

The deepest problem of the day, and the one with real physics in it.

The spiral center is enormously oversampled — the measured sample-density
spread on this trajectory is **1.8 × 10⁶ : 1** (`dens.max()/dens.min()` in
`cs_recon.main`). The spectrum of AᴴA inherits that spread: maxeig (6.6e9)
lives in the heavily-sampled low-frequency directions; the bulk of the
spectrum, carrying the actual image detail, sits 3–4 orders lower.

Why CG never cared: CG is a Krylov method — it builds a polynomial in AᴴA
tuned to the actual spectrum and converges per *cluster* of eigenvalues, not
per condition number. 15–20 iterations on this operator was always enough.

Why FISTA/PDHG stalled: a gradient step with α = 1/maxeig advances each
eigendirection by a factor (1 − α·λᵢ) per iteration. Directions with
λᵢ/maxeig ≈ 1e-3 need *thousands* of iterations to move. Measured: after 300
PDHG iterations the recon amplitude was still ~1000× below the CG solution.
The first "working" TV results were λ-independent unconverged garbage that
nevertheless looked like images — the second time in one session that a
silent failure produced plausible pictures (see the eye-over-metrics rule).

**The fix is density compensation, returning in a new role.** Weight the
data term: `min ‖W^(1/2)(Ax − y)‖²` with w = 1/density. Then WAᴴA ≈ flatter
spectrum — the quadrature-weights argument of `Physics_Notes.md` §5, but now
buying *convergence speed* instead of one-shot adjoint accuracy. Density is
estimated with one extra operator pair:

```python
dens = |A(Aᴴ 1)|        # cs_recon.main; Dirichlet-smoothed local density
w = 1 / clip(dens, dens.max()·1e-4, None);  w /= w.mean()
```

This is the first Pipe–Menon iteration with the full PSF as kernel — the
exact construction the deleted DCF module proved *inadequate for direct
recon* (corr 0.38). As a preconditioner the quality bar is far lower, and it
is plenty: weighted FISTA reaches the CG solution scale in ~20 iterations
(corr 0.945 vs unweighted CG-20, 3 s).

**The honest caveat:** with W in the objective, the minimizer is no longer
the unweighted least-squares solution. Scanner noise is i.i.d. across
samples, so *unweighted* LS is the maximum-likelihood estimator; W
de-weights the oversampled center, i.e. throws away some averaging there in
exchange for usable convergence. Every gridding pipeline (Steve, Faraz) makes
the same trade implicitly — DCF-weighted adjoint *is* weighted LS in
disguise. CS just makes it visible. If this ever matters, the alternatives
are: ADMM with inner CG solves (tested — converges, ~2× slower per quality),
or W as a true preconditioner in a primal-dual scheme rather than an
objective change. Parked.

**Architecture dividend:** W is per-sample diagonal weights — exactly the
slot decision D3 reserved for bin weights. Dynamic binned CS later =
`W_total = W_dcf · W_bin` with zero structural change.

## 5. Operator normalization: one scalar kills two bug classes

`cs_recon.main` builds the solver-facing operator as

```python
c = 1/sqrt(maxeig(A_wᴴ A_w))
A_n = Multiply(c·sqrt(w)) ∘ A_raw      # maxeig(A_nᴴA_n) = 1
y_n = c · sqrt(w) · y
```

Scaling A and y by the *same* constant leaves the least-squares minimizer
unchanged (the objective just scales by c²) but pins maxeig to 1. Three
consequences:

1. FISTA's α = 1, so **threshold = λ literally** — the §3 parameterization
   becomes exact instead of derived.
2. PDHG's primal/dual step split becomes sane: before normalization
   ‖A‖ ≈ 8e4 versus ‖G‖ ≈ 3.5, and sigpy's step-size choice left the primal
   step at ~1e-10 (the byte-identical-TV bug had *two* causes; this was the
   second).
3. Every λ, threshold, and tolerance in the code is now in image-coefficient
   units, portable across datasets.

The maxeig itself comes from 30 power iterations (`sigpy.app.MaxEig`), ~10 s,
computed once and reused by `cs_montage.py` via `cs_sweep_metrics.json`.

## 6. The two priors

### L1-wavelet (synthesis form, FISTA)

```python
W = sigpy.linop.Wavelet(img_shape, wave_name="db4")
min_a ½‖(A_n Wᴴ) a − y_n‖² + λ‖a‖₁ ;   x = Wᴴ a       # cs_recon.wavelet_recon
```

Solve for the *coefficients* a, with the operator A_n∘Wᴴ. Because db4 is
orthogonal (WᴴW = I), maxeig(W A_nᴴA_n Wᴴ) = maxeig(A_nᴴA_n) = 1, so the
precomputed normalization carries over — this is why an orthogonal wavelet
was chosen over a redundant frame for the first pass. The prior says: the
image is sparse in a multiscale basis — smooth regions compress into few
coarse coefficients, noise spreads thinly across many fine ones and gets
thresholded away. Failure mode, clearly visible in the sweep sheet at
t ≥ 0.03: **blocky db4 texture** — surviving coefficients reconstruct as
literal wavelet atoms. (Production fix if ever needed: translation-invariant
/ undecimated wavelets, or cycle spinning. Costs the orthogonality shortcut.)

The complex soft-threshold (`sigpy.prox.L1Reg`) shrinks *magnitudes* and
keeps phase: `a ← a/|a| · max(|a| − t, 0)`. MRI images are complex; never
threshold real and imaginary parts independently.

### TV (analysis form, PDHG)

```python
G = sigpy.linop.FiniteDifference(img_shape)
min_x ½‖A_n x − y_n‖² + λ‖Gx‖₁                        # cs_recon.tv_recon
```

`‖Gx‖₁` can't be handled by FISTA (prox of the *composition* has no closed
form), hence PDHG, which dualizes G and only ever needs prox of ‖·‖₁ itself.
The prior says: the image is piecewise constant — true for a phantom, only
roughly true for a lung. Failure mode at large λ: staircasing and blob-ification
(visible at t = 0.1: edges survive but interior texture is gone and the
resolution inserts smear). On this fully-sampled phantom TV never reached the
wavelet's metrics (best TV: lfCV 0.124 at t = 0.003 vs wavelet 0.084); its
real test is the undersampled regime, where TV's edge-preservation
historically shines.

## 7. First sweep results (fully sampled, static, single bin)

Reference: CG-20 baseline SNR 19.5, lowfreq-CV 0.123. Bars: Steve-equiv SNR
28.7, Faraz lowfreq-CV 0.093, ACR true extents 190/190/148 mm. Full numbers:
`cs_sweep_metrics.json`; images: `cs_sweep_sheet.png`,
`montage_cs_t0.003.png`, `montage_cs_t0.01.png`.

| recon | SNR | lowfreq-CV | extents (mm) | eye |
|---|---|---|---|---|
| wavelet t0.003 | 32.8 | **0.084** | 212/175/250 | most natural texture |
| wavelet t0.01 | 47.2 | **0.085** | 180/172/185 | cleaner; block texture onset |
| wavelet t0.03 | 67.9 | 0.087 | 175/175/140 | clearly blocky |
| wavelet t0.1 | 68.1 | 0.092 | 172/175/138 | checkerboard artifacts |
| tv t0.003 | 18.1 | 0.124 | 250/250/250 | ≈ CG, slightly smoothed |
| tv t0.1 | 84.2 | 0.143 | 172/168/128 | blurred to blobs |

Both wavelet candidates beat both bars simultaneously — the handoff thesis
confirmed on its easiest case. λ choice between t0.003 and t0.01 is
**Hooman's visual verdict, pending** — montages are rendered for exactly that.

## 8. How these metrics lie (read before quoting them)

1. **SNR is inflated by background suppression.** Soft-thresholding zeroes
   background coefficients → background σ collapses → SNR = μ_obj/σ_bg
   explodes (t0.1: SNR 68 on a visibly *worse* image). A denoising prior
   partially buys its SNR by cleaning the *measurement* of noise, not the
   signal. Compare SNR only at fixed visual sharpness, and check the
   resolution inserts.
2. **Extents are threshold-fragile.** The 0.25·p99 mask reads 250 mm (whole
   FOV) on CG because background noise clears the threshold; on CS volumes
   the background is clean so extents suddenly mean something — the columns
   are not comparable across regularization strengths. The wavelet z-extent
   shrinking (185 → 140 with t) is real L1 edge erosion though: dim end
   slices die first.
3. **Per-volume normalization in the montages** (p99.5 → [0,1], shared color
   axis — `cs_montage.py`) makes structure and background comparable across
   pipelines but says nothing about absolute intensity; the three volumes are
   in unrelated units at the source.
4. Twice this session, broken recons produced plausible images and metrics
   (λ no-op; unconverged PDHG). The standing rule held: **no image-quality
   conclusion without eyes on the figure.**

## 9. What's next (handoff CS spec, remaining steps)

1. **λ verdict** — Hooman, from the montages. Then freeze t for the wavelet.
2. **Undersampling experiment** (the point of CS): retrospectively keep 1/2,
   1/4, 1/8 of interleaves; recon CS vs CG vs Steve-equiv at each factor.
   Prediction from the capstone: Steve's gridder degrades first (cell
   averages break as cells empty). Verifying that is publishable material.
   Note: undersampling *changes the density landscape* — recompute `dens`
   per mask, never reuse the fully-sampled W.
3. **Binned/dynamic CS** — W_bin into the same weight slot, temporal
   regularizer across bins. Needs dynamic data or simulated bin weights.
4. Parked: ADMM-with-inner-CG as the statistically cleaner solver; undecimated
   wavelets if block texture survives at the chosen t.

## 10. Quick reproduction

```bash
cd workspace/helpers/recon
../.venv/bin/python cs_recon.py   ../../data/v3_fov250/recon_io --max-iter 100
../.venv/bin/python cs_montage.py ../../data/v3_fov250/recon_io --t-rel 0.003
```

~16 s per wavelet recon, ~22 s per TV recon, M4 Pro CPU. `selftest.py` is
untouched by all of this (operator module unchanged) but run it anyway after
any future operator edit.
