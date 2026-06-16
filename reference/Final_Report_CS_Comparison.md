# Final Report — CS Reconstruction Comparison (ours vs BART vs Lustig vs Steve vs Faraz)

**Date:** 2026-06-15 · **Data:** `data/v3_fov250/recon_io` (ACR phantom, static gas,
single bin, 1 coil, 424 320 spiral samples) · **Status:** comparison phase CLOSED;
this is the final report before the next step (temporal 4D / undersampling).

> **For the presentation agent:** this document is self-contained. Every script
> (purpose / inputs / outputs / run command) is in §7; every figure (what it shows
> / the slide takeaway / FINAL-vs-SUPERSEDED status) is in §8. Read §6 (findings)
> and §8 for slide content; use §5 for the narrative arc. **Do NOT present figures
> tagged SUPERSEDED as conclusions — they encode intermediate, retracted reads.**

---

## 1. TL;DR (one slide)

We benchmarked our compressed-sensing (CS) reconstruction against an independent
implementation (**BART `pics`**, Uecker/Lustig C code) and our personal MATLAB
baseline (**Lustig SparseMRI**), on identical ACR-phantom spiral data, and placed
all pipelines (plus **Steve** and **Faraz** gridding recons) on one ruler.

**Headline result:** after correcting a slice-matching error, **all pipelines
resolve the real ACR resolution-insert grid comparably; our CS is as sharp as
BART and sharper than Lustig.** The independence test (BART) confirms our recon's
geometry and does not beat it. Two scalar-metric traps and one slice-alignment
trap were caught **by eye** mid-analysis and corrected.

---

## 2. Project context

ASAP = spiral hyperpolarized Xe-129 lung MRI recon. Goal chain: understand Steve's
implementation → compare with independent implementations on the same raw data →
build CS → temporal 4D. This session executed the **CS-vs-independent comparison**
(roadmap step 2). All work is in `workspace/` (personal); upstream code is
read-only. Benign phantom signal-processing — no human/PHI data.

---

## 3. Pipelines compared

| Pipeline | Code | Operator | Solver | Prior | Lineage |
|----------|------|----------|--------|-------|---------|
| **ours** | `helpers/recon/cs_recon.py` | FINUFFT (our wrappers) | FISTA (wavelet) / PDHG (TV) / CG | db4 ℓ1-wavelet, TV | sigpy (Lustig-lab) |
| **BART** | `bart pics` (built from source) | BART nuFFT | ADMM / FISTA | ℓ1-wavelet, TV | Uecker lab — **independence test** |
| **Lustig** | `codes/2025_CS/sparseMRI_v0.2` (MATLAB) | Fessler IRT NUFFT | fnlCg (nonlinear-CG, smoothed-ℓ1) | TV+ℓ1 on image (2D wavelet only) | Lustig — **personal baseline** |
| **Steve** | `helpers/recon/steve_kernel_numpy.py` | Gaussian gridding kernel | one-shot gridding | none (filtered) | original author |
| **Faraz** | `codes/.../faraz_recon.mat` (MATLAB) | KB gridding + iterative DCF | gridding | none | parallel implementer |

**Why BART is the independence test:** our CS is built on sigpy (Frank Ong, Lustig
lab); the MATLAB Lustig comparison shares that lineage (personal baseline, partly
circular). BART is a separate codebase → the real independence check.

---

## 4. The objective (all CS pipelines minimize the same thing)

min_x ‖W^½(Ax − y)‖² + λ R(x),  R ∈ {ℓ1-wavelet, TV}

- **A** = non-Cartesian NUFFT operator; **y** = spiral k-space; **W** = DCF weights.
- **ours:** exact ℓ1 prox (FISTA/PDHG), DCF as data-term preconditioner (decision D6).
- **Lustig (`fnlCg.m`, the `phi` = `Phi(x)`):** same objective but *smoothed* ℓ1
  (`(|·|²+l1Smooth)^{p/2}`) + nonlinear-CG; DCF only seeds the init; λ₁=`xfmWeight`,
  λ₂=`TVWeight` (Hooman's old run = 0.01/0.01, both terms, identity transform).
- **BART:** ADMM/FISTA, λ via `-R W:7:0:λ` / `-R T:7:0:λ`.

---

## 5. Narrative arc (for the story slides)

1. **Fixed the metrics first.** The existing `cg_tune.metrics` had two measured
   blind spots (inflated SNR when a prior zeroes the background; full-FOV extent
   from a noise-sensitive threshold). Wrote a NEW `metrics_v2.py` (corner-ROI noise,
   half-max extent, `edge_sharp` guard) — did not touch the old metric, so prior
   results stay reproducible.
2. **Built BART from source** (no homebrew formula; recipe in §9) and ran the first
   ours-vs-BART sweep.
3. **Fixed a BART divergence:** BART wavelet diverged to ~1e32 under default FISTA;
   forced ADMM (`-m`) → stable.
4. **Swept the Lustig baseline** (it had been one fixed point): pure-TV `TVWeight`
   sweep, scored with `metrics_v2`, so it matched BART's λ sweep.
5. **Assembled TV three-way + wavelet two-way** on one ruler.
6. **The texture saga (caught by eye, twice):** a cross-hatch on BART recons looked
   like an artifact. A scalar (HF-band energy) first *inflated* it (a code bug), then
   *dismissed* it (band conflates mottle with grid). **Then Hooman identified it as
   the REAL ACR resolution-insert grid.**
7. **The slice-matching revelation (Hooman):** the phantoms are slightly DIFFERENT
   SIZES across pipelines, so a flip-only orientation aligned the bulk but left the
   thin insert at different z. Fixed-z comparison was invalid — it faked "BART
   resolves, ours blurs." Added **z-affine registration** (scale+shift).
8. **Final, slice-matched result:** all pipelines resolve the grid; ours is sharp.
   Built the **six-way registered montage** (Faraz/Steve/our CG/our wav/BART/Lustig).

---

## 6. Key findings (the meat — one slide each)

**F1 — Independence test passes.** BART independently reproduces the phantom
geometry (orientation corr ≈ 0.92; matched extent). Two unrelated codebases agree.

**F2 — At matched slices, all pipelines resolve the real ACR resolution grid;
ours is as sharp as BART, sharper than Lustig.** *(Figure: `grid_slice_zoom_z52.png`,
`zreg_sixway_montage.png` — FINAL.)*

**F3 — Phantom size differs up to ~23% across pipelines** (z-scale rel. our CG:
Faraz 1.16, BART 1.04, Steve 1.01, our wav 0.98, Lustig 0.94). Real, pipeline-
dependent FOV/scaling difference. **Fixed-z slice comparison is invalid without
z-affine registration.** *(Figure: `zreg_sixway_montage.png`.)*

**F4 — TV is commoditized.** ours/BART/Lustig TV all land ~0.18 lfCV; differences
are small. The CS prior that matters is the wavelet, where Lustig can't compete
(2D-only operator).

**F5 — Scalar image metrics are unreliable here; the eye is the instrument.**
SNR (background-zeroing), lfCV (rewards blur), edge_sharp (rewards texture), and
HF-band energy (conflates grid with mottle) each misled at least once. Every
correct turn this session came from a figure, not a number. *(3rd documented
instance — memory `eye_vs_metric`.)*

**F6 — Lustig's nonlinear-CG stalls** (`obj` ≈ 1094 flat across all 15 iters):
final image ≈ DCF-seeded init. Consistent with the spiral conditioning trap.

---

## 7. Code catalog (everything written this session)

All run with `workspace/helpers/.venv/bin/python` (arm64, finufft 2.5.1, sigpy
0.1.27) unless noted. `<recon_io>` = `data/v3_fov250/recon_io`. Originals untouched.

| Script | Purpose | Inputs | Outputs | Run |
|--------|---------|--------|---------|-----|
| `recon/metrics_v2.py` | Fixed quality metrics (corner-ROI SNR, half-max extent, edge_sharp). Drop-in for `cg_tune.metrics`. | a recon volume (ndarray); CLI: `<recon_io>` (demos old-vs-new on CG-20) | dict{snr,cv,lowfreq_cv,extent_mm,edge_sharp,noise,signal,bg_collapsed} | `metrics_v2.py <recon_io>` |
| `recon/bartio.py` | BART `.cfl/.hdr` read/write (complex64, column-major). No toolbox dep. | `writecfl(name,arr)`; `readcfl(name)` | `.cfl`+`.hdr` files / ndarray | imported |
| `recon/bart_compare.py` | Build BART inputs from recon_io; sweep `pics` (wavelet+TV × λ); score ours + BART with metrics_v2. | `<recon_io>`, `--bart`, `--iters`, `--lambdas` | `bart/{traj,ksp,sens,bart_W_l*,bart_T_l*}.cfl/hdr`, `bart/bart_vs_ours_montage.png`, `bart/bart_compare_metrics.json` | `bart_compare.py <recon_io> --bart <bart>` |
| `recon/tv_threeway.py` | ours/BART/Lustig **TV** on one ruler (metrics_v2). | `<recon_io>`, `--t-rel`, `--iters` | `tv_threeway_montage.png`, `tv_threeway_metrics.json` | `tv_threeway.py <recon_io>` |
| `recon/wavelet_twoway.py` | ours vs BART **wavelet** montage + the texture proof figure. | `<recon_io>`, `--bart`, `--t-rel` | `wavelet_twoway_montage.png`, `wavelet_twoway_metrics.json`, `texture_compare.png` | `wavelet_twoway.py <recon_io> --bart <bart>` |
| `recon/slice_matched_compare.py` | Multi-slice montage, orient (flip-only) to ours. **Superseded by z-affine.** | `<recon_io>`, `--our-t`, `--bart-l`, `--lustig-tv`, `--slices` | `slice_matched_montage.png` | `slice_matched_compare.py <recon_io>` |
| `recon/resolution_sweep.py` | Per-pipeline λ sweep at one z slice (find best by eye). **z NOT affine-matched.** | `<recon_io>`, `--bart`, `--z` | `resolution_sweep_z<z>.png` | `resolution_sweep.py <recon_io> --z 50` |
| `recon/z_register_compare.py` | **z-affine (scale+shift) registration** + fine-z montage (3 pipelines). The methodological fix. | `<recon_io>`, `--our-t`, `--bart-l`, `--lustig-tv`, `--z0`, `--z1` | `zreg_finez_montage.png`; prints fitted scale/shift | `z_register_compare.py <recon_io> --bart <bart>` |
| `recon/zreg_sixway_montage.py` | **FINAL figure:** 6 pipelines (Faraz/Steve/our CG/our wav/BART/Lustig), all z-affine registered to our CG. | `<recon_io>`, `--bart`, `--z0`, `--z1` | `zreg_sixway_montage.png`; prints per-pipeline z-scale | `zreg_sixway_montage.py <recon_io> --bart <bart>` |
| `lustig_oneshot/run_cs_sweep.m` | MATLAB: fnlCg pure-TV sweep over `TVWeight` (xfmWeight=0). | `(in_mat, out_mat, codes_root, tvweights[])` | `lustig_tv_sweep.mat` (`gas_final[numW,100³]`, `tvW`) | via driver |
| `lustig_oneshot/run_lustig_sweep.py` | Driver: reuse ACR_test.mat → MATLAB sweep → score with metrics_v2. | `<recon_io>`, `--tvweights`, `--matlab` | `lustig/lustig_tv_sweep.mat`, `lustig/lustig_tv_sweep_metrics.json` | `.venv_lustig/bin/python run_lustig_sweep.py <recon_io>` |

**Helper functions worth naming for slides:** `orient_to_ours` (flip+1 transpose,
BART), `best_orient_full` / `orient_any` (full 48-orientation, Lustig/Faraz),
`fit_z_affine` + `resample_z` (z scale+shift registration), `grid_hf` (high-freq
band energy — **known unreliable**, see F5).

---

## 8. Figure catalog (every figure, with status)

`<recon_io>` = `data/v3_fov250/recon_io`. **Status: FINAL** = use in presentation;
**PROCESS** = shows the method/journey, fine to show as "how we got there";
**SUPERSEDED** = encodes a retracted conclusion, do NOT present as a result.

| Figure | Shows | Slide takeaway | Status |
|--------|-------|----------------|--------|
| `zreg_sixway_montage.png` | 6 pipelines, z-affine registered, fine-z columns | All resolve the ACR grid at matched z; ours sharp; sizes differ ~23% | **FINAL** |
| `grid_slice_zoom_z52.png` | ours/BART/Lustig zoomed at matched grid slice | ours & Lustig sharp; BART softer here — ours not blurring | **FINAL** |
| `zreg_finez_montage.png` | ours/BART/Lustig fine-z AFTER z-affine | grid lines up in same column once registered | **FINAL** |
| `bart/bart_compare_metrics.json` + `bart/bart_vs_ours_montage.png` | first ours-vs-BART λ sweep (wavelet ADMM-fixed, TV) | independence test; metric numbers | PROCESS |
| `lustig/lustig_tv_sweep_metrics.json` | Lustig TV sweep, metrics_v2 | Lustig TV ~0.18 lfCV, NCG stalls | PROCESS |
| `tv_threeway_montage.png` | ours/BART/Lustig TV at fixed z=50 | TV wash — BUT fixed-z, not size-matched | SUPERSEDED (use registered) |
| `wavelet_twoway_montage.png` | ours vs BART wavelet at fixed z=50 | looked like "ours clean, BART grid" — that was the z-mismatch | SUPERSEDED |
| `texture_compare.png` | ours±DCF vs BART, texture | "BART grid" — later shown to be real ACR structure at a mismatched slice | SUPERSEDED |
| `slice_matched_compare`→`slice_matched_montage.png` | 5 slices, flip-only orientation | structures roughly line up, but no z-scale | SUPERSEDED (flip-only) |
| `resolution_sweep_z50.png` | per-pipeline λ sweep at fixed z=50 | useful λ-vs-detail trend, but z=50 not size-matched | PROCESS (caveat) |
| `finez_montage.png` | fine-z, flip-only orientation | first proof grid is in all three (~z52) | PROCESS |

**Prior-session figures still on disk** (Steve-vs-Faraz / earlier CS): `montage_all.png`,
`montage_{steve,cg,faraz}.png`, `montage_lustig_v3_vs_ours.png`, `cs_sweep_sheet.png`,
`cg_tune_sheet.png`, `baseline_comparison.png` — context, not this comparison's conclusions.

---

## 9. Reproduce from scratch

**Environment:** `workspace/helpers/.venv` (arm64; finufft 2.5.1, sigpy 0.1.27).
Lustig DCF step needs `helpers/lustig_oneshot/.venv_lustig` (torch+torchkbnufft).
MATLAB R2025a on PATH.

**BART build (macOS, from source — no working homebrew formula):**
```bash
brew install fftw openblas gcc make
git clone --depth 1 https://github.com/mrirecon/bart ~/bin/bart-src
ln -sf "$(brew --prefix openblas)/include/cblas.h" \
       "$(brew --prefix openblas)/include/cblas_openblas.h"   # BART expects this name
cd ~/bin/bart-src && gmake clean
CPATH="$(brew --prefix fftw)/include" LIBRARY_PATH="$(brew --prefix fftw)/lib" \
  gmake -j4 CC=gcc-15 FFTW_BASE=$(brew --prefix fftw) \
            BLAS_BASE=$(brew --prefix openblas) OPENBLAS=1 PNG=0
# binary -> ~/bin/bart-src/bart   (mac `make` is 3.81 -> must use gmake; mac `gcc` is clang -> gcc-15)
```

**Run order (produces all FINAL artifacts):**
```bash
cd workspace/helpers/recon
../.venv/bin/python bart_compare.py ../../data/v3_fov250/recon_io --bart ~/bin/bart-src/bart
../lustig_oneshot/.venv_lustig/bin/python ../lustig_oneshot/run_lustig_sweep.py \
    ../../data/v3_fov250/recon_io
../.venv/bin/python z_register_compare.py ../../data/v3_fov250/recon_io --bart ~/bin/bart-src/bart
../.venv/bin/python zreg_sixway_montage.py ../../data/v3_fov250/recon_io --bart ~/bin/bart-src/bart
```

**BART conventions that matter:** trajectory `bart_traj = traj_rad·N/(2π)`, dims
`(3,M,1)`; ksp `(1,M,1)`; single-coil `sens = ones(N,N,N,1)`; **wavelet MUST use
ADMM (`-m`)** (FISTA diverges).

---

## 10. Open threads / next steps

1. **Why phantom sizes differ ~23%** across pipelines (FOV / oversampling
   convention per implementation). Must pin before any quantitative
   extent/resolution comparison. *(Highest priority — affects everything quant.)*
2. **A resolution-targeted metric** (line-pair contrast at the insert slice) so a
   number finally tracks the eye instead of fighting it (F5).
3. **BART per-sample DCF** — `-p` did not wire noncart density in this build;
   would make BART's weighting match ours.
4. **Temporal 4D / undersampling** (the actual next project phase): diaphragm
   binning → per-bin weights; undersampling returns here.

---

## 11. Companion docs
- `reference/BART_Comparison.md` — detailed chronology, the texture/slice saga with
  all the dead-ends (DCF/Toeplitz/kernel/fftmod ruled out).
- `reference/Lustig_CS_Baseline.md`, `reference/Lustig_CS_Tuning.md` — Lustig pipeline.
- `reference/CS_Implementation.md` — our CS theory (objective, λ traps, DCF precond).
- `helpers/recon/AGENTS.md` — contracts + pitfalls (read before touching recon code).
- Memories: `eye_vs_metric` (scalar traps, 3×), `slice_matching_zaffine` (the ~7-23%
  size / z-affine gotcha).
