# helpers/recon — FINUFFT Baseline + CS Layer

> TL;DR: Hooman's own recon pipeline. Validated finufft operator pair
> (`asap_recon.py`) + CG method-of-record + CS layer (`cs_recon.py`).
> Run everything with `../.venv/bin/python`. `README.md` has file roles;
> `../../reference/CS_Implementation.md` has the deep theory.

## Purpose

Owns: the independent reconstruction used to arbitrate Steve-vs-Faraz and to
carry the CS phase (handoff 2026-06-11 spec). Does not own: Steve's code
(repo root, read-only), Faraz's MATLAB (external), scanner data
(`../../data/`, gitignored).

## Entry Points

| Task | Command |
|------|---------|
| Sanity after any operator change | `../.venv/bin/python selftest.py` |
| Reproduce Steve-vs-ours baseline | `compare_baseline.py ../../data/v3_fov250/recon_io` |
| CS λ sweep + contact sheet | `cs_recon.py <recon_io> --max-iter 100` |
| CS vs Steve vs Faraz montage | `cs_montage.py <recon_io> --t-rel 0.003` |
| Produce inputs from a new .dat | `dump_inputs.py meas.dat gp_traj.npy out/` |
| **Fixed quality metrics (use these)** | `metrics_v2.py <recon_io>` (corner-ROI SNR, half-max extent, edge_sharp) |
| **ours vs BART pics sweep** | `bart_compare.py <recon_io> --bart ~/bin/bart-src/bart` |
| **Lustig TV λ sweep (matched to BART)** | `../lustig_oneshot/run_lustig_sweep.py <recon_io>` |
| **z-affine register before slice compare** | `z_register_compare.py <recon_io> --bart <bart>` |
| **FINAL 6-way registered montage** | `zreg_sixway_montage.py <recon_io> --bart <bart>` |

> **CS-comparison summary + code/figure catalog (for a blind agent / presentation):**
> `../../reference/Final_Report_CS_Comparison.md`. BART cfl I/O: `bartio.py`.
> TV three-way: `tv_threeway.py`. Wavelet two-way: `wavelet_twoway.py`.
> Per-pipeline λ sweep at a slice: `resolution_sweep.py`. ALWAYS z-affine register
> (sizes differ ~23%) — flip-only montages (`slice_matched_compare.py`,
> `tv_threeway`, `wavelet_twoway`) are superseded for slice claims.

## Contracts

- Python = `../.venv/bin/python` (arm64, pins finufft 2.5.1 + sigpy 0.1.27).
  Never the conda base.
- BART scripts (`bart_compare`, `z_register_compare`, `zreg_sixway_montage`,
  `wavelet_twoway`, `resolution_sweep`) need the BART binary at
  `~/bin/bart-src/bart` — built FROM SOURCE (no homebrew formula). Build recipe:
  `../../reference/Final_Report_CS_Comparison.md` §9 (gmake, gcc-15, OPENBLAS=1,
  PNG=0, cblas_openblas.h symlink). Pass `--bart <path>` if elsewhere.
- `selftest.py` must pass after any edit to `asap_recon.py`.
- λ/regularization is parameterized as a **soft-threshold in coefficient
  units** (relative to p99 of |W·x_cg|), never in objective-function units —
  objective-unit λ was a measured silent no-op (see Pitfalls).
- Gradient-type solvers (FISTA/PDHG) **require** the DCF weights `w = 1/|AAᴴ1|`
  in the data term. For undersampling experiments, recompute `dens` per
  sampling mask — never reuse the fully-sampled weights.
- Bins enter only as per-sample weights (decision D3); `recon()`'s API does
  not change for dynamic work. DCF and bin weights multiply into the same W.
- Image-quality conclusions require eyes on a figure (standing rule): twice
  in one session broken recons produced plausible images *and* metrics.

## Pitfalls

- **λ no-op trap**: FISTA's effective threshold is α·λ with α = 1/maxeig(AᴴA)
  (maxeig ≈ 6.6e9 here). λ scaled to data magnitudes lands orders below the
  coefficient scale and silently does nothing — TV output was byte-identical
  across 4 decades of λ.
- **Conditioning trap**: unweighted AᴴA spectrum spans the sample-density
  spread (1.8e6 on this trajectory). CG converges in ~20 iters regardless
  (Krylov); FISTA/PDHG stall ~1000× below solution amplitude and the
  unconverged output still looks like an image.
- **SNR metric blind spot**: sparsity priors zero the background, collapsing
  σ_bg and inflating SNR (t=0.1 wavelet: SNR 68 on a visibly worse image).
  Extents from `cg_tune.metrics` are threshold-fragile the same way.
- **`cs_montage.py` input coupling**: it reads maxeig/t_ref from
  `cs_sweep_metrics.json`; rerunning `cs_recon.py` rewrites that file and
  silently changes later montages.
- **Display units**: the three pipelines output unrelated absolute scales;
  montages normalize each volume by p99.5 onto a shared [0,1] axis —
  structure comparable, absolute intensity not.
- `phase_corrected_real` in `faraz_zoom_check.py` is a documented tombstone —
  do not revive (local-phase reference creates hollow blotches and fools CV).

## Status (2026-06-15, updated)

CS layer built and swept; wavelet t0.003–0.01 beats both handoff bars. λ verdict
**still pending** Hooman's eye.
- **(1) metrics — FIXED in NEW module `metrics_v2.py`** (did NOT edit
  `cg_tune.metrics`; old sweep JSONs must stay reproducible). Corner-ROI noise
  (SNR no longer inflates silently; `bg_collapsed` flag), half-max extent
  (CG-20 now 172/168/128 mm, was the full-FOV 250/250/250 bug), new `edge_sharp`
  guard so SNR can't be won by blurring. `metrics_v2.py <recon_io>` demos
  old-vs-new on CG-20. New code should import metrics from here.
- **(2) BART `pics` — RAN** (the real independence test). `bart_compare.py` +
  `bartio.py`. BART built from source at `/Users/hoomham/bin/bart-src/bart`.
  BART reproduces the phantom (extent/orientation cc≈0.92); our wavelet lfCV
  0.157 < BART floor ~0.177, BUT a grid/moiré texture on ALL BART recons (wav+tv)
  confounds the CV head-to-head — see `../../reference/BART_Comparison.md` (run
  table, the wavelet-needs-ADMM fix, the open texture artifact + next steps).
  MATLAB Lustig comparison already done (ours wins; `../lustig_oneshot/`).
- (3) temporal 4D with diaphragm binning — undersampling returns there.

Sibling tools this phase: `metrics_v2.py` (fixed metrics), `bart_compare.py`
(our CS vs BART pics sweep), `bartio.py` (cfl/hdr I/O), `lustig_compare.py`.

## Pitfalls (added 2026-06-15)
- **BART wavelet diverges under FISTA here** (~1e32). pics `-R W` must run with
  `-m` (ADMM); `bart_compare.run_pics` forces it. TV is fine on the default.
- **BART traj units:** `bart_traj = traj_rad * N/(2π)` (N=100), dims (3,M,1).
- **The "BART grid" is the REAL ACR resolution insert (z≈50-54), NOT an artifact.**
  MANDATORY: phantoms are different SIZES, so flip-orientation aligns the bulk but
  NOT the thin insert slice. Fitted z-affine to ours: BART scale 1.02/+1.5slice,
  Lustig 0.93 (7% smaller). Fixed-z compare is INVALID — it faked "BART resolves,
  ours blurs". After z-affine reg (`z_register_compare.py`) the grid aligns and
  at the matched slice ALL THREE resolve it; ours wav t0.003 as sharp or sharper
  than BART. Every prior winner call RETRACTED. Always z-affine register
  (scale+shift) before slice compare. (`zreg_finez_montage.png`,
  `grid_slice_zoom_z52.png`, `BART_Comparison.md` FINAL VERDICT).
- **metrics_v2 is blind to resolution loss; judge the insert slice by EYE.** Best-λ
  is a resolution call at matched z, not min-lfCV.
- **`metrics_v2` ≠ `cg_tune.metrics`:** corner-air noise ROI assumes a centered
  object (true for ACR). Don't mix the two metrics' numbers in one table.
- **`edge_sharp` guards blur, NOT texture.** Noise/texture (BART moiré, Lustig
  grain) inflates it — don't read a higher edge_sharp as finer detail.
- **Three-way TV is a wash** (ours/BART/Lustig all ~0.18 lfCV under metrics_v2);
  differentiation is the WAVELET prior (ours 0.157 < BART 0.177; Lustig's
  Wavelet is 2D-only, can't play). See `tv_threeway.py` +
  `../../reference/BART_Comparison.md`. Lustig TV sweep: `lustig_oneshot/
  run_lustig_sweep.py` (pure TV, xfmWeight=0, metrics_v2).
- **Orientation:** BART → `orient_to_ours` (flips+1 transpose, cc~0.92). Lustig
  needs the FULL 48-orientation search (`tv_threeway.best_orient_full`) — its
  per-slice rot90 mixes axes beyond the BART set. Metrics are orientation-
  invariant; only the montage PLANE needs the match.

## Navigation

- Parent: `../../CLAUDE.md` (workspace rules)
- Theory: `../../reference/CS_Implementation.md`, `../../reference/Physics_Notes.md`
- File-by-file roles: `./README.md`
- Comparison truth: `../../reference/Recon_Comparison_StaticGas.md`
