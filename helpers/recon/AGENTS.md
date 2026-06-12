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

## Contracts

- Python = `../.venv/bin/python` (arm64, pins finufft 2.5.1 + sigpy 0.1.27).
  Never the conda base.
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

## Status (2026-06-12)

CS layer built and swept; wavelet t0.003–0.01 beats both handoff bars
(SNR > 28.7, lowfreq-CV < 0.093). λ verdict pending Hooman's eye. Next:
undersampling experiment (1/2, 1/4, 1/8 interleaves).

## Navigation

- Parent: `../../CLAUDE.md` (workspace rules)
- Theory: `../../reference/CS_Implementation.md`, `../../reference/Physics_Notes.md`
- File-by-file roles: `./README.md`
- Comparison truth: `../../reference/Recon_Comparison_StaticGas.md`
