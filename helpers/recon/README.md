# FINUFFT Recon Baseline

> **Restructure 2026-06-24:** the CS operator library moved to `2026_XeCS_Recon`.
> `asap_recon.py`, `cg_tune.py`, `selftest.py` (+ the CS/4D/comparison code) now live
> there; the ASAP scripts below import them via the `xecs_recon.pth` in `../.venv`.
> Originals preserved in `../_delete/`. See `AGENTS.md`.

Own thin pipeline + library NUFFT operator (handoff D1/D3). Native arm64 —
run with this project's `.venv`. (Historical note: the old Intel anaconda base
was x86_64 under Rosetta, which is why this venv exists; since 2026-06-11 the
system conda base is native arm64 Miniforge, but keep using the venv — it pins
the working finufft/sigpy stack.)

## Setup (done 2026-06-10)

```bash
cd workspace/helpers
/opt/homebrew/bin/python3.11 -m venv .venv      # native arm64
.venv/bin/pip install finufft sigpy numpy scipy matplotlib
```

## Files

Stayers (here):

| File | Role |
|------|------|
| `dump_inputs.py` | **No-GPU input production**: runs Steve's own `raw.py`/`traj` loaders on a `.dat` + trajectory `.npy`, writes his exact `trajx/y/z.npy`, `acq.npy`, `bins.npy` + `meta.json`. Needs `pymapvbvd` (installed) |
| `dump_inputs_dyn.py` | dynamic-recon variant of the input dump |
| `steve_kernel_numpy.py` | Faithful CPU reimplementation of `cudarecon`/`cudarenorm` (single bin/channel): same filter, box, Gaussian, knorm kluge, F-order reshape, FFT, crop. ~1.6 s / 120k samples on 153³. Validated corr 0.95 vs synthetic truth |
| `compare_baseline.py` | The arbiter experiment: our recons vs Steve — uses GPU `savedbin0.npy` if present, else computes Steve-equivalent via the numpy kernel. Flip-search alignment, slice figure. Imports `asap_recon` from XeCS |
| `cs_montage.py` | CS vs Steve vs Faraz montage (cross-repo). Imports `asap_recon`/`cs_recon` from XeCS |
| `faraz_montage.py`, `faraz_zoom_check.py` | Faraz figure tools + ×1.205 zoom-bug doc |
| `convert_calib.py` | calibration format converter |

Moved to `2026_XeCS_Recon` (originals in `../_delete/recon/`):

| File | New home |
|------|----------|
| `asap_recon.py` | XeCS `recon/` — `recon(traj, data, sample_weights=None, method='cg'\|'adjoint')`; FINUFFT type1+type2; Steve-grid-units → radians. CG = method of record (iters=20, lam=0); DCF variant deleted 2026-06-11 |
| `cg_tune.py` | XeCS `recon/` — λ×iters sweep + legacy `metrics()` (in `cs_recon` import chain) |
| `selftest.py`, `selftest_4d.py` | XeCS `recon/` — adjointness + quality-ordering self-tests |
| `cs_recon.py`, `cs_recon_4d.py`, `binning.py`, `surrogates.py` | XeCS `recon/` |
| `cine_4d.py`, `export_4d.py`, `slice_video.py`, `nav_movie.py`, `surrogate_compare.py`, `diaphragm_check.py`, `kernel_check.py` | XeCS `pipeline/` |
| `bart_compare.py`, `bartio.py`, `lustig_compare.py`, `wavelet_twoway.py`, `tv_threeway.py`, `slice_matched_compare.py`, `z_register_compare.py`, `zreg_sixway_montage.py`, `resolution_sweep.py`, `metrics_v2.py` | XeCS `workspace/compare/` |

## Conventions (read before touching)

- Input trajectory in **Steve grid units** (`k·MS/IS + MS/2`); module converts to radians internally. MS=240, IS=100 defaults from `gtypes.py`.
- Default `isign=-1` matches Steve's forward-FFT convention voxel-for-voxel; textbook recon is `isign=+1` (conjugate/flip twin).
- Trajectory auto-tiled to data length (Steve's `kidx = idx % nuniq`).
- Zeros in data (spike/exclusion masking) contribute nothing — no special-casing.
- Bins = `sample_weights` vectors, one `recon()` call per bin. No bin machinery.

## Status / next

- [x] Operator pair validated (selftest PASS, 2026-06-10)
- [x] No-GPU path complete: `dump_inputs.py` (inputs via Steve's own loaders)
      + `steve_kernel_numpy.py` (Steve-equivalent output). Cloud CUDA now
      optional — only for one-time bit-faithfulness certification
- [ ] **Arbiter run on real phantom data** — blocked only on a `.dat` file +
      its gas trajectory `.npy`:
      `dump_inputs.py meas.dat gp_traj.npy out/ && compare_baseline.py out/`
- [ ] Confounder-neutralized comparison (gplb=0 via meta.json edit, killpts)
      per `reference/Recon_Comparison_StaticGas.md` protocol
- [ ] One-time Colab run to certify numpy kernel vs GPU savedbin0.npy
- [x] CS layer built: `cs_recon.py` (sigpy L1-wavelet FISTA / TV PDHG on the
      same finufft operator, DCF-preconditioned — unweighted A^H A stalls
      gradient solvers, density spread 1.8e6). First sweep 2026-06-12:
      wavelet t0.003-0.01 beats both bars (SNR 32.8-47.2 > 28.7,
      lowfreqCV 0.084-0.085 < 0.093). λ choice pending Hooman's visual
      verdict on `cs_sweep_sheet.png`
- [ ] CS vs independent Lustig-lineage CS (BART `pics`) — after λ verdict.
      Undersampling test shelved 2026-06-12; returns with 4D diaphragm
      binning (temporal phase)
