# FINUFFT Recon Baseline

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

| File | Role |
|------|------|
| `asap_recon.py` | Module: `recon(traj, data, sample_weights=None, method='cg'\|'adjoint')`; FINUFFT type1+type2; Steve-grid-units → radians conversion. CG = method of record (iters=20, lam=0 per sweep); DCF variant deleted 2026-06-11 |
| `cg_tune.py` | λ×iters sweep (2026-06-11). Verdict: Tikhonov λ is a no-op on fully-sampled data; gplb filter accounts for ~2 SNR pts of Steve's lead; the rest is bias–variance (his gridder smooths) → smoothing regularizer (CS layer) is the real knob |
| `selftest.py` | Synthetic validation: adjointness dot-product test (machine precision), quality ordering adjoint < +DCF < CG. No scanner data needed |
| `dump_inputs.py` | **No-GPU input production**: runs Steve's own `raw.py`/`traj` loaders on a `.dat` + trajectory `.npy`, writes his exact `trajx/y/z.npy`, `acq.npy`, `bins.npy` + `meta.json`. Needs `pymapvbvd` (installed) |
| `steve_kernel_numpy.py` | Faithful CPU reimplementation of `cudarecon`/`cudarenorm` (single bin/channel): same filter, box, Gaussian, knorm kluge, F-order reshape, FFT, crop. ~1.6 s / 120k samples on 153³. Validated corr 0.95 vs synthetic truth |
| `compare_baseline.py` | The arbiter experiment: our recons vs Steve — uses GPU `savedbin0.npy` if present, else computes Steve-equivalent via the numpy kernel. Flip-search alignment, slice figure |

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
- [ ] CS layer (sigpy L1-wavelet/TV) on top of the same operator
