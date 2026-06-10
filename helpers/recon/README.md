# FINUFFT Recon Baseline

Own thin pipeline + library NUFFT operator (handoff D1/D3). Native arm64 —
**do not run with the conda base python** (it is x86_64 under Rosetta).

## Setup (done 2026-06-10)

```bash
cd workspace/helpers
/opt/homebrew/bin/python3.11 -m venv .venv      # native arm64
.venv/bin/pip install finufft sigpy numpy scipy matplotlib
```

## Files

| File | Role |
|------|------|
| `asap_recon.py` | Module: `recon(traj, data, sample_weights=None, method=...)`; adjoint / Pipe–Menon DCF / CG via FINUFFT type1+type2; Steve-grid-units → radians conversion |
| `selftest.py` | Synthetic validation: adjointness dot-product test (machine precision), quality ordering adjoint < +DCF < CG. No scanner data needed |
| `compare_baseline.py` | The arbiter experiment: recon Steve's npy dumps, compare vs `savedbin0.npy` over all axis-flips, write side-by-side slices |

## Conventions (read before touching)

- Input trajectory in **Steve grid units** (`k·MS/IS + MS/2`); module converts to radians internally. MS=240, IS=100 defaults from `gtypes.py`.
- Default `isign=-1` matches Steve's forward-FFT convention voxel-for-voxel; textbook recon is `isign=+1` (conjugate/flip twin).
- Trajectory auto-tiled to data length (Steve's `kidx = idx % nuniq`).
- Zeros in data (spike/exclusion masking) contribute nothing — no special-casing.
- Bins = `sample_weights` vectors, one `recon()` call per bin. No bin machinery.

## Status / next

- [x] Operator pair validated (selftest PASS, 2026-06-10)
- [ ] Arbiter run on real phantom data — needs Steve's npy dumps
      (`trajx/y/z.npy`, `acq.npy` from one `dyn_recon` run, results.py:258-262)
      and his `savedbin0.npy` output in the same folder
- [ ] Confounder-neutralized comparison (gplb=0, killpts) per
      `reference/Recon_Comparison_StaticGas.md` protocol
- [ ] CS layer (sigpy L1-wavelet/TV) on top of the same operator
