# lustig_oneshot — one-command Lustig CS on recon_io data

Replaces the old 3-stage dance (run `spiral_gpt_ACR_fovdes` in MATLAB → upload
`ACR_data.mat` to Drive → run `spiral3d_frames_mat_hoom.ipynb` in Colab →
download `ACR_test.mat` → run `spiral3d_cs_3D_hoom.m`) with **one command** that
reads a `recon_io` folder directly and reproduces the old analysis exactly.

## Run

```bash
.venv_lustig/bin/python run_lustig.py ../../data/v3_fov250/recon_io
```

Outputs to `<recon_io>/lustig/`:
- `ACR_test.mat` — ktrajs/kdatas/kcomps (exact notebook format)
- `lustig_cs.mat` — `gas(15,100,100,100)`, the 15 fnlCg iterations

Flags: `--out-dir DIR`, `--skip-dcf` (reuse existing ACR_test, skip torch step),
`--matlab /path/to/matlab`.

## Pieces

| File | Role |
|------|------|
| `build_acrtest.py` | recon_io npy → ACR_test.mat. Exact replica of `spiral3d_frames_mat_hoom.ipynb`: recenter grid-index traj (−MS/2), max-radius normalize to [−π,π], **torchkbnufft pipe DCF** at `im_size=(100,100,100)`. Runs in `.venv_lustig`. |
| `run_cs.m` | Headless exact copy of `spiral3d_cs_3D_hoom.m`. NUFFT3D + fnlCg, TVWeight=xfmWeight=0.01, 15 iters, per-slice rot90, saves `gas`. Auto-adds IRT + sparseMRI paths. |
| `run_lustig.py` | Driver: build → `matlab -batch run_cs` → metrics. |
| `.venv_lustig/` | torch + torchkbnufft (arm64 CPU) for the DCF step only. Not the recon `.venv`. |

## Exactness notes (why output reproduces tv01g01)

- Trajectory: recon_io stores grid-index (`g = k·MS/IS + MS/2`); the notebook's
  `k/|k|.max()·π` normalization is scale-invariant, so only the **center offset**
  matters → subtract MS/2 (=120 here). Geometry then matches both the old Lustig
  run and our finufft pipeline.
- **DCF is torchkbnufft pipe** — not MATLAB `voronoidens` (different algorithm).
  This is the one piece that forces the torch venv; it's load-bearing because
  fnlCg barely moves from the DCF-gridded init.
- `im_size=(100,100,100)` — the FIXED size. The `(256,64,256)` cardiac-leftover
  bug from the 2026-05 docs is **not** in this notebook version.

## Caveat carried from analysis

This Lustig pipeline (TV + identity-L1, stalled NCG) produces a softer recon
than our finufft+FISTA wavelet CS; metrics on v3_fov250 match the old-data run
(final lowfreq-CV ≈ 0.115 vs our 0.084). Reproduced here for the comparison set,
not because it wins. See `reference/Lustig_CS_Baseline.md`.
