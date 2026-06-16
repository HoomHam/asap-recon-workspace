# helpers/lustig_oneshot — one-command Lustig MATLAB CS

> TL;DR: reproduces Hooman's pre-project Lustig CS (sparseMRI + IRT) on a
> `recon_io` folder in **one command**, replacing the old
> MATLAB→Colab→Drive→MATLAB dance. `run_lustig.py <recon_io>`. README has the
> walkthrough; `../../reference/Lustig_CS_Tuning.md` has the parameter/sweep
> reference; `../../reference/Lustig_CS_Baseline.md` has the verdict.

## Purpose

Owns: the exact, scripted reproduction of the Lustig CS comparison baseline.
Does not own: the CS *method-of-record* (that's `../recon/`), the Lustig source
toolboxes (`../../codes/2025_CS/`, read-only), scanner data (`../../data/`).

## Entry Points

| Task | Command |
|------|---------|
| Full one-shot on a recon_io folder | `.venv_lustig/bin/python run_lustig.py ../../data/v3_fov250/recon_io` |
| Rebuild ACR_test only (DCF step) | `.venv_lustig/bin/python build_acrtest.py <recon_io> <out.mat>` |
| Reuse existing ACR_test, skip torch | `run_lustig.py <recon_io> --skip-dcf` |
| Montage result vs our CS | `../recon/lustig_compare.py <recon_io> <recon_io>/lustig/lustig_cs.mat` |
| **Pure-TV λ sweep (matched to BART, metrics_v2)** | `.venv_lustig/bin/python run_lustig_sweep.py <recon_io> --tvweights 1e-3 3e-3 1e-2 3e-2 1e-1` → `run_cs_sweep.m` → `<recon_io>/lustig/lustig_tv_sweep.mat` + metrics json |

## Contracts

- DCF step uses **`.venv_lustig`** (torch + torchkbnufft, arm64 CPU) — NOT the
  recon `../recon/../.venv` (finufft/sigpy). Kept separate on purpose; gitignored.
- The DCF MUST stay **torchkbnufft pipe** at `im_size=(100,100,100)` to match the
  old analysis. MATLAB `voronoidens` is a different algorithm — not a drop-in for
  exact reproduction (fine only as a sensitivity test).
- `run_cs.m` is a **byte-for-byte** headless copy of `spiral3d_cs_3D_hoom.m`
  (TV+identity-L1, 0.01/0.01, 15 iters). Do not "improve" it — its job is
  faithful reproduction. Tuning experiments go in new scripts, not here.
- Trajectory exactness hinges on the **center offset**: recon_io is grid-index
  (`g=k·MS/IS+MS/2`); subtract MS/2 before the notebook's scale-invariant
  max-radius normalization. Done in `build_acrtest.py`.

## Pitfalls

- `tv01g01.mat` in `../../codes/2025-09-24_ACR/` was the **old 2025-08-16 ACR**
  scan (425 984 samples), not v3_fov250 (424 320). Use this tool's output for
  same-data comparisons.
- `run_cs.m` saves v7 (`save(out_mat,'gas')`, ~120 MB < 2 GB) so `scipy.io.loadmat`
  reads it. An earlier `-v7.3` save needed h5py; the compare/driver handle both.
- MATLAB must be on PATH as `matlab` (or pass `--matlab`); `run_cs.m` auto-adds
  IRT (`setup.m`) + sparseMRI paths from `codes_root`.
- Verdict: this pipeline is **softer** than our finufft+FISTA CS (stalled NCG,
  DCF only in init). Reproduced for the comparison set, not because it wins.

## Navigation

- Parent: `../../CLAUDE.md` (workspace rules)
- Our CS method-of-record: `../recon/AGENTS.md`
- Params + what to sweep: `../../reference/Lustig_CS_Tuning.md`
- Pipeline + verdict: `../../reference/Lustig_CS_Baseline.md`
