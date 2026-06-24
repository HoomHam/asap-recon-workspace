# helpers/recon — ASAP-Side Recon Tools (non-CS)

> TL;DR: The ASAP-anchored scripts — Steve data prep, FINUFFT-vs-Steve baseline,
> Faraz figures, and the cross-repo CS montage. The CS **operator library + 4D
> pipeline + CS comparison moved to the `2026_XeCS_Recon` repo (2026-06-24).** Run
> with `../.venv/bin/python`.

## Restructure note (2026-06-24)
The CS code that used to live here is now in `2026_XeCS_Recon/` (operators in
`recon/`, 4D pipeline in `pipeline/`, comparison in `workspace/compare/`, Lustig in
`workspace/lustig/`). The originals are preserved in `../_delete/` (safe to remove).
The scripts below that need CS operators (`asap_recon`, `cs_recon`, `cg_tune`) import
them from XeCS via `xecs_recon.pth` in `../.venv/lib/python3.11/site-packages/`.

## Purpose
Owns: the ASAP-anchored analysis — independent FINUFFT recon to arbitrate
Steve-vs-Faraz, Steve data prep, Faraz figures, and the cross-repo CS-vs-Steve-vs-Faraz
montage. Does not own: the CS operators (now `2026_XeCS_Recon`), Steve's code (repo
root, read-only), Faraz's MATLAB (`../../faraz/`), scanner data (`../../data/`).

## Files here (stayers)

| File | Role | CS import? |
|------|------|-----------|
| `dump_inputs.py` | produce recon_io from a `.dat` (Steve's loaders) | no |
| `dump_inputs_dyn.py` | dynamic recon_io dump | no |
| `compare_baseline.py` | FINUFFT vs Steve baseline check | `asap_recon` (XeCS) |
| `steve_kernel_numpy.py` | Steve GPU kernel validation (numpy reimpl) | no |
| `faraz_montage.py` | Faraz figure tools | `asap_recon` (XeCS) |
| `faraz_zoom_check.py` | Faraz ×1.205 zoom-bug documentation | `asap_recon` (XeCS) |
| `convert_calib.py` | calibration format converter | no |
| `cs_montage.py` | our CS vs Steve vs Faraz (cross-repo, ASAP side) | `asap_recon`, `cs_recon` (XeCS) |

## Entry Points

| Task | Command |
|------|---------|
| Sanity for the CS operator (now in XeCS) | `../../../2026_XeCS_Recon/workspace/.venv/bin/python ../../../2026_XeCS_Recon/recon/selftest.py` |
| Reproduce Steve-vs-ours baseline | `../.venv/bin/python compare_baseline.py ../../data/v3_fov250/recon_io` |
| CS vs Steve vs Faraz montage | `../.venv/bin/python cs_montage.py <recon_io> --t-rel 0.003` |
| Produce inputs from a new .dat | `../.venv/bin/python dump_inputs.py meas.dat gp_traj.npy out/` |

## Contracts
- Python = `../.venv/bin/python` (arm64, finufft 2.5.1 + sigpy 0.1.27). Never conda base.
- CS operators resolve from `2026_XeCS_Recon/recon` via `xecs_recon.pth` — if
  `import asap_recon` fails, that `.pth` is missing from `../.venv` site-packages.
- Image-quality conclusions require eyes on a figure (standing rule): twice in one
  session broken recons produced plausible images *and* metrics.

## Pitfalls
- **`cs_montage.py` input coupling**: reads maxeig/t_ref from `cs_sweep_metrics.json`
  (written by XeCS `cs_recon.py`); rerunning the sweep silently changes later montages.
- `phase_corrected_real` in `faraz_zoom_check.py` is a documented tombstone — do not
  revive (local-phase reference creates hollow blotches and fools CV).
- Operators no longer live here — edit them in `2026_XeCS_Recon/recon/`, then rerun
  its `selftest.py`. Editing the copies in `../_delete/` does nothing.

## Navigation
- Parent: `../../CLAUDE.md` (workspace rules)
- CS library: `2026_XeCS_Recon/recon/AGENTS.md`
- CS comparison narrative: `../../reference/Final_Report_CS_Comparison.md`
- Comparison truth (Steve vs Faraz): `../../reference/Recon_Comparison_StaticGas.md`
- File-by-file roles: `./README.md`
