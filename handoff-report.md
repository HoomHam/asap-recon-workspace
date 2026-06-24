# ASAP Workspace Handoff — non-CS (Faraz compare, Steve-vs-ours)

**Date:** 2026-06-24 · **Workstream:** ASAP-side non-CS workspace — Faraz fork,
Steve-vs-ours baseline, cross-impl comparison figures.
**Previous handoff (full prior detail, incl. CS):** `handoffs/handoff-workspace-2026-06-24.md`.

> **CS work decoupled to `2026_XeCS_Recon` on 2026-06-24** (operators, 4D pipeline,
> BART/Lustig comparison). **CS open loops now live there**, not here:
> `2026_XeCS_Recon/handoff-report.md` (library) + `2026_XeCS_Recon/workspace/handoff-report.md`
> (comparison). This handoff covers only the ASAP-side non-CS work. Steve/Tyger cloud
> pipeline → root `../handoff-report.md`. Migration record:
> `reference/start_here_to_decouple.md`, `reference/Repo_Architecture.md` (memory `asap_xecs_migration`).

---

## What this workspace still owns
- **`faraz/`** — Faraz's MATLAB fork (was `codes/2023_Faraz_Recon_HH/`). Entry
  `faraz/spiral_human_20240227.m`; gridding `faraz/gridrecon_fa_20230113.m`.
- **`helpers/recon/`** (stayers) — `compare_baseline.py` (FINUFFT vs Steve baseline),
  `steve_kernel_numpy.py` (Steve GPU kernel CPU reimpl), `dump_inputs*.py` (Steve data
  prep), `faraz_montage.py` / `faraz_zoom_check.py`, `cs_montage.py` (CS-vs-Steve-vs-Faraz,
  ASAP side). These import CS operators from XeCS via `xecs_recon.pth` in `helpers/.venv`.
- **`reference/`** — Steve-vs-Faraz truth: `Recon_Comparison_StaticGas.md`,
  `Recon_Overview_Steve.md`, `Recon_Overview_Faraz.md`. (CS theory docs are mirrored to
  XeCS; the comparison report `Final_Report_CS_Comparison.md` stays here.)

## Established, do not redo
- **Steve vs Faraz static comparison — FINALIZED** (`reference/Recon_Comparison_StaticGas.md`,
  single source of truth: magnitude-vs-real combine, coil-combine, de-apodization).
- Preserved CS originals in `helpers/_delete/` (byte-verified vs XeCS) — safe to
  `git rm` once XeCS is in daily use.

## Open loops (ASAP non-CS only)
| # | Item | Where it stands |
|---|------|-----------------|
| R3 | **Steve rebuttal to `for_steve.md`** | Pending — tell Steve about the DIAPHRAGM fix + ask for review of the comparison. |
| R4 | **Tell Faraz about the ×1.205 zoom bug** | Faraz doesn't know his recon has a scale error (documented in `helpers/recon/faraz_zoom_check.py`). |
| B1 | **Arbiter run: ours vs Steve on real phantom** | `helpers/recon/compare_baseline.py` — blocked on a `.dat` + gas trajectory `.npy`; see `helpers/recon/README.md`. |
| D1 | **Auto_Faraz_Recon.md instruction doc** | Missing. Write in `instructions/`, style of `instructions/Auto_Steve_Recon.md`. Entry `faraz/spiral_human_20240227.m`. |

## Hard constraints
- Never `git commit`/`push` in repo ROOT (`2026_ASAP_Recon/`). Git only in this
  `workspace/` (or in the `2026_XeCS_Recon` repo). 
- Never modify Faraz's code or Steve's upstream code.
- ASAP stayers run with `helpers/.venv/bin/python`; CS operators resolve from XeCS via
  `xecs_recon.pth` (if `import asap_recon` fails, that `.pth` is missing).
- ALWAYS z-affine register before any cross-pipeline slice comparison.
- `cs_montage.py` reads maxeig/t_ref from `cs_sweep_metrics.json` (written by XeCS
  `cs_recon.py`) — rerunning the XeCS sweep silently changes later montages.

## Key paths
| Item | Path |
|------|------|
| Faraz MATLAB fork | `faraz/` (entry `spiral_human_20240227.m`, map `faraz/AGENTS.md`) |
| ASAP stayers | `helpers/recon/` (`compare_baseline`, `steve_kernel_numpy`, `cs_montage`, `faraz_*`, `dump_inputs*`) |
| Steve-vs-Faraz truth | `reference/Recon_Comparison_StaticGas.md` |
| Comparison report | `reference/Final_Report_CS_Comparison.md` |
| Preserved CS originals | `helpers/_delete/` |
| Kento C kernel | `codes/kasap.c` (vs `../asap/asap.c`) |
| Memories | `eye_vs_metric`, `slice_matching_zaffine`, `workspace_system`, `asap_xecs_migration` |

## Suggested skills
- `/handoff update <note>` when the Faraz/Steve items resolve or the arbiter run lands.
- `/intent-layer-maintenance` — AGENTS.md paths shifted in the decouple.
- CS open loops → `2026_XeCS_Recon` handoffs (root + workspace).
