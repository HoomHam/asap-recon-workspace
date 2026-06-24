# ASAP Workspace Handoff — CS recon (now in 2026_XeCS_Recon)

**Date:** 2026-06-24 · **Workstream:** Hooman's CS recon (finufft/CS, 4D dynamic, BART/Lustig compare).
**Previous handoff (full detail):** `workspace/handoffs/handoff-workspace-2026-06-24.md`.

> **Repo decouple: ✅ DONE 2026-06-24.** All CS code + CS workspace moved to the new
> repo **`2026_XeCS_Recon`** (operators `recon/`, 4D `pipeline/`, comparison
> `workspace/compare/`, Lustig `workspace/lustig/`). The ASAP-side originals are
> preserved (not deleted) in `helpers/_delete/`. This file now tracks the CS open loops
> — **most now land in `2026_XeCS_Recon`.** Migration record (reference):
> `reference/start_here_to_decouple.md`, `reference/Repo_Architecture.md`. Execution
> detail: archived handoff above + `2026_XeCS_Recon/workspace/handoffs/XeCS_Handoff_2026-06-24.md`.
> Decouple implementation walkthrough intentionally not repeated here.

---

## Where the CS work lives now
- Operators (`asap_recon`, `cs_recon`, `cs_recon_4d`, `binning`, `surrogates`,
  `cg_tune`, `selftest*`): `2026_XeCS_Recon/recon/`
- 4D pipeline (`cine_4d`, `export_4d`, `slice_video`, `nav_movie`, …): `2026_XeCS_Recon/pipeline/`
- Comparison (BART/Lustig/z-reg/`metrics_v2`): `2026_XeCS_Recon/workspace/compare/`
- Run with `2026_XeCS_Recon/workspace/.venv/bin/python` (py3.11.15; finufft 2.5.1,
  sigpy 0.1.27). `xecs_recon.pth` puts `recon/` on the path from any cwd. `selftest.py`
  + `selftest_4d.py` PASS.
- ASAP-side stayers (`helpers/recon/compare_baseline.py`, `cs_montage.py`, `faraz_*`)
  import operators from XeCS via an `xecs_recon.pth` in `helpers/.venv`.

## Established, do not redo
- **3D+t dynamic CS — IMPLEMENTATION FINAL.** wavelet_xyz (db4) + circular TV_t, PDHG,
  DCF-preconditioned, single-coil. As-built: `reference/Dynamic_4D_CS_Implementation.md`.
  `selftest_4d` = ALL PASS.
  - **Hard-won diaphragm conventions (do NOT re-break):** SI axis = axis2 (corr 0.94);
    `win_ilv=20`, `smooth_win=5`, no median filter. Surrogate auto corr-based selection
    (`prefer="auto"`); on 025JC lo-edge wins (corr 0.83), period 3.52s. `prefer="hi"` =
    DISPLAY line only (clips out of FOV at deep inspiration → not for binning). 027JC nav
    too low-quality → use signal surrogate. Display: sagittal=axis0 no-rot,
    coronal=axis1 90°CW, axial=axis2 90°CCW+fliplr.
- **Static CS comparison — FINALIZED.** `reference/Final_Report_CS_Comparison.md` +
  `reference/BART_Comparison.md`. Verdict: all pipelines resolve the ACR grid; **ours as
  sharp as BART, sharper than Lustig** (after z-affine registration; ~23% size spread).
  BART built from source `~/bin/bart-src/bart`; BART wavelet MUST use ADMM (`-m`).
  RETRACTED (don't resurrect): "ours-wavelet wins on lfCV", "BART moiré is an artifact",
  "DCF removes the grid".

## Open loops (carry forward — still unaddressed; land in XeCS unless noted)

| # | Item | Where it stands |
|---|------|-----------------|
| O1 | **λ verdict (ours wavelet t0.003 vs t0.01)** | Hooman's eye — RESOLUTION call at the insert slice, not min-lfCV. |
| O7 | **λ_t visual pick** | View `2026_XeCS_Recon/workspace/outputs/025JC_sweep_lt/` coronal MP4s; pick the λ_t that keeps inter-bin motion while suppressing per-bin streaks. (Sweep done: lam_t ∈ {0.003,0.01,0.02,0.05,0.1,0.2}, current default 0.05.) |
| ⚠️B0 | **B0 correction check (CS)** | Confirm the CS forward/adjoint uses the **correct-resolution** b-matrix (B0/off-res + coil phase), not the low-res b from the navigator loop. `cs_recon.py` vs `results.py:calcb` / `tyger_recon.py:_diaphragm_navigator`. **Do not lose this.** |
| O8 | **027JC full pipeline** | Nav quality too low for diaphragm surrogate (82/363 valid) → run the signal-surrogate cine instead. |
| O2 | **BART per-sample DCF** | `pics -p` did not wire noncart density. |
| O3 | **Lustig wavelet** | Not possible — `FWT2_PO` is 2D on a 3D volume; TV is the only fair 3-way. |
| O4 | **b-map stage (`calcb` port)** | 027JC/025JC nch=1 so 4D didn't need it; unported for multi-coil. |
| O9 | **data/ cleanup** | Sweep dirs + nav_movie dirs mixed into raw `data/v3_dyn_025JC/recon_io_dyn/`; raw acq/traj must stay, result files should move to `outputs/`. `data/` gitignored. |

### Deferred comparison tasks (→ XeCS workspace)
- **Task A — SUPER-validate every comparison script** (manuscript-grade): operator adjoint
  (`selftest`), `metrics_v2` synthetic-phantom unit check, `bartio` round-trip,
  orientation/`fit_z_affine` recovers a known transform, BART `bart_traj=traj_rad·N/(2π)`
  via PSF + `-m` convergence, Lustig `run_cs_sweep.m` vs `run_cs.m` diff, cross-pipeline
  ACR extent ~190/190/148 mm. Deliverable: `Validation_Report.md` in XeCS workspace.
- **Task B — Explain the ~23% phantom-size difference** — trace each pipeline's
  trajectory-units→grid mapping; PSF test (A) gives per-pipeline FOV directly.

## Instruction docs TODO (the new ask)
Write one `Auto_<Impl>_Recon.md` per main implementation, in the **style of
`workspace/instructions/Auto_Steve_Recon.md`** (build map → call graph → data provenance
→ bash run guide → comparison playbook). Needed in **`2026_XeCS_Recon/workspace/instructions/`**:
- `Auto_CS_Static_Recon.md` — entry `2026_XeCS_Recon/recon/cs_recon.py`
- `Auto_CS_4D_Recon.md` — entry `2026_XeCS_Recon/pipeline/cine_4d.py`
- `Auto_BART_Recon.md` — entry `2026_XeCS_Recon/workspace/compare/bart_compare.py` + `~/bin/bart-src/bart`
- `Auto_Lustig_Recon.md` — entry `2026_XeCS_Recon/workspace/lustig/run_lustig.py`

(ASAP side: `Auto_Faraz_Recon.md` for `workspace/faraz/spiral_human_20240227.m` — see root handoff.)
Also: recreate XeCS `workspace/lustig/.venv_lustig` (`torch torchkbnufft`) if the Lustig flow is needed there.

## Hard constraints (carry forward)
- Never `git commit`/`push` in repo ROOT (`2026_ASAP_Recon/`). Git only in `workspace/`
  or in the `2026_XeCS_Recon` repo. Fork branch `HoomHam/asap_recon` is separate — OK to push.
- Never modify Faraz's code, Steve's upstream code, or `run_cs.m`.
- ALWAYS z-affine register before any cross-pipeline slice comparison.
- Read `2026_XeCS_Recon/recon/AGENTS.md` + `workspace/compare/AGENTS.md` pitfalls and
  `reference/Dynamic_4D_CS_Implementation.md` before touching CS code.
- Trajectory = `fa_spiral_dyn_fancy_v3_20240130_{gp,dp}.npy`; display flips are plot-only.

## Key paths (post-restructure)
- CS library: `2026_XeCS_Recon/recon/` · 4D: `…/pipeline/` · compare: `…/workspace/compare/`
- CS venv: `2026_XeCS_Recon/workspace/.venv/bin/python`
- 4D as-built / theory: `reference/Dynamic_4D_CS_Implementation.md`, `reference/4D_CS_Theory_Limitations.md`
- Comparison report: `reference/Final_Report_CS_Comparison.md`, `reference/BART_Comparison.md`
- Migration record: `reference/start_here_to_decouple.md`, `reference/Repo_Architecture.md`
- Faraz MATLAB fork: `workspace/faraz/` · Preserved CS originals: `helpers/_delete/`
- Data dumps (shared, gitignored): `data/v3_dyn_025JC/recon_io_dyn/`, `data/v3_dyn/recon_io_dyn/`,
  `data/v3_fov250/recon_io/` · BART: `~/bin/bart-src/bart`
- Memories: `eye_vs_metric`, `slice_matching_zaffine`, `workspace_system`, `asap_xecs_migration`

## Suggested skills
- `/handoff update <note>` when λ/λ_t picks made, B0 check done, or instruction docs written.
- `/code-review` on the 4D CS modules before manuscript.
- `/intent-layer-maintenance` — AGENTS.md paths moved in the decouple.
