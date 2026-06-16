# ASAP Recon — Workspace Handoff

**Date:** 2026-06-16 · **Branch context:** workspace git (personal, `HoomHam/asap-recon-workspace`), NOT the main repo. Latest push `f31f5a8`.

> **This session jumped ahead to the LAST item on the list** — Task C, the 4D /
> temporal CS recon — and built it end to end. The comparison-phase tasks (A
> super-validation, B the ~23% size question) were **deliberately deferred** and are
> still open; everything we had for them is kept verbatim below. Pick up whichever
> workstream you want next.

## ⭐ What we did this session — 4D dynamic CS (old Task C, now BUILT)

**Read:** `reference/Dynamic_4D_CS_Implementation.md` — the AS-BUILT doc (modules,
data facts, the 4D objective as coded, the full diaphragm method + every hard-won
decision, run recipe, outputs). This handoff does NOT duplicate it.

One-line summary: extended our static CS operator to a B-bin respiratory cine
(per-bin wavelet + circular temporal-TV PDHG, single coil), with three respiratory
surrogates (signal / pneumotach / diaphragm-via-CS-nav). Validated on real human
data (027JC signal cine; 025JC surrogates). 10 new modules in `helpers/recon/`
(`dump_inputs_dyn, binning, cs_recon_4d, selftest_4d, surrogates, cine_4d,
diaphragm_check, kernel_check, nav_movie, surrogate_compare`). `selftest_4d` = ALL PASS.

Hard-won gotchas now written into the implementation doc (do not re-break):
SI axis = **axis2** (corr 0.94); diaphragm nav `win_ilv=20`/`smooth_win=5`, **no
median filter**; diaphragm **curve** = edge metric with **correlation-based boundary
selection** (clean, period 3.5 s) — NOT hardcoded hi; the **hi edge is the anatomical
dome but clips out of FOV at inspiration** so it's only the nav-movie *display* line;
all orientation flips are **display-only, never on the matrix**.

### NEXT (4D CS workstream), in order
1. **Test + finalize the diaphragm surrogate** (the immediate follow-up). Confirm the
   corr-selected edge curve is robust across patients; decide the final binning
   surrogate (corr-selected edge vs centroid — both clean); apply the fixes to **027JC**
   (its diaphragm was run before the axis/window/selection fixes). Lock `diaphragm_check`
   + `nav_movie` conventions.
2. **Temporal parameter tuning** — sweep `lam_t_rel` (e.g. 0.02/0.05/0.1/0.2) on the
   joint cine; pick by eye (memory `eye_vs_metric`). Also revisit `lam_s_rel`, `bindt`
   (soft-bin width), and B (bins).
3. **Run the full cines for all three surrogates** on 027JC and 025JC (`cine_4d.py
   --surrogate signal|pneumo|diaphragm --stage both`); compare the three binnings.
4. **Natural extensions** (deferred in the impl doc): low-rank temporal option (swap/add
   to temporal-TV); dissolved-phase cine (RBC/TP, `results.py:303`); retrospective
   undersampling experiment (recompute DCF per mask — `recon/AGENTS.md` contract).
5. Graphify is **stale** (pre-4D, Jun 15 21:27) and not auto-updating — re-run
   `/graphify` to refresh the knowledge graph over the new modules.

---

## Comparison phase (PRIOR session) — FINALIZED, do not redo

> **Read:** `reference/Final_Report_CS_Comparison.md` (narrative, per-script I/O §7,
> figure catalog §8, reproduce + BART build §9, open threads §10).

1. **`metrics_v2.py`** — fixed quality metrics (corner-ROI SNR, half-max extent,
   `edge_sharp`). Replaces `cg_tune.metrics` WITHOUT editing it. Caveat: **blind to
   resolution loss** (open issue).
2. **BART built from source** at `~/bin/bart-src/bart` (recipe = Final Report §9).
   `bart_compare.py` runs ours-vs-BART λ sweep. BART **wavelet MUST use ADMM (`-m`)**.
3. **Lustig TV λ-sweep** (`lustig_oneshot/run_lustig_sweep.py` → `run_cs_sweep.m`),
   scored with metrics_v2. Lustig's NCG stalls (obj flat).
4. **Slice-matching solved.** Phantoms differ in SIZE across pipelines (z-scale rel.
   our CG: Faraz 1.16, BART 1.04, Steve 1.01, our wav 0.98, Lustig 0.94 → **~23%
   spread**). Fixed-z compare was INVALID; fixed with **z-affine registration**
   (`z_register_compare.py`).
5. **FINAL verdict (slice-matched):** all pipelines resolve the real ACR grid; **ours
   as sharp as BART, sharper than Lustig.** Figures: `zreg_sixway_montage.png`,
   `grid_slice_zoom_z52.png`, `zreg_finez_montage.png`.
6. **Final report + Intent Layer** written/updated.

### RETRACTED (don't resurrect in any report/manuscript)
- "ours-wavelet wins on lfCV (0.157<0.177)" — blur erasing real resolution.
- "BART has a grid/moiré artifact" — it's the REAL ACR resolution insert; BART resolves it.
- "DCF removes the grid" / "grid is unweighted-LS signature" — disproven. Trail in `BART_Comparison.md`.

---

## Carried over / NOT finalized (open loops — still open)

| # | Item | Where it stands |
|---|------|-----------------|
| O1 | **λ verdict (ours wavelet t0.003 vs t0.01)** | still Hooman's eye; a RESOLUTION call at the insert slice |
| O2 | **BART per-sample DCF** | `pics -p` did not wire noncart density; would match ours' weighting |
| O3 | **Lustig wavelet** | not possible — `FWT2_PO` (2D) on a 3D volume. Wavelet stays ours-vs-BART |
| O4 | **b-map stage** (Steve's `calcb` port) | for general multi-coil. **NOTE: 027JC/025JC are nch=1**, so 4D did NOT need it; still unported for multi-coil data |
| O5 | Steve's rebuttal to `for_steve.md`; tell Faraz about ×1.205 zoom bug | pending |
| O6 | `recon/*.py` token bloat | now larger with the 4D modules; consider splitting `helpers/recon/compare/` + `helpers/recon/dyn4d/` child nodes |

---

## DEFERRED comparison tasks (still the plan for that workstream)

### Task A — SUPER-validate every comparison script (presentation/manuscript-grade)
Be certain the code behind every reported number/figure is correct. Scripts in
`helpers/recon/` + `helpers/lustig_oneshot/` (catalog: Final Report §7). Validate:
1. **Operator:** `../.venv/bin/python selftest.py` passes (adjoint `<Ax,y>=<x,Aᴴy>`).
2. **metrics_v2 unit-check:** synthetic phantoms with KNOWN snr/extent/edge; `bg_collapsed`
   fires on zeroed bg; `edge_sharp` drops under blur. (write `test_metrics_v2.py`.)
3. **bartio round-trip:** `writecfl`→`readcfl` identity (cross-check BART `python/cfl.py`).
4. **Orientation/registration:** feed `best_orient_full`/`fit_z_affine` a KNOWN
   transformed volume → must recover the inverse. (most error-prone; drove the verdict.)
5. **BART invocation:** `bart_traj=traj_rad·N/(2π)` via point-source PSF centered; `-m`
   ADMM convergence.
6. **Lustig reproduction:** diff `run_cs_sweep.m` vs byte-for-byte `run_cs.m` (only the
   TVWeight loop + `xfmWeight=0`).
7. **Cross-pipeline sanity:** every recon hits ACR extent ~190/190/148 mm; flag misses.
Deliverable: `reference/Validation_Report.md`. Consider parallel audit subagents.

### Task B — explain the ~23% phantom-size difference
Each pipeline reconstructs at a slightly different physical size (z-scale Faraz 1.16 …
Lustig 0.94). Root cause = FOV / k-space-scaling / oversampling convention. Trace each
pipeline's trajectory-units→grid mapping (ours `grid_to_radians`; BART [-N/2,N/2];
Lustig max-radius norm; Steve `g=k·MS/IS+MS/2`; Faraz `KSpaceCoor` 1/mm). A point-source
PSF test (Task A.5) gives per-pipeline FOV directly. Pin before any quantitative claim.

---

## Hard constraints (carry forward)
- Run ours with `../.venv/bin/python` (arm64; finufft 2.5.1, sigpy 0.1.27). Lustig DCF
  step uses `lustig_oneshot/.venv_lustig` (torch). NEVER conda base (x86/Rosetta).
- NEVER git commit/push in repo ROOT. All git in `workspace/` only (pushed: `f31f5a8`).
- Never modify Faraz's code, Steve's upstream code, or `run_cs.m` (byte-for-byte original).
- Read `recon/AGENTS.md` pitfalls before touching recon code — z-affine / eye-over-metric
  / BART-needs-ADMM; and `Dynamic_4D_CS_Implementation.md` for the 4D/diaphragm gotchas.
- ALWAYS z-affine register before any cross-pipeline slice comparison (~23% size spread).
- Correct dynamic trajectory = `fa_spiral_dyn_fancy_v3_20240130_{gp,dp}.npy` (date in
  protocol name); display orientation flips are plot-only, never on the matrix.

## Suggested skills
- `/graphify` — refresh the stale knowledge graph over the new 4D modules.
- `/code-review` on the new 4D modules (`cs_recon_4d.py`, `surrogates.py`, `cine_4d.py`)
  before they back a manuscript, and on the comparison scripts (Task A).
- `/intent-layer-maintenance` if `helpers/recon/` is split into `dyn4d/` + `compare/` nodes.
- `/handoff update <note>` as diaphragm finalization / temporal tuning lands.

## Key paths
- 4D AS-BUILT (read first for the new work): `reference/Dynamic_4D_CS_Implementation.md`
- Comparison report: `reference/Final_Report_CS_Comparison.md`; dead-ends: `BART_Comparison.md`
- Code: `helpers/recon/` (+ `lustig_oneshot/`); node docs: their `AGENTS.md`
- 4D dumps/outputs: `data/v3_dyn/recon_io_dyn/` (027JC), `data/v3_dyn_025JC/recon_io_dyn/` (025JC)
- Comparison figures/JSON: `data/v3_fov250/recon_io/` (+ `bart/`, `lustig/`)
- BART binary: `~/bin/bart-src/bart`
- Memories: `eye_vs_metric` (scalar traps ×3), `slice_matching_zaffine` (z-affine gotcha)
