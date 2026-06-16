# ASAP Recon — Workspace Handoff (CS comparison phase CLOSED)

**Date:** 2026-06-15 · **Branch context:** workspace git (personal), NOT the main repo.
**Phase closed this session:** CS-vs-independent comparison (ours / BART / Lustig /
Steve / Faraz) on `data/v3_fov250/recon_io` (ACR phantom, static gas, single bin).
**Next round (this handoff's focus):** (1) SUPER-validate every script we wrote, then
(2) explain the ~23% phantom-size difference, then (3) the final 4D / temporal recon.

> **Read first:** `reference/Final_Report_CS_Comparison.md` — the presentation-ready
> report. It has the narrative arc, a per-script I/O catalog (§7), a figure catalog
> tagged FINAL/PROCESS/SUPERSEDED (§8), reproduce-from-scratch + BART build recipe
> (§9), and open threads (§10). This handoff does NOT duplicate it — it points at it.

---

## What was FINALIZED this session (do not redo)

1. **`metrics_v2.py`** — fixed quality metrics (corner-ROI SNR, half-max extent,
   `edge_sharp`). Replaces `cg_tune.metrics` WITHOUT editing it (old JSONs stay
   reproducible). Caveat carried forward: it is **blind to resolution loss** — see
   open issue. 
2. **BART built from source** at `~/bin/bart-src/bart` (recipe = Final Report §9).
   `bart_compare.py` runs the ours-vs-BART λ sweep. BART **wavelet MUST use ADMM
   (`-m`)** — FISTA diverged.
3. **Lustig TV λ-sweep** (`lustig_oneshot/run_lustig_sweep.py` → `run_cs_sweep.m`),
   scored with metrics_v2, matched to BART's sweep. Lustig's NCG stalls (obj flat).
4. **Slice-matching solved.** The phantoms are different SIZES across pipelines
   (z-scale rel. our CG: Faraz 1.16, BART 1.04, Steve 1.01, our wav 0.98, Lustig
   0.94 → **~23% spread**). Flip-only orientation aligns the bulk but NOT the thin
   ACR resolution-insert slice → fixed-z comparison was INVALID and faked a "BART
   resolves, ours blurs" result. Fixed with **z-affine registration (scale+shift)**:
   `z_register_compare.py` (`fit_z_affine`/`resample_z`).
5. **FINAL verdict (slice-matched):** at the matched grid slice ALL pipelines resolve
   the real ACR resolution grid; **ours is as sharp as BART, sharper than Lustig.**
   Independence test (BART) confirms our geometry, does not beat us. Figures:
   `zreg_sixway_montage.png`, `grid_slice_zoom_z52.png`, `zreg_finez_montage.png`.
6. **Final report + Intent Layer** written/updated (Final_Report, BART_Comparison,
   recon `AGENTS.md`, lustig `AGENTS.md`, workspace `CLAUDE.md` index).

### RETRACTED this session (don't resurrect in any report/manuscript)
- "ours-wavelet wins on lfCV (0.157<0.177)" — that was blur erasing real resolution.
- "BART has a grid/moiré artifact" — it's the REAL ACR resolution insert; BART resolves it.
- "DCF removes the grid" / "grid is unweighted-LS signature" — both disproven (ours-no-DCF
  is grid-free; the grid was a mismatched slice). Full dead-end trail in `BART_Comparison.md`.

---

## Carried over / NOT finalized (open loops)

| # | Item | Where it stands |
|---|------|-----------------|
| O1 | **λ verdict (ours wavelet t0.003 vs t0.01)** | still Hooman's eye; now a RESOLUTION call at the insert slice, not min-lfCV |
| O2 | **BART per-sample DCF** | `pics -p` did not wire noncart density in this build; would make BART weighting match ours |
| O3 | **Lustig wavelet** | not possible — its `Wavelet` is `FWT2_PO` (2D) on a 3D volume. Wavelet stays ours-vs-BART |
| O4 | **b-map stage** (Steve's `calcb` port) | needed before multi-coil / 4D — still unported |
| O5 | Steve's rebuttal to `for_steve.md`; tell Faraz about ×1.205 zoom bug | pending (pre-existing loops) |
| O6 | `recon/*.py` ≈ 16k tokens | approaching 20k Intent-Layer child-node threshold; future split → `helpers/recon/compare/` |

---

## NEXT SESSION — do these, in order

### Task A (FIRST) — SUPER-validate every script (presentation/manuscript-grade)
Goal: be certain the code behind every reported number/figure is correct, so reports
and a future manuscript are solid. The 10 new scripts are in `helpers/recon/` and
`helpers/lustig_oneshot/` (catalog: Final Report §7). Validate, don't assume:

1. **Operator correctness:** `../.venv/bin/python selftest.py` must pass (adjoint test,
   `<Ax,y>=<x,Aᴴy>`). This underpins ALL ours recons.
2. **metrics_v2 unit-check:** synthetic phantoms with KNOWN snr/extent/edge → confirm
   the metric returns them; confirm `bg_collapsed` fires on a zeroed background; confirm
   `edge_sharp` drops under deliberate Gaussian blur. (No test file exists yet — write
   `helpers/recon/test_metrics_v2.py`.)
3. **bartio round-trip:** `writecfl`→`readcfl` of a random complex array is identity
   (also cross-check against BART's own `~/bin/bart-src/python/cfl.py` if present).
4. **Orientation/registration correctness:** feed `best_orient_full`/`fit_z_affine` a
   KNOWN transformed copy of a volume (apply a flip+z-scale, then recover it) → must
   return the inverse. This is the most error-prone code and it drove the verdict.
5. **BART invocation:** confirm `bart_traj = traj_rad·N/(2π)` is right by gridding a
   delta/point-source trajectory and checking the PSF is centered (also answers O-size
   question below). Confirm `-m` ADMM convergence (objective decreasing).
6. **Lustig reproduction:** `run_cs.m` is byte-for-byte the original — confirm
   `run_cs_sweep.m` differs ONLY by the TVWeight loop + `xfmWeight=0` (diff them).
7. **Cross-pipeline sanity:** every recon should hit ACR extent ~190/190/148 mm after
   metrics_v2; flag any that don't.

Deliverable: a `reference/Validation_Report.md` (pass/fail per script + any bugs found
+ fixes). Consider spawning parallel `Explore`/general-purpose subagents per script to
audit in parallel (the user explicitly wants thoroughness here).

### Task B — explain the ~23% phantom-size difference
Each pipeline reconstructs the phantom at a slightly different physical size (z-scale
spread Faraz 1.16 … Lustig 0.94). Root cause = FOV / k-space-scaling / oversampling
convention per pipeline. This MUST be pinned before any quantitative extent/resolution
claim in a manuscript. Approach: for each pipeline trace the trajectory-units → grid
mapping (ours: `asap_recon.grid_to_radians`, `bart_traj=traj_rad·N/(2π)`; Lustig:
`build_acrtest.py` max-radius normalization; BART: its [-N/2,N/2] convention; Steve:
grid-index `g=k·MS/IS+MS/2`; Faraz: his `KSpaceCoor` 1/mm scaling). Compare the implied
mm/voxel. A point-source PSF test (Task A.5) gives the per-pipeline FOV directly.

### Task C — final 4D / temporal recon (the real next phase)
Architecture already decided (D3): bins enter as per-sample weights `W_b`; `recon()` API
unchanged. Needs: dynamic (multi-bin) data — current `recon_io` is single static bin;
diaphragm binning → bin-center weights (`bins.npy` exists). Add a temporal regularizer
(temporal-TV and/or low-rank across bins). **Undersampling returns here** (recompute DCF
per sampling mask — contract in `recon/AGENTS.md`). b-map stage (O4) becomes relevant.
BART `pics` has first-class temporal dims if we want the independent 4D check too.

---

## Hard constraints (carry forward)
- Run ours with `../.venv/bin/python` (arm64; finufft 2.5.1, sigpy 0.1.27). Lustig DCF
  step uses `lustig_oneshot/.venv_lustig` (torch). NEVER conda base (x86/Rosetta).
- NEVER git commit/push in repo ROOT. All git in `workspace/` only.
- Never modify Faraz's code, Steve's upstream code, or `run_cs.m` (byte-for-byte original).
- Read `recon/AGENTS.md` pitfalls before touching recon code — esp. the z-affine /
  eye-over-metric / BART-needs-ADMM entries.
- ALWAYS z-affine register before any cross-pipeline slice comparison (~23% size spread).

## Suggested skills
- `/code-review` on the new comparison scripts (Task A) — fresh eyes on the
  orientation/registration + BART invocation math before they back a manuscript.
- `/intent-layer-maintenance` after Task A if scripts move into a `compare/` child node.
- `/handoff update <note>` as Task A validation lands.

## Key paths
- Report (read first): `reference/Final_Report_CS_Comparison.md`
- Chronology + dead-ends: `reference/BART_Comparison.md`
- Code: `helpers/recon/` (+ `lustig_oneshot/`); node docs: their `AGENTS.md`
- Figures + metrics JSON: `data/v3_fov250/recon_io/` (and `bart/`, `lustig/` subdirs)
- BART binary: `~/bin/bart-src/bart`
- Memories: `eye_vs_metric` (scalar traps ×3), `slice_matching_zaffine` (z-affine gotcha)
