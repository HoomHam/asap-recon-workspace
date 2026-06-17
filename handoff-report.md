# ASAP Recon — Workspace Handoff

**Date:** 2026-06-17 · **Branch context:** workspace git (`HoomHam/asap-recon-workspace`). Latest push `6174438`.

> **Two active workstreams:**
> 1. **3D+t dynamic CS** — implementation now FINAL (this session). λ_t sweep started, more to come.
> 2. **Static CS comparison** (Steve/Faraz/BART/Lustig) — super-validation deferred, still open.
> Pick up whichever is next.

---

## ⭐ 3D+t Dynamic CS — IMPLEMENTATION FINAL

The wavelet_xyz + TV_t pipeline is **locked**. Do not re-architect it without reason.
Read `reference/Dynamic_4D_CS_Implementation.md` for full as-built details.

### What the implementation is

- **B-bin respiratory cine** (B=16 default) from 3-min free-breathing Xe-129 acquisition
- **Objective:** per-bin spatial wavelet (db4) + circular temporal TV across bins, PDHG solver,
  DCF-preconditioned, single coil (nch=1 for 027JC/025JC → no SENSE needed)
- **Three surrogates:** signal / pneumotach / diaphragm-via-CS-nav
- **10 modules** in `helpers/recon/`: `dump_inputs_dyn`, `binning`, `cs_recon_4d`, `selftest_4d`,
  `surrogates`, `cine_4d`, `diaphragm_check`, `kernel_check`, `nav_movie`, `surrogate_compare`
- **New tools added this session:** `slice_video.py` (all-slice grid video, 3 orientations),
  `export_4d.py` (NIfTI + MATLAB export from cine_joint.npy)
- `selftest_4d` = ALL PASS

### Hard-won diaphragm conventions (do NOT re-break)

- SI axis = **axis2** (corr 0.94 vs signal); `win_ilv=20`, `smooth_win=5`, no median filter
- Surrogate curve / binning → **auto corr-based selection** (`prefer="auto"` default):
  on 025JC → lo edge wins (corr 0.83 vs hi 0.63), period 3.52s, 276/349 valid windows
- `prefer="hi"` = anatomical diaphragm dome = the **nav_movie cyan-dashed DISPLAY line** only —
  clips out of FOV at deep inspiration on 025JC (40% valid, 2× harmonic period) → NOT for binning
- On **027JC**: nav quality too low for diaphragm surrogate (82/363 valid, drift curve) → use signal
- Orientation display: sagittal=axis0 no-rot, coronal=axis1 90°CW, axial=axis2 90°CCW+fliplr

### Output structure (locked)

**`outputs/025JC/`** — canonical tracked outputs, committed to git:

| File | Contents |
|------|----------|
| `diaphragm_check.png` | surrogate QA (auto/lo, period 3.52s) |
| `nav_montage.png` | nav grid (coronal, N=80) |
| `nav_coronal.gif/.mp4` | coronal nav movie with diaphragm overlay |
| `axial_joint.gif/.mp4` | axis2, 90°CCW, fliplr |
| `coronal_joint.gif/.mp4` | axis1, 90°CW |
| `sagittal_joint.gif/.mp4` | axis0, no rotation |
| `cine_4d_025JC_diaphragm.nii.gz` | local only (55MB, gitignored) |
| `cine_4d_025JC_diaphragm.mat` | local only (55MB, gitignored) |

`outputs/archive/` — stale axis-numbered files, gitignored.

**Data folder note:** `data/v3_dyn_025JC/recon_io_dyn/` still contains some non-raw result
outputs (sweep dirs `sweep_lt*/`, nav_movie dirs, diaphragm curves). These are not yet cleaned
up — the full data dir cleanup is deferred (raw acq/traj files must stay, result files should
eventually move to `outputs/`). `data/` is gitignored so nothing there is at risk.

---

## ⭐ λ_t Sweep — FIRST PASS DONE, MORE TO COME

### What was run

`lam_t ∈ {0.003, 0.01, 0.02, 0.05, 0.1, 0.2}` on 025JC diaphragm cine (B=16, lam_s=0.01).
Baseline computed once, 6 joint sweeps from shared baseline. Coronal slice_video generated
for each for visual comparison.

Results at `data/v3_dyn_025JC/recon_io_dyn/sweep_lt{val}/` (gitignored).
Copies for viewing: `outputs/025JC_sweep_lt/coronal_lt{val}_joint.mp4` + `montage_lt{val}.png`.

| lam_t | MP4 size | Interpretation |
|-------|----------|---------------|
| 0.003 | 800 KB | near-baseline; streaky |
| 0.01 | 661 KB | mild smoothing |
| 0.02 | 563 KB | moderate |
| 0.05 | 422 KB | current working default |
| 0.1 | 261 KB | heavy smoothing |
| 0.2 | 242 KB | likely over-smoothed |

**Hooman's visual pick is PENDING** — view the coronal MP4s and choose the λ_t that
preserves motion between bins while suppressing per-bin streaks.

### Plan going forward

More sweeps will run to converge on optimal recon. The eventual comparison plan:

- **Intra-4D:** sweep lam_t, lam_s, B (bins), surrogate choice (signal vs diaphragm)
- **Inter-implementation:** compare Steve / Faraz / BART / ours-static / ours-4D once parameters
  are locked. This will require z-affine registration (same ~23% size spread likely applies).
- Goal: pick best-optimized 4D recon then formally compare to baselines. Not yet started.

---

## ⭐ CS Theory Discussion (this session)

A detailed theoretical discussion was had and **written up** in:
`reference/4D_CS_Theory_Limitations.md`

Key points (read the doc for full depth):
1. **Spatial (wavelet) vs temporal (TV) are handled separately** — orthogonal domains,
   independent knobs, correct physics.
2. **Temporal wavelet fails** — 16 bins too few, respiratory motion is smooth (not wavelet-sparse).
   TV is the right prior for slowly-varying signals.
3. **4D spatial treatment doesn't win** — anisotropic motion, incommensurable units, 16<<100.
4. **Two known limitations of current implementation:**
   - Temporal TV penalizes real respiratory motion at fixed voxels (motion-at-fixed-voxel problem)
   - 16 bins is coarse; intra-bin motion already blurred before CS
5. **We accept these for now** (first-pass demonstration). The doc lists 4 concrete tests needed
   before any publication claim (motion preservation curve, bin-count tradeoff, L+S comparison,
   motion-compensated TV).
6. **Next principled step:** L+S (low-rank + sparse) — handles motion without registration.

---

## Comparison phase (PRIOR session) — FINALIZED, do not redo

> **Read:** `reference/Final_Report_CS_Comparison.md` (narrative, per-script I/O §7,
> figure catalog §8, reproduce + BART build §9, open threads §10).

1. **`metrics_v2.py`** — fixed quality metrics (corner-ROI SNR, half-max extent, `edge_sharp`).
2. **BART built from source** at `~/bin/bart-src/bart`. BART **wavelet MUST use ADMM (`-m`)**.
3. **Lustig TV λ-sweep** via one-shot; NCG stalls (obj flat).
4. **Slice-matching solved** via z-affine registration (`z_register_compare.py`); ~23% size spread.
5. **FINAL verdict:** all pipelines resolve ACR grid; **ours as sharp as BART, sharper than Lustig.**
6. Final report + Intent Layer written/updated.

### RETRACTED (don't resurrect)
- "ours-wavelet wins on lfCV" — blur erasing resolution.
- "BART has moiré artifact" — it's the real ACR insert; BART resolves it.
- "DCF removes the grid" — disproven.

---

## Carried over / NOT finalized (open loops)

| # | Item | Where it stands |
|---|------|-----------------|
| O1 | **λ verdict (ours wavelet t0.003 vs t0.01)** | still Hooman's eye; RESOLUTION call at insert slice |
| O2 | **BART per-sample DCF** | `pics -p` did not wire noncart density |
| O3 | **Lustig wavelet** | not possible — `FWT2_PO` 2D on 3D volume |
| O4 | **b-map stage** (`calcb` port) | 027JC/025JC nch=1 so 4D didn't need it; unported for multi-coil |
| O5 | Steve's rebuttal to `for_steve.md`; tell Faraz ×1.205 zoom bug | pending |
| O6 | `recon/*.py` token bloat | consider splitting `dyn4d/` + `compare/` child nodes |
| O7 | **Hooman's λ_t visual pick** | view `outputs/025JC_sweep_lt/` coronal MP4s |
| O8 | **027JC full pipeline** | nav quality too low for diaphragm; run signal surrogate cine instead |
| O9 | **data/ cleanup** | sweep dirs + nav_movie dirs still in `data/v3_dyn_025JC/recon_io_dyn/` |

---

## DEFERRED comparison tasks

### Task A — SUPER-validate every comparison script (presentation/manuscript-grade)
1. **Operator:** `selftest.py` passes (adjoint `<Ax,y>=<x,Aᴴy>`).
2. **metrics_v2 unit-check:** synthetic phantoms with KNOWN snr/extent/edge.
3. **bartio round-trip:** `writecfl`→`readcfl` identity.
4. **Orientation/registration:** `best_orient_full`/`fit_z_affine` on KNOWN transform → recovers inverse.
5. **BART invocation:** `bart_traj=traj_rad·N/(2π)` via PSF; `-m` ADMM convergence.
6. **Lustig reproduction:** diff `run_cs_sweep.m` vs `run_cs.m`.
7. **Cross-pipeline sanity:** every recon hits ACR extent ~190/190/148 mm.
Deliverable: `reference/Validation_Report.md`.

### Task B — explain the ~23% phantom-size difference
Trace each pipeline's trajectory-units→grid mapping. PSF test (Task A.5) gives per-pipeline FOV directly.

---

## Hard constraints (carry forward)
- Run ours with `workspace/helpers/.venv/bin/python` (arm64; finufft 2.5.1, sigpy 0.1.27).
- NEVER git commit/push in repo ROOT. All git in `workspace/` only.
- Never modify Faraz's code, Steve's upstream code, or `run_cs.m`.
- Read `recon/AGENTS.md` pitfalls + `Dynamic_4D_CS_Implementation.md` before touching recon code.
- ALWAYS z-affine register before any cross-pipeline slice comparison.
- Trajectory = `fa_spiral_dyn_fancy_v3_20240130_{gp,dp}.npy`; display flips are plot-only.
- 025JC source `.dat` at `data/xe/human/2024-11-13/025JC/` (needed if dump must be regenerated).

## Suggested skills
- `/handoff update <note>` when λ_t pick is made or 027JC signal cine runs.
- `/code-review` on 4D modules before manuscript.
- `/graphify` — stale (pre-4D).
- `/intent-layer-maintenance` if `helpers/recon/` is split.

## Key paths
- 4D AS-BUILT: `reference/Dynamic_4D_CS_Implementation.md`
- 4D theory + limitations: `reference/4D_CS_Theory_Limitations.md` ← NEW, read for CS discussion
- Comparison report: `reference/Final_Report_CS_Comparison.md`
- Code: `helpers/recon/` (+ `lustig_oneshot/`)
- 025JC dump: `data/v3_dyn_025JC/recon_io_dyn/`
- 025JC λ_t sweep outputs: `outputs/025JC_sweep_lt/`
- 025JC canonical outputs: `outputs/025JC/`
- 027JC dump: `data/v3_dyn/recon_io_dyn/`
- Comparison figures: `data/v3_fov250/recon_io/`
- BART binary: `~/bin/bart-src/bart`
- Memories: `eye_vs_metric`, `slice_matching_zaffine`
