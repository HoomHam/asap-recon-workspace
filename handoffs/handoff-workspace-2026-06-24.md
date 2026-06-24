# ASAP Recon — Workspace Handoff

**Date:** 2026-06-24 · **Branch context:** workspace git (`HoomHam/asap-recon-workspace`) + fork `HoomHam/asap_recon` branch `diaphragm-recon` (root-side). Latest workspace push `dea9ee9`, fork head `db80f16`.

> **Three workstreams + one structural decision:**
> 1. **Root / Tyger GPU recon (Steve's side)** — DIAPHRAGM navigator ported, figures reworked. Fork `diaphragm-recon`, image public at ghcr.io.
> 2. **3D+t dynamic CS** — implementation FINAL, λ_t sweep first pass done. More sweeps to come; λ_t visual pick pending.
> 3. **Static CS comparison** (Steve/Faraz/BART/Lustig) — super-validation deferred, still open.
> 4. **⭐ NEW (2026-06-24): Repo restructure** — workspace will split into two repos. See `workspace/Repo_Architecture.md` and the section below.
> Pick up whichever is next.

---

## ⭐ REPO ARCHITECTURE DECISION (2026-06-24)

**Full plan:** `workspace/Repo_Architecture.md` — read this first before any restructuring.

### The Model

Hooman runs a two-repo pattern per project (main repo + workspace). The CS implementation has outgrown the workspace concept and needs its own home.

```
                   NON-CS                        CS
MAIN REPO    2026_ASAP_Recon (read-only)   2026_XeCS_Recon (NEW — git init needed)
WORKSPACE    ASAP workspace (trimmed)       XeCS workspace (NEW)
DATA         Shared — both repos point at the same recon_io/ via config.json
```

### What goes where (summary)

**Stays in ASAP workspace:**
- `pipeline/` — Steve/Tyger orchestration (asap_run.py, post_process.py)
- `helpers/` — ASAP-specific: dump_inputs, compare_baseline, steve_kernel_numpy, faraz_*, cs_montage
- `faraz/` — Faraz's MATLAB fork reference copy (move from `codes/2023_Faraz_Recon_HH/`)
- `reference/` — Final_Report_CS_Comparison, BART_Comparison, Lustig_CS_*, Recon_Overview_Steve, StaticGas

**Moves to 2026_XeCS_Recon (main repo):**
- `helpers/recon/asap_recon.py`, `cs_recon.py`, `cs_recon_4d.py`, `binning.py`, `surrogates.py`, `selftest*.py`
- `helpers/recon/cine_4d.py`, `export_4d.py`, `slice_video.py`, `nav_movie.py`, `diaphragm_check.py`, `kernel_check.py`, `surrogate_compare.py`

**Moves to 2026_XeCS_Recon / workspace:**
- `helpers/recon/bart_compare.py`, `bartio.py`, `wavelet_twoway.py`, `tv_threeway.py`, `zreg_sixway_montage.py`, `slice_matched_compare.py`, `z_register_compare.py`, `metrics_v2.py`, `resolution_sweep.py`, `cg_tune.py`, `lustig_compare.py`
- `helpers/lustig_oneshot/` → XeCS `workspace/lustig/`
- `reference/CS_Implementation.md`, `Dynamic_4D_CS_Implementation.md`, `4D_CS_Theory_Limitations.md`, `Physics_Notes.md`, `Compressed_Sensing_Dynamic_Imaging*.md`
- `outputs/025JC/`, `outputs/025JC_sweep_lt/` (CS outputs)

**Archived within ASAP workspace:**
- `codes/2025_CS/`, `codes/2025_Xe129_CS/`, `codes/2025-09-24_ACR/` → `archive/`

### Cross-repo import pattern
ASAP workspace scripts that need CS operators: `sys.path` → `2026_XeCS_Recon/recon/`
XeCS compare scripts that need data: `config.json` in XeCS `workspace/data/` pointing at ASAP `recon_io/` path

### Status
**✅ EXECUTED 2026-06-24.** `2026_XeCS_Recon` created (`git init`, commit `c20b0d4`,
45 files); ASAP workspace trimmed (commit `8a07d2d`). Self-tests PASS.

- **Cross-repo import = `.pth`, not `sys.path` edits.** `xecs_recon.pth` in both venvs
  → `2026_XeCS_Recon/recon`. ASAP stayers (`compare_baseline`, `cs_montage`,
  `faraz_*`) resolve operators from XeCS; zero file edits.
- **Nothing deleted** (Hooman's rule): the CS originals + `lustig_oneshot/` moved to
  `helpers/_delete/` (byte-verified vs XeCS first), not `git rm`. Remove later.
- **Deviations:** `cg_tune.py` → XeCS `recon/` (not `compare/` — it's in the
  `cs_recon` import chain, else `selftest_4d` breaks); CS reference docs kept in BOTH
  repos (cross-linked from ASAP docs); `outputs/025JC*` NOT moved (stayed in ASAP).
- **Detail:** `2026_XeCS_Recon/workspace/handoffs/XeCS_Handoff_2026-06-24.md`,
  root `start_here_to_decouple.md` → "Execution Result".
- **Open loops below: paths still pre-restructure** for CS items — CS work now lives
  in `2026_XeCS_Recon` (operators `recon/`, pipeline `pipeline/`, compare
  `workspace/compare/`). Faraz fork now at `workspace/faraz/`.
- **Phase 14 TODO:** instruction docs `Auto_CS_Static/4D/BART/Lustig_Recon.md` (XeCS)
  + `Auto_Faraz_Recon.md` (ASAP); recreate XeCS `.venv_lustig`.

---

## ⭐ ROOT-SIDE: 2026-06-19 — DIAPHRAGM Navigator on Tyger

**Scope:** Fork `HoomHam/asap_recon`, branch `diaphragm-recon` (head `db80f16`, pushed).
All edits live on the fork only; main MEDCAP repo untouched.

### Fork & CI/CD

- **Fork:** `github.com/HoomHam/asap_recon`, remote `hooman`, branch `diaphragm-recon`
- **CI workflow** (`.github/workflows/build-image.yml`): builds amd64 container, pushes to `ghcr.io/hoomham/xe-tyger-recon`
  - Mutable tag: `diaphragm-recon` (rotates with each push)
  - Immutable tag: `:<sha>` (safe for caching in Tyger)
  - Package is **public** — Tyger can pull
- **mrd-fork pinned:** `requirements.txt` → commit `b6b6d18` (matches local writer schema; otherwise `RuntimeError: Invalid schema`)
- **Current image:** `ghcr.io/hoomham/xe-tyger-recon:db80f16ecc5f53a2af6a44f834c571d619ff09a5`
  - Use the `:<sha>` tag in Tyger codespec (mutable tag cached locally)

### tyger_recon.py Changes

- **Ported Steve's DIAPHRAGM navigator** from `main.py:calcLVcb`
  - Low-res navigator loop → diaphragm z via 25–75% dropoff parabola → Savgol smooth → `ilvbin[DIAPHRAGM]`
  - Exports navigator arrays into output.mrd: `nav_coronal`, `nav_diaphragm_z`, `nav_time`, `nav_volume`, `nav_ilvtime`, `nav_volmeastime`
  - Also written to `recon.mat`

- **Diaphragm-vs-apex tracking fix**
  - Root cause: `dyn_usimg_recon` stores GPDYN z-flipped; for Siemens-via-MRD data, diaphragm sits at low z
  - Steve's "walk from high-z end" edge-finder locked **apex** (less-mobile, tracked range ~5 vs ~10 for diaphragm)
  - **Fix:** flip navigator GPDYN z back (`getimg(GPDYN)[::-1]`) so diaphragm sits at high z and finder tracks it
  - **Verified:** line now on inferior lung edge, moves with breathing ✓

### workspace/pipeline/post_process.py — Figure Rework

- **Slice montages** (axial/coronal/sagittal) now **per-bin videos**, 10-wide grid:
  - Lung background trimmed, fixed color axis (vmin=EE bin, vmax=EI bin)
  - Orientation: axial flipud, coronal slice-order reversed, sagittal rot90 ccw
  - Last row dropped if < 4 tiles
  - `gas_montage.png` disabled
  
- **navigator.gif:** coronal projection + tracked diaphragm line (no display flip; nav_coronal now exported apex-up from tyger_recon)

- **resp_traces.png (NEW):** real navigator waveform
  - Green navigator volume (tracked): solid, breath-by-breath, follows FID envelope ✓
  - Green diaphragm-z measurements: **noisy + spreads out in last ~15%** of acquisition (signal decay region) — cloud of points after t≈0.85
  - Image-derived per-bin cyan: demoted to faint (intentional)
  - **FUTURE:** filter/reject nav-z scatter in decay region, sign-align with pneumotach (up=inspiration consistently)

### Latest Good Output

**Location:** `workspace/outputs/25JC/20260619-053136/`
- Contains all figures: videos, montages, traces, navigator gif
- Image tag in run: `db80f16`

### Instruction Doc Pattern (NEW)

Wrote `instructions/Auto_Steve_Recon.md` — **model for future docs** in this project:
- Build map (file roles, data flow)
- Call graph (entry points, GPU dispatch)
- MRD provenance (how arrays get into output)
- Bash run guide (exact commands, var substitution)
- CS-comparison lookup playbook

**Recommendation:** instruction docs are the primary fast-lookup tool; `/graphify` only for broad exploratory cross-repo sweeps.
We want the same style for: BART comparison, Lustig CS, 4D/static CS.

---

## ⭐ WORKSPACE-SIDE: 3D+t Dynamic CS — IMPLEMENTATION FINAL

The wavelet_xyz + TV_t pipeline is **locked**. Do not re-architect it without reason.
Read `reference/Dynamic_4D_CS_Implementation.md` for full as-built details.

### What the Implementation Is

- **B-bin respiratory cine** (B=16 default) from 3-min free-breathing Xe-129 acquisition
- **Objective:** per-bin spatial wavelet (db4) + circular temporal TV across bins, PDHG solver, DCF-preconditioned, single coil (nch=1 for 027JC/025JC → no SENSE needed)
- **Three surrogates:** signal / pneumotach / diaphragm-via-CS-nav
- **10 modules** in `helpers/recon/`: `dump_inputs_dyn`, `binning`, `cs_recon_4d`, `selftest_4d`, `surrogates`, `cine_4d`, `diaphragm_check`, `kernel_check`, `nav_movie`, `surrogate_compare`
- **New tools added this session:** `slice_video.py` (all-slice grid video, 3 orientations), `export_4d.py` (NIfTI + MATLAB export from cine_joint.npy)
- `selftest_4d` = ALL PASS

### Hard-Won Diaphragm Conventions (do NOT re-break)

- SI axis = **axis2** (corr 0.94 vs signal); `win_ilv=20`, `smooth_win=5`, no median filter
- Surrogate curve / binning → **auto corr-based selection** (`prefer="auto"` default):
  on 025JC → lo edge wins (corr 0.83 vs hi 0.63), period 3.52s, 276/349 valid windows
- `prefer="hi"` = anatomical diaphragm dome = the **nav_movie cyan-dashed DISPLAY line** only —
  clips out of FOV at deep inspiration on 025JC (40% valid, 2× harmonic period) → NOT for binning
- On **027JC**: nav quality too low for diaphragm surrogate (82/363 valid, drift curve) → use signal
- Orientation display: sagittal=axis0 no-rot, coronal=axis1 90°CW, axial=axis2 90°CCW+fliplr

### Output Structure (Locked)

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

**Data folder note:** `data/v3_dyn_025JC/recon_io_dyn/` still contains some non-raw result outputs (sweep dirs `sweep_lt*/`, nav_movie dirs, diaphragm curves). These are not yet cleaned up — full data dir cleanup deferred (raw acq/traj must stay, result files should eventually move to `outputs/`). `data/` is gitignored.

---

## ⭐ λ_t Sweep — FIRST PASS DONE, MORE TO COME

### What Was Run

`lam_t ∈ {0.003, 0.01, 0.02, 0.05, 0.1, 0.2}` on 025JC diaphragm cine (B=16, lam_s=0.01).
Baseline computed once, 6 joint sweeps from shared baseline. Coronal slice_video generated for each for visual comparison.

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

**Hooman's visual pick is PENDING** — view the coronal MP4s and choose the λ_t that preserves motion between bins while suppressing per-bin streaks.

### Plan Going Forward

More sweeps will run to converge on optimal recon. The eventual comparison plan:

- **Intra-4D:** sweep lam_t, lam_s, B (bins), surrogate choice (signal vs diaphragm)
- **Inter-implementation:** compare Steve / Faraz / BART / ours-static / ours-4D once parameters are locked. This will require z-affine registration (same ~23% size spread likely applies).
- Goal: pick best-optimized 4D recon then formally compare to baselines. Not yet started.

---

## ⭐ CS Theory Discussion (this session)

Detailed theoretical discussion written up in: `reference/4D_CS_Theory_Limitations.md`

Key points (read the doc for full depth):
1. **Spatial (wavelet) vs temporal (TV) are handled separately** — orthogonal domains, independent knobs, correct physics.
2. **Temporal wavelet fails** — 16 bins too few, respiratory motion is smooth (not wavelet-sparse). TV is the right prior for slowly-varying signals.
3. **4D spatial treatment doesn't win** — anisotropic motion, incommensurable units, 16<<100.
4. **Two known limitations of current implementation:**
   - Temporal TV penalizes real respiratory motion at fixed voxels (motion-at-fixed-voxel problem)
   - 16 bins is coarse; intra-bin motion already blurred before CS
5. **We accept these for now** (first-pass demonstration). The doc lists 4 concrete tests needed before any publication claim (motion preservation curve, bin-count tradeoff, L+S comparison, motion-compensated TV).
6. **Next principled step:** L+S (low-rank + sparse) — handles motion without registration.

---

## ⚠️ CRITICAL: B0 Correction Check (CS Implementation)

**DO NOT LOSE THIS — Hooman will investigate during CS phase.**

- Steve's DIAPHRAGM navigator runs `g_res.calcb(...)` at **low res** (`g.MS = g.IS + 4`, e.g. MS=104) to build the b-matrix for per-usimg navigator recon, then recomputes b at **full res** before the binned recon. See `tyger_recon.py:_diaphragm_navigator` and `results.py:calcb`.
- The b-matrix carries the **B0 / off-resonance + coil phase** correction (`self.b[ich,:,:,:]`, applied in `dyn_usimg_recon` / `dyn_recon`).
- **TODO when continuing CS** (`workspace/helpers/recon/cs_recon.py`): confirm the CS forward/adjoint operator uses the **correct-resolution B0/b** — i.e. did we ever apply a low-res B0 where a full-res one was needed (or skip B0 entirely)? This was flagged off the navigator work; not yet investigated.

---

## Comparison Phase (PRIOR session) — FINALIZED, do not redo

> **Read:** `reference/Final_Report_CS_Comparison.md` (narrative, per-script I/O §7, figure catalog §8, reproduce + BART build §9, open threads §10).

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

## Carried Over / NOT Finalized (Open Loops)

> Paths below are **pre-restructure**. After `2026_XeCS_Recon` is created, update accordingly.
> "Lands In" column = which repo owns this work after restructure.

| # | Item | Where It Stands | Lands In |
|---|------|-----------------|----------|
| R1 | **PR: `diaphragm-recon → MEDCAP:dev`** | Optional — lets Kento pull the navigator fix upstream | Root fork |
| R2 | **Verify 023LL + 016PG `resp_traces.png`** | 016PG EE/EI bins flipped vs 023LL — confirm diaphragm tracking looks correct | ASAP workspace |
| O1 | **λ verdict (ours wavelet t0.003 vs t0.01)** | Hooman's eye; RESOLUTION call at insert slice | XeCS workspace |
| O2 | **BART per-sample DCF** | `pics -p` did not wire noncart density | XeCS workspace |
| O3 | **Lustig wavelet** | Not possible — `FWT2_PO` 2D on 3D volume; TV is only fair 3-way | XeCS workspace |
| O4 | **b-map stage** (`calcb` port) | 027JC/025JC nch=1 so 4D didn't need it; unported for multi-coil | XeCS main repo |
| O5 | **Steve rebuttal + Faraz zoom bug** | `for_steve.md` rebuttal pending; Faraz doesn't know ×1.205 zoom bug | ASAP workspace |
| O6 | ~~`recon/*.py` token bloat~~ | **Superseded by repo restructure plan** — `Repo_Architecture.md` handles the split | — |
| O7 | **Hooman's λ_t visual pick** | View `outputs/025JC_sweep_lt/` coronal MP4s, pick winner | XeCS workspace |
| O8 | **027JC full pipeline** | Nav quality too low for diaphragm; run signal surrogate cine instead | XeCS workspace |
| O9 | **data/ cleanup** | Sweep dirs + nav_movie dirs in `data/v3_dyn_025JC/recon_io_dyn/` | XeCS workspace |
| O10 | **resp_traces.png noise cleanup** | Filter/reject nav-z scatter in decay region; sign-align with pneumotach | ASAP workspace |
| ⚠️ | **B0 correction check (CS)** | Verify CS uses correct-resolution b-matrix (not low-res from navigator loop) — see `helpers/recon/cs_recon.py` vs `results.py:calcb` | XeCS main repo |

### Deferred Comparison Tasks (post-restructure → XeCS workspace)

**Task A — SUPER-validate every comparison script** (presentation/manuscript-grade)
1. Operator: `selftest.py` passes adjoint `<Ax,y>=<x,Aᴴy>`
2. `metrics_v2` unit-check: synthetic phantoms with KNOWN snr/extent/edge
3. `bartio` round-trip: `writecfl`→`readcfl` identity
4. Orientation/registration: `best_orient_full`/`fit_z_affine` on KNOWN transform → recovers inverse
5. BART invocation: `bart_traj=traj_rad·N/(2π)` via PSF; `-m` ADMM convergence
6. Lustig reproduction: diff `run_cs_sweep.m` vs `run_cs.m`
7. Cross-pipeline sanity: every recon hits ACR extent ~190/190/148 mm
Deliverable: `reference/Validation_Report.md` (in XeCS workspace)

**Task B — Explain ~23% phantom-size difference**
Trace each pipeline's trajectory-units→grid mapping. PSF test (Task A.5) gives per-pipeline FOV directly.

---

## Hard Constraints (Carry Forward)

- **Pre-restructure:** run CS with `workspace/helpers/.venv/bin/python` (arm64; finufft 2.5.1, sigpy 0.1.27). After restructure: venv moves to XeCS repo.
- NEVER git commit/push in repo ROOT (`2026_ASAP_Recon/`). All git in `workspace/` or XeCS repo only. Fork branch (`HoomHam/asap_recon`) is a separate repo — OK to push there.
- Never modify Faraz's code, Steve's upstream code, or `run_cs.m`.
- Read `helpers/recon/AGENTS.md` pitfalls + `Dynamic_4D_CS_Implementation.md` before touching CS recon code.
- ALWAYS z-affine register before any cross-pipeline slice comparison.
- Trajectory = `fa_spiral_dyn_fancy_v3_20240130_{gp,dp}.npy`; display flips are plot-only.
- 025JC source `.dat` at `data/xe/human/2024-11-13/025JC/` (needed if dump must be regenerated).

## Suggested Skills

- `/handoff update <note>` when λ_t pick made, restructure executed, 027JC signal cine runs, or resp_traces cleaned.
- `/code-review` on 4D modules before manuscript.
- `/graphify` — stale (pre-4D, pre-restructure).
- `/intent-layer-maintenance` after restructure executes (AGENTS.md files need new paths).

## Key Paths (pre-restructure — update after 2026_XeCS_Recon created)

### Root-side (fork)
- Fork branch: `HoomHam/asap_recon` `diaphragm-recon` (head `db80f16`)
- Tyger entrypoint: `tyger_recon.py` (navigator, b-matrix, export)
- Post-process: `workspace/pipeline/post_process.py` (figures, videos)
- Instruction doc: `workspace/instructions/Auto_Steve_Recon.md` ✅
- Latest Steve output: `workspace/outputs/25JC/20260619-053136/`
- Latest human subjects: `workspace/outputs/023LL/20260620-041321/`, `workspace/outputs/016PG/20260620-042143/`
- Image: `ghcr.io/hoomham/xe-tyger-recon:db80f16ecc5f53a2af6a44f834c571d619ff09a5`

### Workspace-side (CS — pre-restructure paths)
- **Repo architecture plan:** `workspace/Repo_Architecture.md` ← READ THIS FIRST
- 4D AS-BUILT: `workspace/reference/Dynamic_4D_CS_Implementation.md`
- 4D theory + limitations: `workspace/reference/4D_CS_Theory_Limitations.md`
- Comparison report: `workspace/reference/Final_Report_CS_Comparison.md`
- CS code: `workspace/helpers/recon/` (+ `helpers/lustig_oneshot/`)
- 025JC dump: `workspace/data/v3_dyn_025JC/recon_io_dyn/`
- 025JC λ_t sweep outputs: `workspace/outputs/025JC_sweep_lt/`
- 025JC canonical outputs: `workspace/outputs/025JC/`
- 027JC dump: `workspace/data/v3_dyn/recon_io_dyn/`
- Comparison figures: `workspace/data/v3_fov250/recon_io/`
- BART binary: `~/bin/bart-src/bart`
- Memories: `eye_vs_metric`, `slice_matching_zaffine`, `workspace_system`
