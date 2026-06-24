# Repo Architecture Plan — ASAP Recon + XeCS

**Date:** 2026-06-20  
**Status:** Planning — not yet executed

---

## The Model (Two-by-Two)

```
                   NON-CS                          CS
              ─────────────────────────────────────────────────
MAIN REPO  │  2026_ASAP_Recon                 2026_XeCS_Recon
           │  (Steve/Kento, read-only)        (CS operators, Hooman's code)
           │
WORKSPACE  │  ASAP workspace                  XeCS workspace
           │  Steve pipeline, Faraz,          BART, Lustig, CS comparison,
           │  Steve-vs-CS comparison          4D pipeline, CS sweeps
           │
DATA       │  Shared (already lives outside — /workspace/data/ is a convenience copy)
```

Each workspace can **import** from the other's main repo via `sys.path` / `pip install -e`.
Neither workspace needs to own the data — just point at the right paths via config.

---

## What Goes Where

### 2026_ASAP_Recon / workspace (trimmed — non-CS application layer)

**Primary question answered here: "Steve vs Faraz" and "CS vs non-CS"**

```
workspace/
├── pipeline/              ← Steve/Tyger orchestration (asap_run.py, post_process.py, param_gui.py)
├── helpers/               ← non-CS analysis tools
│   ├── dump_inputs*.py    ←   Steve data prep (tied to Steve's loaders)
│   ├── compare_baseline.py ←  FINUFFT vs Steve (non-CS baseline check)
│   ├── steve_kernel_numpy.py ← Steve GPU kernel validation (CPU numpy reimpl)
│   ├── faraz_montage.py   ←   Faraz figure tools
│   ├── faraz_zoom_check.py ←  Faraz zoom bug documentation
│   ├── convert_calib.py   ←   calib format converter
│   └── cs_montage.py      ←   our CS vs Steve vs Faraz (cross-repo comparison)
├── faraz/                 ← Faraz's MATLAB fork reference copy (was codes/2023_Faraz_Recon_HH/)
├── data/                  ← scanner data (gitignored)
├── outputs/               ← Steve/Tyger pipeline outputs
├── instructions/          ← Auto_Steve_Recon.md ✅, Auto_Faraz_Recon.md ❌
├── reference/             ← ASAP-specific docs:
│                          ←   Recon_Overview_Steve.md, Recon_Overview_Faraz.md
│                          ←   Recon_Comparison_StaticGas.md
│                          ←   Final_Report_CS_Comparison.md (cross-impl, belongs here)
│                          ←   BART_Comparison.md, Lustig_CS_*.md
│                          ←   Tyger_Setup.md
├── handoffs/
└── archive/               ← + codes/2025_CS/, codes/2025_Xe129_CS/ (archived)
```

---

### 2026_XeCS_Recon (new repo — CS operators library + CS workspace)

**CS library:** general FINUFFT operator + CS solvers, reusable on any xenon dataset.  
**CS workspace primary question:** "Which CS wins? (ours vs BART vs Lustig)"

```
2026_XeCS_Recon/
├── recon/                 ← CS operator library (was helpers/recon/)
│   ├── asap_recon.py      ←   FINUFFT spiral operator
│   ├── cs_recon.py        ←   static 3D CS (CG/FISTA/PDHG)
│   ├── cs_recon_4d.py     ←   4D temporal CS (circular TV, PDHG)
│   ├── binning.py         ←   soft respiratory binning
│   ├── surrogates.py      ←   FID/pneumo/diaphragm surrogates
│   ├── selftest.py        ←   operator self-test (synthetic)
│   └── selftest_4d.py     ←   4D operator self-test
├── pipeline/              ← 4D CS runnable pipelines
│   ├── cine_4d.py         ←   4D cine pipeline (load→surrogate→bins→solve→montage)
│   ├── export_4d.py       ←   NIfTI + MATLAB export from cine_joint.npy
│   ├── slice_video.py     ←   all-slice cine video (MP4+GIF)
│   ├── nav_movie.py       ←   navigator coronal cine with diaphragm overlay
│   ├── surrogate_compare.py ← respiratory surrogate QA
│   ├── diaphragm_check.py ←   diaphragm surrogate validation
│   └── kernel_check.py    ←   diaphragm detector 1D profile
└── workspace/             ← XeCS application layer
    ├── compare/           ←   CS comparison scripts
    │   ├── bart_compare.py ←    BART pics vs ours
    │   ├── bartio.py       ←    BART .cfl/.hdr I/O
    │   ├── lustig_compare.py ←  Lustig vs ours
    │   ├── wavelet_twoway.py ←  ours vs BART wavelet
    │   ├── tv_threeway.py  ←    ours vs BART vs Lustig TV
    │   ├── slice_matched_compare.py ← z-affine registered slices
    │   ├── z_register_compare.py ← z-affine registration
    │   ├── zreg_sixway_montage.py ← FINAL 6-way montage
    │   ├── resolution_sweep.py ← per-pipeline λ sweep at one slice
    │   ├── cg_tune.py      ←    CG λ × iterations sweep
    │   └── metrics_v2.py   ←    FIXED quality metrics (corner-ROI SNR, half-max, edge_sharp)
    ├── lustig/            ←   Lustig MATLAB one-shot (was helpers/lustig_oneshot/)
    ├── outputs/           ←   CS sweep outputs, 4D cine outputs (025JC etc.)
    ├── data/              ←   config.json pointing to ASAP recon_io dumps
    ├── instructions/      ←   Auto_CS_Static_Recon.md ❌
    │                      ←   Auto_CS_4D_Recon.md ❌
    │                      ←   Auto_BART_Recon.md ❌
    │                      ←   Auto_Lustig_Recon.md ❌
    ├── reference/         ←   CS_Implementation.md, Dynamic_4D_CS_Implementation.md
    │                      ←   4D_CS_Theory_Limitations.md, Physics_Notes.md
    │                      ←   Compressed_Sensing_Dynamic_Imaging*.md
    └── handoffs/
```

---

## Cross-Repo Calls

```
ASAP workspace helpers
  └─ cs_montage.py
       sys.path → 2026_XeCS_Recon/recon/        (imports asap_recon, cs_recon)
       data path → ASAP workspace/data/           (recon_io dumps, Steve output)

XeCS workspace compare
  └─ wavelet_twoway.py, tv_threeway.py, etc.
       imports XeCS recon/ directly (same repo)
       sys.path → 2026_ASAP_Recon/               (Steve's loaders if needed)
       data path → ASAP workspace/data/           (recon_io dumps via config.json)
```

Both workspaces read from the same `recon_io/` path. Neither owns the data.

---

## File Migration Table

| File (current location) | Destination |
|---|---|
| `helpers/recon/asap_recon.py` | XeCS `recon/` |
| `helpers/recon/cs_recon.py` | XeCS `recon/` |
| `helpers/recon/cs_recon_4d.py` | XeCS `recon/` |
| `helpers/recon/binning.py` | XeCS `recon/` |
| `helpers/recon/surrogates.py` | XeCS `recon/` |
| `helpers/recon/selftest*.py` | XeCS `recon/` |
| `helpers/recon/cine_4d.py` | XeCS `pipeline/` |
| `helpers/recon/export_4d.py` | XeCS `pipeline/` |
| `helpers/recon/slice_video.py` | XeCS `pipeline/` |
| `helpers/recon/nav_movie.py` | XeCS `pipeline/` |
| `helpers/recon/surrogate_compare.py` | XeCS `pipeline/` |
| `helpers/recon/diaphragm_check.py` | XeCS `pipeline/` |
| `helpers/recon/kernel_check.py` | XeCS `pipeline/` |
| `helpers/recon/bart_compare.py` | XeCS `workspace/compare/` |
| `helpers/recon/bartio.py` | XeCS `workspace/compare/` |
| `helpers/recon/lustig_compare.py` | XeCS `workspace/compare/` |
| `helpers/recon/wavelet_twoway.py` | XeCS `workspace/compare/` |
| `helpers/recon/tv_threeway.py` | XeCS `workspace/compare/` |
| `helpers/recon/slice_matched_compare.py` | XeCS `workspace/compare/` |
| `helpers/recon/z_register_compare.py` | XeCS `workspace/compare/` |
| `helpers/recon/zreg_sixway_montage.py` | XeCS `workspace/compare/` |
| `helpers/recon/resolution_sweep.py` | XeCS `workspace/compare/` |
| `helpers/recon/cg_tune.py` | XeCS `workspace/compare/` |
| `helpers/recon/metrics_v2.py` | XeCS `workspace/compare/` |
| `helpers/lustig_oneshot/` | XeCS `workspace/lustig/` |
| `helpers/recon/dump_inputs*.py` | ASAP workspace `helpers/` (stays) |
| `helpers/recon/compare_baseline.py` | ASAP workspace `helpers/` (stays) |
| `helpers/recon/steve_kernel_numpy.py` | ASAP workspace `helpers/` (stays) |
| `helpers/recon/faraz_*.py` | ASAP workspace `helpers/` (stays) |
| `helpers/recon/convert_calib.py` | ASAP workspace `helpers/` (stays) |
| `helpers/recon/cs_montage.py` | ASAP workspace `helpers/` (stays) |
| `codes/2023_Faraz_Recon_HH/` | ASAP workspace `faraz/` |
| `codes/2025_CS/` | ASAP workspace `archive/` |
| `codes/2025_Xe129_CS/` | ASAP workspace `archive/` |
| `codes/2025-09-24_ACR/` | ASAP workspace `archive/` |
| `reference/CS_Implementation.md` | XeCS `workspace/reference/` |
| `reference/Dynamic_4D_CS_Implementation.md` | XeCS `workspace/reference/` |
| `reference/4D_CS_Theory_Limitations.md` | XeCS `workspace/reference/` |
| `reference/Physics_Notes.md` | XeCS `workspace/reference/` |
| `reference/Compressed_Sensing_Dynamic_Imaging*.md` | XeCS `workspace/reference/` |
| `reference/Final_Report_CS_Comparison.md` | ASAP workspace `reference/` (stays) |
| `reference/BART_Comparison.md` | ASAP workspace `reference/` (stays) |
| `reference/Lustig_CS_*.md` | ASAP workspace `reference/` (stays) |
| `reference/Recon_Overview_Steve.md` | ASAP workspace `reference/` (stays) |
| `reference/Recon_Comparison_StaticGas.md` | ASAP workspace `reference/` (stays) |

---

## Instruction Docs Needed (6 implementations)

### ASAP workspace `instructions/`
| Doc | Entry point | Status |
|-----|-------------|--------|
| `Auto_Steve_Recon.md` | `pipeline/asap_run.py` → Tyger fork image | ✅ EXISTS |
| `Auto_Faraz_Recon.md` | `faraz/spiral_human_20240227.m` | ❌ needed |

### XeCS workspace `instructions/`
| Doc | Entry point | Status |
|-----|-------------|--------|
| `Auto_CS_Static_Recon.md` | `2026_XeCS_Recon/recon/cs_recon.py` | ❌ needed |
| `Auto_CS_4D_Recon.md` | `2026_XeCS_Recon/pipeline/cine_4d.py` | ❌ needed |
| `Auto_BART_Recon.md` | `workspace/compare/bart_compare.py` + `~/bin/bart-src/bart` | ❌ needed |
| `Auto_Lustig_Recon.md` | `workspace/lustig/run_lustig.py` | ❌ needed |

---

## Execution Order (when ready to execute)

1. `mkdir` + `git init` at `~/Hooman/Work/Codes/2026_XeCS_Recon/`
2. Create folder structure (recon/, pipeline/, workspace/{compare,lustig,outputs,data,instructions,reference,handoffs})
3. Copy XeCS files from ASAP workspace (copy not `git mv` — separate repos)
4. Commit XeCS initial state with CLAUDE.md + README
5. Remove moved files from ASAP workspace + commit
6. Internal ASAP moves: `codes/2023_Faraz_Recon_HH/` → `faraz/`, `codes/2025_*/` → `archive/`
7. Update sys.path references in comparison scripts
8. Move arm64 `.venv` to XeCS workspace (operators live there now)
9. Write `workspace/data/config.json` in XeCS pointing at ASAP recon_io path
10. Run `selftest.py` + `selftest_4d.py` from XeCS to confirm clean imports
11. Write instruction docs in both workspaces

---

## Hard Rules Going Forward

- ASAP workspace: only things tied to Steve's code or cross-impl comparison (Steve as root)
- XeCS workspace: only CS-internal work — operators, 4D pipeline, CS vs BART vs Lustig
- Never put data in either repo — config.json points to the real data location
- XeCS is the library; ASAP workspace is the consumer (import via sys.path or pip install -e)
