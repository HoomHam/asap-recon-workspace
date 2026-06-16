# ASAP Recon — Session Handoff
**Date:** 2026-06-08  
**Next focus:** Study major differences between Steve's (Python/CUDA) and Faraz's (MATLAB) implementations

---

## What Was Built This Session

### 1. Repository Setup

- Cloned `https://github.com/MEDCAP/asap_recon` → `/Users/hoomham/Hooman/Work/Codes/2026_ASAP_Recon/`
- Main repo is **read-only mirror** — git ops forbidden here (enforced in `CLAUDE.md` and `.git/info/exclude`)
- `workspace/` folder created with **its own independent git repo** (`workspace/.git`)

### 2. Folder Structure

```
2026_ASAP_Recon/               ← upstream MEDCAP mirror (NO git ops)
├── CLAUDE.md                  ← root Intent Layer node (git-excluded)
├── main.py / raw.py / recon.py / results.py / gtypes.py
├── asap/asap.c                ← Steve's C kernel
└── workspace/                 ← Hooman's workspace (own git, commit freely)
    ├── CLAUDE.md              ← workspace Intent Layer root
    ├── handoffs/              ← session archives
    ├── reference/
    │   ├── Recon_Overview.md      ← Steve's code map
    │   ├── Recon_Comparison.md    ← 10-point Steve vs Faraz diff
    │   └── Obsidian_Index.md      ← index of Obsidian notes
    ├── archive/
    ├── helpers/
    └── codes/
        ├── AGENTS.md              ← codes/ Intent Layer node
        ├── kasap.c                ← Kento's cloud C kernel
        └── 2023_Faraz_Recon_HH/   ← Hooman's MATLAB fork (active)
            ├── AGENTS.md          ← MATLAB fork Intent Layer node
            └── [all Faraz .m files]
```

### 3. External Paths (read-only references)

| Resource | Path |
|----------|------|
| Faraz's original MATLAB recon | `/Users/hoomham/Hooman/Work/Codes/2023_Faraz_Recon/` |
| Obsidian ASAP notes vault | `~/Library/Mobile Documents/iCloud~md~obsidian/Documents/Action/MRI/ASAP Recon/` |

### 4. Intent Layer (all nodes pass validation)

4 nodes created, all `PASSED with warnings` (no errors):

| Node | File | Purpose |
|------|------|---------|
| Main root | `CLAUDE.md` | Project identity, git rules, code map, pitfalls |
| Workspace root | `workspace/CLAUDE.md` | Session startup, reference docs, entry points |
| codes/ index | `workspace/codes/AGENTS.md` | What lives in codes/, contracts |
| MATLAB fork | `workspace/codes/2023_Faraz_Recon_HH/AGENTS.md` | Full file map, key params, pitfalls |

### 5. Reference Documents Created

- `workspace/reference/Recon_Overview.md` — Steve's data flow, gridding params, GPU dispatch
- `workspace/reference/Recon_Comparison.md` — compressed 10-point Steve vs Faraz diff
- `workspace/reference/Obsidian_Index.md` — index of 14 Obsidian notes with when-to-read guidance

---

## Project Context

**ASAP Recon** = Accelerated Spiral Acquisition and Processing for hyperpolarized Xe-129 lung MRI (non-Cartesian spiral k-space reconstruction).

| Role | Person |
|------|--------|
| Original author | Steve Kadlecek — Python/CUDA (`recon.py`, `raw.py`, `results.py`, `main.py`) |
| Cloud GPU port | Kento — adapting for MR cloud infrastructure (`asap/asap.c` → `kasap.c`) |
| Hooman | Learning, comparing, CS integration — workspace only |
| Faraz | Independent MATLAB implementation, same raw data |
| Pauli | Collaborator on recon validation |

**Goal chain:**
1. ✅ Setup complete
2. **← NEXT: Study major differences between Steve and Faraz implementations**
3. Reconstruct from raw Siemens `.dat` data
4. Build compressed sensing (CS) framework

---

## What To Do Next Session

### Primary Task: Deep-Dive Implementation Comparison

**Start here:** `workspace/reference/Recon_Comparison.md` — already has the 10-point high-level summary.

**Go deeper on these 4 key diffs (highest impact on image quality):**

1. **Kernel: Gaussian (Steve) vs Kaiser-Bessel (Faraz)**
   - Steve's active kernel: `recon.py` → `cudarecon()` — look for `exp(-dsq/kdist0sq)`
   - Steve's commented-out KB path: `recon.py` → `bessi0()` function (defined but not wired in)
   - Faraz's KB: `workspace/codes/2023_Faraz_Recon_HH/createKBkernel.m` (6 lines) + `gridrecon_fa_20230113.m`
   - **Question to answer:** What would it take to swap Gaussian → KB in Steve's Python path?

2. **Density compensation: knorm (Steve) vs iterative DCF (Faraz)**
   - Steve: `recon.py` → `cudarenorm()` — per-cell k/knorm division
   - Faraz: `workspace/codes/2023_Faraz_Recon_HH/iterative_dcf_fa_20190910.m` — 5 CG iterations
   - Reference: `~/Library/.../ASAP Recon/ASAP Faraz vs Steve.md` section "2) Density compensation"
   - **Question to answer:** How many DCF iterations are needed? What's the convergence criterion?

3. **Coil combine: per-coil complex (Steve) vs adaptive Roemer magnitude (Faraz)**
   - Steve: `results.py` → `getimg()` returns per-coil or linear combine
   - Faraz: `workspace/codes/2023_Faraz_Recon_HH/combinecoils_fa.m`
   - **This explains the visual background difference** (zero-mean vs Rician-positive)

4. **Trajectory handling differences**
   - Steve: `raw.py` → `traj.load()` — infers `npts` from |k|² periodicity
   - Faraz: `loadtrajectory3D.m` — loads from calibrated `.mat`, returns coords in 1/mm
   - Reference: `~/Library/.../ASAP Recon/ASAP Trajectory.md`

### Suggested Workflow

1. Read `workspace/reference/Recon_Comparison.md` first (context)
2. Read deep-dive Obsidian note: `~/Library/.../ASAP Recon/ASAP Faraz vs Steve.md` (full 10-point)
3. Side-by-side code reading:
   - `recon.py` (Steve's gridding) alongside `workspace/codes/2023_Faraz_Recon_HH/gridrecon_fa_20230113.m`
   - Focus on kernel, DCF, neighborhood size
4. Write findings to `workspace/reference/Recon_Comparison.md` (extend the existing doc)
5. After understanding: plan what to change in `workspace/codes/2023_Faraz_Recon_HH/` to port quality improvements

---

## Key Files for Next Session

| File | Why |
|------|-----|
| `recon.py` | Steve's gridding kernel — understand Gaussian + see commented KB |
| `raw.py` | Steve's trajectory loading — understand npts inference |
| `results.py` | Steve's image output + coil path |
| `workspace/codes/2023_Faraz_Recon_HH/gridrecon_fa_20230113.m` | Faraz's gridding |
| `workspace/codes/2023_Faraz_Recon_HH/iterative_dcf_fa_20190910.m` | Faraz's DCF |
| `workspace/codes/2023_Faraz_Recon_HH/createKBkernel.m` | Faraz's KB kernel |
| `workspace/reference/Recon_Comparison.md` | Existing summary — extend it |
| Obsidian: `ASAP Faraz vs Steve.md` | Deep analysis notes (already written) |
| Obsidian: `ASAP Recon Pipeline.md` | Ideal pipeline plan (future target) |

---

## Pitfalls / Gotchas to Remember

- **Git boundary:** NEVER run git commit/push commands in `2026_ASAP_Recon/` root (update with pull is allowed). All git ops in `workspace/` only.
- **MATLAB filename rule:** Any `.m` file rename breaks all callers — MATLAB resolves functions by filename.
- **Stale lookup tables:** `grid_lookup_20220418.mat` caches neighbor lookup — delete if trajectory/FOV changes.
- **Hardcoded paths in Faraz's script:** `spiral_human_20240227.m` has Windows paths (`C:\Users\faraz\...`) — must update before running.
- **Steve's background looks negative:** That's complex output (zero-mean). Take `np.abs()` to match MATLAB Rician output.
- `kasap.c` (Kento) and `asap/asap.c` (Steve) differ — don't assume equivalence.

---

## Suggested Skills for Next Session

- `/handoff update` — after completing comparison study, update this doc with findings
- `intent-layer:intent-layer-maintenance` — if reference docs get substantially updated
- No other special skills needed for code reading work

---

## Workspace Git State

**Remote:** `git@github.com:HoomHam/asap-recon-workspace.git`  
**Branch:** `main`  
**Last commit:** `a429655` — "add: codes folder (MATLAB + C implementations)"  
**Status:** Pushed to GitHub, clean working tree  
**Size:** `.git/` = 1.5MB (large .mat files removed from history via git-filter-repo)  
**.gitignore:** Excludes `.DS_Store`, `*.mat`
