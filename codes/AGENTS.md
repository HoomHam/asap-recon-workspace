# codes/ — Hooman's Working Code Copies

## Purpose

Active code Hooman is modifying. NOT read-only mirrors.
Commit freely — this is inside the workspace git.

## Contents

| Item | Origin | Status | Intent |
|------|--------|--------|--------|
| `2023_Faraz_Recon_HH/` | Copied from `/Users/hoomham/Hooman/Work/Codes/2023_Faraz_Recon/` | **Active — Hooman modifying** | Hooman's MATLAB recon: understand + improve Faraz's pipeline |
| `kasap.c` | Kento's cloud GPU C kernel | Reference copy | Study Kento's C adaptation of Steve's ASAP algorithm |

## `2023_Faraz_Recon_HH/` — Hooman's MATLAB Fork

**Original author:** Faraz  
**Hooman's goal:** Learn the pipeline, then improve it (KB kernel already used, iterative DCF, coil combine)

Key files to work in:

| File | Role | Modification priority |
|------|------|-----------------------|
| `spiral_human_20240227.m` | Main study script — entry point | Rename/adapt for Hooman's data paths |
| `gridrecon_fa_20230113.m` | Core gridding (KB kernel + DCF) | Understand first, modify later |
| `fa_spiral_dyn_recon.m` | Dynamic recon orchestrator | Understand first |
| `grid_lookup_20230113.m` | Nearest-neighbor lookup table | Likely keep as-is |
| `iterative_dcf_fa_20190910.m` | Iterative DCF — quality path | Likely keep as-is |
| `combinecoils_fa.m` | Adaptive Roemer coil combine | Likely keep as-is |
| `createKBkernel.m` | KB kernel builder | Likely keep as-is |

**Do NOT rename files without updating all callers** — MATLAB uses filename-as-function-name.

## `kasap.c` — Kento's C Kernel

C implementation of ASAP gridding for cloud GPU deployment.
Compare with `../../asap/asap.c` (Steve's version) to understand what Kento changed.
Do not modify unless explicitly working on C kernel.

## Code Map

| Looking for... | Go to |
|----------------|-------|
| Entry point for MATLAB recon | `2023_Faraz_Recon_HH/spiral_human_20240227.m` |
| Gridding + DCF | `2023_Faraz_Recon_HH/gridrecon_fa_20230113.m` |
| Kento's C kernel | `kasap.c` (compare with `../../asap/asap.c`) |
| Full file map for MATLAB fork | `2023_Faraz_Recon_HH/AGENTS.md` |

## Entry Points

| Task | Start Here |
|------|------------|
| Run MATLAB recon on new data | `2023_Faraz_Recon_HH/spiral_human_20240227.m` |
| Compare Kento vs Steve C kernels | `kasap.c` vs `../../asap/asap.c` |

## Contracts

- `2023_Faraz_Recon_HH/` is a fork — modifications tracked by workspace git
- `kasap.c` is reference only — do not modify unless actively working on C kernel
- Original sources are read-only: `../../` (Steve) and `/Users/hoomham/Hooman/Work/Codes/2023_Faraz_Recon/` (Faraz)

## Pitfalls

- MATLAB filename = function name — renaming any `.m` file breaks all callers silently
- `kasap.c` and `../../asap/asap.c` differ — do not assume equivalence without diffing
