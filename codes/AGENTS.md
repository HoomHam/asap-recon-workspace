# codes/ — Reference C Kernel

## Purpose

Holds Kento's cloud-GPU C kernel reference copy. Commit freely (inside workspace git).

> **Restructure note (2026-06-24):** Faraz's MATLAB fork moved to `../faraz/`.
> Old code folders (`2025_CS/`, `2025_Xe129_CS/`, `2025-09-24_ACR/`) moved to
> `../archive/codes/`. CS operator library + 4D pipeline moved out to the new
> `2026_XeCS_Recon` repo. Only `kasap.c` remains here (root CLAUDE.md points at
> `workspace/codes/kasap.c`).

## Contents

| Item | Origin | Status | Intent |
|------|--------|--------|--------|
| `kasap.c` | Kento's cloud GPU C kernel | Reference copy | Study Kento's C adaptation of Steve's ASAP algorithm |

## `kasap.c` — Kento's C Kernel

C implementation of ASAP gridding for cloud GPU deployment.
Compare with `../../asap/asap.c` (Steve's version) to understand what Kento changed.
Do not modify unless explicitly working on the C kernel.

## Code Map

| Looking for... | Go to |
|----------------|-------|
| Kento's C kernel | `kasap.c` (compare with `../../asap/asap.c`) |
| Faraz MATLAB fork | `../faraz/` (entry `spiral_human_20240227.m`, map `../faraz/AGENTS.md`) |
| Archived old code | `../archive/codes/` |
| CS operators + 4D pipeline | `2026_XeCS_Recon/` repo (see `../CLAUDE.md`) |

## Contracts

- `kasap.c` is reference only — do not modify unless actively working on C kernel
- Original sources read-only: `../../` (Steve), `/Users/hoomham/Hooman/Work/Codes/2023_Faraz_Recon/` (Faraz)

## Pitfalls

- `kasap.c` and `../../asap/asap.c` differ — do not assume equivalence without diffing
