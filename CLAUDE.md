# Workspace — ASAP Recon (Hooman's)

> Project nature: benign MRI recon **code comparison** (phantom data, signal
> processing). Hooman's model preference for this project: **Fable 5**.
> See root `CLAUDE.md` → Session Notes.

## Git: ALLOWED HERE
Commit, push, pull freely. This `.git` is independent of the main repo.
Never go up to `../` and run git commands there.

## Purpose
Personal working space alongside Steve/Kento's ASAP recon code.
- Understand the algorithm
- Compare with Faraz's implementation
- Build toward CS (compressed sensing) integration
- Raw data reconstruction experiments

## Token Efficiency Rules for Agents
- Read files on demand — never preload entire codebase
- Reference `../CLAUDE.md` for project identity and main code map
- For reconstruction physics: ask, don't assume
- Handoff doc is truth for session state — always read it first

## Folder Rules

| Folder | What goes here | File naming |
|--------|---------------|-------------|
| `handoffs/` | End-of-session state | `ASAP_Handoff_YYYY-MM-DD.md` |
| `reference/` | Notes on live code, concepts, comparisons | `Domain_Description.md` |
| `archive/` | Old notes, superseded analysis | `Domain_Description_YYYY-MM.md` |
| `helpers/` | Hooman's scripts — loaders, plotters, CS tools | subfolder by type |

## Active Reference Docs

| Doc | Contents | Status |
|-----|---------|--------|
| `reference/Recon_Comparison_StaticGas.md` | Steve vs Faraz full diff — static phantom, single bin, gas only. Scoring, theory notes, compute cost. **Single source of truth for this comparison.** | ✅ 2026-06-10 |
| `reference/Recon_Overview_Steve.md` | Code map — Steve's pipeline: file roles, data flow, entry points, GPU dispatch, magic numbers | ✅ 2026-06-10 |
| `reference/Recon_Overview_Faraz.md` | Code map — Faraz's MATLAB pipeline: file roles, data flow, entry points, compute model, magic numbers | ✅ 2026-06-10 |
| `reference/Physics_Notes.md` | Educational: non-Cartesian recon from FID + 2D DFT up — signal eq, DCF, gridding, kernels, coil combine, trajectory calibration, CS outlook; every concept mapped to Steve/Faraz code | ✅ 2026-06-10 |

## Archive

Original reference docs (created 2026-06-07) were accidentally deleted before the 2026-06-10 session; revived from ScreenPipe screenshots and archived:

| Doc | What it was | Notes |
|-----|------------|-------|
| `archive/Obsidian_index_2026-06.md` | Index of Obsidian ASAP Recon vault notes | Revived from screenshots; a few rows marked `[unverified]` |
| `archive/Recon_Comparison_2026-06.md` | Old Steve-vs-Faraz summary (LLM-generated) | Revived; superseded by `reference/Recon_Comparison_StaticGas.md`, which corrects its errors (magnitude-vs-real combine, coil-combine claims, de-apodization) |

Obsidian deep-dive comparison (`Action/MRI/ASAP Recon/ASAP Faraz vs Steve.md`) also superseded by `Recon_Comparison_StaticGas.md` for static-gas scope.

## People Quick Ref
- **Steve** — original author, check his code for ground truth
- **Kento** — cloud GPU adaptation, may have changed GPU dispatch
- **Faraz** — parallel impl, same raw data, different codebase
- **Hooman** — that's you, learning + CS integration

## Session Startup
1. Read latest `handoffs/ASAP_Handoff_*.md` (most recent date)
2. Check `reference/` for relevant domain doc
3. Main code lives at `../` — read-only, no git ops there
