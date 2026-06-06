# Workspace — ASAP Recon (Hooman's)

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
*(populate as created)*

| Doc | Contents |
|-----|---------|
| `reference/Recon_Overview.md` | High-level algorithm map — Steve's pipeline |
| `reference/Recon_Comparison.md` | Steve vs Faraz implementation diff |
| `reference/Physics_Notes.md` | k-space, gridding, Xe-129 gas exchange concepts |

## People Quick Ref
- **Steve** — original author, check his code for ground truth
- **Kento** — cloud GPU adaptation, may have changed GPU dispatch
- **Faraz** — parallel impl, same raw data, different codebase
- **Hooman** — that's you, learning + CS integration

## Session Startup
1. Read latest `handoffs/ASAP_Handoff_*.md` (most recent date)
2. Check `reference/` for relevant domain doc
3. Main code lives at `../` — read-only, no git ops there
