# Workspace — ASAP Recon (Hooman's)

**Canon** (adopted 2026-07-12): ledger/facts/cards/tree in `canon/` — end sessions
with `/leave`. Cross-project retractions: `~/Hooman/Work/BELIEFS.md`.
(Root `CLAUDE.md` is Kento-repo-tracked — canon note lives here, not there.)

> Project nature: benign MRI recon **code comparison + dynamic recon** (phantom data:
> ACR phantom, piston-cylinder with HP gas). No PHI. Current focus: dynamic imaging
> of piston cylinder filled/evacuated with hyperpolarized gas. Hooman's model
> preference: **Fable 5**. See root `CLAUDE.md` → Session Notes.

## Git: ALLOWED HERE
Commit, push, pull freely. This `.git` is independent of the main repo.
Never go up to `../` and run git commands there.

## Purpose
Personal working space alongside Steve/Kento's ASAP recon code.
- Understand the algorithm
- Compare with Faraz's implementation
- Build toward CS (compressed sensing) integration
- Raw data reconstruction experiments

## Restructure note (2026-06-24)
The CS work split out into its own repo: **`2026_XeCS_Recon`** (operators in
`recon/`, 4D pipeline in `pipeline/`, comparison in `workspace/compare/`, Lustig in
`workspace/lustig/`). This workspace is now the non-CS / Steve-anchored side. The
ASAP scripts that still need CS operators import them from XeCS via an
`xecs_recon.pth` in `helpers/.venv`. The old CS originals are preserved (not deleted)
in `helpers/_delete/` — safe to remove once XeCS is in daily use. CS theory docs are
mirrored into `2026_XeCS_Recon/workspace/reference/` (XeCS is canonical for CS now).

## Token Efficiency Rules for Agents
- Read files on demand — never preload entire codebase
- Reference `../CLAUDE.md` for project identity and main code map
- For reconstruction physics: ask, don't assume
- Handoff doc is truth for session state — always read it first

## Folder Rules

| Folder | What goes here | File naming |
|--------|---------------|-------------|
| `handoffs/` | End-of-session state | `handoff-YYYY-MM-DD[-N].md` (workspace) / root kept as `handoff-2026-*.md` |
| `reference/` | Notes on live code, concepts, comparisons | `Domain_Description.md` |
| `archive/` | Old notes, superseded analysis (+ `archive/codes/` old code) | `Domain_Description_YYYY-MM.md` |
| `helpers/` | Hooman's non-CS scripts — Steve loaders, baseline, Faraz figures | subfolder by type |
| `faraz/` | Faraz's MATLAB fork reference (was `codes/2023_Faraz_Recon_HH/`) | — |
| `codes/` | `kasap.c` (Kento's C kernel) only — rest archived/moved | — |

`helpers/recon/` has its own node: `helpers/recon/AGENTS.md` (the ASAP-side stayers —
Steve data prep, FINUFFT-vs-Steve baseline, Faraz figures, cross-repo CS montage).
The CS operator library + 4D pipeline + Lustig one-shot moved to `2026_XeCS_Recon`
(see `2026_XeCS_Recon/recon/AGENTS.md`, `workspace/compare/AGENTS.md`).

## Active Reference Docs

| Doc | Contents | Status |
|-----|---------|--------|
| `reference/Recon_Comparison_StaticGas.md` | Steve vs Faraz full diff — static phantom, single bin, gas only. Scoring, theory notes, compute cost. **Single source of truth for this comparison.** | ✅ 2026-06-10 |
| `reference/Recon_Overview_Steve.md` | Code map — Steve's pipeline: file roles, data flow, entry points, GPU dispatch, magic numbers | ✅ 2026-06-10 |
| `reference/Recon_Overview_Faraz.md` | Code map — Faraz's MATLAB pipeline: file roles, data flow, entry points, compute model, magic numbers | ✅ 2026-06-10 |
| `reference/Physics_Notes.md` | Educational: non-Cartesian recon from FID + 2D DFT up — signal eq, DCF, gridding, kernels, coil combine, trajectory calibration, CS outlook; every concept mapped to Steve/Faraz code | ✅ 2026-06-10 |
| `reference/CS_Implementation.md` | Educational: CS layer on the finufft operator — objective, λ-as-threshold parameterization (two measured failures), DCF-as-preconditioner (CG vs gradient solvers), wavelet/TV priors, metric blind spots, first sweep results | ✅ 2026-06-12 |
| `reference/Lustig_CS_Baseline.md` | Hooman's pre-project MATLAB Lustig CS pipeline — density framework mismatch, why Steve's densities fail, diff vs our CS, **one-shot tool + same-data v3 verdict (ours wins)** | ✅ 2026-06-15 |
| `reference/Lustig_CS_Tuning.md` | Lustig CS parameter/conditioning reference — the 3 run scripts decoded, every NUFFT3D/TV/wavelet/fnlCg knob, init, conditioning caveats (DCF-in-init-only, 2D-wavelet-on-3D bug), priority-ranked sweep list | ✅ 2026-06-15 |
| `reference/Final_Report_CS_Comparison.md` | **Presentation-ready final report** — ours vs BART vs Lustig vs Steve vs Faraz. Narrative arc, code catalog (I/O per script), figure catalog (FINAL vs superseded), findings, reproduce-from-scratch, open threads. **Read this first for the comparison phase.** | ✅ 2026-06-15 |
| `reference/BART_Comparison.md` | Detailed chronology of the BART comparison + the texture/slice-matching saga with all dead-ends ruled out | ✅ 2026-06-15 |

## Archive

Original reference docs (created 2026-06-07) were accidentally deleted before the 2026-06-10 session; revived from ScreenPipe screenshots and archived:

| Doc | What it was | Notes |
|-----|------------|-------|
| `archive/Obsidian_index_2026-06.md` | Index of Obsidian ASAP Recon vault notes | Revived from screenshots; a few rows marked `[unverified]` |
| `archive/Recon_Comparison_2026-06.md` | Old Steve-vs-Faraz summary (LLM-generated) | Revived; superseded by `reference/Recon_Comparison_StaticGas.md`, which corrects its errors (magnitude-vs-real combine, coil-combine claims, de-apodization) |
| `archive/CS-Lustig-ASAP-{Pipeline,Status,Adaptation}_2026-05.md`, `archive/ASAP-Lustig-Reconstruction_2026-05.md` | Obsidian wiki pages on the pre-project Lustig CS pipeline (2026-05-22 Cowork audit) | Superseded by `reference/Lustig_CS_Baseline.md`; keep for full stage-by-stage detail |

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
