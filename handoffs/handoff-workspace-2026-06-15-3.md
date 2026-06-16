# ASAP Recon — Workspace Session Handoff

**Date:** 2026-06-15 (workspace/CS workstream — NOT the root Tyger handoff)
**Scope:** Hooman's CS recon work only. Root `./handoff-report.md` is a separate
Tyger/main-repo handoff — do not mix.
**Phase closed this session:** Lustig MATLAB CS reproduced as a one-command tool,
compared against our CS on the **same** v3 data — **ours wins**. Parameter/conditioning
reference for the Lustig code written.
**Next phase:** sweep the Lustig input space to find ITS optimum, compare our-optimum
vs Lustig-optimum head to head.

## Project nature (read first)
Benign engineering: comparison/reimplementation of MRI recon code on **physical
phantom data** (ACR). No human data, no PHI. Model: Hooman's project default is
**Fable 5** (CLAUDE.md). This session ran on **Opus 4.8** (switched mid-session;
Fable 5 access was lost from the picker — account/plan-tier issue, not repo). If a
fresh session is not on Fable 5, mention it.

## What got done this session

1. **One-command Lustig CS pipeline** — `workspace/helpers/lustig_oneshot/`
   (`run_lustig.py <recon_io>`). Collapses the old MATLAB→Colab→Drive→MATLAB dance:
   reads recon_io npy directly → builds `ACR_test.mat` with exact torchkbnufft pipe
   DCF → headless MATLAB `fnlCg` → `gas(15,100,100,100)`. Node: `lustig_oneshot/AGENTS.md`.
   Ran clean on v3_fov250 (~32 s DCF + ~4 min CS).
2. **Same-data comparison** — Lustig CS vs our CS, both on v3_fov250. **Our CS wins**
   (lowfreq-CV 0.084/0.085 vs Lustig 0.115; Lustig SNR 99 is the background-zeroing
   blind-spot, not real gain). Lustig is visibly softer. Montage:
   `data/v3_fov250/recon_io/montage_lustig_v3_vs_ours.png`. Tool:
   `helpers/recon/lustig_compare.py`.
3. **Corrected a stale belief:** the `(256,64,256)` DCF size bug from the 2026-05 docs
   is NOT in the real workflow — the notebook used sets `im_size=(100,100,100)`.
4. **Found the earlier compare was cross-data:** `tv01g01.mat` (and the other two
   saved runs) were the **old 2025-08-16 ACR** scan (425 984 samples), not v3
   (424 320). The one-shot fixes that.
5. **Decoded the 3 Lustig run scripts** and wrote a full parameter/conditioning
   reference (see new docs). Key findings: `_Wavelet_Adop` is misnamed (XFM=1, both
   weights 0 → pure least-squares + edge taper, no sparsity); `_Wavelet` has a real
   bug (`FWT2_PO` is 2D applied to a 3D volume, forces 128³); DCF used only in init
   so **final ≈ DCF-gridded init** (stalled NCG).
6. **Intent-layer patched** (commit `ba8e490`): new `lustig_oneshot/AGENTS.md`, updated
   `helpers/recon/AGENTS.md` status, `workspace/CLAUDE.md` doc table.

## New / key paths

| What | Path |
|------|------|
| One-command Lustig CS | `workspace/helpers/lustig_oneshot/` (`run_lustig.py`, `build_acrtest.py`, `run_cs.m`, `README.md`, `AGENTS.md`) |
| Torch venv for DCF step | `workspace/helpers/lustig_oneshot/.venv_lustig` (gitignored; torch+torchkbnufft arm64) |
| Our-CS-vs-Lustig montage tool | `workspace/helpers/recon/lustig_compare.py` |
| Lustig pipeline + verdict | `workspace/reference/Lustig_CS_Baseline.md` |
| **Lustig params + what to sweep** | `workspace/reference/Lustig_CS_Tuning.md` ← read first next session |
| Lustig run scripts (3 variants) | `workspace/codes/2025-09-24_ACR/spiral3d_cs_3D_hoom*.m` |
| Saved Lustig runs (OLD ACR data) | `workspace/codes/2025-09-24_ACR/{tv01g01,wlet_ph1,wlet_ph1_adopt}.mat` |
| Lustig source toolboxes | `workspace/codes/2025_CS/` (IRT + sparseMRI_v0.2) |
| v3 same-data outputs | `workspace/data/v3_fov250/recon_io/lustig/` (`ACR_test.mat`, `lustig_cs.mat`) |
| Our CS code | `workspace/helpers/recon/cs_recon.py` (+ `AGENTS.md`, `CS_Implementation.md`) |

All commits pushed to `HoomHam/asap-recon-workspace` through `ba8e490`.

## Next session — DO THIS FIRST

Goal: find the **Lustig optimum** by sweeping its input space, then compare
our-optimum vs Lustig-optimum (eye + metrics), same v3 data.

1. Read `reference/Lustig_CS_Tuning.md` — it has the priority-ranked sweep list.
   The headline order: **(a)** fix the wavelet to a real 3D transform (or slice-loop
   `FWT2_PO`) before trusting any wavelet result; **(b)** add DCF preconditioning to
   the *objective* (not just init) or the solver stays stalled and sweeps just buy
   smoothing; **(c)** then sweep `xfmWeight` ∈ {0.003,0.01,0.03,0.1}; **(d)** `TVWeight`;
   **(e)** `imSize` (128³ mandatory for wavelet); iterations LAST (metrics monotonic,
   no sweet spot — confirmed).
2. Write a NEW sweep script (don't edit the good helpers). Likely drives `run_cs.m`
   variants via the one-shot, or a parameterized `run_cs.m` that takes the knobs.
   Keep `run_lustig.py`/`run_cs.m` as the faithful-reproduction baseline; experiments
   go in new files.
3. Compare via `lustig_compare.py` (already undoes the MATLAB per-slice rot90; shared
   p99.5 axis). Metrics via `cg_tune.metrics` — but mind the SNR blind spot
   (background zeroing inflates it); trust lowfreq-CV + the eye.
4. **Eye-over-metrics rule (standing):** send Hooman the montage for a visual verdict
   before claiming any image-quality change.

## Carried-over open loops (do not lose)
1. **λ verdict for OUR CS — t0.003 vs t0.01, Hooman's eye, still pending.** Montages:
   `recon_io/montage_cs_t0.003.png` vs `montage_cs_t0.01.png`. Step-1 "fix metrics +
   freeze λ" never happened — still gates a "final" number (qualitative verdict won't
   change).
2. **BART `pics`** — the real independence test (MATLAB-Lustig done was personal-baseline;
   sigpy is Lustig-lab so sigpy-vs-sigpy circular). `brew install bart`; trajectory
   grid-index → BART units.
3. Tell **Faraz** about the ×1.205 `resizing` zoom bug (`zoom_check_metrics.json`;
   his `~/Hooman/Work/Analysis/2025-09-24_ACR/spiral_gpt_ACR_20250915.m:445-452`).
4. **Steve's rebuttal** to `workspace/for_steve.md` — pending.
5. **b-map stage** (Steve's `calcb` port) — before multi-coil / 4D.
6. Parked solver alts: ADMM-with-inner-CG (cleaner than DCF-weighted objective);
   undecimated wavelets.
7. Temporal 4D (diaphragm binning, bins as W_b weights, temporal TV/low-rank) — the
   eventual phase where undersampling returns.

## Hard constraints
- Our CS: run with `workspace/helpers/.venv/bin/python` (finufft 2.5.1, sigpy 0.1.27, arm64).
- Lustig DCF step: `helpers/lustig_oneshot/.venv_lustig` (torch+torchkbnufft) — separate venv.
- DCF for exact Lustig repro MUST stay torchkbnufft pipe at `(100,100,100)`; `voronoidens`
  is a different algorithm (sensitivity test only).
- NEVER git commit/push in repo root (`../`). All git in `workspace/` (`HoomHam/asap-recon-workspace`).
- Root `CLAUDE.md` Code Map is intentionally NOT updated with the new Lustig tools
  (main repo = no commits); navigation lives in `workspace/CLAUDE.md`.
- Never modify Faraz's code or the Lustig source toolboxes; scanner data stays in
  gitignored `workspace/data/`.
- `selftest.py` must pass after any operator change (operator untouched this session).
- Read `helpers/recon/AGENTS.md` + `lustig_oneshot/AGENTS.md` pitfalls before touching recon code.

## Two-handoff note (process)
This repo now keeps two handoffs: **root** `./handoff-report.md` = Tyger/main-repo;
**workspace** `./workspace/handoff-report.md` = this CS workstream. Archives share
`workspace/handoffs/`; workspace ones carry the `-workspace-` token. The `/handoff`
skill was updated this session to ask root-or-workspace before acting.

## Suggested skills
- `/handoff update <…>` (workspace) as sweep milestones land.
- `/code-review` on the new sweep script + `cs_recon.py` before the BART session.
- `/intent-layer-maintenance` only if drift accumulates — this session's drift was patched manually.
