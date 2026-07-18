# ASAP Recon — Ledger

Append-only session history. Decisions, RETRACTIONS, BRANCH lines.
Never edit old entries. Format per canon law (~/.claude/CLAUDE.md).

---

## 2026-07-12 — adopted into canon

- Project opted into canon (full) by /adopt. Retro-fill harvested from root
  handoff-report.md + workspace handoffs (Jun 2026).
- Repo law (figured out & confirmed this session): root repo has TWO remotes —
  `origin` = MEDCAP/asap_recon (Kento, PULL-ONLY, never push) and
  `hooman` = HoomHam/asap_recon (own fork). Dev branch `diaphragm-recon`
  tracks the fork, 14 commits ahead of origin/main. workspace/ = separate
  repo (HoomHam/asap-recon-workspace), hidden from root git via
  .git/info/exclude.
- Stale root graphify-out/ (21M, Jun 15) deleted per System handoff step 2.

## 2026-07-17 · session 2

- Traced Steve's navigator binning end-to-end: 3 surrogates (SIGNAL/PNEUMOTACH/DIAPHRAGM),
  `bin()` at raw.py:153 (limb-split → fold → searchsorted RANK → equal-count phase [0,1)),
  ×nbins → cyclic Gaussian soft-weight in gridder (recon.py:104-111, nbins=16). asap.c has
  no binning (trajectory generator only).
- Wrote piston-framed reference doc `~/Downloads/Steve_Diaphragm_Binning_Piston.md` (zero
  biology in prose) + explanatory figures (amplitude-rank → representative cycle).
- Established the DIAPHRAGM 4D stack has NO native time axis (bins are amplitude-rank).
  Recipe to derive per-bin time = EE-detect → within-cycle τ → soft-weighted CIRCULAR mean
  (caught the phase-wrap bug). Parked plan: workspace/helpers/recon/bin_time_PLAN.md. Not built.
- XeCS relayed a directive to force the navigator window to 26 interleaves (=1 Thomson set,
  832=32×26). MEASURED on real v3 data (v3_dyn / v3_fov250 / v3_dyn_025JC): nuniqueilvs=832,
  nusimg=32 → ilvperusimg = 26 ALREADY. Never 20. No functional fix needed.
- Added diagnostic-only edit to root tyger_recon.py `_diaphragm_navigator`: comment +
  THOMSON_SET=26 warn-if-mismatch + nav-frame print; hoisted ilvperusimg, removed redundant
  recompute. Zero behavior change. Left LOCAL (uncommitted) — fork push is Hooman's.
- Wrote a cross-repo prompt for the other session to diff its nav impl vs this, update its
  docs, and fix only if a conventional-ASAP (non-CS) path carries a wrong inherited window.
- RECONCILE: large untracked backlog in workspace/ from June sessions (pipeline/, autorun
  archives, outputs/{016PG,023LL,025JC,25JC,piston}, several helpers + reference md) never
  carded/committed. Flagged for a dedicated reconcile pass — NOT retro-carded here (provenance
  uncertain; avoid fabricating cards).

- DECIDED: all generated figures go to workspace/outputs/<task>/ (rule; saved to project memory figures_to_outputs.md).
- DECIDED: keep the tyger_recon.py nav diagnostic as a future-ASAP alignment guard.
- DECIDED: root tyger_recon.py change stays local; Hooman commits/pushes the fork himself.
- BRANCH: diaphragm-binning, parent: B8 (DIAPHRAGM navigator). Investigation of the binning
  algorithm + nav window + 4D time vector. Parked (bin_time.py unbuilt).
