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

## 2026-08-27 · session 3

- Confirmed the batch recon cohort on Ext: 84 session folders in
  /Volumes/HoomHamExt/AIkill_Dynamic (renamed from Dynamic), all with complete s+d Tyger
  outputs; 78/84 also have p (5 subjects missing proton entirely); 3 sessions' d skipped
  dissolved (RBC/TP params absent: 2024-01-18_001JM, 01-22_007RA, 01-31_008CR).
- Built helpers/atlas/dis_atlas.py (C26): per-subject QC atlas — 60 IDs, dates sorted;
  per ID a gp panel page + dp panel page (orientation-major blocks, per date one 10-slice
  row) then per date 6 full all-slice pages (cor/sag/ax × gp/dp); one mp4 per ID×kind
  cycling the 16 bins ×5. Iterated 4 layout rounds with Hooman on 043AS sample.
- Cross-date column alignment solved deterministically (no registration): per-session
  slice picking at shared apex→base fractions of each session's OWN gas extent + fixed-size
  in-plane crop centered on each session's lung bbox. Verified visually (heart notch /
  mediastinum land in same column across 4 repeats of 043AS).
- Slice extent from GAS mask only (≥30 voxels above 0.15·max of bin-mean) — dp too noisy
  to self-select; extent genuinely reaches Z=99 (basal signal), not artifact.
- Output: workspace/outputs/aikill_atlas/dissolved_atlas.pdf (614 pages, 188 MB, NOT in
  git) + 119 videos (1.0 GB) moved to /Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/
  aikill_atlas/videos/ per big-output law, symlinked back at outputs/aikill_atlas/videos.

- DECIDED: atlas slice/crop selection always derives from gas_phase, never dp.
- DECIDED: sessions whose d lacks dissolved fall back to s-folder gas; dp rows/pages skipped.
- DECIDED: atlas videos live on Ext (regenerable — acceptable single-copy risk).
- BRANCH: aikill-atlas, parent: B10 (batch cohort). Cohort QC visualization. Delivered;
  open for per-ID audit findings.
- session 2026-08-31: C27 rt_prepost_fig.py — RT pre/post coronal figures (workspace/outputs/rt_prepost/). No new branch.

## 2026-09-03 · session 4

- Tyger access check for 005JJ. Service-principal cert `tyger-sp.pem` had expired 2026-08-01
  (AADSTS700027); device-code login attempted, never claimed. Kento emailed a new pem
  (CN=TygerTEP-CLIlogin, valid to 2028-08-04) → installed, old kept as
  `tyger-sp.pem.expired-2026-08-01`, `tyger login -f LOGIN_FILE.yml` → role owner.
- 005JJ (2024-05-13, single-coil 67 MB dat, source `/Volumes/HoomHamExt/_5t_images_roundtrip/
  Images/2024-05-13/005JJ`) converted with SIGNAL binning and run on Tyger gpunp as run 558
  using the OLD image `db80f16` (numba-working; guard not needed for nch=1). Succeeded, GPU
  wall 92 s, log clean (16 bins gp+dp, RBC/TP phase solved every bin).
- Download is the weak link: ~1 Mbps effective. `run exec --logs` streaming blew the 10-min
  Bash cap at 52 MiB; `tyger buffer read <id> > file` (dop 32) died at 44 MiB with
  "context deadline exceeded"; `tyger buffer read <id> -o output.mrd -p 4` in background got
  the full 192 MB in ~10 min. plot_recon montage OK (lungs in all 16 bins, expiration dip
  bins 5–8). Output byte-size equals the July run; content differs from byte 30289 (header,
  not inspected).
- Results on Ext (mirror law): /Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/
  tyger_check_2026-09-03/2024-05-13_005JJ_s/ (output.mrd, input.mrd, tyger.log, pngs,
  codespec_db80f16.yml). No p/d, no post_process, by request ("nothing else").
- Left untouched: pipeline/recon_codespec.yml still points at broken d136eb1 (F1 unchanged);
  the July pipeline/post_process.py + codespec + June/July untracked backlog remain
  uncommitted (not this session's). 2026-08-31 session (C27) had appended ledger/cards but
  never committed — included in this commit.

- DECIDED: single-coil subjects can run TODAY on image db80f16; F1 only blocks nch=8 (guard) subjects.
- DECIDED: dyn_recon.submit() must stop relying on `run exec --logs` streaming for the output; switch to run + `buffer read -o -p 4` with retry before the batch resumes.
- BRANCH: tyger-access-check, parent: B6 (tyger). Closed same session (delivered).
