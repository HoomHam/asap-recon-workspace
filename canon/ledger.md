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

## 2026-09-06 · misplaced session (no ASAP work)
- A session opened in this folder for a Mac clock question and then did ~/Hooman/System work (Jon Snow, vault sweep, Brain/Canon mirrors). Nothing ASAP changed. Record lives in System's ledger session 29. The transcript the exit hook drops into workspace/archive/ (cff52ba0) is a duplicate of System's copy — safe to delete.

## 2026-09-12/13 · session 5
- Inventory of v2/v3 dynamics on Ext roundtrip: 27 not reconstructed → fixed the three blockers
  in the fork (Hooman-approved push to `hooman/diaphragm-recon`): `472fbc9` numba/scipy/numpy
  pins (F42), curve_fit-failure fallback, optional dp trajectory; `40d23a4` keeps dissolved as
  unsplit complex + `rbc_tp_separated` tag. Hooman's nav diagnostic + gtypes.py left uncommitted.
- Reworked `pipeline/dyn_recon.py` submit (buffer write / run create / poll / buffer read -p 4,
  retries, runs dir on Ext, `--codespec`, MPLBACKEND=Agg) + `todo_batch.sh` (FORCE=1). Hit and fixed
  `--ttl` HTTP 500 (F44) and a 40-min plt.show() hang (F45). 21 sessions reconstructed on Tyger;
  unusable (no coil signal): 019WR + three tiny 000LL; duplicates skipped (SD03, 000LL 01-17) and
  bigmac 000LL output quarantined on Ext.
- Full audit of `_5t_images_roundtrip/Images` (136 folders): 104 unique sessions done (96 real /
  64 subjects + 8 000-series), 4 failed, 4 duplicate, 20 no trajectory (2022 seqs, F51), 2 not dynamic.
- Provenance: Steve original (branch main 3303276) now kept as plain non-git snapshot
  `Codes/2026_Steve_Recon` (renamed by Hooman, worktree detached); diffs in outputs/steve_vs_fork_diff.
  DIAPHRAGM z-flip verified as a real deviation from Steve (F49).
- SNR_Table_All.xlsx: 11 rows greened (94 total); own measured SNR cols I–S (gas/DP max, min, at EI,
  bins, method, note, split validity) from `helpers/snr_calc.py`; symlinked into XeCS + PCA data/.
- Found Steve's RBC/TP split unreliable (F47): phase stop at ΣaTP=0 jump, 60/86 sessions off target,
  noise gain 1/|sinΔφ| with corr −0.99. DP SNR switched to whitened magnitude (F48). Wrote
  `2steve/04_RBCTP_Split_PhaseStop_Conditioning.md` (not sent) + learning note
  `reference/DP_RBCTP_Split_Noise_Diagnostics.md`. Flagged atlas dp pages (C26) suspect.

- DECIDED: dissolved SNR = covariance-whitened magnitude for all sessions; never report per-part RBC/TP SNR.
- DECIDED: new SNR table columns = max / min / at EI (EI = max lung-volume bin, DIAPHRAGM bins); end-exp dropped from table (kept in CSV).
- DECIDED: 2026_Steve_Recon stays a plain folder (no git) — reference for diffing only.
- DECIDED: bigmac 2023-11-02 000LL output quarantined (byte-duplicate of 2023-11-03), not deleted.
- DECIDED: SNR table shared to XeCS/PCA by symlink (source of truth stays ASAP workspace/data).
- RETRACTION: F1 + F41 → F42 (numba pins fixed the pipeline; db80f16 no longer needed).
- RETRACTION: F40 → F43 (download speed varies; submit() no longer streams).
- RETRACTION: ledger 2026-08-27 "d skipped dissolved (RBC/TP params absent: 001JM, 007RA, 008CR)" → F46 (gas-only acquisitions).
- BRANCH: leftover-recon (B14), parent: B10. Closed (delivered).
- BRANCH: steve-provenance (B15), parent: B8. Closed (z-flip verified, snapshot + diffs).
- BRANCH: snr-table (B16), parent: B14. Closed (delivered).
- BRANCH: rbctp-split-audit (B17), parent: B16. OPEN ★ (note to Steve unsent; atlas dp pages suspect).
- BRANCH: tyger-access-check (B13) — closed ✅ (superseded by B14).

## 2026-09-14 · session 6
- Read-only survey for Hooman (requested via System session system-f3): documented the Siemens .dat → MRD → Tyger
  path ahead of the 2026-09-15 Spinhance × Oxford gas-imaging meeting (raw-data exchange, their modified ISMRMRD).
  Note: `notes/Data_SiemensToMRD_Pipeline_2026-09-14.md` (DRAFT, file:line cited). No code changed, nothing run.
- Findings: our MRD is the MEDCAP MRD v2 fork with a private NdArray+meta layout, unreadable by stock v2 or v1
  tools; gas/dissolved split and Xe constants live in raw.py, not the file (F53). Scanner measured as Avanto
  1.494 T, syngo MR D13, VD multi-RAID .dat (F54, n=1). No "MRD1" term exists anywhere in repo/canon/vault.
- Also noted: dyn_recon default `--codespec` still → broken d136eb1 (C13); `reference/Tyger_Setup.md` stale
  (shows .dat as --input, old basefolder).
- CONFLICT (A5, Hooman to decide): root + workspace CLAUDE.md say "phantom data, no PHI", but the cohort pipeline
  runs on human sessions (.dat headers + session folders carry PHI; IRB cloud letter covers Tyger, not an outside group).
- DECIDED: the pipeline note stays in this repo (`workspace/notes/`), per Hooman.
- BRANCH: siemens-mrd-doc (B18), parent: B6 (tyger). Closed (delivered). ★ stays B17.
- Model: Opus 5 (not Fable 5).
