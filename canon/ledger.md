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

## 2026-09-16 · session 7
- Hooman asked for a way to rank the 104 gas recons by intra-lung "patchiness" (COPD-like bright/dark
  patches with edges, NOT noise, NOT coil shading, lung size must not win). Built `helpers/structure_rank.py`
  (C30): per coronal slice at the max-SNR bin, R = I / local mean (σ10) removes shading/gravity ramp, σ3
  smoothing kills noise + vessels, G = rim-weighted mean |∇R| minus background noise floor, area-weighted
  over ALL slices (per-pixel mean ⇒ size cancels). Second column P = perimeter(bright 5σ)/perimeter(lung
  extent) = fragmentation. Extent = 5σ mask closed r=8 per midline side (bridges defects ≤16 px; never
  bridges the mediastinum). Dead ends recorded in the docstring: 2σ extent adds a partial-volume band
  whose ramp reads as a rim edge; zeroing outside the mask + σ3 smoothing fakes a rim; r=5 closing
  bridged the mediastinum. Outputs `outputs/structure_rank/` (panel_trachea_slice.png = all sessions,
  trachea coronal slice, sorted; panel_maxarea_slice.png; scatter_snr_vs_G.png; structure_table.csv).
  Hooman: "the panel would do and I pick myself" — panels delivered, no further metric tuning.
- Trachea slice = among slices with ≥40 % of max lung area, max signal in the midline column band
  2–25 rows below the apex; 16/103 fell back to the max-area slice (*). 2024-05-13_002ZS absent
  everywhere: d/output.mrd exists but no recon.mat (post_process never ran).
- Cross-session (XeCS session, k0 classification of all 201 ASAP twix): audited which twix each AIkill
  recon used. dyn_recon rule (largest spiral-dyn .dat) re-applied to the roundtrip Images folders =
  XeCS k0 "primary" for 104/104; the 25 leftover logs name the same file 25/25 (F55).
  `outputs/twix_audit_2026-09-16/twix_audit.csv`.
- DECIDED (Hooman, relayed by the XeCS session): four sessions with a GOOD second free-breathing dynamic get
  re-reconstructed with BOTH dynamics merged into one input: 2023-11-02_025VP (MID02741+MID02743),
  2024-01-31_008CR (MID00194+MID00199), 2024-03-06_030DN (MID00545+MID00543), 2024-08-12_013VM
  (MID00408+MID00410). Built `pipeline/merge_dyn.py` (C31) + `dyn_recon.py --merge/--ref` (C9) +
  `pipeline/merged_batch.sh` (C32). Why the merge is recon-legal (F56): time axis = index·TR (MDH
  timestamps unused), cal block = first numspec lines (file B's stripped), spiral arm = gas ilv index mod
  nuniqueilvs (832 v3 / 640 v2) so file A is trimmed to whole arm cycles (<14 s of dose-out tail).
  Validated on 030DN through raw.load_from_arr: gas-dissolved pattern, 832 cycle, RBC/TP fit converged,
  and the loader's own low-SNR filter masks the 15 s dead head of file B + dose-out tail. p binning NOT
  supported on merged input (pneumotach is one wall-clock record) → s+d only. Outputs go to
  `AIkill_Dynamic/<date>_<id>_merged/{s,d}` (originals untouched). Batch launched 20:55 (nohup), image
  40d23a4, logs `pipeline_runs/logs/<session>_merged.log`.
  NOTE 008CR: MID00194 is acquired BEFORE MID00199 (the recon's primary), so the merged head file =
  MID00194 (its cal block is the one used); 8-ch gas-only.
- STANDING RULE requested by Hooman (via XeCS session): before reconning any session, check XeCS
  `workspace/outputs/calspec/k0_sessions.csv` column `extra_free` and merge the good extra dynamics with
  `--merge`. Recorded here + cards; the workspace/CLAUDE.md line is NOT written yet — a peer session
  cannot authorise a CLAUDE.md edit, Hooman to confirm in this session (A6-style: Hooman flips it).
- DECIDED (Hooman, directly, 23:05): rule added to workspace/CLAUDE.md, SCOPED — the 104 cohort is already
  screened, the four merges above are the only ones; the rule applies to future / unscreened data only.
- RECEIVED (written by the XeCS session, not this one — lane note): `notes/calspec_package_2026-09-16/` (18 MB
  COPY; source of truth = XeCS `workspace/outputs/calspec/`). README.md → session_table.csv (104 rows: recon
  twix, breath-hold grade/T1,eff, extra_free_dynamics, spec_block, pooled spectrum numbers, phase-8/16 validity),
  k0_classify.csv, k0_sessions.csv, k0_contact_sheet.png, cohort_map.csv (HC/HC_OLD/EBV/LTX/RT),
  Recon_Input_Rule_2026-09-16.md (the standing rule text). Also sent to the PCA-registration session.
- BRANCH: structure-rank (B19), parent: B16 (snr-table). Closed (delivered).
- RESULT (23:15, all 4 merged sessions done, s+d, image 40d23a4; laptop crash 21:12–23:00 did not touch the
  nohup batch; 025VP/008CR failed once — B without cal block / no dp trajectory — fixed in merge_dyn.py and
  relaunched). Gas SNR at EI (snr_calc, d): 030DN 24.9→27.7, 013VM 27.4→30.1, 008CR 29.4→31.9,
  **025VP 19.1→14.7 (WORSE)**. DP: 030DN 15.8→18.2, 013VM 11.9→14.8, 025VP 10.9→8.4. Why 025VP loses:
  MID02741 is a weak short dose (k0 SNR 53 vs 126 for MID02743); Steve's bins average interleaves with equal
  weight, so a dose 2.4× weaker dilutes the mean signal while noise adds → SNR drops. Merge helps only when
  the two doses are comparable. Recommendation: keep the ORIGINAL 025VP recon; `_merged` kept on Ext for
  the record (bonus: it has an RBC/TP fit, the original had no cal block). `outputs/merged_dyn_2026-09-16/`
  (snr_merged_vs_orig.csv, panel_orig_vs_merged_gas_d.png, 030DN k0 seam plot).
- DECIDED (Hooman, 23:25): 025VP keeps its ORIGINAL recon (MID02743 only). Ext folder renamed
  `2023-11-02_025VP_merged_REJECTED/` + README marker. Rule in workspace/CLAUDE.md now names the 3 merges
  and the comparable-dose condition.
- BRANCH: merged-dyn (B20), parent: B14 (leftover-recon). Closed (delivered). ★ → B20.
- 2026-09-17 01:50 (B19 reopened): Hooman wants the two panels as VIDEOS over all 16 bins and the
  CARINA (bifurcation) as the governing slice. `carina_slice()` = Y template (trachea column × two
  diagonal arms anchored at apex/midline) on the max-over-bins volume; slices with ≥30 % max lung area;
  '*' = column < 5σ (9/107 weak). Dead ends: absolute-brightness airway mask grabs bright lung bases;
  lateral-width-of-airway-component is fooled by disconnected hilar patches. tiles.npz now holds
  (bins, z, x) per session; `panel_video()` → mp4 + gif via ffmpeg, frames kept in `*_frames/`.
  Outputs: panel_bifurcation_slice.{png,mp4,gif}, panel_maxarea_slice.{png,mp4,gif}. Old
  panel_trachea_slice.png deleted (superseded). Panels now include the 3 adopted `_merged` sessions
  (107 rows); `REJECTED` folders excluded.
- 02:10 per-COHORT panels (Hooman: healthies separate, pre/post adjacent, longitudinal adjacent). Cohort
  from XeCS `cohort_map.csv` (sid → HC 8 sessions / HC_OLD 4 / EBV 27 / LTX 53 / RT 14, all 107 mapped).
  Inside a cohort: one block per subject, sessions date-ascending (merged right after its original),
  blocks ordered by the subject's max G_struct. `outputs/structure_rank/cohorts/panel_<cohort>_
  {bifurcation,maxarea}.{png,mp4,gif}` (+ `_frames/`). Flags: `--cohort-map`, `--cohorts-only`.
- /leave 2026-09-17 03:30. Reconciled: `outputs/archive/` had no index row (added, untagged). Left dirty,
  not mine: handoffs/handoff-workspace-2026-06-*.md, pipeline/recon_codespec.yml, helpers/build_status_tab.py
  (untracked but carded C14), archive/*_auto_*.jsonl exports, canon/.nudge. Root repo: nothing committed
  (root is Kento's read-only repo; root handoff-report.md regenerated, untracked there by design).
- Model: Fable 5.1.
- 2026-09-17 ~04:00 session (Fable 5.1): Hooman order — CT overlay for the RT subjects. Ext
  `Work/CT/Data/{01BB,02ZS,03PM}` have Siemens RT 4D-CT (8 phases + Average, 3 mm). Converted to NIfTI
  (C33, Ext mirror path `Codes/2026_ASAP_Recon/ct_overlay/<subj>/ct/`). Asked PCA session for registered
  stacks: 03PM has legacy final (gas only), 01BB/02ZS none → elastix one-shot best; all 16×100³, identity
  affine, orientation arbitrary; FOV 350 (raw.py) ⇒ 3.5 mm. 01BB CT 2023-03-27 = Tyger 2023-03-27_02BB
  (same day; NOT 2023-10-10_002BB). Rigid only (Hooman: no deforming). C34 orientation search launched.
- BRANCH: ct-overlay (B21), parent: B14 (leftover-recon). Open. ★ → B21.
- RESULT (04:40): C34 done for 3 subjects. Orientation F58 (S=−Z, P=−Y, L=+X) — unanimous across subjects.
  Rigid (full-res Mattes MI, ROI = lung⊕25 mm): 01BB rot [6.4,1.7,−1.5]° trans [15,13,11] mm, gas-in-CT-lung
  0.96 / lung covered 0.38 (bad lung, apices unventilated); 02ZS rot [3.8,−0.7,3.0]° trans [17,−3,3] mm,
  0.79/0.77 (left lung largely unventilated, gas slightly below CT diaphragm at right base = inflation
  mismatch, rigid cannot fix); 03PM rot ≤1.7° trans [−8,0,−2] mm, 0.93/0.74 — trachea/carina line up in
  sagittal. Similarity-scale diag 0.981/1.000/0.968 ⇒ FOV 350 holds. Dead ends: (1) sitk registration
  returns CompositeTransform — downcast via GetNthTransform(0).Downcast(); (2) lung mask by
  border-touching-air deleted lungs when pharynx/mouth in FOV; (3) area-threshold slice range let sinus
  air in — use contiguous run around max-area slice. Figures: outputs/ct_overlay/<subj>/; big files on Ext.
  Reply sent to PCA session with F58.
- 04:55 CT-phase check (C34 `phases`, phase_scores.csv): EI ↔ CT phase chosen by MAX CT lung volume (not the 0%
  label). 03PM phase00 = max (2.78 L, MI best) ✓; 01BB phase00 = max (4.62 L, MI flat ±0.005 across phases) ✓;
  02ZS: 0% is NOT end-inhale on this scan (2.18 L) — phase37 = max (2.35 L), MI −0.34 vs −0.30 → 02ZS
  re-registered to phase37 (gas-in-lung 0.79→0.85; base-of-lung spill under CT diaphragm reduced). phase00
  result kept as `outputs/ct_overlay/02ZS_phase00/` + Ext `mri_on_ct_phase00.nii.gz`, `rigid_phase00.tfm`.
  Overlay clip = tissue(HU>−500) ∪ lung⊕15 mm (01BB immobilization bag encloses air gaps, fill-holes body
  mask let noise specks through).
- 05:55 COHORT RESULT (C35): 29 keys (17 LTX subjects, 3 with two 2023 sessions; 9 EBV pre/post) all ran
  (exit 0, ~8 min/key, 3-wide). CT pick = max lung-volume ORIGINAL axial ≤2.5 mm series (all soft-kernel
  1–2 mm); EBV 002JM/007RA/009JT/010AJ used the pre-made NIFTI/CT.nii.gz. Fit table
  `outputs/ct_overlay/cohort_fit_table.csv`; one page per key `cohort_summary.pdf` (53 MB, on Ext, symlinked). Eye review of the
  numeric outliers: 012SR OK (patient tilt −13°); 010AJ OK — gas only in upper lobes (basal non-
  ventilation), z-shift −30 mm follows the trachea; 037GD OK but low SNR (noise blobs, EI bin 2 by count
  — bins 2/7/8 equal by top-signal); 002JM = UNREGISTERED stack (no elastix result), noisy but fits;
  038RL POOR — CT anatomy distorted (mediastinal shift, small high left lung) vs a symmetric MRI, rigid
  cannot reconcile (gas-in-lung 0.64) — needs Hooman's read (different disease state at CT?). Hooman
  approved the 16-bin videos (04:40) before the batch. S1/S2 EBV folders skipped (subject unknown).

## 2026-09-17 · session 8 (Fable 5.1; 03:35 → 19:10)
(the "2026-09-17 ~04:00", RESULT 04:40, 04:55, COHORT RESULT 05:55 bullets above, appended live under the
session-7 header, belong to THIS session — left in place, append-only.)
- CT overlay delivered: 4D-CT (RT) + clinical CT (LTX 2023, EBV pre/post) → NIfTI; rigid-only xenon gas EI
  overlay (Mattes MI, no deformation) for 32 keys; axial all-slice + 1 mm coronal montages; 16-bin videos.
- Tyger grid ↔ scanner orientation established (F58) by 48-combo brute force, unanimous over 3 subjects;
  PCA session recorded it as their F63.
- CT inhale phase/series chosen by max lung-mask volume (02ZS 4D-CT 0% ≠ end-inhale → phase37).
- Eye review of numeric outliers: 038RL poor (distorted CT anatomy); 010AJ/037GD/012SR/002JM OK with caveats.
- Evening: Hooman asked which subject a 3-tile xenon coronal screenshot is. Template matching against the 104
  stacks was NOT discriminative (NCC 0.56–0.64 flat). ScreenPipe placed the screenshot in
  `Workshop_EBV_Dynamic_HH_2025_Finaled.pptx` (Presentations/2025_Workshop), slide text "COPD, Male, 68,
  RUL Blocked (EBV), FEV1 > 32 % improve, EBV-PRE/POST". Deck file not locatable from the sandbox (mdfind/
  find empty); Hooman took over the lookup. Lesson: for "which subject is this figure" go to ScreenPipe
  first, image matching second.
- DECIDED (Hooman): 16-bin video format approved (6 axial + 6 coronal tiles, EI-fixed colour scale) — then
  batch for all.
- DECIDED: rigid only; CT phase/series = max lung volume; big outputs on Ext mirror path, montages in repo,
  53 MB cohort_summary.pdf on Ext + symlink.
- BRANCH: ct-overlay (B21), parent: B14. Open (038RL read, S1/S2 EBV subject ids, workshop-slide subject
  pending Hooman). ★ → B21.

## 2026-09-24 · session 9 (Fable 5.1; ~18:30 → 20:10)
- Hooman asked for a critical read of the 129Xe 1-point Dixon literature (Duke canon, consortium, critics) +
  an opinion on what the field converged on. Delivered `reference/OnePointDixon_Review.md` (tags learn,
  literature): 14 papers read in full by the session, ~25 secondary papers extracted by two sub-agents
  (dossiers condensed into §8), Duke public pipeline `xenon-gas-exchange-consortium` decomposition code
  inspected, Steve's `raw.py`/`results.py` phase model re-read against it.
- Verdict (short): convergence is on acquisition + reporting conventions, not on validation of the
  RBC/membrane decomposition. Imaging RBC:M equals the spectroscopic ratio BY CONSTRUCTION in 1-point
  Dixon, so "agrees with spectroscopy" is not evidence. Unclosed physics: RF/exchange phase intercept
  (TE90 ≠ 1/4Δf, disease-dependent), phase evolution across the readout (consortium: "not rigorously
  examined"), third dissolved resonance (~100° from barrier-1). Confounders (lung volume r≈−0.97, sex
  ±15–20 %, Hb ≤20–40 %, age −0.05 RBC:M/decade) exceed site-to-site scatter. Robust outputs = total
  dissolved, spectroscopic RBC:M, RBC shift, oscillation amplitude.
- SUSPECT (F59, not yet tested): Steve's `TEeff` (raw.py:389) is the phase-vs-TE INTERCEPT converted to
  time, minus/plus killpts terms; the 2πΔf·TE term of the real echo (`self.TE`) is absent, yet
  results.py:324 uses 2π·f·TEeff as the split basis. At 1.5 T (F54) that is ≈0.12°/µs of TE. If the
  fitted `RBCphase[0]−TPphase[0]` differs from 2πΔf·TEeff by ≈2πΔf·TE, the near-singular split in B17
  (Δφ 6–39°) is a code artefact, not an acquisition property. One-line test per session from the cal block.
- Noted for later: our 1.5 T (T2*≈2 ms) is the friendlier field for a Collier-style multi-echo
  (ΔTE 0.7 ms) spiral decomposition on ASAP — no ratio prior, contamination modelled not subtracted.
- Duke pipeline facts recorded: `dixon_decomposition` = global rotation to arctan2(RBC:M,1) + gas-image
  phase map subtraction; Hb/volume/T2* corrections are global scalars; contamination phase hard-coded
  (`optimized_conta_phase = 49.9`).
- Lane note: no scripts written; no root writes. Vault topic note `Literature/Topics/Dixon Methods in HP
  Gas MRI.md` NOT edited (Hooman's/Librarian's lane) — a "see also" pointer to the review is suggested.
- BRANCH: dixon-lit-review (B22), parent: B17 (rbctp-split-audit). Open (F59 test pending). ★ stays B21.
- 20:40 F59 VERIFIED (SUSPECT → VALID): ran Steve's spectral fit locally (CPU) on 041WF/004DS/043AS via
  `helpers/_f59_check.py` (scratch tier). `TEeff` encodes the phase difference EXTRAPOLATED TO t=0
  (intercept 32°/43°, the Kaushik-2016 "RF intercept"); the split then uses Δφ_model = 20.5°/32.3°/≈9°
  (noise gain 2.9/1.9/6.7). The phase difference actually measured at the first spectral sample
  (t = TE + 2·60 µs) is 95.6°/102.6°/≈88° — the gap equals 2πΔf·TE (74°/66°/77°). So our TE = 0.62 ms
  puts RBC/TP ≈ in quadrature already; Steve's "it's just math, no TE90 needed" is right in principle,
  the code just uses the wrong time point. One-line fix in results.py:324–325 (use RBCphase[0]−TPphase[0]).
  043AS shows phase wrap (−632°) in the polyfit path → unwrap/complex-fit needed. Also traj.BW hard-coded
  10 µs vs header dtdyn 5 µs (only the killpts term, ≈1°). Root untouched; goes to Steve via 2steve note
  (to write: 05_RBCTP_Split_TE_Term.md) when Hooman says so.
- 21:05 Wrote `2steve/05_RBCTP_Split_TE_Term.md` (+ README index row, figs `2steve/fig/05_*.png`,
  `05b_*_218ppm_badfits.png`; data `outputs/f59_te_term_2026-09-24/`). Extra finding while making the
  figure: on 218-ppm (RBC-centred) sessions 2023-04-07 030DN, 2023-05-04 034LR, 2023-08-24 030DN the
  2-peak picker locks onto wrong peaks (Δf = 4403 / 991 / 1460 Hz vs the physical ~330 Hz; ratios
  0.53 / 0.66 / 1.64 meaningless) → note 05 asks for a Δf/fRBC gate. Not sent to Steve (Hooman decides).
- 21:40 Wrote `2steve/06_Calcb_Static_Structure.md` (+ README row; figs `2steve/fig/06_*.png`; data
  `outputs/calcb_imprint_2026-09-24/`). Mechanism from code: one cycle-average b (10σ mask + global
  quadratic outside; nch=8 → image-derived magnitude Σ|b| ∈ [1,8]) multiplies all 16 bins, gas takes
  real(). Saved-recon check (real images): negative lung-mask voxels 8.4→0→7.5 % across the SIGNAL cycle
  (041WF 2024-09-03), deep (<−3σ) ≤0.4 % → attenuation at the fixed average-lung boundary, not
  inversion. 13 sessions are 8-channel (list in note). Σ|b| map for 042DR being computed with the XeCS
  numpy replica at MS=104 (`helpers/_calcb_8ch_weight.py`, background) — append numbers to note §C.
  Our solution stated in the note: hybrid_b + imag_qc (XeCS, note 01), magnitude for videos, Σ|b|
  normalisation for 8-ch, EE∪EI navigator mask for b's support.
- 22:05 042DR (8 ch) Σ|b| map done (numpy replica, MS=104, ~20 min): W inside lung min/p5/median/p95/max =
  1.97/2.65/3.66/5.26/6.95, background median 4.51 → the 8-channel bins carry a fixed ×2 (p5–p95)
  lung-shaped gain, bright apical/peripheral rim, darker core, noise weighted above lung. Appended to
  `2steve/06` §C; fig `2steve/fig/06_calcb_W_map_042DR_8ch.png`. Notes 05+06 complete, unsent.
- 22:20 /leave. Session scope was literature + diagnosis; no root code changed, nothing sent to Steve.
- DECIDED (Hooman): 8-channel Σ|b| normalisation deferred — single-channel fixes first.
- DECIDED (Hooman): implement notes 05 + 06 in a NEW session (context hygiene); order given via /handoff.
- DECIDED: dissolved RBC/TP re-split is doable OFFLINE from saved aRBC/aTP (invertible linear map) —
  no Tyger re-run needed for the 86 sessions; calcb imprint is NOT fixable offline (real() already applied).
- BRANCH: dixon-lit-review (B22), parent: B17. Closed (delivered: reference/OnePointDixon_Review.md,
  2steve/05, 2steve/06, F59 VALID, F60). ★ → B17 (rbctp-split-audit reopened: resplit next). B21 stays open.
