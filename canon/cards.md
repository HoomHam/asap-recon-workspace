# ASAP Recon — Cards

One card per load-bearing script. Origin: H(uman)/C(laude)/A(gent)/R(etro-inferred).
Retro-filled 2026-07-12; all origins R. Quick table regenerated from bodies by /leave.

## Quick table

| ID | Script | Branch | Status |
|----|--------|--------|--------|
| C1 | main.py (root) | stem | WORKS (local, 2 uncommitted fixes) |
| C2 | raw.py (root) | stem | FIXED in fork 472fbc9 (fit fallback, dp optional); Bug B (unbounded LM) open |
| C3 | recon.py (root) | stem | WORKS (numba pinned 0.65.1, F42) |
| C4 | results.py (root) | stem | WORKS · ⚠ RBC/TP split unreliable (F47) · 40d23a4 keeps unsplit dissolved |
| C5 | gtypes.py (root) | stem | WORKS |
| C6 | convert_siemens_to_mrd.py (root) | tyger | WORKS (F16 caveat) |
| C7 | tyger_recon.py (root) | diaphragm | WORKS (fork image 40d23a4, cloud GPU) · z-flip = deviation from Steve (F49) · nav 26-guard local |
| C8 | pipeline/batch_recon.py | batch | SUPERSEDED by C29 todo_batch.sh (F1 retracted) |
| C9 | pipeline/dyn_recon.py | batch | WORKS all nch on 40d23a4 (F42); robust submit (F43 F44), MPLBACKEND=Agg (F45) |
| C10 | pipeline/asap_run.py | auto-steve | WORKS |
| C11 | pipeline/post_process.py | auto-steve | WORKS · writes rbc_tp_separated to recon.mat |
| C12 | pipeline/param_gui.py | auto-steve | WORKS |
| C13 | pipeline/recon_codespec*.yml | tyger | CONFIG — use recon_codespec_40d23a4.yml (F42); base file still d136eb1 |
| C14 | helpers/build_status_tab.py | batch | WORKS |
| C15 | helpers/probe_nch.py | batch | WORKS (nch=8 only, F8) |
| C16 | helpers/check_frbc.py | batch | WORKS (F8 caveat) |
| C17 | helpers/diag_fit.py | batch | WORKS |
| C18 | helpers/recon/compare_baseline.py | arbiter | WORKS (blocked: missing .dat+traj) |
| C19 | helpers/recon/steve_kernel_numpy.py | arbiter | WORKS (~1e-6 vs GPU) |
| C20 | helpers/recon/dump_inputs{,_dyn}.py | arbiter | WORKS |
| C21 | helpers/recon/convert_calib.py | arbiter | WORKS |
| C22 | helpers/recon/faraz_montage.py + faraz_zoom_check.py | faraz | WORKS |
| C23 | helpers/recon/cs_montage.py | xecs-bridge | WORKS (coupled to XeCS sweep json) |
| C24 | helpers/recon/diaphragm_bin_demo.py | diaphragm-binning | WORKS (synthetic demo) |
| C25 | helpers/recon/diaphragm_bintime_demo.py | diaphragm-binning | WORKS (synthetic demo) |
| C26 | helpers/atlas/dis_atlas.py | aikill-atlas | gp pages fine; ⚠ dp pages SUSPECT (F47) — SUPERSEDED by C39 rbc_tm_atlas.py 2026-09-25 (still imported for layout helpers) |
| C27 | helpers/rt_prepost_fig.py | main | WORKS (RT pre/post figs 2026-08-31) |
| C28 | helpers/snr_calc.py | snr-table | WORKS (gas + whitened DP SNR, 104 sessions, F48) |
| C29 | pipeline/todo_batch.sh | leftover-recon | WORKS (25-session worklist, FORCE=1) |
| C30 | helpers/structure_rank.py | structure-rank | WORKS (107 sessions incl. 3 merged; carina-slice panels as mp4/gif over 16 bins, F57) |
| C31 | pipeline/merge_dyn.py | merged-dyn | WORKS (030DN validated through raw loader, F56) |
| C32 | pipeline/merged_batch.sh | merged-dyn | DONE 2026-09-16 (4 run, 3 adopted: 030DN 013VM 008CR; 025VP REJECTED — weak 2nd dose) |
| C33 | helpers/ct/ct_to_nifti.py | ct-overlay | WORKS (RT 4D-CT → 27 NIfTI phases on Ext) |
| C34 | helpers/ct/ct_mri_overlay.py | ct-overlay | WORKS (search/refine/overlay/video/phases; F58 orientation; rigid only) |
| C35 | helpers/ct/ct_cohort_batch.py | ct-overlay | DONE 2026-09-17 (29 LTX+EBV keys; 038RL poor fit; S1/S2 skipped) |
| C36 | helpers/resplit_rbctp.py | rbctp-resplit (B23) | WORKS 2026-09-24 (--fit steve: 79 re-split, 10 gated; --fit xecs: 68 re-split, 6 gated, 15 no cal block; stability gate) |
| C37 | helpers/specfit_template.py | rbctp-resplit (B23) | DONE 2026-09-25 — precision fix, accuracy unproven (F71); fallback-only rule proposed, not deployed |
| C38 | helpers/te90_table.py | aikill-atlas (B12) | DONE 2026-09-25 — TE90 per session from the XeCS fits, 88 rows + strip plot (F72) |
| C39 | helpers/atlas/rbc_tm_atlas.py | aikill-atlas (B12) | WORKS 2026-09-25 — gas/RBC/TM atlas; `--root/--tag` = production tree (container split); outputs aikill_atlas/ (raw b, superseded) and aikill_atlas_b44/ (F82) |
| C41 | helpers/atlas/perbin_phaseref.py | perbin-phaseref (B24) | DONE 2026-09-25 — per-bin gas phase reference test on fork f5ac7c5 outputs: no effect (F77), F74 retracted |
| C42 | helpers/atlas/rbc_repro_tests.py | perbin-phaseref (B24) | DONE 2026-09-25 — mirror-bin / cross-subject / cross-date reproducibility of aRBC structure (F75 F76) |
| C43 | helpers/atlas/bphase_fix.py | perbin-phaseref (B24) | WORKS 2026-09-25 — proxy removal of the static fine-scale b phase from the split (F79); 050EH/045VS; --write stores recon_resplit_xecs_bfix.mat |
| C44 | helpers/atlas/bphase_confirm.py | perbin-phaseref (B24) | DONE 2026-09-25 — F79 confirmed against the exported b (Tyger 649/650, fork c8366c3); true-b fix vs proxy |
| C45 | helpers/atlas/bphase_sigma.py | perbin-phaseref (B24) | DONE 2026-09-25 — how much of b's phase to keep: band-resolved static residual vs σ, |dis| floor; σ 2.5 chosen |
| C46 | helpers/atlas/bsmooth_compare.py | perbin-phaseref (B24) | DONE 2026-09-25 — container 19b7365 (smoothed b) vs raw container vs offline true-b: match 0.98, gas −0.7 % rms (F80) |
| C47 | helpers/atlas/bsmooth_videos.py | perbin-phaseref (B24) | DONE 2026-09-25 — raw-b vs smoothed-b atlas videos (3 sessions); inline variant: raw / σ2.5 / σ10 / no-b rows (videos_sigma/) |
| C48 | helpers/atlas/bphase_sigma_cohort.py | perbin-phaseref (B24) | DONE 2026-09-25 — 57 sessions (+4 stragglers): σ = 5 chosen (F81); fallback for container-gated sessions |
| C49 | pipeline/run_prod_b44.sh + prod_b44_todo.txt | prod-b44 (B25) | DONE 2026-09-25 — 61/61, σ 4.4 (fork 04f445b) → Ext AIkill_Dynamic_b44/; 59 split, 2 container-gated (F82) |
| C50 | pipeline/mrd_to_mat.py | prod-b44 (B25) | WORKS 2026-09-25 — output.mrd → recon.mat (gas, |gas|, aRBC/aTP, |dis|, calcb_b, rbc_tp_* meta, nav), no figures |
| C40 | helpers/atlas/split_trust.py | aikill-atlas (B12) | DONE 2026-09-25 — trust grade; `--root/--tag` = production tree → split_trust_b44/ (A18/B23/C13/D5, F82) |

## Cards

### C1 · main.py (root)
- why: local tkinter GUI entry; canonical source of DIAPHRAGM navigator (`calcLVcb`, :279-330) and dataset→traj resolution (`ID_callback`, :390-442)
- origin: R · branch: stem · facts: F31
- in: Siemens .dat + traj .npy via basefolder/datatype/date/subjectID/
- out: on-screen recon (needs CUDA — fails on M4)
- status: WORKS locally, two uncommitted fixes (path `/`, hidden-file filter)

### C2 · raw.py (root)
- why: TWIX load, traj load, binning, two-Lorentzian spectral fit (lorfit :13, loop :313-380)
- origin: R · branch: stem · facts: F4 F5 F9 F13 F18 F42 F46
- in: .dat/MRD arrays, trajectory
- out: k-space, ilvbin, fRBC/fTP ratios
- status: fork 472fbc9 (2026-09-12): non-converged Lorentzian curve_fit → clear fit lists (gas-only / unsplit dissolved) instead of abort; dp trajectory optional (gas-only v3_20230821). Bug B (unbounded warm-start LM, non-physical params) still open. FFT pattern classifier correctly labels 2024-01 8-ch scans gas-only (F46).

### C3 · recon.py (root)
- why: CUDA gridding kernels (cudarecon/cudarenorm) via numba.cuda
- origin: R · branch: stem · facts: F1(retracted) F17 F42
- in: k-space samples, traj, b-matrix → out: gridded volume
- status: WORKS — numba 0.65.1/llvmlite 0.47.0 pinned in requirements.txt (F42); identical to Steve's original

### C4 · results.py (root)
- why: image recon, calcb (B0/coil-phase b-matrix), dyn_recon/dyn_usimg_recon, RBC/TP split (:298-327 Steve original)
- origin: R · branch: stem · facts: F4 F5 F47 F48
- status: WORKS for gas. ⚠ Steve's RBC/TP phase sweep stops at the ΣaTP=0 discontinuity and the one-point Dixon is ill-conditioned (F47) — stored aRBC/aTP invalid in 60/86 split sessions. Fork 40d23a4: when fRBC empty, dissolved kept as unsplit complex + `rbc_tp_separated` attribute (was: July guard dropped DPDYN). Suggested Steve fix: `2steve/04_RBCTP_Split_PhaseStop_Conditioning.md`

### C5 · gtypes.py (root)
- why: global enums (bintype, imgtype) + gvar defaults; basefolder :27 Hooman-specific, intentionally uncommitted
- origin: R · branch: stem · facts: F31
- status: WORKS

### C6 · convert_siemens_to_mrd.py (root)
- why: local .dat → MRD; bakes recon params into MRD header (only channel to headless cloud job)
- origin: R · branch: tyger · facts: F16 F49
- in: Siemens .dat, --binning → out: input.mrd
- status: WORKS (breaks on `._*.dat` shadows; dyn_recon stages symlinks). Dynamic array bit-identical to Steve's mapVBVD read (030DN, F49)

### C7 · tyger_recon.py (root)
- why: Tyger cloud-GPU Docker entrypoint; hosts ported _diaphragm_navigator + MRD nav-array export
- origin: R · branch: diaphragm · facts: F2 F30 F42 F49
- in: input.mrd (params :66-83) → out: output.mrd (dissolved tagged `rbc_tp_separated` since 40d23a4)
- status: WORKS on fork image 40d23a4; NVIDIA GPU (cloud) only. dp trajectory optional (472fbc9). DIAPHRAGM z-flip (db80f16) is a real deviation from Steve's GUI navigator (F49). Hooman's Thomson-26 diagnostic still uncommitted (local).

### C8 · workspace/pipeline/batch_recon.py
- why: batch driver — iterate subjects, call dyn_recon.py --methods s,p,d, green-mark Excel row, save incrementally
- origin: R · branch: batch · facts: F1(retracted) F10 F15
- in: SNR_Table_All.xlsx, source/archive drives
- out: /Volumes/HoomHamExt/Dynamic/<date>_<id>/{s,p,d}/ (now AIkill_Dynamic), updated Excel
- status: SUPERSEDED for the leftover set by C29 todo_batch.sh; default drive paths stale (HoomHam dead)

### C9 · workspace/pipeline/dyn_recon.py
- why: per-subject driver — resolve, _clean_datadir (symlinks around ._*), convert, submit to Tyger, publish, post-process
- origin: R→C (reworked 2026-09-12) · branch: batch · facts: F16 F42 F43 F44 F45
- in: `/Volumes/HoomHamExt/_5t_images_roundtrip/Images/<date>/<id>/` (+ `--data-root`), `--codespec`
- out: `/Volumes/HoomHamExt/AIkill_Dynamic/<date>_<id>/{s,p,d}/` (output.mrd, input.mrd, tyger.log with run/buffer ids, codespec.yml, pngs, recon.mat, fig/); scratch run dirs on Ext `.../2026_ASAP_Recon/pipeline_runs/`
- status: WORKS for nch=1 and nch=8 on image 40d23a4. submit = buffer create (NO --ttl) → buffer write → run create → poll → run logs → buffer read -o -p 4, retried; plot/post under MPLBACKEND=Agg

### C10 · workspace/pipeline/asap_run.py
- why: single-dataset orchestrator, 6 stages resolve→param GUI→convert→submit→publish→post_process ("analyze 25JC with Steve")
- origin: R · branch: auto-steve
- out: workspace/outputs/<dataset>/<timestamp>/
- status: WORKS

### C11 · workspace/pipeline/post_process.py
- why: output.mrd → recon.mat, signal_pneumo.npz, per-bin slice videos, navigator.gif, resp_traces.png
- origin: R→C · branch: auto-steve · facts: F47
- out: recon.mat now also carries `rbc_tp_separated` (1 split, 0 unsplit, -1 untagged pre-40d23a4)
- status: WORKS (resp_traces nav-z scatter in last ~15% un-cleaned; July pneumotach drift fix included). `dissolved_phase_real/imag/magnitude` are Steve's aRBC/aTP — not trustworthy for split sessions (F47)

### C12 · workspace/pipeline/param_gui.py
- why: standalone tkinter param picker (mirrors main.py fields incl. DIAPHRAGM) → params.json
- origin: R · branch: auto-steve
- status: WORKS

### C13 · workspace/pipeline/recon_codespec*.yml
- why: Tyger job spec (image sha); variants `recon_codespec_db80f16.yml`, `_472fbc9.yml`, `_40d23a4.yml`
- origin: R→C · branch: tyger · facts: F2 F42
- status: CONFIG — current = `recon_codespec_40d23a4.yml`; base `recon_codespec.yml` still points at d136eb1 (July flip, uncommitted, not this session's)

### C14 · workspace/helpers/build_status_tab.py
- why: regenerate Excel "Recon Status" tab from drives (blocker + HYPERLINK per ungreened subject)
- origin: R · branch: batch · facts: F15
- status: WORKS (conda base python, PYTHONPATH=.) — tab now stale vs 2026-09-13 greens

### C15 · workspace/helpers/probe_nch.py
- why: CPU-convert each subject to MRD, report nch/numspec/fRBC, predict guard; submits nothing
- origin: R · branch: batch · facts: F8
- status: WORKS — faithful ONLY for nch=8

### C16 · workspace/helpers/check_frbc.py
- why: guard predictor (fRBC presence)
- origin: R · branch: batch · facts: F8
- status: WORKS, same nch=1 caveat

### C17 · workspace/helpers/diag_fit.py
- why: monkeypatch raw.curve_fit, swallow failures, print per-iTE SNR/peak-sep/linewidths — the tool that found Bug B
- origin: R · branch: batch · facts: F5 F34
- status: WORKS (CPU, conda base)

### C18 · workspace/helpers/recon/compare_baseline.py
- why: arbiter — FINUFFT recon vs Steve baseline on real phantom data
- origin: R · branch: arbiter · facts: F21 F23
- in: recon_io/ (.dat + gas traj .npy)
- status: WORKS, blocked on missing .dat+traj (open loop B1)

### C19 · workspace/helpers/recon/steve_kernel_numpy.py
- why: faithful CPU float64 reimplementation of Steve's cudarecon/cudarenorm
- origin: R · branch: arbiter
- status: WORKS (~1e-6 vs GPU; Colab certification optional/pending)

### C20 · workspace/helpers/recon/dump_inputs.py + dump_inputs_dyn.py
- why: dump Steve-format .npy inputs from any .dat via his own loaders (CPU); _dyn for dynamic
- origin: R · branch: arbiter
- status: WORKS (--fov flag; traj class hardcodes 350)

### C21 · workspace/helpers/recon/convert_calib.py
- why: cal-struct .mat → Steve-format trajectory .npy (k in 1/mm)
- origin: R · branch: arbiter
- status: WORKS

### C22 · workspace/helpers/recon/faraz_montage.py + faraz_zoom_check.py
- why: Faraz comparison figures + verification of ×1.205 zoom bug
- origin: R · branch: faraz · facts: F22
- status: WORKS (phase_corrected_real kept as intentionally-unused cautionary trap)

### C23 · workspace/helpers/recon/cs_montage.py
- why: CS-vs-Steve-vs-Faraz montage (ASAP side); imports CS operators from XeCS via xecs_recon.pth; reads maxeig/t_ref from XeCS cs_sweep_metrics.json
- origin: R · branch: xecs-bridge · facts: F35
- status: WORKS but coupled — rerunning XeCS sweep silently changes later montages

## Scratch / reference (uncarded)
`read_mapvbvd.py` (root, standalone loader, UNKNOWN) · `workspace/codes/kasap.c` (Kento reference, F20) · `asap/asap.c` (Steve reference, F20) · `helpers/_delete/` (byte-verified CS originals, moved to XeCS 2026-06-24) · `helpers/calib/` (duplicate .npy pair) · `pipeline/runs/`, `batch_recon.log`, `__pycache__/` · 2026-09 scratch tier: `helpers/_roundtrip_audit.py` (roundtrip image-tree audit → outputs/roundtrip_audit_2026-09-13/), `helpers/_zorder_check.py` (F49 check → outputs/zorder_check_2026-09-12/), `helpers/_rbctp_fig.py` (F47 evidence figure → outputs/snr_2026-09-13/04_rbctp_split_evidence.png), `helpers/_f59_check.py` (reruns Steve's spectral fit on d/input.mrd, compares split basis Δφ vs measured → outputs/f59_te_term_2026-09-24/, F59), `helpers/_calcb_8ch_weight.py` (calcb Σ|b| map via XeCS numpy replica at MS=104 → outputs/calcb_imprint_2026-09-24/) · `Codes/2026_Steve_Recon/` (plain snapshot of Steve main 3303276, diff reference only).

Note: root CLAUDE.md still lists cs_recon.py / cs_recon_4d.py under helpers/recon — STALE, they moved to 2026_XeCS_Recon in the 2026-06-24 decouple.

### C24 · helpers/recon/diaphragm_bin_demo.py
- why: reproduce Steve's `raw.py:153 bin()` verbatim on a synthetic navigator curve; 4-panel
  figure showing amplitude-rank → representative breath/stroke cycle (limb split, fold, equal-count)
- origin: A · branch: diaphragm-binning · facts: F37, F38
- in: none (self-contained synthetic navigator)
- out: workspace/outputs/diaphragm_binning/diaphragm_binning.png (A: z(t) colored by slope;
  B: the fold; C: representative cycle bin-vs-amplitude; D: equal-count + breath pooling)
- see: bin_time_PLAN.md (same dir)

### C25 · helpers/recon/diaphragm_bintime_demo.py
- why: derive the REAL per-bin time vector for the position-rank 4D stack — EE-detect →
  within-cycle τ → soft-weighted CIRCULAR mean per bin (matches recon.py gridder weights)
- origin: A · branch: diaphragm-binning · facts: F37, F38
- in: none (synthetic; swap real nav_volume/ilvtime to use)
- out: workspace/outputs/diaphragm_binning/diaphragm_bintime.png (C: bin→time monotonic but
  non-uniform; D: real time gaps vs uniform — clusters at turns). Console prints per-bin phase+sec.
- note: plan + recipe in workspace/helpers/recon/bin_time_PLAN.md; bin_time.py NOT yet built

### C26 · helpers/atlas/dis_atlas.py
- why: subject-level QC atlas of the AIkill_Dynamic batch recons — group 84 sessions by
  subject ID (dates sorted), gas (gp) + dissolved (dp) at the highest gas-signal bin
- origin: H · branch: main · facts: F46 F47
- in: /Volumes/HoomHamExt/AIkill_Dynamic/*/d/recon.mat (gas_phase +
  dissolved_phase_magnitude, (16,Z,Y,X)); s/recon.mat gas fallback for the 3 gas-only sessions
  (2024-01-18_001JM, 01-22_007RA, 01-31_008CR — gas-only ACQUISITIONS, F46, not skipped dissolved)
- out: workspace/outputs/aikill_atlas/dissolved_atlas.pdf (per ID: gp panel page + dp panel
  page [orientation-major blocks cor/sag/ax, per date one 10-slice row], then per date 6
  full all-slice pages cor/sag/ax × gp/dp) + videos/<ID>_{gp,dp}.mp4 (16 bins ×5, 5 fps)
- note: orientation follows pipeline/post_process.py (axial flipud, coronal reversed,
  sagittal rot90 ccw); slice extent from GAS mask (≥30 voxels above 0.15·max, bin-mean);
  cross-date column alignment = per-session slices at shared apex→base fractions + fixed-size
  crop centered on each session's lung bbox (shift+scale, no registration); black bg, no gaps
- status: ⚠ dp pages/videos SUSPECT (F47): `dissolved_phase_magnitude` = |aRBC+i·aTP| is noise-amplified by Steve's split for split sessions. Header stamped; marker `outputs/aikill_atlas/SUSPECT_F47.md`. Gas pages fine. Covers 84 sessions (pre-2026-09-12 cohort).

### C27 · helpers/rt_prepost_fig.py
- why: RT-study figure — per subject, 10 coronal lung slices per visit, pre-RT row over
  post-RT row, for the 9 RT subjects in SNR_Table_All (Hooman ask 2026-08-31)
- origin: H · branch: main · facts: —
- in: /Volumes/HoomHamExt/AIkill_Dynamic/<date>_<id>/s/recon.mat (gas_phase (16,Z,Y,X));
  visit list + SNR hardcoded from workspace/data/SNR_Table_All.xlsx Study=RT rows
- out: workspace/outputs/rt_prepost/rt_<ID>.png (per subject), rt_all_pairs.png
  (002ZS/003PM/004DS stacked), rt_all_singles.png (001BB/005JJ/006MM/008TP/009ML)
- note: EI bin (brightest = end-inspiration), per-visit norm to 99.5 pct of lung voxels,
  slice picks = 10 even fractions over lung Y-extent trimmed 8% each end, union Z-X crop
  across a subject's visits; coronal orientation matches post_process.py. Excluded: 007IT
  (container curve_fit fail, no recon — now reconstructed 2026-09-12), 001BB pre 2023-03-27
  (not on drive; its post folder is named 002BB — ID mismatch vs Excel 001BB)

### C28 · helpers/snr_calc.py
- why: own image-domain SNR for every AIkill_Dynamic session — gas + total dissolved, per bin,
  with max / min / end-inspiration summaries for SNR_Table_All.xlsx (Hooman ask 2026-09-13)
- origin: C (method approved by H) · branch: snr-table · facts: F47 F48 F50 F52
- in: `/Volumes/HoomHamExt/AIkill_Dynamic/<date>_<id>/<binning>/recon.mat` (default `d`: gas_phase,
  dissolved_phase_real/imag, rbc_tp_separated); `--only`, `--binning`, `--out`
- out: `workspace/outputs/snr_2026-09-13/snr_table.csv` — one row per session: EI/exp bins
  (max/min lung-mask volume), gas/DP SNR at EI/exp/max/min + bins, per-bin SNR and lung-voxel
  strings, σ, DP method; diagnostic RBC/TP per-part SNR (unreliable). Excel cols I–S written from it.
- method: σ = std in 8 corner 10³ cubes minus 3-voxel-dilated lung; per-bin mask gas>5σ,
  blobs ≥500 vox; DP = √max(vᵀC⁻¹v − 2, 0) mean over mask (C = background covariance of
  real/imag; E[q]=2 verified)
- status: WORKS (104 sessions, 0 errors, ~1 s/session)

### C29 · pipeline/todo_batch.sh
- why: detached worklist driver for the 25 leftover v2/v3 sessions (inventory 2026-09-12)
- origin: C · branch: leftover-recon · facts: F42 F45
- in: `<codespec.yml> [filter-regex]`, env `FORCE=1`; roundtrip Images + Ext `staged_src/` symlink dirs
- out: per-session logs + `SUMMARY.txt` in `/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/pipeline_runs/logs/`; outputs via C9
- status: WORKS (run with nohup; bigmac 2023-11-02 000LL commented out — duplicate)

### C30 · helpers/structure_rank.py
- why: Hooman wanted the recons with the most intra-lung patchiness (COPD-like bright/dark patches with
  edges) at good SNR, size-normalised, with one panel of coronal trachea slices so his eye picks
- origin: C · branch: structure-rank (B19) · facts: F57
- in: `AIkill_Dynamic/<session>/d/recon.mat` gas_phase; `outputs/snr_2026-09-13/snr_table.csv` (max-SNR bin);
  `--only`, `--sort <col>`, `--figs-only` (rebuild panels from csv + tiles.npz without recompute)
- out: `outputs/structure_rank/{structure_table.csv, panel_bifurcation_slice.{png,mp4,gif},
  panel_maxarea_slice.{png,mp4,gif}, *_frames/, scatter_snr_vs_G.png, tiles.npz (bins,z,x per session)}`;
  `--figs-only --sort <col> --fps N --no-video` rebuild without recompute; per-cohort panels
  (`cohorts/panel_<HC|HC_OLD|EBV|LTX|RT>_{bifurcation,maxarea}.*`, subject blocks date-ordered) from
  `--cohort-map notes/calspec_package_2026-09-16/cohort_map.csv`, `--cohorts-only`
- method: per coronal slice at max-SNR bin: extent = 5σ mask closed r=8 per midline side + fill; rim weight
  0→1 over 3–8 px; R = I/localmean(σ10); G = w-mean |∇ gaussian(R, σ3)| − background floor; P = perimeter
  (bright)/perimeter(extent); session = area-weighted mean over slices. Carina slice (2026-09-17) = Y
  template (trachea column × diagonal arms, anchored at apex/midline) on max-over-bins volume, slices
  ≥30 % max area; * = column < 5σ. Videos: one frame per bin, fixed grey per session, ffmpeg mp4+gif
- status: WORKS (103/104; 002ZS has no recon.mat). Dead ends in docstring. Needs scikit-image (venv).

### C31 · pipeline/merge_dyn.py
- why: sessions where the sequence was stopped and restarted have a second good free-breathing dynamic;
  Hooman wants both doses in ONE binned recon
- origin: C (Hooman's decision via XeCS session) · branch: merged-dyn (B20) · facts: F56
- in: `--data-dir Images/<date>/<id> --mids MIDa,MIDb[,…] --gp-traj --dp-traj [--ref MID] --binning S|D --out input.mrd`
- out: input.mrd (MEDCAP MRD via root converter, monkeypatched reader — root untouched), merge_report.json
  (lines, numspec, MDH t0/t1, gap, trim), merge_k0.png (|k0| gas across seam)
- method: A keeps cal block, trimmed to whole arm cycles (nuniquesmp/npts gas ilv); B's cal block stripped;
  concat along lines; header equality asserted (TR/TE/dwell/numspec)
- status: WORKS on 030DN (27460 lines = 229 s imaging). p binning unsupported.

### C32 · pipeline/merged_batch.sh
- why: run the 4 merged sessions unattended (harness kills long foreground jobs)
- origin: C · branch: merged-dyn (B20) · facts: F56
- in: WORK list inside (date id mids ref); `SPEC` env (default recon_codespec_40d23a4.yml)
- out: `AIkill_Dynamic/<date>_<id>_merged/{s,d}`; logs `pipeline_runs/logs/<session>_merged.log`, SUMMARY
- status: launched 2026-09-16 20:55 via nohup

### C33 · helpers/ct/ct_to_nifti.py
- why: Hooman has RT-planning 4D-CT for 01BB/02ZS/03PM and wants the xenon MRI overlaid on it; step 1 =
  DICOM → one NIfTI per respiratory phase
- origin: H (2026-09-17 order) · branch: ct-overlay (B21) · facts: —
- in: `/Volumes/HoomHamExt/Work/CT/Data/<subj>/CT/<series>/*.dcm` (Siemens RT chest, 512², 0.977 mm, 3 mm, HFS)
- out: `/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/ct_overlay/<subj>/ct/{phase00..phase87,average,pbv}.nii.gz`
  + index.json (HU int16; SimpleITK series reader, LPS geometry in header)
- status: WORKS — 27 volumes, 819 MB (Ext, mirror path per big-output rule)

### C34 · helpers/ct/ct_mri_overlay.py
- why: rigid (NO deformation — Hooman: diseased lungs have missing regions, never stretch to fill CT lung)
  overlay of registered xenon gas EI bin on CT phase00 (end-inhale); axial all slices + interpolated coronal
- origin: H · branch: ct-overlay (B21) · facts: (pending — orientation of Tyger (Z,Y,X) grid vs LPS)
- in: PCA-session stacks (01BB/02ZS elastix one-shot `registered_elastix.npy`, 03PM legacy
  `Analysis/2024-12-06_003PM/reg/5000/registered.nii`), EI bins 8/8/6; CT phase00 from C33; FOV 350 mm ⇒ 3.5 mm
- out: Ext `ct_overlay/<subj>/{mri_ei_native,mri_on_ct,ct_lungmask}.nii.gz, rigid.tfm, axial_png/, coronal_png/`;
  repo `outputs/ct_overlay/<subj>/{orientation_scores.csv, register.json, axial_montage.png, axial.gif,
  coronal_montage.png, coronal.gif, summary_3plane.png}`
- method: `search` = 6 axis perms × 8 flips, each rigid Mattes-MI on 3.5 mm CT (ROI = lung ⊕ 25 mm), pick
  lowest MI + cross-subject consensus; `refine` = full-res rigid from it; similarity (7-DOF) run only as a
  scale diagnostic of the 350 mm FOV; CT lung mask = air < −400 HU ∩ per-slice filled body (pharynx-safe)
- status: WORKS 2026-09-17 — same orientation won for all 3 subjects (F58); `video` = 16-bin mp4/gif (all bins through the one rigid transform, 6 axial + 6 coronal tiles); specs also from Ext `ct_overlay/manifest.json` (C35); rigid fits: rot ≤6.4°, trans ≤17 mm,
  scale diag 0.97–1.00. `phases` subcommand = CT-phase match diagnostic (phase_scores.csv); `refine <s> phase=phaseNN` overrides the CT phase (02ZS → phase37, max lung volume). v1 lung-mask bug
  (mouth-connected air deleted lungs) fixed with per-slice body fill.

### C35 · helpers/ct/ct_cohort_batch.py
- why: Hooman: "do it for all" — extend the CT overlay to the EBV (pre/post) and LTX (2023) clinical CTs
- origin: H · branch: ct-overlay (B21) · facts: F58
- in: `Work/CT/CT EBV/*`, `Work/CT/CT LTX/*` (DICOM ORIGINAL axial ≤2.5 mm, or given NIFTI/CT.nii.gz); PCA elastix
  one-shot stacks (fallback unregistered stack.npy); COHORT table inside = CT folder → xenon session(s)
- out: Ext `ct_overlay/<ct_subject>/ct/<series_tag>.nii.gz` + index.json (lung volume per candidate, `_chosen`);
  Ext `ct_overlay/manifest.json` (key `<ct_subject>__<session>`); then C34 outputs per key
- method: `prep` converts every candidate, picks max lung-mask volume (ties → thicker, softer kernel);
  `run <key>` = C34 refine (F58 orientation fixed) + overlay + video; parallel via `keys | xargs -P 3`
- status: DONE 2026-09-17 05:55 — 29/29 keys ran; 1 poor fit (038RL, distorted CT anatomy), 1 unregistered input (002JM);
  S1/S2 EBV folders skipped (subject unknown). cohort_fit_table.csv + cohort_summary.pdf in outputs/ct_overlay/

### C36 · helpers/resplit_rbctp.py
- why: Hooman's directive (session 9 /handoff): re-split the RBC/TP maps of the 86 split sessions offline with the
  k0 basis angle (F59) instead of re-running Tyger; prove the B17 near-singularity was the code artefact
- origin: H (order) / C (design + code) · branch: rbctp-resplit (B23, parent B17) · facts: F59 F61 F62 F63 F64
- in: Ext `AIkill_Dynamic/<key>/d/{input.mrd, recon.mat, tyger.log}`; root `raw.py`/`gtypes.py` imported unchanged
  (CPU spectral fit); `helpers/snr_calc.py` masks (corner box, lung_mask, EXCL_DILATE) for comparability with note 04
- out: `outputs/resplit_2026-09-24/{fits,rows,logs,fig}/`, `resplit_summary.csv`; Ext `d/recon_resplit.mat`
  (aRBC, aTP float32 (nbins,z,y,x) + dphi_old/new, ph_new, R_new, gate flags, corr/negTP per bin) — never touches recon.mat
- method: stages fit (parallel subprocesses) → split (invert old basis from saved aRBC/aTP, masked complex sum,
  0.001-rad sweep, min |R−target| with both sums > 0; gate loose/strict) → summary → fig / cohort
- gotchas: killpts term is −13° not −1° (F61); a gated session's stale .mat from an earlier run must be deleted by
  hand (042DR was); `_merged` and `_REJECTED` folders are included by discovery — filter by key when tabulating
- status: DONE 2026-09-24 — 79 re-split, 10 gated; 045VS figure shown; Hooman's eye verdict + patch go pending

### C37 · helpers/specfit_template.py
- why: Hooman: "try your fix for RBC-low and membrane infiltration, see if it stays stable on the good data, find a threshold" — the M3
  two-membrane decomposition flips on RBC-weak blocks (F69); the split needs only one complex coefficient per line
- origin: H (order) / C (design + code) · branch: rbctp-resplit (B23) · facts: F65 F66 F68 F69 (+ result rows to come)
- in: XeCS FID cache (Ext `2026_XeCS_Recon/calspec/cache/<tag>.npz`), `notes/calspec_package_2026-09-16/{resplit_inputs_2026-09-24.csv, specfit.py}`
- out: `outputs/resplit_2026-09-24/template_study{,_adaptive}.csv`, `template_params.json`, `fig/template_study{,_adaptive}.png`, log
- method: cohort line shapes per protocol (median of stable M3 rows: RBC/mem1/mem2 shift + FWHM, a2/a1, −105° relative phase) →
  basis T_gas, T_rbc, T_mem(lumped) at the block's located gas frequency; complex lstsq for 3 coefficients on a grid over
  RBC/mem shifts (±1.5 ppm; adaptive: ±2 ppm × width scales 0.7–1.4 × a2/a1 0.8–2.5); dphi_k0 and lumped/scalar ratio
  from the coefficients at t_k0 = 10 µs; stability = drop-2 vs drop-3, odd/even reps, first/second half; same tests on specfit
- status: fixed shapes: stability 2–4× better than specfit, agreement −1.6° median MAD 4.3° on stable rows, but >10° off on
  19/74 with high residual (frozen shapes misfit); adaptive grid running 2026-09-25 02:40

### C38 · helpers/te90_table.py
- why: Hooman 2026-09-25: "I want a comparison of TE90 in all subjects" — no new fitting, arithmetic on the cal-block fits
- origin: H (order) / C (code) · branch: aikill-atlas (B12) · facts: F59 F61 F63 F72
- in: `outputs/resplit_2026-09-24/resplit_summary_xecs.csv` (φ_k0 = dphi_new, Δf = df_xecs, carrier, gates);
  `notes/calspec_package_2026-09-16/resplit_inputs_2026-09-24.csv` (t_k0_us, TE_us — the time φ_k0 refers to: 610 µs v2 / 630 µs v3)
- out: `outputs/te90_2026-09-25/te90_table.csv` (88 rows: TE, t_k0, Δf, φ_k0 ± rep sd, φ at TE, TE90 ± sd, ΔTE, F_lump, atlas_set flag),
  `te90_strip.png` (TE90 + φ_k0 per session, carrier colour, gated = x, 2023|2024 divider)
- method: φ(t) = φ_k0 + 360·Δf·(t − t_k0) (RBC − TP advancing with t, F59); TE90 = t_k0 + (90° − φ_k0)/(360·Δf);
  TE90_sd = rep-to-rep sd of φ_k0 / (360·Δf). `atlas_set` = carrier ON_RBC ∧ status resplit (61) — the C39 input list
- gotcha: TE_ms in the summary CSV is 0.62 for all rows (MRD header) while the calspec twix says 600 µs for v2 — the table
  uses the calspec numbers because φ_k0 was evaluated at that t_k0; TE90 is a REPORT number, the split uses the measured φ_k0

### C39 · helpers/atlas/rbc_tm_atlas.py
- why: the B12 atlas dp pages were SUSPECT (F47, Steve's split); Hooman 2026-09-25: regenerate from the re-split maps,
  RBC-carrier sessions only ("we hit the RBC and in a few the middle — do not mix them"), ignore sessions without inputs
- origin: H (order + scope) / C (code) · branch: aikill-atlas (B12) · facts: F47 F66 F70 F72
- in: `outputs/te90_2026-09-25/te90_table.csv` atlas_set == yes (61 keys); Ext `AIkill_Dynamic/<key>/d/recon.mat` (gas_phase) +
  `d/recon_resplit_xecs.mat` (aRBC, aTP, dphi_new_deg, F_lump, insp_bin, carrier, gate); imports layout helpers from C26
- out: `outputs/aikill_atlas/rbc_tm_atlas.pdf` (per ID: panel pages gp / rbc / tm, then per date full-slice pages cor/sag/ax × 3 kinds;
  last page = appendix listing every Ext session NOT included and why), `videos_rbctm/<ID>_{gp,rbc,tm}.mp4` (16 bins × 5; 934 MB → symlink to Ext `Work/Codes/2026_ASAP_Recon/aikill_atlas/videos_rbctm/`, same pattern as C26 videos/),
  `rbc_tm_atlas_sessions.csv` (bin shown, φ_k0, F_lump, ratio, vmax per kind), `rbc_tm_atlas.log`
- method: C26 layout unchanged (gas-picked bin = max gas total signal, gas-mask extent, shared crop window per ID, 10 slices at
  apex→base fractions); RBC/TM = max(aRBC,0) / max(aTP,0), window 0 … 99.5th percentile inside the gas mask at the shown bin;
  titles carry φ_k0, F_lump, bin shown + insp_bin. LUMPED scale (TM = |mem1+mem2|)
- status: DONE 2026-09-25 02:04 — 61 sessions / 51 IDs, 703 pages, 203 MB (under the 250 MB laptop rule), 153 mp4; appendix lists 43 excluded. Eye-flag: 2024-12-05_043AS RBC is noise-only (stable=no, F_lump=1.00 → M2 fallback?) — Hooman's verdict pending

### C41 · helpers/atlas/perbin_phaseref.py
- why: Hooman "let's try your fix": per-bin phase reference from the bin's own complex gas image (consortium-style) to
  remove the F74 static RBC imprint
- origin: H (go) / C (design + code) · branch: perbin-phaseref (B24, parent B12) · facts: F74 F77
- in: fork f5ac7c5 output.mrd (gas_phase_image, gas_phase_complex, dissolved aRBC+i·aTP + rbc_tp_* meta) on Ext
  `Work/Codes/2026_ASAP_Recon/tyger_gascplx_2026-09-25/<key>/d/` (045VS, 043AS 02-21, 000LL 2023-01-31; spec
  `pipeline/recon_codespec_f5ac7c5.yml`, submit `pipeline/run_gascplx_2026-09-25.sh` from stored d/input.mrd)
- out: `outputs/perbin_phaseref_2026-09-25/<key>_s{3,6}.json`, `fig/<key>_s3.png`
- method: un-rotate the container split (z = z_rot·e^{−i ph_b}), ψ_b = phase of lung-masked Gaussian-smoothed G_b
  (σ 3 / 6 vox; variants gas, gas 2σ, self-reference 3σ as the bound), z·e^{−iψ_b}, re-sweep to the target, re-split;
  F74 static test + out-of-wedge before/after
- result: ψ_b in the lung median 0.6–1.7°, p90 2–4° → aRBC unchanged (corr 0.99–1.00), imprint index unchanged; even
  the self-reference bound leaves it → F74 mechanism retracted; F77

### C42 · helpers/atlas/rbc_repro_tests.py
- why: after C41 failed, decide whether the aRBC fine structure is noise, a fixed artifact, or signal — three
  reproducibility tests (mirror bins, different subjects at fixed voxels, same subject across dates)
- origin: C (design) after Hooman's "I saw it in RBC separately" · branch: perbin-phaseref (B24) · facts: F75 F76
- in: Ext `AIkill_Dynamic/<key>/d/{recon.mat, recon_resplit_xecs.mat}`, atlas set (61); split_trust.csv grades for the
  cross-subject subset (A/B)
- out: `outputs/split_trust_2026-09-25/{mirror_bin_repro, cross_subject_fixedvoxel, cross_date_same_subject}.csv`,
  `fig/2024-09-10_045VS_mirrorbins.png`; also `static_test_maps.csv`, `static_test_basal.csv` (inline, same session)
- result: mirror hp corr aRBC 0.55 (signal, F75); cross-subject 0.06 (no fixed pattern); cross-date 0.31 (subject-
  specific, F76). First run was inline (ledger s11); file saved verbatim afterwards, not re-run
- addendum 06:10 (inline, not in the file): fine-scale (< 1.5 vox) corr between bins 0 and 7 inside/outside the calcb
  support → `fine_static_support.csv` (F79): aRBC 0.84 / aTM 0.77 / |dis| 0.47 / gas 0.33 / background 0

### C43 · helpers/atlas/bphase_fix.py
- why: Hooman on 050EH_rbc.mp4: "there is non-moving structure that is the same in all frames, what is that?" — the
  frame-invariant voxel-scale speckle inside a sharp lung boundary (F79) = calcb b's fine phase; remove it offline
- origin: H (observation) / C (diagnosis + code) · branch: perbin-phaseref (B24) · facts: F78 F79
- in: Ext `AIkill_Dynamic/<key>/d/{recon.mat, recon_resplit_xecs.mat}` (common frame z_rot = aRBC·e^{iφ}+aTP)
- out: `outputs/bphase_fix_2026-09-25/<key>.json`, `fig/<key>.png` (before/after bins 0/7/12, δ map, MIN over bins,
  fine components); `--write` → Ext `d/recon_resplit_xecs_bfix.mat` (aRBC, aTP, delta_deg)
- method: δ(x) = angle(mean_b z_rot) − smooth(angle, σ 2.5 vox, lung-weighted); z' = z_rot·e^{−iδ}; re-sweep to the
  spectroscopic target; re-split. Metric = corr of the < 1.5-vox component between bins 0 and 7 (static texture)
- caveat: proxy — removes ALL static fine phase, real or not (unmeasurable at per-voxel RBC SNR ~2 anyway); the
  proper fix is in calcb (low-pass b's phase inside the support = 2steve/06 §4 hybrid_b) — pending b export c8366c3
- result: 050EH static corr 0.82 → −0.10, oow 16.7 → 10.1 %; 045VS 0.93 → 0.32, 5.3 → 2.3 %; |dis| unchanged

### C44 · helpers/atlas/bphase_confirm.py
- why: F79 was inferred from the maps alone; confirm with b itself and compare the proxy fix with the true one
- origin: C · branch: perbin-phaseref (B24) · facts: F79
- in: fork c8366c3 output.mrd (`calcb_b` item + images + rbc_tp meta) on Ext `Work/Codes/2026_ASAP_Recon/tyger_bexport_2026-09-25/<key>/d/`
  (050EH, 045VS; spec `pipeline/recon_codespec_c8366c3.yml`, `pipeline/run_bexport_2026-09-25.sh`)
- out: `outputs/bphase_fix_2026-09-25/<key>_confirm.json`, `fig/<key>_confirm.png`
- method: fine phase of b (σ 2.5 lung-weighted smoothing removed) vs the proxy δ (bin-mean z) and vs the observed fine
  aRBC (prediction sin δ_b·|z|); three splits: container, proxy fix, true-b fix (z·e^{−iδ_b}); static-texture metric + oow
- result: prediction corr 0.91 / 0.95; proxy vs b fine phase corr 0.68 / 0.92 inside the support, 0.18 / 0.50 outside
  (polynomial region has no fine b phase — the proxy removes other static fine phase there); true-b fix reaches the
  |dis| floor; proxy overshoots slightly (removes real static fine phase too). Maps proxy vs true-b corr 0.95–0.98

### C45 · helpers/atlas/bphase_sigma.py
- why: Hooman: "a low-pass is still a fixed structure on every bin — what is the REAL correct thing?" → measure the
  static component of the split maps per spatial band as a function of the smoothing width, against the |dis| floor
- origin: H (question) / C (code) · branch: perbin-phaseref (B24) · facts: F79 F80
- in: fork c8366c3 output.mrd with `calcb_b` (Ext tyger_bexport_2026-09-25/{050EH,045VS}/d/)
- out: `outputs/bphase_fix_2026-09-25/<key>_sigma.csv`, `fig/<key>_sigma.png`
- method: for σ ∈ {0, 1.5, 2.5, 4, 6, 10, quadratic-everywhere}: remove (phase_b − lowpass_σ) from every bin, re-sweep,
  re-split; corr(bin0, bin7) of aRBC / aTM in bands < 1.5, 1.5–6, > 6 vox vs the same for |dis| (cannot carry phase)
- result: σ 2.5 puts aRBC at/below the |dis| floor in every band on both sessions; σ ≥ 4 no gain; quadratic everywhere
  worse on negatives → fork patch default σ = 2.5

### C46 · helpers/atlas/bsmooth_compare.py
- why: prove the fork fix (19b7365) does in the container what the offline true-b fix did, and that the gas image is
  not harmed
- origin: C · branch: perbin-phaseref (B24) · facts: F79 F80
- in: Ext `tyger_bexport_2026-09-25/<key>/d/output.mrd` (raw b, c8366c3) + `tyger_bsmooth_2026-09-25/<key>/d/output.mrd`
  (19b7365; spec `pipeline/recon_codespec_19b7365.yml`, `pipeline/run_bsmooth_2026-09-25.sh`); 050EH, 045VS (000LL
  reconstructed too, no raw-b twin to compare)
- out: `outputs/bphase_fix_2026-09-25/<key>_bsmooth.json`, `fig/<key>_bsmooth_frames.png` (aRBC bins 0/4/7/11 raw vs
  smoothed, coronal + sagittal)
- result: F80 numbers. Frames: sharp lung boundary and fine speckle gone; 050EH RBC = smooth blobs (its true SNR)

### C47 · helpers/atlas/bsmooth_videos.py
- why: Hooman: "can I see the videos for all three" — the atlas-format bin videos with raw-b and smoothed-b rows
- origin: H · branch: B24 · facts: F80
- in: Ext AIkill_Dynamic (raw maps) + tyger_bsmooth_2026-09-25 output.mrd (19b7365); reuses C39 Session/IDGroup/render_video
- out: `outputs/bphase_fix_2026-09-25/videos/<ID>_{gp,rbc,tm}.mp4` (rows raw b / smooth b); inline follow-up (ledger 08:10,
  not in the file): `videos_sigma/<ID>_*.mp4` rows raw b / σ 2.5 / σ 10 / no b from the raw-b export (045VS, 050EH, 000LL)
- note: Hooman first read the two-row video as "still there" — the rows were not obvious; label reads `cor raw b rbc` etc.

### C48 · helpers/atlas/bphase_sigma_cohort.py
- why: Hooman: "do 1, 2.5, 5 and 10 and decide" — on the cohort, not on two sessions
- origin: H (order) / C (code) · branch: B24 · facts: F80 (+ result row to come)
- in: Ext `Work/Codes/2026_ASAP_Recon/tyger_bexport_2026-09-25/<key>/d/output.mrd` (raw b export, c8366c3) for the 61
  atlas sessions — produced by `pipeline/bexport_batch.sh` (xargs −P 3 over `pipeline/bexport_todo.txt`, restart-safe
  skip on existing ≥ 300 MB output; detached with nohup setsid because Claude background tasks die at 10 min)
- out: `outputs/bphase_fix_2026-09-25/cohort_sigma.csv`, `cohort_sigma.png`
- method: per session/σ the C44 offline split; static bands vs |dis| floor, out-of-wedge, lung RBC / bg sd, corr vs raw
- status: DONE 17:20 — 57 sessions in the first pass (F81 table); 4 stragglers (039CP 03-11, 036RL 09-27 truncated by
  round cuts; 020JS, 051VM left unsplit by the container gate → fallback to the offline angle/target + one raw-b split)
- batch notes: Claude background tasks die at 10 min and `nohup setsid` dies with the tool call — the batch ran as a
  Monitor command in 30-min rounds (restart-safe skip ≥ 300 MB; a cut download > 300 MB slipped through once →
  EOFError caught by the sweep); hotspot (172.20.10.x) kills 400 MB downloads; Ext dropped twice (bad cable, Hooman)
  → mount guard added to run_bexport; 2 in flight after the drop

### C49 · pipeline/run_prod_b44.sh (+ prod_b44_todo.txt)
- why: Hooman 17:45: run production with σ 4.4 for every session with clean spectroscopy/fit (61) before discussing the rest
- origin: H (order, σ) / C (script) · branch: prod-b44 (B25) · facts: F79 F80 F81
- in: Ext `AIkill_Dynamic/<key>/d/input.mrd` (stored DIAPHRAGM inputs); spec `recon_codespec_04f445b.yml`
- out: Ext `AIkill_Dynamic_b44/<key>/d/{output.mrd, tyger.log, codespec.yml}` then `recon.mat` via C50
- method: as run_bexport (mount guard; skip when output exists, > 300 MB AND fully readable by the mrd reader — the
  earlier size-only skip let a truncated file through); run as Monitor rounds, 2 in flight (bad cable)
- status: launched after CI 36192122624

### C50 · pipeline/mrd_to_mat.py
- why: post_process.py makes gifs/figures (~1 min/session); the atlas and trust tools need only recon.mat
- origin: C · branch: prod-b44 (B25) · facts: —
- in: `<session>/d/output.mrd` (fork ≥ c8366c3 items) · out: `<session>/d/recon.mat` (gas_phase, gas_phase_magnitude,
  dissolved_phase_real/imag = aRBC/aTP when rbc_tp_separated = 1, dissolved_phase_magnitude, calcb_b, rbc_tp_* strings, nav_*)
- note: the OLD recon.mat convention (post_process) has the same names; readers must check `rbc_tp_separated`

### C40 · helpers/atlas/split_trust.py
- why: Hooman on the RBC/TM atlas: "I still see the fake lung structures you thought came from the static b maps — check
  that; and give me some marker of how much to trust the separation"
- origin: H (question) / C (tests + code) · branch: aikill-atlas (B12) · facts: F47 F66 F73 F74 (+2steve/03, /06)
- in: Ext `AIkill_Dynamic/<key>/d/{recon.mat (gas_phase), recon_resplit_xecs.mat}`; `outputs/te90_2026-09-25/te90_table.csv`
  (61 keys + fit numbers); masks from `helpers/snr_calc.py` (same as the re-split)
- out: `outputs/split_trust_2026-09-25/split_trust.csv` (61 rows), `cohort_trust.png`, `fig/<key>.png` (EI+EE coronal:
  gas, |z|/σ, phase-in-wedge, out-of-wedge + calcb support contour, resolvable), `fig/<key>_static.png` (EI vs EE gas /
  |dissolved| / RBC fraction, coronal + sagittal), `run.log`; grades merged into `aikill_atlas/rbc_tm_atlas_sessions.csv`
- method: z = aRBC·e^{iφ_k0} + aTP (stored basis); PHYSICAL wedge = aRBC ≥ 0 ∧ aTP ≥ 0; out-of-wedge = error. Noise per
  voxel σ(x) = σ_corner·A(x)/⟨A⟩_corner with A the un-divided gridding envelope (2steve/03; fitted radially per session,
  L² ≈ 8–10k vs kernel theory 5066 — both carried, conclusions hold under both); noise-only prediction of out-of-wedge
  from each voxel's phase-SNR; regions = inside/outside calcb support (|F_avg| > 10·mean|noise| proxy) + rim band.
  Static test: global z-shift EE→EI from gas x-corr; corr(EI,EE) of RBC fraction at fixed vs shifted voxels, |dissolved|
  as the anatomy yardstick; imprint index = (fixed−shifted)_rbcfrac − (fixed−shifted)_|dis|. Resolvable voxel: σφ ≤ φ_k0/4.
  Grade = worst of {gain ≤1.35/1.7, rep-sd ≤4/8°, resolvable ≥60/35 %, excess ≤2/5 %}; two C's → D.
- gotchas: corner σ alone underestimates central noise 2–4× (envelope) → the first pass called 99 % of voxels resolvable;
  the z-shift is one rigid number per session (apex moves less than the base) — the yardstick absorbs most of that
- status: DONE 03:20; 0 failures; verdicts in F73/F74

## Outputs index

| Output dir (workspace/outputs/) | From | Facts | Status |
|---|---|---|---|
| te90_2026-09-25/ | C38 | F72 | VALID (te90_table.csv 88 rows, te90_strip.png) |
| aikill_atlas_b44/ | C39 (--root) | F82 | VALID — PRODUCTION atlas (σ 4.4): rbc_tm_atlas.pdf, rbc_tm_atlas_sessions.csv (+ grades), videos_rbctm → Ext Work/Codes/2026_ASAP_Recon/aikill_atlas_b44/ |
| split_trust_b44/ | C40 (--root) | F82 | VALID — PRODUCTION trust table (split_trust.csv, fine_static_b44.csv, cohort_trust.png, fig/) |
| split_trust_2026-09-25/ | C40 C42 | F73 F74(R) F75 F76 | VALID (split_trust.csv 61, cohort_trust.png, fig/<key>{,_static}.png ×61; static_test_{maps,basal}.csv, mirror_bin_repro.csv, cross_subject_fixedvoxel.csv, cross_date_same_subject.csv, fig/…_mirrorbins.png). imprint_index column = weak statistic, F74 retracted |
| perbin_phaseref_2026-09-25/ | C41 | F77 | VALID (3 sessions × σ 3/6 json + fig) |
| bphase_fix_2026-09-25/ | C43 C44 C45 C46 | F79 F80 | VALID (050EH, 045VS: <key>.json, _confirm, _sigma.csv, _bsmooth.json + figs) |
| ct_overlay/ | C33 C34 C35 | F58 | VALID — RT 3 + LTX 20 + EBV 9 keys (cohort_fit_table.csv, cohort_summary.pdf, bins_video.mp4 per key; 038RL poor fit) (per subj: summary_3plane, axial/coronal montages, register.json, orientation_scores.csv, phase_scores.csv; gifs + per-slice PNGs + NIfTI on Ext `Codes/2026_ASAP_Recon/ct_overlay/`) |
| snr_2026-09-13/ | C28, _rbctp_fig.py | F47 F48 F52 | VALID (snr_table.csv, rbc_tp_solve_check.csv, 04_rbctp_split_evidence.png, test CSVs) |
| f59_te_term_2026-09-24/ | _f59_check.py (scratch) | F59 | VALID (f59.json 208/218-ppm sessions, f59_oddsessions.json bad fits, 05_rbctp_split_te_term.png, 05b_…_218ppm_badfits.png → 2steve/05) |
| calcb_imprint_2026-09-24/ | _calcb_8ch_weight.py, inline recon.mat pass (scratch) | F47 | VALID (neg_imprint.json + 06_neg_imprint_*.png per-bin negative fraction; W_2024-01-17_042DR.{json,npy,log} + 06_W_map_*.png Σ\|b\| static gain, 8 ch → 2steve/06) |
| structure_rank/ | C30 | F57 | VALID (panels sorted by G_struct; eye-pick tool, not a verdict) |
| twix_audit_2026-09-16/ | inline script (ledger s7) | F55 | VALID (twix_audit.csv, 104 rows) |
| merged_dyn_2026-09-16/ | C31 C32 + snr_calc | F56 | VALID (snr_merged_vs_orig.csv, panel_orig_vs_merged_gas_d.png; 025VP merge hurts) |
| archive/ | ? (pre-canon) | ? | untagged — verify before trust (reconciled 2026-09-17) |
| roundtrip_audit_2026-09-13/ | _roundtrip_audit.py | F51 | VALID (roundtrip_audit.csv, 136 folders) |
| steve_vs_fork_diff/ | git diff main/dev/diaphragm-recon | F49 | VALID (4 .diff files) |
| zorder_check_2026-09-12/ | _zorder_check.py | F49 | VALID (030DN edge-finder figure) |
| leftover_recon_2026-09-12/ | copies of AIkill_Dynamic montages | F47 | VALID gas; dp montages magnitude-of-split caveat |
| aikill_atlas/ | C26 C39 | F47 F72 | dissolved_atlas.pdf + videos/: gp VALID, ⚠ dp pages SUSPECT (SUSPECT_F47.md) — superseded by rbc_tm_atlas.pdf + videos_rbctm/ + rbc_tm_atlas_sessions.csv (C39, 61 ON_RBC sessions, 2026-09-25) |
| rt_prepost/ | C27 | — | VALID |
| diaphragm_binning/ | C24 C25 | F37 F38 | VALID |
| 016PG/ · 023LL/ · 025JC/ · 25JC/ · 025JC_sweep_lt/ · piston/ | ? (June backlog, pre-canon) | ? | untagged — verify before trust (reconciled 2026-09-13) |
| resplit_2026-09-24/ | C36 C37 | F59 F61–F71 | VALID — fits/ (89 Steve refits), rows{,_xecs}/, resplit_summary{,_xecs}.csv (79 / 68 re-split), fit_xcheck.csv (Steve vs XeCS, 70), template_study{,_adaptive}.csv + template_params.json (C37), logs; fig/: cohort + per-session before/after (steve, _xecs), tyger_old_vs_new_6651591.png, template_study{,_adaptive}.png; big .mat on Ext `AIkill_Dynamic/<key>/d/recon_resplit{,_xecs}.mat` (79 / 68 × 128 MB) and Tyger validation runs Ext `Work/Codes/2026_ASAP_Recon/tyger_specfit_2026-09-24/{045VS,041WF}/d/` |
