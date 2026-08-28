# ASAP Recon — Cards

One card per load-bearing script. Origin: H(uman)/C(laude)/A(gent)/R(etro-inferred).
Retro-filled 2026-07-12; all origins R. Quick table regenerated from bodies by /leave.

## Quick table

| ID | Script | Branch | Status |
|----|--------|--------|--------|
| C1 | main.py (root) | stem | WORKS (local, 2 uncommitted fixes) |
| C2 | raw.py (root) | stem | BROKEN (Bug A, Bug B, 042DR crash) |
| C3 | recon.py (root) | stem | WORKS on old numba only (F1) |
| C4 | results.py (root) | stem | WORKS |
| C5 | gtypes.py (root) | stem | WORKS |
| C6 | convert_siemens_to_mrd.py (root) | tyger | WORKS (F16 caveat) |
| C7 | tyger_recon.py (root) | diaphragm | WORKS (fork image, cloud GPU only) · +nav 26-guard (F37, local) |
| C8 | pipeline/batch_recon.py | batch | WORKS, blocked by F1 |
| C9 | pipeline/dyn_recon.py | batch | WORKS, blocked by F1 |
| C10 | pipeline/asap_run.py | auto-steve | WORKS |
| C11 | pipeline/post_process.py | auto-steve | WORKS |
| C12 | pipeline/param_gui.py | auto-steve | WORKS |
| C13 | pipeline/recon_codespec.yml | tyger | CONFIG (sha bump pending, F1) |
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
| C26 | helpers/atlas/dis_atlas.py | aikill-atlas | WORKS (atlas delivered 2026-08-27) |

## Cards

### C1 · main.py (root)
- why: local tkinter GUI entry; canonical source of DIAPHRAGM navigator (`calcLVcb`, :279-330) and dataset→traj resolution (`ID_callback`, :390-442)
- origin: R · branch: stem · facts: F31
- in: Siemens .dat + traj .npy via basefolder/datatype/date/subjectID/
- out: on-screen recon (needs CUDA — fails on M4)
- status: WORKS locally, two uncommitted fixes (path `/`, hidden-file filter)

### C2 · raw.py (root)
- why: TWIX load, traj load, binning, two-Lorentzian spectral fit (lorfit :13, loop :313-380)
- origin: R · branch: stem · facts: F4 F5 F9 F13 F18
- in: .dat/MRD arrays, trajectory
- out: k-space, ilvbin, fRBC/fTP ratios
- status: BROKEN — Bug A (multi-coil fit failure aborts run), Bug B (unbounded warm-start LM escapes to non-physical params), 042DR `dp_array is None` crash

### C3 · recon.py (root)
- why: CUDA gridding kernels (cudarecon/cudarenorm) via numba.cuda
- origin: R · branch: stem · facts: F1 F17
- in: k-space samples, traj, b-matrix → out: gridded volume
- status: WORKS on older numba; image d136eb1's numba rejects kernels

### C4 · results.py (root)
- why: image recon, calcb (B0/coil-phase b-matrix), dyn_recon/dyn_usimg_recon, dissolved-phase guard
- origin: R · branch: stem · facts: F4 F5
- status: WORKS

### C5 · gtypes.py (root)
- why: global enums (bintype, imgtype) + gvar defaults; basefolder :27 Hooman-specific, intentionally uncommitted
- origin: R · branch: stem · facts: F31
- status: WORKS

### C6 · convert_siemens_to_mrd.py (root)
- why: local .dat → MRD; bakes recon params into MRD header (only channel to headless cloud job)
- origin: R · branch: tyger · facts: F16
- in: Siemens .dat, --binning → out: input.mrd
- status: WORKS (breaks on `._*.dat` shadows)

### C7 · tyger_recon.py (root)
- why: Tyger cloud-GPU Docker entrypoint; hosts ported _diaphragm_navigator + MRD nav-array export
- origin: R · branch: diaphragm · facts: F2 F30
- in: input.mrd (params :66-83) → out: output.mrd
- status: WORKS on fork image; NVIDIA GPU (cloud) only

### C8 · workspace/pipeline/batch_recon.py
- why: batch driver — iterate subjects, call dyn_recon.py --methods s,p,d, green-mark Excel row, save incrementally
- origin: R · branch: batch · facts: F1 F10 F15
- in: SNR_Table_All.xlsx, source/archive drives
- out: /Volumes/HoomHamExt/Dynamic/<date>_<id>/{s,p,d}/, updated Excel
- status: WORKS (uncommitted done_already()/pneumotach edits); blocked by F1

### C9 · workspace/pipeline/dyn_recon.py
- why: per-subject driver — resolve, _clean_datadir (symlinks around ._*), convert, submit to Tyger, post-process
- origin: R · branch: batch · facts: F1 F16
- status: WORKS; blocked by F1 for guarded subjects

### C10 · workspace/pipeline/asap_run.py
- why: single-dataset orchestrator, 6 stages resolve→param GUI→convert→submit→publish→post_process ("analyze 25JC with Steve")
- origin: R · branch: auto-steve
- out: workspace/outputs/<dataset>/<timestamp>/
- status: WORKS

### C11 · workspace/pipeline/post_process.py
- why: output.mrd → recon.mat, signal_pneumo.npz, per-bin slice videos, navigator.gif, resp_traces.png
- origin: R · branch: auto-steve
- status: WORKS (resp_traces nav-z scatter in last ~15% un-cleaned)

### C12 · workspace/pipeline/param_gui.py
- why: standalone tkinter param picker (mirrors main.py fields incl. DIAPHRAGM) → params.json
- origin: R · branch: auto-steve
- status: WORKS

### C13 · workspace/pipeline/recon_codespec.yml
- why: Tyger job spec (image tag); pinned to d136eb1
- origin: R · branch: tyger · facts: F1 F2
- status: CONFIG — needs sha bump after numba repin

### C14 · workspace/helpers/build_status_tab.py
- why: regenerate Excel "Recon Status" tab from drives (blocker + HYPERLINK per ungreened subject)
- origin: R · branch: batch · facts: F15
- status: WORKS (conda base python, PYTHONPATH=.)

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
`read_mapvbvd.py` (root, standalone loader, UNKNOWN) · `workspace/codes/kasap.c` (Kento reference, F20) · `asap/asap.c` (Steve reference, F20) · `helpers/_delete/` (byte-verified CS originals, moved to XeCS 2026-06-24) · `helpers/calib/` (duplicate .npy pair) · `pipeline/runs/`, `batch_recon.log`, `__pycache__/`.

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
- origin: H · branch: main · facts: —
- in: /Volumes/HoomHamExt/AIkill_Dynamic/*/d/recon.mat (gas_phase +
  dissolved_phase_magnitude, (16,Z,Y,X)); s/recon.mat gas fallback for the 3 sessions
  whose d skipped dissolved (2024-01-18_001JM, 01-22_007RA, 01-31_008CR)
- out: workspace/outputs/aikill_atlas/dissolved_atlas.pdf (per ID: gp panel page + dp panel
  page [orientation-major blocks cor/sag/ax, per date one 10-slice row], then per date 6
  full all-slice pages cor/sag/ax × gp/dp) + videos/<ID>_{gp,dp}.mp4 (16 bins ×5, 5 fps)
- note: orientation follows pipeline/post_process.py (axial flipud, coronal reversed,
  sagittal rot90 ccw); slice extent from GAS mask (≥30 voxels above 0.15·max, bin-mean);
  cross-date column alignment = per-session slices at shared apex→base fractions + fixed-size
  crop centered on each session's lung bbox (shift+scale, no registration); black bg, no gaps
