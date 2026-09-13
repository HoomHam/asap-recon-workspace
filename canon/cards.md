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
| C26 | helpers/atlas/dis_atlas.py | aikill-atlas | ⚠ dp pages SUSPECT (F47) — regenerate with whitened \|z\| |
| C27 | helpers/rt_prepost_fig.py | main | WORKS (RT pre/post figs 2026-08-31) |
| C28 | helpers/snr_calc.py | snr-table | WORKS (gas + whitened DP SNR, 104 sessions, F48) |
| C29 | pipeline/todo_batch.sh | leftover-recon | WORKS (25-session worklist, FORCE=1) |

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
`read_mapvbvd.py` (root, standalone loader, UNKNOWN) · `workspace/codes/kasap.c` (Kento reference, F20) · `asap/asap.c` (Steve reference, F20) · `helpers/_delete/` (byte-verified CS originals, moved to XeCS 2026-06-24) · `helpers/calib/` (duplicate .npy pair) · `pipeline/runs/`, `batch_recon.log`, `__pycache__/` · 2026-09 scratch tier: `helpers/_roundtrip_audit.py` (roundtrip image-tree audit → outputs/roundtrip_audit_2026-09-13/), `helpers/_zorder_check.py` (F49 check → outputs/zorder_check_2026-09-12/), `helpers/_rbctp_fig.py` (F47 evidence figure → outputs/snr_2026-09-13/04_rbctp_split_evidence.png) · `Codes/2026_Steve_Recon/` (plain snapshot of Steve main 3303276, diff reference only).

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

## Outputs index

| Output dir (workspace/outputs/) | From | Facts | Status |
|---|---|---|---|
| snr_2026-09-13/ | C28, _rbctp_fig.py | F47 F48 F52 | VALID (snr_table.csv, rbc_tp_solve_check.csv, 04_rbctp_split_evidence.png, test CSVs) |
| roundtrip_audit_2026-09-13/ | _roundtrip_audit.py | F51 | VALID (roundtrip_audit.csv, 136 folders) |
| steve_vs_fork_diff/ | git diff main/dev/diaphragm-recon | F49 | VALID (4 .diff files) |
| zorder_check_2026-09-12/ | _zorder_check.py | F49 | VALID (030DN edge-finder figure) |
| leftover_recon_2026-09-12/ | copies of AIkill_Dynamic montages | F47 | VALID gas; dp montages magnitude-of-split caveat |
| aikill_atlas/ | C26 | F47 | ⚠ dp pages SUSPECT (SUSPECT_F47.md) |
| rt_prepost/ | C27 | — | VALID |
| diaphragm_binning/ | C24 C25 | F37 F38 | VALID |
| 016PG/ · 023LL/ · 025JC/ · 25JC/ · 025JC_sweep_lt/ · piston/ | ? (June backlog, pre-canon) | ? | untagged — verify before trust (reconciled 2026-09-13) |
