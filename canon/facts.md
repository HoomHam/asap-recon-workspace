# ASAP Recon — Facts (belief table)

Flip Status, never delete rows. VALID / RETRACTED / SUSPECT.
Retro-filled 2026-07-12 from root handoff-report.md + workspace handoffs (Jun 2026).

| ID | Claim | Status | Since | Note |
|----|-------|--------|-------|------|
| F1 | Pipeline DOWN: newer numba in fork image `d136eb1` rejects `recon.py` CUDA kernels; every guarded recon blocked until numpy/scipy/numba pinned in requirements.txt | RETRACTED→F42 | 2026-07-10 | root handoff |
| F2 | Always use immutable `:<sha>` container tag in codespec, never mutable `diaphragm-recon` tag — Tyger nodes cache the mutable tag | VALID | 2026-06-19 | |
| F3 | requirements.txt must pin mrd-fork to `b6b6d18` or run dies with `RuntimeError: Invalid schema` | VALID | 2026-06-19 | |
| F4 | Guard predicate is NOT "multi-coil" — it is "no spectral fit succeeded"; `025VP` is nch=1 with numspec=0 and guard still fires | VALID | 2026-07-10 | |
| F5 | `len(fRBC) > 0` NOT sufficient to trust dissolved-phase output — guard never checks fit params are physical (Bug B) | VALID | 2026-07-10 | |
| F6 | Do NOT trust nch=1 "guard silent" CPU verdict as proof a subject reconstructs — fit runs and diverges differently per scipy version (the `007IT` mistake) | VALID | 2026-07-10 | |
| F7 | Container's older scipy is the honest one; loosening maxfev would make `007IT` "run" and silently emit a WRONG dissolved-phase image | VALID | 2026-07-10 | |
| F8 | `probe_nch.py`/`check_frbc.py` are faithful oracles for nch=8 (fit never runs) but NOT for nch=1 (fit runs, can diverge) | VALID | 2026-07-10 | |
| F9 | Multi-coil spectral-fit results discarded at raw.py:362, but any one of 120 discarded curve_fit failures still aborts whole recon (Bug A); one-line fix `if(numspec>0 and self.nch==1)` | VALID | 2026-07-10 | |
| F10 | "55/62 done" subject count | RETRACTED | 2026-07-10 | real number = 83 green, matches 83 dirs on HoomHamExt |
| F11 | "`005DS` (row 46) is a permanent fail" | RETRACTED | 2026-07-10 | plain guard case (nch=8, fRBC=0), reconstructs gas-only |
| F12 | The 26 NO DYNAMIC SCAN subjects were never acquired with spiral-dyn sequence, unrecoverable; remove from denominator (incl. row 3 `001SW`, highest-SNR red herring) | VALID | 2026-07-10 | |
| F13 | "missing `_dp.npy` → gas-only, warned not fatal" (pipeline/AGENTS.md claim) | RETRACTED | 2026-07-10 | FALSE: raw.py multiplies dp_array (None) unconditionally and crashes (`042DR`) |
| F14 | exfat/fskit volumes wedge to EPERM at volume root (both HoomHam & HooMain); NOT TCC — reproduces with FDA; fix `diskutil unmount force … && diskutil mount diskNsN` | VALID | 2026-07-10 | |
| F15 | Pipeline needs `mrd` AND `openpyxl`; only conda base python has `mrd` (arm64 .venv does not) — launch batch with conda base python | VALID | 2026-07-10 | |
| F16 | AppleDouble `._*.dat` shadow files on exfat break convert_siemens_to_mrd (picks 2nd-largest .dat as breath-hold ref); any new raw-folder reader must filter `._*` | VALID | 2026-07-10 | |
| F17 | recon.py has Kaiser-Bessel kernel (`bessi0`) but COMMENTED OUT — active gridding kernel is Gaussian | VALID | — | |
| F18 | raw.py infers npts from |k|² periodicity, not header — fragile if trajectory changes | VALID | — | |
| F19 | Steve's output background looks negative because complex/zero-mean; take np.abs() for MATLAB-like (Rician) visuals | VALID | 2026-06-08 | |
| F20 | `asap/asap.c` (Steve) and `workspace/codes/kasap.c` (Kento) are different — never assume equivalence | VALID | — | |
| F21 | On spiral operator, gradient solvers (FISTA/PDHG) stall ~1000× without DCF preconditioning while CG converges; unconverged output still looks like an image | VALID | 2026-06-12 | |
| F22 | Faraz recon ×1.205 magnified; root cause = `resizing` k-rescale in his analysis script (not his KB+DCF engine); correct fix = discard |k|>Nyquist or grid larger, NEVER rescale k | VALID | 2026-06-11 | |
| F23 | Tikhonov λI is a no-op on fully-sampled data (AᴴA well-conditioned; penalizes amplitude not roughness) — gap to Steve's SNR closes only via smoothing/CS regularizer | VALID | 2026-06-11 | |
| F24 | Trust Hooman's eye over scalar metrics: interior/low-freq-CV SNR misled twice in one day; SNR inflates when priors zero the background | VALID | 2026-06-12 | |
| F25 | "ours-wavelet wins on lfCV" | RETRACTED | 2026-06-24 | was blur erasing resolution — do not resurrect |
| F26 | "BART has a moiré artifact" | RETRACTED | 2026-06-24 | it is the real ACR insert; BART resolves it |
| F27 | "DCF removes the grid" | RETRACTED | 2026-06-24 | disproven |
| F28 | "(256,64,256) DCF-size bug in real Lustig workflow" (2026-05 docs) | RETRACTED | 2026-06-15 | notebook uses im_size=(100,100,100) |
| F29 | Lustig `_Wavelet` variant has real bug: FWT2_PO (2D) applied to 3D volume forcing 128³; `_Wavelet_Adop` misnamed (weights 0 → pure least-squares); DCF only in init so final ≈ DCF-gridded init (stalled NCG) | VALID | 2026-06-15 | |
| F30 | "DIAPHRAGM binning silently falls back to SIGNAL (tyger_recon.py lacks navigator path)" | RETRACTED | 2026-06-19 | FIXED: ported main.py calcLVcb + flipped GPDYN z (`[::-1]`) so edge-finder tracks diaphragm not apex |
| F31 | Local edits to main.py/raw.py/tyger_recon.py/convert_siemens_to_mrd.py live on fork branch only; root pull can wipe uncommitted ones | VALID | 2026-06-12 | |
| F32 | MATLAB resolves by filename — renaming any Faraz .m breaks callers; `grid_lookup_20220418.mat` = stale neighbor cache, delete if traj/FOV changes; `spiral_human_20240227.m` has hardcoded Windows paths | VALID | 2026-06-08 | |
| F33 | mapvbvd end-of-file UserWarning (byte-offset overrun) during .dat read is benign | VALID | 2026-06-18 | |
| F34 | Do NOT submit a Tyger GPU job to test a hypothesis testable on CPU (one run burned confirming what diag_fit.py explained free) | VALID | 2026-07-10 | |
| F35 | UNVERIFIED: whether CS forward/adjoint operator applies correct-resolution B0/b-matrix (vs low-res one from navigator loop) | SUSPECT | 2026-06-19 | open caveat, never checked; CS side now lives in XeCS |
| F36 | Root repo: `origin` = MEDCAP (Kento) PULL-ONLY; pushes go to `hooman` fork (HoomHam/asap_recon) on branch `diaphragm-recon`; workspace/ = own repo, hidden via .git/info/exclude | VALID | 2026-07-12 | confirmed this session |
| F37 | ASAP tyger navigator undersampled-image = 26 interleaves (= 1 complete Thomson set of v3 traj; `ilvperusimg = nuniqueilvs/nusimg = 832/32`), block-aligned from ilv 0; nusimg=32 hardcoded (convert:183, tyger reader default). killpts trims SAMPLE axis only → interleave block phase preserved. MEASURED on v3_dyn/v3_fov250/v3_dyn_025JC (N=424320, npts=510, nuniqueilvs=832). It was NEVER 20 — the XeCS "20 bad default" did not apply here | VALID | 2026-07-17 | measured, not assumed |
| F38 | DIAPHRAGM `bin()` (raw.py:153) is equal-COUNT amplitude rank (percentile), not equal-mL: respects inflation ORDERING (short strokes never forced into deep bins) but bin width in mL varies (narrow at turns, wide mid-stroke) → residual iso-inflation blur. True same-size binning needs fixed-mL-width position bins (accept unequal/sparse counts). Design change, not current code | VALID | 2026-07-17 | Hooman's iso-inflation point |
| F39 | Tyger service-principal cert (`tyger-sp.pem`) expires YEARLY; Kento issues the new one. Current pem valid to 2028-08-04. Failure signature: `AADSTS700027 ... certificate ... expired`. Check: `openssl x509 -in tyger-sp.pem -noout -enddate` | VALID | 2026-09-03 | renewed this session |
| F40 | Tyger output download is ~1 Mbps from this laptop; `run exec --logs` streaming and default `buffer read` (dop 32) both fail on 192 MB outputs. Reliable = `tyger buffer read <id> -o file -p 4` (~10 min). Output buffer id is on the "buffer read" lines of tyger.log; TTL 1 hr after job end | RETRACTED→F43 | 2026-09-03 | measured run 558; speed varies |
| F41 | Image `db80f16` still runs on Tyger (run 558, 2026-09-03): single-coil subjects reconstructable now; F1 (numba) only blocks nch=8 subjects needing the dissolved guard | RETRACTED→F42 | 2026-09-03 | |
| F42 | Fork image `472fbc9` pins numpy 2.4.6 / scipy 1.17.1 / numba 0.65.1 / llvmlite 0.47.0 (= db80f16's working set); numba 0.66 `TypeError: Signature mismatch` in calcb gone — 8-ch sessions (005DS, 017AK, 023DB, 004LR, 029CK, 02BB, 032WS) reconstruct on Tyger. Current image `40d23a4` (+ unsplit-dissolved tag) | VALID | 2026-09-13 | runs 559–6xx |
| F43 | Tyger download speed is NOT fixed: ~1 Mbps on 09-03, 192–209 MB in 7–12 s on 09-12/13. dyn_recon.submit() now = buffer create → buffer write -i → run create -b → poll run show → run logs → buffer read -o -p 4, all retried (no stdout streaming) | VALID | 2026-09-13 | |
| F44 | `tyger buffer create --ttl <D.HH:MM:SS>` returns HTTP 500 on the spinhance Tyger server; plain `buffer create` works — never pass --ttl | VALID | 2026-09-13 | reproduced |
| F45 | `tyger_deploy/plot_recon.py` ends with `plt.show()`; on the macosx matplotlib backend an unattended batch BLOCKS until the window is closed (40-min stall, 000LL 2024-02-02). dyn_recon runs plot/post with `MPLBACKEND=Agg` | VALID | 2026-09-13 | |
| F46 | 2024-01 8-ch sessions 042DR (01-17), 001JM (01-18), 007RA (01-22), 008CR (01-31) are gas-only ACQUISITIONS (raw.py FFT pattern `gas-only`). No successful cohort run ever dropped dissolved; ledger 2026-08-27 attribution "RBC/TP params absent" was wrong | VALID | 2026-09-13 | tyger.log pattern lines |
| F47 | Steve's RBC/TP split (results.py 298–327) is unreliable: stop test `(R−target)(lastR−target)<0` fires at the ΣaTP=0 discontinuity and R uses whole-volume sums → logged R off target in 60/86 split sessions (19 valid, 7 partial). One-point Dixon is ill-conditioned: corr(aRBC,aTP)=−cosΔφ, per-part noise gain 1/\|sinΔφ\| (045VS Δφ≈6°, ×9). `dissolved_phase_real/imag` AND `dissolved_phase_magnitude`=\|aRBC+i·aTP\| are not trustworthy for split sessions | VALID | 2026-09-13 | outputs/snr_2026-09-13/rbc_tp_solve_check.csv; 2steve/04 |
| F48 | Total dissolved SNR done right = covariance-whitened magnitude: q=vᵀC⁻¹v=\|z\|²/σ² (C = background covariance of stored real/imag), E[q]=2 in noise (measured 2.00–2.09, split and unsplit), SNR = mean_lung √max(q−2,0); phase-independent | VALID | 2026-09-13 | helpers/snr_calc.py |
| F49 | Fork DIAPHRAGM z-flip (`db80f16`) is a real deviation from Steve, not an MRD artifact: Steve mapVBVD raw == Kento MRD raw bit-identical (030DN); Steve's unflipped edge-finder tracks upper (apex-like) edge z≈36, 7.7 px excursion vs flipped diaphragm z≈67, 16.8 px | VALID | 2026-09-12 | n=1 (030DN), helpers/_zorder_check.py |
| F50 | SNR_Table_All.xlsx col E "Ave. SNR" = Hooman's hand table (static, ranking only), not a measured recon SNR; measured gas/DP SNR = cols I–N | VALID | 2026-09-13 | XeCS notes/SNR_Table_Hooman_2026-09-09.md |
| F51 | 20 roundtrip image sessions on 2022 sequences (fancy_20220421/0422/0602/0812/1103; 15 in SNR table) have NO trajectory file anywhere on Ext or laptop → not reconstructable until trajectories are generated | VALID | 2026-09-13 | roundtrip_audit.csv |
| F52 | Gas max-SNR bin ≠ end-inspiration (max lung-volume) bin in 45/94 sessions (27 before EI, 18 after) — early-inhale gas delivery; DIAPHRAGM bins are phase-folded (0/15 = end-exp, EI ≈ 6–9) | VALID | 2026-09-13 | full-precision fields |
| F53 | Our MRD = MRD v2 yardl binary via MEDCAP mrd-fork `b6b6d18` (NOT ISMRMRD v1 HDF5): whole acquisitions/trajectories/pneumotach as NdArrays tagged by fork-only `meta` keys, empty `encoding`, no Acquisition records, recon params in userParameters → stock MRD v2 (arrays have no meta; schema string check) and v1 tools cannot read it; gas/dissolved split + Xe freq/BW/FOV live in raw.py, not the file | VALID | 2026-09-14 | notes/Data_SiemensToMRD_Pipeline_2026-09-14.md; fork schema from installed mrd/protocols.py vs ismrmrd/mrd main model |
| F54 | Scanner = Siemens MAGNETOM Avanto 1.494 T, syngo MR D13; .dat is VD multi-RAID container (nMeas=1, MID/FID match filename) | VALID | 2026-09-14 | n=1 dynamic .dat (2024-01 session); cohort-wide not checked |
