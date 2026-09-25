# Fork patch 2026-09-24 — RBC/TP split: specfit basis angle + ratio, stop criterion, gas magnitude export

**Status: APPLIED 2026-09-25 01:48 as fork commit 6651591 (branch diaphragm-recon, pushed to `hooman`; image ghcr.io/hoomham/xe-tyger-recon:6651591…; spec pipeline/recon_codespec_6651591.yml). Validated on 045VS + 041WF (F70). v3.** `fork_patch_2026-09-24.diff` (this folder) was built from copies
of the root working tree (`diaphragm-recon` @ 40d23a4 plus Hooman's uncommitted `tyger_recon.py`
edits) and passes `git apply --check` against that working tree and against the index. Apply only
after Hooman's explicit go: `git apply workspace/notes/fork_patch_2026-09-24.diff` in the root, then
rebuild the Tyger image and re-run. Files: `raw.py`, `results.py`, `tyger_recon.py`, **new `specfit.py`**
(vendored verbatim from `notes/calspec_package_2026-09-16/specfit.py`, XeCS 2026-09-24; re-vendor if XeCS
updates it — the CSV and the Tyger path must run identical code).

v1 of this patch (same evening) kept Steve's two-Lorentzian spectrum fit and only corrected the
basis angle (fitted `RBCphase[0] − TPphase[0]` moved to k0) plus a Δf/fRBC gate. v2 replaces the fit
itself after the cross-session comparison with XeCS (ledger session 10, F65–F69): Steve's fit measures
the right lumped phase when it converges, but the peak picker and the unbounded LM fail on 8/89
sessions and the RBC shoulder makes any unconstrained two-line fit multimodal.

## Rationale (2steve notes 04 + 05, facts F59–F69)

1. **Spectral fit = `specfit.fit_block`** (`raw.py`, replaces lines 315–438 of the old block): time-domain
   complex fit of the pooled cal-block FID with a chemical-shift prior (gas L / RBC L / membrane Voigt
   [M2] or two Lorentzians [M3], RBC 218 ± 4 ppm and membrane windows from the FITTED gas line, bounded
   trf). Input = the kept FIDs (Steve's first-100 > 2× last-100 filter, no `exp(−n/200)` apodisation),
   `first_sample = killpts` so phases stay referenced to ADC sample 0 = TE, `t_k0 = killpts·dtdyn`.
   Per-rep chunks (`nsmpperusimg/npts` lines) give the within-block phase SD. Returns directly what
   the split needs: `dphiRBCTP` = RBC − lumped-membrane phase AT THE IMAGE k0 (the old `2πf·TEeff`
   was the t = 0 intercept, missing `2πΔf·TE` ≈ 74°, F59), `RBCTPratio` = a_RBC/|mem1+mem2| (the ratio a
   two-component k0 split reproduces by construction, F66; the literature's scalar ratio =
   `RBCTPratio / F_lump`), `fRBC`, `fTP` (lumped centroid). Legacy fields (`TEeff`, `RBCphase`, …) are
   filled for old callers; `sspect/sspectfit` for the GUI come from the model.
2. **Gate** = specfit's `valid` (SNR ≥ 20, RBC/mem ≥ 0.05, per-rep SD ≤ 30°), plus a stability refit
   (one more leading sample dropped, gas line fixed): reject when the k0 angle moves by > 10° (F69;
   catches 042DR at 31°), plus |sin Δφ| ≥ 0.3. On rejection the params stay empty and the existing
   fork fallback keeps the dissolved image unsplit (`rbc_tp_separated=0`) with the reason in the log.
   Same rule as the offline re-split (`resplit_rbctp.py --fit xecs`, `STAB_MAX_DEG`).
3. **Sample handling**: killpts = 2 is kept — ADC sample 0 is a ~1.8× transient and sample 1 overshoots
   by median 7 % (> 10 % in 39/94 blocks); keeping sample 1 moves the M3 decomposition and the angle by
   ~6° (F68). The offline CSV must use the same drop (requested from XeCS).
4. **Stop criterion** (`results.py`): masked complex sum → sweep ph at 0.001 rad, take the phase
   minimising |R − target| among phases with both masked sums positive (the old first-sign-change
   test on whole-volume sums fired at the ΣaTP = 0 pole in 60/86 sessions). The split is linear, so
   the sums at every ph follow from one masked complex sum: no per-phase volume work.
5. **Provenance in the output** (`tyger_recon.py`): `rbc_tp_dphi_deg`, `rbc_tp_sin_dphi`,
   `rbc_tp_df_hz`, `rbc_tp_target`, per-bin `rbc_tp_ph_rad` / `rbc_tp_R`, plus `rbc_tp_model_used`,
   `rbc_tp_ratio_scalar`, `rbc_tp_F_lump`, `rbc_tp_snr_diss`, `rbc_tp_dphi_rep_sd_deg` next to
   `rbc_tp_separated`; `results.rbc_tp_split` holds the same dict. Maps are on the LUMPED scale.
6. **Gas magnitude** (`results.py` + `tyger_recon.py`): single-channel `|F·b|` per bin exported as
   a second NdArrayFloat item `gas_phase_magnitude` alongside `real(F·b)` (note 06: `real()`
   attenuates the moving rim outside b's fixed mask). nch > 1 untouched (Σ|b| deferred).

## Before applying — check

- Downstream readers of `output.mrd` (`dyn_recon.py` / XeCS `post_process.py`) must key items by
  meta name and ignore the new `gas_phase_magnitude` item (a reader that takes "the first
  NdArrayFloat" would still get the gas image, since it is appended after it).
- The marginal class (v2 2023 sessions, |sin Δφ| 0.4–0.7: carrier between the lines, F63) passes
  the gate; the 0.3 threshold is a choice, not a measurement.
- Smoke test of the vendored path on `input.mrd` (045VS 54.3°, ratio 0.393; 004DS 55.3°; 042DR 59.9°)
  matches `specfit` on the XeCS cache with the same drop; differences vs the first CSV are the
  sample-1 issue above. Not runnable end-to-end locally (CUDA) — a Tyger rerun of 045VS is the test.
- Nothing here changes the gridding, calcb, or the 8-channel path.
