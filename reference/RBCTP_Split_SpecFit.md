# RBC/TP split — the spectral fit that feeds it, and why (2026-09-24/25)

**Status:** VALID; the fit is live in fork commit `6651591` (`diaphragm-recon`), validated on Tyger (F70).
Companion docs: `notes/fork_patch_2026-09-24.md` (what changed in the fork), `../2steve/04, 05, 07` (what Steve
should change and why), `reference/OnePointDixon_Review.md` (the field), `reference/DP_RBCTP_Split_Noise_Diagnostics.md`.
Facts F59–F70. Scripts: `helpers/resplit_rbctp.py` (C36), `helpers/specfit_template.py` (C37), `helpers/_fit_xcheck.py`.

## 1. What a one-point split needs from the spectrum
The dissolved image at a voxel is `z = a_RBC·e^{iφ_RBC} + a_TP·e^{iφ_TP}` at the time the image's k0 was
sampled (t_k0 = TE + killpts·dt_img). The 2×2 split needs exactly three numbers per session:
1. **Δφ_k0 = φ_RBC − φ_TP at t_k0** — the basis angle. Only the difference matters; the global phase sweep
   absorbs the common term.
2. **the whole-lung RBC/TP target** that fixes the global phase (Kaushik 2016; F60: the imaging ratio then
   equals the spectroscopic one by construction).
3. **conditioning**: noise gain 1/|sin Δφ|; corr(aRBC, aTP) = −cos Δφ in pure noise.

Steve's code took Δφ from `2π·f·TEeff`, the phase-vs-TE intercept at t = 0 (F59): missing 2πΔf·TE ≈ 74°,
basis 6–32° apart on data that sit 45–90° apart → noise gain 2–10, TP negative in the lung in 54/86 sessions.
The k0 angle = fitted first-sample phase difference + 2πΔf·killpts·(dt_img − dt_spec) = −13° (F61, not −1°).

## 2. The fit now used: `specfit.py` (from the XeCS calspec work, vendored verbatim)
Time-domain complex fit of the pooled cal-block FID, S(t) = Σ a e^{iφ} e^{2πift} e^{−πw_L t}[Gaussian],
t from ADC sample 0 = TE. Lines: gas L / RBC L / membrane Voigt (M2) or L+L (M3, Robertson 2017), bounded
`least_squares` trf, chemical-shift prior referenced to the FITTED gas line (RBC 218 ± 4 ppm, mem1 198–205,
mem2 193.5–198, widths ≥ 6/8 ppm, non-overlapping windows). Samples 0 and 1 dropped (F68). Output at t_k0:
`dphi_k0_deg` (RBC vs lumped membrane mem1+mem2), `ratio_lumped` = a_RBC/|mem1+mem2|, `ratio_scalar` =
a_RBC/(a1+a2), `F_lump`, Δf, SNR, per-rep phase SD, stability (angle change when one more sample is dropped),
`valid`/`reason`. M3 unless degenerate (a1/a2 ∉ [0.3, 3], |f1−f2| < 2 ppm, rep-SD blow-up) → M2.

### Why it beats Steve's frequency-domain two-Lorentzian fit (`raw.py:314–392` original) — the four reasons
1. **Never picks the wrong peak.** The gas line is a fitted component and the reference; the carrier position
   is irrelevant. Steve's 4-peak picker locked onto noise on 7/89 (v2 2023 blocks with no dissolved signal) and
   the unbounded LM returned a negative ratio once (042DR). Five of note-04's 19 "valid" sessions were such
   garbage fits whose ratio happened to match a meaningless target (F64).
2. **Bounded, with a prior.** RBC is a shoulder on the membrane line at RBC/mem 0.1–0.3; an unconstrained
   two-dissolved-line fit is multimodal (4/5 high-SNR blocks collapsed without the prior, F67).
3. **No apodisation, no baseline hack.** Steve multiplies each FID by exp(−n/200) (12 ms) before the FFT —
   an extra 27 Hz (1.5 ppm) Lorentzian broadening on every line (gas 2×) that biases the a·w widths; plus an
   edge-mean baseline subtraction that cannot remove the gas tail. The time-domain fit works on the raw FID.
4. **Carries its own validity** (SNR, rep SD, stability, model), so the recon can refuse to split instead of
   splitting garbage. Steve's fit returns a number no matter what.

What Steve's fit gets right: when it converges, its RBC − TP phase at the first sample equals the M3
LUMPED-membrane phase within +3° median (MAD 8°, n = 70) — the two-line fit is the right *estimator* for a
two-component split; M2's single-Voigt phase is 13–16° off because the two membrane lines sit −105° apart (F65).

## 3. The ratio scale (F66) — read this before comparing to literature
The split returns a_RBC and |mem1+mem2| by construction, so its target must be the LUMPED ratio
(×1.5–1.8 the scalar ratio a_RBC/(a1+a2) that Robertson/Duke/Bier report). Every aRBC/aTP map from the fork
and from `recon_resplit_xecs.mat` is on the lumped scale; divide the map ratio by the session's `F_lump`
(median 1.48; in the .mat and in the output.mrd meta `rbc_tp_F_lump`) for the literature scale. Steve's a·w
ratio sits ×1.15 below the scalar and ×1.6 below the lumped ratio.

## 4. Conventions that bit us
- Phase origin: `specfit` phases are at ADC sample 0 = TE (verified by construction; its early docstring said
  otherwise). Steve's are at his first kept sample (TE + killpts·dt_spec = TE + 120 µs). Image k0 = TE + 10 µs.
- Sign: both e^{+2πift}, RBC the higher frequency, Δφ = RBC − TP, advancing with t.
- ADC samples 0 AND 1 are transients (s0 ≈ 1.8×, s1 +7 % median, +10 % in 39/94 blocks) → killpts = 2 is
  right; the offline CSV and the container must drop the same samples (F68).
- v2 2023 protocol: carrier between the lines, TE 600 µs, weaker pulse → Δφ_k0 ≈ 35°, noise gain ≈ 1.7.
  Physics, not code (F63).
- MRD cal lines == twix cache sample-for-sample, line-for-line (checked 045VS): no conversion offset.

## 5. Where things are
| what | where |
|---|---|
| fit code (vendored) | root `specfit.py` (fork 6651591); copy `notes/calspec_package_2026-09-16/specfit.py` |
| per-session inputs (75 blocks) | `notes/calspec_package_2026-09-16/resplit_inputs_2026-09-24.csv` |
| FID cache (read-only, XeCS) | Ext `/Volumes/HoomHamExt/Work/Codes/2026_XeCS_Recon/calspec/cache/<tag>.npz` |
| offline re-split | `helpers/resplit_rbctp.py --fit xecs` → `outputs/resplit_2026-09-24/resplit_summary_xecs.csv`, Ext `d/recon_resplit_xecs.mat` (68) |
| Steve-fit re-split (reference) | `--fit steve` → `resplit_summary.csv`, Ext `d/recon_resplit.mat` (79) |
| cross-checks | `outputs/resplit_2026-09-24/{fit_xcheck.csv, template_study*.csv}`, `fig/` |
| Tyger validation | Ext `Work/Codes/2026_ASAP_Recon/tyger_specfit_2026-09-24/{045VS,041WF}/d/` |
| new-image spec | `pipeline/recon_codespec_6651591.yml` |

## 6. Open: RBC-weak blocks and the template solve (C37)
Even with the prior, blocks with scalar RBC/mem ≲ 0.15 move 10–30° with the line subset (042DR, 047SS, 040RP…):
the 16-parameter M3 decomposition is what flips. `helpers/specfit_template.py` tests the alternative: freeze
the line SHAPES at the cohort values (per protocol) and solve one complex coefficient per line by linear least
squares on a global grid over shifts (and, `--adaptive`, widths and the mem2/mem1 ratio) — no initial guess,
no local minima. First result (fixed shapes): 2–4× more stable than specfit under line-subset and sample
drops, agrees with specfit to −1.6° median / MAD 4° on stable blocks, but disagrees > 10° on 19/74 where its
residual is high (frozen shapes misfit those subjects). The adaptive grid is the next test; the decision rule
(when to switch) follows from it. Results appended below when done.

### 6b. Result (2026-09-25 02:50, F71)
- Adaptive grid (widths × a2/a1 × ±2 ppm shifts) did not lower the residual (median 0.08 either way) and made the
  flagged blocks LESS stable (odd/even p90 4.8° → 14.1°): the 10–30° template-vs-specfit disagreements on ~1/5 of
  the blocks are not a line-shape misfit; they sit in the data's own uncertainty (noise, wash-in, model error) and
  no estimator can be called right there without ground truth.
- Verdict: the fixed-shape template is a **precision** fix, not a proven **accuracy** fix. Rule proposed (not
  deployed): specfit M3 when its stability passes; fixed-shape template only as a flagged fallback where specfit is
  unstable AND template odd/even + halves stability < 10° AND relative residual < 0.2. Rescues 042DR, 021JM, 023DB;
  047SS 10-23, 040RP 12-11, 037GD stay unsplit.
- Next test for accuracy: simulation ground truth (synthetic FIDs from the cohort M3 parameters, varied RBC/mem,
  SNR, membrane phase, wash-in; bias and variance of both estimators). Not run.
