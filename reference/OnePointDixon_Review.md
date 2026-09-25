# One-point Dixon for hyperpolarized 129Xe gas-exchange MRI — a critical literature review

**Status:** DRAFT v1, 2026-09-24 (Fable 5.1 session). Tags: `learn`, `literature`.
**Ask (Hooman):** read the 1-point Dixon papers (Duke canon, consortium, critics), report, give an opinion on
the approach the field converged on.
**Sources read in full:** Kaushik 2016 MRM; Robertson 2017 MRM; Wang 2018 MRM (3T); Wang 2018 Thorax;
Hahn 2018 MRM; Ruppert 2019 MRM (Rizi lab, HH co-author); Niedbalski 2021 MRM (consortium position paper,
gas-exchange + calibration sections); Willmering 2021 MRM; Niedbalski 2022 MRM; Bechtel 2023 MRM;
Collier 2021 MRM; Mummy 2024 Front Med; Leewiwatwong 2025 MRM (multi-site reference values); abstracts
of Wang 2017 Med Phys, Garrison 2023 Radiol CTI, Leewiwatwong 2026 MRM, García Delgado 2026 MRM.
Secondary set (Qing 2014, Norquay 2015/2017, Plummer 2023/2024, Collier 2024, Lu 2024, Mummy 2025,
Costelle 2025, Niedbalski 2019/2020/2023, Kern 2019, Willmering 2025 phantom, Marshall 2021 review, Teague 2024,
Xie 2019, Ruppert 2016 CSSR, Kaushik 2013/2014, Shim 2026 editorial) extracted by two sub-agents; see §8.
(`Fain_2017_Hahn` in the PDF library is neonatal 1H UTE, not xenon — ignored.)
Also inspected: Duke public pipeline `TeamXenonDuke/xenon-gas-exchange-consortium` (main, 2026-09-24),
Steve's `raw.py` / `results.py` in this repo.
Local PDFs: `~/Hooman/Work/Literature/PDFs/`. Plain-text extracts were session-scratch (not kept).

---

## 1. What the field converged on (the "consortium 1-point Dixon")

| Item | Converged value (3T) | Source |
|---|---|---|
| Sequence | interleaved 3D radial, gas view then dissolved view, 1000+1000 projections, 64 pts/arm, 0.64 ms readout, 16 s | Niedbalski 2021 Table 6 |
| TR (same-frequency) / flips | 15 ms; 0.5° gas / 20° dissolved (TR90,equiv ≈ 249 ms) | consortium; Ruppert 2019 |
| Fast variant | 5.4 ms / 12° (2100 proj, 11 s) or 8.5 ms / 15° (1200 proj, 10 s) | Niedbalski 2022; Leewiwatwong 2025 |
| Dissolved RF | 1-lobe Hanning-windowed sinc 0.65–0.69 ms (3T); 3-lobe 1.15–1.25 ms (1.5T); centred 218 ppm → **208 ppm since 2025** | consortium; Bechtel 2023; Leewiwatwong 2025 |
| TE | "TE90" from calibration; expected 0.45–0.50 ms at 3T (0.8–1.1 ms at 1.5T); fallback 0.47 (1.0) ms | consortium Eq. 4 |
| Calibration | 500 dissolved FIDs at TR 15 ms + 20 gas FIDs; discard first 100; fit next 67 in time domain (Lorentzian RBC/gas, Voigt membrane) → RBC:M, Δf, phases, flip | consortium §3.3; Bier 2019 |
| Decomposition | global zero-order phase rotation so that Re/Im (or their ratio inside the lung mask) match spectroscopic RBC:M; then per-voxel subtraction of the gas-image phase map ("B0 correction") | Kaushik 2016; Duke pipeline `img_utils.dixon_decomposition` |
| Gas-contamination | keep <10 % of dissolved signal by pulse tuning; if not, Hahn dual-echo or Willmering spectroscopy-informed subtraction; acceptable threshold ≈ 9 ± 3 % | Wang 2018; Willmering 2021; Leewiwatwong 2026 |
| Post-hoc scalings | T2* (global, from linewidth: 0.99 ms M, 1.04 ms RBC), flip-angle/TR (MOXE), Hb (MOXE, Hb0 = 14 g/dL), lung-volume, 218→208 ppm factor 0.89 | Leewiwatwong 2025; Bechtel 2023; Duke pipeline |
| Quantification | voxel-wise M/gas, RBC/gas; colour bins of 1 SD (Box-Cox) of a pooled healthy 18–30 y reference; report defect/low/high % | Wang 2017; Leewiwatwong 2025 |
| Reference (208 ppm, 3T, n = 50, 3 vendors) | vent 0.71 ± 0.14; M/gas 0.97 ± 0.27 ×10⁻²; RBC/gas 0.48 ± 0.20 ×10⁻²; RBC:M 0.49 ± 0.11 (≈0.55 at 218 ppm) | Leewiwatwong 2025 |

Duke pipeline decomposition, as shipped (checked in source):
```python
diffphase     = correct_b0(image_gas, mask)            # gas-image phase, zero-mean inside lung
desired_angle = arctan2(rbc_m_ratio, 1.0)              # assumes RBC ⟂ membrane at k0
current_angle = angle(sum(image_dissolved[mask]))
image_dixon   = image_dissolved * exp(1j*(desired_angle - current_angle)) * exp(-1j*diffphase)
image_rbc, image_membrane = imag(image_dixon), real(image_dixon)   # sign-flipped if mean < 0
```
No per-voxel Δφ, no readout-phase model, no T2* difference between compartments, no third peak.

## 2. How it got here (2013 → 2026)

- **2013–2014 (Duke, 1.5T).** Total dissolved-phase 3D radial imaging (Kaushik 2013), then spectroscopy showing RBC:barrier 3.3× lower in IPF, r = 0.89 with DLCO (Kaushik 2014). The decomposition problem is born.
- **2014 (UVA).** Qing et al. do the "proper" thing: hierarchical IDEAL, 3 TEs (longest 3.98 ms), needs very high polarization.
- **2016 (Kaushik).** 1-point Dixon in humans: one sub-ms TE, phase-shift + ratio prior + gas-phase B0 map. Honest limitations section: (i) the two lines are ~140 Hz wide and overlap 30–40 %, (ii) the Dixon condition holds only at k0, (iii) TE90 is not 1/(4Δf): predicted 735 µs, measured 410–470 µs healthy and 160–285 µs IPF; phase-vs-TE intercepts of 33° (healthy) and 55° (IPF) attributed to "phase evolution during the RF pulse" that "must be investigated by simulating the Bloch equations". That investigation never appeared in the literature I read.
- **2017 (Wang, Med Phys).** Binning maps against a healthy reference (n = 13): mean ± SD bins, defect/low/high %. Healthy: barrier 4.9 ± 1.5 ×10⁻³, RBC 2.6 ± 1.0 ×10⁻³; 15 % of healthy RBC voxels already "low".
- **2017 (Robertson).** Time-domain fitting shows **three** dissolved resonances (RBC, barrier 1 ≈ 201 ppm, barrier 2 ≈ 196 ppm, ~100° out of phase with each other). 2-peak fits mis-assign barrier signal to RBC; RBC:barrier changes from 65 % → 47 % reduction in IPF; the previously reported −2.4 ppm RBC shift in IPF shrinks to −0.7 ppm. Paper ends: the third peak "necessitates investigating the impact on the 1-point Dixon methods… separate imaging… may be improved by multi-echo… IDEAL".
- **2018.** Wang moves to 3T: T2* ≈ 1.1 ms, TE90 = 0.47 ± 0.02 ms (vs 0.37 predicted), gas contamination 3 %, RBC image SNR **2.9 ± 1.5**, spectrum appended in-breath so RBC:barrier is measured at the imaging TR/flip. Hahn (Wisconsin) shows gas contamination mimics disease and adds a second echo to remove it. Ruppert (Rizi lab) establishes TR–flip equivalence: the dissolved contrast is set by RRF = (1 − cos α)/TR; different groups' ratios are only comparable after conversion.
- **2021.** Consortium position paper standardizes the acquisition and lists the open problems in its own words: "high degree of undersampling (~15 % of Nyquist), excitation pulse design, disregard of local phase variations, and chemical-shift-induced phase evolution during the radial read-out"; "the extent to which this impacts the derived quantitative imaging metrics has not been rigorously examined". Willmering (Cincinnati) gives a contamination correction that needs no extra echo. Collier (Sheffield) publishes the four-echo EPSI alternative.
- **2022–2024.** Fast Dixon via TR–flip equivalence (Niedbalski). Hb adjustment (Bechtel): RBC:M moves up to 20 % across the normal Hb range, 40 % in one subject with Hb 10.2; the 218-ppm pulse under-excites the membrane (×0.89 ± 0.11) → recommendation to excite at 208 ppm. Age/sex/BMI dependence (Plummer, Mummy, Collier): RBC:M −0.05 per decade, males +0.17. Lung-volume dependence (Garrison, UVA/IDEAL): M/gas vs volume r = −0.97; COPD-vs-healthy differences in M/gas and RBC/gas become non-significant after volume correction.
- **2025–2026.** Multi-vendor reference distributions (Leewiwatwong): site-to-site means within ±8 %, Cohen's d ≤ 0.36; sex differences of +14 % M/gas (F) and +21 % RBC/gas (M) survive Hb correction. Leewiwatwong 2026: contamination acceptable up to ~9 %; both correction methods need a **one-time empirical phase calibration** (Hahn ≈ 230°, Willmering ≈ 49°, hard-coded as `optimized_conta_phase = 49.9` in the pipeline).

So: the field converged on the *acquisition and reporting conventions* (2016–2021) and has since been adding global scalar corrections (T2*, TR–flip, Hb, volume, excitation offset, contamination, bias field) on top of an unchanged decomposition.

## 3. The decomposition, written out

At k0 the dissolved signal is
S = a_M e^{iφ_M} + a_R e^{iφ_R}, with Δφ = φ_R − φ_M = 2πΔf·TE_eff + φ_RF(pulse, exchange, B0).

1-point Dixon needs Δφ = 90° and a global rotation that is fixed by the prior a_R/a_M = RBC:M (spectroscopy). Two consequences that the literature under-states:

1. **The global ratio is not a validation.** The imaging RBC:M is forced to equal the spectroscopic RBC:M inside the mask by construction (Kaushik Eq. 3; Duke `desired_angle`). Every paper that reports "imaging ratio agreed with spectroscopy" for 1-point Dixon is reporting the constraint, not a measurement. Only regional deviations from the global ratio carry information. In Collier's four-echo method the agreement (r = 0.90, bias 10.9 % explained by TR/flip) is a genuine check because nothing was imposed.
2. **Conditioning.** For a general Δφ the 2×2 inversion has noise gain 1/|sin Δφ| and corr(a_R, a_M) = −cos Δφ (HH's derivation in `DP_RBCTP_Split_Noise_Diagnostics.md`). Duke sidesteps this by *assuming* Δφ = 90° and taking Re/Im: no noise amplification, but any true departure from quadrature becomes a **bias** (cross-talk of membrane into RBC ∝ cos Δφ) rather than variance. Steve's `results.py` does the opposite: it models Δφ explicitly, so a wrong or small Δφ produces the amplified, anti-correlated maps HH found (60/86 sessions invalid). Neither pipeline measures Δφ per voxel.

Where Δφ actually comes from is the least-understood part of the whole method:

| Observation | Value | Source |
|---|---|---|
| Naïve TE90 = 1/(4Δf) at 1.5T | 735 µs | Kaushik 2016 |
| Measured TE90, healthy / IPF, 1.5T | 429 ± 30 / 240 ± 60 µs | Kaushik 2016 Table 1 |
| Phase-vs-TE intercept, healthy / IPF | 33 ± 4° / 55 ± 10° | Kaushik 2016 |
| Predicted vs measured TE90 at 3T | 0.37 vs 0.47 ± 0.02 ms | Wang 2018 |
| Barrier-1 vs barrier-2 phase (same excitation) | ~100° apart although only 83 Hz apart | Robertson 2017 |
| Barrier-1 phase, IPF vs healthy | +20° | Robertson 2017 |

The intercept is not just the half-pulse delay (at 1.5T a 1.2 ms pulse would give ~73°, not 33°); it depends on disease and on which barrier sub-peak dominates. The consortium formula TE90 = TE + (90° − Δφ_meas)/(360°·Δf) absorbs whatever it is into a per-subject calibration. That is pragmatic, but it means "TE90" is a fitted quantity of a 2-peak model applied to a 3-peak spectrum, and the phase it targets is only defined for the whole-lung average.

**Phase evolution during the readout.** Δf ≈ 690 Hz at 3T; with 10 µs dwell the RBC–membrane phase advances ~2.5° per sample (Collier quotes ~10°/point at 1.5T with 16 kHz BW). Wang 2018: <100° over the ~32 samples that carry dissolved signal before it hits the noise floor. So the two compartments are in quadrature at k0 and rotate away from it across k-space; the outer k-space of the "RBC" image is contaminated by membrane signal at a phase that depends on |k|. This is a chemical-shift blurring/cross-talk term that the consortium acknowledges and nobody has quantified for the binning metrics. Brodsky-style k-space chemical-shift modelling (used by Collier, Kammerman, Zanette) is the fix; it costs nothing in acquisition time but needs the model.

**Effective resolution.** Dissolved signal is at the noise floor after ~0.3 ms of a 0.64 ms readout (Wang 2018). Nominal 6.25 mm voxels (64³ over 40 cm) are decorative; the membrane/RBC PSF is set by T2* ≈ 1 ms and the density compensation, probably 15–20 mm FWHM. Add 1000 projections ≈ 15 % Nyquist. The reference distributions therefore describe a strongly blurred, streak-contaminated quantity, and "defect %" is a function of the PSF as much as of physiology. This is also why 1.5T (T2* ≈ 2 ms) and 3T references are not interchangeable (Leewiwatwong 2025 says so).

## 4. Catalogue of critiques and their status

| # | Issue | Evidence | Addressed? |
|---|---|---|---|
| 1 | Gas-phase contamination mimics disease | Hahn 2018 (artifacts ≈ disease; correction removes them); Willmering 2021 (9.5 ± 4.8 % typical at 3T, 5.9–17 % across coils); Leewiwatwong 2026 (limit ≈ 9 %) | Mostly. Both corrections need an empirical platform-specific phase (230° / 49°) — a symptom that the contamination model is incomplete (RF-amp nonlinearity, receive-chain phase). GE could not meet <10 % in the 2025 reference study and used the dual echo. |
| 2 | Quadrature only at k0; phase evolves through readout | Kaushik 2016; Wang 2018 (<100°); consortium 2021 ("not rigorously examined"); Collier 2021 (10°/pt) | **Open.** No published quantification of the bias on M/gas, RBC/gas, defect %. |
| 3 | "TE90" ≠ 1/(4Δf); intercept unexplained, disease-dependent | Kaushik 2016; Wang 2018 | **Open.** Handled per-subject by calibration; physics never resolved. |
| 4 | Third dissolved resonance; barrier 1 and 2 ~100° apart | Robertson 2017 | **Ignored** by the imaging pipeline (2-peak Dixon; Voigt membrane in the fit is a lineshape patch, not a phase model). Robertson himself suggested multi-echo. |
| 5 | B0 correction uses the gas-image phase | Kaushik 2016; Hahn 2018 (calls it a surrogate assuming spatially invariant RF phase) | Partly. Assumes the gas voxel (airway-weighted, different PSF, ~40× lower flip, single resonance with its own susceptibility shift) shares the dissolved voxel's phase. Reasonable at first order; unverified. |
| 6 | Global RBC:M prior; regional ratio assumed to be modulated only by data phase | Kaushik 2016; Hahn 2018 (spatial invariance of flip and T2* "only valid for gross effects") | **Structural.** Not fixable without more echoes. |
| 7 | The prior itself depends on the fitting model | Robertson 2017 (2- vs 3-peak: RBC:barrier differs ~30 % in healthy); Leewiwatwong 2025 vs Collier (RBC shift 218.2 vs 216.9 ppm from Voigt/time-domain vs Lorentzian/frequency-domain) | **Open.** The consortium chose one fit; cross-group ratios are not comparable without re-fitting raw FIDs. |
| 8 | Unequal excitation at 218 ppm | Bechtel 2023 (membrane at 89 %; empirical factor 0.89 ± 0.11) | Fixed going forward (208 ppm). All 2016–2024 disease thresholds carry the bias; the "constant" has a 12 % SD. |
| 9 | Ratios depend on TR and flip | Ruppert 2019; Niedbalski 2022; Bechtel 2023 | Addressed by TR90,equiv; but every cross-group conversion runs through MOXE with healthy literature geometry (δ/d = 0.185, tX = 1.6 s). |
| 10 | Hb dependence | Bechtel 2023 (≤20 % RBC:M in normal Hb range; 40 % at Hb 10.2) | Correction exists; validated cross-sectionally on n = 18 with a CR of 0.09; sensitive to δ/d (±8–15 % in ILD). Reference distributions in 2025 are *not* Hb-adjusted. |
| 11 | Sex, age, BMI | Mummy 2024 (RBC:M −0.05/decade, M +0.17, BMI −0.07/10); Plummer 2023; Collier 2024; Leewiwatwong 2025 (F: M/gas +14 %, RBC/gas −17 %, persists after Hb) | Partly. Consortium recommends one 18–30 y mixed-sex reference for all adults and "interpret with age-dependent norms". Any 60-y-old woman will have "RBC defects" by construction. |
| 12 | Lung-volume dependence | Garrison 2023 (M/gas r = −0.97 with volume; COPD differences vanish after correction); Qing 2014b (RV→TLC: TP/gas ×2, RBC/gas ×3.5); Collier 2024 (EIVt→TLC: RBC/gas −33 %, M/gas −39 %, RBC:M only −9 %); Kern 2019 | Partly. Duke pipeline now has `vol_correction`; the 2025 reference admits inflation could not be verified at any site. The single largest confounder in the literature I read. RBC:M is the one ratio that is inflation-robust. |
| 12b | Calibration spectrum samples a moving target | Ruppert 2016 (RBC amplitude swings ~15 % peak-to-peak at ~1 Hz; 12 % at TLC vs 8 % at RV); Norquay 2015 (plasma↔RBC exchange τ ≈ 12 ms, i.e. faster than TR) | Consortium averages 67 FIDs ≈ 1 s (≈ one cardiac cycle) — adequate for the mean, but the imaging data average a different cardiac/TR history than the spectrum. |
| 12c | Gas and dissolved images have different sampling | Niedbalski 2023 (gas ~fully sampled vs dissolved ~16 %; "could artificially elevate membrane or RBC signal" after division); RBC SSIM between two same-day scans 0.14 | Not addressed; the ratio maps divide two images with different PSF and streak structure. |
| 13 | RBC chemical shift moves with sO2 | Norquay 2017 (5.5 ppm over sO2 0.1–1; ~1 ppm drop in a 35-s breath-hold); Kaushik 2014/Robertson 2017 (IPF shift) | Small effect on TE90 phase at TE ≈ 0.45 ms (1 ppm ≈ 0.6°), but it is a physiological variable buried inside the calibration. |
| 14 | SNR | Wang 2018 RBC SNR 2.9 ± 1.5; Niedbalski 2022 RBC 8 ± 3; reference study excludes RBC SNR < 5 | Not addressed. At SNR 5 the noise CV (20 %) is half the reference CV (42 %); bin membership of a voxel is coin-flip-grade for the RBC map. Site-specific SNR alters defect %. |
| 15 | Undersampling / PSF | consortium ("~15 % of Nyquist"); Lu 2024 ("may require revisiting the … recommended 1-point Dixon acquisition using only 1000 radial views"); Plummer 2023 (dissolved gridded at kernel sharpness 0.1 vs 0.3 for gas, "reduced the effective resolution"); Plummer 2024 (decay-modelled CS: 3× lower RMSE in simulation, global ratio and Dixon angle unchanged) | Partly (fast Dixon 2×; CS recon in Cincinnati). Chemical-shift phase along the readout is still not in any forward model. |
| 15b | Multi-vendor deployment cost | Mummy 2025 XeCITE: GE contamination "could not be reliably suppressed below 20 %" → second echo; hardest part was "harmonizing the reading and k-space scaling of radial trajectories"; fast Dixon "not compatible with the Hahn multi-echo gas suppression"; 5 of 8 sites failed first qualification on gas-exchange SNR/off-resonance | Being worked through; shows how much of the standard is scanner plumbing. |
| 15c | QA phantom | Willmering 2025: acetone/DMSO + Fe(acac)₃ at 212/194 ppm, T1 1.25 s, field-independent; validated with vendor 3-echo mDixon, not 1-point | The harmonization phantom for 1-point Dixon was itself decomposed with multi-echo. |
| 16 | Reference-relative, unitless output | Wang 2017 onward | By design. 2.3 % of healthy voxels are "defect" whatever the physiology. |
| 17 | Repeatability | Garrison 2023 ICC 0.88 M/gas, 0.71 RBC/gas; Hahn 2019; CR RBC:M 0.09 | Acceptable for group studies; individual-patient thresholds are within the repeatability band. |

## 5. The alternatives and why they lost

- **Hierarchical IDEAL (UVA, Qing 2014).** 3 TEs up to 3.98 ms at 1.5T, image-domain decomposition; needs polarization the field did not have; at 3T the TEs exceed T2*. Won the physics, lost the SNR budget.
- **Four-echo flyback 3D radial EPSI (Sheffield, Collier 2021/2024).** ΔTE 0.7 ms at 1.5T (effective-NSA optimum), k-space chemical-shift model (Brodsky/Wiesinger) with per-compartment T2* and gas as a fourth species, so contamination is *modelled* not subtracted, no ratio prior needed for the decomposition, and cardiogenic oscillations come for free from the k0 samples every 40 ms. Sensitivity test: replacing subject-specific frequencies/T2* by population means changes ratios by only 2–3.5 % — so the calibration spectrum is nearly optional. Cost: 2 cm isotropic voxels, 332 projections, 40° flip (2021); 2024 version is standalone (frequency-tailored pulse, 934 projections, TR 15 ms / 22°, TEs 0.57–2.67 ms, 20 dummy FIDs provide the calibration), n = 62 normative equations, intra-session ICC RBC:M 0.98. Proposed 3T variant: ΔTE 0.35 ms, two interleaves with two echoes each. Collier 2024 states the obvious gap: a regional comparison of 1-point Dixon and multi-echo maps in the same subjects "is warranted" and does not exist. This is the method I would pick if SNR permitted.
- **Model-based IDEAL with R2* maps (Wisconsin, Kammerman 2020).** Dissolved R2* as an IPF biomarker; same family as Collier.
- **Spiral IDEAL + MOXE (Toronto, Zanette 2018/2019).** Parametric gas-exchange maps; spiral trajectory — directly relevant to ASAP.
- **3D CSI (Wuhan).** Full spectrum per voxel in ~9 s at coarser resolution; consortium mentions it as "likely to provide improved separation".
- **CSSR/XTC spectroscopy (UVA, Rizi).** Global but physically modelled (MOXE); the consortium and Duke now lean on MOXE for every correction while imaging with the method MOXE cannot describe (2-peak, single TE).

They lost for non-physics reasons: a single-echo sequence is trivial to port across vendors; RBC SNR at 3T is already marginal so extra echoes were unaffordable; the Duke pipeline, ISMRMRD converters, reference values, FDA (ventilation) and Polarean support created a network standard. The consortium paper is explicit that 1-point Dixon was chosen for "relative simplicity and robustness", and that multi-point methods "have not been as well-published or broadly disseminated".

## 6. My verdict

1. **It is a good engineering compromise for 3T and a weak measurement of RBC vs membrane.** Given T2* ≈ 1 ms, RBC SNR ≈ 5, and a 10–16 s breath-hold, a single sub-ms echo with a global ratio prior is close to the only thing that produces two maps at all. The field's convergence is a convergence on *feasibility and standardization*, not on validation. The multi-site reference paper demonstrates precision (sites agree to ±8 %) and says nothing about accuracy, because there is no ground truth for regional RBC vs membrane xenon and the only global check is imposed by the method itself.
2. **The physics debt is concentrated in three places nobody closed:** the RF/exchange-induced phase intercept (why TE90 is 30–60 % shorter than 1/(4Δf) and shorter still in fibrosis), the phase evolution across the readout (k-space cross-talk; the consortium's own "not rigorously examined"), and the third resonance whose phase differs by ~100° from its neighbour and shifts by 20° in IPF. Together they mean the "membrane" image is a projection of two vectors of unknown, disease-dependent relative phase onto an axis defined by a 2-peak fit. I would not rule out that part of the "elevated membrane uptake in ILD" signature is a projection effect — it is at least a hypothesis that no paper tests.
3. **The correction stack is drifting toward a tower of global scalars** (T2*, TR–flip via MOXE, Hb via MOXE, lung volume, 218→208 ppm, contamination phase, seven N4 passes, Box-Cox bins), each calibrated on small healthy cohorts and each applied uniformly to a voxel-wise ratio image whose noise and PSF are never modelled. Every scalar makes the healthy reference tighter and the machinery more self-consistent; none makes a voxel more accurate. The genuinely robust outputs are the ones that need no decomposition: total dissolved/gas, the spectroscopic RBC:M (repeatable, r² ≈ 0.84 with DLCO), RBC shift, and RBC oscillation amplitude.
4. **Biggest confounders outrank the decomposition errors.** Lung volume (r ≈ −0.97), sex (+14/−17 %), Hb (≤20–40 %) and age (−0.05 RBC:M per decade) each move the metrics more than the site-to-site scatter, and the standard reference is a mixed-sex 18–30 y cohort at an unverified inflation level. Consortium 1-point Dixon "defect %" in an older woman with mild anaemia is not interpretable without the covariate models that only exist as single-site prototypes.
5. **What I would keep, what I would change.** Keep the interleaved radial/spiral gas–dissolved acquisition, the in-breath spectrum, 208 ppm excitation, TR–flip standardization, and the reference-distribution reporting for ventilation and total dissolved. Change the decomposition: at minimum a k-space chemical-shift model with the measured Δf and per-compartment T2* (two echoes if SNR allows; Collier's 3T design), report per-subject Δφ and conditioning, and stop calling a ratio-constrained global agreement "validation".

## 7. Consequences for ASAP / this repo

- **Total dissolved is safe; the RBC/TP split is not** — consistent with `2steve/04`. Steve's `results.py` models Δφ explicitly (`phaseRBC − phaseTP = 2πΔf·TEeff`) and inverts the 2×2; Duke assumes quadrature and takes Re/Im. Steve's version is more honest and more fragile.
- **Check before anything else (hypothesis, not yet verified — F59):** in `raw.py:386–392`, `deltaphase` is the *intercept* difference of the phase-vs-TE fits (extrapolation to TE = 0) and `TEeff = deltaphase/(2πΔf) − killpts/spectBW + killpts/BW`. The 2πΔf·TE term of the actual echo time (`self.TE` from `alTE`) does not appear in `TEeff`, yet `results.py:324–325` uses `2π·f·TEeff` as the basis phases. Our data are 1.5 T (F54), so Δf ≈ 340 Hz and the missing term is ≈ 0.12° per µs of TE: ≈ 40° at TE 0.3 ms, ≈ 110° at TE 0.9 ms. HH measured Δφ_model = 6°, 11°, 39° from the anti-correlation of the split maps. Test: compare `RBCphase[0] − TPphase[0]` (fitted phase difference at the first spectral sample, which sits at the same TE as the image k0) with `2π·(fRBC−fTP)·TEeff`, per session, from the calibration block. If they differ by ≈ 2πΔf·TE, the basis vectors are wrong and the ill-conditioning is a code artefact, not an acquisition property — and the fix is one line. If they agree, HH's original reading (acquisition far from quadrature) stands and only a TE change or multi-echo helps. Either way, log Δφ at k0 into the output.
- **1.5 T is the friendlier field for doing it properly.** T2* ≈ 2 ms (Qing, Collier: TP 2.2 ms, RBC 1.9 ms), so Collier's four-echo design (ΔTE 0.7 ms, first three echoes inside ~2 ms) applies unchanged, and ASAP's spiral is the trajectory Zanette/Doganay used for spiral IDEAL. A two- to four-echo spiral with k-space chemical-shift modelling (gas, M, RBC; frequencies and T2* from the existing calibration block) would give a decomposition that needs no ratio prior and models gas contamination instead of subtracting it. Worth a feasibility branch once the dynamic-phantom work is stable.
- **If a 1-point split is kept for the cohort:** (i) gate on |sin Δφ| ≥ 0.3 (already proposed in `2steve/04`), (ii) subtract gas contamination from k-space using the in-breath spectrum (Willmering: works retroactively, no extra echo), (iii) apply the gas-phase B0 map as Duke does, (iv) report total DP as the whitened magnitude (E[q] = 2 test), (v) never compare RBC/TP across sessions acquired at different TR/flip without the TR90 conversion.

## 8. Secondary-source extractions (two sub-agent dossiers, condensed; numbers as reported)

### 8A. Duke-lineage papers

**Kaushik 2013 JAP (total dissolved, no Dixon).** 1.5 T GE; 1.2 ms 3-lobe sinc on RBC (+3832 Hz), TE/TR 0.932/7.5 ms, 32³, 1001 rays, 0.5°/22°. Gas contamination "~20 % of the dissolved-phase excitation" was accepted. Ratio normalization removes B1 "to first order". Transfer CV correlates with TLC (r = 0.77 supine) because 1 L was 11–17 % of TLC.

**Kaushik 2014 JAP (spectroscopy in IPF).** 3-Lorentzian complex fit, area = amplitude × FWHM. RBC:barrier 0.55 ± 0.13 HV vs 0.16 ± 0.03 IPF; r = 0.89 with DLCO; same-session repeat 6.6 %, inter-day −8 % (two subjects −24 %, −16 %). RBC:gas vs inflation r = −0.78. Already in 2014 the paper warns that the ratio depends on TR/flip, on equal excitation (6 % lower flip on RBC), on hematocrit ("not measured") and on inflation — the things "fixed" in 2019–2023.

**Wang 2018 Thorax.** 1.5 T, TE90 ≈ 0.9 ms, 6.3 mm. HV (n = 13, 33.6 ± 15.7 y): barrier 0.52 ± 0.13 ×10⁻², RBC 0.28 ± 0.06 ×10⁻², RBC/barrier 0.58 ± 0.13, BarrierHigh 3 %, RBCLow 14.6 ± 6.9 %. IPF: barrier 188 ± 36 % of reference, RBC 72 ± 25 %, RBC/barrier 0.21. DLCO r: barrier −0.75, RBC 0.72, RBC/barrier 0.94; CT fibrosis score r = 0.09–0.27 (all n.s.). Repeat over 7.7 months: within-subject CV 16 % barrier, 18 % RBC, 11 % RBC/barrier ("upper bound"). Bin definitions differ from later papers (BarrierHigh ≥ 3 SD; RBCLow > 1 SD). No B0/phase/contamination description.

**Niedbalski 2019 MRM (keyhole decay mapping).** 40–50 % global HP-magnetization attenuation over a 3600-view scan; C1 mapped by 2-key keyhole. The standard Dixon pipeline applies no spatial decay correction across its 1000+ views.

**Niedbalski 2020 JAP (oscillation keyhole on Dixon data).** Mixed 1.5 T/3 T retrospective cohort; the only explicit statement of the Dixon phase step in the Duke set: "images were phase shifted such that the ratio of the mean signals of the imaginary and real channels matched the spectroscopic RBC/barrier ratio" (one global phase, no B0 map mentioned). "Nominal ~3 mm … true resolution is ~6 × 6 × 6 mm³"; "~15 % sampling at the edge of k-space". Oscillation: HV 8.7 ± 4.1 %, PAH 4.2 ± 2.8 %, IPF 13.8 ± 3.8 %; image vs spectroscopy R = 0.85. Limitation in their words: "phase evolution continues throughout the radial readout, although the effect … has yet to be rigorously examined."

**Niedbalski 2023 NMR Biomed (single breath-hold).** Dedicated GX at TR 2.7 ms / 12°, 1900 proj (TR90 249 ms). Single-breath vs dedicated: ICC mem/gas 0.97, RBC/gas 0.99 (means), but SSIM RBC **0.14**, Dice low-RBC 0.60. "The true resolution is significantly poorer than nominal"; gas nearly fully sampled vs dissolved ~16 % "could potentially lead to artificially elevated membrane or RBC signal … as a result of scaling by the gas signal."

**Lu 2024 MRM (oscillation phantoms).** Digital phantom: osc CV 0.27/0.39/0.57 for 5000/2000/1000 views; SNR ≥ 3 needed; healthy RBC-transfer SNR 7.9 ± 5.5, CTEPH 2.3 ± 1.0. One healthy volunteer's "defect" went 29.1 % → 0.3 % by changing the key radius alone. States the 1000-view consortium acquisition "suffers from significant radial undersampling artifacts" and should be revisited (fast Dixon 2200 views). Reference SD for oscillation went 9.0 % → 3.8 % between papers by re-curating the cohort.

**Costelle 2025 JAP (circuit model).** Osc = 3.73/V′C + 6.76 (r²adj 0.53); PVR r² = 0.27 vs RHC; capillary blood volume estimated from RBC-transfer images assuming linearity with RBC_mean/RBC_ref, young reference, no Hb correction; authors: "cannot yet be used clinically".

**Mummy 2025 JMRI (XeCITE multi-centre COPD trial, terminated early).** 7 + 1 sites, 3 vendors, all 3 T. Calibration 520 FIDs / TE 0.45 / 20° / dwell 39 µs sets frequency, voltage and TE for Dixon. GE sites: second dissolved echo at 3.8 ms because contamination "could not be reliably suppressed below 20 %" (target < 10 %). Qualification: 2 sites passed first try, 5 second, 1 third; failures were gas-exchange SNR (off-resonance excitation ×3, T/R switching ×2), wrong resolution ×2. "The most challenging aspect … was harmonizing the reading and k-space scaling of radial trajectories across the three major scanner platforms." Fast Dixon incompatible with the Hahn dual-echo; recommends comparing Hahn vs Willmering (done in Leewiwatwong 2026) and moving to the 10-s/2200-view scan. Dynamic spectroscopy needs RBC SNR > 100 (DE 75–100 mL). No Dixon phase/B0 description beyond citing Kaushik.

**Shim/Woods/Mugler 2026 AJRCCM editorial.** Treats gas-exchange metrics as near-clinical; no methodological caveats.

Agent's cross-paper notes: TE drifted 0.932 → "≈0.9" (1.5 T) → 0.45–0.50 (3 T) with 1/(4Δf) = 0.73 ms at 1.5 T, so TE90 is empirical and pulse-definition dependent; excitation centre RBC (218) → "dissolved" → 208 ppm; ratio provenance changed (Lorentzian AUC → time-domain Voigt); bin definitions changed four times; reference cohorts changed five times; third peak absent from every fit model.

### 8B. Critics and alternatives

**Collier 2021 MRM (Sheffield four-echo).** See §5. Extra numbers: contamination in the dissolved interleave = 0.15 ± 0.02 of the gas interleave (≈ 0.1° off-resonance flip); gas-contamination image vs ventilation image r = 0.93 ± 0.03 pixelwise; T2* TP 2.18/2.4 ms, RBC 1.89/2.3 ms (HV/IPF); RBC shift 216.9 ± 0.6 (HV) vs 215.7 ± 0.9 ppm (IPF). Explicit: 1-point Dixon "would be ill-conditioned for the reliable detection of a second dissolved peak in the TP signal"; the 3-echo IDEAL of Qing assumed equal T2*, no contamination and instantaneous k-space ("severe image blurring" when done in image space).

**Collier 2024 MRM (normative, n = 62, 20–69 y).** RBC:M(M) = −0.00362·age + 0.60, RBC:M(F) = −0.00170·age + 0.44 (34 % vs 21 % decline over 20–70 y); M:Gas no age dependence; T2*_RBC rises with age (ρ = 0.60); RBC shift depends on sex (p = 0.001). TLC vs end-inspiratory tidal: RBC:M −8.6 %, RBC:Gas −33 %, M:Gas −39 %. Intra-session ICC RBC:M 0.983 (BA −8.3/+5.9 %), RBC:Gas 0.748. Adopting consortium TR/FA "had a negative impact on the SNR of the RBC k-space center signal".

**Qing 2014a JMRI (UVA 3-echo IDEAL).** 1.5 T Siemens; excitation midway between TP and RBC; TEs 0.74/2.36/3.98 ms; gas at 2 TEs for a B0 map; TR 19 ms, 23°/0.4°, 11 s; 7.6 × 7.6 × 17 mm³; hierarchical IDEAL with one T2*; error from unequal T2* ≤ 3–5 %. FID vs image RBC fraction agree (0.27/0.24 etc.). Repeatability 5–12 %, inflation blamed. Notes that "the signal corresponds to a range of exchange times" set by FA/TR.

**Qing 2014b NMR Biomed (CSSR vs imaging, n = 45).** RV → TLC: TP/GP ×2, RBC/GP ×3.5, RBC/TP ×1.5, septal thickness only 7–10 % (n.s.). CSSR RBC/TP healthy 0.29 ± 0.04; COPD 0.15 ± 0.05; asthma bimodal. Imaging − CSSR RBC/TP = 0.037 ± 0.036. Steady state at 23°/19 ms ≈ CSSR delay of 26 ms (6 µm wall) to 100 ms (20 µm wall): the exchange time an image "represents" depends on the pathology. Airway gas is in the spectrum but not in the ratio images.

**Norquay 2017 MRM.** RBC–plasma separation 20.4 → 25.5 ppm over sO2 0–1, exponential, field-independent; in vivo RBC shift falls ~1 ppm in 35 s (sO2 0.87 → 0.80) and oscillates 180° out of phase with the RBC amplitude. Quote: Dixon methods "rely upon a fixed chemical shift difference … drift in the 129Xe-RBC peak position would affect the intensity of the 129Xe-RBC phase image." Arithmetic: 1 ppm at 3 T over TE 0.45 ms ≈ 6°; the full sO2 range ≈ 28°.

**Norquay 2015 MRM.** T1 2.2 → 7.8 s over sO2; RBC ↔ plasma exchange τ ≈ 12 ms (ka 0.022, kb 0.062 ms⁻¹); plasma/RBC magnetization ratio 0.36 at Hct 0.48, 0.56 at Hct 0.39 — RBC:M is Hb-dependent by construction, and the compartments are exchange-coupled on the TR timescale.

**Xie 2019 NMR Biomed (Wuhan single-breath CSSR + DWI).** Global RBC fraction η 0.48 (young) → 0.33 (older healthy) → 0.28 (COPD); RBC/TP 0.47 → 0.30 → 0.24; SVRd 196 → 133 → 77 cm⁻¹. Age moves the "prior" as much as disease.

**Plummer 2023 MRM (Cincinnati, age-dependent reference, n = 30, 5–68 y).** Consortium Dixon at 3 T Philips with Willmering contamination removal; dissolved gridded at kernel sharpness 0.1 (gas 0.3) — "reduced the effective resolution but provided sufficient SNR". All four metrics significantly skewed; age explains 37 % of RBC:M median variance, 18 % of RBC median (peak ≈ 30 y), 9 % of membrane; simulated healthy medians vary "up to 50 %" over 5–80 y. Linear binning "assumes … a normal distribution … may not always be true"; RBC:M "more robust to differences in lung inflation".

**Plummer 2024 MRM (decay-modelled CS recon).** Forward model with RF/T1/T2* decay, wavelet + TV, ADMM; simulation RMSE 0.115 → 0.041; higher SNR and sharpness on gas-exchange data; median dissolved:gas and real/imaginary angle unchanged; 2× undersampling feasible. Gridding "ill-suited for undersampling"; "<20 % sampling". Chemical-shift phase along the 2.1 ms readout still not modelled.

**Willmering 2025 JMR (phantom).** See §4 row 15c.

**Teague 2024 Med Phys (UVA IDEAL, asthma).** Controls young/old RBC/gas 0.37 vs 0.30 %; same-day repeatability < 5 %; RBC/tissue "did not correlate with any other measures".

**Kern 2019 MRM (Hannover, localized CSSR vs IDEAL).** RBC/gas is the least reproducible metric under either method (CV 13–21 %); IDEAL "effectively extrapolated to an echo time of zero assuming a common T2*"; the gas field map is assumed to reflect dissolved off-resonance; lung inflation "one major source of variability".

**Marshall 2021 review (Sheffield).** Single-point method "suffers from contamination … though a technique for removal … has been reported"; gas images from multi-resonant imaging and dedicated ventilation images "are not interchangeable"; dissolved-phase standardization "remains work to be done".

**Ruppert 2016 MRM (CSSR pulsations, UVA).** RBC peak-to-peak pulsation 15.1 ± 3.8 % at ~1 Hz, 12 % at TLC vs 8 % at RV, up to 50 % variation between breath-holds; τ-independent → capillary volume, not velocity. A short calibration spectrum carries a cardiac-phase bias of that order unless whole cycles are averaged.

Agent's synthesis: critics converge on (1) single-TE physics (phase evolution across the readout; only k0 is clean), (2) TE90 and the RBC:M prior are not constants (sO2, cardiac phase, age, sex, Hct), (3) lung volume dominates variance and RBC:M is the only inflation-robust ratio, (4) age/sex references are mandatory, (5) compartment-specific T2* and frequencies belong in the decomposition. Best-evidenced: readout phase evolution (analytic) and the sO2/sex dependence of the RBC shift (two in-vitro calibrations at two fields + in vivo). Weakest objection in magnitude: unequal T2* (≤ 3–5 %). Still open: no correction for chemical-shift phase across the readout; single global angle and global RBC:M from a spectrum with a different cardiac/TR history; global T2*; third peak unresolvable at one TE; no published regional comparison of 1-point vs multi-echo in the same subjects.
