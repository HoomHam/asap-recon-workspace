# Respiratory binning for the s/d/p dynamic recons — decision draft (2026-10-26, DRAFT for Hooman, nothing applied)

Status: DRAFT. Coordinated by the ASAP session; audited by PCA (2026-pca-registration-04) and XeCS (2026-xecs-recon-e3).
Goal (Hooman): the binning that gives (1) representative binned data — one physiological state per bin, (2) balanced
interleave count per bin, (3) realistic signal amplitude vs time relative to the breath-by-breath dynamics; for s/d/p.

## 1. Candidates

| id | method | origin | one line |
|---|---|---|---|
| A | phase fraction | Steve SIGNAL (raw.py) | bin = (t − t_EE)/T_breath, EE = strict k0 min ±50 ilv; every breath a scaled copy |
| B | absolute time since EE | PCA | bin = (t − t_EE)/T_ref; tail → last bin (B) or dropped (Bx) |
| C | periodic line-constrained k-means | Ismail 2025 Med Phys 52:4318 §2.3 | (phase, amplitude) points of all cycles, 16-means with periodic phase, centroids on the mean breath curve; unbalanced counts, missing projections IMPUTED |
| D | amplitude rank, inhale/exhale split | Steve DIAPHRAGM (raw.bin) | equal counts by rank; time smeared; needs nav for d |
| E | time-line k-means | ASAP variant of C | x = abs time/T_ref instead of phase |
| F | **two-anchor** | XeCS | inspiration on an ABSOLUTE clock (bins 0–5, τ_I = 0.4 × median T), exhale + pause FRACTIONAL (bins 6–15) |

Kernels: `hard`, `soft1` = Steve's exp(−d²/2) cyclic (production, w_eff = 2.5 × count), `pou05` = σ 0.5 bin + partition
of unity (each interleave's weights sum to 1; XeCS recon/binning.py soft_membership).

## 2. Physics that decides it (XeCS, from `<key>_breaths_raw.csv`, 5 v3 sessions)

- t_B (bolus peak after EE) = 1.34 / 1.34 / 1.44 / 1.37 / 1.43 s, slope vs period ≈ 0 → fixed absolute time.
- t_C (nav max) absorbs ~25 % of the period variation; exhale + pause absorbs ~75 % (r 0.70–0.93). tC/T ≈ 0.35–0.43.
- So: phase bins (slope 1 everywhere) smear B; absolute time (slope 0) is right on inspiration, wrong past EI
  (long breaths pile the pause into the last bin, short ones never reach it). Cardiac-like: systole fixed, diastole stretches.
- Ismail's amplitude axis is a geometry surrogate (CT projection intensity ∝ lung volume); xenon k0 is bolus + decay +
  dose drift (−36…+24 %/min) → clustering on (phase, k0) sorts by dose, wrong axis. Its imputation of missing
  projections across respiratory states is not acceptable for gridding/CS. Its usable idea (constrain the representative
  curve, respect the wrap) is inside F.

## 3. Offline prototype (ASAP, `helpers/_binning_eval.py`, `outputs/binning_eval_2026-10-26/`)

Exact Tyger EE + auto-exclusions (§12 of BinTiming_Facts), 5 sessions, k0 (s) and nav volume (d, nav's own EE).
Truth = median over breaths of amplitude at absolute time after EE; bin time = kernel-weighted, EE wrap folded per bin.
Median over sessions, k0 signal:

| method/kernel | count CV | min/max count | amp RMS /(B−A) | B/A recovered | t_B error s | within-bin amp sd | within-bin time sd s |
|---|---|---|---|---|---|---|---|
| A phase soft1 (production) | 0.022 | 0.95 | 0.075 | 0.81 | +0.20 | 0.168 | 0.40 |
| A phase pou05 | 0.024 | 0.95 | 0.070 | 0.91 | +0.20 | 0.131 | 0.31 |
| B abstime pou05 | 0.066 | 0.77 | 0.035 | 0.93 | +0.16 | 0.115 | 0.22 |
| B abstime hard | 0.234 | 0.29 | 0.027 | 0.93 | +0.06 | 0.103 | 0.16 |
| C Ismail hard | 0.297 | 0.38 | 0.058 | 0.97 | +0.18 | 0.060 | 0.29 |
| D amp-rank on k0 | 0.25–0.32 | 0.00–0.31 | 0.14 | — | — | 0.04 | 0.5–0.7 |
| **F two-anchor pou05** | **0.043** | **0.87** | **0.034** | **0.92** | +0.16 | 0.121 | 0.23 |
| F two-anchor hard | 0.054 | 0.84 | 0.032 | 0.93 | +0.06 | 0.106 | 0.19 |

Per session, B/A recovered (bin / per-breath): F pou05 015HK 4.20/4.58, 014KS 3.94/4.68, 039CP 4.55/4.92,
016PG 6.56/6.51, 025JC 2.65/2.91 — vs A soft1 3.85, 3.24, 3.97, 4.70, 2.42 on the same denominators.
On the navigator (d): F pou05 amp RMS 0.025 vs A soft1 0.036, D amp-rank 0.052 (D over-recovers B/A 1.44: rank bins
select the deepest breaths' extremes, a state but not the typical breath); counts F 0.94 min/max.

Reading: F is best on every criterion jointly; B equals it on fidelity but loses count balance; A (production) loses
~20 % of B/A and has the widest within-bin time spread; C/D are amplitude-homogeneous but unbalanced or time-blind.
Residual F deficit (~8 % B/A) = kernel σ 0.5 on 0.3 s inspiratory bins + bin 0 straddling EE; hard bins give 0.93.

## 4. Proposed approach (F, to be ratified)

1. EE anchors: Steve's strict k0 minima (±N = 50 ilv, header `sigminN`), unchanged. Recon-free, works for s and p;
   for d the same k0 anchors are used (nav max is not needed as an anchor).
2. τ_I = split × median period, split = 0.4 (session-median tC/T 0.35–0.43; sensitivity table in §5). Nav-free.
3. φ = 0.4 (t − t_EE)/τ_I on 0 ≤ t − t_EE < τ_I (6 bins × ≈0.3 s absolute: A, rise, B ≈ bin 4, C ≈ bin 6);
   φ = 0.4 + 0.6 (t − t_EE − τ_I)/(T_k − τ_I) after (10 fractional bins, equal count per breath by construction).
4. Exclusions: breaths with T_k < 1.5 τ_I; auto low-signal sets as today; interleaves outside first/last EE → −1.
5. Kernel: σ = 0.5 bin + partition of unity (needs a fork change: per-interleave normalisation of binwt; bindist0sq → 0.5
   and a precomputed Σ_b w_b,i). Keep bin centres at b/16 (bin 0 = EE, Steve's convention).
6. d variant (F′, XeCS): exhale limb binned by nav-amplitude fraction instead of time fraction; test only, expect
   agreement within kernel noise. p: same φ on the pneumotach EE clock once the alignment is trusted (not now).
7. Export with the recon (bin-export patch): ilvbin, ee_times, τ_I, per-bin w_eff, mean time, nav z.

Fork footprint: raw.py SIGNAL block (≈15 lines: φ instead of tee/T), recon.py kernel width + normalisation (≈10 lines),
converter flags `binmode=twoanchor`, `split`, `binsigma`. Old behaviour stays the default until validated.

## 5. Validation and tests before it enters the fork (XeCS T1–T8, adopted; ASAP additions T9–T10)

Offline first (ilv tables), then ONE Tyger run on 015HK + 016PG, then the 5.
- T1 kernel-weighted k0 B/A per bin vs per-breath smoothed median within 5 % (now −19 % A, −28 % 016PG); same for C/B, nextEE/C.
- T2 within-bin time sd ≤ 0.3 s on bins 0–6 in all 5; nav-z sd on exhale bins ≤ diaphragm scheme's.
- T3 count balance: max/min hard count ≤ 1.5 over 16 bins (B fails on 016PG; report A/B/D/F side by side).
- T4 Σ_FOV gas_phase(bin)/k0_wmean(bin) flat ±2 % with exact membership (§12 test) — must hold with partition-of-unity weights.
- T5 stress cases 016PG (irregular) and 014KS (dead window): a bin populated by < 60 % of breaths is flagged, never imputed.
- T6 image: diaphragm 10–90 % rise distance per bin (Ismail's metric) on d/ vs s/ vs F; bin-6 (EI) image at least as sharp as diaphragm bin 7/8.
- T7 downstream: PCA's per-bin global budget (m_A, ρJ) from F-binned images vs the per-breath k0 budget agrees within the hold-vs-tail decay band.
- T8 no imputation, no borrowed interleaves; excluded = zero weight.
- T9 split sensitivity (0.35 / 0.40 / 0.45 / session tC/T) changes B/A recovered by < 3 % and count CV by < 0.02.
- T10 regression: A (production) reproduced bit-exactly when `binmode` is absent (default path untouched).

## 6. Scores (1–5: representativeness · count balance · amplitude-vs-time · irregular breathers · s/d/p applicability)

| | A | B | C | D | F |
|---|---|---|---|---|---|
| XeCS | 3·5·2·2·5 | 3·2·3·2·5 | 3·2·3·3·2 | 3·5·2·3·3 | 5·4·4·4·4 |
| PCA (geometry-weighted audit, before F existed) | 2·5·2·2·4 | 3·3·4·2·4 | 4·2·3·3·3 (Cn) | 4·5·3(4 at σ .5)·4·3 | — |
| ASAP (from §3 numbers) | 3·5·2·2·5 | 4·3·4·2·5 | 3·2·4·3·2 | 3·5·2·3·3 | 4·5·4·4·4 |

## 7. Open points for Hooman

- Ratify F (two-anchor) as the candidate; split fixed 0.4 vs per-session tC/T; σ 0.5 + partition of unity vs Steve's σ 1.
- Whether the d/ recon moves to F too (one binning for s/d) or keeps the diaphragm rank bins as a geometry reference.
- Go for the fork patch + one Tyger validation run (015HK, 016PG).

## 8. PCA audit (received after §3; `2026_PCA_registration/.../outputs/binning_2026-10-26/scheme_metrics.csv`)

Independent emulator (reproduces ilvbin_tyger to 0.000 bins), metric = geometry blur (kernel-weighted within-bin SD of
diaphragm z, vox) + k0 B/A recovered + hard-count CV. Findings: (1) D (nav rank) has the lowest geometry blur at every
kernel; A has the worst; (2) **the kernel is the biggest single lever**: σ 1 → 0.5 → hard moves D's B/A recovered
0.90 → 1.00 → 1.03 and blur 1.72 → 0.94 → 0.60 vox (015HK), B 0.91 → 1.00 → 1.02; cost w_eff/hard 2.5 → 1.25 → 1
(noise × √2 at σ 0.5); (3) Ismail with nav amplitude (Cn) = D's geometry with worse balance — no gain; with k0 in the
feature it bins on the measured quantity — reject; (4) 014KS D count CV 0.20 = the interpolated dead window (fix:
NaN, not interp). PCA proposal: D for s/d/p + σ 0.5 + dead-window NaN fix; per-breath k0 (unbinned) for content/time;
an absolute-time signal stack only as an optional display. Note: PCA wrote "bindist0sq 2 → 8" for σ 0.5; the kernel is
exp(−d²/bindist0sq) so σ² = bindist0sq/2 → **σ 0.5 ⇔ bindist0sq 0.5**, not 8 (8 would be σ 2).

## 9. Geometry blur added to the ASAP prototype (`outputs/binning_eval_2026-10-26/geometry_blur.txt`)

Mean within-bin SD of nav z (vox) | k0 B/A recovered, per session, kernel soft1 / pou05 / hard:

| | 015HK | 014KS | 039CP | 016PG | 025JC |
|---|---|---|---|---|---|
| A phase | 2.00/1.41/1.23 · .84/.91/.93 | 2.02/1.66/1.54 · .69/.76/.79 | 1.48/1.11/1.00 · .81/.91/.93 | 2.91/2.59/2.50 · .72/.89/.91 | 1.27/1.01/0.93 · .83/.91/.93 |
| F two-anchor | 1.95/1.34/1.14 · .85/.92/.93 | 1.91/1.49/1.34 · .77/.84/.87 | 1.49/1.12/1.00 · .82/.92/.95 | 2.91/2.58/2.47 · .84/1.01/1.03 | 1.27/1.02/0.94 · .85/.91/.93 |
| D nav rank | 1.73/0.94/0.61 · .86/.96/.99 | 1.35/0.78/0.53 · .81/.90/.93 | 1.27/0.71/0.45 · .84/.96/1.00 | 2.05/1.22/0.81 · .75/.95/1.19 | 0.99/0.57/0.39 · .87/.94/.98 |

Reading: on geometry D is 2× better than any time axis at the same kernel; F ≈ A ≈ B on geometry (the time axes are
equivalent there); on signal-vs-time F ≥ B > A (§3), and D is time-blind on the plateau bins (§13 of BinTiming).
τ_I fixed at 2.0 s cohort-wide (F2) breaks count balance on short-period subjects (015HK CV 0.24, 025JC 0.31) →
keep τ_I = 0.4 × median period (`split_sensitivity.txt`: B/A within 1 % for split 0.30–0.45).

## 10. Consolidated recommendation (ASAP, for Hooman to ratify)

Two purposes, two axes, one kernel change:
1. **Kernel for all binnings: σ 0.5 bin + partition of unity** (`bindist0sq` 2 → 0.5 in recon.py plus per-interleave
   normalisation Σ_b w = 1). Biggest gain on every criterion for every scheme; noise × √2 (SNR budget: gas fine,
   dissolved to be checked per session). Both auditors converge on this.
2. **d/ = D (Steve's nav-rank, inhale/exhale split) stays the geometry stack** (registration, PCA's lane), with the
   dead-window fix (nav_volume → NaN inside auto-excluded sets, so bins are −1 like the kernel treats them).
3. **s/ = F two-anchor becomes the signal-time stack** (Step-1 A→B→C dynamics): best amplitude-vs-time among
   balanced schemes, B/A recovered 0.92 (pou05) vs 0.81 production, counts CV ≤ 0.09, nav-free. Per-breath k0 stays
   the content reference (unbinned). Steve's phase (A) retired for s/ once T1–T10 pass; stays reachable as default.
4. p/: untouched until the pneumotach alignment is validated.
5. Ismail 2025: not adopted; cited as the single-signal CT case; its constraint idea is inside F, its imputation is
   excluded (T8).
Fork footprint: raw.py SIGNAL block (φ two-anchor, ≈15 lines), recon.py kernel (≈10 lines), navigator NaN rule
(2 lines), converter flags `binmode`, `split`, `binsigma`; bin-export patch for validation. Validation run: 015HK +
016PG first (s with F, d with D, both σ 0.5), then the 5.

## 11. XeCS audit of the prototype + final arms (`outputs/binning_eval_2026-10-26/final_arms.txt`)

XeCS code audit of `bins_twoanchor` + `pou05`: PASS. Two evaluation fixes applied: (1) the EE-wrap fold now uses the
bin's nominal coordinate (data-dependent median flipped on 016PG: bin 0 reported +9 s); (2) truth B/A = median of the
Savgol-smoothed per-breath picks (`B_over_A_smsm`, 4.34/4.41/4.78/5.17/2.76) instead of extreme-value picks (3 % high).
Added arms: per-breath exhale weight normalisation w ← w·clip((T_med−τ_I)/(T_k−τ_I), 0.5, 2) on members with
t−t_EE ≥ τ_I (`_pbn`); τ_I fixed at 1.8 s cohort-wide (`F_1.8s`).
Median over 5 sessions, k0 (count CV | min/max | amp RMS | B/A ÷ B/A_smsm | t_B err s | within-bin time sd s):

| arm | hard | pou05 (σ 0.5, partition of unity) |
|---|---|---|
| A phase | 0.025 · 0.94 · 0.065 · 0.979 · +0.06 · 0.29 | 0.024 · 0.95 · 0.070 · 0.956 · +0.20 · 0.31 |
| B abstime | 0.234 · 0.29 · 0.027 · 0.981 · +0.06 · 0.10 | 0.066 · 0.77 · 0.031 · 0.981 · +0.16 · 0.18 |
| F two-anchor 0.4 | 0.054 · 0.84 · 0.032 · 0.978 · +0.06 · 0.19 | 0.043 · 0.87 · 0.034 · 0.960 · +0.16 · 0.23 |
| **F 0.4 + pbn** | **0.029 · 0.92 · 0.029 · 0.974 · +0.06 · 0.19** | 0.024 · 0.94 · 0.031 · 0.959 · +0.16 · 0.24 |
| F 1.8 s | 0.146 · 0.70 · 0.028 · 0.975 | 0.143 · 0.73 · 0.032 · 0.951 — count balance lost on short periods → dropped |

Per session (hard, B/A ÷ smsm): A 0.98/0.83/0.96/1.15/0.98; F 0.99/0.93/0.98/1.30/0.98; F+pbn 0.98/0.92/0.96/**1.04**/0.97
(015HK/014KS/039CP/016PG/025JC). Per-breath normalisation is what fixes 016PG (long pauses decayed below the median
A pooled into bin 15 → B/A over-shoot 1.30 → 1.04). 014KS 0.92 is the scan's own t_B drift (sd 0.46 s) + dose rise,
not a binning degree of freedom (XeCS). Production A at σ 1 (soft1) sits at 0.81–0.86 on the raw truth (§3).

## 12. E′ — one axis for s/ and d/ (`eprime_geometry.txt`, `eprime_pbnorm.txt`)

XeCS objection to §10: s/ on F and d/ on D puts Step-1 numbers on different interleave sets (a d/ bin 7 ≠ an s/ bin 7),
and D's EE bins pool 1–2 s of real time where the dissolved state is time-driven. PCA: F's geometry blur equals A's
(known property, record it). Both asked for E′ = two-anchor with the exhale limb on navigator-amplitude RANK
(within-exhale rank of nav_volume pooled over breaths, equal count; inspiration unchanged on the absolute clock).

Geometry blur = within-bin SD of nav z (vox), split by limb (insp bins 0–5 | exhale 6–15); medians over 5 sessions, hard / pou05:

| arm | navz sd insp | navz sd exh | B/A ÷ smsm | count CV | time sd insp | time sd exh |
|---|---|---|---|---|---|---|
| A phase | 1.32 / 1.58 | 1.18 / 1.32 | 0.979 / 0.956 | 0.025 | 0.31 / 0.33 | 0.35 / 0.37 |
| F two-anchor | 1.09 / 1.27 | 1.12 / 1.26 | 0.978 / 0.960 | 0.054 | 0.09 / 0.15 | 0.25 / 0.28 |
| F + pbn | 1.09 / 1.26 | 1.11 / 1.25 | 0.974 / 0.959 | 0.029 | 0.09 / 0.15 | 0.25 / 0.29 |
| **E′** | 1.09 / 1.28 | **0.57 / 0.76** | **1.004 / 0.977** | **0.024 / 0.014** | 0.12 / 0.18 | 0.36 / 0.38 |
| E′ + pbn | 1.10 / 1.28 | 0.57 / 0.74 | 0.998 / 0.971 | 0.052 / 0.043 | 0.12 / 0.18 | 0.36 / 0.38 |
| D nav rank | 0.56 / 0.84 | 0.52 / 0.78 | 1.034 / 0.991 | 0.001 | **1.00 / 1.05** | 0.44 / 0.48 |

Reading:
- E′ matches D on exhale geometry (0.57 vs 0.52 vox) and F on inspiratory timing (0.12 s), with the best k0
  fidelity (1.00) and balanced counts (CV 0.02). The remaining geometry gap to D is on the inspiratory limb only
  (1.09 vs 0.56 vox): at a fixed time after EE the diaphragm sits at different heights in different breaths
  (slope 0.25). That is the price of a time axis where B lives; PCA asks that it be recorded as a property.
- D is time-blind on BOTH limbs: inspiratory within-bin time sd 1.0 s (rank bins mix early and late inhale across
  breaths), exhale 0.44 s. For dissolved A/B states this is the wrong axis; its geometry advantage on inspiration is
  bought with 1 s of time mixing.
- 016PG under E′: B/A over-shoots (1.44 hard; 1.25 with pbn) for the same reason as F without pbn (lowest-volume
  exhale rank = most-decayed pause ends). pbn helps less here than for F (1.04) because rank already sorts those
  into bin 15 regardless of weight. Known limit on 4–12 s breathers: flag the bin, never impute (T5).
- Bin 0 under Steve's centre convention straddles the wrap: E′ bin 0 = first 0.3 s of inspiration + lowest 3 % of
  exhale volume (both ≈ A; t = −0.1 s). XeCS tweak 6(b) (centres at (b+½)/16) would separate them; test, not blocker.

## 13. Frozen rule (for PCA's emulator, test (a)) — ASAP prototype as coded in `helpers/_binning_eval.py`

- EE list = Tyger's (`<key>_ee_tyger.csv`); breath k: t_EE,k ≤ t < t_EE,k+1; tee = t − t_EE,k; T_k = t_EE,k+1 − t_EE,k;
  outside [first, last EE) → −1; Tyger auto-excluded interleaves → −1 (after binning).
- τ_I = 0.4 × median(diff(EE)) per session (all gaps, global, NOT per breath). Breaths with T_k < 1.5 τ_I → −1.
- F: φ = 0.4·tee/τ_I for tee < τ_I; else 0.4 + 0.6·(tee − τ_I)/(T_k − τ_I) (equal FRACTION of each breath's
  exhale, not rank). φ clipped to [0, 1). Inspiration therefore spans 6.4 bins (0.4 × 16), exhale 9.6 — PCA's 6/10
  emulation is a different rule (slightly). Pause handled by the fraction (it is the end of the exhale limb).
- E′: inspiration as F; exhale φ = 0.4 + 0.6·(r + ½)/N, r = descending rank of nav_volume among ALL kept exhale
  members (tee ≥ τ_I, finite nav) pooled over breaths; NaN nav → −1.
- Bin coordinate c = 16 φ; hard bin = floor(c + ½) mod 16 (Steve's centres at b/16, bin 0 straddles the wrap);
  kernel w = exp(−d²/(2σ²)), d cyclic in bins, σ = 1 (today, bindist0sq 2) or 0.5 (bindist0sq 0.5); `pou` divides
  each interleave's weights by their sum.
- pbn: w_i ← w_i · clip((T_med − τ_I)/(T_k − τ_I), 0.5, 2) for members with tee ≥ τ_I only.
- Truth for T1: median over breaths of `B_over_A_smsm` (breaths with T > 1.5 × median dropped).

## 14. Revised recommendation (supersedes §10 points 2–3)

1. Kernel σ 0.5 + partition of unity for all binnings (both auditors; noise ×1.41 vs today, PCA measured
   0.531/0.751/1.000 relative bin-mean noise for σ 1 / 0.5 / hard).
2. **One axis for s/ and d/: E′** (inspiration absolute clock, exhale nav rank) so Step-1 per-bin numbers from d/ masks
   and s/ signal are the same interleave sets. Needs the navigator in every recon (one low-res pass, ~1 min GPU;
   `export_nav` plumbing from the bin-export patch). Fallback if E′ fails T1–T6 on the Tyger run: F (+pbn) for s/ and
   D* (dead-window NaN fix, PCA: 014KS count CV 0.20 → 0.017, free) for d/, with T11 (never ratio a d/ bin to an
   s/ bin; recon type in every column name) and T12 (D's EE bins reported with their real-time sd next to every
   dissolved number).
3. The F/E′ stack is the signal-time stack; registration / Jacobian input stays the geometry stack (D* or E′ once
   T6 + PCA's cascade bench say so). Record: time-axis bins have A-level geometry blur on inspiration (one time, mixed
   diaphragm heights across breaths).
4. Per-breath exhale normalisation: in for F; for E′ optional (helps 016PG less); decide on the Tyger run.
5. p/ untouched; Ismail cited, not adopted (one-signal CT case; constraint idea inside E′, imputation excluded).
6. Order of work after Hooman's go: fork patch (raw.py φ rule + `binmode`/`split`/`binsigma`/`pbn` flags, recon.py
   kernel width + per-interleave normalisation, navigator NaN rule + `export_nav`, bin export) → test (a) offline
   match to 0.000 bins → Tyger 015HK + 016PG (s and d with E′, σ 0.5) → T1–T12 → the 5 → PCA cascade bench.

## 15. Amendments to the frozen rule (XeCS endorsement of E′, 2026-10-26 late; adopted into §13/§14 pending Hooman)

- **split = 0.375 (= 6/16)**, not 0.4: the EE→EI junction then falls exactly on a bin boundary (bin 6 centre with
  centres at b/16), so no bin is half absolute clock / half nav rank. Inside the measured insensitivity (< 1 % on B/A
  for 0.30–0.45). τ_I = 0.375 × median period.
- **A for Step-1 tables = bin 0** (first ≈ 0.3 s after EE on the absolute clock) or the per-breath k0_A — never
  min-over-bins. E′'s exhale rank sorts the most-decayed pause ends into bin 15; min-over-bins picks a sub-median A
  and over-shoots B/A (016PG 1.44 / 1.25). With A = bin 0 the over-shoot is gone by construction; bin 15 keeps its
  honest label "end of exhale, geometry-sorted, decay-mixed". Same rule for D* in the fallback.
- Exhale nav rank uses nav_volume AFTER the dead-window NaN fix (never the interpolated values).
- The inspiratory geometry gap (E′ 1.09 vs D 0.56 vox) is the breath-to-breath tidal-volume spread at a fixed time
  after EE, same order as the diaphragm's own travel inside a 0.3 s bin (~20 mm/s × 0.3 s ≈ 6 mm ≈ 1.7 vox); a
  time-homogeneous bin is not meaningfully blurrier than its own duration. Recorded as a property.
- Status: XeCS endorses E′ as the one axis for s/ and d/; PCA endorses with its three conditions (§12) and is
  emulating E′ for test (a). Over to Hooman.

## 16. PCA independent emulation of the frozen rule (test (a), offline half) — PASS

PCA `_bin_schemes.py` re-implemented F and E′ from §13 alone. Medians over 5, mine vs PCA (hard/pou05): E′ count CV
0.024 vs 0.026; E′ navz sd insp 1.09/1.28 vs 1.02/1.24, exh 0.57/0.76 vs 0.59/0.79; E′ B/A 1.004/0.977 vs
1.009/0.967; F count CV per session 0.025/0.096/0.071/0.054/0.020 vs 0.024/0.094/0.060/0.042/0.032. Only definitional
mismatch: inspiratory time sd (PCA 0.27 s includes the wrap-straddling bin 0). Rule is reproducible; the 0.000-bin
match is now the export's job. PCA position, adopted: (1) E′ = the ONE axis for Step-1 accounting (s/ and d/ same
interleave sets); (2) registration geometry stack stays **D\* until PCA's cascade bench** (markers / Jacobian / leak,
015HK, E′ vs D* at σ 0.5) says otherwise — E′ inspiratory bins carry 1.5–2× D*'s diaphragm spread, which is the motion
registration must resolve; §14.3 stays "D* or E′ once the bench says so"; (3) 016PG stays a flagged case under every
time-axis scheme (B/A 1.56 pou05 in PCA's emulation with min-over-bins; gone with A = bin 0, §15); (4) (b+½)/16 bin
centres tested on the Tyger run, not a blocker. Draft ready for Hooman from all three sessions.

## 17. Hooman's decision + the ≥ 50-session statistics (2026-10-26 night)

Hooman (via PCA): "to decide on F, I need many many more subjects added to the stat first. If the stats stay the same
we go for it, at least 50 datasets need to show us the same." Addendum: "the bolus injection can be slightly different
per study." → No fork change. Export extended to every old-tree session (108 with input.mrd + navigator-bearing
d/recon.mat; CPU 10 s each; no GPU needed, the navigator exists for all). `helpers/_session_index.py` →
`outputs/bin_timing_2026-10-26/session_index.csv`, one row per session: cohort (cohort_map), protocol (TR), seq_roster,
EE count, exclusions (incl. mid-scan), period stats, t_B / t_C medians, slopes vs period, exhale share, τ_I per session
and the fraction of peaks inside it, injection profile (rise t50/t90, early-arrival fraction at 0.25/0.5 s), B/A, C/B,
k0 trend, `layout_ok`.
Excluded from the pool (`layout_ok` = 0): 9 mixed-ADC sessions (gas/diss lines not alternating in input.mrd; raw-k0
export invalid as written; `_EXCLUDED_mixed_adc.txt`), 7 TR 22.3 ms sessions + 000LL 2024-02-02 (TR 25.9, 8-ch):
different sequences, dissolved/gas at EE ≈ 1 → line slicing wrong. Pool = 91 sessions (v3 62, v2 29).

Medians by protocol × cohort (t_B s | t_B/T | slope t_B vs T | t_C/T | exhale slope | frac t_B < 0.375 T | t50 s |
early frac 0.5 s | B/A | n irregular | n mid-scan exclusion):

| | n | t_B | t_B/T | slope | t_C/T | exh slope | B inside τ_I | t50 | early | B/A | irreg | mid |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| v3 HC | 4 | 1.38 | 0.30 | +0.11 | 0.39 | +0.74 | 0.98 | 0.77 | 0.17 | 4.39 | 3 | 1 |
| v3 EBV | 18 | 1.26 | 0.32 | +0.18 | 0.34 | +0.74 | 0.90 | 0.81 | 0.15 | 3.27 | 7 | 0 |
| v3 LTX | 31 | 1.39 | 0.35 | +0.17 | 0.38 | +0.73 | 0.79 | 0.83 | 0.14 | 3.62 | 13 | 1 |
| v3 RT | 9 | 1.22 | 0.38 | +0.28 | 0.41 | +0.78 | 0.43 | 0.80 | 0.15 | 2.63 | 6 | 0 |
| v2 LTX | 16 | 1.35 | 0.33 | +0.18 | 0.40 | +0.63 | 0.87 | 0.81 | 0.16 | 3.73 | 4 | 5 |
| v2 EBV | 4 | 1.32 | 0.23 | +0.06 | 0.33 | +0.83 | 1.00 | 0.83 | 0.12 | 3.77 | 3 | 2 |
| v2 RT | 5 | 1.28 | 0.33 | +0.21 | 0.37 | +0.64 | 0.96 | 0.71 | 0.24 | 3.53 | 4 | 2 |
| v2 HC / HC_OLD | 2 + 2 | 1.31 / 1.45 | 0.30 / 0.29 | +0.25 / +0.31 | 0.42 / 0.34 | +0.39 | 0.67 / 0.99 | 0.74 / 1.02 | 0.21 / 0.03 | 5.97 / 5.46 | 0 | 0 |

What holds on 62 v3 sessions (the 5-session physics): t_B absolute 1.2–1.4 s in every cohort; the exhale + pause
absorbs 73–78 % of the period variation (v2 63–83 %); injection profile uniform (t50 0.7–0.9 s, early fraction
0.12–0.24) except the HC_OLD pair (0.03) and the non-v2/v3 sequences — so delivery looks consistent within v2/v3 at the
k0 level; per-session columns are there for Hooman's clustering by date / setup.
What does NOT hold: the t_B slope vs period is +0.17–0.28 in patients (HC +0.11), not ≈ 0, and the fraction of peaks
inside a 0.375·T inspiratory limb is 0.79 (LTX) and 0.43 (RT: short periods, t_B/T 0.38). **τ_I must be per session
from the session's own t_B / t_C** (e.g. τ_I = max(0.375·T_med, 1.25·t_B,med), or τ_I = t_C,med with the nav), never a
cohort constant — exactly Hooman's addendum. The two-anchor structure itself (fixed inspiratory clock, stretched
exhale) is confirmed at n = 62; the split parameter is not.
Irregular-breather fraction (period IQR/median > 0.3 or any breath > 1.5× median): v3 29/62; mid-scan exclusions 2/62
(v2 9/29). Scheme emulation A/B/F/D*/E′ over the 91 = PCA's `_bin_schemes.py` (joins on key), result pending.
Pooled v3 (62, layout_ok) medians [IQR]: t_B 1.35 s [1.16–1.50], t_B/T 0.345 [0.30–0.37], slope t_B vs T +0.18
[+0.07..+0.28], t_C/T 0.375 [0.34–0.40], exhale slope +0.74 [0.62–0.81], T 4.0 s [3.3–4.7], period CV 0.11, t50
0.82 s [0.72–0.88], early fraction 0.14 [0.11–0.21], B/A 3.3 [2.6–4.2], k0 trend −10 %/min [−33..0], excluded
interleaves 12 % [7–16]. Only **19 of 62** sessions keep ≥ 95 % of their peaks inside a 0.375·T inspiratory limb →
the split parameter is per session or the scheme fails on half the cohort.

## 18. Scheme emulation on the 91 (PCA `_bin_schemes.py`, `2026_PCA_registration/.../outputs/binning_2026-10-26/cohort_table.md`)

v3 medians by cohort (EBV 18 / HC 4 / LTX 31 / RT 9); count CV · B/A recovered · navz sd insp · navz sd exh (vox) · insp time sd (s):
- TODAY (σ 1): A phase .022–.041 · 0.83–0.88 · 1.15–2.33 · 1.08–1.83 · 0.66–0.97; D* .022–.035 · 0.80–0.90 · 0.84–1.79 ·
  0.81–1.41 · 1.0–1.3; E′ .024–.031 · 0.79–0.89 · 1.06–2.03 · 0.78–1.41 · 0.7–1.0.
- HARD: A .022–.041 · 0.93–0.98 · 0.88–1.53 · 0.89–1.30 · 0.40–0.70; D* .022–.035 · 0.90–1.03 · 0.41–0.68 · 0.38–0.50 ·
  0.86–1.16; F(0.4T) .026–.051 · 0.92–0.97 · 0.76–1.12 · 0.89–1.34 · 0.21–0.26; E′(0.4T) .024–.031 · 0.93–0.995 ·
  0.74–1.08 · 0.41–0.60 · 0.21–0.25.
- Per-session τ_I = max(0.375·T_med, 1.25·t_B,med) with the bin split still fixed 6.4/9.6: fidelity unchanged
  (0.91–1.00) but balance breaks on fast breathers (RT CV 0.17) — τ_I/T grows while inspiration keeps 6.4 bins.
- **E′ad (PCA fix): per-session τ_I AND n_insp = round(16·τ_I/T_med) clamped 4–10, exhale gets the rest**:
  CV .032/.052/.033/.045 · B/A 0.92/0.99/0.97/0.98 · navz insp 0.93/1.07/0.99/0.74 · exh 0.47/0.58/0.63/0.45 ·
  insp time sd 0.21–0.25 (EBV/HC/LTX/RT). Balance of D*, inspiratory timing of F, exhale geometry of D*, on all four
  v3 cohorts; same picture on v2 (29). pou05: B/A 0.88–0.97, geometry +0.1–0.2 vox, as on the 5.

Readings: (1) ranking on 62 v3 + 29 v2 = the 5-session audit — Hooman's "stats stay the same" condition is met for the
STRUCTURE; (2) kernel σ is still the first lever (A at σ 1 0.83–0.88 → hard 0.93–0.98); (3) the split must be per
session AND the bin allocation must follow it (E′ad), otherwise balance fails on RT; (4) D* remains the registration
geometry stack (time-axis inspiratory blur 0.74–1.1 vox vs D* 0.41–0.68 on every cohort). Two sessions failing in
PCA's loader (042DR, 001JM) = the mixed-ADC list, consistent.

## 19. Frozen rule, v2 (supersedes §13/§15 parameters; pending Hooman)

E′ad: EE = Steve's k0 strict minima (sigminN); τ_I = max(0.375·T_med, 1.25·t_B,med) per session (nav-free; t_B from
the per-breath k0 peak), n_insp = round(16·τ_I/T_med) clamped [4, 10]; inspiration bins 0..n_insp−1 on the absolute
clock (φ = (n_insp/16)·(t−t_EE)/τ_I), exhale bins n_insp..15 on nav-volume rank pooled over kept breaths (NaN nav
after the dead-window fix → −1); breaths with T < 1.5·τ_I → −1; kernel σ 0.5 + partition of unity; A = bin 0; bin
centres b/16 (Steve), (b+½)/16 to be tested on the Tyger run. Export n_insp and τ_I with the recon. Registration
stack stays D* until PCA's cascade bench. Everything else as §14/§15 (T1–T12, no imputation, p/ untouched).

## 20. Mix-binning applied (Hooman's vote 10-27) + the pneumotach clock finding

- Hooman: against adopting E′ad (navigator fails on low-SNR sessions; bolus-based criteria not future-proof with
  proportional delivery; keep studying s/d/p) → MIX applied ONCE to 015HK / 016PG / 014KS as a fourth stack
  ("mix-binning"), four-way comparison s/d/p/mix, nothing adopted. Fork commits (diaphragm-recon → hooman):
  508b329 MIX mode + binsigma/binpou knobs + bin export; efde661 converter accepts MIX; 36f67fb int-params fix
  (508b329 broke every binned recon: float params → float kernel indices; caught by the first run = T10).
  Output tree `Ext/AIkill_Dynamic_mix/<key>/{m,p}`; `pipeline/run_mix3.sh`; verification `helpers/_mix_check.py`.
- **Pneumotach clock (PCA, `pneumo_clock_fit.csv`)**: matching pneumotach-volume minima to the Tyger k0 EE list gives
  t_k0 = a·t_pn + b with a ≈ 0.989–0.991 and b ≈ −26…−32 s on 015HK/014KS/016PG/039CP (residual 0.1 s, 23–29 of 24–29
  EE matched); at the coded zero shift only 0–9 breaths match. Rank bins on the pneumotach volume give k0 B/A 1.4–2.2
  as coded vs 3.9–5.1 re-aligned (diaphragm rank on the same interleaves 3.7–4.9). So the pneumotach SIGNAL is fine;
  the first-packet alignment (converter zeroes the pneumo clock at packet 0 ≡ gas interleave 0, raw.py:184-202) misses
  by ~30 s and the pneumo clock runs ~1 % fast vs iilv·TR. Every old-tree p/ recon (93) is mis-binned; a p rerun
  as coded reproduces it. No Steve knob exists (gvar.bindt is read from the header but unused).
  Proposed fork fix (Hooman's call): per-session linear alignment of the pneumotach clock to the gas lines by the
  minima-to-minima fit (or cross-correlation of volume vs k0), applied in raw.load_pneumotach_from_arrays before
  cropping, with pneumo_clock_a / pneumo_clock_b exported; then p = a second amplitude-rank stack (volume instead
  of diaphragm height) and the four-way is meaningful. p reruns HELD until decided.

## 21. MIX stacks delivered (2026-10-27 03:04–03:12 Tyger, image 36f67fb) — `outputs/binning_eval_2026-10-26/mix_check.txt`

| | 015HK | 016PG | 014KS |
|---|---|---|---|
| τ_I / n_insp | 1.68 s / 7 | 2.86 s / 6 | 1.92 s / 6 |
| binned / total ilv | 5569 / 7017 | 10654 / 11664 | 7124 / 8945 |
| hard count CV (min/max) | 0.028 (0.94) | 0.059 (0.85) | 0.032 (0.91) |
| exported ilvbin vs offline rule | 0.0000 | 0.0000 | 0.0000 |
| Σ_FOV gas / k0_wmean per bin | 0.986–1.017 | 0.983–1.017 | 0.972–1.022 |
| median t after EE, insp bins 1→6 | 0.23→1.44 s | 0.47→2.98 s | 0.32→1.94 s |
| nav z sd, insp / exhale bins (vox) | 0.9–1.5 / 0.3–0.8 | 1.9–3.0 / 0.4–1.0 | 0.9–1.4 / 0.3–0.8 |

Tree `Ext/AIkill_Dynamic_mix/<key>/m/` (recon.mat with bin-export keys + nav + calcb_b + rbc_tp_*). s = `AIkill_Dynamic_b44d_sig/<key>/s`,
d = `AIkill_Dynamic_b44d/<key>/d` on the same image lineage; p held (§20). Four-way image comparison = PCA's loader.
Not adopted; one-off per Hooman's vote.

## 22. Four-way image read (PCA, `2026_PCA_registration/.../outputs/step1_global_2026-10-26/{fourway.png, fourway_bins.csv, summary_mix.csv}`)

EI-mask real-sum gas B/A, d / s / mix vs per-breath k0: 015HK 3.95 / 3.86 / 4.40 vs 4.33; 014KS 3.86 / 2.96 / 4.06 vs
4.41; 016PG 3.68 / 3.06 / 6.04 vs 5.20 (mix over-shoot on the irregular breather: dose-weighted exhale bins, as the
emulator predicted). Mix t_B 1.44 / 1.61 / 1.91 s (bins 6 / 5 / 4 of n_insp 7 / 6 / 6), C at bins 7 / 6 / 6; C/B 0.95 /
0.96 / 0.78. Dissolved: D/G shape binning-independent — D/G C/A 0.48 / 0.65 / 0.80 (mix) vs 0.49 / 0.65 / 0.83 (d);
mask J 1.44 / 1.31 / 1.55 vs 1.42 / 1.30 / 1.50; dissolved B/A 2.09 / 2.80 / 5.42 vs 1.95 / 2.57 / 3.43. Noise:
background σ per bin mix/d median 1.10 / 1.11 / 1.15 (max 1.18 / 1.36 / 1.40) — below the √2 the kernel alone
predicts (partition of unity helps); gas SNR at A 8.0 / 8.3 / 2.1 vs d 9.8 / 10.9 / 3.8 (016PG mix bin 0 at SNR 2 =
the price). Sharpness proxy (|∇G| on the lung-boundary shell) d ≈ s ≈ mix within 1–5 % — it does not resolve the
0.5-vs-1-vox inspiratory diaphragm spread; only the registration cascade bench can (PCA's lane, on Hooman's go).
Old p (clock broken) for content only: 1.68 / 1.54 / 1.21. Nothing adopted.
Reading for Hooman: on regular breathers mix recovers the per-breath B/A that s loses (014KS 4.06 vs 2.96) at ~10 %
noise; on 4–12 s breathers it over-shoots and bin 0 starves; dissolved D/G and mask J do not depend on the binning.
- Addendum (XeCS, 10-27): the 1 % rate is Tyger's synthetic ilvtime (17.200 ms) vs the MDH gas-interleave spacing
  17.3715 ms (a = 0.99013), not the pneumotach clock; offset = T_cal 26.73 s (30 cal reps) + window term. Converter
  patch model = known constants; PCA's fit = check. All "t_B in seconds" so far are on the synthetic grid (×1.00997
  true). XeCS pilot: p as shipped has phase concentration 0.23–0.52 vs 0.83–0.90 on the true clock; content peak
  p − d = +2.9 / −0.3 / +7.2 bins (random per session), s − d = +0.6…+0.7 bins.

## 23. XeCS d/p/s anchor-offset test on 53 clock-aligned sessions (Hooman's direct order; XeCS `outputs/dps_offset/`)

- **p/ stacks (XeCS F432)**: pneumotach offset Δ = 26.8 s [26.7, 27.1] (32.5–32.7 s inside the 2024-07-29..09-10
  window, 9 sessions) ≈ 6.9 breaths. Breath-phase concentration of the p bins R = 0.36 [0.22, 0.49] as shipped vs
  0.95 on the true clock (d 0.89, s 0.97); 47/53 sessions R < 0.6; p/ image content fundamental 0.16 of the mean vs
  0.47 for d/ and s/ (3× smeared); p peak − d peak +0.84 [−3.3, +4.4] bins, uncorrelated with Δ/T mod 1 (rank binning
  + in-scan drift scramble it, not a plain shift). Every session recovers on the true clock → the converter patch with
  the known constants (T_cal, window b, 17.200/17.3715) is sufficient; no free fit.
  **Rule: no old-tree p/ image enters any table until p is rerun on the fixed converter.**
- **Anchors (XeCS F433, 0.258 s/bin)**: s bin 0 (strict k0 min) sits +0.170 s [+0.064, +0.304] = +0.70 bins AFTER
  inspiration onset, constant across bins (the strict ±50-ilv minimum lands late on the flat trough); k0 peak − flow
  reversal −0.167 s; nav EI turnaround −0.030 s, nav EE turnaround −0.102 s; d inhale/exhale boundary at rank bin
  7.5 [7.0, 7.8]. In the images the s/ content peak is +1.20 bins [+0.70, +1.93] after the d/ peak (30/53 > 1 bin):
  an anchor/axis difference present before any physiology; b44d vs old-tree agree to 0.01 bin for d and s.
  **Rule: any same-bin d-vs-s number carries the ~1.2-bin offset.** Applies to MIX too (same k0-min anchor as s →
  MIX bin 0 is ~0.17 s after inspiration onset; τ_I and t_B are measured from that late anchor).

## 24. Pneumotach clock patch APPLIED (Hooman's order, 2026-10-27) — fork commit 0209c1a

Model with no free parameters: the converter reads the MDH timestamps (2.5 ms units) and writes the true time of every
gas interleave since the first acquired line (`pneumo_gas_true_s`); raw.py maps the pneumotach clock (+ record-start
offset, −5.8 s inside 2024-07-29..09-10, else 0) onto the synthetic ilvtime by piecewise-linear interpolation between
the gas interleaves BEFORE cropping. Header knobs `pneumoclock` (default 0 = old behaviour) and `pneumo_offset_s`.
PCA's minima-to-EE fit is vendored (`pneumoclock.py`) as the CHECK only; its result is exported with the recon
(`pneumo_clock_chk_a/b/n_matched/n_ee/resid_sd_s/ok`, plus `pneumo_clock_T_cal_s`, `pneumo_clock_ilv_true_ms`,
`pneumo_time_aligned`, `pneumo_vol_aligned`).
Offline on 015HK: T_cal 26.733 s (XeCS 26.7325), true spacing 17.3714 ms (XeCS 17.3715); check BEFORE mapping a 0.991
b −32.4 s 21/29 matched (= PCA's table); AFTER mapping a 1.0003 b −0.10 s 29/29 matched, resid 0.07 s — the −0.10 s
is the k0-min vs navigator-turnaround anchor offset of §23, not a clock residual. p reruns for 015HK / 016PG / 014KS
on image 0209c1a → `Ext/AIkill_Dynamic_mix/<key>/p/` (run_mix3.sh, MODES=p, --pneumoclock).
- XeCS F434 (direct, 10-27): d vs s with the F433 anchor removed — a per-session constant shift (median 1.20 bins =
  0.29 s) makes d/s 0.97–1.01 in every bin (rms 0.032 vs 0.170 raw), B/A d 2.72 = s 2.74; mapping through pressure-
  phase fraction does worse (+0.48 bin, ±10 %; content follows absolute time since EE, F431). On the existing recons d
  and s are the same global curve once the anchor is removed. XeCS per-session clock constants: XeCS
  outputs/cohort/clock_align.json; the fork patch measures T_cal and the spacing from the MDH timestamps instead.
- XeCS exceptions to the record-start rule (clock_align.json, 53 sessions) — for the p reruns beyond the three:
  REAL offsets (apply the measured offset, not the rule): 2024-04-25_042DR −1.7…−2.3 s; 2024-10-10_019LR −0.3…−0.5 s;
  2024-02-27_009JT ≈ 0 ± 0.5 (surrogates disagree). Cal block ≠ 30 reps (handled by the MDH reader): 005JJ 22 reps
  19.60 s, 009ML 25 reps 22.28 s. Window term confirmed on all 9 in-window sessions (−5.64…−5.95 s). Weak-fit
  sessions (006KL, 040RP, 025VP 04-22, 012EB, 009JT 10-08, 026CH, 027JC): rule fine; if the check flags them, the
  k0 surrogate (−0.05…−0.16 s) is the better reference. The −0.10 s check residual = −(0.17 − 0.06): k0 min vs
  pressure crossing (+0.17) minus volume-min vs pressure-crossing (~0.06).
  TODO for a p batch: `--pneumo-offset-s` per session from clock_align.json for the three real-offset sessions.
- First p rerun on 0209c1a: 015HK mapped (T_cal 26.733 s, spacing 17.3714 ms) but the vendored check module was not
  in the Docker image; 016PG / 014KS fell back to the old clock because the raw file ends on an unpaired gas line
  (map one entry longer than ilvtime). Fixed in fork commit 40c248b (Dockerfile COPY, converter truncates to the
  paired count, raw.py tolerates +1); p reruns repeated on image 40c248b.

## 25. p stacks DELIVERED on image 40c248b (2026-10-27 01:03–01:10) — `outputs/binning_eval_2026-10-26/p_check.txt`

| | 015HK | 016PG | 014KS |
|---|---|---|---|
| clock mapped / T_cal / true spacing | yes / 26.733 s / 17.3714 ms | yes / 26.730 s / 17.3716 ms | yes / 26.730 s / 17.3714 ms |
| check after mapping: a / b / matched / resid | 1.0003 / −0.10 s / 29/29 / 0.07 s | 1.0001 / +0.42 s / 24/24 / 0.13 s | 0.9967 / +0.57 s / 24/28 / 0.13 s |
| binned / total ilv | 7017 / 7017 | 11664 / 11664 | 8945 / 8945 |
| hard count CV | 0.001 | 0.000 | 0.000 |
| Σ_FOV gas / k0_wmean per bin | 0.96–1.05 | 0.99–1.01 | 0.89–1.21 * |

\* 014KS: the p bins are assigned to every interleave including the auto-excluded dead window (zero-signal
interleaves land in bins 0–3 by volume rank); the kernel skips them in the image but my check's k0 weighting does not
(use the `excluded` flag from `<key>_ilvbin_tyger.csv` when reading p/ilvbin). 014KS check b +0.57 s / a 0.9967 is
the dead-window breaths' minima mismatch, still ok (24/28). The +0.4–0.6 s b on 016PG/014KS vs −0.10 on 015HK is
volume-minimum-vs-k0-minimum timing per subject (plateau breathers), not clock.
Tree `Ext/AIkill_Dynamic_mix/<key>/p/` with ilvbin (bin_mode 0), ee_times, bin_weff, pneumo_clock_*,
pneumo_time_aligned / pneumo_vol_aligned, calcb_b, rbc_tp_*. The four-way (s/d/p/mix) is now complete on the three
sessions for PCA's loader next session. Old-tree p/ remains invalid (rule §23).
