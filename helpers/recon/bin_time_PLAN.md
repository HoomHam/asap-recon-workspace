# bin_time — PLAN & NOTES (parked 2026-07-17)

Status: **NOT built.** Plan + reference demos only. Picking up other problems first.

## Problem

DIAPHRAGM (and PNEUMOTACH) 4D recon stack `rspace(nbins, x, y, z)` is ordered by
**amplitude rank** (`raw.py:153 bin()`), not time. Each bin pools interleaves from
all real breaths across the whole scan → no absolute acquisition time survives.
The pipeline emits **no per-bin time vector today**:
- `results.dyn_recon` returns only `rspace(nbins,…)`.
- navigator dict (`tyger_recon._diaphragm_navigator`) has `nav_time` = per *frame*
  absolute time, nothing per *bin*.

So the 4D "time axis" is currently just bin index 0..15, which is **monotonic with
breath phase but NOT proportional to time**. Any dynamic fit (gas wash-in/out τ,
ventilation rate) on bin index is distorted.

## Key distinction by bintype (do not forget)

| bintype | ilvbin meaning | time vector |
|---------|----------------|-------------|
| SIGNAL | within-breath **time fraction** `(t−EE_k)/(EE_{k+1}−EE_k)` (`raw.py:476-486`) | **trivial**: `t_bin[b] = (b/nbins)·T̄`, uniform. Already time-binned. |
| DIAPHRAGM / PNEUMOTACH | amplitude **rank**, equal-count (`bin()`) | **must average** within-breath time per bin (this doc) |

## The recipe (verified on synthetic curve — see demos)

1. **EE troughs** on `ilvvol[DIA]` via the N-window two-pass minima finder
   (port `raw.py:449-466`; the DIAPHRAGM path currently skips it). Phase=0 anchors.
2. **Per interleave** within-breath elapsed time
   `τ_i = (t_i − EE_before)/(EE_after − EE_before) ∈ [0,1)`.
3. **Soft-weighted average per bin using the SAME gridder weights**
   `w(b,i) = exp(−bindist(b, binctr_i)² / bindist0sq)`, cyclic (`recon.py:104-111`,
   `bindist0sq = 2`, `binctr = ilvbin·nbins`). Consistency with what physically
   entered image bin b.
4. **CIRCULAR mean, not linear** — τ wraps at EE (bin 0 and bin nbins−1 are cyclic
   neighbors). Linear mean corrupts wrap bins (first run gave bin0 = 0.497s, wrong).
   → angle = 2π·τ, average sin/cos, `atan2`, unwrap across (phase-ordered) bins,
   zero-offset at bin 0. ×`T̄` (mean EE-to-EE period) for seconds.

Output: `t_bin[b]` (cycles) + `t_bin_s[b]` (seconds) + per-bin circ-std.

## Caveats (bite hard)

- **Non-uniform in time.** Equal-count clusters bins where diaphragm dwells (near
  EE/EI); spreads where it moves fast. Demo: gaps 0.12–0.30s vs 0.198s uniform.
  MUST use real `t_bin` for any τ/rate fit, never bin index.
- **Phase clock, not real seconds.** Period drifts breath-to-breath; you averaged a
  *normalized* cycle. Report dimensionless phase primary; seconds via `T̄` = nominal.
- **Inhale/exhale asymmetry is real** — don't force symmetry (demo: inhale span
  ~1.56s, exhale ~1.46s).
- **Per-bin spread ≠ 0** — circ-std tells tightness; report it alongside t_bin.

## Reference implementations (in this dir)

Scripts live here (`helpers/recon/`); figures write to
`workspace/outputs/diaphragm_binning/` (rule: figures go to outputs/).
- `diaphragm_bin_demo.py` — reproduces `bin()` verbatim; 4-panel figure of the
  amplitude-rank → representative-cycle mapping. Fig: `outputs/diaphragm_binning/diaphragm_binning.png`.
- `diaphragm_bintime_demo.py` — full time-vector recipe above (EE detect, τ,
  soft-weighted **circular** mean per bin). Fig: `outputs/diaphragm_binning/diaphragm_bintime.png`.
  Both run on a synthetic navigator; swap in real `nav_volume`/`ilvtime` to use.

## TODO when resumed

- [ ] Write `bin_time.py`: `def bin_time(ilvvol, ilvtime, ilvbin, nbins, bindist0sq=2, N=?) -> (t_bin_cyc, t_bin_s, t_bin_std, Tbar)`.
- [ ] Pull real `nav_volume` (= `ilvvol[DIA]`) + `nav_ilvtime` from the navigator
      dict; pick N for the minima window to match real breath length in interleaves.
- [ ] Decide where the time vector rides in the output MRD (per-bin meta / sidecar).
- [ ] Validate on a real dataset vs SIGNAL binning (SIGNAL gives a free uniform-time
      cross-check).
- [ ] Sanity: does soft-weight averaging match the actual per-bin interleave set the
      gridder used? (same binctr, same bindist0sq — yes, but confirm nbins wrap.)
