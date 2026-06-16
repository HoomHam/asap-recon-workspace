# ASAP Recon — Session Handoff

**Date:** 2026-06-11 (updated; covers 06-10 → 06-11 sessions)
**State:** Comparison phase complete. Own FINUFFT recon built, validated, and used as arbiter on real phantom data. Both big Faraz questions (geometry, uniformity) resolved.
**Next focus:** regularize CG → fix/delete DCF → CS layer.

---

## Current state (what exists and works)

### Reference docs (committed, `workspace/reference/`)
- `Recon_Comparison_StaticGas.md` — **single source of truth** for Steve-vs-Faraz: 10-step diff, scoring, theory notes, compute-cost analysis, comparison-protocol rules.
- `Recon_Overview_Steve.md`, `Recon_Overview_Faraz.md` — per-pipeline code maps with magic numbers.
- `Physics_Notes.md` — educational non-Cartesian primer, every concept mapped to both codebases.
- `workspace/for_steve.md` — report for Steve (implementation analysis + GPU-necessity claim + 7 advocate questions). **Shared/being shared with Steve; expect his rebuttal — respond using the comparison doc + benchmarks.**

### Recon code (committed, `workspace/helpers/recon/`, run via `workspace/helpers/.venv/bin/python`)
- `asap_recon.py` — `recon(traj, data, sample_weights=None, method='adjoint'|'adjoint_dcf'|'cg')`. FINUFFT type1/2 as Aᴴ/A (adjointness 6e-15). Bins = weight vectors, never architecture.
- `dump_inputs.py` — Steve's npy dumps from any `.dat` via his own loaders, CPU-only, `--fov` flag (his traj class hardcodes 350).
- `steve_kernel_numpy.py` — faithful CPU `cudarecon`/`cudarenorm` (float64; ~1e-6 vs GPU expected, certification pending).
- `convert_calib.py` — cal-struct `.mat` → Steve-format trajectory `.npy` (k in 1/mm).
- `faraz_montage.py`, `faraz_zoom_check.py` — comparison figures + zoom-bug verification (montage layout = Faraz's 6×8, slices 16..63).
- `selftest.py` — synthetic validation, run after any operator change.

### Test data (`workspace/data/v3_fov250/`, gitignored, local)
ACR phantom: `.dat` (512×832, gas-only, 1ch) + v3 FOV250 calibration + Faraz's recon (`recon_io/faraz/faraz_recon.mat`, axes (x,z,y) display-permuted!) + his analysis script copy. Outputs in `recon_io/`: `montage_all.png`, `montage_zoomfix.png`, `zoom_check_metrics.json`.

### Environment
Native arm64 everywhere (Miniforge base since 2026-06-11; project venv `workspace/helpers/.venv` = python3.11 + finufft 2.5.1 + sigpy — keep using the venv, it pins the stack).

## Established results (don't re-derive)

1. **CG ≈ Steve-equiv at corr 0.984** on the phantom — Steve's normalized-average gridder ≈ least-squares inverse for fully-sampled static data.
2. **Geometry: ours/Steve correct** (object 178×180×150 mm vs ACR true 190⌀×148); **Faraz ×1.205 magnified** — root cause his `resizing` block (`~/Hooman/Work/Analysis/2025-09-24_ACR/spiral_gpt_ACR_20250915.m:445-452`, flag line 16): `alpha = ((matsize/2-1)/fov)/kmax_meas = 0.8298` rescales k. Proven: α-emulated CG matches his volume at scale exactly 1.00; corrected-Faraz vs steve-equiv volume corr 0.98. His engine (KB+DCF) is fine; the rescale is in his *analysis script*, dataset-dependent. **Correct fix: discard |k|>Nyquist samples or grid larger — never rescale k.** (His 80-grid clips the trajectory corners at 0.188 > 0.16 1/mm; our IS=100 has headroom — that's why he "needed" it.)
3. **Uniformity:** his interior is genuinely slightly flatter in low frequencies (lowfreq-CV 0.093 vs our 0.110), pixel noise identical (~0.09). Partly his 3.125 mm voxels. NOT a display artifact: for 1-ch 1-rep data his combine reduces to magnitude, so |.| vs |.| is the fair comparison.
4. **Performance (M4 Pro CPU, measured):** one Steve-kernel pass of 10M samples @240³ = 22.4 s; FINUFFT adjoint same data = 0.3 s. 16-bin/8-ch projections: Steve-arch verbatim ~48 min; one-pass restructure ~5 min; FINUFFT CG ~20 min. GPU unnecessary at this scale — claim documented in `for_steve.md`.

## Hard-won lessons (cost us real time — respect them)

- **Send figures, trust Hooman's eye over scalar metrics.** Interior-CV misled twice in one day (missed lowfreq shading; then "improved" while suppression artifacts wrecked the image). Decompose low-freq vs noise; visual verdict gates any image-quality claim.
- **Local phased-real-part is a trap** — `phase_corrected_real` in `faraz_zoom_check.py` kept unused as cautionary reference. A correct version must follow Steve's `calcb`: phase from full-data recon + polynomial fill in low-signal voxels.
- Faraz's mat volumes are display-permuted (x,z,y); both codebases use *forward* FFT k→image (flips vs any ifft-based recon); CG extent-threshold metric breaks on unregularized noise floors (the 250/250/250 artifact in the metrics table).

## Next steps (in order)

1. ~~Regularize CG~~ **DONE 2026-06-11, negative result that settles the roadmap** (`cg_tune.py`, `cg_tune_metrics.json`): Tikhonov λ swept over 4 decades × iters {10,15,20,30} — SNR moves only 19.6→20.9. λI is a no-op on fully-sampled data (AᴴA well-conditioned; it penalizes amplitude, not roughness). Steve's gplb=300 filter applied to our data: +2 SNR (19.5→21.5). Remaining gap to his 28.7 = **bias–variance**: his kernel-regression gridder smooths (biased, low variance), CG is the unbiased LS estimate. The knob that closes it is a smoothing regularizer — i.e. the CS layer itself. CG defaults locked: iters=20, lam=0.
2. ~~DCF fix-or-delete~~ **DELETED 2026-06-11** — `pipe_menon_dcf` + `adjoint_dcf` removed (CG supersedes; tombstone comment in `asap_recon.py` says how to rebuild properly if a fast DCF preview is ever wanted). `recon()` default method is now `'cg'`.
3. **CS layer — now the direct next step** — sigpy L1-wavelet/TV with the FINUFFT operator; soft-bin weights as diagonal W in the data term; temporal regularization across bins. Per the sweep finding, this is also what delivers noise parity with Steve's filtered gridder (expect to beat his SNR at equal or better resolution).
4. **Proper b-map output stage (optional, multi-channel future)** — port Steve's `calcb` phase-fill approach; needed before multi-coil data and fair real-part comparisons.
5. **Confounder-neutralized comparison** (paper-time): `gplb: 0` in `recon_io/meta.json`, killpts both ways on Faraz's side.
6. **One-time Colab certification** of `steve_kernel_numpy` vs GPU `savedbin0.npy` (optional; closes the ~1e-6 float question in `for_steve.md` §4c).
7. **Tell Faraz** about the resizing bug (numbers in Established results #2 and `zoom_check_metrics.json`).
8. **Respond to Steve's advocacy** when it arrives.

Reproduce-everything quick start:
```bash
cd workspace/helpers/recon
../.venv/bin/python selftest.py                                          # operator sanity
../.venv/bin/python compare_baseline.py  ../../data/v3_fov250/recon_io   # arbiter
../.venv/bin/python faraz_zoom_check.py  ../../data/v3_fov250/recon_io   # zoom-fix figure
```

## Rules that always apply

- **Project nature:** benign phantom-data MRI recon code comparison — ordinary signal-processing work (see root `CLAUDE.md` Session Notes). Hooman's model preference for this project: **Fable 5**.

- **NEVER git commit/push in repo root** (`2026_ASAP_Recon/`); pull OK. All git in `workspace/` (own repo, `HoomHam/asap-recon-workspace`).
- Never modify Faraz's codebases (`2023_Faraz_Recon/`, `workspace/data/.../faraz/`) — changes go in our `workspace/helpers/` only.
- Scanner data stays local (`workspace/data/` is gitignored).

## Key files

| File | Why |
|------|-----|
| `workspace/reference/Recon_Comparison_StaticGas.md` | Comparison truth + protocol rules |
| `workspace/reference/Physics_Notes.md` | Concept→code map; §10 fastest orientation |
| `workspace/for_steve.md` | The claim Steve will attack; keep consistent with it |
| `workspace/helpers/recon/README.md` | Scaffold usage + status checklist |
| `workspace/data/v3_fov250/recon_io/zoom_check_metrics.json` | All measured numbers |
| `workspace/archive/fable_5_asap_06-10-2026.txt` | Dialogue export (06-10 session) |

## Suggested skills

- `/handoff update` after each next-step lands
- `intent-layer:intent-layer-maintenance` if reference docs change substantially
