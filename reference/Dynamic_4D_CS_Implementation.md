# 4D Temporal CS for Dynamic Xe-129 Ventilation — AS-BUILT

> This is the **real implementation** (not the plan). Everything below is built,
> run, and validated on real human dynamic data (027JC, 025JC). All code lives in
> `workspace/helpers/recon/`. Nothing is committed to the main repo.

## What it does

Reconstructs a 3-min free-breathing hyperpolarized Xe-129 spiral acquisition into a
**B-bin respiratory cine** using **our** static-CS operator (FINUFFT + sigpy,
DCF-preconditioned, L1-wavelet) extended to 4D with a **circular temporal-TV**
coupling. Bins come from three respiratory surrogates (signal / pneumotach /
diaphragm). Bins enter only as **per-sample weights** `w_b = dcf · m_b` (handoff
decision D3); the only cross-bin coupling is the temporal regularizer.

## Data facts (both patients, same `fa_spiral_dyn_fancy_v3_20240130` sequence)

- **Single coil** (`nch=1`) → no SENSE / `calcb` port needed.
- gas-dissolved acquisition (`ilvperTR=2`); we reconstruct the **gas** cine.
- `MS=240`, `IS=100`, `npts=510`, `killpts=2`, `TR≈17.2 ms`, 832 unique interleaves.
- 027JC: 7272 gas interleaves (644 excluded). 025JC: 6980 (1166 excluded).
- Correct trajectory = `fa_spiral_dyn_fancy_v3_20240130_{gp,dp}.npy` (date in the
  protocol name; 425984 rows = 512×832, killpts→510×832=424320 = proven `meta`).

## Modules (`workspace/helpers/recon/`)

| File | Role | Key params / notes |
|------|------|--------------------|
| `dump_inputs_dyn.py` | `.dat` → dump (all coils, traj, ilvtime, signal+pneumo surrogates, exclude mask, meta). Reads TWIX directly via mapvbvd (copies `_read_twix`/`_parse_pneumotach` to avoid the `mrd` dep). Runs Steve's CPU `load_from_arr`. | `--datdir --gp-traj --dp-traj --pneumotach --seqname --out` |
| `binning.py` | surrogate `ilvvol`→ soft per-bin membership `M(B,nilv)`. Ports `raw.bin` insp/exp mirror → circular phase φ∈[0,1) → circular-Gaussian membership (`sigma_bins`), partition-of-unity normalized. `tile_to_samples` → per-sample. | B=16, sigma_bins=0.75 |
| `cs_recon_4d.py` | 4D operators + solvers. `Finufft4D` (B copies of the 3D op), `temporal_diff_op` (circular Δ along bin axis = `@TVOPDt` port, built from `Circshift−Identity`), `spatial_wavelet_op` (db4, axes=(1,2,3)), `dcf_and_norm`. `recon_4d_baseline` (per-bin wavelet, Stage 3) + `recon_4d_joint` (wavelet+temporal-TV PDHG, Stage 4). **MS is threaded** (don't hardcode 240). | lam_s_rel, lam_t_rel rel. to p99 baseline coeffs |
| `selftest_4d.py` | adjointness of all 3 ops + λ_t coupling sanity. **All PASS.** | `../.venv/bin/python selftest_4d.py` |
| `surrogates.py` | three surrogates → per-interleave `ilvvol`. signal/pneumo = load from dump. **diaphragm = CS-nav** (see below). | — |
| `cine_4d.py` | driver: dump → bin → Stage-3 baseline + Stage-4 joint → montages. | `--surrogate --bins --stage --max-iter --lam-s-rel --lam-t-rel` |
| `diaphragm_check.py` | diaphragm curve QA: raw + smoothed + signal/pneumo/diaph overlay. | `--metric edge --win-ilv 20 --smooth-win 5` |
| `kernel_check.py` | shows the 1D S-I profile (“kernel”) + centroid vs half-max edges. | — |
| `nav_movie.py` | per-window nav-image movie (mp4/gif/montage), diaphragm line overlay. | `--view coronal|sagittal --win-ilv 20` |
| `surrogate_compare.py` | overlay signal/pneumo(/diaphragm) + breathing-period cross-check. | — |

## The 4D objective (as implemented)

For `X = (B,N,N,N)`, single coil:

    min_X  Σ_b ‖ √(c·w_dcf·m_b) (A x_b − y) ‖²
         + λ_s · Σ_b ‖ W_db4 x_b ‖₁        (per-bin spatial wavelet — the 3D winner)
         + λ_t · ‖ D_t X ‖₁                 (circular temporal TV across bins)

- `w_dcf`, `c` from `dcf_and_norm` (= `cs_recon.py` math: `w=1/|A Aᴴ 1|`, `c=1/√maxeig`).
- Joint solve = sigpy PDHG, `G=Vstack([W,D])`, `proxg=Stack([L1,L1])`, warm-started
  from the baseline. DCF baked into A ⇒ PDHG does not stall.
- Cost on CPU finufft: ~adjoint 0.18 s / forward 0.10 s at M=3.7M, N=100. Full
  B=16 baseline+joint ≈ 10–12 min/patient.

## Respiratory surrogates → `ilvvol[0,1]`

1. **signal** — `Σ|FID[:8]|` per interleave (Steve's `raw.py:429`), dumped directly.
2. **pneumotach** — vendor binary → savgol → integrate P→volume (`raw.py:180`), dumped.
3. **diaphragm** — our CS-nav (the hard-won part, see next section).

All three → one `binning.py` soft binner. signal & pneumo agree (~3.6 s, clean).

## Diaphragm method (`surrogates.diaphragm_curve`) — read this, it's subtle

Per window of `win_ilv` **consecutive** interleaves, CG-recon a low-res nav image,
measure an S-I position, smooth, interpolate onto all interleaves.

**Hard-won settled choices (each was a bug we hit):**
- **S-I axis = `axis2`** (`SI_AXIS=2`). Found empirically: per-axis lung centroid vs
  the signal surrogate gave axis2 corr **0.94** (axis0/1: 0.11/0.32). The 1D kernel
  = `Σ|img|` over the other two axes → profile along axis2 (this is the
  "sagittal-collapsed" / full-signal profile; identical whether you think coronal or
  sagittal).
- **`win_ilv=20`** (~10 windows/breath here) and **`smooth_win=5`** savgol. Bigger
  windows (40) undersample the breath; a savgol window **longer than the breath
  over-smooths** (the 11-window default flattened it). **No median filter** — it
  deleted the real sharp inspiration excursions.
- **Curve metric = `edge` with correlation-based boundary selection.** `_halfmax_edges`
  returns both lung boundaries (lo, hi); FOV-pinned values are rejected; the boundary
  with the **higher |corr| to the signal surrogate wins** (on 025JC: lo 0.83 vs hi
  0.63 → lo). This is the **“beautiful” clean curve** (period 3.52 s, smoothed hugs
  raw). Do NOT hardcode hi.
- **hi vs lo nuance:** the **hi (inferior) edge is anatomically the diaphragm dome**,
  but it **clips out of FOV at deep inspiration** → sparse/sharp/noisy curve (period
  reads a 2× harmonic). The lo (apex) edge stays in-FOV → clean. So:
  - **surrogate CURVE / binning** → corr-selected edge (clean).
  - **nav_movie DISPLAY line** → hi edge (shows where the dome is on the image).
  These are two different jobs; keep them separate.

## Orientation — display-only, NEVER the matrix

The recon matrix and every saved array stay in **raw axis2 coordinates**. `nav_movie`
applies only matplotlib display ops — `.T` (S-I vertical) and `invert_yaxis()`
(inferior at bottom) — applied identically everywhere. No `np.flip`/rotation on data.
`--view {coronal,sagittal}` changes only which axis is summed for the 2D picture
(coronal=sum axis1, sagittal=sum axis0); the diaphragm number is identical either way.
nav_movie final overlay = **only the cyan-dashed hi (diaphragm) line**.

## How to run (per patient)

```bash
cd workspace/helpers/recon
PY=../.venv/bin/python    # arm64 venv: finufft 2.5.1, sigpy 0.1.27
# 0. dump
$PY dump_inputs_dyn.py --datdir <027JC dir> \
   --gp-traj <..._v3_20240130_gp.npy> --dp-traj <..._dp.npy> \
   --pneumotach <pneumotach_file> --seqname fancy_v3_20240130 [--no-ref] \
   --out ../../data/<v3_dyn[_025JC]>/recon_io_dyn
# 1. surrogate QA
$PY diaphragm_check.py <dump> --metric edge --win-ilv 20 --smooth-win 5
$PY surrogate_compare.py <dump>            # signal vs pneumo period agreement
$PY nav_movie.py <dump> --view coronal --win-ilv 20   # diaphragm-line movie
# 2. cine (per surrogate)
$PY cine_4d.py <dump> --surrogate signal|pneumo|diaphragm --bins 16 --stage both \
   --lam-s-rel 0.01 --lam-t-rel 0.05
# selftest anytime: $PY selftest_4d.py   (must print ALL PASS)
```

## Outputs

- 027JC `workspace/data/v3_dyn/recon_io_dyn/`: dump + `cine_signal/{cine_baseline.npy,
  cine_joint.npy, montage_baseline.png, montage_joint.png}` (joint = cleaner, streaks
  suppressed, motion preserved).
- 025JC `workspace/data/v3_dyn_025JC/recon_io_dyn/`: dump + `diaphragm_check.png`,
  `kernel_check.png`, `nav_movie_{coronal,sagittal}/`, `ilvvol_{signal,pneumo,diaphragm}.npy`,
  `diaphragm_pos_{raw,smooth}.npy`.

## Status / done

- ✅ Stage 0 (dump + ops selftest), Stage 2 (binning), Stage 3 (baseline cine),
  Stage 4 (joint temporal-TV cine) — validated on 027JC signal.
- ✅ Three surrogates working; signal/pneumo cross-validated; diaphragm CS-nav curve
  beautiful (corr-selected edge).

## Open / next

- Run full cines for pneumo + diaphragm (027JC, 025JC).
- `lam_t` sweep on the joint cine (currently 0.05 rel).
- Apply diaphragm fixes everywhere (027JC diaphragm rerun).
- Later: low-rank temporal option; dissolved phase; retrospective undersampling
  (recompute DCF per mask).

## Files referenced (read-only ground truth)
`raw.py` (`bin`:153, signal:429, pneumo:180), `results.py` (`calcb`:61,
`dyn_usimg_recon`:193, `dyn_recon`:243), `main.py:277` (`calcLVcb`), upstream
trajectory loader `raw.traj`; old MATLAB temporal op `workspace/codes/2025_Xe129_CS/@TVOPDt`.
