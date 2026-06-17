# 4D CS Theory, Design Choices, and Known Limitations

**Date:** 2026-06-17  
**Context:** Written after completing the lam_t sweep (lam_t ∈ {0.003, 0.01, 0.02, 0.05, 0.1, 0.2})
on 025JC diaphragm surrogate. These are the theoretical questions that came up and the honest
answers — including where the current implementation is a known approximation.

**Read next:** `Dynamic_4D_CS_Implementation.md` for the as-built code details.  
**Open tests:** See §6 — these should be run before any publication claim.

---

## 1. Why spatial and temporal CS are handled differently

The objective is:

```
min_X  Σ_b ‖√(w_dcf · m_b) (A x_b − y)‖²
     + λ_s · Σ_b ‖W_db4 x_b‖₁       ← spatial: per-bin wavelet
     + λ_t · ‖D_t X‖₁               ← temporal: circular TV across bins
```

**Spatial (λ_s, wavelet db4 per bin):**
- Denoises each individual bin image independently
- Exploits sparsity of lung structure in wavelet domain (edges → sparse coefficients)
- Tuned from static CS work; λ_s = 0.01 rel is settled
- Applied identically to every bin, no coupling between bins

**Temporal (λ_t, circular TV across bins):**
- Couples the 16 bins via `D_t X` (circular difference along bin axis)
- Penalizes bin-to-bin voxel difference → smooths the cine in the respiratory direction
- Only meaningful in 4D — has no equivalent in static recon
- New knob, requires visual sweep to tune

**Why separate:** Spatial and temporal sparsity live in orthogonal domains. Wavelet captures
intra-bin structure (edges, tissue boundaries). TV captures inter-bin smoothness (respiratory
motion continuity). They compete: λ_t too high → bins look identical (motion erased); λ_s too
high → each bin blurry/blocky. Decoupling them allows independent tuning without coupled artifacts.

---

## 2. Why temporal wavelet does not work

Wavelets require the signal to be sparse in that domain. On the temporal (bin) axis this fails
for two reasons:

**Not enough samples.** 16 bins → 2-3 wavelet levels max before hitting the signal length.
The resulting transform is barely more than a DCT of 16 numbers. Not sparse — just a different
basis.

**Wrong signal model.** Respiratory motion is a smooth oscillation across bins. Smooth signals
are NOT sparse in wavelet domain (they ARE sparse in Fourier). Penalizing ‖W_t X‖₁ would
destroy the motion curve, not the noise.

**Why TV works instead:** TV penalizes Σ|differences|. Smooth respiratory motion → small
differences between adjacent bins → small TV penalty → motion preserved. Random undersampling
streaks → large bin-to-bin changes → penalized. TV is the right prior for slowly-varying signals.

**The right tool for frequency-domain temporal regularization** would be **low-rank** — respiratory
motion lives in a very low-dimensional subspace (1-2 principal components). L+S (low-rank + sparse)
exploits that naturally. That is the principled next step.

---

## 3. Why treating the cine as a 4D spatial volume does not obviously win

The idea: treat (x, y, z, b) as a 4D spatial volume, apply 4D wavelet or 4D TV jointly.

**Problems:**

- **The 4th axis is not spatially coherent.** Spatial dims have smooth isotropic structure.
  The bin axis has respiratory motion — directional, anisotropic (diaphragm translates along z,
  not x or y). A 4D isotropic wavelet imposes wrong physics.

- **Anisotropic "voxel size".** Spatial = 2.4 mm. Temporal = 1/16 of a breath cycle. These
  units are incommensurable without explicit normalization. Mixing them in a joint operator
  conflates physical scales.

- **16 is too few.** A 4D db4 wavelet on (100,100,100,16) collapses to 1 level in the bin
  direction immediately. Haar-like pair only. No multiresolution benefit.

**Where 4D spatial works:** k-t SPARSE, k-t BLAST — many time frames (hundreds), long enough
temporal axis, registered coordinate frame. That regime is fundamentally different from 16 bins.

---

## 4. Known limitations of the current implementation (honest accounting)

### 4a. Temporal TV penalizes real motion

`‖D_t X‖₁ = Σ_b Σ_voxel |x_b[voxel] - x_{b+1}[voxel]|`

This assumes the same voxel = same anatomy across bins. That is only true if the lung does not
translate or deform between bins. It does — the diaphragm moves ~2-3 cm from expiration to
inspiration. At a fixed voxel:
- End-expiration: lung tissue
- End-inspiration: air (or different anatomy)

Temporal TV then penalizes the real respiratory deformation as if it were noise. The λ_t sweep
reveals this directly:
- `lam_t = 0.2`: diaphragm barely moves between bins — motion suppressed
- `lam_t = 0.003`: motion preserved but per-bin streaks return
- The "sweet spot" chosen visually is a compromise, not a clean solution

### 4b. 16 bins is coarse

- `ilvvol_diaphragm` has ~276 measured positions across ~100 s scan
- These collapse to 16 bins via soft membership
- Each bin averages ~17 interleaves — intra-bin motion already blurred before CS
- The 16-bin cine cannot resolve sub-bin motion regardless of λ_t

### 4c. Why we accepted both limitations for now

Not because they are solved — because the alternatives cost more than they currently buy:

| Fix | Cost | Prerequisite |
|-----|------|-------------|
| Motion-compensated TV | Deformable registration across bins | Needs the cine first (chicken-and-egg) |
| L+S low-rank | New operator + solver | ~1 week implementation |
| k-t SPARSE | Raw k-t data, not binned | Different acquisition protocol |
| More bins (B=32) | Per-bin SNR collapses | More data or higher undersampling tolerance |

The current output shows respiratory motion and is cleaner than no regularization.
That is the bar for a first-pass demonstration.

---

## 5. What the lam_t sweep tells us

Ran lam_t ∈ {0.003, 0.01, 0.02, 0.05, 0.1, 0.2} on 025JC diaphragm cine, B=16, lam_s=0.01.
Coronal slice_video generated for each. MP4 file sizes as a proxy:

| lam_t | MP4 size | Interpretation |
|-------|----------|---------------|
| 0.003 | 800 KB | near-baseline; streaks visible |
| 0.01 | 661 KB | mild temporal smoothing |
| 0.02 | 563 KB | moderate |
| 0.05 | 422 KB | current working default |
| 0.1 | 261 KB | heavy smoothing |
| 0.2 | 242 KB | likely over-smoothed; motion risk |

**Decision pending:** visual read by Hooman to pick the optimum. Eye-over-metrics rule applies.

---

## 6. Tests needed before publication claims (open)

These are the quantitative/qualitative validations that would make this publishable:

**Test 1 — Motion preservation across bins**  
Measure diaphragm position (same edge metric from `surrogates.py`) across the 16 bins of the
reconstructed cine. Plot vs bin number. Compare lam_t=0.003 vs lam_t=0.2. The curve should
match the original `ilvvol_diaphragm` shape. Heavy TV will flatten it — that flattening is
quantifiable suppression of real motion.

**Test 2 — Bin count tradeoff**  
Run B ∈ {8, 16, 32}. Plot per-bin SNR and diaphragm excursion vs B. Shows the
SNR-vs-temporal-resolution tradeoff explicitly.

**Test 3 — L+S comparison**  
Implement low-rank + sparse decomposition on the same data. Compare to TV temporal on the same
patient. Demonstrates whether the motion-suppression problem in §4a actually matters quantitatively.

**Test 4 — Motion-compensated TV (later)**  
Register adjacent bins (ANTs or simple affine), apply TV in deformed frame. Closes the §4a gap
properly. Chicken-and-egg: needs the unregistered cine first, which we now have.

---

## 7. Conceptual hierarchy (summary)

```
Current implementation:         wavelet_xyz + TV_t        (decoupled, accepted approximation)
Principled extension:           L+S (low-rank_t + sparse_xyz)
Full solution:                  motion-compensated L+S
k-t methods (different regime): 4D spatial on many frames (hundreds, not 16)
```

The current implementation is at the bottom of this hierarchy. It works for demonstration.
Moving up requires the tests in §6 to justify the cost.
