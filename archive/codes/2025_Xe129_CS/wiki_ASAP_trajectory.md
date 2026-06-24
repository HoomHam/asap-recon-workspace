# ASAP Trajectory — Technical Reference

**For:** Claude Code CLI sessions working on CS adaptation  
**Source:** kasap.c (2024_Steve_dynamiccode), ASAP Trajectory.md (Obsidian)

---

## What ASAP Is

ASAP = **As Spiral As Possible**. A 3D non-Cartesian k-space trajectory. Not purely radial — each spoke is twisted as fast as hardware allows, making it a 3D spiral-radial hybrid. This improves incoherence for CS compared to simple radial.

---

## Trajectory Structure

Three nested levels:

```
NREPS repetitions
  └── NI interleaves per rep
        └── NPTS readout points per interleaf
```

Total k-space points per acquisition = `NI × NPTS × NREPS`

### Example parameters (proton ACR phantom, 2026):
| Parameter | Value | Meaning |
|-----------|-------|---------|
| NI | 26 | interleaves per repetition |
| NPTS | 512 | readout samples per interleaf |
| NREPS | 32 | total repetitions |
| dt | 10 µs | dwell time |
| at | NPTS × dt = 5.12 ms | acquisition time per interleaf |
| t0 | 40 µs | dead time before readout starts |
| FOV | 350 mm | |
| ms | 160 | matrix size |
| gamma | 42.577 MHz/T | proton (use 11.777 MHz/T for Xe129) |

> **Note:** Confirm NI, NPTS, NREPS, FOV, ms for each specific scan from the sequence parameters.

---

## How the Trajectory Is Generated (kasap.c)

1. **Spoke directions:** `NI` unit vectors (`initv[]`) distributed on the sphere — either Fibonacci (fast) or Thompson-optimized (better uniformity via electrostatic repulsion).

2. **Radial law `kr(t)`:** Monotonic ramp from 0 to `kmax = ms/(2×FOV)`. Shape controlled by power-law parameter `n`. Hardware limits enforced: MAXG = 40 mT/m, MAXS = 150 T/m/s.

3. **Spiral twist:** At each time point, an angular velocity `w = min(wG, wS)` is derived from gradient/slew constraints. A rotation axis `rotvec` precesses around x-hat, and k-vectors rotate around it — creating the spiral character.

4. **Gradient ramp-down:** Last samples tapered to g=0 within slew limit.

5. **Per-rep rotation:** All spokes rotated 90° around a rep-specific axis (`reprot[irep]`), itself drawn from a second hedgehog of `NREPS` directions. This covers k-space across reps.

---

## Output Arrays from kasap.c

Two flat arrays per rep (packed as NPTS × NI, column-major):

- `kx[lin], ky[lin], kz[lin]` — k-space coordinates (units: cycles/m or 1/FOV, check normalization)
- `gx[lin], gy[lin], gz[lin]` — gradients in T/m

**Indexing:** `lin = ir + j * NPTS` where `ir` = sample index, `j` = interleaf index.

Output text files from main():
- `kspacetraj.txt` — header: `NI NREPS NPTS`, then `kx ky kz` per line
- `fancytraj.txt` — same structure for gradients

---

## Key Pipeline Steps Before CS

From the ASAP Recon Pipeline (Obsidian):

1. **Trajectory calibration:** Two-slice phase method for k(t), gradient delay correction (Δx,Δy,Δz), Savitzky-Golay smoothing of phase derivative, enforce k(0)≈0.
2. **Raw data conditioning:** Noise spike repair, interleaf DC offset subtraction, per-leaf gain normalization.
3. **DCF:** Iterative (3–5 CG iterations), warm-started from analytical spiral DCF. Kernel: Kaiser-Bessel (k=5, os=2-3).
4. **Gridding:** NUFFT with z-padding, apodization.

---

## What CS Reconstruction Needs from ASAP

Three arrays (matching `spiral3d_cs_Dt.m` / `spiral3d_cs_DzDt.m`):

| Variable | Content | Size (3D static case) |
|----------|---------|----------------------|
| `ktrajs` | normalized k-space coordinates | [Ntotal × 3] or [3 × Ntotal] |
| `kdatas` | complex raw k-space signal | [Ntotal] or [Ncoils × Ntotal] |
| `kcomps` | density compensation weights | [Ntotal] |

Where `Ntotal = NI × NPTS × NREPS` (all reps combined for 3D static recon).

> **TODO:** Confirm exact raw data format — see `wiki_CS_adaptation.md` for details.
