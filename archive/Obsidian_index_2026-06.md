# Obsidian ASAP Recon Notes — Index

> **Revival note (2026-06-10):** Original file (created 2026-06-07 10:41) was accidentally deleted. Reconstructed from ScreenPipe screenshots (deleted after revival). Rows the screenshots cut off are marked `[unverified]`.

**Vault path:** `~/Library/Mobile Documents/iCloud~md~obsidian/Documents/Action/MRI/ASAP Recon/`

All files are primary source — read them directly, do not duplicate content into workspace.

## Files

| File | Contents | When to read |
|------|----------|--------------|
| `ASAP Faraz vs Steve.md` | 10-point deep comparison: kernel, DCF, coil combine, orientation, B0 — with practical recommendations | Comparing implementations, debugging visual differences |
| `ASAP Recon Pipeline.md` | Ideal pipeline: traj calibration, DCF, gridding, coil combine, QC, CS-ready architecture | Planning improvements, CS integration |
| `ASAP Recon Faraz Approach.md` | Full function-by-function map of Faraz's MATLAB pipeline | Understanding Faraz's code |
| `ASAP Trajectory.md` | Trajectory math: ASAP spiral design, k-space coverage, calibration approach | Trajectory debugging, calibration work |
| `ASAP Recon Check w Pauli.md` | Discussion with Pauli on recon validation | Cross-check questions, Pauli collab context |
| `3D PSF Incoherence Measurement.md` | PSF measurement approach for CS incoherence validation | CS framework input |
| `Measuring Coherence in K-space.md` | [contents column cut off in screenshot] | CS / sampling design |
| `Other Coherence Matrices.md` | Additional coherence matrix approaches | CS framework |
| `Trajectory & Calibration Improvements.md` | Calibration improvement plan: phase method, DC offset, gradient delay | [when-to-read cut off] |
| `build_calibration_from_xyz — Script Manual (ASAP Calibration).md` | Manual (ASAP Calibration) [row partially visible] | Calibration |
| `load_rawdata_20250816 Manual.md` | How to load raw data with new loader | Raw data loading |
| `CTFT, DTFT and DFT.md` | Math background: Fourier transforms for MRI | Concepts refresher |

*[unverified — present in vault but not visible in screenshots: `Code Manual.md`, `build_calibration_from_xyz Manual.md`, `build_calibration_from_xyz Script.md`; original index may or may not have listed them]*

## Key Cross-References (Obsidian links)

- Trajectory ↔ Calibration: `ASAP Trajectory` ↔ `Trajectory & Calibration Improvements`
- CS work: `3D PSF Incoherence` ↔ `Measuring Coherence in K-space` ↔ `Other Coherence Matrices`
- Comparison: `ASAP Faraz vs Steve` ↔ `ASAP Recon Faraz Approach`
