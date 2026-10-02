# ASAP side — how the two cines were made (no sparsity prior)

Data: 006KL 2024-03-12 MID174, 1H free breathing, 13000 lines = 15.6 passes of the same 832 interleaves, body coil 2 channels.
Code: `2026_ASAP_Recon/workspace/helpers/_h1c*.py`; one command per arm: `_h1c_pipeline.py real <arm> ...`.

## Common to both arms
1. **Channels.** Thermal noise covariance from differences along the readout in the stopped tail of each line; whitening; only
   the principal virtual channel is used (the second one adds about 4 % at high k and needs a fitted sensitivity).
2. **Forward model.** Image on the full body field (LR 176 x AP 120 x SI 160 voxels of 3.5 mm = 616 x 420 x 560 mm), restricted
   to a body support mask; NUFFT at the measured trajectory on a 2x finer time grid, then the **receiver filter** (estimated
   from the data, `shared/receiver_filter_v1.npz`) and decimation. Without the filter the model is wrong by 8-10 % of the
   signal power at mid and high k, because the trajectory moves about 1.2/350 mm^-1 per 5 us sample and the receiver band
   only covers +-150 mm along the gradient while the body reaches +-300 mm.
3. **Breathing model.** 16 node images on the shared cyclic coordinate (node b at b/16, the agreed phase definition); each
   line uses its own coordinate by linear interpolation between the two neighbouring nodes, so nothing is averaged inside a bin.
4. **Estimator.** Linear least squares, conjugate gradients, 20 iterations, data weights = density compensation to the power 0.75.
5. **Temporal coupling along the motion.** First pass: weakly coupled reconstruction (plain cyclic coupling 0.05, 12 iterations).
   From it a displacement field by demons registration of the end-expiration nodes (15, 0, 1) to the end-inspiration nodes (7, 8),
   and one amplitude per node (1D search). Second pass: quadratic penalty on x_b - W_b x_(b+1), W_b = trilinear warp by the
   amplitude difference times the field. Coupling therefore shares data between phases without flattening the diaphragm.
6. **Spatial regularisation (labelled "prior").** Tikhonov 0.01 weighted by 1/p(r), p = (first-pass node magnitude smoothed by
   a 9 mm Gaussian, relative to its 90th percentile, clipped to 0.1..1) squared: a Gaussian prior whose variance follows the
   low-resolution image ("dark stays dark at 9 mm"). One fixed reweighting, never iterated; quadratic, so the result is linear
   in the data. XeCS ruled it inside "no sparsity prior" on the condition that every arm has a **twin without it**
   (uniform Tikhonov 0.003): `*_noprior`.
7. First pass, motion field and prior are built from the scored training lines only. The display cine uses all steady-state lines.

## The two arms
- **sharp**: motion coupling 0.5, no post-filter.
- **snr**: motion coupling 1.0, then the readout filter exp(-t/1.5 ms) applied to the image as a radial k-space weight
  (same filter as the Tyger rerun of 2026-10-01, expressed in the image).

## Tried and not used (negative results)
second channel with fitted sensitivity (4 % at high k) · phase constraint (XeCS; fails, fat-water phase) · self-navigation of each
line's breathing coordinate (no gain on the phantom) · MP-PCA across phases (keeps 14 of 16 components: nothing here is thermal
noise) · smoothing the body outside the 350 mm box (clearly worse) · stronger coupling than 1 (low k starts to pay).
