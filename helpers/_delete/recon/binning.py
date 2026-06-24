"""Soft respiratory binning: per-interleave surrogate -> per-bin membership weights.

One binner for all three surrogates (signal / pneumotach / diaphragm). Takes the
continuous per-interleave volume `ilvvol in [0,1]` that the surrogate produced and
returns a soft membership matrix M (B, nilv): M[b, i] = how much interleave i
belongs to bin b. These multiply into the recon weight slot (w_b = dcf * M[b]),
exactly handoff decision D3.

Phase, not amplitude. We port raw.bin's insp/exp mirroring (raw.py:153-177) so the
breath unwraps into a monotonic phase coordinate phi in [0,1): expiration occupies
the lower half, inspiration is mirrored into the upper half, so phi runs over a
full cycle (EE -> EI -> EE). phi is treated as CIRCULAR to match the circular
temporal-TV operator (bin B wraps to bin 1) used in the 4D solver.

Soft membership is a circular Gaussian around each bin center, then normalized
across bins to a partition of unity (sum_b M[b,i] = 1) so no interleave is over- or
under-counted. Excluded interleaves get zero membership in every bin.
"""

import numpy as np


def phase_from_surrogate(ilvvol, exclude=None):
    """Continuous breathing phase phi in [0,1) per interleave, via insp/exp mirror.

    Mirrors raw.bin (raw.py:153-177): samples on the falling edge (inspiration in
    Steve's sign convention) are reflected above the expiration range, then the
    rank of the combined value is the phase. NaN / excluded -> NaN.
    """
    v = np.asarray(ilvvol, dtype=float).copy()
    bad = ~np.isfinite(v)
    if exclude is not None:
        bad = bad | np.asarray(exclude, dtype=bool)
    good = ~bad
    vg = v[good]
    if vg.size < 3:
        phi = np.full(v.shape, np.nan)
        return phi

    # breath direction on the GOOD samples (1 = falling edge, Steve's "inspiration")
    breathdir = np.zeros(vg.size)
    breathdir[1:] = (vg[1:] < vg[:-1]).astype(float)
    if vg[1] < vg[0]:
        breathdir[0] = 1
    iv = vg[breathdir == 1]
    ev = vg[breathdir == 0]
    maxiv = 0.0 if iv.size == 0 else iv.max() * 1.00001
    maxev = 0.0 if ev.size == 0 else ev.max()
    vm = vg.copy()
    vm[breathdir == 1] = maxiv + maxev - vm[breathdir == 1]

    # rank -> phase in [0,1)
    order = np.argsort(vm, kind="mergesort")
    rank = np.empty(vm.size)
    rank[order] = np.arange(vm.size)
    phi_g = rank / vm.size

    phi = np.full(v.shape, np.nan)
    phi[good] = phi_g
    return phi


def _circ_dist(a, b):
    """Circular distance on the unit interval, in [0, 0.5]."""
    d = np.abs(a - b)
    return np.minimum(d, 1.0 - d)


def soft_membership(phi, n_bins=16, sigma_bins=0.75, normalize=True):
    """Circular-Gaussian soft membership M (B, nilv) from per-interleave phase.

    sigma_bins: Gaussian width in units of bin spacing (1/B). ~0.75 gives
    neighboring-bin overlap without washing the cine out. NaN phase -> 0 in all bins.
    normalize: scale so sum_b M[b,i] = 1 on valid interleaves (partition of unity).
    """
    phi = np.asarray(phi, dtype=float)
    nilv = phi.size
    valid = np.isfinite(phi)
    centers = (np.arange(n_bins) + 0.5) / n_bins
    sigma = sigma_bins / n_bins

    M = np.zeros((n_bins, nilv))
    pv = phi[valid]
    for b in range(n_bins):
        d = _circ_dist(pv, centers[b])
        M[b, valid] = np.exp(-(d ** 2) / (2.0 * sigma ** 2))
    if normalize:
        colsum = M.sum(axis=0)
        nz = colsum > 0
        M[:, nz] /= colsum[nz]
    return M


def tile_to_samples(M_ilv, npts):
    """(B, nilv) per-interleave membership -> (B, npts*nilv) per-sample.

    Matches acq_dyn.npy layout: reshape(getimg[:, ich, :], npts*nilv, order='F')
    puts each interleave's npts samples in a contiguous block, so np.repeat with
    npts along the interleave axis reproduces it.
    """
    return np.repeat(M_ilv, npts, axis=1)


def membership_from_surrogate(ilvvol, n_bins=16, sigma_bins=0.75, exclude=None):
    """Convenience: surrogate volume -> (phi, M_ilv)."""
    phi = phase_from_surrogate(ilvvol, exclude=exclude)
    M = soft_membership(phi, n_bins=n_bins, sigma_bins=sigma_bins)
    return phi, M
