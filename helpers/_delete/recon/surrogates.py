"""Respiratory surrogates -> per-interleave volume ilvvol in [0,1] on ilvtime.

Three surrogates, one output shape. SIGNAL and PNEUMOTACH are already computed by
Steve's loader and saved by dump_inputs_dyn.py, so they are plain loads. DIAPHRAGM
is the new one: it reconstructs a navigator image per short window of consecutive
interleaves (CS Use B) and reads the diaphragm edge position out of it -- a sharper
input image than Steve's low-res gridded dyn_usimg_recon means a cleaner edge.

Diaphragm pipeline (port of main.calcLVcb + dyn_usimg_recon, results.py:193 / main.py:277):
  window = `win_ilv` consecutive gas interleaves -> our recon (CG by default, or
  wavelet CS) at low res -> 1D drop-off profile along the S-I axis -> 25%/75%
  crossings -> quadratic fit -> sub-pixel 50% crossing = diaphragm position. savgol
  smooth across windows, interpolate onto every interleave.
"""

import json
import os

import numpy as np
from scipy.signal import savgol_filter, medfilt

import asap_recon as ar
from cs_recon import FinufftForward, wavelet_recon
import sigpy as sp


def load_signal(dump_dir):
    return np.load(os.path.join(dump_dir, "ilvvol_signal.npy"))


def load_pneumo(dump_dir):
    return np.load(os.path.join(dump_dir, "ilvvol_pneumo.npy"))


# Superior-inferior (diaphragm-motion) axis of the recon volume. Determined
# empirically: per-window lung centroid along each axis vs the signal surrogate
# gave axis 2 period == signal period and corr +0.94 (axis 0/1: corr 0.11/0.32).
# So the diaphragm moves along axis 2, NOT axis 0 as an early version assumed.
SI_AXIS = 2


def _si_profile(img, axis=SI_AXIS):
    """1D lung signal profile along the S-I axis (sum over the other two axes)."""
    other = tuple(a for a in range(3) if a != axis)
    return np.sum(np.abs(img), axis=other)


def _centroid_position(img, axis=SI_AXIS):
    """Signal centre-of-mass along the S-I axis -- robust at low SNR (never fails).

    Whole-lung profile (not a single boundary), median-thresholded so the noise
    floor doesn't drag the centroid to the FOV centre. This is the measurement that
    correlated +0.94 with the signal surrogate; the edge fit is kept for reference.
    """
    prof = _si_profile(img, axis).astype(float)
    prof = prof - np.median(prof)
    prof[prof < 0] = 0
    s = prof.sum()
    if s <= 0:
        return None
    return float((np.arange(img.shape[axis]) * prof).sum() / s)


def _halfmax_edges(img, axis=SI_AXIS):
    """Both lung boundaries = sub-pixel half-max crossings of the S-I profile.

    Background-subtracted profile, linear-interpolated crossings of 0.5*max on the
    rising (superior) and falling (inferior) sides. The diaphragm is whichever of
    the two moves with breathing -- diaphragm_curve picks the higher-variance side.
    Returns (lo, hi) in axis-voxel units, or (None, None).
    """
    p = _si_profile(img, axis).astype(float)
    p = p - np.median(p)
    p[p < 0] = 0
    if p.max() <= 0:
        return None, None
    thr = 0.5 * p.max()
    above = np.where(p > thr)[0]
    if above.size < 1:
        return None, None
    i, j = above[0], above[-1]

    def cross(a, b):  # linear interp where p crosses thr between idx a (<thr) and b
        if a < 0 or p[b] == p[a]:
            return float(b)
        return float(a + (thr - p[a]) / (p[b] - p[a]))
    lo = cross(i - 1, i) if i > 0 else float(i)
    hi = cross(j + 1, j) if j < len(p) - 1 else float(j)
    return lo, hi


def _edge_position(img, axis=SI_AXIS):
    """Sub-pixel diaphragm edge along the S-I axis (port main.py:290-306)."""
    dropoff = _si_profile(img, axis)
    nav_N = img.shape[axis]
    lo, hi = dropoff.min(), dropoff.max()
    if hi <= lo:
        return None
    p1 = next((p for p in range(nav_N - 1, 0, -1)
               if dropoff[p] > lo + 0.25 * (hi - lo)), 0)
    p2 = next((p for p in range(p1, 0, -1)
               if dropoff[p] > lo + 0.75 * (hi - lo)), 0)
    if p1 - p2 < 2:
        return None
    p = np.polyfit(np.arange(p2, p1), dropoff[p2:p1], 2)
    mid = (hi + lo) / 2
    disc = p[1] ** 2 - 4 * p[0] * (p[2] - mid)
    if disc < 0 or p[0] == 0:
        return None
    det = np.sqrt(disc)
    for m in ((-p[1] + det) / (2 * p[0]), (-p[1] - det) / (2 * p[0])):
        if p2 < m < p1:
            return float(m)
    return None


def diaphragm_surrogate(dump_dir, **kw):
    """Convenience: just the per-interleave ilvvol (see diaphragm_curve)."""
    return diaphragm_curve(dump_dir, **kw)["ilvvol"]


def diaphragm_curve(dump_dir, win_ilv=20, nav_N=64, method="cg", metric="edge",
                    prefer="auto", smooth_win=5, cg_iters=20, nav_max_iter=30, verbose=True):
    """Diaphragm motion curve via per-window navigator recons (CS Use B).

    method : 'cg' (default; fast, iterative -- already sharper than Steve's grid)
             'wavelet' (true CS nav; computes a per-window DCF, slower).
    metric : 'edge' (default; half-max lung boundary -- see prefer for which one)
             'centroid' (S-I signal centre-of-mass: mid-lung, robust at low SNR)
    prefer : 'auto' (default; pick the boundary with higher |corr| to signal surrogate.
             On 025JC this is the LO/apex edge, corr 0.83 vs hi 0.63. LO tracks
             breathing cleanly and stays in-FOV. This is the binning surrogate.)
             'hi'  (inferior/dome edge = anatomical diaphragm dome -- this is what
             the nav_movie cyan-dashed line displays. Clips out of FOV at deep
             inspiration on 025JC -> only 40% valid windows, 2x-harmonic period.
             Good for DISPLAY, not for binning surrogate.)
             'lo'  (superior/apex edge -- force lo regardless of corr)
    win_ilv : interleaves per nav window. Smaller = more windows/breath (finer
             temporal sampling) but lower nav SNR. 20 -> ~10 windows/breath here.
    smooth_win : savgol window (samples) for the displayed/binning curve. MUST be
             shorter than the breathing period in windows or it averages the breath
             away -- the original 11 over-smoothed. 5 keeps the oscillation.
    Returns dict:
      times    : (n_edge,) window-center time (s) of each measured edge
      pos_raw  : (n_edge,) raw sub-pixel diaphragm position (nav-voxels, axis 0)
      pos_smooth : (n_edge,) savgol-smoothed positions
      ilvvol   : (ntotalilvs,) [0,1], NaN where no window covered it (binning input)
      ilvtime  : (ntotalilvs,) interleave times
    """
    meta = json.load(open(os.path.join(dump_dir, "meta.json")))
    npts, nilv, MS = meta["npts"], meta["ntotalilvs"], meta["MS"]
    traj = np.stack([np.load(os.path.join(dump_dir, f"traj{a}.npy")) for a in "xyz"], 1)
    y = np.ascontiguousarray(np.load(os.path.join(dump_dir, "acq_dyn.npy"))[0],
                             dtype=np.complex128)
    ilvtime = np.load(os.path.join(dump_dir, "ilvtime.npy"))

    M = npts * nilv
    traj_tiled = ar.tile_traj(traj.astype(float), M)   # grid units, per-sample
    n_win = nilv // win_ilv

    times, cen, elo, ehi = [], [], [], []
    for w in range(n_win):
        i0, i1 = w * win_ilv, (w + 1) * win_ilv
        sl = slice(i0 * npts, i1 * npts)
        yw = y[sl]
        if np.count_nonzero(yw) < 0.5 * len(yw):       # mostly-excluded window
            continue
        tg = traj_tiled[sl]
        if method == "wavelet":
            tr = ar.grid_to_radians(tg, MS)
            A_raw = FinufftForward(tr, nav_N)
            dens = np.abs(A_raw(A_raw.H(np.ones(len(yw), dtype=complex))))
            wd = 1.0 / np.clip(dens, dens.max() * 1e-4, None); wd /= wd.mean()
            Aw = sp.linop.Multiply((len(yw),), np.sqrt(wd)) * A_raw
            L = sp.app.MaxEig(Aw.H * Aw, dtype=np.complex128, show_pbar=False).run()
            c = 1.0 / np.sqrt(L)
            A = sp.linop.Multiply((len(yw),), c * np.sqrt(wd)) * A_raw
            img_cg = ar.recon(tg, yw, method="cg", MS=MS, IS=nav_N, cg_iters=15)
            W = sp.linop.Wavelet((nav_N,) * 3, wave_name="db4")
            lam = 0.01 * np.percentile(np.abs(W(img_cg)), 99)
            img = wavelet_recon(A, c * np.sqrt(wd) * yw, lam, nav_max_iter)
        else:
            img = ar.recon(tg, yw, method="cg", MS=MS, IS=nav_N, cg_iters=cg_iters)
        c0 = _centroid_position(img)
        lo, hi = _halfmax_edges(img)
        times.append(float(np.mean(ilvtime[i0:i1])))
        cen.append(c0 if c0 is not None else np.nan)
        elo.append(lo if lo is not None else np.nan)
        ehi.append(hi if hi is not None else np.nan)
        if verbose and (w % 25 == 0):
            print(f"  nav window {w+1}/{n_win}")

    times = np.asarray(times); cen = np.asarray(cen)
    elo = np.asarray(elo); ehi = np.asarray(ehi)
    # reject edges pinned at the FOV boundary (half-max crossing failed / lung ran
    # off the image) -- those blow-ups otherwise masquerade as huge "motion"
    bb = 1.0
    elo[(elo <= bb) | (elo >= nav_N - 1 - bb)] = np.nan
    ehi[(ehi <= bb) | (ehi >= nav_N - 1 - bb)] = np.nan

    if metric == "centroid":
        pos = cen
    elif metric == "edge_quad":
        pos = cen  # legacy quad fit unused now; centroid fallback
    else:  # 'edge'
        if prefer == "hi":
            pos = ehi
            if verbose:
                print("edge curve: using hi (inferior/dome) boundary [anatomical diaphragm]")
        elif prefer == "lo":
            pos = elo
            if verbose:
                print("edge curve: using lo (superior/apex) boundary")
        else:  # 'auto': legacy corr-based selection
            sig = None
            sp_ = os.path.join(dump_dir, "ilvvol_signal.npy")
            if os.path.exists(sp_):
                sig = np.interp(times, ilvtime, np.load(sp_))
            def corr_to_sig(e):
                if sig is None:
                    return np.nanstd(e)
                m = np.isfinite(e) & np.isfinite(sig)
                if m.sum() < 6:
                    return -1.0
                a, b = e[m] - e[m].mean(), sig[m] - sig[m].mean()
                return abs(float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-30)))
            clo, chi = corr_to_sig(elo), corr_to_sig(ehi)
            pos = elo if clo > chi else ehi
            if verbose:
                print(f"edge curve: |corr w/ signal| lo {clo:.2f} vs hi {chi:.2f} -> "
                      f"using {'lo' if clo > chi else 'hi'} boundary")
    ok = np.isfinite(pos)
    times, pos = times[ok], pos[ok]
    if len(pos) < 6:
        raise SystemExit(f"diaphragm: only {len(pos)} valid positions -- check nav recon")
    times = np.asarray(times); pos_raw = np.asarray(pos)
    w = min(smooth_win, len(pos_raw)); w = w if w % 2 else w - 1
    pos_smooth = savgol_filter(pos_raw, w, 2) if w >= 3 else pos_raw.copy()
    ilvvol = np.interp(ilvtime, times, pos_smooth, left=np.nan, right=np.nan)
    good = np.isfinite(ilvvol)
    lo, hi = np.nanmin(ilvvol), np.nanmax(ilvvol)
    ilvvol[good] = (ilvvol[good] - lo) / (hi - lo + 1e-30)
    print(f"diaphragm: {len(pos_raw)} nav edges over {n_win} windows; "
          f"covered {int(good.sum())}/{nilv} interleaves")
    return {"times": times, "pos_raw": pos_raw, "pos_smooth": pos_smooth,
            "ilvvol": ilvvol, "ilvtime": ilvtime}
