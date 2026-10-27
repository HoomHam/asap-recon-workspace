"""Scratch (2026-10-26, Hooman's binning decision): offline evaluation of respiratory binning candidates on the
k0 signal (SIGNAL / s) and the diaphragm navigator volume (DIAPHRAGM / d) of the 5 v3 sessions. No recon.

Candidates (bins -> 16; -1 = unused):
  A  phase      Steve SIGNAL: (t - t_EE,k) / T_k                                  (raw.py)
  B  abstime    PCA: (t - t_EE,k) / T_ref, T_ref = median period, tail -> bin 15 ('B' keeps tail; 'Bx' drops it)
  C  ismail     periodic line-constrained k-means in (phase, amplitude), centroids on the mean breath curve
                (Ismail 2025 Med Phys 52:4318 §2.3, eqs 2-4)
  D  amprank    Steve DIAPHRAGM rule applied to the signal: amplitude rank, inhale/exhale split (raw.bin)
  E  tline      like C but x = absolute time since EE / T_ref (tail wrapped), amplitude = signal
Kernel: HARD (one bin) and SOFT (Steve's exp(-d^2/2), cyclic) where a continuous bin coordinate exists (A, B, D, E);
C is hard only (cluster id).

Metrics per (session, signal, method, kernel) — Hooman's three criteria:
  count_cv        CV of effective count per bin (balance)            lower = better
  amp_rms         RMS of (bin amp - truth(t_bin)) / (B - A), truth = median over breaths of amp at that absolute
                  time after EE; t_bin = kernel-weighted mean time of the bin (realistic amplitude vs time)
  BA_bin / BA_pb  max/min bin amp vs per-breath median peak/trough (amplitude dynamics recovered)
  tB_bin / tB_pb  time of the max bin vs per-breath median t_peak
  amp_wsd         mean over bins of within-bin weighted amplitude sd / (B - A) (representativeness: one state per bin)
  t_wsd           mean over bins of within-bin weighted time sd (s), bins 2-13 only
  unused          fraction of interleaves with no bin
usage: python helpers/_binning_eval.py <key> ...   -> outputs/binning_eval_2026-10-26/{summary.csv, <key>_<sig>_<method>.csv}
"""
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from _bin_timing import steve_bin, steve_ee_minima   # noqa: E402  (verbatim raw.py rules, verified == root)

BT = os.path.join(os.path.dirname(__file__), '..', 'outputs', 'bin_timing_2026-10-26')
OUT = os.path.join(os.path.dirname(__file__), '..', 'outputs', 'binning_eval_2026-10-26')
NB = 16


def rd(p):
    return list(csv.DictReader(open(p)))


def col(rows, name):
    return np.array([float(r[name]) for r in rows])


# ---------------------------------------------------------------- breath bookkeeping
def breath_index(t, ee):
    k = np.searchsorted(ee, t, side='right') - 1
    inside = (k >= 0) & (k < len(ee) - 1)
    k = np.where(inside, k, -1)
    tee = np.full(len(t), np.nan); per = np.full(len(t), np.nan)
    tee[inside] = t[inside] - ee[k[inside]]
    per[inside] = ee[k[inside] + 1] - ee[k[inside]]
    return k, tee, per


# ---------------------------------------------------------------- binning candidates
def bins_phase(t, sig, ee):
    k, tee, per = breath_index(t, ee)
    return np.where(k >= 0, tee / per, -1.0)


def bins_abstime(t, sig, ee, keep_tail=True):
    k, tee, per = breath_index(t, ee)
    tref = np.median(np.diff(ee))
    f = tee / tref
    if keep_tail:
        f = np.minimum(f, (NB - 0.5) / NB - 1e-6)      # tail -> bin 15 (plateau)
    else:
        f = np.where(f >= 1.0, -1.0, f)
    return np.where(k >= 0, f, -1.0)


def bins_amprank(t, sig, ee):
    v = sig.astype(float).copy()
    return steve_bin(v)


def line_kmeans(x, a, periodic=True, iters=60):
    """Ismail §2.3: k-means on (x, a) with periodic x in [0,1), centroids projected onto the mean curve
    y(x) = mean amplitude in each of 16 x-bins (interpolated to a continuous line). Returns cluster id (0..15,
    ordered by centroid x) per point and the centroids."""
    ok = np.isfinite(x) & np.isfinite(a)
    xs, as_ = x[ok], a[ok]
    xc = (np.arange(NB) + 0.5) / NB
    yc = np.array([as_[(xs >= j / NB) & (xs < (j + 1) / NB)].mean() if np.any((xs >= j / NB) & (xs < (j + 1) / NB)) else np.nan
                   for j in range(NB)])
    good = np.isfinite(yc)
    # continuous periodic line, 1600 points
    xl = (np.arange(NB * 100) + 0.5) / (NB * 100)
    if periodic:
        yl = np.interp(xl, np.concatenate([xc[good] - 1, xc[good], xc[good] + 1]), np.tile(yc[good], 3))
    else:
        yl = np.interp(xl, xc[good], yc[good])
    cx, cy = xc.copy(), np.interp(xc, xl, yl)
    lab = np.zeros(len(xs), int)
    for it in range(iters):
        dx = np.abs(xs[:, None] - cx[None, :])
        if periodic:
            dx = np.minimum(dx, 1 - dx)
        d2 = dx ** 2 + (as_[:, None] - cy[None, :]) ** 2
        new = d2.argmin(1)
        if it and np.array_equal(new, lab):
            break
        lab = new
        for kk in range(NB):
            m = lab == kk
            if not m.any():
                continue
            # periodic mean of x: unwrap around the current centroid
            xx = xs[m]
            if periodic:
                xx = cx[kk] + ((xx - cx[kk] + 0.5) % 1.0 - 0.5)
            mx, my = xx.mean() % 1.0 if periodic else xx.mean(), as_[m].mean()
            # project onto the line (nearest line point)
            dl = np.abs(xl - mx)
            if periodic:
                dl = np.minimum(dl, 1 - dl)
            j = np.argmin(dl ** 2 + (yl - my) ** 2)
            cx[kk], cy[kk] = xl[j], yl[j]
    order = np.argsort(cx)
    rank = np.empty(NB, int); rank[order] = np.arange(NB)
    out = np.full(len(x), -1.0)
    out[ok] = rank[lab]
    return out, cx[order], cy[order]


def bins_ismail(t, sig, ee):
    k, tee, per = breath_index(t, ee)
    ph = np.where(k >= 0, tee / per, np.nan)
    a = (sig - np.nanmin(sig)) / (np.nanmax(sig) - np.nanmin(sig))
    lab, cx, cy = line_kmeans(ph, a, periodic=True)
    return np.where(lab >= 0, lab / NB, -1.0)          # cluster id as bin coordinate (hard use only)


def bins_tline(t, sig, ee):
    k, tee, per = breath_index(t, ee)
    tref = np.median(np.diff(ee))
    x = np.where(k >= 0, np.minimum(tee / tref, 1 - 1e-6), np.nan)
    a = (sig - np.nanmin(sig)) / (np.nanmax(sig) - np.nanmin(sig))
    lab, cx, cy = line_kmeans(x, a, periodic=True)
    return np.where(lab >= 0, lab / NB, -1.0)


def bins_twoanchor(t, sig, ee, split=0.4):
    """XeCS proposal F: inspiration on an ABSOLUTE clock (bins 0..split*16), exhale+pause FRACTIONAL (rest).
    tau_I = split * median period (nav-free; == median t_C when split = median t_C/T). Breaths with T < 1.5 tau_I
    excluded. phi = split*(t-tA)/tau_I for t-tA < tau_I, else split + (1-split)*(t-tA-tau_I)/(T-tau_I)."""
    k, tee, per = breath_index(t, ee)
    tau = split * np.median(np.diff(ee))
    insp = tee < tau
    phi = np.where(insp, split * tee / tau, split + (1 - split) * (tee - tau) / np.maximum(per - tau, 1e-6))
    phi = np.where((k >= 0) & (per >= 1.5 * tau), np.clip(phi, 0, 1 - 1e-9), -1.0)
    return phi


METHODS = {'A_phase': bins_phase, 'B_abstime': bins_abstime, 'Bx_abstime_drop': lambda t, s, e: bins_abstime(t, s, e, False),
           'C_ismail': bins_ismail, 'D_amprank': bins_amprank, 'E_tline': bins_tline, 'F_twoanchor': bins_twoanchor}
SOFT_OK = {'A_phase', 'B_abstime', 'Bx_abstime_drop', 'D_amprank', 'F_twoanchor'}
KERNELS = ('hard', 'soft1', 'pou05')      # Steve sigma=1 (bindist0sq=2) ; XeCS sigma=0.5 + partition of unity


# ---------------------------------------------------------------- metrics
def weights(binfrac, kernel):
    n = len(binfrac); ok = binfrac >= 0
    c = binfrac * NB
    w = np.zeros((n, NB))
    for ib in range(NB):
        d = np.minimum.reduce([np.abs(ib - c), np.abs(ib + NB - c), np.abs(ib - NB - c)])
        if kernel == 'soft1':
            w[:, ib] = np.where(ok, np.exp(-d * d / 2.0), 0.0)
        elif kernel == 'pou05':
            w[:, ib] = np.where(ok, np.exp(-d * d / (2 * 0.5 ** 2)), 0.0)
        else:
            w[:, ib] = ok & (np.floor(c + 0.5).astype(int) % NB == ib)
    if kernel == 'pou05':
        w = w / np.maximum(w.sum(1, keepdims=True), 1e-12)
    return w


def wstat(x, w):
    m = np.isfinite(x); ww = w * m[:, None]; s = ww.sum(0) + 1e-12
    mean = (ww * np.nan_to_num(x)[:, None]).sum(0) / s
    sd = np.sqrt(np.maximum((ww * (np.nan_to_num(x)[:, None] - mean) ** 2).sum(0) / s, 0))
    return mean, sd


def truth_curve(t, sig, ee):
    """median over breaths of amplitude at absolute time after EE (0.05 s grid), per-breath A, B, t_B."""
    k, tee, per = breath_index(t, ee)
    grid = np.arange(0, np.nanmax(per) + 0.05, 0.05)
    curves, A, B, tB = [], [], [], []
    for kk in range(len(ee) - 1):
        m = k == kk
        if m.sum() < 10:
            continue
        tt, ss = tee[m], sig[m]
        A.append(ss[0]); B.append(ss.max()); tB.append(tt[ss.argmax()])
        curves.append(np.interp(grid, tt, ss, right=np.nan))
    curves = np.array(curves)
    med = np.nanmedian(curves, 0)
    return grid, med, np.median(A), np.median(B), np.median(tB)


def evaluate(key, signame, t, sig, ee, excl, BA_pb_override=None, pb_norm=False, split=0.4):
    """BA_pb_override: per-breath B/A truth from smoothed picks (breaths_raw B_over_A_smsm median) instead of the
    extreme-value picks of truth_curve (XeCS audit: raw picks are ~3 % high). pb_norm: per-breath weight
    normalisation on the exhale bins (members with t after EE >= tau_I): w *= (T_med - tau)/(T_k - tau), clipped
    [0.5, 2], so a 12 s breath does not out-vote a 4 s one (XeCS (ii))."""
    grid, med, A_pb, B_pb, tB_pb = truth_curve(t, sig, ee)
    BA_truth = BA_pb_override if BA_pb_override else B_pb / max(A_pb, 1e-9)
    k, tee, per = breath_index(t, ee)
    tmed = np.median(np.diff(ee)); tau = split * tmed
    rows = []
    for mname, fn in METHODS.items():
        try:
            bf = fn(t, sig, ee)
        except Exception as e:
            print(f'   {mname}: FAILED {e!r}'); continue
        bf = np.where(excl == 1, -1.0, bf)                     # auto-excluded interleaves never reach an image
        for kernel in (KERNELS if mname in SOFT_OK else ('hard',)):
            soft = kernel != 'hard'
            w = weights(bf, kernel)
            if pb_norm:
                fac = np.clip((tmed - tau) / np.maximum(per - tau, 1e-6), 0.5, 2.0)
                fac = np.where((k >= 0) & (tee >= tau), fac, 1.0)
                w = w * fac[:, None]
            weff = w.sum(0)
            am, asd = wstat(sig, w)
            # bin time: fold the EE wrap per bin around the bin's own weighted-median BIN COORDINATE (not the raw
            # phase: for two-anchor / abs-time bins the phase of bin-0 members spreads 0.02-0.08 vs 0.97-0.99 on
            # irregular breathers and the median lands mid-way — XeCS audit). Members > half a cycle ahead count as
            # negative time to the next EE; members > half a cycle behind count as t + previous period.
            bfc = np.where(bf >= 0, bf, np.nan)
            nxt = per - tee
            prevper = np.full(len(t), np.nan); kk = k >= 1; prevper[kk] = ee[k[kk]] - ee[k[kk] - 1]
            tfold = np.zeros((len(t), NB))
            for ib in range(NB):
                # deterministic anchor: the bin's nominal coordinate ib/NB (no data-dependent median, which flips
                # when the pre-EE members out-weigh the post-EE ones, e.g. 016PG long plateaus)
                dph = bfc - ib / NB
                tfold[:, ib] = np.where(dph > 0.5, -nxt, np.where(dph < -0.5, tee + prevper, tee))
            mw = np.isfinite(tfold) & (w > 0); ws = (w * mw).sum(0) + 1e-12
            tm = (np.nan_to_num(tfold) * w * mw).sum(0) / ws
            tsd = np.sqrt(np.maximum(((np.nan_to_num(tfold) - tm) ** 2 * w * mw).sum(0) / ws, 0))
            tref = np.median(np.diff(ee))
            truth_at = np.interp(np.clip(tm, 0, grid[-1]), grid, med)
            live = weff > 0.5
            amp_rms = np.sqrt(np.nanmean(((am[live] - truth_at[live]) / (B_pb - A_pb)) ** 2))
            ib_max = int(np.argmax(np.where(live, am, -np.inf)))
            am_min = am[live].min()
            rows.append(dict(key=key, signal=signame, method=mname, kernel=kernel,
                             count_cv=weff.std() / weff.mean(), count_min_over_max=weff.min() / weff.max(),
                             amp_rms=amp_rms, BA_bin=am[live].max() / max(am_min, 1e-9), BA_pb=BA_truth,
                             BA_pb_raw=B_pb / max(A_pb, 1e-9),
                             tB_bin=tm[ib_max], tB_pb=tB_pb, amp_wsd=np.mean(asd[live]) / (B_pb - A_pb),
                             t_wsd=np.mean(tsd[live]), unused=np.mean(bf < 0), tref=tref, empty_bins=int((~live).sum())))
            with open(f'{OUT}/{key}_{signame}_{mname}_{kernel}.csv', 'w') as f:
                f.write('bin,w_eff,amp_wmean,amp_wsd,t_wmean_s,t_wsd_s,truth_amp_at_t\n')
                for b in range(NB):
                    f.write(f'{b},{weff[b]:.1f},{am[b]:.6g},{asd[b]:.6g},{tm[b]:.3f},{tsd[b]:.3f},{truth_at[b]:.6g}\n')
    return rows


def main(keys):
    os.makedirs(OUT, exist_ok=True)
    allrows = []
    for key in keys:
        ty = rd(f'{BT}/{key}_ilvbin_tyger.csv')
        t = col(ty, 't_s'); excl = col(ty, 'excluded').astype(int)
        ee_sig = np.loadtxt(f'{BT}/{key}_ee_tyger.csv', skiprows=1)
        k0 = col(rd(f'{BT}/{key}_ilv_raw.csv'), 'k0_gas_sum8')
        dia = rd(f'{BT}/{key}_ilvbin_diaphragm_tyger.csv')
        navv = col(dia, 'nav_volume')
        # navigator EE: Steve's strict-minimum rule on the navigator volume itself (+-50 ilv), NaN -> +inf
        vv = np.where(np.isfinite(navv), navv, np.inf)
        ee_nav = t[steve_ee_minima(vv, t, 50)]
        print(f'== {key}: k0 EE {len(ee_sig)}, nav EE {len(ee_nav)}, excluded {excl.sum()}')
        allrows += evaluate(key, 's_k0', t, k0, ee_sig, excl)
        allrows += evaluate(key, 'd_nav', t, navv, ee_nav, excl | (~np.isfinite(navv)).astype(int))
    with open(f'{OUT}/summary.csv', 'w') as f:
        wr = csv.DictWriter(f, fieldnames=list(allrows[0].keys())); wr.writeheader(); wr.writerows(allrows)
    # compact print: per signal x method x kernel, median over sessions
    import collections
    agg = collections.defaultdict(list)
    for r in allrows:
        agg[(r['signal'], r['method'], r['kernel'])].append(r)
    print(f'\n{"signal":6s} {"method":16s} {"kern":5s} | count_cv  min/max | amp_rms | BA_bin/BA_pb | tB_bin-tB_pb(s) | amp_wsd | t_wsd(s) | unused')
    for (s, m, kk), rs in sorted(agg.items()):
        g = lambda f: np.median([r[f] for r in rs])
        ba = np.median([r['BA_bin'] / r['BA_pb'] for r in rs]); dt = np.median([r['tB_bin'] - r['tB_pb'] for r in rs])
        print(f'{s:6s} {m:16s} {kk:5s} | {g("count_cv"):7.3f}  {g("count_min_over_max"):7.2f} | {g("amp_rms"):7.3f} | '
              f'{ba:12.2f} | {dt:15.2f} | {g("amp_wsd"):7.3f} | {g("t_wsd"):8.2f} | {g("unused"):6.2f}')


if __name__ == '__main__':
    main(sys.argv[1:])
