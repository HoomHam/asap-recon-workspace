"""Scratch: gas-k0 breathing period of the 000HH dynamic twix (never k0-classified — 2022 NO_TRAJ).

Reuses XeCS helpers/asap_k0_classify.py (read_k0, classify) + calspec_extract.detect_cal, same
gas-parity / cal-block logic as the 201-twix cohort table. Adds a sliding-window period so a
fast->slow switch inside ONE file shows up.

Usage: helpers/.venv/bin/python helpers/_hh_k0_scan.py <twix.dat> [...]
Writes: outputs/hh_fastslow_2026-10-03/<tag>_k0.png + k0_hh.csv
"""
import csv
import pathlib
import sys

import numpy as np

XECS = pathlib.Path.home() / 'Hooman/Work/Codes/2026_XeCS_Recon/workspace/helpers'
sys.path.insert(0, str(XECS))
from asap_k0_classify import read_k0, classify     # noqa: E402
from calspec_extract import detect_cal             # noqa: E402

OUT = pathlib.Path(__file__).resolve().parents[1] / 'outputs' / 'hh_fastslow_2026-10-03'


def window_period(t, k0, win=40.0, step=10.0):
    from scipy.signal import find_peaks
    dt = float(np.median(np.diff(t)))
    sm = np.convolve(k0, np.ones(3) / 3.0, mode='same')
    rng = max(float(k0.max() - k0.min()), 1e-30)
    pk, _ = find_peaks(sm, distance=max(int(1.2 / dt), 2), prominence=0.15 * rng)
    tp = t[pk]
    rows = []
    for a in np.arange(t[0], t[-1] - win + 1e-9, step):
        p = tp[(tp >= a) & (tp < a + win)]
        rows.append((a - t[0] + win / 2, float(np.median(np.diff(p))) if p.size >= 3 else np.nan))
    return pk, np.array(rows)


def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    OUT.mkdir(parents=True, exist_ok=True)
    out = []
    for f in map(pathlib.Path, sys.argv[1:]):
        tag = f'{f.parent.parent.name}_{f.parent.name}__{f.name.split("_")[1]}'
        k, ts, nlin, ncha, proto = read_k0(f)
        n_cal, nrep, t_fid, _, _ = detect_cal(ts, nlin)
        kk = k[:, :nrep]
        medl = np.median(kk[:, n_cal:], axis=1)
        par = 0 if medl[0::2].mean() > medl[1::2].mean() else 1
        k0 = kk[par::2, n_cal:].mean(axis=0).astype(float)
        t = t_fid[n_cal:nrep, 0]
        c = classify(t, k0)
        pk, wp = window_period(t, k0)
        row = dict(tag=tag, protocol=proto, ncha=ncha, nlin=nlin, nrep=int(nrep), n_cal=int(n_cal), **c,
                   win_period_min=np.nanmin(wp[:, 1]) if wp.size else np.nan,
                   win_period_max=np.nanmax(wp[:, 1]) if wp.size else np.nan)
        out.append(row)
        print({k_: row[k_] for k_ in ('tag', 'protocol', 'ncha', 'cls', 'dur_s', 'snr', 'n_peaks', 'period_s',
                                      'period_cv', 'win_period_min', 'win_period_max')}, flush=True)
        fig, ax = plt.subplots(2, 1, figsize=(11, 5), sharex=True)
        ax[0].plot(t - t[0], k0, lw=0.7)
        ax[0].plot(t[pk] - t[0], k0[pk], 'r.', ms=4)
        ax[0].set_ylabel('gas k0')
        ax[0].set_title(f'{tag}  {proto}  {c["cls"]}  period {c.get("period_s")} s  (cv {c.get("period_cv")})')
        if wp.size:
            ax[1].plot(wp[:, 0], wp[:, 1], 'o-')
        ax[1].set_ylabel('period, 40 s window [s]')
        ax[1].set_xlabel('time [s]')
        fig.tight_layout()
        fig.savefig(OUT / f'{tag}_k0.png', dpi=110)
        plt.close(fig)
    keys = list(dict.fromkeys(k_ for r in out for k_ in r))
    with open(OUT / 'k0_hh.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(out)
    print('wrote', OUT)


if __name__ == '__main__':
    main()
