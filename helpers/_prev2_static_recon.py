"""Scratch: pre-v2 (2022 fancy) trajectory check — time-averaged gas recon of a raw twix with a candidate
calibrated trajectory, next to Faraz's own recon of the same session.

Why a custom reader: 20220812 acqOrder 7 ('GP-DP-DP-Spectra') files interleave 1024-sample spectra
(1 per 31 scans) with 512-sample imaging scans; mapVBVD (Python AND MATLAB) sets NCol=1024 for all
scans and returns garbage for the 512 ones. Here we walk the MDHs ourselves and keep only nsamp==512.

Gas selection: imaging scans in time order, phase-of-cycle (mod ncyc) with the largest median |k0|
= gas. Arm index (--arm mdh, default, = Faraz fa_spiral_dyn_recon 2022 branch): Line + nLin*(Rep mod 32)
from the MDH counters; (--arm count): (gas-scan counter + shift) mod narms (the v2/ASAP rule).

Usage: helpers/.venv/bin/python helpers/_prev2_static_recon.py <twix> <tag> [--traj gp.npy] [--ncyc 3]
       [--faraz img.mat] [--shift 0]
Writes: outputs/prev2_calib_2026-10-03/static_<tag>.png (+ .npy)
"""
import argparse
import pathlib
import struct

import numpy as np

WS = pathlib.Path(__file__).resolve().parents[1]
OUT = WS / 'outputs' / 'prev2_calib_2026-10-03'
TRAJ = WS / 'data/xe/human/traj/fa_spiral_dyn_fancy_v2_20230131_gp.npy'
FOV, N = 350.0, 80


def read_scans(path, nsamp_img=512):
    """Imaging scans (nscan, nch, nsamp) complex64 in acquisition order + MDH timestamps (s), Line, Rep."""
    data, ts, lin, rep = [], [], [], []
    with open(path, 'rb') as fh:
        nmeas = struct.unpack('<II', fh.read(8))[1]
        fh.seek(8 + 152 * (nmeas - 1) + 8)
        off = struct.unpack('<Q', fh.read(8))[0]
        fh.seek(off)
        pos = off + struct.unpack('<I', fh.read(4))[0]
        while True:
            fh.seek(pos)
            h = fh.read(192)
            if len(h) < 192:
                break
            dma = struct.unpack('<I', h[0:4])[0] & 0x1FFFFFF
            ev = struct.unpack('<Q', h[40:48])[0]
            ns, nc = struct.unpack('<HH', h[48:52])
            if ev & 1 or dma == 0:
                break
            if ns == nsamp_img and nc > 0:
                buf = fh.read(dma - 192)
                a = np.frombuffer(buf, dtype=np.uint8).reshape(nc, 32 + 8 * ns)[:, 32:]
                data.append(a.copy().view('<c8').reshape(nc, ns))
                ts.append(struct.unpack('<I', h[12:16])[0] * 2.5e-3)
                ln, _, _, _, _, _, rp = struct.unpack('<7H', h[52:66])
                lin.append(ln)
                rep.append(rp)
            pos += dma
    return np.stack(data), np.array(ts), np.array(lin), np.array(rep)


def pipe_dcf(kr, n_it=12):
    import finufft
    w = np.ones(kr.shape[0], np.complex128)
    for _ in range(n_it):
        g = finufft.nufft3d1(kr[:, 0], kr[:, 1], kr[:, 2], w, (2 * N,) * 3, eps=1e-4)
        back = finufft.nufft3d2(kr[:, 0], kr[:, 1], kr[:, 2], g, eps=1e-4)
        w = w / np.maximum(np.abs(back), 1e-12)
    return np.real(w)


def recon(kd, traj):
    """kd: (narms, nch, nsamp) arm-averaged k-space; traj (narms*nsamp, 3) cycles/mm. SoS image (N,N,N)."""
    import finufft
    kr = (traj * FOV) * 2 * np.pi / N                      # radians on an N grid
    keep = np.ones(kr.shape[0], bool)
    keep.reshape(kd.shape[0], -1)[:, :2] = False           # killpts = 2 (as raw.py)
    kr = np.ascontiguousarray(kr[keep])
    w = pipe_dcf(kr)
    img = 0
    for c in range(kd.shape[1]):
        y = kd[:, c, :].reshape(-1)[keep] * w
        im = finufft.nufft3d1(kr[:, 0], kr[:, 1], kr[:, 2], y.astype(np.complex128), (N,) * 3, eps=1e-6)
        img = img + np.abs(im) ** 2
    return np.sqrt(img)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('twix')
    ap.add_argument('tag')
    ap.add_argument('--traj', default=str(TRAJ))
    ap.add_argument('--ncyc', type=int, default=3)
    ap.add_argument('--shift', type=int, default=0)
    ap.add_argument('--narms', type=int, default=640)
    ap.add_argument('--faraz', default='')
    ap.add_argument('--arm', default='mdh', choices=['mdh', 'count'])
    ap.add_argument('--nlin', type=int, default=20)
    ap.add_argument('--gaperiod', type=int, default=32)
    ap.add_argument('--window', default='', help='a0:n gas scans')
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    d, ts, lin, rep = read_scans(a.twix)
    k0 = np.sqrt((np.abs(d[:, :, 1]) ** 2).sum(1))
    med = [np.median(k0[i::a.ncyc]) for i in range(a.ncyc)]
    gpar = int(np.argmax(med))
    g = d[gpar::a.ncyc]
    glin, grep_ = lin[gpar::a.ncyc], rep[gpar::a.ncyc]
    if a.window:
        a0, n = map(int, a.window.split(':'))
        g, glin, grep_ = g[a0:a0 + n], glin[a0:a0 + n], grep_[a0:a0 + n]
    print(f'{a.tag}: {d.shape[0]} imaging scans, nch {d.shape[1]}; phase medians {np.round(med, 6)} -> gas phase {gpar}; '
          f'{g.shape[0]} gas scans', flush=True)
    traj = np.load(a.traj).reshape(-1, 3)
    ns = d.shape[2]
    narms = traj.shape[0] // ns
    if a.arm == 'mdh':
        arm = (glin + a.nlin * (grep_ % a.gaperiod) + a.shift) % narms
    else:
        arm = (np.arange(g.shape[0]) + a.shift) % narms
    kd = np.zeros((narms, d.shape[1], ns), np.complex128)
    cnt = np.bincount(arm, minlength=narms)
    for i in range(narms):
        kd[i] = g[arm == i].mean(0) if cnt[i] else 0
    img = recon(kd, traj)
    np.save(OUT / f'static_{a.tag}.npy', img.astype(np.float32))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows = [('ours: ' + pathlib.Path(a.traj).stem + f' arm={a.arm} shift{a.shift} win {a.window or "all"}', img)]
    if a.faraz:
        import scipy.io as sio
        fz = sio.loadmat(a.faraz, variable_names=['img_gp'])['img_gp'].mean(-1)
        rows.append(('Faraz img_gp mean over bins', fz))
    fig, ax = plt.subplots(len(rows), 3, figsize=(10, 3.4 * len(rows)), squeeze=False)
    for r, (lab, v) in enumerate(rows):
        c = np.array(v.shape) // 2
        for j in range(3):
            ax[r, j].imshow(np.take(v, c[j], axis=j), cmap='gray', vmin=0, vmax=np.percentile(v, 99.8))
            ax[r, j].set_title(f'{lab}\naxis{j} mid', fontsize=7)
            ax[r, j].axis('off')
    fig.suptitle(f'{a.tag}: time-averaged gas, {g.shape[0]} gas scans', fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / f'static_{a.tag}.png', dpi=90)
    # crude sharpness: normalized gradient energy inside signal
    gm = np.sqrt(sum(x ** 2 for x in np.gradient(img)))
    print(f'{a.tag}: sharpness {float((gm ** 2).sum() / (img ** 2).sum()):.4f}  wrote {OUT / f"static_{a.tag}.png"}')


if __name__ == '__main__':
    main()
