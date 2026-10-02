"""Scratch (2026-10-01): respiratory surrogate for 1H dynamics, packed as a synthetic pneumotach.

1H k0 carries breathing at only ~1%, under a much larger fixed per-interleave pattern
(the first samples are not exactly at k=0), so Steve's SIGNAL binning bins the pattern.
Here: k0 (early samples, coil-combined) / mean-per-interleave pattern -> high-pass
(7 s moving median) -> smooth -> volume-like trace V(t), t = line * TR from scan start.
Steve's PNEUMOTACH path integrates pressure (cumsum) into volume, so we write dV
in the vendor packet format the converter parses (0xA6 0x20, ms stamp, uint16 P).

usage: _proton_resp_pneumo.py <dat> <out_pneumotach_file> <diag_png> [--flip]
"""
import struct
import sys
import numpy as np
import mapvbvd
from scipy.ndimage import median_filter
from scipy.signal import savgol_filter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

NILV, FS_OUT = 832, 200.0


def k0_resp(dat):
    tw = mapvbvd.mapVBVD(dat, quiet=True)
    tw = tw[-1] if isinstance(tw, list) else tw
    tw.image.flagRemoveOS = False
    tw.image.squeeze = True
    tr = float(tw.hdr.MeasYaps[('alTR', '0')]) * 1e-6
    raw = tw.image.unsorted()
    raw = raw if raw.ndim == 3 else raw[:, None, :]
    k0 = np.sqrt((np.abs(raw[2:6].mean(0)) ** 2).sum(0))          # (lines,)
    il = np.arange(k0.size) % NILV
    ss = np.arange(k0.size) >= NILV                                # pattern from steady state only
    pat = np.array([k0[ss & (il == i)].mean() for i in range(NILV)])
    r = k0 / pat[il]
    win = int(round(7.0 / tr)) | 1
    r = r - median_filter(r, size=win, mode='nearest')            # removes approach-to-SS + drift
    r = savgol_filter(r, int(round(0.5 / tr)) | 1, 2)
    return np.arange(k0.size) * tr, r, tr


def write_pneumo(path, t, v):
    tt = np.arange(0.0, t[-1] + 1.0, 1 / FS_OUT)
    vv = np.interp(tt, t, v)
    vv = (vv - vv.min()) / (vv.max() - vv.min())
    dv = np.diff(vv, prepend=vv[0])
    c = 19.0 / np.abs(dv).max()                                    # P decodes as -20..70: keep dv*c within +-19
    # quantize the cumulative sum so integration does not drift
    q = np.round((np.cumsum(dv * c) + 20.0 * np.arange(1, dv.size + 1)) * 65535 / 90.0)
    u = np.diff(q, prepend=0).astype(np.int64)
    u = np.clip(u, 0, 65535)
    # the vendor parser finds packets by scanning for A6 20 anywhere, so no payload
    # byte pair may form that marker; nudge such values (pressure carry keeps the sum)
    def has_marker(bs):
        return any(bs[i] == 0xA6 and bs[i + 1] == 0x20 for i in range(len(bs) - 1))
    pkt = bytearray()
    carry = 0
    for ti, ui in zip(tt, u):
        ms = int(round(ti * 1000))
        while has_marker(struct.pack('<I', ms)):
            ms += 1
        uu = int(ui) + carry
        carry = 0
        while has_marker(struct.pack('<H', uu)) or (uu & 0xFF) == 0x20 or (uu >> 8) == 0xA6:
            uu += 1
            carry -= 1
        b = bytearray(54)
        b[0], b[1] = 0xA6, 0x20
        b[2:6] = struct.pack('<I', ms)
        b[32:34] = struct.pack('<H', uu)
        pkt += b
    pkt += bytes([0xA6, 0x20]) + bytes(52 + 54)    # pad: parser skips the last packet start
    with open(path, 'wb') as f:
        f.write(bytes(pkt))
    return tt, vv


if __name__ == '__main__':
    dat, out, png = sys.argv[1:4]
    t, r, tr = k0_resp(dat)
    if '--flip' in sys.argv:
        r = -r
    tt, vv = write_pneumo(out, t, r)
    f = np.fft.rfftfreq(r.size, tr)
    P = np.abs(np.fft.rfft((r - r.mean()) * np.hanning(r.size))) ** 2
    b = (f > 0.08) & (f < 0.8)
    print(f'TR {tr*1e3:.1f} ms, {r.size} lines, resp peak {f[b][P[b].argmax()]:.3f} Hz, wrote {out}')
    fig, ax = plt.subplots(2, 1, figsize=(12, 6))
    ax[0].plot(t, r, lw=.6); ax[0].set_xlim(0, 60); ax[0].set_title('k0 respiratory surrogate (first 60 s)')
    ax[1].plot(tt, vv, lw=.6); ax[1].set_title('normalized volume written as pneumotach (full scan)')
    ax[1].set_xlabel('s from scan start')
    fig.tight_layout(); fig.savefig(png, dpi=80)
