"""Scratch (2026-10-01): image-based ("d") respiratory surrogate for a 1H dynamic -> synthetic pneumotach.

Uses the low-k navigator frames from _proton_selfnav.py (ALL frames, incl. first pass),
the dome ROI centre from its .npz, tracks the dome edge (+-45 mm window, gradient xcorr),
smooths by 1 frame, volume = -dz (dome descends on inspiration), and writes it in the
vendor pneumotach format via _proton_resp_pneumo.write_pneumo.

usage: _proton_dome_pneumo.py <frames.npy> <selfnav.npz> <out_pneumotach> <out_npz> [ILV=26] [TR=0.0172] [DX=11]
"""
import os, sys
import numpy as np
from scipy.ndimage import gaussian_filter1d
sys.path.insert(0, os.path.dirname(__file__))
from _proton_resp_pneumo import write_pneumo

fr_f, nav_f, out_p, out_npz = sys.argv[1:5]
ILV = int(sys.argv[5]) if len(sys.argv) > 5 else 26
TR = float(sys.argv[6]) if len(sys.argv) > 6 else 0.0172
DX = float(sys.argv[7]) if len(sys.argv) > 7 else 11.0
SI, UP, HALF = 2, 20, 3

fr = np.load(fr_f, mmap_mode='r')
nav = np.load(nav_f)
c = nav['roi_centre']
N = fr.shape[1]
sl = [slice(max(ci - HALF, 0), ci + HALF + 1) for ci in c]
sl[SI] = slice(None)
other = tuple(a for a in range(3) if a != SI)
prof = np.array([np.asarray(f[tuple(sl)]).mean(axis=other) for f in fr])
x = np.arange(N)
xf = np.linspace(0, N - 1, N * UP)
P = np.array([np.interp(xf, x, p) for p in prof])
ref = P[int(np.ceil(832 / ILV)):].mean(0)                    # steady-state reference
g0 = np.gradient(gaussian_filter1d(ref, UP))
zmm = xf * DX
e = np.argmax(np.abs(g0) * ((zmm > 50) & (zmm < 190)))
W = (np.abs(np.arange(len(xf)) - e) <= int(45 / DX * UP)).astype(float)
dz = []
for p in P:
    g = np.gradient(gaussian_filter1d(p, UP)) * W
    cc = np.correlate(g, g0 * W, 'full')
    L = len(g)
    dz.append((np.argmax(cc[L - 1 - 3 * UP:L + 3 * UP]) - 3 * UP) / UP * DX)
dz = gaussian_filter1d(np.array(dz) - np.median(dz), 1.0)
t = (np.arange(len(dz)) + 0.5) * ILV * TR
tt, vv = write_pneumo(out_p, t, -dz)
np.savez(out_npz, t=t, dz=dz, edge_mm=zmm[e])
print(f'dome edge {zmm[e]:.0f} mm, {len(dz)} frames, dz 5-95% {np.percentile(dz,95)-np.percentile(dz,5):.1f} mm, wrote {out_p}')
