"""Scratch (h1 challenge): score a 16-phase reconstruction of the phantom against its truth (item 7 of the kit).

score(mag, truth=...)   mag (16, 100, 100, 100) magnitude in the Tyger frame. One global scale is fitted first
                        (least squares on the tissue mask), so nobody is charged for units.
  err_box / err_lung    rms(|x| - |truth|) / rms(|truth|) in the whole box / in lung_core, median over phases
  lung_bias             mean |x| in lung_core / mean |truth| in lung_core - 1
  vessel_contrast       (mean in vessel voxels - mean in parenchyma voxels) of x over the same of the truth
  dome                  excursion and 10-90 % width as in the kit, next to the truth's own values
usage: _h1c_phantom_score.py <label>=<file.npy> [...]
"""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _h1c as H
from _h1c_score import dome, KIT2

TRUTH = os.path.join(H.EXT, 'phantom', 'truth_v1.npz')


def score(mag, truth_file=TRUTH, key='truth_ksphere'):
    kit = np.load(KIT2); T = np.abs(np.load(truth_file)[key]).astype(np.float32)
    lung, tis = kit['lung_core'], kit['tissue']
    ves = lung[None] & (T > 1.6 * np.median(T[:, lung], axis=1)[:, None, None, None])
    parn = lung[None] & (T < 1.15 * np.median(T[:, lung], axis=1)[:, None, None, None])
    s = (mag[:, tis] * T[:, tis]).sum() / (mag[:, tis] ** 2).sum()
    x = mag * s
    eb, el, lb, vc, pos, wid, tpos, twid = [], [], [], [], [], [], [], []
    for b in range(16):
        d = x[b] - T[b]
        eb.append(np.sqrt((d ** 2).mean() / (T[b] ** 2).mean())); el.append(np.sqrt((d[lung] ** 2).mean() / (T[b][lung] ** 2).mean()))
        lb.append(x[b][lung].mean() / T[b][lung].mean() - 1)
        vc.append((x[b][ves[b]].mean() - x[b][parn[b]].mean()) / (T[b][ves[b]].mean() - T[b][parn[b]].mean()))
        p, w = dome(x[b]); pos.append(p); wid.append(w * H.DX); p, w = dome(T[b]); tpos.append(p); twid.append(w * H.DX)
    pos, tpos = np.array(pos), np.array(tpos)
    return dict(scale=float(s), err_box=float(np.median(eb)), err_lung=float(np.median(el)), lung_bias=float(np.median(lb)), vessel_contrast=float(np.median(vc)),
                dome_excursion_mm=float((np.nanmax(pos) - np.nanmin(pos)) * H.DX), truth_excursion_mm=float((tpos.max() - tpos.min()) * H.DX),
                dome_width_mm=float(np.nanmean(wid)), truth_width_mm=float(np.mean(twid)),
                dome_pos_rms_err_mm=float(np.sqrt(np.nanmean((pos - tpos) ** 2)) * H.DX), vessel_voxels=int(ves.sum() / 16))


def show(lab, m):
    print(f"{lab:26s} err box {m['err_box']:.3f}  err lung {m['err_lung']:.3f}  lung bias {m['lung_bias']:+.3f}  vessel contrast {m['vessel_contrast']:.3f}  "
          f"dome exc {m['dome_excursion_mm']:5.2f} (truth {m['truth_excursion_mm']:.2f}) mm  width {m['dome_width_mm']:5.2f} (truth {m['truth_width_mm']:.2f}) mm  dome pos err {m['dome_pos_rms_err_mm']:.2f} mm", flush=True)


if __name__ == '__main__':
    for a in sys.argv[1:]:
        lab, f = a.split('=', 1)
        show(lab, score(np.abs(np.load(f)).astype(np.float32)))
