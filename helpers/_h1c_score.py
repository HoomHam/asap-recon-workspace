"""Scratch (h1 challenge): judging-kit v2 scores for one ASAP arm.

metrics(a, b)        items 1a, 1b, 2, 3, 4, 5 from complex half images (16, 100, 100, 100) in the Tyger frame
run_arm(tag, ...)    three reconstructions (half A, half B, scored training set) with the 4D least-squares engine,
                     scores incl. held-out lines (6) and held-out interleaves (6b); saves crops + json
usage: _h1c_score.py <tag> <lam_t> <lam_s> <iterations> [exp_filter_ms]
"""
import os, sys, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _h1c as H, _h1c_dyn as DY

KIT2 = '/Volumes/HoomHamExt/Work/Codes/2026_XeCS_Recon/h1challenge/shared/judging_kit_v2.npz'
COL = (slice(45, 82), slice(50, 60), slice(26, 38))


def dome(mag):
    """mag (SI, AP, LR) -> (half-crossing position, 10-90 % width) in voxels along SI inside the kit column."""
    p = mag[COL].mean((1, 2))
    ih = int(np.argmax(p)); hi = p[ih]; lo = np.median(p[:18])
    out = []
    for lev in (0.1, 0.5, 0.9):
        t = lo + lev * (hi - lo)
        i = ih
        while i > 0 and p[i - 1] > t:
            i -= 1
        if i == 0:
            return np.nan, np.nan
        out.append((i - 1) + (t - p[i - 1]) / (p[i] - p[i - 1]))
    return out[1], out[2] - out[0]


def metrics(a, b, kit=None):
    kit = np.load(KIT2) if kit is None else kit
    S = (a + b) / 2; D = (a - b) / 2
    res = {k: [] for k in ('noise_lung_over_tissue', 'lung_over_tissue', 'tissue_snr', 'artefact_bg_over_tissue', 'dome_pos_vox', 'dome_width_mm')}
    for ph in range(S.shape[0]):
        P = lambda m: (np.abs(S[ph][m]) ** 2).mean() - (np.abs(D[ph][m]) ** 2).mean()
        sd = lambda m: np.sqrt((np.abs(D[ph][m]) ** 2).mean() / 2)
        T = np.sqrt(max(P(kit['tissue']), 0))
        res['noise_lung_over_tissue'].append(sd(kit['lung_core']) / T)
        res['lung_over_tissue'].append(np.sqrt(max(P(kit['lung_core']), 0)) / T)
        res['tissue_snr'].append(T / sd(kit['tissue']))
        res['artefact_bg_over_tissue'].append(np.sqrt(max(P(kit['background']), 0)) / T)
        pos, wid = dome(np.abs(S[ph]))
        res['dome_pos_vox'].append(pos); res['dome_width_mm'].append(wid * H.DX)
    out = {k: float(np.nanmedian(v)) for k, v in res.items() if k not in ('dome_pos_vox', 'dome_width_mm')}
    pos = np.array(res['dome_pos_vox'])
    out['dome_excursion_mm'] = float((np.nanmax(pos) - np.nanmin(pos)) * H.DX)
    out['dome_width_mm'] = float(np.nanmean(res['dome_width_mm']))
    out['dome_pos_vox'] = [float(v) for v in pos]
    return out


def run_arm(tag, lam_t, lam_s, it, exp_ms=None, d=None, save=True, recon_fn=None, moco=None, post=None):
    kit = np.load(KIT2)
    d = DY.Dyn() if d is None else d
    tr = kit['scored_train_lines']
    ho_il = d.ho_ilv_lines
    win = None if exp_ms is None else np.exp(-(np.arange(H.KILL, H.NPTS) * H.DWELL) / (exp_ms * 1e-3)).astype(np.float32)
    f0 = recon_fn if recon_fn is not None else (lambda lines: d.recon(lines, lam_s=lam_s, lam_t=lam_t, it=it, win=win, verbose=False, moco=moco))
    f = f0 if post is None else (lambda lines: post(f0(lines)))
    t0 = time.time()
    xa = DY.crop(f(np.intersect1d(tr, kit['half_a_lines']))); xb = DY.crop(f(np.intersect1d(tr, kit['half_b_lines'])))
    m = metrics(xa, xb, kit)
    x = f(tr)
    U = d.fwd_all(x)
    m['heldout_lines'] = d.chi2(x, kit['heldout_lines'], U).tolist()
    m['heldout_interleaves'] = d.chi2(x, ho_il, U).tolist()
    m['train'] = d.chi2(x, tr[::6], U).tolist()
    m.update(tag=tag, lam_t=lam_t, lam_s=lam_s, it=it, exp_ms=exp_ms, seconds=time.time() - t0)
    if save:
        os.makedirs(os.path.join(H.EXT, 'arms'), exist_ok=True)
        np.save(os.path.join(H.EXT, 'arms', f'{tag}_scored_crop.npy'), DY.crop(x).astype(np.complex64))
        np.save(os.path.join(H.EXT, 'arms', f'{tag}_halfA_crop.npy'), xa.astype(np.complex64)); np.save(os.path.join(H.EXT, 'arms', f'{tag}_halfB_crop.npy'), xb.astype(np.complex64))
        os.makedirs(os.path.join(H.OUT, 'scores'), exist_ok=True)
        json.dump(m, open(os.path.join(H.OUT, 'scores', f'{tag}.json'), 'w'), indent=1)
    return m


def show(m):
    ho, hi = np.array(m['heldout_lines']), np.array(m['heldout_interleaves'])
    tr = np.array(m['train'])
    print(f"{m['tag']:22s} noise/tis {m['noise_lung_over_tissue']:.4f}  lung/tis {m['lung_over_tissue']:.3f}  tisSNR {m['tissue_snr']:6.1f}  bg/tis {m['artefact_bg_over_tissue']:.3f}  "
          f"dome exc {m['dome_excursion_mm']:5.2f} mm width {m['dome_width_mm']:5.2f} mm | held-out lines low/mid/high {ho[0,0]:9.0f} {ho[0,1]:7.1f} {ho[0,2]:6.2f} | "
          f"held-out ilv {hi[0,0]:9.0f} {hi[0,1]:7.1f} {hi[0,2]:6.2f} | train {tr[0,0]:9.0f} {tr[0,1]:7.1f} {tr[0,2]:6.2f}  ({m['seconds']:.0f} s)", flush=True)


if __name__ == '__main__':
    tag, lam_t, lam_s, it = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
    show(run_arm(tag, lam_t, lam_s, it, float(sys.argv[5]) if len(sys.argv) > 5 else None))
