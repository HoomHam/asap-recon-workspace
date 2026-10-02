"""Scratch (h1 challenge): ASAP's kit code on both sides' delivered halves (items 1-5), item 3 on the shared region
(kit background mask inside ASAP's undilated support), plus a narrow dome column for the edge width."""
import os, sys, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _h1c as H, _h1c_dyn as DY, _h1c_score as SC
kit = dict(np.load(SC.KIT2)); sup = DY.crop(np.load(os.path.join(H.EXT, 'support_g.npy')).astype(np.float32)) > 0.5
shared = kit['background'] & sup; print('shared background voxels', int(shared.sum()))
kit2 = dict(kit); kit2['background'] = shared
rows = []
for side, d in (('asap', os.path.join(H.EXT, 'deliver', 'real')), ('xecs', '/Volumes/HoomHamExt/Work/Codes/2026_XeCS_Recon/h1challenge/deliver/real')):
    arms = sorted({f[len(side) + 1:-len('_halfA.npy')] for f in os.listdir(d) if f.startswith(side + '_') and f.endswith('_halfA.npy')})
    for arm in arms:
        a = np.load(os.path.join(d, f'{side}_{arm}_halfA.npy')); b = np.load(os.path.join(d, f'{side}_{arm}_halfB.npy'))
        m = SC.metrics(a, b, kit2)
        S = np.abs((a + b) / 2); SC_COL = SC.COL
        SC.COL = (slice(45, 82), slice(53, 57), slice(30, 34)); w = np.nanmean([SC.dome(S[p])[1] for p in range(16)]) * H.DX; SC.COL = SC_COL
        rows.append(dict(side=side, arm=arm, narrow_dome_width_mm=float(w), **{k: v for k, v in m.items() if k != 'dome_pos_vox'}))
        print(f"{side} {arm:14s} lung noise/tis {m['noise_lung_over_tissue']:.4f}  lung/tis {m['lung_over_tissue']:.3f}  tissue SNR {m['tissue_snr']:6.1f}  shared bg/tis {m['artefact_bg_over_tissue']:.3f}  "
              f"excursion {m['dome_excursion_mm']:5.2f} mm  10-90 width kit column {m['dome_width_mm']:5.2f} mm, narrow column {w:5.2f} mm", flush=True)
json.dump(rows, open(os.path.join(H.EXT, 'deliver', 'kit_scores_by_asap.json'), 'w'), indent=1)
