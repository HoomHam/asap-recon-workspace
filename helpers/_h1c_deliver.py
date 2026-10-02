"""Scratch (h1 challenge): assemble the ASAP delivery from the finished pipeline runs.

deliver/<src>/asap_<arm>_display.npy   (16, 100, 100, 100) float32 magnitude, Tyger frame, all steady-state lines
deliver/<src>/asap_<arm>_{halfA,halfB,scored}.npy  complex64 (kit items)
deliver/scores_asap.json + scores_asap.md  one row per (source, arm)
src = real | v1 | v1b ; arm = sharp | sharp_noprior | snr | snr_noprior
"""
import os, sys, json, shutil
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _h1c as H

A = os.path.join(H.EXT, 'arms'); DL = os.path.join(H.EXT, 'deliver'); SC = os.path.join(H.OUT, 'scores')
ARMS = ('sharp', 'sharp_noprior', 'snr', 'snr_noprior')
rows = []
for src in ('real', 'v1', 'v1b'):
    for arm in ARMS:
        fj = os.path.join(SC, f'{src}_{arm}.json')
        if not os.path.exists(fj):
            continue
        od = os.path.join(DL, src); os.makedirs(od, exist_ok=True)
        m = json.load(open(fj))
        for part in ('halfA', 'halfB', 'scored'):
            shutil.copyfile(os.path.join(A, f'{src}_{arm}_{part}_crop.npy'), os.path.join(od, f'asap_{arm}_{part}.npy'))
        fd = os.path.join(A, f'{src}_{arm}_display_crop.npy')
        if os.path.exists(fd):
            np.save(os.path.join(od, f'asap_{arm}_display.npy'), np.abs(np.load(fd)).astype(np.float32))
        ho, hi = np.array(m['heldout_lines'])[0], np.array(m['heldout_interleaves'])[0]
        r = dict(src=src, arm=arm, noise_lung_over_tissue=m['noise_lung_over_tissue'], lung_over_tissue=m['lung_over_tissue'], tissue_snr=m['tissue_snr'],
                 background_over_tissue=m['artefact_bg_over_tissue'], dome_excursion_mm=m['dome_excursion_mm'], dome_width_mm=m['dome_width_mm'],
                 heldout_lines=ho.tolist(), heldout_interleaves=hi.tolist(), settings={k: m[k] for k in ('lt', 'lam_s', 'prior', 'wpow', 'it', 'apod_ms')})
        if 'phantom' in m:
            r['phantom'] = {k: m['phantom'][k] for k in ('err_box', 'err_lung', 'lung_bias', 'vessel_contrast', 'dome_excursion_mm', 'truth_excursion_mm', 'dome_pos_rms_err_mm')}
        rows.append(r)
os.makedirs(DL, exist_ok=True)
json.dump(rows, open(os.path.join(DL, 'scores_asap.json'), 'w'), indent=1)
with open(os.path.join(DL, 'scores_asap.md'), 'w') as f:
    f.write('# ASAP arms, kit v2 scores (thermal units for the held-out items; samples 2-50 / 50-200 / 200-512)\n\n')
    f.write('| data | arm | lung noise / tissue | lung / tissue | background / tissue | dome excursion (mm) | dome 10-90 width (mm) | held-out lines | held-out interleaves |\n|---|---|---|---|---|---|---|---|---|\n')
    for r in rows:
        f.write(f"| {r['src']} | {r['arm']} | {r['noise_lung_over_tissue']:.4f} | {r['lung_over_tissue']:.3f} | {r['background_over_tissue']:.3f} | {r['dome_excursion_mm']:.1f} | {r['dome_width_mm']:.1f} | "
                f"{r['heldout_lines'][0]:.0f} / {r['heldout_lines'][1]:.0f} / {r['heldout_lines'][2]:.1f} | {r['heldout_interleaves'][0]:.0f} / {r['heldout_interleaves'][1]:.0f} / {r['heldout_interleaves'][2]:.1f} |\n")
    ph = [r for r in rows if 'phantom' in r]
    if ph:
        f.write('\n## Phantom, against the truth (k-sphere limited)\n\n| phantom | arm | error box | error lung | lung bias | vessel contrast | dome excursion (truth) mm | dome position error mm |\n|---|---|---|---|---|---|---|---|\n')
        for r in ph:
            p = r['phantom']
            f.write(f"| {r['src']} | {r['arm']} | {p['err_box']:.3f} | {p['err_lung']:.3f} | {p['lung_bias']:+.3f} | {p['vessel_contrast']:.2f} | {p['dome_excursion_mm']:.2f} ({p['truth_excursion_mm']:.2f}) | {p['dome_pos_rms_err_mm']:.2f} |\n")
print(open(os.path.join(DL, 'scores_asap.md')).read())
