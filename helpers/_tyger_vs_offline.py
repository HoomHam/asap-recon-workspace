#!/usr/bin/env python3
"""Scratch (2026-09-25, B23): compare the Tyger recon made with the patched image (commit 6651591: specfit +
k0 split) against the OFFLINE re-split of the old recon (resplit_rbctp.py --fit xecs) for the same session.
Expected: same basis angle and target -> aRBC/aTP maps equal up to float32 round-off of the old export and
the 0.001-rad phase grid. Also prints the rbc_tp_* meta the new image writes.
usage: python _tyger_vs_offline.py <key> [<new_root>]
"""
import sys
import numpy as np
import scipy.io as sio
import mrd
from pathlib import Path

key = sys.argv[1]
NEW = Path(sys.argv[2] if len(sys.argv) > 2 else '/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/tyger_specfit_2026-09-24') / key / 'd'
OLD = Path('/Volumes/HoomHamExt/AIkill_Dynamic') / key / 'd'

# meta from the new output.mrd
meta = {}
with mrd.BinaryMrdReader(str(NEW / 'output.mrd')) as r:
    r.read_header()
    for it in r.read_data():
        if isinstance(it, mrd.StreamItem.NdArrayComplexFloat) and it.value.meta.get('dissolved_phase_image'):
            meta = {k: v[0].value if v else '' for k, v in it.value.meta.items()}
            dp_new = it.value.data
        elif isinstance(it, mrd.StreamItem.NdArrayFloat):
            for k in ('gas_phase_image', 'gas_phase_magnitude'):
                if it.value.meta.get(k):
                    meta[k + '_shape'] = it.value.data.shape
print('new output.mrd meta:')
for k, v in meta.items():
    print(f'  {k}: {v}')

off = sio.loadmat(OLD / 'recon_resplit_xecs.mat')
aR, aT = off['aRBC'].astype(np.float64), off['aTP'].astype(np.float64)
print(f"offline: dphi_new {float(np.ravel(off['dphi_new_deg'])[0]):.1f} deg, target {float(np.ravel(off['target'])[0]):.4f}, ph_new per bin {np.round(off['ph_new'].ravel(), 3)}")
nR, nT = dp_new.real.astype(np.float64), dp_new.imag.astype(np.float64)
print('shapes new', dp_new.shape, 'offline', aR.shape)
if nR.shape == aR.shape:
    for name, a, b in (('aRBC', nR, aR), ('aTP', nT, aT)):
        rel = np.linalg.norm(a - b) / np.linalg.norm(b)
        c = np.corrcoef(a.ravel(), b.ravel())[0, 1]
        print(f'{name}: rel RMS diff {rel:.4f}, corr {c:.6f}, max|diff|/max|offline| {np.abs(a - b).max() / np.abs(b).max():.4f}')
    # per-bin phase check: the new run's global phase vs the offline one (both absolute on the same rotated image)
    if 'rbc_tp_ph_rad' in meta:
        ph_new = np.array([float(v) for v in meta['rbc_tp_ph_rad'].split()])
        print('per-bin ph: tyger', np.round(ph_new, 3))
        print('per-bin |dphi(tyger - offline)| deg:', np.round(np.degrees(np.abs(np.angle(np.exp(1j * (ph_new - off['ph_new'].ravel()))))), 2))
# gas image: real(F*b) must be unchanged by the patch
g_old = sio.loadmat(OLD / 'recon.mat', variable_names=['gas_phase'])['gas_phase'].astype(np.float64)
if (NEW / 'recon.mat').is_file():
    g_new = sio.loadmat(NEW / 'recon.mat', variable_names=['gas_phase'])['gas_phase'].astype(np.float64)
    print(f'gas_phase new vs old: rel RMS diff {np.linalg.norm(g_new - g_old) / np.linalg.norm(g_old):.2e}')
