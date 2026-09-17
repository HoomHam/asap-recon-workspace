#!/usr/bin/env python
"""CT overlay for the EBV + LTX cohorts (extends ct_mri_overlay.py to clinical CTs).

Usage:
    ct_cohort_batch.py prep  [ct_subject ...]   # DICOM/NIfTI -> NIfTI candidates, lung volume per series, pick max, write manifest
    ct_cohort_batch.py merge                    # merge per-subject manifest_part.json -> manifest.json
    ct_cohort_batch.py list                     # print manifest keys
    ct_cohort_batch.py run   <key> [...]        # refine (fixed orientation F58) + overlay + video for manifest keys
    ct_cohort_batch.py keys                     # one key per line (for xargs -P)

CT series choice: ORIGINAL axial, 0.5 <= thk <= 2.5 mm (fallback <= 3.5), >= 80 slices; every candidate is
converted, the one with the LARGEST lung-mask volume wins (= deepest inspiration, same rule as the RT phase
pick); ties within 3 % -> thicker slices, then softer kernel. Clinical CTs are single breath-hold inspiration
scans, so no phase tag: the chosen series tag is used as `phase`.

MRI: PCA-session elastix one-shot stack (registered_elastix.npy); if missing -> unregistered stack.npy
(prov says so). EI bin = 'auto' (max ventilated volume). Sessions per CT: LTX "(2023)" folders -> every 2023
xenon session of that subject (042DR has none -> 2024-01-17); EBV pre/post -> matching-date session
(006KL: one MRI for both CTs). S1/S2 EBV folders: subject unknown -> skipped.
"""
import glob
import json
import os
import shutil
import subprocess
import sys

import numpy as np
import SimpleITK as sitk

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ct_mri_overlay as ov  # noqa: E402

CT_ROOT = '/Volumes/HoomHamExt/Work/CT'
PCA = '/Volumes/HoomHamExt/Work/Codes/2026_PCA_registration'
ORIENT = ('012', '--+')  # F58

# ct_subject key -> (cohort folder, CT folder name, [xenon sessions])
COHORT = {
    # LTX
    '004LR': ('CT LTX', '004LR (2023)', ['2023-03-10_004LR']),
    '011BA': ('CT LTX', '011BA (2023)', ['2023-05-02_011BA']),
    '012SR': ('CT LTX', '012SR (2023)', ['2023-05-02_012SR']),
    '017AK': ('CT LTX', '017AK (2023)', ['2023-03-08_017AK']),
    '021JM': ('CT LTX', '021JM (2023)', ['2023-04-18_021JM']),
    '023DB': ('CT LTX', '023DB (2023)', ['2023-03-09_023DB', '2023-04-26_023DB']),
    '025VP': ('CT LTX', '025VP (2023)', ['2023-11-02_025VP']),
    '029CK': ('CT LTX', '029CK (2023)', ['2023-03-24_029CK']),
    '030DN': ('CT LTX', '030DN (2023)', ['2023-04-07_030DN', '2023-08-24_030DN']),
    '032WS': ('CT LTX', '032WS (2023)', ['2023-03-31_032WS', '2023-09-25_032WS']),
    '036RL': ('CT LTX', '036RL (2023)', ['2023-11-02_036RL']),
    '037GD': ('CT LTX', '037GD (2023)', ['2023-11-13_037GD']),
    '038RL': ('CT LTX', '038RL (2023)', ['2023-11-16_038RL']),
    '039CP': ('CT LTX', '039CP (2023)', ['2023-11-21_039CP']),
    '040RP': ('CT LTX', '040RP (2023)', ['2023-12-11_040RP']),
    '041WF': ('CT LTX', '041WF (2023)', ['2023-12-19_041WF']),
    '042DR': ('CT LTX', '042DR (2023)', ['2024-01-17_042DR']),
    # EBV
    '002JM_pre': ('CT EBV', '002JM_pre', ['2023-10-25_002JM']),
    '003KH_pre': ('CT EBV', '003KH_pre', ['2023-10-31_003KH']),
    '003KH_post': ('CT EBV', '003KH_post', ['2023-12-15_003KH']),
    '005DS_pre': ('CT EBV', '005DS_Pre', ['2023-11-17_005DS']),
    '006KL_pre': ('CT EBV', '006KL_pre', ['2024-03-12_006KL']),
    '006KL_post': ('CT EBV', '006KL_post', ['2024-03-12_006KL']),
    '007RA_pre': ('CT EBV', '007RA_pre', ['2024-01-22_007RA']),
    '009JT_pre': ('CT EBV', '009JT_pre', ['2024-02-27_009JT']),
    '010AJ_pre': ('CT EBV', '010AJ_pre', ['2024-03-28_010AJ']),
}


def _kernel_num(k):
    import re
    m = re.search(r'(\d+)', str(k))
    return int(m.group(1)) if m else 99


def dicom_candidates(folder, thk_max=2.5):
    import pydicom
    cands = []
    for root, dirs, files in os.walk(folder):
        fs = [f for f in files if not f.startswith('.')]
        if len(fs) < 80:
            continue
        try:
            h = pydicom.dcmread(os.path.join(root, fs[0]), stop_before_pixels=True)
        except Exception:
            continue
        itype = '/'.join(getattr(h, 'ImageType', []))
        iop = [round(float(x)) for x in getattr(h, 'ImageOrientationPatient', [0] * 6)]
        thk = float(getattr(h, 'SliceThickness', 0) or 0)
        if 'ORIGINAL' not in itype or iop != [1, 0, 0, 0, 1, 0] or not (0.5 <= thk <= thk_max):
            continue
        k = getattr(h, 'ConvolutionKernel', '?')
        k = k[0] if isinstance(k, (list, pydicom.multival.MultiValue)) else k
        cands.append(dict(dir=root, n=len(fs), series=str(getattr(h, 'SeriesNumber', '?')), kernel=str(k), thk=thk,
                          tag=f's{getattr(h, "SeriesNumber", "x")}_{str(k).replace(",", "").replace(" ", "")}_{thk:g}mm'))
    if not cands and thk_max < 3.5:
        return dicom_candidates(folder, 3.5)
    return cands


def convert(series_dir, out_path):
    r = sitk.ImageSeriesReader()
    ids = r.GetGDCMSeriesIDs(series_dir)
    files = r.GetGDCMSeriesFileNames(series_dir, ids[0])
    r.SetFileNames(files)
    img = sitk.Cast(r.Execute(), sitk.sitkInt16)
    sitk.WriteImage(img, out_path, True)
    return img


def lung_volume_L(img):
    lung = ov.ct_lung_mask(sitk.Cast(img, sitk.sitkFloat32))
    return float(lung.sum() * np.prod(img.GetSpacing()) / 1e6)


def prep(subjects):
    for cs in subjects:
        manifest = {}
        coh, folder, sessions = COHORT[cs]
        src = os.path.join(CT_ROOT, coh, folder)
        ct_dir = os.path.join(ov.EXT, cs, 'ct'); os.makedirs(ct_dir, exist_ok=True)
        index = {}
        nif = [p for p in (os.path.join(src, 'NIFTI', 'CT.nii.gz'), os.path.join(src, 'CT.nii.gz')) if os.path.exists(p)]
        if nif:
            out = os.path.join(ct_dir, 'nifti.nii.gz')
            if not os.path.exists(out):
                shutil.copy(nif[0], out)
            img = sitk.ReadImage(out)
            index['nifti'] = dict(src=nif[0], size=list(img.GetSize()), spacing=[round(x, 3) for x in img.GetSpacing()],
                                  lung_L=round(lung_volume_L(img), 3), kernel='?', thk=float(img.GetSpacing()[2]))
            print(cs, 'nifti', index['nifti'])
        else:
            for c in dicom_candidates(src):
                out = os.path.join(ct_dir, f'{c["tag"]}.nii.gz')
                img = sitk.ReadImage(out) if os.path.exists(out) else convert(c['dir'], out)
                index[c['tag']] = dict(src=os.path.relpath(c['dir'], src), n=c['n'], kernel=c['kernel'], thk=c['thk'],
                                       size=list(img.GetSize()), spacing=[round(x, 3) for x in img.GetSpacing()],
                                       lung_L=round(lung_volume_L(img), 3))
                print(cs, c['tag'], index[c['tag']]['lung_L'], 'L')
        if not index:
            print(cs, 'NO CT CANDIDATES'); continue
        vmax = max(v['lung_L'] for v in index.values())
        best = sorted([t for t, v in index.items() if v['lung_L'] >= 0.97 * vmax],
                      key=lambda t: (-index[t]['thk'], _kernel_num(index[t]['kernel'])))[0]
        index['_chosen'] = best
        json.dump(index, open(os.path.join(ct_dir, 'index.json'), 'w'), indent=2)
        for sess in sessions:
            el = f'{PCA}/elastix_oneshot_2026-09-16/{sess}/registered_elastix.npy'
            st = f'{PCA}/xecs_registrations_2026-09-15/{sess}/stack.npy'
            mri, prov = (el, 'elastix one-shot PCAMetric2 (PCA session)') if os.path.exists(el) else (st, 'UNREGISTERED Tyger stack (no elastix result)')
            key = f'{cs}__{sess}'
            manifest[key] = dict(mri=mri, ei='auto', phase=best, ct_subject=cs, cohort=coh.split()[-1], session=sess,
                                 prov=f'{prov}; CT {coh}/{folder} series {best} (lung {index[best]["lung_L"]} L, max of {len(index) - 1})')
            print('  manifest', key, '->', best)
        json.dump(manifest, open(os.path.join(ov.EXT, cs, 'manifest_part.json'), 'w'), indent=2)


def merge():
    manifest = json.load(open(ov.MANIFEST)) if os.path.exists(ov.MANIFEST) else {}
    for cs in COHORT:
        part = os.path.join(ov.EXT, cs, 'manifest_part.json')
        if os.path.exists(part):
            manifest.update(json.load(open(part)))
    json.dump(manifest, open(ov.MANIFEST, 'w'), indent=2)
    print(len(manifest), 'keys in manifest')


def run(keys):
    perm, flips = ORIENT
    for key in keys:
        print('=== RUN', key, flush=True)
        ov.cmd_refine([key], (perm, flips))
        ov.cmd_overlay([key])
        ov.cmd_video([key])


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'prep':
        prep(sys.argv[2:] or list(COHORT))
    elif cmd == 'merge':
        merge()
    elif cmd in ('list', 'keys'):
        m = json.load(open(ov.MANIFEST))
        for k in m:
            print(k if cmd == 'keys' else f'{k}: {m[k]["phase"]} | {m[k]["prov"][:60]}')
    elif cmd == 'run':
        run(sys.argv[2:])
