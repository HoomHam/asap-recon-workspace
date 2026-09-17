#!/usr/bin/env python
"""Convert RT 4D-CT DICOM series -> NIfTI (one file per respiratory phase).

Input : /Volumes/HoomHamExt/Work/CT/Data/<subj>/CT/<series_dir>/*.dcm
Output: /Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/ct_overlay/<subj>/ct/<tag>.nii.gz
        + index.json (tag -> series description, n slices, spacing, origin, direction)

Geometry: SimpleITK ImageSeriesReader sorts slices by ImagePositionPatient and
writes the LPS direction cosines into the NIfTI (nibabel reads it back as RAS).
Voxel values are Hounsfield units (rescale slope/intercept applied by the reader).
"""
import glob
import json
import os
import re
import sys

import SimpleITK as sitk

SRC = '/Volumes/HoomHamExt/Work/CT/Data'
DST = '/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/ct_overlay'
SUBJECTS = ['01BB', '02ZS', '03PM']


def tag_for(desc: str) -> str:
    """Short tag from SeriesDescription: 'phase00', 'phase12', ..., 'average', 'pbv'."""
    m = re.search(r'(\d+)%', desc)
    if m:
        return f'phase{int(m.group(1)):02d}'
    if 'Average' in desc:
        return 'average'
    if 'PBV' in desc:
        return 'pbv'
    return re.sub(r'\W+', '_', desc).strip('_').lower()


def convert_series(series_dir: str, out_path: str) -> dict:
    reader = sitk.ImageSeriesReader()
    ids = reader.GetGDCMSeriesIDs(series_dir)
    if not ids:
        raise RuntimeError(f'no DICOM series in {series_dir}')
    if len(ids) > 1:
        print(f'  WARNING {len(ids)} series UIDs in {series_dir}, taking first')
    files = reader.GetGDCMSeriesFileNames(series_dir, ids[0])
    reader.SetFileNames(files)
    reader.MetaDataDictionaryArrayUpdateOn()
    reader.LoadPrivateTagsOff()
    img = reader.Execute()
    img = sitk.Cast(img, sitk.sitkInt16)
    sitk.WriteImage(img, out_path, useCompression=True)
    desc = reader.GetMetaData(0, '0008|103e').strip() if reader.HasMetaDataKey(0, '0008|103e') else ''
    return {
        'series_dir': os.path.basename(series_dir),
        'description': desc,
        'n_slices': img.GetSize()[2],
        'size': list(img.GetSize()),
        'spacing': [round(s, 5) for s in img.GetSpacing()],
        'origin': [round(o, 3) for o in img.GetOrigin()],
        'direction': [round(d, 4) for d in img.GetDirection()],
        'file': os.path.basename(out_path),
    }


def main(subjects):
    for subj in subjects:
        out_dir = os.path.join(DST, subj, 'ct')
        os.makedirs(out_dir, exist_ok=True)
        index = {}
        for series_dir in sorted(glob.glob(os.path.join(SRC, subj, 'CT', '*'))):
            if not os.path.isdir(series_dir):
                continue
            first = next((f for f in os.listdir(series_dir) if not f.startswith('.')), None)
            if first is None:
                continue
            d = sitk.ReadImage(os.path.join(series_dir, first))
            desc = d.GetMetaData('0008|103e').strip() if d.HasMetaDataKey('0008|103e') else os.path.basename(series_dir)
            tag = tag_for(desc)
            out_path = os.path.join(out_dir, f'{tag}.nii.gz')
            if os.path.exists(out_path):
                print(f'{subj} {tag}: exists, skip')
                with open(os.path.join(out_dir, 'index.json')) as fh:
                    index[tag] = json.load(fh).get(tag, {'file': f'{tag}.nii.gz'})
                continue
            print(f'{subj} {tag}: converting "{desc}"')
            index[tag] = convert_series(series_dir, out_path)
            print(f'   -> {index[tag]["size"]} spacing {index[tag]["spacing"]}')
        with open(os.path.join(out_dir, 'index.json'), 'w') as fh:
            json.dump(index, fh, indent=2)


if __name__ == '__main__':
    main(sys.argv[1:] or SUBJECTS)
