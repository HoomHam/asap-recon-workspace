"""Scratch: header diff of every 2022 (pre-v2) spiral-dyn twix vs the measured-calibration twix and v2/v3
references — which calibration matches which sequence version (ADC/gradient/protocol parameters).

Usage: helpers/.venv/bin/python helpers/_prev2_headers.py
Writes: outputs/prev2_calib_2026-10-03/headers.csv (all scalar MeasYaps keys that vary across files)
"""
import csv
import glob
import pathlib
import warnings

import mapvbvd

warnings.filterwarnings('ignore')
IMG = '/Volumes/HoomHamExt/_5t_images_roundtrip/Images'
CAL = '/Volumes/HoomHamExt/AIdiff_FMIG_Drive_Data/calibrations'
HH = '/Volumes/HoomHamExt/Work/Images/MRI/Human'
OUT = pathlib.Path(__file__).resolve().parents[1] / 'outputs' / 'prev2_calib_2026-10-03'

FILES = sorted(glob.glob(f'{IMG}/2022-*/*/meas_*fa_spiral_dyn_fancy_2022*.dat')) \
    + [f'{IMG}/2023-11-09/000HH/meas_MID02845_FID05090_fa_spiral_dyn_fancy_20221103.dat'] \
    + sorted(glob.glob(f'{HH}/2022-10-2[67]_000HH/meas_*spiral_dyn*.dat')) \
    + sorted(glob.glob(f'{CAL}/*/meas_*_X.dat')) \
    + [f'{IMG}/2023-03-10/004LR/meas_MID00124_FID02516_fa_spiral_dyn_fancy_v2_20230131.dat',
       f'{IMG}/2024-03-06/030DN/meas_MID00543_FID08467_fa_spiral_dyn_fancy_v3_20240130.dat']


def flat(y):
    out = {}
    for k, v in y.items():
        key = '.'.join(map(str, k))
        if isinstance(v, (int, float, str)):
            out[key] = v
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for f in FILES:
        p = pathlib.Path(f)
        try:
            tw = mapvbvd.mapVBVD(f, quiet=True)
            tw = tw[-1] if isinstance(tw, list) else tw
            y = flat(tw.hdr['MeasYaps'])
            im = tw.image
            r = dict(file=p.name, session=f'{p.parent.parent.name}/{p.parent.name}', size_mb=round(p.stat().st_size / 1e6),
                     NCol=int(im.NCol), NLin=int(im.NLin), NRep=int(im.NRep), NCha=int(im.NCha), NSet=int(im.NSet), **y)
        except Exception as e:  # noqa: BLE001
            r = dict(file=p.name, session=f'{p.parent.parent.name}/{p.parent.name}', err=f'{type(e).__name__}: {e}'[:120])
        rows.append(r)
        print(r['session'], r['file'][:60], r.get('err', ''), flush=True)
    keys = list(dict.fromkeys(k for r in rows for k in r))
    fixed = {'file', 'session', 'size_mb', 'err'}
    vary = [k for k in keys if k in fixed or len({str(r.get(k)) for r in rows}) > 1]
    vary = [k for k in vary if not any(s in k for s in ('Coil', 'lRxChannel', 'Adj', 'ulVersion', 'tSequenceVariant',
                                                       'sProtConsistencyInfo', 'lSequenceID', 'Patient', 'sAutoAlign'))]
    with open(OUT / 'headers.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=vary, extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)
    print('wrote', OUT / 'headers.csv', len(rows), 'files', len(vary), 'varying keys')


if __name__ == '__main__':
    main()
