#!/usr/bin/env python3
"""
Merge two (or more) free-breathing spiral-dyn twix files of one session into ONE
dynamic acquisition for the Xe-129 recon, so both doses feed one binned image.

Why this is legal for the recon (checked against raw.py / results.py, 2026-09-16):
  * time axis is index * TR (raw.py: ilvtime = arange(ntotalilvs) * TR) — MDH
    timestamps are never used, so the wall-clock gap between the two acquisitions
    simply collapses. SIGNAL and DIAPHRAGM binning are phase-per-breath (EE-to-EE)
    and read the breathing from the data itself, so a collapsed gap costs at most the
    one "breath" that straddles the seam.
  * the calibration (spectroscopy) block is the first `numspec` lines of the dynamic
    array and is used once (RBC/TP fit). It must exist ONLY at the head of the merged
    array -> file B's cal block is stripped.
  * gas/dissolved interleaves alternate by line parity (ilvperTR = 1/2/3, detected
    from the FFT of sample 0), and the spiral arm of gas interleave i is i mod
    nuniqueilvs (nuniquesmp / npts, e.g. 832 for v3, 640 for v2). Every file restarts
    the sequence at arm 0 after its cal block, so file A is TRIMMED to a whole number
    of arm cycles (loses < 1 cycle = < 14 s of the dose-out tail) before B is appended.
  * amplitude differences between doses do not matter: noise normalisation is global,
    EE detection is local-minimum based, DIAPHRAGM tracks an edge position.
  * PNEUMOTACH binning is NOT supported for merged input (the pneumotach trace is one
    wall-clock record; it would need cropping + re-stitching on the collapsed axis).
    Run s and d only.

Usage (standalone: writes input.mrd + merge_report.json + merge_k0.png):
  python merge_dyn.py --data-dir <Images/date/id> --mids MID00543,MID00545 \
      --gp-traj <gp.npy> --dp-traj <dp.npy> [--ref MID00515] --binning DIAPHRAGM \
      --out <run_dir>/input.mrd
Called by dyn_recon.py --merge.
"""
import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
import convert_siemens_to_mrd as conv          # noqa: E402  (root, read-only)
import raw as rawmod                           # noqa: E402


def _find_dat(datadir, mid):
    hits = [f for f in os.listdir(datadir)
            if f.endswith('.dat') and not f.startswith('._') and f'_{mid}_' in f]
    if len(hits) != 1:
        raise SystemExit(f'{mid}: expected exactly one .dat in {datadir}, found {hits}')
    return str(Path(datadir) / hits[0])


def _numspec_lines(meta, traj, npts_raw):
    """Replicates raw.py line 256 in the tyger path (npts = raw npts - killpts)."""
    npts = npts_raw - traj.killpts
    return int(meta['numspec_raw'] / 20 * traj.nsmpperusimg / npts + .1) if meta['numspec_raw'] else 0


def _ilv_per_tr(img):
    """raw.py lines 263-287: FFT of sample 0 over lines tells gas / gas-dp / gas-dp-dp."""
    pat = np.abs(np.fft.fft(img[0, 0, :]))
    n = len(pat)
    if max(pat[(n // 3 - 5):(n // 3 + 5)]) > pat[0] / 5:
        return 3
    if max(pat[(n // 2 - 5):(n // 2 + 5)]) > pat[0] / 3:
        return 2
    return 1


def _timestamps(dat_file):
    """MDH timestamps (s) for the report only — the recon never sees them."""
    import mapvbvd
    tw = mapvbvd.mapVBVD(dat_file, quiet=True)
    tw = tw[-1] if isinstance(tw, list) else tw
    ts = np.asarray(tw.image.timestamp).ravel().astype(float) * 2.5e-3
    return float(ts[0]), float(ts[-1])


def build_merged(datadir, mids, gp_traj, dp_traj, out_mrd, binning='DIAPHRAGM',
                 ref_mid=None, pneumotach_file=None, params=None, report_dir=None):
    datadir = Path(datadir)
    files = [_find_dat(datadir, m) for m in mids]
    traj = rawmod.traj()
    # dp trajectory is optional (v3_20230821 has gp only): go through load_traj_from_array,
    # which accepts None; load_traj_from_npfile does not
    traj.load_traj_from_array(np.load(gp_traj), np.load(dp_traj) if dp_traj else None, 32)

    parts, report = [], dict(mids=mids, files=[os.path.basename(f) for f in files], parts=[])
    meta0 = None
    for i, f in enumerate(files):
        arr, meta = conv._read_twix(f)                     # (npts, nch, nlines)
        if meta0 is None:
            meta0 = meta
        nspec = _numspec_lines(meta, traj, arr.shape[0])
        t0, t1 = _timestamps(f)
        info = dict(mid=mids[i], lines=int(arr.shape[2]), numspec=nspec,
                    mdh_t0=t0, mdh_t1=t1, dur_s=round(t1 - t0, 2))
        if i == 0:
            head, img = arr[..., :nspec], arr[..., nspec:]
            p = _ilv_per_tr(img)
            cycle = int(traj.nuniquesmp / (arr.shape[0] - traj.killpts))   # gas interleaves per arm cycle
            ngas = img.shape[2] // p
            keep_gas = (ngas // cycle) * cycle
            img = img[..., :keep_gas * p]
            info.update(ilvperTR=p, cycle_gas_ilv=cycle, gas_ilv=ngas, kept_gas_ilv=keep_gas,
                        trimmed_lines=int(arr.shape[2] - nspec - img.shape[2]),
                        trimmed_s=round((ngas - keep_gas) * meta['TR_us'] * 1e-6, 2))
            parts += [head, img]
        else:
            img = arr[..., nspec:]                          # drop B's cal block
            gap = t0 - report['parts'][-1]['mdh_t1']
            info.update(gap_before_s=round(gap, 1), stripped_cal_lines=nspec)
            parts.append(img)
            # timing must match; numspec may differ (025VP MID02743 has NO cal block at all,
            # numspec_raw 0 -> nothing to strip; the head file's block is the one used)
            for k in ('TR_us', 'TE_us', 'dtdyn_ns'):
                if meta[k] != meta0[k]:
                    raise SystemExit(f'{mids[i]}: header {k} {meta[k]} != {meta0[k]} of {mids[0]}')
            if meta['numspec_raw'] != meta0['numspec_raw']:
                info['numspec_note'] = f'numspec_raw {meta["numspec_raw"]} vs head {meta0["numspec_raw"]}'
        report['parts'].append(info)
    merged = np.concatenate(parts, axis=2)
    report['merged_lines'] = int(merged.shape[2])
    report['merged_imaging_s'] = round((merged.shape[2] - report['parts'][0]['numspec'])
                                       / report['parts'][0]['ilvperTR'] * meta0['TR_us'] * 1e-6, 1)

    ref_file = _find_dat(datadir, ref_mid) if ref_mid else None
    report['reference'] = os.path.basename(ref_file) if ref_file else None

    # hand the merged array to the root converter without touching the root code:
    # patch its file classifier + reader for this one call
    orig_classify, orig_read = conv._classify_dat_files, conv._read_twix

    def _classify(dat_dir, dat_files=None, seqname=None):
        return files[0], ref_file

    def _read(dat_file):
        if dat_file == files[0]:
            return merged, meta0
        return orig_read(dat_file)

    conv._classify_dat_files, conv._read_twix = _classify, _read
    try:
        conv.convert_siemens_to_mrd(str(datadir), str(out_mrd), gp_traj, dp_traj,
                                    pneumotach_file=pneumotach_file, killpts=2,
                                    params=dict(params or {}, binning=binning))
    finally:
        conv._classify_dat_files, conv._read_twix = orig_classify, orig_read

    report_dir = Path(report_dir) if report_dir else Path(out_mrd).parent
    report_dir.mkdir(parents=True, exist_ok=True)
    with open(report_dir / 'merge_report.json', 'w') as f:
        json.dump(report, f, indent=1)
    _k0_figure(merged, report, meta0, report_dir / 'merge_k0.png')
    return out_mrd, report


def _k0_figure(merged, report, meta, path):
    """|k0| of every gas line across the merged array, seams marked — eyeball check."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    p0 = report['parts'][0]
    p = p0['ilvperTR']
    img = merged[..., p0['numspec']:]
    gas = np.abs(img[:8, :, ::p]).sum(axis=(0, 1))
    t = np.arange(gas.size) * meta['TR_us'] * 1e-6
    fig, ax = plt.subplots(figsize=(14, 3.2), dpi=110)
    ax.plot(t, gas, lw=0.5)
    seam = p0['kept_gas_ilv']
    for part in report['parts'][1:]:
        ax.axvline(seam * meta['TR_us'] * 1e-6, color='r', ls='--', lw=1)
        ax.text(seam * meta['TR_us'] * 1e-6, gas.max(), f" {part['mid']} (gap {part['gap_before_s']} s)",
                color='r', fontsize=8, va='top')
        seam += (part['lines'] - part['stripped_cal_lines']) // p
    ax.set_xlabel('recon time axis  (index * TR, gap collapsed)  [s]')
    ax.set_ylabel('|k0| gas (sum first 8 pts)')
    ax.set_title(f"merged {' + '.join(report['mids'])}: {report['merged_imaging_s']} s imaging, "
                 f"A trimmed {p0['trimmed_s']} s to whole arm cycles ({p0['cycle_gas_ilv']} ilv)", fontsize=9)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data-dir', required=True)
    ap.add_argument('--mids', required=True, help='comma-separated, acquisition order, e.g. MID00543,MID00545')
    ap.add_argument('--gp-traj', required=True)
    ap.add_argument('--dp-traj', default=None)
    ap.add_argument('--ref', default=None, help='breath-hold MID to store as reference_acquisition')
    ap.add_argument('--pneumotach', default=None)
    ap.add_argument('--binning', default='DIAPHRAGM', choices=['SIGNAL', 'DIAPHRAGM'])
    ap.add_argument('--out', required=True, help='output input.mrd path')
    args = ap.parse_args()
    _, rep = build_merged(args.data_dir, args.mids.split(','), args.gp_traj, args.dp_traj, args.out,
                          binning=args.binning, ref_mid=args.ref, pneumotach_file=args.pneumotach)
    print(json.dumps(rep, indent=1))


if __name__ == '__main__':
    main()
