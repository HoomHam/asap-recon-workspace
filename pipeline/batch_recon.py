#!/usr/bin/env python3
"""Batch driver for dyn_recon.py over a list of (excel_row, date, id) subjects.

For each subject: run all three recons (s,p,d) via dyn_recon.py, and on success
green-highlight (FF92D050, cols B-F) that subject's row in SNR_Table_All.xlsx and
save the workbook immediately (so the batch is resumable / crash-safe).
Failures are logged and the batch continues.
"""
import subprocess, sys, os, datetime
from pathlib import Path
import openpyxl
from openpyxl.styles import PatternFill

PIPELINE = Path(__file__).resolve().parent
REPO     = PIPELINE.parent.parent
PY       = sys.executable
DYN      = PIPELINE / 'dyn_recon.py'
XLSX     = REPO / 'workspace' / 'data' / 'SNR_Table_All.xlsx'
OUT_ROOT  = '/Volumes/HoomHamExt/Dynamic'
DATA_ROOT = '/Volumes/HoomHam/Images/human'   # source raw folders (dyn_recon default)
GREEN    = PatternFill(start_color='FF92D050', end_color='FF92D050', fill_type='solid')

# (excel_row, date, id). Selection = ungreened rows whose folder has a whitelisted
# spiral-dyn seq (v2_20230131 / v3_20230821 / v3_20240130), SNR-desc. 2022 subjects
# excluded — they use old seq fa_spiral_dyn_fancy_20220812 (not in dyn_recon whitelist).
# 007IT (row 17, top SNR) skipped — GPU curve_fit non-convergence, data-dependent.
WORKLIST = [
    (37, '2023-03-09', '023DB'),
    (46, '2024-02-02', '005DS'),
    (53, '2023-03-10', '004LR'),
    (57, '2023-03-08', '017AK'),
    (77, '2023-03-24', '029CK'),
]


def mark_green(row):
    wb = openpyxl.load_workbook(XLSX)
    ws = wb['Sheet1']
    for c in range(2, 7):          # cols B..F
        ws.cell(row, c).fill = GREEN
    wb.save(XLSX)


def has_pneumotach(date, sid):
    """Mirror dyn_recon.resolve(): a pneumotach file is any name containing
    'pneumotach' (case-insensitive), excluding AppleDouble and real .dat."""
    src = Path(DATA_ROOT) / date / sid
    if not src.is_dir():
        return False
    return any('pneumotach' in f.lower() and not f.startswith('._')
               and not f.endswith('.dat') for f in os.listdir(src))


def done_already(date, sid):
    """Complete = every *runnable* method produced output.mrd. SIGNAL(s) and
    DIAPHRAGM(d) always run; PNEUMOTACH(p) only when a pneumotach file exists
    (dyn_recon skips it otherwise), so p is required only then."""
    base = Path(OUT_ROOT) / f'{date}_{sid}'
    required = ['s', 'd'] + (['p'] if has_pneumotach(date, sid) else [])
    # recon.mat alone also counts (output.mrd pruned from finished trees, 2026-10-21)
    return all((base / m / 'output.mrd').is_file() or (base / m / 'recon.mat').is_file() for m in required)


def main():
    results = []
    for i, (row, date, sid) in enumerate(WORKLIST, 1):
        tag = f'{date}_{sid}'
        print(f'\n########## [{i}/{len(WORKLIST)}] {tag} (row {row}) ##########', flush=True)
        if done_already(date, sid):
            print(f'[batch] {tag} already complete — marking green, skipping recon', flush=True)
            mark_green(row); results.append((tag, 'skipped-done')); continue
        t0 = datetime.datetime.now()
        rc = subprocess.run(
            [PY, str(DYN), '--date', date, '--id', sid, '--methods', 's,p,d', '--keep-run']
        ).returncode
        dt = (datetime.datetime.now() - t0).total_seconds() / 60
        if done_already(date, sid):          # complete regardless of rc (salvage may finish it)
            mark_green(row)
            import glob as _g
            for rd in _g.glob(str(PIPELINE / 'runs' / f'{date}_{sid}_*')):
                __import__('shutil').rmtree(rd, ignore_errors=True)
            print(f'[batch] {tag} OK ({dt:.1f} min) — row {row} greened', flush=True)
            results.append((tag, f'ok {dt:.1f}min'))
        else:
            print(f'[batch] {tag} FAILED (rc={rc}, complete={done_already(date, sid)})', flush=True)
            results.append((tag, f'FAIL rc={rc}'))
    print('\n########## BATCH SUMMARY ##########')
    for tag, st in results:
        print(f'  {tag:22s} {st}')


if __name__ == '__main__':
    main()
