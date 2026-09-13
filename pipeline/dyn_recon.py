#!/usr/bin/env python3
"""
Dynamic Xe-129 recon automation — one command from a raw-data folder to the
dynamic image + all figures, dropped in an output folder named <date>_<id>.

Two ways to point it at the data:

    # resolve <data-root>/<date>/<id>  ->  <out-root>/<date>_<id>
    python dyn_recon.py --date 2024-09-10 --id 008TP

    # explicit folders
    python dyn_recon.py --data-dir /Volumes/HoomHamExt/_5t_images_roundtrip/Images/2024-09-10/008TP \
                        --out-dir  /Volumes/HoomHamExt/AIkill_Dynamic/2024-09-10_008TP

What it does (mirrors the manual run, nothing in the root repo is modified):
  1. resolve    find the raw folder, the *largest* dynamic spiral .dat, its
                gas/dissolved trajectory, and the pneumotach file (if present)
  2. convert    convert_siemens_to_mrd -> input.mrd (params baked into header)
  3. submit     buffer write -> tyger run create -> poll -> buffer read -o -p 4
                (no stdout streaming: the link is ~1-3 Mbps and `run exec --logs`
                dies on 190 MB outputs; every transfer step retries)
  4. plot       plot_recon.py   -> output_gp.png (+ output_dp.png)
  5. post       post_process.py -> recon.mat, signal_pneumo.npz, fig/*.gif, resp_traces.png
  6. publish    copy every artifact into <out-dir>

Root repo stays read-only: this script only *calls* the root scripts; it never
edits or commits there.  It lives in workspace/ where git is allowed.
"""
import argparse
import datetime as _dt
import json
import os
import shutil
import subprocess
import sys
import time
from glob import glob
from pathlib import Path

# --- repo layout -------------------------------------------------------------
PIPELINE_DIR = Path(__file__).resolve().parent            # workspace/pipeline
WORKSPACE    = PIPELINE_DIR.parent                        # workspace
REPO_ROOT    = WORKSPACE.parent                           # repo root (read-only)

CONVERT      = REPO_ROOT / 'convert_siemens_to_mrd.py'
PLOT_SCRIPT  = REPO_ROOT / 'tyger_deploy' / 'plot_recon.py'
POST_SCRIPT  = PIPELINE_DIR / 'post_process.py'
TYGER_SPEC   = PIPELINE_DIR / 'recon_codespec.yml'        # fork image w/ DIAPHRAGM
TRAJ_DIR     = WORKSPACE / 'data' / 'xe' / 'human' / 'traj'
# scratch run dirs hold ~1 GB each (input.mrd + output.mrd) -> Ext, not the laptop
# (big-output law, mirror path of this project)
RUNS_DIR     = Path('/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/pipeline_runs')

TYGER        = os.path.expanduser('~/bin/tyger')
LOGIN_FILE   = '/Users/hoomham/Hooman/Work/Spinhance/Tyger/LOGIN_FILE.yml'
# NOTE: never pass --ttl to `tyger buffer create` — this server answers 500 (2026-09-12);
# default-TTL buffers live for days, enough for a slow download + salvage
TRIES       = 4              # retries per transfer step
POLL_S       = 30             # run status poll interval

# respiratory binning methods -> output subfolder letter
BINNING_LETTER = {'SIGNAL': 's', 'PNEUMOTACH': 'p', 'DIAPHRAGM': 'd'}
LETTER_BINNING = {v: k for k, v in BINNING_LETTER.items()}

# defaults for the --date/--id resolve mode (HoomHam died; data now on Ext)
DEFAULT_DATA_ROOT = '/Volumes/HoomHamExt/_5t_images_roundtrip/Images'
DEFAULT_OUT_ROOT  = '/Volumes/HoomHamExt/AIkill_Dynamic'

PY = sys.executable  # same interpreter that launched this script


def _die(msg, code=1):
    print(f'\n[dyn_recon] ERROR: {msg}', file=sys.stderr)
    sys.exit(code)


def _log(msg):
    stamp = _dt.datetime.now().strftime('%H:%M:%S')
    print(f'[dyn_recon {stamp}] {msg}', flush=True)


# --- stage 1: resolve --------------------------------------------------------
def _known_seqnames():
    """Seqnames we have trajectories for (from the traj seqnames.txt)."""
    seqfile = TRAJ_DIR / 'seqnames.txt'
    if not seqfile.is_file():
        _die(f'seqnames.txt not found: {seqfile}')
    seqs = []
    for line in seqfile.read_text().splitlines():
        parts = line.split()
        if len(parts) >= 2:
            seqs.append(parts[1])       # (fileformat, seqname)
    return seqs


def _real_dats(datadir):
    """.dat files, skipping macOS AppleDouble '._' shadow files."""
    return [f for f in os.listdir(datadir)
            if f.endswith('.dat') and not f.startswith('._')]


def resolve(args):
    # -- data dir --
    if args.data_dir:
        datadir = Path(args.data_dir)
    else:
        if not (args.date and args.id):
            _die('give either --data-dir or both --date and --id')
        datadir = Path(args.data_root) / args.date / args.id
    if not datadir.is_dir():
        _die(f'data folder not found: {datadir}')

    # -- date/id (for the output folder name) --
    date = args.date
    subj_id = args.id
    if not (date and subj_id):
        # infer <date>/<id> from the tail of an explicit --data-dir
        parts = datadir.resolve().parts
        subj_id = subj_id or parts[-1]
        date = date or (parts[-2] if len(parts) >= 2 else 'unknown')

    # -- dynamic .dat: largest spiral-dyn file --
    seqnames = _known_seqnames()
    dats = _real_dats(datadir)
    if not dats:
        _die(f'no .dat files in {datadir}')
    cand = [f for f in dats if any(s in f for s in seqnames)]
    if not cand:
        _die(f'no spiral-dyn .dat (known seqnames {seqnames}) in {datadir}\n'
             f'  found: {dats}')
    dyn = max(cand, key=lambda f: os.path.getsize(datadir / f))
    seqname = args.seqname or next(s for s in seqnames if s in dyn)

    # -- trajectory --
    gp = TRAJ_DIR / f'{seqname}_gp.npy'
    dp = TRAJ_DIR / f'{seqname}_dp.npy'
    gp = str(gp) if gp.is_file() else None
    dp = str(dp) if dp.is_file() else None
    if gp is None:
        _die(f'gas-phase trajectory missing: {TRAJ_DIR}/{seqname}_gp.npy')

    # -- pneumotach (optional) --
    # Most files are named `pneumotach_<id>_...`, but some are `<id>_..._pneumotach`
    # (007RA) or have typo/duplicate suffixes (008CR `_pneumotache`, `_2_pneumotach`).
    # Match 'pneumotach' anywhere (case-insensitive); skip AppleDouble + real .dat.
    pneumo = next((f for f in sorted(os.listdir(datadir))
                   if 'pneumotach' in f.lower()
                   and not f.startswith('._') and not f.endswith('.dat')), None)
    pneumo = str(datadir / pneumo) if pneumo else None

    _log(f'data dir   : {datadir}')
    _log(f'date / id  : {date} / {subj_id}')
    _log(f'dynamic    : {dyn} ({os.path.getsize(datadir / dyn)} bytes)')
    _log(f'seqname    : {seqname}')
    _log(f'gp / dp    : {os.path.basename(gp)} / {os.path.basename(dp) if dp else "(none)"}')
    _log(f'pneumotach : {os.path.basename(pneumo) if pneumo else "(none)"}')

    # -- output base dir (per-method results go in <base>/{s,p,d}) --
    if args.out_dir:
        outbase = Path(args.out_dir)
    else:
        outbase = Path(args.out_root) / f'{date}_{subj_id}'
    _log(f'output base: {outbase}  (methods -> s/p/d subfolders)')

    return dict(datadir=str(datadir), seqname=seqname, gp=gp, dp=dp,
                pneumo=pneumo, outbase=outbase, date=date, id=subj_id)


# --- stage 2: convert --------------------------------------------------------
def _clean_datadir(datadir, run_dir):
    """convert_siemens_to_mrd globs the -i dir for .dat and picks the 2nd-largest
    as the breath-hold reference. macOS AppleDouble '._*.dat' shadows (4 KB) slip
    into that sort and get mistaken for the reference -> mapVBVD IndexError.
    Stage a symlink dir holding only real .dat files so convert never sees them."""
    real = [f for f in os.listdir(datadir)
            if f.endswith('.dat') and not f.startswith('._')]
    junk = [f for f in os.listdir(datadir)
            if f.endswith('.dat') and f.startswith('._')]
    if not junk:
        return str(datadir)               # nothing to hide; use dir as-is
    clean = run_dir / 'src'
    clean.mkdir(parents=True, exist_ok=True)
    for f in real:
        link = clean / f
        if not link.exists():
            os.symlink(Path(datadir) / f, link)
    _log(f'staged {len(real)} real .dat (hid {len(junk)} AppleDouble) -> {clean}')
    return str(clean)


def convert(ds, run_dir, binning):
    input_mrd = run_dir / 'input.mrd'
    cmd = [PY, str(CONVERT),
           '-i', _clean_datadir(ds['datadir'], run_dir), '-o', str(input_mrd),
           '--seqname', ds['seqname'],
           '--gp-traj', ds['gp'],
           '--binning', binning]
    if ds['dp']:
        cmd += ['--dp-traj', ds['dp']]
    if ds['pneumo']:
        cmd += ['--pneumotach', ds['pneumo']]
    _log(f'converting ({binning}) -> {input_mrd}')
    if subprocess.run(cmd).returncode != 0:
        _die('conversion failed')
    if not input_mrd.is_file() or input_mrd.stat().st_size == 0:
        _die('conversion produced no MRD output')
    return input_mrd


# --- stage 3: submit ---------------------------------------------------------
def _check_login():
    if not os.path.exists(TYGER):
        _die(f'tyger CLI not found at {TYGER}')
    r = subprocess.run([TYGER, 'login', 'status'], capture_output=True, text=True)
    if r.returncode != 0 or 'logged in' not in (r.stdout + r.stderr).lower():
        _die(f'not logged in to Tyger. Run:  {TYGER} login -f {LOGIN_FILE}')


def _tyger(*argv):
    return subprocess.run([TYGER, *argv], capture_output=True, text=True)


def _retry(step, what):
    """Run step() until it returns a truthy value; back off between tries."""
    for i in range(1, TRIES + 1):
        try:
            out = step()
            if out:
                return out
        except RuntimeError as e:
            _log(f'{what}: try {i}/{TRIES} failed: {e}')
        if i < TRIES:
            time.sleep(30 * i)
    _die(f'{what} failed after {TRIES} tries')


def _checked(r, what):
    if r.returncode != 0:
        raise RuntimeError(f'{what} exit {r.returncode}: {(r.stderr or r.stdout).strip()[-400:]}')
    return r.stdout.strip()


def _upload(input_mrd):
    """Fresh buffer per try: a half-written buffer can't be rewritten."""
    buf = _checked(_tyger('buffer', 'create'), 'buffer create')
    _checked(_tyger('buffer', 'write', buf, '-i', str(input_mrd)), 'buffer write')
    return buf


def _download(buf, output_mrd):
    if output_mrd.exists():
        output_mrd.unlink()
    _checked(_tyger('buffer', 'read', buf, '-o', str(output_mrd), '-p', '4'), 'buffer read')
    if not output_mrd.is_file() or output_mrd.stat().st_size == 0:
        raise RuntimeError('buffer read produced an empty output.mrd')
    return output_mrd


def submit(input_mrd, run_dir, spec):
    _check_login()
    output_mrd = run_dir / 'output.mrd'
    log_path = run_dir / 'tyger.log'
    shutil.copy2(spec, run_dir / 'codespec.yml')          # provenance: image sha

    mb = input_mrd.stat().st_size / 1e6
    t0 = time.time()
    _log(f'uploading input.mrd ({mb:.0f} MB) ...')
    in_buf = _retry(lambda: _upload(input_mrd), 'upload')
    _log(f'uploaded in {time.time() - t0:.0f} s -> input buffer {in_buf}')
    out_buf = _retry(lambda: _checked(_tyger('buffer', 'create'),
                                      'buffer create'), 'output buffer create')

    run_id = _retry(lambda: _checked(_tyger('run', 'create', '-f', str(spec),
                                            '-b', f'input={in_buf}',
                                            '-b', f'output={out_buf}'), 'run create'),
                    'run create')
    # buffer ids first, so a killed session can salvage with `tyger buffer read`
    with open(log_path, 'w') as flog:
        flog.write(f'# run={run_id} image_spec={spec.name} '
                   f'input_buffer={in_buf} output_buffer={out_buf}\n')
    _log(f'run {run_id} created (output buffer {out_buf}); polling every {POLL_S} s')

    status, last, t0 = None, None, time.time()
    while status not in ('Succeeded', 'Failed', 'Canceled', 'TimedOut'):
        time.sleep(POLL_S)
        r = _tyger('run', 'show', run_id)
        if r.returncode != 0:
            continue                                  # transient API blip; keep polling
        status = json.loads(r.stdout).get('status')
        if status != last:
            _log(f'run {run_id}: {status} ({time.time() - t0:.0f} s)')
            last = status

    logs = _tyger('run', 'logs', run_id)
    with open(log_path, 'a') as flog:
        flog.write(logs.stdout)
        flog.write(logs.stderr)
        flog.write(f'# final status={status}\n')
    if status != 'Succeeded':
        tail = '\n'.join(log_path.read_text().splitlines()[-15:])
        _die(f'tyger run {run_id} {status}. Last log:\n{tail}')

    t0 = time.time()
    _log(f'downloading output buffer {out_buf} ...')
    _retry(lambda: _download(out_buf, output_mrd), 'download')
    _log(f'recon done -> output.mrd ({output_mrd.stat().st_size / 1e6:.0f} MB, '
         f'{time.time() - t0:.0f} s download)')
    return output_mrd


# --- stage 4+5: plot + post-process + publish --------------------------------
def publish(output_mrd, run_dir, target_dir, args):
    _log('rendering montages (plot_recon.py)...')
    plot_cmd = [PY, str(PLOT_SCRIPT), str(output_mrd)]
    if args.slice is not None:
        plot_cmd += ['--slice', str(args.slice)]
    # plot_recon.py ends with plt.show(); on the macosx backend that blocks until the
    # window is closed, which stalls an unattended batch. Agg: files only, show() no-op.
    headless = {**os.environ, 'MPLBACKEND': 'Agg'}
    subprocess.run(plot_cmd, cwd=run_dir, env=headless)  # writes *_gp.png / *_dp.png beside output.mrd

    target_dir.mkdir(parents=True, exist_ok=True)
    for name in ('output.mrd', 'input.mrd', 'tyger.log', 'codespec.yml'):
        src = run_dir / name
        if src.is_file():
            shutil.copy2(src, target_dir / name)
    for png in run_dir.glob('*.png'):
        shutil.copy2(png, target_dir / png.name)

    _log('post-processing (mat + GIFs)...')
    subprocess.run([PY, str(POST_SCRIPT), str(target_dir),
                    '--input-mrd', str(target_dir / 'input.mrd')])

    _log(f'DONE. Results -> {target_dir}')
    return target_dir


# --- main --------------------------------------------------------------------
def main():
    global RUNS_DIR
    ap = argparse.ArgumentParser(
        description='Automated dynamic Xe-129 recon: raw folder -> images + figures.')
    ap.add_argument('--date', help='study date, e.g. 2024-09-10 (resolve mode)')
    ap.add_argument('--id',   help='subject id, e.g. 008TP (resolve mode)')
    ap.add_argument('--data-root', default=DEFAULT_DATA_ROOT,
                    help=f'root that holds <date>/<id> (default {DEFAULT_DATA_ROOT})')
    ap.add_argument('--out-root',  default=DEFAULT_OUT_ROOT,
                    help=f'root for <date>_<id> output (default {DEFAULT_OUT_ROOT})')
    ap.add_argument('--data-dir', help='explicit raw folder (overrides --date/--id)')
    ap.add_argument('--out-dir',  help='explicit output folder (overrides --out-root)')
    ap.add_argument('--methods', default='s,p,d',
                    help="binning methods to run, as letters s(ignal)/p(neumotach)/"
                         "d(iaphragm), comma-separated (default 's,p,d' = all three)")
    ap.add_argument('--seqname', help='override auto-detected trajectory seqname')
    ap.add_argument('--codespec', default=str(TYGER_SPEC),
                    help=f'Tyger run spec (pins the image sha; default {TYGER_SPEC.name})')
    ap.add_argument('--runs-dir', default=str(RUNS_DIR),
                    help=f'scratch run dirs (default {RUNS_DIR})')
    ap.add_argument('--slice', type=int, default=None, help='montage slice index')
    ap.add_argument('--keep-run', action='store_true',
                    help='keep the scratch run dirs (default: delete after publish)')
    ap.add_argument('--force', action='store_true',
                    help='re-run methods even if their output.mrd already exists')
    args = ap.parse_args()
    RUNS_DIR = Path(args.runs_dir)
    spec = Path(args.codespec).resolve()
    if not spec.is_file():
        _die(f'codespec not found: {spec}')

    # parse & validate requested methods
    letters = [x.strip().lower() for x in args.methods.split(',') if x.strip()]
    bad = [x for x in letters if x not in LETTER_BINNING]
    if bad:
        _die(f'unknown method letter(s) {bad} — use s/p/d')
    binnings = [LETTER_BINNING[x] for x in letters]

    ds = resolve(args)

    # pneumotach binning needs a pneumotach file — drop it (with a warning) if absent
    if 'PNEUMOTACH' in binnings and not ds['pneumo']:
        _log('WARNING: no pneumotach file — skipping the PNEUMOTACH (p) recon')
        binnings = [b for b in binnings if b != 'PNEUMOTACH']

    _check_login()  # fail fast before any conversion
    _log(f'codespec   : {spec}')
    _log(f'running methods: {", ".join(binnings)}')

    for binning in binnings:
        letter = BINNING_LETTER[binning]
        target = ds['outbase'] / letter
        if (target / 'output.mrd').is_file() and not args.force:
            _log(f'==================== {binning} -> {target}  (already done, skip) ====================')
            continue

        # salvage: a prior (killed) attempt may have left a *finished* recon in a run
        # dir — publish that instead of paying for another Tyger recon. Require a sane
        # size (normal output.mrd is ~190-205 MB); 0-byte / truncated files from a
        # download killed mid-stream must NOT be trusted.
        MIN_MRD = 100 * 1024 * 1024
        prior = sorted(RUNS_DIR.glob(f'{ds["date"]}_{ds["id"]}_{letter}_*'))
        prior = [p for p in prior
                 if (p / 'output.mrd').is_file() and (p / 'output.mrd').stat().st_size > MIN_MRD]
        if prior and not args.force:
            run_dir = prior[-1]
            _log(f'==================== {binning} -> {target}  (salvage {run_dir.name}) ====================')
            publish(run_dir / 'output.mrd', run_dir, target, args)
            continue

        _log(f'==================== {binning} -> {target} ====================')
        ts = _dt.datetime.now().strftime('%Y%m%d-%H%M%S')
        run_dir = RUNS_DIR / f'{ds["date"]}_{ds["id"]}_{letter}_{ts}'
        run_dir.mkdir(parents=True, exist_ok=True)
        _log(f'run dir    : {run_dir}')

        input_mrd = convert(ds, run_dir, binning)
        output_mrd = submit(input_mrd, run_dir, spec)
        publish(output_mrd, run_dir, target, args)

        if not args.keep_run:
            shutil.rmtree(run_dir, ignore_errors=True)

    _log(f'ALL DONE -> {ds["outbase"]}  ({", ".join(BINNING_LETTER[b] for b in binnings)})')


if __name__ == '__main__':
    main()
