#!/usr/bin/env python3
"""
ASAP recon pipeline — one command to run Steve's reconstruction on a dataset
through the Tyger cloud GPU.

    python asap_run.py 025JC              # GUI param picker, then full run
    python asap_run.py 025JC --no-gui     # use gvar defaults, no GUI
    python asap_run.py 025JC --params p.yml
    python asap_run.py 025JC --slice 50   # montage slice index for plot_recon

Stages (all reuse existing root-repo code; nothing in the root repo is modified):
  1. resolve  — find the dataset dir + its gas/dissolved trajectory (mirrors
                main.py ID_callback)
  2. params   — param_gui.py (default) / --params file / gvar defaults
  3. convert  — convert_siemens_to_mrd(): .dat -> input.mrd with params baked
                into the MRD header
  4. submit   — tyger run exec: stream input.mrd to the cloud GPU job, stream
                output.mrd back (this is the cloud path; main.py instead runs
                tyger_recon.py locally, which needs CUDA and fails on the Mac)
  5. publish  — plot_recon.py montages, copied with the run artifacts into
                workspace/outputs/<dataset>/<timestamp>/

The root repo is read-only here: NEVER git commit / push from the repo root.
"""
import argparse
import datetime as _dt
import json
import os
import re
import shutil
import subprocess
import sys
from glob import glob
from pathlib import Path

# --- repo layout -------------------------------------------------------------
PIPELINE_DIR = Path(__file__).resolve().parent          # workspace/pipeline
WORKSPACE    = PIPELINE_DIR.parent                       # workspace
REPO_ROOT    = WORKSPACE.parent                          # repo root (read-only)
DATA_BASE    = WORKSPACE / 'data' / 'xe'                 # == gvar.basefolder
TYGER_SPEC   = PIPELINE_DIR / 'recon_codespec.yml'  # fork image w/ DIAPHRAGM binning
PLOT_SCRIPT  = REPO_ROOT / 'tyger_deploy' / 'plot_recon.py'
OUTPUTS_DIR  = WORKSPACE / 'outputs'
RUNS_DIR     = PIPELINE_DIR / 'runs'
PARAM_GUI    = PIPELINE_DIR / 'param_gui.py'
TYGER        = os.path.expanduser('~/bin/tyger')
LOGIN_FILE   = '/Users/hoomham/Hooman/Work/Spinhance/Tyger/LOGIN_FILE.yml'

# root repo on path for convert_siemens_to_mrd + gtypes
sys.path.insert(0, str(REPO_ROOT))


def _die(msg, code=1):
    print(f'\n[asap_run] ERROR: {msg}', file=sys.stderr)
    sys.exit(code)


# --- stage 1: resolve dataset ------------------------------------------------
def resolve_dataset(token, datatype='human'):
    """Mirror main.py ID_callback: locate the subject dir, its dynamic .dat
    file(s), the matching gas/dissolved trajectory .npy, and any pneumotach."""
    base = DATA_BASE / datatype
    if not base.is_dir():
        _die(f'data base not found: {base}')

    matches = sorted(p for p in glob(str(base / '*' / f'*{token}*')) if os.path.isdir(p))
    if not matches:
        _die(f'no subject dir matching "{token}" under {base}')
    if len(matches) > 1:
        lst = '\n  '.join(matches)
        _die(f'"{token}" is ambiguous ({len(matches)} matches):\n  {lst}\n'
             f'Re-run with a more specific token (e.g. include the date).')
    datadir = Path(matches[0])

    trajdir = base / 'traj'
    seqfile = trajdir / 'seqnames.txt'
    if not seqfile.is_file():
        _die(f'seqnames.txt not found: {seqfile}')
    seqs = []
    for line in seqfile.read_text().splitlines():
        parts = line.split()
        if len(parts) >= 2:
            seqs.append((parts[0], parts[1]))  # (fileformat, seqname)

    rawfiles = [f for f in os.listdir(datadir) if not f.startswith('.')]
    seqname, dynfiles = None, []
    for _fmt, seq in seqs:
        ds = [str(datadir / f) for f in rawfiles
              if seq in f and f.endswith('.dat')]
        if ds:
            seqname, dynfiles = seq, sorted(ds, key=os.path.getsize)
            break
    if not seqname:
        # fall back to Bruker-style rawdata.job0 naming
        ds = [str(datadir / f) for f in rawfiles if 'rawdata.job0' in f]
        if ds:
            dynfiles = sorted(ds, key=os.path.getsize)
        else:
            _die(f'no recognised raw data in {datadir} '
                 f'(looked for seqnames {[s for _, s in seqs]} or rawdata.job0)')

    gp = trajdir / f'{seqname}_gp.npy' if seqname else None
    dp = trajdir / f'{seqname}_dp.npy' if seqname else None
    gp = str(gp) if gp and gp.is_file() else None
    dp = str(dp) if dp and dp.is_file() else None
    if gp is None:
        _die(f'gas-phase trajectory missing for seqname "{seqname}" in {trajdir}')
    if dp is None:
        print(f'[asap_run] WARNING: no dissolved-phase trajectory for "{seqname}"'
              f' — reconstructing gas-phase only.', file=sys.stderr)

    pneumo = next((str(datadir / f) for f in rawfiles if 'pneumotach' in f), None)

    print(f'[asap_run] dataset   : {datadir}')
    print(f'[asap_run] seqname    : {seqname}')
    print(f'[asap_run] dyn .dat   : {[os.path.basename(d) for d in dynfiles]}')
    print(f'[asap_run] gp / dp    : {os.path.basename(gp)} / '
          f'{os.path.basename(dp) if dp else "(none)"}')
    print(f'[asap_run] pneumotach : {os.path.basename(pneumo) if pneumo else "(none)"}')
    return dict(datadir=str(datadir), dynfiles=dynfiles, seqname=seqname,
                gp=gp, dp=dp, pneumo=pneumo)


# --- stage 2: params ---------------------------------------------------------
PARAM_KEYS = ('MS', 'IS', 'nbins', 'griddx', 'bindt', 'gplb', 'dplb',
              'freqfilter', 'binning')


def _gvar_defaults():
    from gtypes import gvar
    g = gvar()
    return dict(MS=g.MS, IS=g.IS, nbins=g.nbins, griddx=g.griddx, bindt=g.bindt,
                gplb=g.gplb, dplb=g.dplb, freqfilter=g.freqfilter, binning='SIGNAL')


def _load_param_file(path):
    p = Path(path)
    if not p.is_file():
        _die(f'params file not found: {p}')
    text = p.read_text()
    if p.suffix.lower() in ('.yml', '.yaml'):
        try:
            import yaml
            return yaml.safe_load(text) or {}
        except ImportError:
            _die('pyyaml not installed — give a .json params file or `pip install pyyaml`')
    return json.loads(text)


def get_params(args, run_dir):
    defaults = _gvar_defaults()
    if args.params:
        params = {**defaults, **_load_param_file(args.params)}
        print(f'[asap_run] params from {args.params}')
        return params
    if args.no_gui:
        print('[asap_run] params: gvar defaults (--no-gui)')
        return defaults

    # launch the standalone GUI; it writes params.json and exits 0 on "Run"
    out = run_dir / 'params.json'
    preset = run_dir / '_defaults.json'
    preset.write_text(json.dumps(defaults))
    cmd = [sys.executable, str(PARAM_GUI), '--out', str(out),
           '--preset', str(preset), '--dataset', args.dataset]
    print('[asap_run] opening param GUI — set values, click "Run recon"...')
    rc = subprocess.run(cmd).returncode
    if rc != 0 or not out.is_file():
        _die('param GUI cancelled — aborting run.')
    params = {**defaults, **json.loads(out.read_text())}
    return params


# --- stage 3: convert --------------------------------------------------------
def convert(ds, params, run_dir):
    from convert_siemens_to_mrd import convert_siemens_to_mrd
    input_mrd = run_dir / 'input.mrd'
    print(f'[asap_run] converting -> {input_mrd}')
    convert_siemens_to_mrd(
        ds['datadir'], str(input_mrd),
        gp_traj_file=ds['gp'], dp_traj_file=ds['dp'],
        pneumotach_file=ds['pneumo'], params=params,
        dat_files=ds['dynfiles'])
    if not input_mrd.is_file() or input_mrd.stat().st_size == 0:
        _die('conversion produced no MRD output')
    return input_mrd


# --- stage 4: submit to tyger ------------------------------------------------
def _check_tyger_login():
    if not os.path.exists(TYGER):
        _die(f'tyger CLI not found at {TYGER}')
    r = subprocess.run([TYGER, 'login', 'status'], capture_output=True, text=True)
    if r.returncode != 0 or 'logged in' not in (r.stdout + r.stderr).lower():
        _die(f'not logged in to Tyger. Run:  {TYGER} login {LOGIN_FILE}')
    print(f'[asap_run] tyger: {r.stdout.strip() or r.stderr.strip()}')


def submit_to_tyger(input_mrd, run_dir):
    _check_tyger_login()
    output_mrd = run_dir / 'output.mrd'
    log_path = run_dir / 'tyger.log'
    print('[asap_run] submitting to Tyger (cloud GPU) — streaming...')
    with open(input_mrd, 'rb') as fin, open(output_mrd, 'wb') as fout, \
            open(log_path, 'w') as flog:
        proc = subprocess.run(
            [TYGER, 'run', 'exec', '-f', str(TYGER_SPEC), '--logs'],
            stdin=fin, stdout=fout, stderr=flog)
    log_text = log_path.read_text()
    m = re.search(r'run (\d+)', log_text)
    run_id = m.group(1) if m else None
    if proc.returncode != 0:
        tail = '\n'.join(log_text.splitlines()[-15:])
        _die(f'tyger run failed (exit {proc.returncode}). Last log lines:\n{tail}\n'
             f'If the job finished but download timed out, the output buffer lives'
             f' ~1 hr: {TYGER} buffer read <outputBufferId> -o {output_mrd}')
    if not output_mrd.is_file() or output_mrd.stat().st_size == 0:
        _die(f'tyger returned empty output.mrd (run {run_id}). '
             f'Recover within 1 hr: {TYGER} buffer read <outputBufferId> -o {output_mrd}')
    print(f'[asap_run] tyger run {run_id} done -> {output_mrd} '
          f'({output_mrd.stat().st_size} bytes)')
    return output_mrd, run_id


# --- stage 5: publish --------------------------------------------------------
def publish(output_mrd, params, run_id, ds, args, run_dir):
    plot_cmd = [sys.executable, str(PLOT_SCRIPT), str(output_mrd)]
    if args.slice is not None:
        plot_cmd += ['--slice', str(args.slice)]
    print('[asap_run] rendering montages (plot_recon.py)...')
    subprocess.run(plot_cmd, cwd=run_dir)  # writes *_gp.png / *_dp.png beside output.mrd

    ts = _dt.datetime.now().strftime('%Y%m%d-%H%M%S')
    out = OUTPUTS_DIR / args.dataset / ts
    out.mkdir(parents=True, exist_ok=True)
    for png in run_dir.glob('*.png'):
        shutil.copy2(png, out / png.name)
    shutil.copy2(output_mrd, out / 'output.mrd')
    (out / 'params.json').write_text(json.dumps(params, indent=2))
    if (run_dir / 'tyger.log').is_file():
        shutil.copy2(run_dir / 'tyger.log', out / 'tyger.log')
    manifest = dict(dataset=args.dataset, datadir=ds['datadir'], seqname=ds['seqname'],
                    tyger_run_id=run_id, params=params, timestamp=ts)
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2))

    # stage 6: post-process — .mat + fig/ GIFs
    POST_SCRIPT = PIPELINE_DIR / 'post_process.py'
    print('[asap_run] post-processing (mat + GIFs)...')
    post_cmd = [sys.executable, str(POST_SCRIPT), str(out),
                '--input-mrd', str(run_dir / 'input.mrd')]
    subprocess.run(post_cmd)

    print(f'\n[asap_run] DONE. Results -> {out}')
    return out


# --- main --------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Run Steve's ASAP recon on Tyger.")
    ap.add_argument('dataset', help='subject token, e.g. 025JC or 25JC')
    ap.add_argument('--datatype', default='human', help='species folder (default human)')
    ap.add_argument('--no-gui', action='store_true', help='use gvar defaults, skip GUI')
    ap.add_argument('--params', metavar='FILE', help='params .json/.yml (skips GUI)')
    ap.add_argument('--slice', type=int, default=None, help='montage slice index')
    args = ap.parse_args()

    ds = resolve_dataset(args.dataset, args.datatype)

    ts = _dt.datetime.now().strftime('%Y%m%d-%H%M%S')
    run_dir = RUNS_DIR / f'{args.dataset}_{ts}'
    run_dir.mkdir(parents=True, exist_ok=True)
    print(f'[asap_run] run dir   : {run_dir}')

    params = get_params(args, run_dir)
    print(f'[asap_run] params    : {params}')
    input_mrd = convert(ds, params, run_dir)
    output_mrd, run_id = submit_to_tyger(input_mrd, run_dir)
    publish(output_mrd, params, run_id, ds, args, run_dir)


if __name__ == '__main__':
    main()
