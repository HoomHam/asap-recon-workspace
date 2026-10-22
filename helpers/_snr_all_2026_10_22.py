"""Scratch (2026-10-22, Hooman: "SNR for all and show me the table"): snr_calc.session_snr over every reconstructed Xe
session in all trees, DIAPHRAGM (d) binning, one row per session.

Primary recon per session = production settings (deapod): AIkill_Dynamic_b44d if present, else AIkill_Dynamic_2022
(all its d/ folders ran 6b01071 or ee3c91f = b44 + deapod), else AIkill_Dynamic (pre-b44 build, Sept). The
AIkill_Dynamic value sits next to it where both exist (same column as SNR_Table_All.xlsx I–S, 2026-09-13).
Alternates are skipped: _allcoil/_autodrop/_drop1 (coil-drop variants, canonical name kept by _pick_best_2022),
_win (s-only k0 re-cuts), _MID142 (weak twin of 000ID), _dptraj (008LL test), *_REJECTED; *_merged replaces its
single-file original. Diffusion b0/b and 006RL A/B are separate dynamics -> separate rows.
AIkill_Dynamic rows reuse outputs/snr_2026-09-13/snr_table.csv (that tree is unchanged since); the rest is computed.

writes outputs/snr_2026-10-22/{snr_folders.csv, snr_sessions.csv, snr_sessions.xlsx}
"""
import csv
import re
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import snr_calc as S

EXT = Path('/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon')
WS = Path(__file__).resolve().parents[1]
OUT = WS / 'outputs' / 'snr_2026-10-22'
OLD_CSV = WS / 'outputs' / 'snr_2026-09-13' / 'snr_table.csv'
OUTCOME = WS / 'outputs' / 'leftover_recon_2026-10-21' / 'undone43_outcome.csv'
RX = re.compile(r'^(\d{4}-\d{2}-\d{2})_([0-9A-Za-z-]+?)(?:_(.*))?$')
SKIP = ('allcoil', 'autodrop', 'drop1', 'win', 'MID142', 'dptraj', 'REJECTED')
FLAGS = {'2022-02-02_unk-MID16_merged': 'NOT A LUNG (static bright point)',
         '2022-11-04_000HH': 'noise (weak scan)',
         '2022-06-02_008LL': 'SUSPECT: dissolved-only set-up?',
         '2023-01-24_000LL': 'weak', '2022-03-10_006RL_A_merged': 'weak'}


def canonical(tree):
    """key -> (date, subject, variant) for the folders that represent a session's dynamic in this tree."""
    keys = {}
    for p in sorted((EXT / tree).iterdir()):
        m = RX.match(p.name)
        if not m or not (p / 'd' / 'recon.mat').is_file():
            continue
        date, sid, suf = m.groups()
        if suf and any(s in suf for s in SKIP):
            continue
        keys[p.name] = (date, sid, suf or '')
    for k, (date, sid, suf) in list(keys.items()):      # *_merged replaces the single-file original
        if suf == 'merged' and f'{date}_{sid}' in keys:
            del keys[f'{date}_{sid}']
    return keys


def run_one(arg):
    tree, key = arg
    S.DYN = EXT / tree
    try:
        r = S.session_snr(key, 'd')
    except Exception as e:                              # keep going; record the failure
        r = dict(key=key, binning='d', error=f'{type(e).__name__}: {e}')
    r['tree'] = tree
    return r


def spec_of(tree, key):
    f = EXT / tree / key / 'd' / 'codespec.yml'
    m = re.search(r'([0-9a-f]{7})[0-9a-f]{33}', f.read_text()) if f.is_file() else None
    return m.group(1) if m else ('pre-b44 (Sept)' if tree == 'AIkill_Dynamic' else '?')


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    trees = {t: canonical(t) for t in ('AIkill_Dynamic_b44d', 'AIkill_Dynamic_2022', 'AIkill_Dynamic')}
    old = {r['key']: r for r in csv.DictReader(open(OLD_CSV))} if OLD_CSV.is_file() else {}
    rows, todo = {}, []
    for t, keys in trees.items():
        for k in keys:
            if t == 'AIkill_Dynamic' and k in old and not old[k].get('error'):
                rows[(t, k)] = dict(old[k], tree=t, reused='2026-09-13')
            else:
                todo.append((t, k))
    print(f'{len(rows)} reused, {len(todo)} to compute', flush=True)
    with ProcessPoolExecutor(int(sys.argv[1]) if len(sys.argv) > 1 else 4) as ex:
        for i, r in enumerate(ex.map(run_one, todo), 1):
            rows[(r['tree'], r['key'])] = r
            print(f'[{i}/{len(todo)}] {r["tree"]}/{r["key"]}: '
                  + (f'ERROR {r["error"]}' if r.get('error') else
                     f'gas EI/EE {r["gas_snr_insp"]:.1f}/{r["gas_snr_exp"]:.1f}  dp {r["dp_snr_insp"]:.1f}'), flush=True)

    fields = ['tree', 'key', 'spec'] + [f for f in ('insp_bin', 'exp_bin', 'gas_snr_insp', 'gas_snr_exp', 'dp_snr_insp',
                                                    'dp_snr_exp', 'gas_snr_max', 'gas_snr_min', 'dp_snr_max', 'dp_snr_min',
                                                    'dp_method', 'lung_vox_insp', 'lung_vox_exp', 'bins', 'error')]
    fmt = lambda v: f'{v:.1f}' if isinstance(v, float) else ('' if v is None else v)
    with open(OUT / 'snr_folders.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields + ['reused'], extrasaction='ignore')
        w.writeheader()
        for (t, k), r in sorted(rows.items()):
            w.writerow({kk: fmt(v) for kk, v in dict(r, spec=spec_of(t, k)).items()})

    notes = {}
    if OUTCOME.is_file():
        for r in csv.DictReader(open(OUTCOME)):
            for k in re.split(r'[ +/]+', r['output_key(s)']):
                if k.startswith('20'):
                    notes[k.split('/')[-1]] = r['quality_note']
    sess = {}
    for t in ('AIkill_Dynamic', 'AIkill_Dynamic_2022', 'AIkill_Dynamic_b44d'):   # later trees win -> production last
        for k, (date, sid, suf) in trees[t].items():
            variant = suf.replace('merged', '').strip('_') if suf else ''
            sk = (date, sid, variant)
            ent = sess.setdefault(sk, {})
            if t == 'AIkill_Dynamic':
                ent['old'] = rows.get((t, k))
            ent.update(primary=rows.get((t, k)), ptree=t, pkey=k)
    out = []
    for (date, sid, variant), e in sorted(sess.items()):
        p, o = e['primary'] or {}, e.get('old') or {}
        g = lambda r, f: (fmt(float(r[f])) if r.get(f) not in (None, '', 'nan') and str(r.get(f)) != 'nan' else '')
        out.append(dict(date=date, subject=sid, variant=variant, recon=e['ptree'].replace('AIkill_Dynamic', 'AIkill'),
                        spec=spec_of(e['ptree'], e['pkey']), merged='y' if 'merged' in e['pkey'] else '',
                        gas_EI=g(p, 'gas_snr_insp'), gas_EE=g(p, 'gas_snr_exp'), gas_max=g(p, 'gas_snr_max'),
                        gas_min=g(p, 'gas_snr_min'), dp_EI=g(p, 'dp_snr_insp'), dp_EE=g(p, 'dp_snr_exp'),
                        dp=('none' if 'none' in str(p.get('dp_method', '')) else
                            'split' if str(p.get('dp_method', '')).startswith('split') else
                            'unsplit' if p.get('dp_method') else ''),
                        lung_vox_EI=p.get('lung_vox_insp', ''), gas_EI_sept=g(o, 'gas_snr_insp') if e['ptree'] != 'AIkill_Dynamic' else '',
                        dp_EI_sept=g(o, 'dp_snr_insp') if e['ptree'] != 'AIkill_Dynamic' else '',
                        note=FLAGS.get(e['pkey']) or ('' if notes.get(e['pkey'], 'clean') in ('clean', 'ok')
                                                      else notes[e['pkey']]),
                        error=p.get('error', '')))
    cols = list(out[0])
    with open(OUT / 'snr_sessions.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(out)
    import openpyxl
    from openpyxl.styles import Font
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = 'SNR 2026-10-22'
    ws.append(cols)
    for c in ws[1]:
        c.font = Font(bold=True)
    for r in out:
        ws.append([float(r[c]) if c in ('gas_EI', 'gas_EE', 'gas_max', 'gas_min', 'dp_EI', 'dp_EE', 'gas_EI_sept',
                                        'dp_EI_sept') and r[c] != '' else r[c] for c in cols])
    ws.freeze_panes = 'A2'
    wb.save(OUT / 'snr_sessions.xlsx')
    print(f'{len(out)} session rows -> {OUT}/snr_sessions.csv + .xlsx; folders -> snr_folders.csv')


if __name__ == '__main__':
    main()
