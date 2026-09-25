#!/usr/bin/env python
"""Gas + RBC + tissue/membrane (TM) atlas from the XeCS-fit re-split maps.

Replaces the dp pages of dis_atlas.py (SUSPECT, F47): dissolved maps here come from
<session>/d/recon_resplit_xecs.mat (aRBC, aTP = k0-basis split with the XeCS specfit
angle, LUMPED scale: TM = |mem1+mem2|; divide RBC/TM by F_lump for the literature ratio).
Gas comes from <session>/d/recon.mat (unchanged by the re-split, F70).

Session set = outputs/te90_2026-09-25/te90_table.csv rows with atlas_set == yes
(carrier ON_RBC and re-split not gated) -- Hooman's order 2026-09-25: RBC-carrier
sessions only, no mixing with the BETWEEN / ON_MEM / BELOW_MEM / ABOVE_RBC carriers,
sessions without inputs ignored. Every other Ext session is listed on the appendix page.

Layout per subject ID (dates sorted) = dis_atlas.py: panel pages gp / rbc / tm
(orientation-major rows, 10 slices at shared apex->base fractions, gas-picked bin,
gas-mask extent + fixed crop window), then per date full-slice pages cor/sag/ax x
gp/rbc/tm, and one mp4 per ID per kind cycling the 16 bins.

Outputs: workspace/outputs/aikill_atlas/rbc_tm_atlas.pdf, videos_rbctm/<ID>_{gp,rbc,tm}.mp4,
         rbc_tm_atlas_sessions.csv (what went in, angle, F_lump, bin shown)
Usage: rbc_tm_atlas.py [--only ID] [--no-video] [--no-full] [--root /Volumes/HoomHamExt/AIkill_Dynamic_b44 [--tag b44]]
  --root: production tree with the container split in recon.mat (C49/C50); outputs go to outputs/aikill_atlas_<tag>/
"""
from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io as sio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import imageio.v2 as imageio
from PIL import Image, ImageDraw

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dis_atlas import (ORIENTS, orient, montage, THRESH, MINVOX, PAD, NSEL, DPI,
                       FPS, REPEATS, UPSCALE)

SRC = Path('/Volumes/HoomHamExt/AIkill_Dynamic')
WS = Path(__file__).resolve().parents[2]
OUT = WS / 'outputs' / 'aikill_atlas'
VID = OUT / 'videos_rbctm'
TABLE = WS / 'outputs' / 'te90_2026-09-25' / 'te90_table.csv'
KINDS = ('gp', 'rbc', 'tm')
SRC_NOTE = ''
LABEL = {'gp': 'gas', 'rbc': 'RBC', 'tm': 'TM (membrane, lumped)'}


class Session:
    def __init__(self, key, vols, meta):
        self.key, self.date, self.sid = key, key[:10], key[11:]
        self.vols, self.meta = vols, meta
        gas = vols['gp']
        sig = gas.reshape(gas.shape[0], -1).sum(axis=1)
        self.best = int(np.argmax(sig))
        mean_vol = gas.mean(axis=0)
        self.mask = mean_vol > THRESH * (mean_vol.max() + 1e-9)
        self.vmaxes = {}
        for k, v in vols.items():
            sel = v[self.best][self.mask] if self.mask.any() else v[self.best]
            self.vmaxes[k] = float(np.percentile(sel, 99.5) + 1e-9)
        self.extent = {}
        for axis, _ in ORIENTS:
            msk = np.moveaxis(self.mask, axis - 1, 0)
            idx = np.where(msk.sum(axis=(1, 2)) >= MINVOX)[0]
            if len(idx) == 0:
                idx = np.where(msk.any(axis=(1, 2)))[0]
            if len(idx) == 0:
                idx = np.arange(msk.shape[0])
            plane = msk[idx].any(axis=0)
            rows = np.where(plane.any(axis=1))[0]
            cols = np.where(plane.any(axis=0))[0]
            r0, r1 = max(rows[0] - PAD, 0), min(rows[-1] + PAD + 1, msk.shape[1])
            c0, c1 = max(cols[0] - PAD, 0), min(cols[-1] + PAD + 1, msk.shape[2])
            self.extent[axis] = (int(idx[0]), int(idx[-1]), r0, r1, c0, c1)

    def vol(self, kind, b=None):
        return self.vols[kind][self.best if b is None else b]

    def vmax(self, kind):
        return self.vmaxes[kind]

    def tag(self):
        m = self.meta
        return (f'φk0={m["dphi"]:.0f}°  F_lump={m["F_lump"]:.2f}  '
                f'bin {self.best} (insp {m["insp_bin"]})')


class IDGroup:
    def __init__(self, sessions):
        self.sessions = sessions
        self.size, self.window, self.sel = {}, {}, {}
        for axis, _ in ORIENTS:
            ex = [s.extent[axis] for s in sessions]
            H = max(e[3] - e[2] for e in ex)
            W = max(e[5] - e[4] for e in ex)
            self.size[axis] = (H, W)
            nr, nc = np.moveaxis(sessions[0].mask, axis - 1, 0).shape[1:]
            fr = np.linspace(0, 1, NSEL)
            if axis == 2:
                fr = fr[::-1]
            for s in sessions:
                lo, hi, r0, r1, c0, c1 = s.extent[axis]
                wr = int(np.clip(round((r0 + r1) / 2 - H / 2), 0, nr - H))
                wc = int(np.clip(round((c0 + c1) / 2 - W / 2), 0, nc - W))
                self.window[(s, axis)] = (wr, wc)
                self.sel[(s, axis)] = np.round(lo + fr * (hi - lo)).astype(int)

    def tiles(self, sess, kind, axis, slices, b=None):
        vol = np.moveaxis(sess.vol(kind, b), axis - 1, 0)
        H, W = self.size[axis]
        r0, c0 = self.window[(sess, axis)]
        return [orient(vol[s, r0:r0 + H, c0:c0 + W], axis) for s in slices]


def panel_rows(idg, kind, b=None):
    rows = []
    for axis, name in ORIENTS:
        for sess in idg.sessions:
            sel = idg.sel[(sess, axis)]
            rows.append((sess, kind, name, montage(
                idg.tiles(sess, kind, axis, sel, b), len(sel)), sess.vmax(kind)))
    return rows


def render_panel_page(pdf, sid, idg, kind):
    rowdata = panel_rows(idg, kind)
    width_in = 13.0
    heights = [img.shape[0] / img.shape[1] for _, _, _, img, _ in rowdata]
    fig_h = sum(heights) * (width_in * 0.93) + 0.7
    fig = plt.figure(figsize=(width_in, min(fig_h, 200)), facecolor='black')
    gs = fig.add_gridspec(len(rowdata), 1, height_ratios=heights,
                          left=0.065, right=1.0, top=1 - 0.45 / fig_h,
                          bottom=0.0, hspace=0.0)
    for i, (sess, kind, name, img, vmax) in enumerate(rowdata):
        ax = fig.add_subplot(gs[i], facecolor='black')
        ax.imshow(img, cmap='gray', vmin=0, vmax=vmax, interpolation='nearest')
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.set_ylabel(f'{name}\n{sess.date}\n{kind}', fontsize=5, rotation=0,
                      ha='right', va='center', labelpad=16, color='white')
    sub = ' · '.join(f'{s.date}: {s.tag()}' for s in idg.sessions)
    fig.suptitle(f'{sid} — {LABEL[kind]}, best gas bin per session\n{sub}',
                 fontsize=8, color='white')
    pdf.savefig(fig, dpi=DPI, facecolor='black')
    plt.close(fig)


def render_full_page(pdf, sid, idg, sess, kind, axis, name):
    lo, hi = sess.extent[axis][0], sess.extent[axis][1]
    full = np.arange(lo, hi + 1)
    if axis == 2:
        full = full[::-1]
    img = montage(idg.tiles(sess, kind, axis, full), NSEL)
    h_over_w = img.shape[0] / img.shape[1]
    fig_h = max(13 * h_over_w + 0.5, 3)
    fig = plt.figure(figsize=(13, fig_h), facecolor='black')
    ax = fig.add_axes([0.0, 0.0, 1.0, 1 - 0.4 / fig_h], facecolor='black')
    ax.imshow(img, cmap='gray', vmin=0, vmax=sess.vmax(kind), interpolation='nearest')
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    fig.suptitle(f'{sid}  {sess.date} — {name} {LABEL[kind]}, slices {lo}–{hi} '
                 f'({sess.tag()})', fontsize=9, color='white')
    pdf.savefig(fig, dpi=DPI, facecolor='black')
    plt.close(fig)


def render_video(sid, idg, kind):
    nbins = idg.sessions[0].vols['gp'].shape[0]
    per_bin = [panel_rows(idg, kind, b) for b in range(nbins)]
    maxw = max(img.shape[1] for rows in per_bin for _, _, _, img, _ in rows)
    lab_w = 100
    frames = []
    for b, rows in enumerate(per_bin):
        parts = []
        for sess, kind, name, img, vmax in rows:
            u8 = (np.clip(img / vmax, 0, 1) * 255).astype(np.uint8)
            padded = np.zeros((u8.shape[0], maxw), dtype=np.uint8)
            padded[:, :u8.shape[1]] = u8
            parts.append(padded)
        canvas = np.vstack(parts)
        im = Image.fromarray(canvas).resize(
            (canvas.shape[1] * UPSCALE, canvas.shape[0] * UPSCALE), Image.NEAREST)
        full = Image.new('L', (im.width + lab_w, im.height + 24), 0)
        full.paste(im, (lab_w, 24))
        d = ImageDraw.Draw(full)
        d.text((4, 6), f'{sid} {kind}  bin {b:02d}/{nbins - 1}', fill=255)
        y = 24
        for sess, kind, name, img, _ in rows:
            d.text((4, y + img.shape[0] * UPSCALE // 2 - 5),
                   f'{name[:3]} {sess.date[2:]} {kind}', fill=200)
            y += img.shape[0] * UPSCALE
        arr = np.asarray(full)
        arr = arr[:arr.shape[0] - arr.shape[0] % 2, :arr.shape[1] - arr.shape[1] % 2]
        frames.append(arr)
    w = imageio.get_writer(VID / f'{sid}_{kind}.mp4', fps=FPS, codec='libx264',
                           quality=8, macro_block_size=None)
    for _ in range(REPEATS):
        for f in frames:
            w.append_data(f)
    w.close()


def load_session_container(key, root):
    """Production tree (C49/C50): recon.mat holds the container split (dissolved_phase_real/imag = aRBC/aTP
    when rbc_tp_separated == '1') plus rbc_tp_* meta. Sessions the container gate left unsplit raise."""
    m = sio.loadmat(root / key / 'd' / 'recon.mat',
                    variable_names=['gas_phase', 'dissolved_phase_real', 'dissolved_phase_imag', 'rbc_tp_separated',
                                    'rbc_tp_dphi_deg', 'rbc_tp_F_lump', 'rbc_tp_ratio_scalar', 'rbc_tp_target'])
    if str(m.get('rbc_tp_separated', [''])[0]).strip() != '1':
        raise RuntimeError('container left this session unsplit (gate)')
    g = m['gas_phase']
    vols = {'gp': g.astype(np.float32),
            'rbc': np.maximum(m['dissolved_phase_real'], 0).astype(np.float32),
            'tm': np.maximum(m['dissolved_phase_imag'], 0).astype(np.float32)}
    f = lambda k, d=float('nan'): float(str(m[k][0]).strip()) if k in m else d
    gas_sig = g.reshape(g.shape[0], -1).sum(axis=1)
    meta = {'dphi': f('rbc_tp_dphi_deg'), 'F_lump': f('rbc_tp_F_lump'), 'insp_bin': int(np.argmax(gas_sig)),
            'ratio_scalar': f('rbc_tp_ratio_scalar'), 'carrier': 'ON_RBC', 'gate': 'container'}
    return Session(key, vols, meta)


def load_session(key):
    p = SRC / key
    g = sio.loadmat(p / 'd' / 'recon.mat', variable_names=['gas_phase'])['gas_phase']
    r = sio.loadmat(p / 'd' / 'recon_resplit_xecs.mat',
                    variable_names=['aRBC', 'aTP', 'dphi_new_deg', 'F_lump', 'insp_bin',
                                    'ratio_scalar', 'carrier', 'gate_strict'])
    vols = {'gp': g.astype(np.float32),
            'rbc': np.maximum(r['aRBC'], 0).astype(np.float32),
            'tm': np.maximum(r['aTP'], 0).astype(np.float32)}
    sc = lambda k: r[k].ravel()[0]
    meta = {'dphi': float(sc('dphi_new_deg')), 'F_lump': float(sc('F_lump')),
            'insp_bin': int(sc('insp_bin')), 'ratio_scalar': float(sc('ratio_scalar')),
            'carrier': str(r['carrier'][0]).strip(), 'gate': str(r['gate_strict'][0]).strip()}
    return Session(key, vols, meta)


def appendix(pdf, table, included, skipped):
    """Text page: which sessions are in, which are out and why."""
    t = table.set_index('key')
    ext = sorted(p.name for p in SRC.iterdir()
                 if (p / 'd' / 'recon.mat').exists() and '_merged' not in p.name
                 and '_REJECTED' not in p.name)
    lines = [f'RBC/TM atlas — {len(included)} sessions in (carrier ON_RBC, XeCS re-split, not gated). '
             f'Source maps: <session>/d/recon_resplit_xecs.mat (2026-09-24). LUMPED scale.', '']
    if SRC_NOTE:
        lines[0] = f'RBC/TM atlas — {len(included)} sessions. {SRC_NOTE} LUMPED scale.'
    lines.append('OUT (Ext sessions not in this atlas):')
    for k in ext:
        if k in included:
            continue
        if k in t.index:
            row = t.loc[k]
            why = f'carrier {row.carrier} ({row.carrier_ppm} ppm)' if row.carrier != 'ON_RBC' \
                else f'gated ({row.status}, stable={row.xecs_stable})'
        else:
            why = 'no cal-block fit / no dissolved recon'
        lines.append(f'  {k:22s} {why}')
    lines += ['', 'Load errors:'] + [f'  {k}: {e}' for k, e in skipped]
    fig = plt.figure(figsize=(13, 0.16 * len(lines) + 1), facecolor='white')
    fig.text(0.02, 0.98, '\n'.join(lines), va='top', family='monospace', fontsize=7)
    pdf.savefig(fig, dpi=DPI)
    plt.close(fig)


def main():
    argv = sys.argv[1:]
    only = argv[argv.index('--only') + 1] if '--only' in argv else None
    do_video, do_full = '--no-video' not in argv, '--no-full' not in argv
    root = Path(argv[argv.index('--root') + 1]) if '--root' in argv else None      # production tree (container split)
    tag = argv[argv.index('--tag') + 1] if '--tag' in argv else ('b44' if root else '')
    global OUT, VID, SRC_NOTE
    if root:
        SRC_NOTE = f'Source: {root}/<key>/d/recon.mat — container split, calcb phase low-pass 4.4 vox (fork 04f445b, 2026-09-25).'
    if tag:
        OUT = OUT.parent / f'aikill_atlas_{tag}'; VID = OUT / 'videos_rbctm'
    OUT.mkdir(parents=True, exist_ok=True)
    VID.mkdir(exist_ok=True)
    table = pd.read_csv(TABLE)
    keys = sorted(table[table.atlas_set.eq('yes')].key)
    groups = defaultdict(list)
    for k in keys:
        groups[k[11:]].append(k)
    if only:
        groups = {only: groups[only]}
    skipped, included, rows = [], [], []
    pdf_path = OUT / (f'rbc_tm_sample_{only}.pdf' if only else 'rbc_tm_atlas.pdf')
    with PdfPages(pdf_path) as pdf:
        for sid in sorted(groups):
            sessions = []
            for key in sorted(groups[sid]):
                try:
                    sessions.append(load_session_container(key, root) if root else load_session(key))
                    included.append(key)
                except Exception as e:
                    skipped.append((key, repr(e)))
            if not sessions:
                continue
            idg = IDGroup(sessions)
            for kind in KINDS:
                render_panel_page(pdf, sid, idg, kind)
            if do_full:
                for sess in sessions:
                    for axis, name in ORIENTS:
                        for kind in KINDS:
                            render_full_page(pdf, sid, idg, sess, kind, axis, name)
            if do_video:
                for kind in KINDS:
                    render_video(sid, idg, kind)
            for s in sessions:
                rows.append({'key': s.key, 'sid': sid, 'bin_shown': s.best, **s.meta,
                             'vmax_gp': s.vmaxes['gp'], 'vmax_rbc': s.vmaxes['rbc'],
                             'vmax_tm': s.vmaxes['tm']})
            print(f'[atlas] {sid}: {len(sessions)} session(s)', flush=True)
        if not only:
            appendix(pdf, table, set(included), skipped)
    pd.DataFrame(rows).to_csv(OUT / ('rbc_tm_sample_sessions.csv' if only
                                     else 'rbc_tm_atlas_sessions.csv'), index=False)
    print(f'[atlas] PDF -> {pdf_path}  ({len(included)} sessions)')
    for s in skipped:
        print('[atlas] skipped', *s)


if __name__ == '__main__':
    sys.exit(main())
