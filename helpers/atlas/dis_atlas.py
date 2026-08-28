#!/usr/bin/env python
"""Gas + dissolved atlas over /Volumes/HoomHamExt/AIkill_Dynamic.

Groups sessions by subject ID (dates sorted). Per ID:
  - one panel page, rows grouped ORIENTATION-major:
      coronal block: for each date -> gp row, dp row  (10 slices each)
      sagittal block: same
      axial block:   same
    Slice picking, masks, bboxes all come from GAS (gp) — dp is lower SNR.
    Slice indices + in-plane bbox are shared across the ID's sessions
    (union extent), so columns align between repeats.
  - then per date, full-slice pages: cor-gp, cor-dp, sag-gp, sag-dp, ax-gp, ax-dp
  - one mp4 per ID: panel cycling the 16 bins, 5 repeats

Source: <session>/d/recon.mat -> gas_phase + dissolved_phase_magnitude (16,Z,Y,X).
If d has no dissolved (RBC/TP params were absent), gas comes from s/recon.mat
and the dp rows/pages are skipped.
Displayed bin (PDF) = highest gas total-signal bin.
Orientation follows workspace/pipeline/post_process.py:
  axial=Z flipud, coronal=Y reversed order, sagittal=X rot90 ccw.
Black background, no gaps — tiles touch.

Outputs: workspace/outputs/aikill_atlas/dissolved_atlas.pdf + videos/<ID>.mp4
Usage: dis_atlas.py [--only ID]
"""
from __future__ import annotations

import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import scipy.io as sio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import imageio.v2 as imageio
from PIL import Image, ImageDraw

SRC = Path('/Volumes/HoomHamExt/AIkill_Dynamic')
OUT = Path(__file__).resolve().parents[2] / 'outputs' / 'aikill_atlas'
VID = OUT / 'videos'
THRESH = 0.15          # gas data mask: > 15% of max
MINVOX = 30            # slice is data-bearing if >= this many mask voxels
PAD = 2                # in-plane bbox pad
NSEL = 10              # slices per row
DPI = 150
FPS = 5
REPEATS = 5
UPSCALE = 2

ORIENTS = ((2, 'coronal'), (3, 'sagittal'), (1, 'axial'))


def orient(tile, axis):
    if axis == 1:
        return np.flipud(tile)
    if axis == 3:
        return np.rot90(tile, 1)
    return tile


class Session:
    def __init__(self, date, gas, dp, src):
        self.date, self.gas, self.dp, self.src = date, gas, dp, src
        sig = gas.reshape(gas.shape[0], -1).sum(axis=1)
        self.best = int(np.argmax(sig))
        mean_vol = gas.mean(axis=0)
        self.mask = mean_vol > THRESH * (mean_vol.max() + 1e-9)
        self.gvmax = float(np.percentile(gas[self.best][self.mask], 99.5)) \
            if self.mask.any() else float(gas[self.best].max() + 1e-9)
        self.dvmax = None
        if dp is not None:
            sel = dp[self.best][self.mask] if self.mask.any() else dp[self.best]
            self.dvmax = float(np.percentile(sel, 99.5) + 1e-9)
        self.extent = {}                    # axis -> (lo, hi, r0, r1, c0, c1)
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
        v = self.gas if kind == 'gp' else self.dp
        return v[self.best if b is None else b]

    def vmax(self, kind):
        return self.gvmax if kind == 'gp' else self.dvmax


class IDGroup:
    """Anatomy-aligned columns across repeats: per-session slice picking at
    shared apex->base fractions, plus a fixed-size in-plane crop window
    centered on each session's own lung bbox (shift+scale alignment,
    no registration)."""

    def __init__(self, sessions):
        self.sessions = sessions
        self.size = {}                      # axis -> (H, W) shared tile size
        self.window = {}                    # (sess, axis) -> (r0, c0)
        self.sel = {}                       # (sess, axis) -> slice indices
        for axis, _ in ORIENTS:
            ex = [s.extent[axis] for s in sessions]
            H = max(e[3] - e[2] for e in ex)
            W = max(e[5] - e[4] for e in ex)
            self.size[axis] = (H, W)
            nr, nc = np.moveaxis(sessions[0].mask, axis - 1, 0).shape[1:]
            fr = np.linspace(0, 1, NSEL)
            if axis == 2:                   # coronal: reverse order
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


def montage(tiles, ncols):
    h, w = tiles[0].shape
    nrows = int(np.ceil(len(tiles) / ncols))
    out = np.zeros((nrows * h, ncols * w), dtype=np.float32)
    for i, t in enumerate(tiles):
        r, c = divmod(i, ncols)
        out[r * h:(r + 1) * h, c * w:(c + 1) * w] = t
    return out


def panel_rows(idg, kind, b=None):
    """[(sess, kind, name, montage, vmax)] orientation-major, date-minor."""
    rows = []
    for axis, name in ORIENTS:
        for sess in idg.sessions:
            if kind == 'dp' and sess.dp is None:
                continue
            sel = idg.sel[(sess, axis)]
            rows.append((sess, kind, name, montage(
                idg.tiles(sess, kind, axis, sel, b), len(sel)),
                sess.vmax(kind)))
    return rows


def render_panel_page(pdf, sid, idg, kind):
    rowdata = panel_rows(idg, kind)
    if not rowdata:
        return
    width_in = 13.0
    heights = [img.shape[0] / img.shape[1] for _, _, _, img, _ in rowdata]
    fig_h = sum(heights) * (width_in * 0.93) + 0.6
    fig = plt.figure(figsize=(width_in, min(fig_h, 200)), facecolor='black')
    gs = fig.add_gridspec(len(rowdata), 1, height_ratios=heights,
                          left=0.065, right=1.0, top=1 - 0.35 / fig_h,
                          bottom=0.0, hspace=0.0)
    for i, (sess, kind, name, img, vmax) in enumerate(rowdata):
        ax = fig.add_subplot(gs[i], facecolor='black')
        ax.imshow(img, cmap='gray', vmin=0, vmax=vmax, interpolation='nearest')
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.set_ylabel(f'{name}\n{sess.date}\n{kind}', fontsize=5, rotation=0,
                      ha='right', va='center', labelpad=16, color='white')
    label = 'gas (gp)' if kind == 'gp' else 'dissolved (dp)'
    fig.suptitle(f'{sid} — {label}, best gas bin per session',
                 fontsize=10, color='white')
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
    ax.imshow(img, cmap='gray', vmin=0, vmax=sess.vmax(kind),
              interpolation='nearest')
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    fig.suptitle(f'{sid}  {sess.date} — {name} {kind}, slices {lo}–{hi} '
                 f'(bin {sess.best}, gas-picked)', fontsize=10, color='white')
    pdf.savefig(fig, dpi=DPI, facecolor='black')
    plt.close(fig)


def render_video(sid, idg, kind):
    nbins = idg.sessions[0].gas.shape[0]
    per_bin = [panel_rows(idg, kind, b) for b in range(nbins)]
    if not per_bin[0]:
        return
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


def load_session(date, p):
    """d/recon.mat primary; if dissolved absent there, gas from s/recon.mat."""
    m = sio.loadmat(p / 'd' / 'recon.mat',
                    variable_names=['gas_phase', 'dissolved_phase_magnitude'])
    gas, dp, src = m.get('gas_phase'), m.get('dissolved_phase_magnitude'), 'd'
    if dp is None:
        src = 's (no dissolved in d)'
        dp = None
        sm = sio.loadmat(p / 's' / 'recon.mat', variable_names=['gas_phase'])
        gas = sm['gas_phase']
    return Session(date, gas.astype(np.float32),
                   None if dp is None else dp.astype(np.float32), src)


def main():
    only = sys.argv[sys.argv.index('--only') + 1] if '--only' in sys.argv else None
    OUT.mkdir(parents=True, exist_ok=True)
    VID.mkdir(exist_ok=True)
    pat = re.compile(r'^(\d{4}-\d{2}-\d{2})_(\w+)$')
    groups = defaultdict(list)
    for p in sorted(SRC.iterdir()):
        m = pat.match(p.name)
        if m and (p / 'd' / 'recon.mat').exists():
            groups[m.group(2)].append((m.group(1), p))
    if only:
        groups = {only: groups[only]}
    skipped = []
    pdf_path = OUT / (f'sample_{only}.pdf' if only else 'dissolved_atlas.pdf')
    with PdfPages(pdf_path) as pdf:
        for sid in sorted(groups):
            sessions = []
            for date, p in sorted(groups[sid]):
                try:
                    sessions.append(load_session(date, p))
                except Exception as e:
                    skipped.append((sid, date, repr(e)))
            if not sessions:
                skipped.append((sid, '-', 'no usable sessions'))
                continue
            idg = IDGroup(sessions)
            render_panel_page(pdf, sid, idg, 'gp')
            render_panel_page(pdf, sid, idg, 'dp')
            for sess in sessions:
                for axis, name in ORIENTS:
                    render_full_page(pdf, sid, idg, sess, 'gp', axis, name)
                    if sess.dp is not None:
                        render_full_page(pdf, sid, idg, sess, 'dp', axis, name)
            render_video(sid, idg, 'gp')
            render_video(sid, idg, 'dp')
            note = '; '.join(s.src for s in sessions if s.src != 'd')
            print(f'[atlas] {sid}: {len(sessions)} session(s)'
                  + (f'  [{note}]' if note else ''), flush=True)
    print(f'[atlas] PDF -> {pdf_path}')
    if skipped:
        print('[atlas] skipped:')
        for s in skipped:
            print('   ', *s)


if __name__ == '__main__':
    sys.exit(main())
