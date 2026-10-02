"""Scratch (2026-10-01): copy an input.mrd, multiplying the dynamic (and reference) acquisitions by a
readout window w(n) = exp(-n*dwell/T) along the sample axis. Everything else (header, trajectories,
pneumotach) is passed through unchanged.

usage: _mrd_apodize.py <in.mrd> <out.mrd> <T_ms> [dwell_us=5]
"""
import sys
import numpy as np
import mrd

src, dst, T = sys.argv[1], sys.argv[2], float(sys.argv[3]) * 1e-3
dw = (float(sys.argv[4]) if len(sys.argv) > 4 else 5.0) * 1e-6

with mrd.BinaryMrdReader(src) as r:
    header = r.read_header()
    items = list(r.read_data())

n_mod = 0
out = []
for it in items:
    if isinstance(it, mrd.StreamItem.NdArrayComplexFloat) and (
            it.value.meta.get('dynamic_acquisition') or it.value.meta.get('reference_acquisition')):
        d = np.asarray(it.value.data)                         # (nch, npts, nilv) channels-first
        w = np.exp(-np.arange(d.shape[1]) * dw / T).astype(np.float32)
        d = np.ascontiguousarray(d * w[None, :, None]).astype(np.complex64)
        it = mrd.StreamItem.NdArrayComplexFloat(mrd.NdArray(data=d, meta=it.value.meta))
        n_mod += 1
    out.append(it)

with mrd.BinaryMrdWriter(dst) as wtr:
    wtr.write_header(header)
    wtr.write_data(out)
print(f'apodized {n_mod} acquisition array(s) with exp(-t/{T*1e3:g} ms) -> {dst}')
