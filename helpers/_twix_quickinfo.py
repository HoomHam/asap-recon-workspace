"""Scratch: cheap twix inventory — ASCCONV header fields (regex, no mapVBVD) + MDH walk (scan count, sample
lengths, first/last timestamp -> real duration). Usage: _twix_quickinfo.py <dat> [...] -> TSV on stdout."""
import re, struct, sys, collections, pathlib
KEYS = {'prot': r'tProtocolName\s*=\s*"([^"]*)"', 'nuc': r'sTXSPEC\.asNucleusInfo\[0\]\.tNucleus\s*=\s*"([^"]*)"',
        'TR_us': r'alTR\[0\]\s*=\s*(\d+)', 'nrep': r'lRepetitions\s*=\s*(\d+)', 'acqOrder': r'sWipMemBlock\.alFree\[9\]\s*=\s*(\d+)',
        'numRes-1': r'sWipMemBlock\.alFree\[8\]\s*=\s*(\d+)', 'FOVro': r'sSliceArray\.asSlice\[0\]\.dReadoutFOV\s*=\s*([\d.]+)',
        'nslc': r'sSliceArray\.lSize\s*=\s*(\d+)', 'thk': r'sSliceArray\.asSlice\[0\]\.dThickness\s*=\s*([\d.]+)'}
print('file\tMB\tprot\tnuc\tTR_ms\tacqOrder\tnscan\tnsamp\tnch\tdur_s\tnslc\tFOV')
for f in sys.argv[1:]:
    p = pathlib.Path(f)
    with open(f, 'rb') as fh:
        nmeas = struct.unpack('<II', fh.read(8))[1]; fh.seek(8 + 152 * (nmeas - 1) + 8)
        off = struct.unpack('<Q', fh.read(8))[0]; fh.seek(off); hl = struct.unpack('<I', fh.read(4))[0]
        hdr = fh.read(hl).decode('latin-1')
        m = hdr.rfind('### ASCCONV BEGIN'); asc = hdr[m:hdr.find('### ASCCONV END', m)] if m >= 0 else hdr
        v = {k: (re.search(r, asc).group(1) if re.search(r, asc) else '') for k, r in KEYS.items()}
        pos = off + hl; ns = collections.Counter(); t0 = t1 = None; nch = 0
        while True:
            fh.seek(pos); h = fh.read(192)
            if len(h) < 192: break
            dma = struct.unpack('<I', h[0:4])[0] & 0x1FFFFFF; ev = struct.unpack('<Q', h[40:48])[0]
            s, c = struct.unpack('<HH', h[48:52])
            if ev & 1 or dma == 0: break
            if s > 0:
                ns[s] += 1; nch = c; ts = struct.unpack('<I', h[12:16])[0] * 2.5e-3
                t0 = ts if t0 is None else t0; t1 = ts
            pos += dma
    dur = round(t1 - t0, 1) if t0 is not None else ''
    print(f"{p.parent.name}/{p.name}\t{p.stat().st_size/1e6:.0f}\t{v['prot']}\t{v['nuc']}\t{int(v['TR_us'])/1000 if v['TR_us'] else ''}\t{v['acqOrder']}\t{sum(ns.values())}\t{dict(ns)}\t{nch}\t{dur}\t{v['nslc']}\t{v['FOVro']}", flush=True)
