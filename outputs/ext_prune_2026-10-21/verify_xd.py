"""Verify <x>/input.mrd == d/input.mrd except the binning string (x = s: SIGNAL, p: PNEUMOTACH; d: DIAPHRAGM).
The first differing byte must be the binning string's 1-byte varint length; everything after the two strings must match.
Line format in: '<size_x> <size_d> <delta> <relpath x/input.mrd>'. Prints OK/FAIL per pair, full streaming compare."""
import sys

ROOT = '/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/'
CH = 64 << 20
WANT = {'s': b'SIGNAL', 'p': b'PNEUMOTACH'}

for line in open(sys.argv[1]):
    rel = line.split()[3]
    x = ROOT + rel
    letter = rel.split('/')[-2]
    d = x[:-len(f'{letter}/input.mrd')] + 'd/input.mrd'
    why, n = '', 0
    with open(x, 'rb') as fx, open(d, 'rb') as fd:
        hx, hd = fx.read(1 << 20), fd.read(1 << 20)
        i = next((k for k in range(min(len(hx), len(hd))) if hx[k] != hd[k]), None)
        if i is None:
            why = 'no diff in head'
        else:
            sx, sd = WANT[letter], b'DIAPHRAGM'
            ok = hx[i:i + 1 + len(sx)] == bytes([len(sx)]) + sx and hd[i:i + 1 + len(sd)] == bytes([len(sd)]) + sd
            if not ok:
                why = f'first diff at {i} is not the binning string'
            else:
                fx.seek(i + 1 + len(sx))
                fd.seek(i + 1 + len(sd))
                while True:
                    a, b = fx.read(CH), fd.read(CH)
                    if a != b:
                        why = f'payload differs after byte ~{n}'
                        break
                    if not a:
                        break
                    n += len(a)
    print(('FAIL' if why else 'OK  '), i, n, rel, why, flush=True)
