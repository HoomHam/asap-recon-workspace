"""Scratch: XeCS's off-resonance candidate for the fixed misfit. Demodulate the residual of the support-limited static
fit (channel 0) by exp(-i 2 pi df t) and grid it: at which df does its energy focus, and where in the body?"""
import os, sys, json, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _h1c as H, _h1c_dyn as DY
G = DY.G
ybar = np.load(os.path.join(H.EXT, 'ybar_train.npy')).astype(np.complex128)[0]
k = H.load()['k']; op = H.Nufft(k, np.arange(H.NILV), n=G)
dcf = np.load(os.path.join(H.EXT, 'dcf_g.npy')); Wt = dcf / dcf.sum()
x = np.load(os.path.join(H.EXT, 'static_g_ls_ch0.npy')).astype(np.complex128); M = np.load(os.path.join(H.EXT, 'support_g.npy'))
r = ybar - op.fwd(x); r[:, :48] = 0
t = (np.arange(H.KILL, H.NPTS) * H.DWELL)[None, :]
X, Y, Z = np.meshgrid(*[(np.arange(n) - n // 2) * H.DX for n in G], indexing='ij')
reg = {'arms (|LR| > 210 mm)': M & (np.abs(X) > 210), 'SI ends (|SI| > 180 mm)': M & (np.abs(Z) > 180) & (np.abs(X) <= 210),
       'central torso': M & (np.abs(Z) <= 180) & (np.abs(X) <= 210), 'outside the body': ~M}
dfs = np.arange(-2000, 2001, 100)
E = {k_: [] for k_ in reg}; tot = []
for df in dfs:
    g = op.adj(Wt * r * np.exp(-2j * np.pi * df * t)); p = np.abs(g) ** 2
    tot.append(p.sum())
    for k_, m in reg.items(): E[k_].append(p[m].mean())
    if df in (-1000, -500, -200, 0, 200, 500, 1000): np.save(os.path.join(H.EXT, f'resid_demod_{df:+d}.npy'), np.abs(g).astype(np.float32))
bg = np.array(E['outside the body'])
print('df (Hz):            ', ' '.join(f'{d:6d}' for d in dfs[::2]))
for k_ in list(reg)[:3]:
    print(f'{k_:24s}', ' '.join(f'{v:6.2f}' for v in (np.array(E[k_]) / bg)[::2]), ' <- mean residual-image power / power outside the body')
# same sweep on the DATA (not the residual): where does each region focus? sharpness = gradient energy / energy in region
print('data, image gradient energy by region (normalised to df = 0):')
S = {k_: [] for k_ in list(reg)[:3]}
for df in dfs:
    g = np.abs(op.adj(Wt * ybar * np.exp(-2j * np.pi * df * t)))
    gr = sum(np.gradient(g, axis=a) ** 2 for a in range(3))
    for k_ in S: S[k_].append(gr[reg[k_]].sum() / (g[reg[k_]] ** 2).sum())
i0 = list(dfs).index(0)
for k_ in S:
    s = np.array(S[k_]) / S[k_][i0]
    print(f'{k_:24s}', ' '.join(f'{v:6.3f}' for v in s[::2]), f' best df {dfs[np.argmax(s)]:+d} Hz')
json.dump(dict(df=dfs.tolist(), resid_power={k_: list(map(float, v)) for k_, v in E.items()}, sharp={k_: list(map(float, v)) for k_, v in S.items()}), open(os.path.join(H.OUT, 'diag_offres.json'), 'w'), indent=1)
