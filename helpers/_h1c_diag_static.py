"""Scratch: which constraint leaves the fixed-pattern misfit? Static fits of the pass-averaged data on grid G:
one channel or two (sensitivity ratio), body support on or off. Residual x13 = noise units of the 13-pass mean."""
import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _h1c as H, _h1c_dyn as DY
G = DY.G
ybar = np.load(os.path.join(H.EXT, 'ybar_train.npy')).astype(np.complex128)
k = H.load()['k']; op = H.Nufft(k, np.arange(H.NILV), n=G)
dcf = np.load(os.path.join(H.EXT, 'dcf_g.npy')); Wt = dcf / dcf.sum()
S = np.load(os.path.join(H.EXT, 'sens_g.npy')).astype(np.complex128); M = np.load(os.path.join(H.EXT, 'support_g.npy'))
bands = ((0, 8), (8, 48), (48, 198), (198, 510))
sig = [[(np.abs(ybar[c][:, a:b]) ** 2).mean() * 13 for a, b in bands] for c in range(2)]
print('signal x13: ch0', np.round(sig[0], 1), 'ch1', np.round(sig[1], 1))
def run(chs, mask, tag, its=(15, 40)):
    Sx = [S[c] if len(chs) == 2 else np.ones(G) for c in chs]
    m = M if mask else 1.0
    rhs = m * sum(np.conj(s_) * op.adj(Wt * ybar[c]) for s_, c in zip(Sx, chs))
    nrm0 = lambda x: m * sum(np.conj(s_) * op.adj(Wt * op.fwd(s_ * x)) for s_ in Sx)
    sc = np.vdot(rhs, nrm0(rhs)).real / np.vdot(rhs, rhs).real
    def cb(i, x, r):
        if (i + 1) in its:
            e = [[(np.abs((op.fwd(s_ * x) - ybar[c])[:, a:b]) ** 2).mean() * 13 for a, b in bands] for s_, c in zip(Sx, chs)]
            print(f'{tag:34s} it {i+1:3d}: ' + ' | '.join(f'ch{c} ' + ' '.join(f'{v:9.2f}' for v in e_) for c, e_ in zip(chs, e)), flush=True)
    return H.cg(lambda x: nrm0(x) + 1e-4 * sc * x, rhs, it=max(its), cb=cb)
run([0], False, 'ch0 alone, no support')
x0 = run([0], True, 'ch0 alone, support')
run([1], True, 'ch1 alone, support')
run([0, 1], False, 'two channels (R), no support')
x = run([0, 1], True, 'two channels (R), support', its=(15, 40, 100))
np.save(os.path.join(H.EXT, 'static_g_ls.npy'), x.astype(np.complex64)); np.save(os.path.join(H.EXT, 'static_g_ls_ch0.npy'), x0.astype(np.complex64))
