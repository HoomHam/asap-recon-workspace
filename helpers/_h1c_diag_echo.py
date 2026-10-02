"""Scratch: is the fixed misfit a shifted copy of k-space (unspoiled magnetisation from earlier TRs)?
C(q) = sum_j r_j conj(X(k_j + o_j + q)) for all q at once = IFFT[conj(x) * A'^H r], A' on the points k_j + o_j.
Hypotheses for the offset o_j: none; end point of the previous interleave; of the previous two."""
import os, sys, numpy as np, finufft
from scipy import ndimage as ndi
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _h1c as H, _h1c_dyn as DY
G = DY.G; G2 = tuple(2 * g for g in G); dx2 = H.DX / 2
ybar = np.load(os.path.join(H.EXT, 'ybar_train.npy')).astype(np.complex128)[0]
k = H.load()['k']; op = H.Nufft(k, np.arange(H.NILV), n=G)
x = np.load(os.path.join(H.EXT, 'static_g_ls_ch0.npy')).astype(np.complex128)
r = ybar - op.fwd(x)
r[:, :48] = 0                                                   # mid / high band misfit only
X = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(x)))
Xp = np.zeros(G2, complex); c = [g // 2 for g in G]
Xp[c[0]:c[0] + G[0], c[1]:c[1] + G[1], c[2]:c[2] + G[2]] = X
xu = np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(Xp))) * 8    # x on the 1.75 mm grid
kend = k[:, -1]                                                  # (832, 3)
hyp = {'no offset': np.zeros((H.NILV, 3)), 'previous interleave end': np.roll(kend, 1, axis=0), 'minus previous end': -np.roll(kend, 1, axis=0),
       'previous two ends': np.roll(kend, 1, axis=0) + np.roll(kend, 2, axis=0), 'own end (same interleave, earlier pass)': kend}
qax = [np.fft.fftshift(np.fft.fftfreq(g, dx2)) for g in G2]
for name, o in hyp.items():
    kk = k[:, H.KILL:] + o[:, None, :]
    kap = (2 * np.pi * dx2 * kk.reshape(-1, 3))
    kap = (kap + np.pi) % (2 * np.pi) - np.pi
    g = finufft.nufft3d1(kap[:, 0].copy(), kap[:, 1].copy(), kap[:, 2].copy(), r.ravel(), G2, isign=+1, eps=1e-4)
    C = np.abs(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(np.conj(xu) * g))))
    med = np.median(C); i = np.unravel_index(np.argmax(C), C.shape)
    q = [qax[a][i[a]] for a in range(3)]
    # second peak away from the first
    C2 = C.copy(); C2[max(i[0] - 6, 0):i[0] + 7, max(i[1] - 6, 0):i[1] + 7, max(i[2] - 6, 0):i[2] + 7] = 0; j = np.unravel_index(np.argmax(C2), C.shape)
    print(f'{name:42s}: peak/median {C.max()/med:7.1f} at q = ({q[0]:+.4f}, {q[1]:+.4f}, {q[2]:+.4f}) cyc/mm |q| {np.linalg.norm(q):.4f} ; 2nd peak/median {C2.max()/med:6.1f} at |q| {np.linalg.norm([qax[a][j[a]] for a in range(3)]):.4f} ; 99.9 pct/median {np.percentile(C, 99.9)/med:.1f}', flush=True)
