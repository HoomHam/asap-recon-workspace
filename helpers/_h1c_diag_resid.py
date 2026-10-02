"""Scratch: what is the misfit of the 4D least-squares model made of? Residual per line (training lines) split into
the part that repeats for every pass of an interleave (fixed pattern: calibration / static model) and the part that
changes from pass to pass (dynamics not in the model + noise). Per readout band, both channels, noise units."""
import os, sys, time, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _h1c as H, _h1c_dyn as DY
d = DY.Dyn()
train = d.kit['scored_train_lines']
t0 = time.time(); x = d.recon(train, lam_s=1e-3, lam_t=0.02, it=12); print('recon %.0f s' % (time.time() - t0))
np.save(os.path.join(H.EXT, 'x_diag_full.npy'), x.astype(np.complex64))
U = d.fwd_all(x)
r = d.predict(x, train, U) - d.v[:1, train, H.KILL:]              # (1, lines, 510)
y = d.v[:1, train, H.KILL:]
il = train % H.NILV; cnt = np.bincount(il, minlength=H.NILV)
bands = ((0, 8), (8, 48), (48, 198), (198, 510))
out = {}
for c in range(1):
    rb = np.zeros((H.NILV, 510), np.complex128); np.add.at(rb, il, r[c]); rb[cnt > 0] /= cnt[cnt > 0][:, None]; rb = rb[cnt > 0]
    for a, b in bands:
        tot = (np.abs(r[c][:, a:b]) ** 2).mean()
        fixed = (np.abs(rb[:, a:b]) ** 2).mean() - 1.0 / cnt[cnt > 0].mean()            # mean-of-13 noise removed
        sig = (np.abs(y[c][:, a:b]) ** 2).mean()
        print(f'ch{c} samples {a+2:3d}-{b+2:3d}: signal {sig:10.1f}  residual {tot:9.2f}  fixed-pattern part {fixed:9.2f}  varying part {tot - fixed:9.2f}  (noise = 1)')
        out[f'ch{c}_{a+2}_{b+2}'] = dict(signal=float(sig), resid=float(tot), fixed=float(fixed))
    np.save(os.path.join(H.EXT, f'resid_fixed_ch{c}.npy'), rb.astype(np.complex64))
json.dump(out, open(os.path.join(H.OUT, 'diag_resid.json'), 'w'), indent=1)
# is the varying part breathing the model missed? correlation of per-line low-band residual power with the surrogate speed
e = (np.abs(r[0][:, 8:48]) ** 2).mean(1)
vol = H.load()['vol']; sp = np.abs(np.gradient(vol))[train]
print('corr(low-band residual power, |d vol/dt|) = %.3f ; corr with vol = %.3f' % (np.corrcoef(e, sp)[0, 1], np.corrcoef(e, vol[train])[0, 1]))
t = train * 0.0172
f = np.fft.rfftfreq(4096, 0.0172)
# spectrum of the per-line residual energy (uniform in time apart from held-out gaps): cardiac line near 1-1.3 Hz?
ee = np.zeros(d.L); ee[train] = e - e.mean()
P = np.abs(np.fft.rfft(ee[H.NILV:H.NILV + 8192] * np.hanning(8192))) ** 2; ff = np.fft.rfftfreq(8192, 0.0172)
top = np.argsort(P[5:])[::-1][:8] + 5
print('strongest lines in residual-energy spectrum (Hz):', np.round(np.sort(ff[top]), 3))
