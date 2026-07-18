"""
What is the REAL averaged time of each DIAPHRAGM bin in the 4D image?

The 4D stack axis is amplitude-RANK (phase-ordered), not time. To place each
bin on a physical time axis you must, per interleave, compute its within-breath
elapsed time (seconds since its own end-expiration, EE), then take the SAME
soft-weighted average the gridder used to build that bin's image.

    t_bin[b] = sum_i w(b,i) * tau_i / sum_i w(b,i)
    w(b,i)   = exp(-bindist(b, binctr_i)^2 / bindist0sq)   (cyclic, matches recon.py)

This reproduces bin() (raw.py:153) + the gridder bin weight (recon.py:104-111)
and shows the bin->time lookup is NONLINEAR (equal-count clusters bins in time
near the EE/EI turning points).
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

def asap_bin(val):
    bins = val.copy()
    v = np.array([x for x in val if not np.isnan(x)])
    breathdir = np.zeros(len(v))
    if v[1] < v[0]:
        breathdir[0] = 1
    for i in range(1, len(v)):
        if v[i] < v[i - 1]:
            breathdir[i] = 1
    iv = np.sort(v[breathdir == 1]); ev = np.sort(v[breathdir == 0])
    maxiv = max(iv) * 1.00001; maxev = max(ev)
    vf = v.copy(); vf[breathdir == 1] = maxiv + maxev - v[breathdir == 1]
    vsort = np.sort(vf)
    out = np.array([np.searchsorted(vsort, x) for x in vf], dtype=float) / len(v)
    return out, breathdir

def find_EE(vol, N=30):
    """Two-pass local-minimum finder, ported from raw.py:449-466. Returns idx."""
    ee = []
    for j in range(N, len(vol) - N - 1):
        is_min = True
        for k in range(j - N, j + N + 1):
            if k != j and vol[k] <= vol[j]:
                is_min = False; break
        if is_min:
            ee.append(j)
    return np.array(ee)

# ---- same synthetic navigator as before ----
rng = np.random.default_rng(0)
TR = 0.02226
t = np.arange(0, 18.0, TR)
period = 3.2
amp = 1.0 + 0.25 * np.sin(2 * np.pi * t / 11.0)
base = 0.15 * np.sin(2 * np.pi * t / 15.0)
z = base + amp * (0.5 - 0.5 * np.cos(2 * np.pi * t / period))
z += 0.01 * rng.standard_normal(len(t))
z = (z - z.min()) / (z.max() - z.min())

nbins = 16
bindist0sq = 2.0
ilvbin, breathdir = asap_bin(z)
binctr = ilvbin * nbins

# ---- per-interleave within-breath time ----
ee = find_EE(z, N=30)
ee_t = t[ee]
# fractional breath phase tau in [0,1): 0 at preceding EE, 1 at next EE
tau = np.full(len(t), np.nan)
sec_since_ee = np.full(len(t), np.nan)     # absolute seconds since preceding EE
for k in range(len(ee) - 1):
    m = (t >= ee_t[k]) & (t < ee_t[k + 1])
    tau[m] = (t[m] - ee_t[k]) / (ee_t[k + 1] - ee_t[k])
    sec_since_ee[m] = t[m] - ee_t[k]
mean_period = np.mean(np.diff(ee_t))

# ---- soft-weighted mean time per output bin (EXACT recon weighting) ----
def cyclic_bindist(ibin, bc):
    return np.minimum.reduce([np.abs(ibin - bc),
                              np.abs(ibin + nbins - bc),
                              np.abs(ibin - nbins - bc)])
valid = ~np.isnan(tau) & (binctr >= 0)
bc = binctr[valid]; tv = tau[valid]
# phase is CYCLIC (wraps at EE) and the bin weight is cyclic too -> a linear
# mean corrupts bins near the wrap. Use the CIRCULAR mean (average the angle).
ang = 2 * np.pi * tv
bin_phase = np.zeros(nbins); bin_phase_std = np.zeros(nbins)
for b in range(nbins):
    w = np.exp(-cyclic_bindist(b, bc) ** 2 / bindist0sq)
    sw = np.sum(w * np.sin(ang)); cw = np.sum(w * np.cos(ang)); W = np.sum(w)
    mean_ang = np.arctan2(sw, cw)
    bin_phase[b] = (mean_ang / (2 * np.pi)) % 1.0          # in [0,1)
    R = np.sqrt(sw ** 2 + cw ** 2) / W                      # resultant length
    bin_phase_std[b] = np.sqrt(max(-2 * np.log(R), 0)) / (2 * np.pi)  # circ std (cycles)
# unwrap for display: bins are phase-ordered, so make phase monotonic across bins
bin_phase_unwrap = np.unwrap(bin_phase * 2 * np.pi) / (2 * np.pi)
bin_phase_unwrap -= bin_phase_unwrap[0]                    # start at 0 (relative clock)
bin_time_s = bin_phase_unwrap * mean_period                # representative-breath seconds

# ============================================================== FIGURE
fig = plt.figure(figsize=(15, 8.5))
gs = GridSpec(2, 2, figure=fig, hspace=0.34, wspace=0.24)
C_RISE, C_FALL = '#2c7fb8', '#d95f0e'

# A: navigator with detected EE troughs -> defines breath phase
axA = fig.add_subplot(gs[0, 0])
axA.plot(t, z, '-', color='0.7', lw=0.8)
axA.scatter(t, z, s=5, c=np.where(breathdir == 1, C_FALL, C_RISE))
axA.scatter(ee_t, z[ee], marker='v', s=70, c='k', zorder=5, label='detected EE (phase=0)')
axA.set_title('A.  Navigator + end-expiration troughs\nphase tau = fraction of breath since preceding EE', fontsize=11)
axA.set_xlabel('time (s)'); axA.set_ylabel('displacement z'); axA.legend(fontsize=8)

# B: per-interleave amplitude-bin vs its within-breath phase (the raw scatter)
axB = fig.add_subplot(gs[0, 1])
axB.scatter(tau[valid], binctr[valid], s=8,
            c=np.where(breathdir[valid] == 1, C_FALL, C_RISE))
axB.set_title('B.  Each interleave: within-breath phase  ->  amplitude bin\nmonotonic but NOT proportional', fontsize=11)
axB.set_xlabel('within-breath phase tau  (0=EE, ~0.4=EI, 1=next EE)')
axB.set_ylabel('bin center (ilvbin x nbins)')

# C: THE ANSWER -- bin index vs its real averaged time
axC = fig.add_subplot(gs[1, 0])
axC.errorbar(range(nbins), bin_phase_unwrap, yerr=bin_phase_std, fmt='o-', color='#238b45',
             capsize=3, label='circular-mean phase +/- circ std')
# uniform-time reference (what you'd WRONGLY assume if bin idx == time)
axC.plot(range(nbins), np.linspace(bin_phase_unwrap[0], bin_phase_unwrap[-1], nbins), '--', color='0.6',
         label='naive: bin index proportional to time')
axC.set_title('C.  REAL averaged time of each bin\n(this is the 4D time vector)', fontsize=11)
axC.set_xlabel('bin index'); axC.set_ylabel('mean within-breath phase (cycles)')
axC.legend(fontsize=8, loc='upper left')
for b in range(nbins):
    axC.annotate('%.2fs' % bin_time_s[b], (b, bin_phase[b]), fontsize=6.5,
                 textcoords='offset points', xytext=(0, 7), ha='center', color='0.3')

# D: bin spacing in time -> shows non-uniform sampling of the breath
axD = fig.add_subplot(gs[1, 1])
dt = np.diff(bin_time_s)
axD.bar(range(1, nbins), dt, color='#41b6c4', edgecolor='white')
axD.axhline(mean_period / nbins, color='r', ls='--', lw=1,
            label='uniform spacing (period/nbins = %.3fs)' % (mean_period / nbins))
axD.set_title('D.  Real time gap between consecutive bins\nbig gaps = fast diaphragm; small = dwell near EE/EI', fontsize=11)
axD.set_xlabel('bin transition (b-1 -> b)'); axD.set_ylabel('delta t (s)')
axD.legend(fontsize=8)

fig.suptitle('DIAPHRAGM 4D time vector: bins are amplitude-rank -> must AVERAGE within-breath time per bin (period=%.2fs)' % mean_period,
             fontsize=12.5, y=0.99)
import os
out = os.path.join(os.path.dirname(__file__), '..', '..', 'outputs', 'diaphragm_binning', 'diaphragm_bintime.png')
fig.savefig(out, dpi=130, bbox_inches='tight')
print('saved', out)
print('mean breath period (s):', round(mean_period, 3), ' nEE:', len(ee))
print('bin :  phase(cyc)   time_s')
for b in range(nbins):
    print('%3d :  %.3f       %.3f' % (b, bin_phase_unwrap[b], bin_time_s[b]))
