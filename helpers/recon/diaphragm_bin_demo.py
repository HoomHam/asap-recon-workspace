"""
Illustrate exactly how ASAP's DIAPHRAGM binning maps interleaves -> one
representative breath cycle. Reproduces raw.py:153 bin() verbatim and applies
it to a synthetic diaphragm-displacement navigator curve.

Answer to the question: placement is by AMPLITUDE RANK, not time. Time enters
only through the sign of the slope (rising vs falling) which splits the breath
into an inhale limb and an exhale limb so the same displacement value gets two
different bins.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# ---------------------------------------------------------------- exact bin()
def asap_bin(val):
    """Verbatim port of raw.py:153-177 bin()."""
    bins = val.copy()
    v = np.array([x for x in val if not np.isnan(x)])
    if len(v) > 2:
        breathdir = np.zeros(len(v))
        if v[1] < v[0]:
            breathdir[0] = 1
        for i in range(1, len(v)):
            if v[i] < v[i - 1]:
                breathdir[i] = 1              # 1 = falling (this sample)
        iv = np.sort(v[breathdir == 1])
        ev = np.sort(v[breathdir == 0])
        maxiv = 0.0 if len(iv) == 0 else max(iv) * 1.00001
        maxev = 0.0 if len(ev) == 0 else max(ev)
        v = v.copy()
        v[breathdir == 1] = maxiv + maxev - v[breathdir == 1]   # FOLD falling limb up
        vsort = np.sort(v)
        vidx = 0
        for i in range(len(bins)):
            if np.isnan(bins[i]):
                continue
            bins[i] = np.searchsorted(vsort, v[vidx])           # RANK = bin
            vidx += 1
        bins[np.isnan(bins)] = -1
        bins /= len(v)
    return bins, breathdir

# ---------------------------------------------------- synthetic navigator curve
# Diaphragm z displacement over several real breaths of varying amplitude/period,
# sampled once per interleave (per gas TR). Rescaled to [0,1] like ilvvol.
rng = np.random.default_rng(0)
TR = 0.02226                      # gas interleave spacing (s), ASAP typical
t = np.arange(0, 18.0, TR)        # ~18 s of acquisition
# breathing: sum of a base sinusoid with slowly drifting amplitude + baseline wander
period = 3.2
amp = 1.0 + 0.25 * np.sin(2 * np.pi * t / 11.0)         # amplitude varies breath-to-breath
base = 0.15 * np.sin(2 * np.pi * t / 15.0)              # slow baseline drift
z = base + amp * (0.5 - 0.5 * np.cos(2 * np.pi * t / period))   # 0=EE (up), peak=EI (down)
z += 0.01 * rng.standard_normal(len(t))                # navigator noise
# rescale to [0,1] (raw.rescale)
z = (z - z.min()) / (z.max() - z.min())

ilvbin, breathdir = asap_bin(z)   # per-interleave fractional bin in [0,1)
nbins = 16
binctr = ilvbin * nbins           # results.py:248  binarr = bins * nbins
rising = breathdir == 0           # exhale->inhale (diaphragm descending, z rising)
falling = breathdir == 1

# label each interleave with which real breath it belongs to (for the pooling panel)
breath_id = np.floor(t / period).astype(int)

# ============================================================== FIGURE
fig = plt.figure(figsize=(15, 9))
gs = GridSpec(2, 2, figure=fig, hspace=0.32, wspace=0.22)
C_RISE, C_FALL = '#2c7fb8', '#d95f0e'

# --- Panel A: navigator displacement vs time, colored by slope sign ----------
axA = fig.add_subplot(gs[0, 0])
axA.plot(t, z, '-', color='0.7', lw=0.8, zorder=1)
axA.scatter(t[rising], z[rising], s=9, c=C_RISE, label='rising  (breathdir=0, inhale limb)', zorder=2)
axA.scatter(t[falling], z[falling], s=9, c=C_FALL, label='falling (breathdir=1, exhale limb)', zorder=2)
axA.set_title('A.  Navigator: diaphragm displacement z(t)\nEACH DOT = one interleave; color = sign of slope', fontsize=11)
axA.set_xlabel('time (s)'); axA.set_ylabel('displacement (0=EE, 1=EI)')
axA.legend(fontsize=8, loc='upper right')
# annotate: same amplitude, two limbs
axA.annotate('same z, two limbs\n-> two different bins', xy=(t[len(t)//3], 0.5),
             xytext=(2.0, 0.82), fontsize=8, color='k',
             arrowprops=dict(arrowstyle='->', color='0.4'))

# --- Panel B: the FOLD -------------------------------------------------------
# original amplitude on x; folded value on y. Falling limb lifted above rising.
axB = fig.add_subplot(gs[0, 1])
v = z.copy()
iv = np.sort(v[falling]); ev = np.sort(v[rising])
maxiv = max(iv) * 1.00001; maxev = max(ev)
vfold = v.copy(); vfold[falling] = maxiv + maxev - v[falling]
axB.scatter(v[rising], vfold[rising], s=9, c=C_RISE, label='rising limb: value unchanged')
axB.scatter(v[falling], vfold[falling], s=9, c=C_FALL, label='falling limb: v -> (maxiv+maxev) - v')
axB.axhline(maxev, color='0.5', ls='--', lw=0.8)
axB.text(0.02, maxev + 0.03, 'maxev  (top of exhale limb)', fontsize=7.5, color='0.4')
axB.set_title('B.  The FOLD: hysteresis loop -> single monotonic axis\ny = value that gets rank-sorted', fontsize=11)
axB.set_xlabel('displacement z (original amplitude)')
axB.set_ylabel('folded value (sorted by this)')
axB.legend(fontsize=8, loc='lower center')

# --- Panel C: representative single breath cycle -----------------------------
# amplitude vs bin center. This IS the representative cycle: exhale limb fills
# low bins, inhale limb fills high bins, wrapping cyclically.
axC = fig.add_subplot(gs[1, 0])
axC.scatter(binctr[rising], z[rising], s=10, c=C_RISE, label='inhale limb -> low bins')
axC.scatter(binctr[falling], z[falling], s=10, c=C_FALL, label='exhale limb -> high bins')
axC.set_title('C.  Representative breath cycle (all breaths pooled)\nx = bin center (0..%d), y = displacement' % nbins, fontsize=11)
axC.set_xlabel('bin index  (ilvbin x nbins)'); axC.set_ylabel('displacement z')
axC.set_xticks(range(0, nbins + 1, 2))
axC.legend(fontsize=8, loc='center right')
axC.annotate('cyclic: bin 0 and bin %d\nare neighbors (wrap-around)' % (nbins - 1),
             xy=(0.5, 0.02), xytext=(4.5, 0.12), fontsize=8,
             arrowprops=dict(arrowstyle='->', color='0.4'))

# --- Panel D: equal-count property + pooling ---------------------------------
axD = fig.add_subplot(gs[1, 1])
# histogram of interleaves per integer bin, and show breath-of-origin mix
edges = np.arange(0, nbins + 1)
disc = np.clip(np.floor(binctr[binctr >= 0]).astype(int), 0, nbins - 1)
ids = breath_id[binctr >= 0]
nbreaths = ids.max() + 1
cmap = plt.cm.viridis(np.linspace(0, 1, nbreaths))
bottom = np.zeros(nbins)
for b in range(nbreaths):
    counts = np.array([np.sum((disc == k) & (ids == b)) for k in range(nbins)])
    axD.bar(range(nbins), counts, bottom=bottom, color=cmap[b], width=0.9,
            edgecolor='white', linewidth=0.3)
    bottom += counts
axD.set_title('D.  Equal-count binning + pooling\nbar height = # interleaves; color = which real breath', fontsize=11)
axD.set_xlabel('discrete bin index'); axD.set_ylabel('# interleaves')
axD.set_xticks(range(0, nbins, 2))
axD.text(0.5, 0.95, 'each bin ~equal population;\nfilled from MANY real breaths',
         transform=axD.transAxes, fontsize=8, va='top',
         bbox=dict(boxstyle='round', fc='white', ec='0.7'))

fig.suptitle('ASAP DIAPHRAGM binning: amplitude-rank -> representative breath cycle  (raw.py:153 bin())',
             fontsize=13, y=0.99)
import os
out = os.path.join(os.path.dirname(__file__), '..', '..', 'outputs', 'diaphragm_binning', 'diaphragm_binning.png')
fig.savefig(out, dpi=130, bbox_inches='tight')
print('saved', out)
print('nbins', nbins, 'ninterleaves', len(t),
      'rising', rising.sum(), 'falling', falling.sum())
