#!/usr/bin/env python
"""TE90 per dynamic session from the XeCS spectral fits (no new fitting).

phi(t) = phi_k0 + 2*pi*df*(t - t_k0)   [RBC - TP, advancing with t; F59/F61]
TE90   = t_k0 + (90 deg - phi_k0) / (360 deg * df)

Inputs
  outputs/resplit_2026-09-24/resplit_summary_xecs.csv   (dphi_new = phi_k0 [deg], df_xecs [Hz], carrier, gates)
  notes/calspec_package_2026-09-16/resplit_inputs_2026-09-24.csv  (t_k0_us, TE_us per session -- the time phi_k0 refers to)
Outputs
  outputs/te90_2026-09-25/te90_table.csv
  outputs/te90_2026-09-25/te90_strip.png
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

WS = Path(__file__).resolve().parents[1]
OUT = WS / 'outputs' / 'te90_2026-09-25'
OUT.mkdir(parents=True, exist_ok=True)

s = pd.read_csv(WS / 'outputs/resplit_2026-09-24/resplit_summary_xecs.csv')
s = s[~s.key.str.contains('_merged|_REJECTED')]
i = pd.read_csv(WS / 'notes/calspec_package_2026-09-16/resplit_inputs_2026-09-24.csv')
i = i[i.primary].drop_duplicates('session')[['session', 'TE_us', 't_k0_us', 'rbc_ppm', 'mem1_ppm', 'mem2_ppm']]

d = s.merge(i, left_on='key', right_on='session', how='left')
d = d[d.dphi_new.notna()].copy()
d['phi_k0_deg'] = d.dphi_new
d['df_hz'] = d.df_xecs
d['t_k0_us'] = d.t_k0_us.fillna(d.TE_ms * 1000 + 10)
d['deg_per_us'] = 360.0 * d.df_hz * 1e-6
d['TE90_us'] = d.t_k0_us + (90.0 - d.phi_k0_deg) / d.deg_per_us
d['dTE_us'] = d.TE90_us - d.TE_us                  # how much later than the acquired TE
d['phi_at_TE_deg'] = d.phi_k0_deg - (d.t_k0_us - d.TE_us) * d.deg_per_us
d['TE90_sd_us'] = d.dphi_rep_sd / d.deg_per_us     # rep-to-rep spread of phi_k0 -> TE90 spread
d['year'] = d.key.str[:4]
d['atlas_set'] = (d.carrier.eq('ON_RBC') & d.status.eq('resplit')).map({True: 'yes', False: 'no'})

cols = ['key', 'year', 'carrier', 'carrier_ppm', 'status', 'xecs_stable', 'atlas_set',
        'TE_us', 't_k0_us', 'df_hz', 'phi_k0_deg', 'dphi_rep_sd', 'phi_at_TE_deg',
        'TE90_us', 'TE90_sd_us', 'dTE_us', 'F_lump', 'ratio_scalar', 'snr_diss']
t = d[cols].sort_values('key').round(1)
t.to_csv(OUT / 'te90_table.csv', index=False)

# ---- strip plot: TE90 per session, sorted by date; carrier + gate coded ----
t = t.reset_index(drop=True)
fig, axs = plt.subplots(2, 1, figsize=(16, 8), sharex=True,
                        gridspec_kw={'height_ratios': [2, 1]})
colors = {'ON_RBC': 'tab:blue', 'BETWEEN': 'tab:orange', 'ON_MEM': 'tab:red',
          'BELOW_MEM': 'tab:purple', 'ABOVE_RBC': 'tab:green'}
x = np.arange(len(t))
for car, c in colors.items():
    m = t.carrier.eq(car)
    if not m.any():
        continue
    gated = m & t.status.ne('resplit')
    ok = m & t.status.eq('resplit')
    axs[0].errorbar(x[ok], t.TE90_us[ok], yerr=t.TE90_sd_us[ok], fmt='o', color=c,
                    ms=5, capsize=2, lw=0.8, label=f'{car} (n={int(m.sum())})')
    axs[0].errorbar(x[gated], t.TE90_us[gated], yerr=t.TE90_sd_us[gated], fmt='x',
                    color=c, ms=7, capsize=2, lw=0.8, mfc='none')
    axs[1].scatter(x[m], t.phi_k0_deg[m], color=c, s=18,
                   marker='o' if True else 'x')
axs[0].plot([], [], 'kx', label='gated (unstable / invalid)')
for ax in axs:
    ax.axvline(x[t.year.eq('2024')].min() - 0.5, color='grey', lw=0.8, ls='--')
axs[0].axhline(t.TE_us.median(), color='grey', lw=0.8, ls=':')
axs[0].text(0.3, t.TE_us.median() + 8, 'acquired TE', fontsize=8, color='grey')
axs[0].set_ylabel('TE90 [µs]')
med = t[t.atlas_set.eq('yes')].groupby('year').TE90_us.median()
axs[0].set_title('TE90 = t_k0 + (90° − φ_k0)/(360·Δf), from XeCS cal-block fits · '
                 + ' · '.join(f'{y} median {v:.0f} µs' for y, v in med.items())
                 + f' · dashed = 2023|2024', fontsize=10)
axs[0].legend(fontsize=8, ncol=3, loc='upper left')
axs[1].axhline(90, color='grey', lw=0.8, ls=':')
axs[1].set_ylabel('φ_k0 = RBC−TP at k0 [deg]')
axs[1].set_xticks(x)
axs[1].set_xticklabels(t.key, rotation=90, fontsize=5.5)
axs[1].set_xlim(-1, len(t))
fig.tight_layout()
fig.savefig(OUT / 'te90_strip.png', dpi=150)
print(t.to_string(index=False))
print('\natlas-set medians:', med.round(0).to_dict())
print('written', OUT)
