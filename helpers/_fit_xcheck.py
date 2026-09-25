#!/usr/bin/env python3
"""Scratch (2026-09-24, B23): cross-check the RBC/TP split inputs from Steve's 2-Lorentzian spectrum
fit (resplit_summary.csv, via raw.load_from_arr) against the XeCS calspec time-domain fits
(2026_XeCS_Recon workspace/outputs/calspec/cohort/subjects_pooled.csv: M2 = gas L / RBC L / mem Voigt,
M3 = gas L / RBC L / mem1 L + mem2 L). Everything is moved to the IMAGE k0 time:
  Steve: phases at his first spectral sample (ADC sample killpts=2, t = TE + 120 us) -> k0 by -2*pi*df*110 us
  XeCS : phases at ADC sample 0 (t = TE; fit_fid uses t = arange(t0, N)*dwell)   -> k0 by +2*pi*df*10 us
For M3 the membrane is LUMPED as the complex sum mem1 + mem2 evaluated at k0 (the image split sees one
non-RBC component), so dphi = angle(RBC / (mem1 + mem2)) and ratio = a_RBC / |mem1 + mem2|.
Writes outputs/resplit_2026-09-24/fit_xcheck.csv and prints cohort statistics.
"""
import csv
import numpy as np
from pathlib import Path

HZ = 17.61                                   # Hz per ppm at 17.612 MHz (XeCS calspec)
ASAP = Path(__file__).resolve().parents[1]
XECS = ASAP.parents[1] / '2026_XeCS_Recon' / 'workspace' / 'outputs' / 'calspec' / 'cohort' / 'subjects_pooled.csv'
A = {r['key']: r for r in csv.DictReader(open(ASAP / 'outputs' / 'resplit_2026-09-24' / 'resplit_summary.csv'))}
X = list(csv.DictReader(open(XECS)))
w = lambda d: (d + 180) % 360 - 180
f = lambda r, k: float(r[k]) if r.get(k) not in (None, '', 'nan') else np.nan
T_K0 = 10e-6                                 # image k0 after ADC start (killpts * dtdyn)
T_STEVE = 120e-6                             # Steve's first spectral sample after ADC start


def comp(x, p, t):
    """M3 component as a complex number at time t (from ADC sample 0)."""
    return f(x, f'M3_{p}_a') * np.exp(1j * np.radians(f(x, f'M3_{p}_phi_deg'))) * \
        np.exp(2j * np.pi * f(x, f'M3_{p}_ppm_from_gas') * HZ * t) * np.exp(-np.pi * f(x, f'M3_{p}_wL_hz') * t)


rows = []
for x in X:
    a = A.get(x['session'])
    if not a:
        continue
    df_x = (f(x, 'rbc_ppm_from_gas') - f(x, 'mem_ppm_from_gas')) * HZ
    d2 = w(f(x, 'rbc_phi_deg') - f(x, 'mem_phi_deg'))
    d2_k0 = w(d2 + 360 * df_x * T_K0)
    mem = comp(x, 'mem1', T_K0) + comp(x, 'mem2', T_K0)
    rbc = comp(x, 'rbc', T_K0)
    d3_k0 = w(np.degrees(np.angle(rbc / mem)))
    r3_lumped = abs(rbc) / abs(mem)
    r3_scalar = f(x, 'M3_ratio_rbc_mem')
    row = dict(session=x['session'], seq=x['seq'], xecs_valid=x['valid'], snr_diss=f(x, 'snr_diss'),
               asap_status=a['status'], asap_gate=a['gate_loose'], cond=a.get('cond_new', ''),
               df_steve=f(a, 'df'), df_xecs=df_x,
               dphi_k0_steve=f(a, 'dphi_new'), dphi_k0_m2=d2_k0, dphi_k0_m3lumped=d3_k0,
               dphi_steve1st_minus_m2at120=w(f(a, 'dphi_new') - f(a, 'killpts_corr_deg') - w(d2 + 360 * df_x * T_STEVE)),
               m3_mem1_minus_mem2_deg=w(f(x, 'M3_mem1_phi_deg') - f(x, 'M3_mem2_phi_deg')),
               m3_a1_over_a2=f(x, 'M3_mem1_a') / f(x, 'M3_mem2_a'),
               ratio_steve=f(a, 'target'), ratio_m2=f(x, 'ratio_rbc_mem'), ratio_m3_scalar=r3_scalar,
               ratio_m3_lumped=r3_lumped, bic_m2=f(x, 'bic'), bic_m3=f(x, 'M3_bic'))
    rows.append(row)

out = ASAP / 'outputs' / 'resplit_2026-09-24' / 'fit_xcheck.csv'
with open(out, 'w', newline='') as fh:
    wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    wr.writeheader()
    for r in rows:
        wr.writerow({k: (f'{v:.4g}' if isinstance(v, float) else v) for k, v in r.items()})
print(f'-> {out}  ({len(rows)} sessions in both tables)')

ok = [r for r in rows if r['asap_status'] == 'resplit']


def st(name, v):
    v = np.array(v, float); v = v[np.isfinite(v)]
    print(f'{name:58s} n={len(v):2d} median {np.median(v):+7.2f}  MAD {np.median(np.abs(v - np.median(v))):5.2f}  '
          f'IQR {np.percentile(v, 25):+.2f}..{np.percentile(v, 75):+.2f}')


st('Steve(1st sample) - M2 moved to 120 us  [deg]', [r['dphi_steve1st_minus_m2at120'] for r in ok])
st('Steve k0 - M2 k0  [deg]', [w(r['dphi_k0_steve'] - r['dphi_k0_m2']) for r in ok])
st('Steve k0 - M3 lumped k0  [deg]', [w(r['dphi_k0_steve'] - r['dphi_k0_m3lumped']) for r in ok])
st('M3 mem1 - mem2 phase  [deg]', [r['m3_mem1_minus_mem2_deg'] for r in ok])
st('df XeCS - df Steve  [Hz]', [r['df_xecs'] - r['df_steve'] for r in ok])
st('ratio M2 / Steve', [r['ratio_m2'] / r['ratio_steve'] for r in ok])
st('ratio M3 scalar (a_rbc/(a1+a2)) / Steve', [r['ratio_m3_scalar'] / r['ratio_steve'] for r in ok])
st('ratio M3 LUMPED (a_rbc/|mem1+mem2|) / Steve', [r['ratio_m3_lumped'] / r['ratio_steve'] for r in ok])
print('BIC prefers M3 in', sum(r['bic_m3'] < r['bic_m2'] for r in ok), 'of', len(ok))
bad = [(r['session'], round(w(r['dphi_k0_steve'] - r['dphi_k0_m3lumped']), 1)) for r in ok
       if abs(w(r['dphi_k0_steve'] - r['dphi_k0_m3lumped'])) > 20]
print('|Steve - M3 lumped| > 20 deg:', bad)
gated = [(r['session'], r['asap_gate'], r['xecs_valid'], round(r['dphi_k0_m3lumped'], 1), round(r['ratio_m3_lumped'], 3))
         for r in rows if r['asap_status'] != 'resplit']
print('ASAP-gated sessions present in XeCS table (gate, xecs_valid, M3 lumped dphi_k0, ratio):')
for g in gated:
    print('  ', g)
