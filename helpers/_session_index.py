"""Scratch (2026-10-26, Hooman: "≥ 50 datasets must show the same stats before F/E′"): one row per session from the
bin_timing exports (k0 per interleave, Tyger EE list, auto exclusions, per-breath A/B/C with nav), with the roster
labels and the per-session INJECTION PROFILE (Hooman: bolus delivery may differ per study).

Columns (outputs/bin_timing_2026-10-26/session_index.csv):
  key, sid, date, cohort (cohort_map.csv), protocol (TR: 14.9 ms v2 / 17.2 ms v3 / other), nilv, n_EE, scan_s,
  T_med_s, T_cv, frac_long_breaths (> 1.5 x median), irregular_flag (period IQR/median > 0.3 or any long breath),
  n_excl_ilv, frac_excl, excl_midscan (an excluded run that ends > 20 s before the last EE),
  tB_med_s, tB_iqr_s, tB_sd_s, tB_over_T, slope_tB_vs_T (s/s, breaths kept), r_tB_T,
  tC_med_s, tC_over_T, slope_tC_vs_T, slope_exh_vs_T (= d(T - tC)/dT), r_exh_T,
  tau_I_0375 (= 0.375 T_med), frac_tB_inside_tauI,
  rise_t50_s / rise_t90_s (time from EE to 50 % / 90 % of (B - A) on the per-breath k0), early_frac_0p5s (k0 at
  EE + 0.5 s minus A over B - A), early_frac_0p25s, BA_smsm_med, CB_smsm_med, k0_trend_pct_per_min (log k0_B vs time),
  nav_ok (nav max available on >= 80 % of breaths).
usage: python helpers/_session_index.py            (all keys with <key>_breaths_raw.csv present)
"""
import csv
import glob
import os
import sys

import numpy as np

BT = os.path.join(os.path.dirname(__file__), '..', 'outputs', 'bin_timing_2026-10-26')
ROSTER = os.path.join(os.path.dirname(__file__), '..', 'notes', 'calspec_package_2026-09-16', 'cohort_map.csv')
SESSTAB = os.path.join(os.path.dirname(__file__), '..', 'notes', 'calspec_package_2026-09-16', 'session_table.csv')
BATCHLOG = os.path.join(BT, '_batch_all108.log')


def rd(p):
    return list(csv.DictReader(open(p)))


def col(rows, name):
    return np.array([float(r[name]) if r[name] not in ('', 'nan') else np.nan for r in rows])


def slope(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 5:
        return np.nan, np.nan
    p = np.polyfit(x[m], y[m], 1)
    r = np.corrcoef(x[m], y[m])[0, 1]
    return p[0], r


def one(key, cohort, seq_roster='', n_EE_log=np.nan):
    ilv = rd(f'{BT}/{key}_ilv_raw.csv'); t = col(ilv, 't_s'); k0 = col(ilv, 'k0_gas_sum8')
    ty = rd(f'{BT}/{key}_ilvbin_tyger.csv'); excl = col(ty, 'excluded').astype(int)
    ee = np.atleast_1d(np.loadtxt(f'{BT}/{key}_ee_tyger.csv', skiprows=1))
    br = rd(f'{BT}/{key}_breaths_raw.csv')
    TR = float(np.median(np.diff(t)))
    proto = 'v2' if abs(TR - 0.0149) < 5e-4 else ('v3' if abs(TR - 0.0172) < 5e-4 else f'TR{TR * 1e3:.1f}')
    per = col(br, 'period_s'); tB = col(br, 't_B_after_EE_s'); tC = col(br, 't_C_after_EE_s'); tEE = col(br, 't_EE_s')
    Tmed = float(np.median(per)); keep = per < 1.5 * Tmed
    # exclusion runs
    runs = []
    i = 0
    while i < len(excl):
        if excl[i]:
            j = i
            while j < len(excl) and excl[j]:
                j += 1
            runs.append((t[i], t[j - 1])); i = j
        else:
            i += 1
    midscan = any(r[1] < ee[-1] - 20 for r in runs) if len(ee) else False
    # breath-level kept set (EE inside an excluded run dropped, long breaths dropped)
    exEE = np.array([bool(excl[int(np.argmin(np.abs(t - te)))]) for te in tEE])
    kb = keep & ~exEE & np.isfinite(tB)
    # injection profile per breath from the raw k0 trace
    t50, t90, e05, e025 = [], [], [], []
    for k in range(len(ee) - 1):
        if k >= len(br) or not kb[k]:
            continue
        m = (t >= ee[k]) & (t < ee[k + 1])
        tt, ss = t[m] - ee[k], k0[m]
        if len(tt) < 10:
            continue
        A = ss[:3].mean(); iB = int(np.argmax(ss)); B = ss[iB]
        if B <= A:
            continue
        f = (ss[:iB + 1] - A) / (B - A)
        t50.append(tt[int(np.argmax(f >= 0.5))]); t90.append(tt[int(np.argmax(f >= 0.9))])
        e05.append(np.interp(0.5, tt, (ss - A) / (B - A))); e025.append(np.interp(0.25, tt, (ss - A) / (B - A)))
    sl_tB, r_tB = slope(per[kb], tB[kb]); sl_tC, r_tC = slope(per[kb], tC[kb]); sl_ex, r_ex = slope(per[kb], per[kb] - tC[kb])
    tau = 0.375 * Tmed
    ba = col(br, 'B_over_A_smsm'); cb = col(br, 'C_over_B_smsm'); kB = col(br, 'k0_B_sm') if 'k0_B_sm' in br[0] else col(br, 'k0_B')
    tr_m = kb & np.isfinite(kB) & (kB > 0)
    trend = np.nan
    if tr_m.sum() >= 5:
        trend = 100 * 60 * np.polyfit(tEE[tr_m], np.log(kB[tr_m]), 1)[0]
    q = lambda v: (np.nanmedian(v) if np.isfinite(v).any() else np.nan)
    iqr = lambda v: (np.nanpercentile(v, 75) - np.nanpercentile(v, 25)) if np.isfinite(v).sum() > 3 else np.nan
    return dict(key=key, sid=key.split('_')[1], date=key.split('_')[0], cohort=cohort, protocol=proto, seq_roster=seq_roster,
                TR_ms=round(TR * 1e3, 2), nilv=len(t), n_EE=len(ee), n_EE_tyger_log=n_EE_log, scan_s=round(float(t[-1]), 1), T_med_s=round(Tmed, 2),
                T_cv=round(float(np.std(per[keep]) / np.mean(per[keep])), 3) if keep.sum() > 2 else np.nan,
                frac_long_breaths=round(float(np.mean(~keep)), 3),
                irregular_flag=int((iqr(per) / Tmed > 0.3) or (~keep).any()),
                n_excl_ilv=int(excl.sum()), frac_excl=round(float(excl.mean()), 3), excl_midscan=int(midscan),
                n_breaths_kept=int(kb.sum()),
                tB_med_s=round(q(tB[kb]), 3), tB_iqr_s=round(iqr(tB[kb]), 3), tB_sd_s=round(float(np.nanstd(tB[kb])), 3),
                tB_over_T=round(q(tB[kb] / per[kb]), 3), slope_tB_vs_T=round(sl_tB, 3), r_tB_T=round(r_tB, 2),
                tC_med_s=round(q(tC[kb]), 3), tC_over_T=round(q(tC[kb] / per[kb]), 3), slope_tC_vs_T=round(sl_tC, 3),
                slope_exh_vs_T=round(sl_ex, 3), r_exh_T=round(r_ex, 2),
                tau_I_0375=round(tau, 3), frac_tB_inside_tauI=round(float(np.mean(tB[kb] < tau)), 3) if kb.any() else np.nan,
                rise_t50_s=round(q(np.array(t50)), 3), rise_t90_s=round(q(np.array(t90)), 3),
                early_frac_0p5s=round(q(np.array(e05)), 3), early_frac_0p25s=round(q(np.array(e025)), 3),
                BA_smsm_med=round(q(ba[kb]), 2), CB_smsm_med=round(q(cb[kb]), 3), k0_trend_pct_per_min=round(trend, 1),
                nav_ok=int(np.mean(np.isfinite(tC[kb])) >= 0.8) if kb.any() else 0)


def main():
    roster = {r['sid']: r['cohort'] for r in rd(ROSTER)}
    seqtab = {r['session']: r['seq'] for r in rd(SESSTAB)} if os.path.exists(SESSTAB) else {}
    eelog = {}; dg = {}
    if os.path.exists(BATCHLOG):
        import re
        for line in open(BATCHLOG):
            m = re.match(r'^(\S+): binning=\S+ nilv \d+ TR [\d.]+ EE (\d+)', line)
            if m:
                eelog[m.group(1)] = int(m.group(2))
            m = re.match(r'^(\S+): nav .*diss/gas sum8 at A ([\d.]+)', line)
            if m:
                dg.setdefault(m.group(1), float(m.group(2)))
    keys = sorted(os.path.basename(p)[:-len('_breaths_raw.csv')] for p in glob.glob(f'{BT}/*_breaths_raw.csv'))
    rows, fails = [], []
    for k in keys:
        try:
            r = one(k, roster.get(k.split('_')[1], 'UNK'), seqtab.get(k, ''), eelog.get(k, np.nan))
            r['diss_over_gas_at_A'] = dg.get(k, np.nan)
            # layout_ok: v2/v3 TR AND the gas/diss line slicing is sane (dissolved k0 << gas k0 at EE); the TR 22.3 /
            # 25.9 ms sequences and the mixed-ADC sessions fail this -> keep out of every pooled statistic
            r['layout_ok'] = int(r['protocol'] in ('v2', 'v3') and dg.get(k, 9) < 0.6)
            rows.append(r)
        except Exception as e:
            fails.append((k, repr(e)))
    with open(f'{BT}/session_index.csv', 'w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print(f'{len(rows)} sessions -> {BT}/session_index.csv; failed {len(fails)}')
    for k, e in fails:
        print('  FAIL', k, e[:150])
    # compact summary by protocol x cohort
    import collections
    g = collections.defaultdict(list)
    for r in rows:
        if r['layout_ok']:
            g[(r['protocol'], r['cohort'])].append(r)
    print(f'layout_ok sessions: {sum(r["layout_ok"] for r in rows)} of {len(rows)} (v3 {sum(1 for r in rows if r["layout_ok"] and r["protocol"]=="v3")}, v2 {sum(1 for r in rows if r["layout_ok"] and r["protocol"]=="v2")})')
    print(f'{"proto":6s} {"cohort":7s} n | tB_med  tB/T  slope_tB  | tC/T  slope_exh | tB<tauI | t50  e0.5 | BA   | irreg midexcl')
    for (p, c), rs in sorted(g.items()):
        q = lambda f: np.nanmedian([r[f] for r in rs if r[f] == r[f]])
        print(f'{p:6s} {c:7s} {len(rs):2d} | {q("tB_med_s"):5.2f} {q("tB_over_T"):5.2f} {q("slope_tB_vs_T"):+6.2f}   | {q("tC_over_T"):5.2f} {q("slope_exh_vs_T"):+6.2f}  | {q("frac_tB_inside_tauI"):5.2f}   | {q("rise_t50_s"):4.2f} {q("early_frac_0p5s"):4.2f} | {q("BA_smsm_med"):4.2f} | {sum(r["irregular_flag"] for r in rs):2d}   {sum(r["excl_midscan"] for r in rs):2d}')


if __name__ == '__main__':
    main()
