"""Scratch (2026-10-22): outputs/snr_2026-10-22/snr_sessions.csv -> snr_sessions.html (local, sortable, no external
assets; subject codes stay off the web). Summary on top: sessions, median gas/DP SNR by year and by recon tree."""
import csv
import html
import statistics as st
from collections import defaultdict
from pathlib import Path

OUT = Path(__file__).resolve().parents[1] / 'outputs' / 'snr_2026-10-22'
rows = list(csv.DictReader(open(OUT / 'snr_sessions.csv')))
cols = ['date', 'subject', 'variant', 'recon', 'spec', 'gas_EI', 'gas_EE', 'gas_max', 'gas_min', 'dp_EI', 'dp_EE', 'dp',
        'lung_vox_EI', 'gas_EI_sept', 'dp_EI_sept', 'note', 'error']
num = {'gas_EI', 'gas_EE', 'gas_max', 'gas_min', 'dp_EI', 'dp_EE', 'lung_vox_EI', 'gas_EI_sept', 'dp_EI_sept'}
gvals = sorted(float(r['gas_EI']) for r in rows if r['gas_EI'])
lo, hi = gvals[len(gvals) // 10], gvals[9 * len(gvals) // 10]


def med(vals):
    v = [float(x) for x in vals if x]
    return f'{st.median(v):.1f} (n={len(v)})' if v else '–'


groups = defaultdict(list)
for r in rows:
    groups[('year', r['date'][:4])].append(r)
    groups[('recon', r['recon'])].append(r)
summ = ''.join(f'<tr><td>{k}</td><td>{html.escape(v)}</td><td>{len(g)}</td><td>{med(x["gas_EI"] for x in g)}</td>'
               f'<td>{med(x["dp_EI"] for x in g)}</td></tr>' for (k, v), g in sorted(groups.items()))


def cell(c, v):
    if c == 'gas_EI' and v:
        t = min(max((float(v) - lo) / (hi - lo + 1e-9), 0), 1)
        return f'<td class="n" style="background:hsla({int(120 * t)},65%,45%,.28)">{v}</td>'
    return f'<td class="{"n" if c in num else ""}">{html.escape(v)}</td>'


body = ''.join('<tr class="{}">{}</tr>'.format('flag' if r['note'] or r['error'] else '',
                                               ''.join(cell(c, r[c]) for c in cols)) for r in rows)
head = ''.join(f'<th data-n="{int(c in num)}">{c}</th>' for c in cols)
doc = f'''<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Xe Dynamic SNR</title><style>
:root{{--bg:#fbfaf7;--fg:#1d1d1b;--mut:#6b6a66;--line:#e2dfd8;--hd:#f0ede6;--flag:#fff3d6}}
@media (prefers-color-scheme:dark){{:root{{--bg:#1b1b1a;--fg:#e9e6df;--mut:#9a978f;--line:#33322f;--hd:#262523;--flag:#3a3222}}}}
body{{background:var(--bg);color:var(--fg);font:13px/1.4 -apple-system,system-ui,sans-serif;margin:16px}}
h1{{font-size:18px;margin:0 0 4px}} p{{color:var(--mut);margin:0 0 12px;max-width:900px}}
table{{border-collapse:collapse}} th,td{{padding:3px 7px;border-bottom:1px solid var(--line);white-space:nowrap}}
th{{background:var(--hd);position:sticky;top:0;cursor:pointer;text-align:left}} td.n{{text-align:right;font-variant-numeric:tabular-nums}}
tr.flag td{{background:var(--flag)}} .wrap{{overflow-x:auto}} .sum{{margin-bottom:16px}}
</style></head><body>
<h1>Xe-129 dynamic SNR — all reconstructed sessions ({len(rows)} rows)</h1>
<p>DIAPHRAGM bins. gas = mean(lung)/σ (corner noise); dp = total dissolved, whitened |z| Rician-corrected (snr_calc.py).
EI/EE = bin with largest/smallest lung mask. Primary recon = production settings (b44d, or the 2022 tree: b44 + deapod);
*_sept = the pre-b44 AIkill_Dynamic build (same as SNR_Table_All I–S). Deapod raises edge noise, so production and Sept values are
not interchangeable. Highlighted rows carry a note. Click a header to sort.</p>
<table class="sum"><tr><th>by</th><th>group</th><th>rows</th><th>median gas EI</th><th>median dp EI</th></tr>{summ}</table>
<div class="wrap"><table id="t"><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>
<script>
document.querySelectorAll('#t th').forEach((th,i)=>th.onclick=()=>{{const tb=document.querySelector('#t tbody'),n=th.dataset.n==='1',
d=th.dataset.d=th.dataset.d==='a'?'d':'a',rs=[...tb.rows];rs.sort((a,b)=>{{let x=a.cells[i].textContent,y=b.cells[i].textContent;
if(n){{x=parseFloat(x);y=parseFloat(y);x=isNaN(x)?-1e9:x;y=isNaN(y)?-1e9:y;return d==='a'?x-y:y-x}}return d==='a'?x.localeCompare(y):y.localeCompare(x)}});
rs.forEach(r=>tb.appendChild(r))}});
</script></body></html>'''
(OUT / 'snr_sessions.html').write_text(doc)
print(OUT / 'snr_sessions.html')
