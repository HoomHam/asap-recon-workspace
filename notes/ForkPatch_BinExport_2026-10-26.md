# Fork patch (PREPARED, not applied): export bin timing + navigator with every binning mode — 2026-10-26

Status: DRAFT, root untouched. Apply on branch `diaphragm-recon` of the `hooman` fork after Hooman's go, build the
image (CI on push), add `pipeline/recon_codespec_<sha>.yml`, then `SPEC=… GO=1 zsh pipeline/run_sig_b44d.sh …`.
Why: the PCA Step-1 loader wants, inside the SIGNAL-binned recon output, (i) a kernel-weighted sub-voxel diaphragm
position per bin and (ii) the per-bin mean real time, without an external period table. Today the navigator only
runs (and `nav_*` are only exported) when `binning == DIAPHRAGM` (`tyger_recon.py:270-274`), and the bin fraction per
interleave (`g_raw.ilvbin[bt]`) and the EE times (`g_raw.volmeasEEtime[bt]`) are never exported.
Steve's-knob check: no existing knob does this (main.py GUI has no export); the arrays already exist in `raw`.

## 1. `tyger_recon.py` — after the binning choice (replace lines 266-283)

```python
    binning = user_string.get('binning', 'SIGNAL')
    # export_nav=1 (MRD user long): run the diaphragm navigator even when binning by SIGNAL / PNEUMOTACH so the
    # output carries nav_* (sub-voxel diaphragm z per 26-ilv frame, interpolated per interleave) alongside the
    # bins actually used. Costs one low-res pass (~1 min GPU). Default off = unchanged behaviour.
    export_nav = int(user_long.get('export_nav', 0))
    nav = None
    if (binning == 'DIAPHRAGM' or export_nav) and not len(g_raw.ilvbin[bintype.DIAPHRAGM]):
        nav = _diaphragm_navigator(g, g_raw, g_traj, g_res)
        print(f'recomputing full-res b (MS={g.MS})...', file=stderr)
        g_res.calcb(g, g_raw, g_traj, g.usegpu)

    if binning == 'DIAPHRAGM' and len(g_raw.ilvbin[bintype.DIAPHRAGM]):
        bt = bintype.DIAPHRAGM
    elif binning == 'PNEUMOTACH' and len(g_raw.ilvbin[bintype.PNEUMOTACH]):
        bt = bintype.PNEUMOTACH
    else:
        bt = bintype.SIGNAL
    print(f'binned dynamic recon with {bt.name} binning...', file=stderr)
    g_res.dyn_recon(g, g_raw, g_traj, g_raw.ilvbin[bt], g.usegpu)

    # per-interleave bin assignment + per-bin timing, with the SAME kernel the gridding used (recon.py:104-111),
    # so a reader can reproduce any per-bin average without an external table.
    bins = np.asarray(g_raw.ilvbin[bt], dtype='float64')          # fraction in [0,1), -1 = not binned
    binexp = _bin_export(g, g_raw, bins, bt, nav)
    if nav is None:
        nav = {}
    nav.update(binexp)
    ...
    _write_results_to_mrd(g_res, header, output, nav=nav if len(nav) else None)
```

## 2. `tyger_recon.py` — new helper (put next to `_diaphragm_navigator`)

```python
def _bin_export(g, g_raw, bins, bt, nav):
    """Per-interleave bins and per-bin kernel-weighted timing / navigator statistics.

    bins: g_raw.ilvbin[bt] (fraction of the cycle in [0,1), -1 = excluded). Kernel = recon.py cudarecon:
    w = exp(-d^2 / bindist0sq) with d the cyclic distance between the integer bin and bins*nbins, bindist0sq = 2.
    Time after EE: SIGNAL / PNEUMOTACH use volmeasEEtime[bt]; DIAPHRAGM has no EE list, so the EE of the
    navigator itself is used (strict minima of ilvvol within +-50 interleaves, raw.py's rule).
    """
    from recon import bindist0sq
    nb = g.nbins
    t = np.asarray(g_raw.ilvtime, dtype='float64')
    ok = bins >= 0
    c = bins * nb
    w = np.zeros((len(bins), nb))
    for ib in range(nb):
        d = np.minimum.reduce([np.abs(ib - c), np.abs(ib + nb - c), np.abs(ib - nb - c)])
        w[:, ib] = np.where(ok, np.exp(-d * d / bindist0sq), 0.0)
    ee = np.asarray(g_raw.volmeasEEtime[bt], dtype='float64')
    if bt == bintype.DIAPHRAGM or ee.size < 2:
        v = np.asarray(g_raw.ilvvol[bintype.DIAPHRAGM], dtype='float64')
        vv = np.where(np.isnan(v), np.inf, v)
        N = 50
        ee = np.array([t[j] for j in range(N, len(t) - N - 1)
                       if np.isfinite(vv[j]) and np.all(vv[j - N:j] > vv[j]) and np.all(vv[j + 1:j + N + 1] > vv[j])])
    tee = np.full(len(t), np.nan)
    if ee.size >= 2:
        k = np.searchsorted(ee, t, side='right') - 1
        inside = (k >= 0) & (k < ee.size - 1)
        tee[inside] = t[inside] - ee[k[inside]]
    def wstat(x):
        m = np.isfinite(x)
        ww = w * m[:, None]
        s = ww.sum(0) + 1e-12
        mean = (ww * np.nan_to_num(x)[:, None]).sum(0) / s
        sd = np.sqrt(np.maximum((ww * (np.nan_to_num(x)[:, None] - mean) ** 2).sum(0) / s, 0.0))
        return mean.astype('float32'), sd.astype('float32')
    out = {'ilvbin': bins.astype('float32'), 'ilvtime': t.astype('float32'),
           'ee_times': ee.astype('float32'), 'bin_weff': w.sum(0).astype('float32'),
           'bin_n_hard': np.array([np.sum(ok & (np.floor(c + 0.5).astype(int) % nb == ib)) for ib in range(nb)], 'float32')}
    out['bin_time_mean_s'], out['bin_time_sd_s'] = wstat(tee)
    if nav is not None and len(g_raw.ilvvol[bintype.DIAPHRAGM]):
        out['bin_nav_volume'], out['bin_nav_volume_sd'] = wstat(np.asarray(g_raw.ilvvol[bintype.DIAPHRAGM], 'float64'))
        fz, ft = nav['nav_diaphragm_z'].astype('float64'), nav['nav_time'].astype('float64')
        good = np.isfinite(fz)
        navz = np.interp(t, ft[good], fz[good], left=np.nan, right=np.nan)
        out['bin_nav_z_vox'], out['bin_nav_z_vox_sd'] = wstat(navz)
    return out
```

`_write_results_to_mrd` already streams every key of `nav` as an NdArrayFloat tagged with its own name, so no change
there. `ilvbin` is the production bin array (what `dyn_recon` gridded); `bin_time_mean_s` is the kernel-weighted
mean of (t − t_EE) over the interleaves of each bin (SIGNAL: exactly b/16 × period on average).

## 3. `workspace/pipeline/mrd_to_mat.py`

```python
NAV_KEYS = ('nav_coronal', 'nav_diaphragm_z', 'nav_time', 'nav_volume', 'nav_ilvtime', 'nav_volmeastime',
            'ilvbin', 'ilvtime', 'ee_times', 'bin_weff', 'bin_n_hard', 'bin_time_mean_s', 'bin_time_sd_s',
            'bin_nav_volume', 'bin_nav_volume_sd', 'bin_nav_z_vox', 'bin_nav_z_vox_sd')
```
and `d = Path(sess_dir) / 'd'` → take the sub-folder letter from the path (the SIGNAL tree writes `…/s/`).

## 4. Converter

`convert_siemens_to_mrd.py` must write `export_nav` (user long) when asked; `dyn_recon.py` gets `--export-nav`
passing it through (same plumbing as `--sigmin-n`, commit 262b16f).

## 5. Offline equivalent available today (no image build)

`workspace/outputs/bin_timing_2026-10-26/<study>_{ilv,bins,breaths}.csv` from `helpers/_bin_timing_export.py`
give the same per-interleave / per-bin numbers (SIGNAL bins re-derived with Steve's EE rule from the stored k0
signal; nav z from the b44d d/ recon of the same session). Numerically identical to what the patch would export,
except that the patch reads `ilvbin` straight from the recon instead of re-deriving it.
