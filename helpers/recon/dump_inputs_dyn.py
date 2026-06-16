"""Dump the DYNAMIC (multi-coil, all-interleave) inputs for the 4D CS pipeline.

Unlike dump_inputs.py (static, channel-0 only, single bin) this dumps every
coil and every interleave of the 3-min free-breathing acquisition, plus the
per-interleave respiratory surrogates Steve's loader already computes (SIGNAL,
and PNEUMOTACH if a pneumotach file is given). The diaphragm surrogate is built
separately (surrogates.py) because it needs a recon.

Reads the Siemens .dat directly via mapvbvd (the `_read_twix` / `_parse_pneumotach`
bodies are copied verbatim from convert_siemens_to_mrd.py so we do NOT need the
`mrd` package). The .dat -> arrays -> raw.load_from_arr path is exactly what
tyger_recon.reconstruct_from_mrd does, minus the MRD round-trip.

Outputs (into <outdir>):
    trajx/y/z.npy   gas-phase trajectory, grid-index units (rescaled to MS)
    acq_dyn.npy     (nch, npts*ntotalilvs) complex64, F-order per coil, mask applied
    ilvtime.npy     (ntotalilvs,) seconds, time of each gas interleave
    ilvvol_signal.npy   (ntotalilvs,) normalized signal surrogate (raw.py:429)
    ilvvol_pneumo.npy   (ntotalilvs,) normalized pneumotach surrogate (if pneumo given)
    exclude_ilv.npy (ntotalilvs,) bool, interleaves dropped by low-SNR exclusion
    meta.json       npts, nuniquesmp, ntotalilvs, nch, MS, IS, TR, bindt, nbins,
                    ilvperTR, killpts, datfile

Usage:
    .venv/bin/python dump_inputs_dyn.py --datdir DIR --gp-traj gp.npy --dp-traj dp.npy \
        [--pneumotach FILE] [--out OUTDIR] [--seqname STR] [--ms 240] [--is 100]
"""

import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")  # raw.py / load path may plot; never block headless
import numpy as np
import mapvbvd

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
sys.path.insert(0, REPO_ROOT)  # Steve's modules live at repo root

from gtypes import gvar, imgtype, graddir, bintype  # noqa: E402
import raw as steve_raw                              # noqa: E402


# --- copied from convert_siemens_to_mrd.py to avoid importing it (it `import mrd`) ---

def _read_twix(dat_file):
    """Read a Siemens .dat -> (raw (npts_full, nch, nilv) complex64, meta dict)."""
    twix = mapvbvd.mapVBVD(dat_file)
    if isinstance(twix, list):
        twix = twix[-1]
    twix.image.flagRemoveOS = False
    twix.image.squeeze = True
    raw = twix.image.unsorted().astype("complex64")
    if raw.ndim == 2:
        raw = raw[:, np.newaxis, :]  # add singleton channel dim -> (npts, 1, nilv)

    meta = dict(TR_us=0.0, TE_us=0.0, DPoff=0.0, dtdyn_ns=0.0,
                numspec_raw=0, dtspec_us=0.0)
    def _try(fn):
        try:
            return fn()
        except Exception:
            return None
    meta["TR_us"] = _try(lambda: float(twix.hdr.MeasYaps[("alTR", "0")])) or meta["TR_us"]
    meta["TE_us"] = _try(lambda: float(twix.hdr.MeasYaps[("alTE", "0")])) or meta["TE_us"]
    meta["DPoff"] = _try(lambda: float(twix.hdr.MeasYaps[("sWipMemBlock", "adFree", "2")])) or meta["DPoff"]
    meta["dtdyn_ns"] = _try(lambda: float(twix.hdr.MeasYaps[("sRXSPEC", "alDwellTime", "0")])) or meta["dtdyn_ns"]
    meta["dtdyn_ns"] = _try(lambda: float(twix.hdr.MeasYaps[("sRXSPEC", "alDwellTime", "1")])) or meta["dtdyn_ns"]
    meta["numspec_raw"] = _try(lambda: int(
        twix.hdr.MeasYaps[("sWipMemBlock", "alFree", "10")] *
        twix.hdr.MeasYaps[("sWipMemBlock", "alFree", "11")] + 0.1)) or meta["numspec_raw"]
    meta["dtspec_us"] = _try(lambda: float(twix.hdr.MeasYaps[("sWipMemBlock", "alFree", "12")])) or meta["dtspec_us"]
    return raw, meta


def _parse_pneumotach(pneumo_file):
    """Vendor pneumotach binary -> float64 (2, N): row0 = t (s, zeroed), row1 = pressure."""
    updatesendsize = 54
    with open(pneumo_file, mode="rb") as f:
        a = f.read()
    idx = [j for j in range(len(a) - updatesendsize) if a[j] == 0xA6 and a[j + 1] == 0x20]
    print(f"found {len(idx)} pressure measurements in {pneumo_file}")
    t = np.zeros(len(idx), dtype="float64")
    P = np.zeros(len(idx), dtype="float64")
    for cnt, j in enumerate(idx):
        t[cnt] = a[j + 5] * 2**24 + a[j + 4] * 2**16 + a[j + 3] * 2**8 + a[j + 2]
        P[cnt] = -20 + 90 * (a[j + 33] * 2**8 + a[j + 32]) / 65535.0
    t -= t[0]
    t /= 1000
    return np.stack((t, P))


def _classify_dat_files(datdir, seqname=None):
    """Largest .dat = dynamic free-breathing; second-largest = breath-hold ref."""
    dats = sorted([os.path.join(datdir, f) for f in os.listdir(datdir)
                   if f.endswith(".dat")], key=os.path.getsize)
    if seqname:
        dats = [f for f in dats if seqname in os.path.basename(f)]
    if not dats:
        sys.exit(f"no .dat files in {datdir}")
    dyn = dats[-1]
    ref = dats[-2] if len(dats) >= 2 else None
    print(f"dynamic:   {os.path.basename(dyn)} ({os.path.getsize(dyn)} bytes)")
    if ref:
        print(f"reference: {os.path.basename(ref)} ({os.path.getsize(ref)} bytes)")
    return dyn, ref


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datdir", required=True, help="folder with the .dat file(s)")
    ap.add_argument("--gp-traj", required=True, help="gas-phase trajectory .npy")
    ap.add_argument("--dp-traj", required=True, help="dissolved-phase trajectory .npy")
    ap.add_argument("--pneumotach", default=None, help="pneumotach binary (optional)")
    ap.add_argument("--seqname", default=None, help="substring filter on .dat names")
    ap.add_argument("--out", default="recon_io_dyn", help="output directory")
    ap.add_argument("--ms", type=int, default=None, help="override MS (default gvar)")
    ap.add_argument("--is", type=int, dest="IS", default=None, help="override IS (default gvar)")
    ap.add_argument("--no-ref", action="store_true", help="skip the breath-hold reference .dat")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    g = gvar()
    if args.ms is not None:
        g.MS = args.ms
    if args.IS is not None:
        g.IS = args.IS
    killpts = 2

    dynfile, reffile = _classify_dat_files(args.datdir, args.seqname)
    if args.no_ref:
        reffile = None

    # --- read raw .dat(s) -> channels-first (nch, npts, nilv), trim killpts on samples ---
    dynraw, m = _read_twix(dynfile)                       # (npts, nch, nilv)
    dyn_arr = np.ascontiguousarray(dynraw.transpose(1, 0, 2))[:, killpts:, :]
    ref_arr = None
    if reffile is not None:
        refraw, _ = _read_twix(reffile)
        ref_arr = np.ascontiguousarray(refraw.transpose(1, 0, 2))[:, killpts:, :]
    print(f"dyn raw channels-first {dyn_arr.shape}  TR_us={m['TR_us']}  numspec_raw={m['numspec_raw']}")

    meta = {"TR": m["TR_us"] * 1e-6, "TE": m["TE_us"] * 1e-6, "DPoff": m["DPoff"],
            "dtdyn": m["dtdyn_ns"] * 1e-9, "dtspec": m["dtspec_us"] * 1e-6,
            "numspec": int(m["numspec_raw"])}

    pneumo_arr = _parse_pneumotach(args.pneumotach) if args.pneumotach else None

    # --- trajectory (gas + dissolved) ---
    trajec = steve_raw.traj()
    trajec.killpts = killpts
    trajec.load_traj_from_npfile(args.gp_traj, args.dp_traj, nusimg=32)
    if trajec.npts == 0:
        sys.exit("trajectory load failed")
    print(f"traj: npts/ilv={trajec.npts}  nuniquesmp={trajec.nuniquesmp}")

    # --- run Steve's CPU loader (noise norm, gas/diss split, spike filter,
    #     exclusion ranges, SIGNAL + PNEUMOTACH surrogates) ---
    r = steve_raw.raw()
    r.load_from_arr(trajec, ref_arr, dyn_arr, pneumo_arr, "mrd_siemens", meta)
    print(f"raw: npts={r.npts} nch={r.nch} ntotalilvs={r.ntotalilvs} TR={r.TR} ilvperTR={r.ilvperTR}")

    trajec.rescale_to_MS(g.MS, g.IS)

    # per-sample exclusion mask (verbatim from results.dyn_recon / dump_inputs.py)
    npts, nilv = r.npts, r.ntotalilvs
    mask = np.ones(npts * nilv)
    exclude_ilv = np.zeros(nilv, dtype=bool)
    for iilv in range(nilv):
        for el in r.excluderanges:
            if iilv * r.TR >= el.start and iilv * r.TR <= el.stop:
                mask[(iilv * npts):((iilv + 1) * npts)] = 0
                exclude_ilv[iilv] = True

    # --- dump ---
    np.save(os.path.join(args.out, "trajx"), trajec.gettraj(imgtype.GPDYN, graddir.X))
    np.save(os.path.join(args.out, "trajy"), trajec.gettraj(imgtype.GPDYN, graddir.Y))
    np.save(os.path.join(args.out, "trajz"), trajec.gettraj(imgtype.GPDYN, graddir.Z))

    gpdyn = r.getimg(imgtype.GPDYN)  # (npts, nch, ntotalilvs)
    acq = np.empty((r.nch, npts * nilv), dtype=np.complex64)
    for ich in range(r.nch):
        acq[ich] = np.reshape(gpdyn[:, ich, :], npts * nilv, order="F") * mask
    np.save(os.path.join(args.out, "acq_dyn"), acq)

    np.save(os.path.join(args.out, "ilvtime"), np.asarray(r.ilvtime, dtype=float))
    np.save(os.path.join(args.out, "exclude_ilv"), exclude_ilv)

    sig = r.ilvvol[bintype.SIGNAL]
    if len(sig):
        np.save(os.path.join(args.out, "ilvvol_signal"), np.asarray(sig, dtype=float))
    pn = r.ilvvol[bintype.PNEUMOTACH]
    if len(pn):
        np.save(os.path.join(args.out, "ilvvol_pneumo"), np.asarray(pn, dtype=float))

    with open(os.path.join(args.out, "meta.json"), "w") as f:
        json.dump({"npts": int(npts), "nuniquesmp": int(trajec.nuniquesmp),
                   "ntotalilvs": int(nilv), "nch": int(r.nch),
                   "MS": int(g.MS), "IS": int(g.IS), "TR": float(r.TR),
                   "bindt": float(g.bindt), "nbins": int(g.nbins),
                   "ilvperTR": int(r.ilvperTR), "killpts": int(killpts),
                   "datfile": os.path.basename(dynfile)}, f, indent=2)

    print(f"wrote dump -> {args.out}")
    print(f"  acq_dyn {acq.shape} complex64  | signal={len(sig)>0} pneumo={len(pn)>0}"
          f"  | excluded {int(exclude_ilv.sum())}/{nilv} interleaves")


if __name__ == "__main__":
    main()
