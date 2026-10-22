"""Scratch: gas-only resubmit of a GP-DP input.mrd that has NO dissolved trajectory.

raw.py detects gas-dissolved (ilvperTR 2) from k0, but with no dissolved_phase_trajectory item the DPDYN gridding runs
on an empty trajectory and the GPU faults (CudaAPIError 700 at cuMemFree; 2022-02-02_2 MID14 single and MID16-24
merged, 2026-10-21; Feb-2022 seq fancy_20220202 has no dissolved calibration anywhere). Here the dynamic array keeps
only the gas phase (the parity with the larger median k0), so the detector sees gas-only and the in-set counter is the
gas arm counter; TR (alTR = per-scan spacing on this set-up, MDH 15.2 ms) is doubled to the gas spacing.
Then dyn_recon.submit + publish for each binning, same as a normal run.

--mode refasdyn instead promotes the reference_acquisition to the dynamic (and drops the reference): with two
equal-size .dat files the root converter's size sort can pick a different dynamic than dyn_recon logs (2022-02-03
000ID: MID00142 reconstructed, MID00143 sat in the reference slot), so the other scan is recoverable from the MRD.

usage: python3 helpers/_mrd_gas_only_resubmit.py <src input.mrd> <target study dir> [--mode gasonly|refasdyn]
       [--spec pipeline/recon_codespec_6b01071.yml]
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import mrd

PIPE = Path(__file__).resolve().parents[1] / 'pipeline'
sys.path.insert(0, str(PIPE))
import dyn_recon as D  # noqa: E402


def ref_as_dyn(src, dst, binning):
    with mrd.BinaryMrdReader(str(src)) as r:
        h = r.read_header()
        items = list(r.read_data())
    ref = next(it for it in items if 'reference_acquisition' in it.value.meta)
    out = []
    for it in items:
        if it is ref:
            continue
        if 'dynamic_acquisition' in it.value.meta:
            it.value.data = ref.value.data
            print(f'dynamic <- reference_acquisition {ref.value.data.shape}')
        out.append(it)
    for p in h.user_parameters.user_parameter_string:
        if p.name == 'binning':
            p.value = binning
    with mrd.BinaryMrdWriter(str(dst)) as w:
        w.write_header(h)
        w.write_data(iter(out))


def gas_only(src, dst, binning):
    with mrd.BinaryMrdReader(str(src)) as r:
        h = r.read_header()
        items = list(r.read_data())
    for it in items:
        a = it.value
        if 'dynamic_acquisition' in a.meta:
            d = a.data                                   # (nch, npts, nlines)
            k0 = np.abs(d[:, 1:4, :]).mean((0, 1))
            par = int(np.argmax([np.median(k0[0::2]), np.median(k0[1::2])]))
            ratio = np.median(k0[1 - par::2]) / np.median(k0[par::2])
            a.data = np.ascontiguousarray(d[:, :, par::2])
            print(f'gas parity {par}: kept {a.data.shape[2]} of {d.shape[2]} lines (other/gas k0 {ratio:.3f})')
    for p in h.user_parameters.user_parameter_double:
        if p.name == 'TR':
            p.value *= 2
            print(f'TR -> {p.value * 1e3:.1f} ms')
    for p in h.user_parameters.user_parameter_string:
        if p.name == 'binning':
            p.value = binning
    with mrd.BinaryMrdWriter(str(dst)) as w:
        w.write_header(h)
        w.write_data(iter(items))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('src')
    ap.add_argument('target')
    ap.add_argument('--mode', default='gasonly', choices=['gasonly', 'refasdyn'])
    ap.add_argument('--spec', default=str(PIPE / 'recon_codespec_6b01071.yml'))
    a = ap.parse_args()
    target = Path(a.target)
    build = gas_only if a.mode == 'gasonly' else ref_as_dyn
    for letter, binning in (('s', 'SIGNAL'), ('d', 'DIAPHRAGM')):
        run_dir = D.RUNS_DIR / f'{target.name}_{a.mode}_{letter}_{time.strftime("%Y%m%d-%H%M%S")}'
        run_dir.mkdir(parents=True, exist_ok=True)
        inp = run_dir / 'input.mrd'
        build(a.src, inp, binning)
        out = D.submit(inp, run_dir, Path(a.spec))
        D.publish(out, run_dir, target / letter, argparse.Namespace(slice=None))


if __name__ == '__main__':
    main()
