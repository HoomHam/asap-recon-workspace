#!/usr/bin/env python
"""Rigid overlay of registered xenon-129 gas MRI (EI bin) on RT 4D-CT (phase 0% = end-inhale).

Usage:
    ct_mri_overlay.py search  [subj ...]   # orientation search: 6 axis perms x 8 flips, rigid MI, low-res
    ct_mri_overlay.py refine  [subj ...]   # full-res rigid refine of chosen orientation (+ similarity diag)
    ct_mri_overlay.py overlay [subj ...]   # axial (all lung slices) + coronal (1 mm interp) overlays
    ct_mri_overlay.py all     [subj ...]
    ct_mri_overlay.py video   [subj ...]   # 16-bin video through the rigid transform (mp4 repo, gif Ext)
    ct_mri_overlay.py phases  [subj ...]   # diagnostic: MI/containment of the EI bin vs every CT phase

MRI stacks are Tyger numpy order (bins, Z, Y, X), 100^3, FOV 350 mm -> 3.5 mm isotropic (raw.py FOV=350).
NIfTI legacy files are (X,Y,Z,T) in nibabel = same memory order after moveaxis(-1,0) (PCA fact F53).
No deformation anywhere: rigid (6 DOF) is the delivered transform; similarity (7 DOF) only as a
scale diagnostic for the FOV assumption.

Outputs
    Ext  : /Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/ct_overlay/<subj>/
           mri_ei_native.nii.gz, mri_on_ct.nii.gz, rigid.tfm, ct_lungmask.nii.gz, axial_png/, coronal_png/
    Repo : workspace/outputs/ct_overlay/<subj>/  orientation_scores.csv, register.json, montages, gifs
"""
import csv
import glob
import itertools
import json
import os
import sys

import numpy as np
import SimpleITK as sitk
from scipy import ndimage as ndi

EXT = '/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/ct_overlay'
REPO_OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'outputs', 'ct_overlay')
PCA = '/Volumes/HoomHamExt/Work/Codes/2026_PCA_registration'
MRI_VOX = 3.5  # mm, 350 mm FOV / 100

SUBJECTS = {
    # subj: (registered stack, EI bin (0-based), CT phase tag, provenance)
    '01BB': (f'{PCA}/elastix_oneshot_2026-09-16/2023-03-27_02BB/registered_elastix.npy', 8, 'phase00',
             'elastix one-shot PCAMetric2 (no legacy final); Tyger 2023-03-27_02BB; CT same day'),
    '02ZS': (f'{PCA}/elastix_oneshot_2026-09-16/2023-11-06_002ZS/registered_elastix.npy', 8, 'phase00',
             'elastix one-shot PCAMetric2 (no legacy final); Tyger 2023-11-06_002ZS; CT 2023-10-23'),
    '03PM': ('/Volumes/HoomHamExt/Work/Analysis/2024-12-06_003PM/reg/5000/registered.nii', 6, 'phase00',
             'legacy v11 final (gas only); Tyger 2023-12-06_003PM; CT 2023-11-20'),
}


MANIFEST = os.path.join(EXT, 'manifest.json')   # written by ct_cohort_batch.py: key -> {mri, ei, phase, prov, ct_subject}


def _spec(subj):
    """(mri_path, ei_bin, ct_phase_tag, provenance) for a key; EI 'auto' = bin with max ventilated volume."""
    if subj in SUBJECTS:
        return SUBJECTS[subj]
    with open(MANIFEST) as fh:
        m = json.load(fh)[subj]
    ei = m['ei']
    if ei == 'auto':
        ei = auto_ei(m['mri'])
    return m['mri'], ei, m['phase'], m['prov']


def _ct_subject(subj):
    if subj in SUBJECTS:
        return subj
    with open(MANIFEST) as fh:
        return json.load(fh)[subj].get('ct_subject', subj)


def load_stack(path):
    if path.endswith('.npy'):
        return np.load(path)
    import nibabel as nib
    return np.moveaxis(np.asarray(nib.load(path).dataobj), -1, 0)


def auto_ei(path):
    stack = np.abs(load_stack(path)).astype(np.float32)
    thr = 0.25 * np.percentile(stack, 99.5)
    vol = [(b > thr).sum() for b in stack]
    return int(np.argmax(vol))


AXES = 'SPL'  # target sitk array order (k=S, j=P, i=L) -> sitk index (i=L, j=P, k=S) == CT LPS


# ---------------------------------------------------------------- loading
def load_mri_ei(path, ei):
    stack = load_stack(path)
    v = np.abs(stack[ei]).astype(np.float32)
    v /= np.percentile(v, 99.5)
    return np.clip(v, 0, 1.5)


def orient_mri(v, perm, flips):
    """v: (Z,Y,X) Tyger order. perm: tuple of v-axes -> (S,P,L). flips: signs per (S,P,L)."""
    a = np.transpose(v, perm)
    for ax, f in enumerate(flips):
        if f < 0:
            a = np.flip(a, axis=ax)
    return np.ascontiguousarray(a)


def mri_to_sitk(a, center_at):
    img = sitk.GetImageFromArray(a)
    img.SetSpacing((MRI_VOX,) * 3)
    m = a > 0.25
    m = ndi.binary_opening(m, iterations=1)
    lab, n = ndi.label(m)
    if n:
        sizes = ndi.sum(m, lab, range(1, n + 1))
        keep = np.isin(lab, 1 + np.where(sizes > 0.05 * sizes.max())[0])
        m = keep
    cz, cy, cx = ndi.center_of_mass(m)             # array (k,j,i)
    cen_phys = np.array([cx, cy, cz]) * MRI_VOX     # (x,y,z) with origin 0
    img.SetOrigin(tuple(np.asarray(center_at) - cen_phys))
    return img, m


def ct_lung_mask(ct):
    """Binary lung+airway mask from HU: air < -400 inside body (components not touching x/y border)."""
    a = sitk.GetArrayFromImage(ct)                # (k,j,i)
    air = a < -400
    # body = per-slice largest soft-tissue component, holes filled (encloses lungs + airways;
    # robust to trachea->pharynx->mouth connecting lung air to outside air)
    body = np.zeros_like(air)
    for k in range(a.shape[0]):
        tis = ndi.binary_opening(a[k] > -500, iterations=2)
        lab, n = ndi.label(tis)
        if n == 0:
            continue
        sizes = ndi.sum(tis, lab, range(1, n + 1))
        big = lab == (1 + int(np.argmax(sizes)))
        body[k] = ndi.binary_fill_holes(big)
    lung = air & body
    lab, n = ndi.label(lung)
    sizes = ndi.sum(lung, lab, range(1, n + 1))
    keep = 1 + np.where(sizes > 0.02 * sizes.max())[0]
    lung = np.isin(lab, keep)
    lung = ndi.binary_closing(lung, iterations=2)
    lung = ndi.binary_fill_holes(lung)
    return lung


def np_to_sitk_like(a, ref, dtype=np.float32):
    img = sitk.GetImageFromArray(a.astype(dtype))
    img.CopyInformation(ref)
    return img


def resample_iso(img, spacing, interp=sitk.sitkLinear, default=-1024.0):
    sp = np.array(img.GetSpacing()); sz = np.array(img.GetSize())
    new_sz = np.round(sz * sp / spacing).astype(int).tolist()
    return sitk.Resample(img, new_sz, sitk.Transform(), interp, img.GetOrigin(), (spacing,) * 3,
                         img.GetDirection(), default, img.GetPixelID())


# ---------------------------------------------------------------- registration
def register_rigid(fixed, moving, fixed_mask, center, iters=200, shrink=(4, 2, 1), sigmas=(4, 2, 0),
                   transform=None, sampling=0.2, seed=1):
    reg = sitk.ImageRegistrationMethod()
    reg.SetMetricAsMattesMutualInformation(numberOfHistogramBins=32)
    reg.SetMetricSamplingStrategy(reg.RANDOM)
    reg.SetMetricSamplingPercentage(sampling, seed)
    reg.SetMetricFixedMask(fixed_mask)
    reg.SetInterpolator(sitk.sitkLinear)
    reg.SetOptimizerAsRegularStepGradientDescent(learningRate=2.0, minStep=1e-3, numberOfIterations=iters,
                                                 relaxationFactor=0.6, gradientMagnitudeTolerance=1e-5)
    reg.SetOptimizerScalesFromPhysicalShift()
    reg.SetShrinkFactorsPerLevel(list(shrink))
    reg.SetSmoothingSigmasPerLevel(list(sigmas))
    reg.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()
    if transform is None:
        transform = sitk.Euler3DTransform()
        transform.SetCenter(tuple(float(c) for c in center))
    reg.SetInitialTransform(transform, inPlace=False)
    out = reg.Execute(fixed, moving)
    return out, reg.GetMetricValue(), reg.GetOptimizerStopConditionDescription()



def as_param(tf, cls):
    """Downcast a registration result (Transform / CompositeTransform) to `cls`."""
    t = tf
    if t.GetName() == 'CompositeTransform':
        t = sitk.CompositeTransform(t).GetNthTransform(0)
    t = t.Downcast() if hasattr(t, 'Downcast') else t
    return cls(t)

def containment(fixed_lung_img, moving, moving_mask_arr, tf):
    """Fraction of MRI ventilated voxels landing inside the CT lung mask."""
    mm = sitk.GetImageFromArray(moving_mask_arr.astype(np.uint8)); mm.CopyInformation(moving)
    res = sitk.Resample(mm, fixed_lung_img, tf, sitk.sitkNearestNeighbor, 0, sitk.sitkUInt8)
    r = sitk.GetArrayFromImage(res).astype(bool)
    l = sitk.GetArrayFromImage(fixed_lung_img).astype(bool)
    return float((r & l).sum() / max(r.sum(), 1)), float((r & l).sum() / max(l.sum(), 1))


def ct_windowed(ct):
    return sitk.Cast(sitk.Clamp(ct, sitk.sitkFloat32, -1000, 200), sitk.sitkFloat32)


def ct_paths(subj, phase):
    return os.path.join(EXT, _ct_subject(subj), 'ct', f'{phase}.nii.gz')


def load_ct_and_masks(subj, phase):
    ct = sitk.ReadImage(ct_paths(subj, phase), sitk.sitkFloat32)
    lung = ct_lung_mask(ct)
    lung_img = np_to_sitk_like(lung, ct, np.uint8)
    lung_img = sitk.Cast(lung_img, sitk.sitkUInt8)
    # registration ROI = lung dilated ~25 mm (includes chest wall / mediastinum edges)
    roi = ndi.binary_dilation(lung, iterations=1, structure=np.ones((3, 3, 3)))
    roi = ndi.binary_dilation(roi, iterations=8)
    roi_img = sitk.Cast(np_to_sitk_like(roi, ct, np.uint8), sitk.sitkUInt8)
    cen = ndi.center_of_mass(lung)                  # (k,j,i)
    center = ct.TransformContinuousIndexToPhysicalPoint((float(cen[2]), float(cen[1]), float(cen[0])))
    return ct, lung_img, roi_img, center


# ---------------------------------------------------------------- commands
def cmd_search(subjects):
    for subj in subjects:
        path, ei, phase, _ = SUBJECTS[subj]
        out_dir = os.path.join(REPO_OUT, subj); os.makedirs(out_dir, exist_ok=True)
        ct, lung_img, roi_img, center = load_ct_and_masks(subj, phase)
        ct_lo = resample_iso(ct_windowed(ct), MRI_VOX)
        lung_lo = resample_iso(lung_img, MRI_VOX, sitk.sitkNearestNeighbor, 0)
        roi_lo = resample_iso(roi_img, MRI_VOX, sitk.sitkNearestNeighbor, 0)
        v = load_mri_ei(path, ei)
        rows = []
        for perm in itertools.permutations(range(3)):
            for flips in itertools.product((1, -1), repeat=3):
                a = orient_mri(v, perm, flips)
                mov, mmask = mri_to_sitk(a, center)
                try:
                    tf, mi, stop = register_rigid(ct_lo, mov, roi_lo, center, iters=120, shrink=(2, 1), sigmas=(2, 0))
                    c_in, c_cov = containment(lung_lo, mov, mmask, tf)
                except Exception as e:  # noqa
                    mi, c_in, c_cov, stop = np.nan, np.nan, np.nan, f'ERR {e}'
                rows.append(dict(perm=''.join(map(str, perm)), flips=''.join('+' if f > 0 else '-' for f in flips),
                                 mi=mi, contain=c_in, cover=c_cov, stop=stop[:40]))
                print(f'{subj} perm={rows[-1]["perm"]} flips={rows[-1]["flips"]} MI={mi:.4f} in={c_in:.3f} cov={c_cov:.3f}')
        rows.sort(key=lambda r: r['mi'] if np.isfinite(r['mi']) else 0)
        with open(os.path.join(out_dir, 'orientation_scores.csv'), 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
        print(f'{subj} BEST by MI: {rows[0]}')


def choose_orientation(subj, override=None):
    if override:
        return override
    with open(os.path.join(REPO_OUT, subj, 'orientation_scores.csv')) as fh:
        rows = list(csv.DictReader(fh))
    rows.sort(key=lambda r: float(r['mi']))
    return rows[0]['perm'], rows[0]['flips']


def parse_orient(perm_s, flips_s):
    return tuple(int(c) for c in perm_s), tuple(1 if c == '+' else -1 for c in flips_s)


def cmd_refine(subjects, orient=None, phase_override=None):
    for subj in subjects:
        path, ei, phase, prov = _spec(subj)
        phase = phase_override or phase
        out_dir = os.path.join(REPO_OUT, subj); os.makedirs(out_dir, exist_ok=True)
        ext_dir = os.path.join(EXT, subj); os.makedirs(ext_dir, exist_ok=True)
        perm_s, flips_s = choose_orientation(subj, orient)
        perm, flips = parse_orient(perm_s, flips_s)
        ct, lung_img, roi_img, center = load_ct_and_masks(subj, phase)
        ctw = ct_windowed(ct)
        if min(ct.GetSpacing()) < 1.4:                              # thin-slice CT: register at 1.5 mm iso
            ctw = resample_iso(ctw, 1.5); roi_img = resample_iso(roi_img, 1.5, sitk.sitkNearestNeighbor, 0)
        v = load_mri_ei(path, ei)
        a = orient_mri(v, perm, flips)
        mov, mmask = mri_to_sitk(a, center)
        sitk.WriteImage(mov, os.path.join(ext_dir, 'mri_ei_native.nii.gz'), True)
        sitk.WriteImage(lung_img, os.path.join(ext_dir, 'ct_lungmask.nii.gz'), True)
        # stage 1: low-res rigid (as in search), stage 2: full-res rigid from it
        ct_lo = resample_iso(ctw, MRI_VOX); roi_lo = resample_iso(roi_img, MRI_VOX, sitk.sitkNearestNeighbor, 0)
        tf1, mi1, _ = register_rigid(ct_lo, mov, roi_lo, center, iters=200, shrink=(2, 1), sigmas=(2, 0))
        tf2, mi2, stop2 = register_rigid(ctw, mov, roi_img, center, iters=150, shrink=(4, 2), sigmas=(3, 1),
                                         transform=as_param(tf1, sitk.Euler3DTransform), sampling=0.05)
        c_in, c_cov = containment(lung_img, mov, mmask, tf2)
        sitk.WriteTransform(tf2, os.path.join(ext_dir, 'rigid.tfm'))
        # scale diagnostic: similarity transform seeded from rigid (NOT used for the overlay)
        sim = sitk.Similarity3DTransform()
        e = as_param(tf2, sitk.Euler3DTransform)
        sim.SetCenter(e.GetCenter())
        sim.SetMatrix(e.GetMatrix()); sim.SetTranslation(e.GetTranslation())
        tf3, mi3, _ = register_rigid(ct_lo, mov, roi_lo, center, iters=200, shrink=(2, 1), sigmas=(2, 0), transform=sim)
        scale = as_param(tf3, sitk.Similarity3DTransform).GetScale()
        # resample MRI onto CT grid (B-spline, no deformation: rigid only)
        mri_on_ct = sitk.Resample(mov, ct, tf2, sitk.sitkBSpline, 0.0, sitk.sitkFloat32)
        sitk.WriteImage(mri_on_ct, os.path.join(ext_dir, 'mri_on_ct.nii.gz'), True)
        e2 = e
        info = dict(subject=subj, mri=path, ei_bin=ei, ct_phase=phase, provenance=prov,
                    mri_voxel_mm=MRI_VOX, orientation=dict(perm=perm_s, flips=flips_s, meaning='MRI (Z,Y,X) axes -> (S,P,L) with signs'),
                    rigid=dict(euler_deg=[float(np.degrees(x)) for x in (e2.GetAngleX(), e2.GetAngleY(), e2.GetAngleZ())],
                               translation_mm=[float(x) for x in e2.GetTranslation()], center=list(e2.GetCenter()),
                               mi_lowres=float(mi1), mi_fullres=float(mi2), stop=stop2,
                               mri_inside_ctlung=c_in, ctlung_covered_by_mri=c_cov),
                    similarity_diag=dict(scale=float(scale), mi=float(mi3),
                                         note='7-DOF fit only to check the 350 mm FOV assumption; overlay uses rigid'))
        with open(os.path.join(out_dir, 'register.json'), 'w') as fh:
            json.dump(info, fh, indent=2)
        print(json.dumps(info['rigid'], indent=1)); print('similarity scale', scale)


# ---------------------------------------------------------------- figures
def _body_mask(C):
    """Per-slice largest soft-tissue component with holes filled (same rule as ct_lung_mask)."""
    body = np.zeros(C.shape, bool)
    for k in range(C.shape[0]):
        tis = ndi.binary_opening(C[k] > -500, iterations=2)
        lab, n = ndi.label(tis)
        if n:
            sizes = ndi.sum(tis, lab, range(1, n + 1))
            body[k] = ndi.binary_fill_holes(lab == (1 + int(np.argmax(sizes))))
    return body


def _contig_run(area, frac):
    """Indices of the contiguous run around argmax(area) where area >= frac*max (drops sinus/pharynx air)."""
    thr = frac * area.max(); i = int(np.argmax(area)); lo = hi = i
    while lo > 0 and area[lo - 1] >= thr:
        lo -= 1
    while hi < len(area) - 1 and area[hi + 1] >= thr:
        hi += 1
    return np.arange(lo, hi + 1)


def _overlay_rgb(ct_slice, mri_slice, lung_slice, vmin=-1100, vmax=400, thr=0.2, alpha=0.55):
    import matplotlib
    g = np.clip((ct_slice - vmin) / (vmax - vmin), 0, 1)
    rgb = np.stack([g, g, g], -1)
    cm = matplotlib.colormaps['jet']
    col = cm(np.clip(mri_slice, 0, 1))[..., :3]
    w = (mri_slice > thr).astype(float) * alpha
    rgb = rgb * (1 - w[..., None]) + col * w[..., None]
    # CT lung contour in cyan
    edge = lung_slice ^ ndi.binary_erosion(lung_slice)
    rgb[edge] = [0.0, 1.0, 1.0]
    return rgb


def _montage(imgs, ncols, labels=None, out=None, title=None):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    n = len(imgs); nrows = int(np.ceil(n / ncols))
    h, w = imgs[0].shape[:2]
    fig, axs = plt.subplots(nrows, ncols, figsize=(ncols * 2.2, nrows * 2.2 * h / w + 0.3))
    axs = np.atleast_1d(axs).ravel()
    for k, ax in enumerate(axs):
        ax.axis('off')
        if k < n:
            ax.imshow(imgs[k]); ax.set_title(labels[k] if labels else '', fontsize=7, pad=1)
    if title:
        fig.suptitle(title, fontsize=9, y=0.995)
    plt.tight_layout(pad=0.2, rect=(0, 0, 1, 0.985)); fig.savefig(out, dpi=110); plt.close(fig)


def cmd_overlay(subjects):
    import imageio.v2 as imageio
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    for subj in subjects:
        path, ei, phase, prov = _spec(subj)
        out_dir = os.path.join(REPO_OUT, subj); ext_dir = os.path.join(EXT, subj)
        with open(os.path.join(out_dir, 'register.json')) as fh:
            phase = json.load(fh).get('ct_phase', phase)          # phase actually registered
        ct = sitk.ReadImage(ct_paths(subj, phase), sitk.sitkFloat32)
        mri = sitk.ReadImage(os.path.join(ext_dir, 'mri_on_ct.nii.gz'))
        lung = sitk.ReadImage(os.path.join(ext_dir, 'ct_lungmask.nii.gz'))
        with open(os.path.join(out_dir, 'register.json')) as fh:
            info = json.load(fh)
        sp = ct.GetSpacing()
        C = sitk.GetArrayFromImage(ct); M = sitk.GetArrayFromImage(mri); L = sitk.GetArrayFromImage(lung).astype(bool)
        # clip MRI overlay to tissue ∪ (lung ⊕ 15 mm): drops noise specks in air/immobilization-bag gaps,
        # keeps gas up to 15 mm outside the CT lung contour (visible misfit) and soft-tissue-embedded signal
        keep = (C > -500) | ndi.binary_dilation(L, iterations=5, structure=np.ones((3, 3, 3)))
        M = M * keep
        # ---------- axial: every CT slice with lung, radiological (anterior up, patient L on image right)
        ks = _contig_run(L.sum(axis=(1, 2)), 0.10)
        k0, k1 = max(ks[0] - 2, 0), min(ks[-1] + 3, C.shape[0])
        Lsel = L.copy(); Lsel[:k0] = False; Lsel[k1:] = False          # lung proper, for coronal/sagittal ranges
        # crop in-plane to body bbox
        body = C > -500
        jj = np.where(body.any(axis=(0, 2)))[0]; ii = np.where(body.any(axis=(0, 1)))[0]
        j0, j1, i0, i1 = jj[0], jj[-1] + 1, ii[0], ii[-1] + 1
        ax_dir = os.path.join(ext_dir, 'axial_png'); os.makedirs(ax_dir, exist_ok=True)
        frames, labels = [], []
        kstep = max(1, int(round(3.0 / sp[2])))                 # ~3 mm between exported axial slices
        for k in range(k0, k1, kstep):
            rgb = _overlay_rgb(C[k, j0:j1, i0:i1], M[k, j0:j1, i0:i1], L[k, j0:j1, i0:i1])
            z_mm = ct.TransformIndexToPhysicalPoint((0, 0, int(k)))[2]
            fn = os.path.join(ax_dir, f'axial_{k:03d}.png')
            plt.imsave(fn, rgb)
            frames.append((rgb * 255).astype(np.uint8)); labels.append(f'k={k} z={z_mm:.0f}')
        imageio.mimsave(os.path.join(ext_dir, 'axial.gif'), frames, duration=0.15, loop=0)
        step = max(1, len(frames) // 40)
        _montage(frames[::step], 8, labels[::step], os.path.join(out_dir, 'axial_montage.png'),
                 title=f'{subj} axial: CT {phase} + Xe gas EI bin {ei} (rigid). cyan = CT lung contour')
        # ---------- coronal: resample both to 1 mm iso (CT 3 mm slices interpolated), slices every 4 mm
        iso = 1.0
        ct_i = resample_iso(ct, iso, sitk.sitkBSpline)
        mri_i = sitk.Resample(mri, ct_i, sitk.Transform(), sitk.sitkLinear, 0.0)
        lung_i = sitk.Resample(lung, ct_i, sitk.Transform(), sitk.sitkNearestNeighbor, 0)
        Ci = sitk.GetArrayFromImage(ct_i); Li = sitk.GetArrayFromImage(lung_i).astype(bool)
        Mi = sitk.GetArrayFromImage(sitk.Resample(np_to_sitk_like(M, ct), ct_i, sitk.Transform(), sitk.sitkLinear, 0.0))
        Lsel_i = sitk.GetArrayFromImage(sitk.Resample(np_to_sitk_like(Lsel, ct, np.uint8), ct_i, sitk.Transform(), sitk.sitkNearestNeighbor, 0)).astype(bool)
        js = _contig_run(Lsel_i.sum(axis=(0, 2)), 0.05)
        jsel = range(max(js[0] - 5, 0), min(js[-1] + 6, Ci.shape[1]), 4)
        kk = np.where(Lsel_i.any(axis=(1, 2)))[0]; kk0, kk1 = max(kk[0] - 15, 0), min(kk[-1] + 16, Ci.shape[0])
        ii = np.where((Ci > -500).any(axis=(0, 1)))[0]; ci0, ci1 = ii[0], ii[-1] + 1
        co_dir = os.path.join(ext_dir, 'coronal_png'); os.makedirs(co_dir, exist_ok=True)
        frames, labels = [], []
        for j in jsel:
            # rows: S at top -> flip k
            rgb = _overlay_rgb(Ci[kk0:kk1, j, ci0:ci1][::-1], Mi[kk0:kk1, j, ci0:ci1][::-1], Li[kk0:kk1, j, ci0:ci1][::-1])
            y_mm = ct_i.TransformIndexToPhysicalPoint((0, int(j), 0))[1]
            plt.imsave(os.path.join(co_dir, f'coronal_{j:03d}.png'), rgb)
            frames.append((rgb * 255).astype(np.uint8)); labels.append(f'y={y_mm:.0f} (P+)')
        imageio.mimsave(os.path.join(ext_dir, 'coronal.gif'), frames, duration=0.2, loop=0)
        _montage(frames, 6, labels, os.path.join(out_dir, 'coronal_montage.png'),
                 title=f'{subj} coronal (1 mm interp from 3 mm axial): CT {phase} + Xe gas EI bin {ei} (rigid)')
        # ---------- 3-plane summary panel
        kc = int(np.round(ndi.center_of_mass(Lsel)[0])); jc = int(np.round(ndi.center_of_mass(Lsel_i)[1]))
        ii_l = np.where(Lsel_i.any(axis=(0, 1)))[0]; ic = int(np.round(0.5 * (ii_l[0] + ndi.center_of_mass(Lsel_i)[2])))  # right lung
        pan = [_overlay_rgb(C[kc, j0:j1, i0:i1], M[kc, j0:j1, i0:i1], L[kc, j0:j1, i0:i1]),
               _overlay_rgb(Ci[kk0:kk1, jc, ci0:ci1][::-1], Mi[kk0:kk1, jc, ci0:ci1][::-1], Li[kk0:kk1, jc, ci0:ci1][::-1]),
               _overlay_rgb(Ci[kk0:kk1, :, ic][::-1], Mi[kk0:kk1, :, ic][::-1], Li[kk0:kk1, :, ic][::-1])]
        fig, axs = plt.subplots(1, 3, figsize=(15, 5.5))
        for a_, im, t in zip(axs, pan, ['axial (A up, L right)', 'coronal (S up, L right)', 'sagittal, right lung (S up, P right)']):
            a_.imshow(im); a_.set_title(t); a_.axis('off')
        r = info['rigid']
        fig.suptitle(f"{subj}: rot {np.round(r['euler_deg'],1)} deg, trans {np.round(r['translation_mm'],1)} mm, "
                     f"MRI-in-CTlung {r['mri_inside_ctlung']:.2f}, CT-lung covered {r['ctlung_covered_by_mri']:.2f}, "
                     f"similarity-scale diag {info['similarity_diag']['scale']:.3f}", fontsize=9)
        plt.tight_layout(); fig.savefig(os.path.join(out_dir, 'summary_3plane.png'), dpi=110); plt.close(fig)
        print(subj, 'overlays done:', out_dir)


def cmd_phases(subjects):
    """Diagnostic: which 4D-CT phase best matches the MRI EI bin (low-res rigid MI per phase, fixed orientation)."""
    for subj in subjects:
        path, ei, _, _ = _spec(subj)
        out_dir = os.path.join(REPO_OUT, subj)
        perm, flips = parse_orient(*choose_orientation(subj))
        v = load_mri_ei(path, ei); a = orient_mri(v, perm, flips)
        rows = []
        for tag in sorted(os.path.basename(f)[:-7] for f in glob.glob(os.path.join(EXT, subj, 'ct', 'phase*.nii.gz')) + [ct_paths(subj, 'average')]):
            ct, lung_img, roi_img, center = load_ct_and_masks(subj, tag)
            ct_lo = resample_iso(ct_windowed(ct), MRI_VOX); roi_lo = resample_iso(roi_img, MRI_VOX, sitk.sitkNearestNeighbor, 0)
            lung_lo = resample_iso(lung_img, MRI_VOX, sitk.sitkNearestNeighbor, 0)
            mov, mmask = mri_to_sitk(a, center)
            tf, mi, _ = register_rigid(ct_lo, mov, roi_lo, center, iters=200, shrink=(2, 1), sigmas=(2, 0))
            c_in, c_cov = containment(lung_lo, mov, mmask, tf)
            lung_l = float(sitk.GetArrayFromImage(lung_img).sum() * np.prod(ct.GetSpacing()) / 1e6)
            rows.append(dict(phase=tag, ct_lung_L=round(lung_l, 3), mi=mi, contain=c_in, cover=c_cov))
            print(subj, rows[-1])
        with open(os.path.join(out_dir, 'phase_scores.csv'), 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)


def _display_grid(ct, lung, spacing=1.5, margin_mm=25):
    """CT resampled to `spacing` iso, cropped to lung bbox + margin (for per-bin resampling)."""
    L = sitk.GetArrayFromImage(lung).astype(bool)
    kk, jj, ii = [np.where(L.any(axis=ax))[0] for ax in ((1, 2), (0, 2), (0, 1))]
    sp = np.array(ct.GetSpacing()); m = np.round(margin_mm / sp).astype(int)
    lo = [max(ii[0] - m[0], 0), max(jj[0] - m[1], 0), max(kk[0] - m[2], 0)]
    hi = [min(ii[-1] + m[0] + 1, ct.GetSize()[0]), min(jj[-1] + m[1] + 1, ct.GetSize()[1]), min(kk[-1] + m[2] + 1, ct.GetSize()[2])]
    crop = ct[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]]
    return resample_iso(crop, spacing, sitk.sitkLinear)


def cmd_video(subjects, fps=4):
    """16-bin video: every bin of the registered stack through the SAME rigid transform, on CT.
    Frame = 6 axial + 6 coronal slices spanning the lung. mp4 in repo outputs, gif on Ext."""
    import imageio.v2 as imageio
    import subprocess
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    for subj in subjects:
        path, ei, phase0, prov = _spec(subj)
        out_dir = os.path.join(REPO_OUT, subj); ext_dir = os.path.join(EXT, subj)
        with open(os.path.join(out_dir, 'register.json')) as fh:
            info = json.load(fh)
        phase = info.get('ct_phase', phase0)
        perm, flips = parse_orient(info['orientation']['perm'], info['orientation']['flips'])
        tf = sitk.ReadTransform(os.path.join(ext_dir, 'rigid.tfm'))
        ct = sitk.ReadImage(ct_paths(subj, phase), sitk.sitkFloat32)
        lung = sitk.ReadImage(os.path.join(ext_dir, 'ct_lungmask.nii.gz'))
        mov_ref = sitk.ReadImage(os.path.join(ext_dir, 'mri_ei_native.nii.gz'))   # carries origin/spacing
        grid = _display_grid(ct, lung)
        lung_g = sitk.Resample(lung, grid, sitk.Transform(), sitk.sitkNearestNeighbor, 0)
        Cg = sitk.GetArrayFromImage(grid); Lg = sitk.GetArrayFromImage(lung_g).astype(bool)
        keep = (Cg > -500) | ndi.binary_dilation(Lg, iterations=10)
        stack = np.abs(load_stack(path)).astype(np.float32)
        norm = np.percentile(stack[ei], 99.5); stack = np.clip(stack / norm, 0, 1.5)
        nb = stack.shape[0]
        ks = _contig_run(Lg.sum(axis=(1, 2)), 0.10); js = _contig_run(Lg.sum(axis=(0, 2)), 0.05)
        ksel = np.linspace(ks[0] + 0.1 * len(ks), ks[-1] - 0.1 * len(ks), 6).astype(int)
        jsel = np.linspace(js[0] + 0.1 * len(js), js[-1] - 0.1 * len(js), 6).astype(int)
        frames = []
        for b in range(nb):
            a = orient_mri(stack[b], perm, flips)
            mov = sitk.GetImageFromArray(np.ascontiguousarray(a)); mov.CopyInformation(mov_ref)
            Mg = sitk.GetArrayFromImage(sitk.Resample(mov, grid, tf, sitk.sitkLinear, 0.0, sitk.sitkFloat32)) * keep
            tiles_ax = [_overlay_rgb(Cg[k], Mg[k], Lg[k]) for k in ksel]
            tiles_co = [_overlay_rgb(Cg[:, j, :][::-1], Mg[:, j, :][::-1], Lg[:, j, :][::-1]) for j in jsel]
            fig, axs = plt.subplots(2, 6, figsize=(16, 2.6 * (tiles_ax[0].shape[0] / tiles_ax[0].shape[1] + tiles_co[0].shape[0] / tiles_co[0].shape[1]) + 0.6))
            for a_, t in zip(axs[0], tiles_ax):
                a_.imshow(t); a_.axis('off')
            for a_, t in zip(axs[1], tiles_co):
                a_.imshow(t); a_.axis('off')
            tag = 'EI' if b == ei else ''
            fig.suptitle(f'{subj}  bin {b:02d}/{nb - 1} {tag}   CT {phase} + Xe gas (rigid, all bins same transform)   top: axial (A up, L right)   bottom: coronal (S up, L right)', fontsize=10)
            plt.tight_layout(pad=0.2, rect=(0, 0, 1, 0.96))
            fig.canvas.draw()
            fr = np.asarray(fig.canvas.buffer_rgba())[..., :3].copy(); plt.close(fig)
            fr = fr[:fr.shape[0] // 2 * 2, :fr.shape[1] // 2 * 2]      # even dims for h264
            frames.append(fr)
        gif = os.path.join(ext_dir, 'bins_video.gif'); imageio.mimsave(gif, frames, duration=1.0 / fps, loop=0)
        mp4 = os.path.join(out_dir, 'bins_video.mp4')
        fr_dir = os.path.join(ext_dir, 'bins_video_frames'); os.makedirs(fr_dir, exist_ok=True)
        for i, fr in enumerate(frames):
            imageio.imwrite(os.path.join(fr_dir, f'bin{i:02d}.png'), fr)
        subprocess.run(['ffmpeg', '-y', '-loglevel', 'error', '-framerate', str(fps), '-i', os.path.join(fr_dir, 'bin%02d.png'),
                        '-c:v', 'libx264', '-pix_fmt', 'yuv420p', '-crf', '20', mp4], check=True)
        print(subj, 'video:', mp4, '| gif:', gif)


if __name__ == '__main__':
    cmd = sys.argv[1]
    subs = sys.argv[2:] or list(SUBJECTS)
    subs = [x for x in subs if '=' not in x]
    orient = None; phase_override = None
    if cmd == 'refine' and len(sys.argv) > 2 and '=' in sys.argv[-1]:  # refine 01BB perm=012 flips=+-+ [phase=phase37]
        kv = dict(x.split('=') for x in sys.argv[3:]); subs = [sys.argv[2]]
        orient = (kv['perm'], kv['flips']) if 'perm' in kv else None; phase_override = kv.get('phase')
    if cmd in ('search', 'all'):
        cmd_search(subs)
    if cmd in ('refine', 'all'):
        cmd_refine(subs, orient, phase_override)
    if cmd in ('overlay', 'all'):
        cmd_overlay(subs)
    if cmd == 'phases':
        cmd_phases(subs)
    if cmd in ('video', 'all'):
        cmd_video(subs)
