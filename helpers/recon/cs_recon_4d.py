"""4D temporal compressed sensing on the validated FINUFFT operator.

Extends the static CS (cs_recon.py) to a B-bin respiratory cine. Single coil
(027JC is nch=1), so the 4D forward operator is just B independent copies of the
proven 3D operator; bins differ ONLY by their per-sample weight w_b = dcf * m_b
(handoff decision D3). The only coupling between bins is the temporal regularizer.

Objective (X has shape (B, N, N, N)):

    min_X  sum_b || sqrt(c*w_dcf*m_b) (A x_b - y) ||^2
         + lam_s * || W_spatial X ||_1     (per-bin db4 wavelet  -- the 3D winner)
         + lam_t * || D_t X ||_1           (circular temporal TV across bins -- NEW)

Stage 3 baseline = drop the temporal term: B independent wavelet_recon calls
(reuses cs_recon.wavelet_recon verbatim). Stage 4 = joint PDHG with both priors.

lam_s / lam_t are passed RELATIVE to the p99 coefficient magnitude of the Stage-3
baseline (same convention as cs_recon.T_RELS; objective-scaled lambda is a measured
no-op -- see AGENTS.md). DCF is baked into A, so PDHG does not stall (AGENTS.md
conditioning trap is about gradient solvers WITHOUT DCF).
"""

import numpy as np
import sigpy as sp

import asap_recon as ar
from cs_recon import FinufftForward, wavelet_recon

N_DEFAULT = 100
EPS = 1e-7


# ----------------------------------------------------------------- DCF / norm

def dcf_and_norm(traj_rad, N, eps=EPS):
    """DCF weights w (1/density), normalization c = 1/sqrt(maxeig(A_w^H A_w)).

    Identical math to cs_recon.py:134-147 but factored out for reuse. Computed
    once on the (single, shared) trajectory -- bin membership does not change it.
    """
    A_raw = FinufftForward(traj_rad, N, eps=eps)
    M = traj_rad.shape[0]
    dens = np.abs(A_raw(A_raw.H(np.ones(M, dtype=complex))))
    w = 1.0 / np.clip(dens, dens.max() * 1e-4, None)
    w /= w.mean()
    Aw = sp.linop.Multiply((M,), np.sqrt(w)) * A_raw
    L = sp.app.MaxEig(Aw.H * Aw, dtype=np.complex128, show_pbar=False).run()
    return w, 1.0 / np.sqrt(L), L


# ----------------------------------------------------------------- 4D operators

class Finufft4D(sp.linop.Linop):
    """A_4D: (B,N,N,N) image stack -> (B,M) samples. Bin-independent, shared traj."""

    def __init__(self, traj_rad, B, n, eps=EPS, scale=1.0):
        self.traj_rad = traj_rad
        self.B = B
        self.n = n
        self.eps = eps
        self.scale = scale
        super().__init__((B, traj_rad.shape[0]), (B, n, n, n))

    def _apply(self, input):
        out = np.empty((self.B, self.traj_rad.shape[0]), dtype=np.complex128)
        for b in range(self.B):
            out[b] = self.scale * ar.forward(self.traj_rad, input[b], eps=self.eps)
        return out

    def _adjoint_linop(self):
        return _Finufft4DAdjoint(self.traj_rad, self.B, self.n, self.eps, self.scale)


class _Finufft4DAdjoint(sp.linop.Linop):
    def __init__(self, traj_rad, B, n, eps, scale):
        self.traj_rad = traj_rad
        self.B = B
        self.n = n
        self.eps = eps
        self.scale = scale
        super().__init__((B, n, n, n), (B, traj_rad.shape[0]))

    def _apply(self, input):
        out = np.empty((self.B, self.n, self.n, self.n), dtype=np.complex128)
        for b in range(self.B):
            out[b] = self.scale * ar.adjoint(self.traj_rad, input[b], self.n, eps=self.eps)
        return out

    def _adjoint_linop(self):
        return Finufft4D(self.traj_rad, self.B, self.n, self.eps, self.scale)


def temporal_diff_op(ishape):
    """Circular first difference along the bin axis 0: (D X)[b] = X[(b+1)%B] - X[b].

    Port of @TVOPDt (workspace/codes/2025_Xe129_CS/@TVOPDt). Built from sigpy
    primitives so the adjoint (-circular divergence) is exact and auto-derived.
    """
    return sp.linop.Circshift(ishape, [-1], axes=[0]) - sp.linop.Identity(ishape)


def spatial_wavelet_op(ishape, wave_name="db4"):
    """Per-bin db4 wavelet (axes 1,2,3) -- treats bin axis 0 as batch. Orthogonal."""
    return sp.linop.Wavelet(ishape, axes=(1, 2, 3), wave_name=wave_name)


# ----------------------------------------------------------------- recons

def _weighted_op_and_data(traj_rad, y, w_dcf, c, m_samp_b, N):
    """Single-bin weighted, normalized operator + data for bin b. sqrt(c*w_dcf*m_b)."""
    M = traj_rad.shape[0]
    sw = c * np.sqrt(w_dcf) * np.sqrt(np.clip(m_samp_b, 0, None))
    A = sp.linop.Multiply((M,), sw) * FinufftForward(traj_rad, N, eps=EPS)
    return A, sw * y


def recon_4d_baseline(traj_grid, y, M_samp, N=N_DEFAULT, MS=ar.MS_DEFAULT,
                      lam_s_rel=0.01, max_iter=60, w_dcf=None, c=None, verbose=True):
    """Stage 3: B independent per-bin wavelet CS (no temporal coupling).

    traj_grid : (Nuniq,3) grid-index trajectory (MS-grid units).
    y         : (M,) complex samples (M = npts*ntotalilvs).
    M_samp    : (B, M) per-sample soft membership (binning.tile_to_samples output).
    MS        : grid-units matrix size of traj_grid (240 for the real dump).
    Returns (B,N,N,N) complex cine and the (w_dcf, c) used (so Stage 4 can reuse).
    """
    B, M = M_samp.shape
    traj_rad = ar.grid_to_radians(ar.tile_traj(np.asarray(traj_grid, float), M), MS)
    if w_dcf is None or c is None:
        w_dcf, c, _ = dcf_and_norm(traj_rad, N)
    y = np.ascontiguousarray(y, dtype=np.complex128)

    # threshold reference from a CG bin (most-populated) -> p99 of its wavelet coeffs
    W = spatial_wavelet_op((1, N, N, N))
    b_ref = int(np.argmax(M_samp.sum(1)))
    img_ref = ar.recon(traj_grid, y, sample_weights=w_dcf * M_samp[b_ref],
                       method="cg", MS=MS, IS=N, cg_iters=20)
    t_ref = float(np.percentile(np.abs(W(img_ref[None])), 99))
    lam_s = lam_s_rel * t_ref

    cine = np.empty((B, N, N, N), dtype=np.complex128)
    for b in range(B):
        A, y_n = _weighted_op_and_data(traj_rad, y, w_dcf, c, M_samp[b], N)
        cine[b] = wavelet_recon(A, y_n, lam_s, max_iter)
        if verbose:
            print(f"  baseline bin {b+1:2d}/{B}  |x| p99 "
                  f"{np.percentile(np.abs(cine[b]),99):.3e}")
    return cine, (w_dcf, c, t_ref)


def recon_4d_joint(traj_grid, y, M_samp, N=N_DEFAULT, MS=ar.MS_DEFAULT,
                   lam_s_rel=0.01, lam_t_rel=0.05, max_iter=80, baseline=None,
                   w_dcf=None, c=None, max_power_iter=15, verbose=True):
    """Stage 4: joint wavelet + circular temporal-TV via PDHG.

    baseline : optional (B,N,N,N) Stage-3 cine used to set both thresholds. If
    None, a baseline is computed first (recommended -- it is also the reference).
    """
    B, M = M_samp.shape
    traj_rad = ar.grid_to_radians(ar.tile_traj(np.asarray(traj_grid, float), M), MS)
    if w_dcf is None or c is None:
        w_dcf, c, _ = dcf_and_norm(traj_rad, N)
    y = np.ascontiguousarray(y, dtype=np.complex128)

    if baseline is None:
        if verbose:
            print("computing Stage-3 baseline for thresholds + reference ...")
        baseline, (w_dcf, c, _) = recon_4d_baseline(
            traj_grid, y, M_samp, N=N, MS=MS, lam_s_rel=lam_s_rel,
            max_iter=max(30, max_iter // 2), w_dcf=w_dcf, c=c, verbose=verbose)

    ishape = (B, N, N, N)
    W = spatial_wavelet_op(ishape)
    D = temporal_diff_op(ishape)
    lam_s = lam_s_rel * float(np.percentile(np.abs(W(baseline)), 99))
    lam_t = lam_t_rel * float(np.percentile(np.abs(D(baseline)), 99))
    if verbose:
        print(f"lam_s={lam_s:.3e} (rel {lam_s_rel})  lam_t={lam_t:.3e} (rel {lam_t_rel})")

    # normalized stacked forward operator: (B,N,N,N) -> (B,M)
    Wmat = (c * np.sqrt(w_dcf))[None, :] * np.sqrt(np.clip(M_samp, 0, None))
    A = sp.linop.Multiply((B, M), Wmat) * Finufft4D(traj_rad, B, N)
    y_n = Wmat * y[None, :]

    G = sp.linop.Vstack([W, D])
    proxg = sp.prox.Stack([sp.prox.L1Reg(W.oshape, lam_s),
                           sp.prox.L1Reg(D.oshape, lam_t)])
    x = sp.app.LinearLeastSquares(A, y_n, x=baseline.astype(np.complex128).copy(),
                                  G=G, proxg=proxg, max_iter=max_iter,
                                  max_power_iter=max_power_iter,
                                  show_pbar=verbose).run()
    return x, baseline
