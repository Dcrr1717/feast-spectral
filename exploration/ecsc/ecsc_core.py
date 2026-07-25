"""Paper E (ECSC + Stark): coupled-channel engine.

ECSC potential (quantum/dense plasma screening):
    V(r) = -(Z/r) e^{-lam r} cos(lam r)
         = -(Z/2r) [ e^{-mu+ r} + e^{-mu- r} ],   mu+- = lam (1 +- i)

i.e. HALF-SUM OF TWO YUKAWAS AT COMPLEX SCREENING PARAMETERS -- so the
complex-scaled matrix elements reuse the Paper D machinery verbatim.
Decay of e^{-mu+- r e^{i th}} requires Re(mu+- e^{i th}) > 0, i.e.
cos th - sin th > 0  =>  THETA < PI/4 (hard analytic cap, asserted).

Class StarkChannelsE = Paper D's StarkChannels + scr='yukawa'|'ecsc'.
scr='yukawa' must reproduce d0_architecture.StarkChannels exactly
(regression-certified in e0 gates).
"""

import sys
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from basis import _q_functions  # noqa: E402

THETA_CAP = np.pi / 4          # ECSC analytic rotation cap
_GATED = False


class StarkChannelsE:
    """Coupled-channel complex-scaled Stark + screened Coulomb
    (Yukawa or ECSC screening)."""

    def __init__(self, nb, lmax, beta, nquad=None):
        self.nb, self.lmax, self.beta = nb, lmax, beta
        xmax = 4.0 * nb + 24.0 * (lmax + 1) + 120.0
        nq = nquad or max(3000, 16 * nb)
        t, w = leggauss(nq)
        self.x = 0.5 * xmax * (t + 1.0)
        self.w = 0.5 * xmax * w
        self.q, self.dq = {}, {}
        for l in range(lmax + 1):
            self.q[l], self.dq[l] = _q_functions(nb, l, self.x)
        self.cert = max(
            float(np.max(np.abs((self.q[l] * self.w) @ self.q[l].T
                                - np.eye(nb))))
            for l in range(lmax + 1))

    def _screen(self, lam, r, theta, scr):
        """Screening factor f(r) with r -> r e^{i theta} applied."""
        if lam == 0.0:
            return 1.0
        z = r * np.exp(1j * theta)
        if scr == "yukawa":
            return np.exp(-lam * z)
        if scr == "ecsc":
            if theta >= THETA_CAP:
                raise ValueError(
                    f"ECSC rotation angle {theta} >= pi/4: screened "
                    "blocks diverge (e^{-lam z} cos(lam z) grows)")
            mup, mum = lam * (1 + 1j), lam * (1 - 1j)
            return 0.5 * (np.exp(-mup * z) + np.exp(-mum * z))
        raise ValueError(f"unknown screening '{scr}'")

    def build(self, Z, lam, F, theta, scr="ecsc"):
        nb, lmax, b = self.nb, self.lmax, self.beta
        N = (lmax + 1) * nb
        H = np.zeros((N, N), dtype=complex)
        ph1, ph2 = np.exp(-1j * theta), np.exp(-2j * theta)
        r = self.x / (2 * b)
        f = self._screen(lam, r, theta, scr)
        for l in range(lmax + 1):
            qw = self.q[l] * self.w
            T = 2 * b * b * ((self.dq[l] * self.w) @ self.dq[l].T)
            if l > 0:
                T += 2 * b * b * l * (l + 1) * ((qw / self.x**2) @ self.q[l].T)
            T = 0.5 * (T + T.T)
            V = -Z * 2 * b * ((qw * f / self.x) @ self.q[l].T)
            Vth = ph1 * 0.5 * (V + V.T)
            i0 = l * nb
            H[i0:i0 + nb, i0:i0 + nb] = ph2 * T + Vth
        for l in range(lmax):
            cl = (l + 1) / np.sqrt((2 * l + 1) * (2 * l + 3))
            R = ((self.q[l] * self.w * r) @ self.q[l + 1].T)
            blk = F * cl * np.exp(1j * theta) * R
            i0, j0 = l * nb, (l + 1) * nb
            H[i0:i0 + nb, j0:j0 + nb] = blk
            H[j0:j0 + nb, i0:i0 + nb] = blk.T
        return H


def track_point(S, F, lam, pred, radius, tol, thetas, scr="ecsc",
                imcut=1e-8):
    """Eigenvalue nearest `pred` within `radius`, theta-stable to `tol`.

    The second angle pairs to the FIRST angle's pick (not to `pred`):
    independent nearest-to-pred picks at each angle can select
    different states that accidentally pair within `tol` (caught by
    the Paper E round-2 audit at the n=2 blue lam=0.125 point)."""
    ev = np.linalg.eigvals(S.build(1.0, lam, F, thetas[0], scr=scr))
    ev = ev[(np.abs(ev - pred) < radius) & (ev.imag < imcut)]
    if not ev.size:
        return None
    e0 = ev[np.argmin(np.abs(ev - pred))]
    ev = np.linalg.eigvals(S.build(1.0, lam, F, thetas[1], scr=scr))
    ev = ev[(np.abs(ev - e0) < radius) & (ev.imag < imcut)]
    if not ev.size:
        return None
    e1 = ev[np.argmin(np.abs(ev - e0))]
    d = abs(e1 - e0)
    return (e1, d) if d <= tol else None


# ---------------- import-time certification gate ------------------------
# Every production run importing this engine passes the hydrogen limit
# or dies (P-E0.0 discipline, enforced here so the claim "a failed
# check aborts every run" is literally true).
def _hydrogen_gate():
    global _GATED
    if _GATED:
        return
    S = StarkChannelsE(60, 2, 0.35)
    ev = np.sort(np.linalg.eigvals(S.build(1.0, 0.0, 0.0, 0.0)).real)
    exact = [-0.5, -0.125, -0.125, -1.0 / 18, -1.0 / 18]
    dev = max(abs(a - b) for a, b in zip(ev[:5], exact))
    if not (dev < 1e-9):
        raise RuntimeError(
            f"ecsc_core hydrogen gate FAILED (dev={dev:.1e}) -- no "
            "production run may proceed on a failed assembly")
    _GATED = True


_hydrogen_gate()
