"""Paper F: coupled-channel complex-scaled Dirac + screening + field.

Each kappa-channel is a RadialDiracBasis on a SHARED Legendre grid
(same nq, xmax for every channel -> identical x,W -> consistent
cross-channel <r> integrals).  The field F z couples channels via the
relativistic dipole coefficients A_L (large) and A_S (small).

Per-channel diagonal block = the shifted complex-symmetric Dirac pencil
of Paper B (cs_core.cs_pencil): [[Vth, c ph K],[c ph K^T, Vth-2c^2]],
Vth = e^{-i th} <f|V_screened|f>, K = <f|O+ f>.  Field-free (F=0) it is
block-diagonal, each block reproducing exact Dirac-Coulomb (Sommerfeld).

Field block between channels kappa', kappa (F z, scales e^{+i th}):
  F e^{+i th} [[A_L(k',k) R, 0],[0, A_S(k',k) R]],  R = <f_k'| r |f_k>.
z is diagonal in the large/small space (commutes with beta), so it
never couples large-to-small; A_S uses the (-k',-k) spin-angular parity.
"""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(HERE))

from dirac import RadialDiracBasis, C_AU, exact_dirac_coulomb  # noqa: E402
from angular import A_L, A_S  # noqa: E402

C = C_AU


class CoupledDiracStark:
    def __init__(self, nb, kappas, Z, beta, m2=1, nq=None, xmax=None,
                 c=C_AU, cert_tol=1e-2):
        # cert_tol is loose because cs_pencil uses ONLY f, V and
        # K=<f|O+f>; the enforced c2/c3 (1/r^2 operator-product certs,
        # unused here) degrade on the Legendre grid at high Z.  The
        # certs that matter -- c1 (orthonormality) and c4 (IBP) -- stay
        # ~1e-12 and are re-checked in self.cert below.
        self.nb, self.kappas, self.Z, self.beta = nb, list(kappas), Z, beta
        self.m2, self.c = m2, c
        self.nchan = len(kappas)
        nq = nq or max(6000, 20 * nb)
        # shared grid: fixed xmax across channels (legendre branch uses
        # x = 0.5 X (t+1), independent of kappa) -> identical grid.
        lmax = max(abs(k) for k in kappas)
        xmax = xmax or (4.0 * nb + 24.0 * (lmax + 1) + 240.0)
        self.B = []
        gam = []
        for k in kappas:
            g = float(np.sqrt(k * k - (Z / c) ** 2))
            gam.append(g)
            B = RadialDiracBasis(nb, k, Z, beta, nquad=nq, xmax=xmax,
                                 c=c, cert_tol=cert_tol, mu=g,
                                 quad="legendre")
            self.B.append(B)
        self.gamma = gam
        self.cert = max(max(B.certs[0], B.certs[3]) for B in self.B)
        # shared grid handles
        self.r = self.B[0].r
        self.W = self.B[0].W
        # dipole radial blocks R[(i,j)] = <f_ki | r | f_kj> (shared grid)
        self._R = {}

    def _dipole_R(self, i, j):
        key = (i, j)
        if key not in self._R:
            fi, fj = self.B[i].f, self.B[j].f
            self._R[key] = np.einsum('ak,bk,k->ab', fi, self.r * fj, self.W)
        return self._R[key]

    def _channel_pencil(self, i, theta, lam):
        """Shifted complex-symmetric Dirac pencil for channel i."""
        B = self.B[i]
        c, Z, r = self.c, self.Z, B.r
        ph = np.exp(-1j * theta)
        f = B.f
        K = np.einsum('ak,bk,k->ab', f, (2 * c * B.h), B.W)   # <f|O+ f>
        if lam == 0.0:
            V = -Z * np.einsum('ak,bk,k->ab', f, (1.0 / r) * f, B.W)
        else:
            scr = np.exp(-lam * r * np.exp(1j * theta))
            V = -Z * np.einsum('ak,k,bk,k->ab', f, scr / r, f, B.W)
        Vth = ph * V
        n = B.nb
        Hs = np.block([[Vth, c * ph * K],
                       [c * ph * K.T, Vth - 2 * c * c * np.eye(n)]])
        return 0.5 * (Hs + Hs.T)

    def build(self, F, theta, lam=0.0):
        nb, nc = self.nb, self.nchan
        N = 2 * nb * nc
        H = np.zeros((N, N), dtype=complex)
        for i in range(nc):
            blk = self._channel_pencil(i, theta, lam)
            H[2 * nb * i:2 * nb * (i + 1),
              2 * nb * i:2 * nb * (i + 1)] = blk
        if F != 0.0:
            phf = np.exp(1j * theta)
            for i in range(nc):
                for j in range(nc):
                    if i == j:
                        continue
                    aL = A_L(self.kappas[i], self.kappas[j], self.m2)
                    aS = A_S(self.kappas[i], self.kappas[j], self.m2)
                    if abs(aL) < 1e-14 and abs(aS) < 1e-14:
                        continue
                    R = self._dipole_R(i, j)
                    blk = np.zeros((2 * nb, 2 * nb), dtype=complex)
                    blk[:nb, :nb] = F * phf * aL * R
                    blk[nb:, nb:] = F * phf * aS * R
                    H[2 * nb * i:2 * nb * (i + 1),
                      2 * nb * j:2 * nb * (j + 1)] = blk
        return 0.5 * (H + H.T)


def channel_sommerfeld(kappa, Z, nmax=6, c=C_AU):
    """Exact Dirac-Coulomb binding energies for a channel (for gates)."""
    out = {}
    for n in range(abs(kappa), nmax + 1):
        try:
            out[n] = exact_dirac_coulomb(n, kappa, Z, c)
        except Exception:
            pass
    return out
