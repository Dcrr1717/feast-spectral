"""Radial Dirac-Coulomb solver on orthonormal generalized-Laguerre bases.

Three discretizations of the radial Dirac operator

    h_kappa = [[ V + c^2,  c*O+ ],          O+/- = -/+ d/dr + kappa/r,
               [ c*O-,     V - c^2 ]]       V = -Z/r,

are provided: NAIVE (independent Laguerre sets for both components), RKB
(restricted kinetic balance) and DKB (dual kinetic balance, Shabaev et al.,
PRL 93, 130405 (2004)), all built on the L2-orthonormal Laguerre functions
phi_j ~ r^(ell+1) e^(-beta r) with ell(ell+1) = kappa(kappa+1).

Operator identities used (order matters!):

    O- O+ u = -u'' + kappa(kappa-1) u / r^2
    O+ O- u = -u'' + kappa(kappa+1) u / r^2

CERTIFIED ASSEMBLY: every build checks, at run time, the internal
identities that the discrete matrices must satisfy exactly (up to
quadrature/roundoff),

    C1:  <phi_i | phi_j>                  = delta_ij
    C2:  <phi | O+O- phi>  (via d2f)      = <O- phi | O- phi>   (= 4c^2 <g|g>)
    C3:  <phi | O-O+ phi>  (via d2f)      = <O+ phi | O+ phi>   (= 4c^2 <h|h>)

C2/C3 certify simultaneously the second-derivative recurrence, the
quadrature near the r->0 singularity, the integration-by-parts boundary
terms, AND the centrifugal coefficients kappa(kappa+-1) -- i.e. exactly the
class of bug that a wrong operator-product expansion introduces.
"""

from __future__ import annotations

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import gammaln
from scipy.linalg import eigh

C_AU = 137.035999084          # 1/alpha, CODATA 2022

_LEGGAUSS_CACHE: dict[int, tuple[np.ndarray, np.ndarray]] = {}


def exact_dirac_coulomb(n: int, kappa: int, Z: float, c: float = C_AU) -> float:
    """Point-nucleus Dirac-Coulomb binding energy E - c^2 (Hartree)."""
    g = np.sqrt(kappa * kappa - (Z / c) ** 2)
    nr = n - abs(kappa)
    return c * c / np.sqrt(1.0 + (Z / c / (nr + g)) ** 2) - c * c


def exact_ladder(kappa: int, Z: float, kmax: int = 12) -> np.ndarray:
    """First kmax exact binding energies of channel kappa (ascending)."""
    n0 = abs(kappa) if kappa < 0 else kappa + 1
    return np.array([exact_dirac_coulomb(n, kappa, Z) for n in range(n0, n0 + kmax)])


def laguerre_funcs(nb: int, ell: int, x: np.ndarray):
    """Orthonormal Laguerre functions q_j(x) with first and second
    derivatives, via the numerically stable orthonormal recurrence."""
    a = 2 * ell + 2
    q = np.zeros((nb, x.size))
    q[0] = np.exp((ell + 1) * np.log(x) - 0.5 * x - 0.5 * gammaln(a + 1))
    if nb > 1:
        q[1] = (a + 1 - x) * q[0] / np.sqrt(1.0 * (1 + a))
    for j in range(1, nb - 1):
        q[j + 1] = ((2 * j + a + 1 - x) * q[j]
                    - np.sqrt(j * (j + a)) * q[j - 1]) / np.sqrt((j + 1.0) * (j + 1 + a))
    # first derivatives (complete them BEFORE any second derivative uses them)
    dq = np.zeros_like(q)
    for j in range(nb):
        dq[j] = ((ell + 1 + j) / x - 0.5) * q[j]
        if j > 0:
            dq[j] -= np.sqrt(j * (j + a)) / x * q[j - 1]
    # second derivatives from the *complete* dq
    d2q = np.zeros_like(q)
    for j in range(nb):
        A = (ell + 1 + j) / x - 0.5
        Ap = -(ell + 1 + j) / x ** 2
        d2q[j] = Ap * q[j] + A * dq[j]
        if j > 0:
            B = np.sqrt(j * (j + a)) / x
            Bp = -np.sqrt(j * (j + a)) / x ** 2
            d2q[j] -= Bp * q[j - 1] + B * dq[j - 1]
    return q, dq, d2q


class RadialDiracBasis:
    """Laguerre machinery for one (kappa, Z, beta) channel, with certified
    quadrature and the DKB/RKB/naive Hamiltonians."""

    def __init__(self, nb: int, kappa: int, Z: float, beta: float,
                 nquad: int | None = None, xmax: float | None = None,
                 c: float = C_AU, cert_tol: float = 1e-7):
        # cert_tol calibration: structural assembly bugs show up as O(0.1-1)
        # relative violations; converged-quadrature noise sits at <~1e-9
        # (the O+ channel of ell=0 has a non-vanishing integrand at r=0 and
        # is the slowest to converge). 1e-7 separates the two cleanly.
        assert kappa != 0
        self.nb, self.kappa, self.Z, self.beta, self.c = nb, kappa, Z, beta, c
        self.ell = kappa if kappa > 0 else -kappa - 1

        nq = nquad or max(4000, 16 * nb)
        if nq not in _LEGGAUSS_CACHE:
            _LEGGAUSS_CACHE[nq] = leggauss(nq)
        t, w0 = _LEGGAUSS_CACHE[nq]
        X = xmax or (4.0 * nb + 24.0 * (self.ell + 1) + 200.0)
        x = 0.5 * X * (t + 1.0)
        w = 0.5 * X * w0

        q, dq, d2q = laguerre_funcs(nb, self.ell, x)
        b = beta
        self.r = x / (2 * b)
        self.W = w / (2 * b)                       # dr-measure weights
        self.f = np.sqrt(2 * b) * q
        self.df = (2 * b) ** 1.5 * dq
        self.d2f = (2 * b) ** 2.5 * d2q
        r = self.r
        self.g = (self.df + (kappa / r) * self.f) / (2 * c)     # (1/2c) O- f
        self.h = (-self.df + (kappa / r) * self.f) / (2 * c)    # (1/2c) O+ f
        # operator products -- ORDER MATTERS:
        self.OmOp = -self.d2f + kappa * (kappa - 1) / r ** 2 * self.f  # O-O+ f
        self.OpOm = -self.d2f + kappa * (kappa + 1) / r ** 2 * self.f  # O+O- f
        self.V = -Z / r

        # ---- certified assembly ------------------------------------
        ME = self._me
        Kg = 4 * c * c * ME(self.g, self.g)      # <O- f | O- f>
        Kh = 4 * c * c * ME(self.h, self.h)      # <O+ f | O+ f>
        c1 = np.max(np.abs(ME(self.f, self.f) - np.eye(nb)))
        c2 = np.max(np.abs(ME(self.f, self.OpOm) - Kg)) / max(1.0, np.max(np.abs(Kg)))
        c3 = np.max(np.abs(ME(self.f, self.OmOp) - Kh)) / max(1.0, np.max(np.abs(Kh)))
        # I1 (IBP identity <f_i|O+ f_j> = <O- f_i|f_j>): catches O+/O- sign
        # errors and boundary-term violations in one shot.
        c4 = np.max(np.abs(ME(self.f, self.h) - ME(self.g, self.f)))
        self.certs = (c1, c2, c3, c4)            # relative to natural scale
        if max(self.certs) > cert_tol:
            raise RuntimeError(f"assembly certification failed: {self.certs}")

    # -- integrals ----------------------------------------------------
    def _me(self, u, v):
        return np.einsum('ik,jk,k->ij', u, v, self.W)

    def _mev(self, u, p, v):
        return np.einsum('ik,jk,k->ij', u, p * v, self.W)

    # -- solvers --------------------------------------------------------
    def _ingap(self, H, S):
        ev = eigh(0.5 * (H + H.T), 0.5 * (S + S.T), eigvals_only=True)
        ev = np.sort(ev) - self.c ** 2
        return ev[(ev > -2 * self.c ** 2 + 1e-6) & (ev < -1e-9)]

    def solve_dkb(self) -> np.ndarray:
        """DKB, Shabaev sign convention (PRL 93, 130405):

            Psi+_j = (  f_j            , +(1/2c) O- f_j )
            Psi-_j = ( -(1/2c) O+ f_j  ,  f_j           )

        The MINUS sign is load-bearing: it flips one component only, which
        changes the span (not a basis-sign congruence). It follows from the
        positronic limit P = c O+ Q / (E - c^2 - V) -> -(1/2c) O+ Q at
        E ~ -c^2, and it is what gives the H-- block the correct
        charge-conjugate nonrelativistic limit (kinetic +(1/2)<O+f|O+f>).
        With s=+1 the positronic limit has the wrong-sign kinetic term and
        spurious states appear deep inside the gap (observed empirically).
        """
        c, V = self.c, self.V
        f, g, h = self.f, self.g, self.h
        ME, MEv = self._me, self._mev
        Sff = ME(f, f)
        Spp = Sff + ME(g, g)
        Spm = -ME(f, h) + ME(g, f)
        Smm = ME(h, h) + Sff
        Hpp = (MEv(f, V, f) + c * c * Sff
               + 4 * c * c * ME(g, g)                   # <f|cO+|g> + <g|cO-|f>
               + MEv(g, V, g) - c * c * ME(g, g))
        Hmm = (MEv(h, V, h) + c * c * ME(h, h)
               - 4 * c * c * ME(h, h)                   # s=-1: <-h|cO+|f>+<f|cO-|-h>
               + MEv(f, V, f) - c * c * Sff)
        Hpm = (-(MEv(f, V, h) + c * c * ME(f, h))       # <f|(V+c^2)|-h>
               + 2 * c * c * ME(f, h)                   # <f|cO+|f> = 2c^2<f|h>
               - 0.5 * ME(g, self.OmOp)                 # <g|cO-|-h> = -(1/2)<g|O-O+f>
               + MEv(g, V, f) - c * c * ME(g, f))
        H = np.block([[Hpp, Hpm], [Hpm.T, Hmm]])
        S = np.block([[Spp, Spm], [Spm.T, Smm]])
        return self._ingap(H, S)

    def solve_rkb(self) -> np.ndarray:
        """RKB: large in {f_j}, small in {O- f_j} with its Gram metric."""
        c, V = self.c, self.V
        f = self.f
        chi = 2 * c * self.g                          # O- f
        ME, MEv = self._me, self._mev
        Sll = ME(f, f)
        Sss = ME(chi, chi)
        Hll = MEv(f, V, f) + c * c * Sll
        # <f_i | c O+ | chi_j> = c <O- f_i | chi_j> = c (Sss)_ij  (by parts)
        Hss = MEv(chi, V, chi) - c * c * Sss
        H = np.block([[Hll, c * Sss], [c * Sss.T, Hss]])
        S = np.block([[Sll, 0 * Sll], [0 * Sll, Sss]])
        return self._ingap(H, S)

    def solve_naive(self) -> np.ndarray:
        """Both components in the SAME Laguerre set (pollution showcase)."""
        c, V = self.c, self.V
        f, df = self.f, self.df
        kappa, r = self.kappa, self.r
        ME, MEv = self._me, self._mev
        Opf = -df + kappa / r * f
        Sff = ME(f, f)
        H = np.block([[MEv(f, V, f) + c * c * Sff, c * ME(f, Opf)],
                      [c * ME(f, Opf).T, MEv(f, V, f) - c * c * Sff]])
        S = np.block([[Sff, 0 * Sff], [0 * Sff, Sff]])
        return self._ingap(H, S)


def validate(kappa: int, Z: float, nb: int = 60, betas=None,
             method: str = "dkb", kcheck: int = 4):
    """Best-over-beta validation of the lowest kcheck states vs exact.

    Returns dict with per-state errors, the variational-violation flag
    (any Galerkin level strictly BELOW its exact counterpart), and certs.
    """
    betas = betas if betas is not None else [0.5 * max(Z, 1.0), 0.9 * max(Z, 1.0), 1.4 * max(Z, 1.0)]
    ladder = exact_ladder(kappa, Z)
    best = None
    for b in betas:
        B = RadialDiracBasis(nb, kappa, Z, b)
        gap = getattr(B, f"solve_{method}")()
        if gap.size < kcheck:
            continue
        errs = gap[:kcheck] - ladder[:kcheck]
        score = abs(errs[0])
        if best is None or score < best["err"][0].__abs__():
            best = dict(beta=b, err=errs, levels=gap[:kcheck], certs=B.certs,
                        viol=bool(np.any(errs < -1e-8 * np.maximum(1.0, np.abs(ladder[:kcheck])))))
    return best
