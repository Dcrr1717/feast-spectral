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

import math

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import gammaln, roots_genlaguerre
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


def laguerre_funcs(nb: int, mu: float, x: np.ndarray):
    """Orthonormal Laguerre functions with REAL prefactor power mu > 0,

        q_j(x) = sqrt(j!/Gamma(j+2mu+1)) x^mu e^{-x/2} L_j^{(2mu)}(x),

    plus first and second derivatives (stable orthonormal recurrence).
    Integer case mu = ell+1 recovers the nonrelativistic basis; the
    cusp-adapted relativistic case uses mu = gamma = sqrt(kappa^2-(Z/c)^2)
    (Grant's L-spinor prefactor), for which the point-nucleus Dirac
    eigenfunctions are finite linear combinations of the q_j.
    """
    a = 2.0 * mu
    q = np.zeros((nb, x.size))
    q[0] = np.exp(mu * np.log(x) - 0.5 * x - 0.5 * gammaln(a + 1))
    if nb > 1:
        q[1] = (a + 1 - x) * q[0] / np.sqrt(1.0 * (1 + a))
    for j in range(1, nb - 1):
        q[j + 1] = ((2 * j + a + 1 - x) * q[j]
                    - np.sqrt(j * (j + a)) * q[j - 1]) / np.sqrt((j + 1.0) * (j + 1 + a))
    # first derivatives (complete them BEFORE any second derivative uses them)
    dq = np.zeros_like(q)
    for j in range(nb):
        dq[j] = ((mu + j) / x - 0.5) * q[j]
        if j > 0:
            dq[j] -= np.sqrt(j * (j + a)) / x * q[j - 1]
    # second derivatives from the *complete* dq
    d2q = np.zeros_like(q)
    for j in range(nb):
        A = (mu + j) / x - 0.5
        Ap = -(mu + j) / x ** 2
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
                 c: float = C_AU, cert_tol: float = 1e-7,
                 mu: float | None = None, quad: str = "legendre"):
        # cert_tol calibration: structural assembly bugs show up as O(0.1-1)
        # relative violations; converged-quadrature noise sits at <~1e-9
        # (the O+ channel of ell=0 has a non-vanishing integrand at r=0 and
        # is the slowest to converge). 1e-7 separates the two cleanly.
        # mu: prefactor power. Default ell+1 (nonrelativistic regularity);
        # pass mu = sqrt(kappa^2-(Z/c)^2) for the cusp-adapted basis. NOTE:
        # the C2/C3 boundary terms scale as r^(2mu-1), so for mu < 1/2
        # (Z ~> 118 at |kappa|=1) those certifications legitimately degrade.
        assert kappa != 0
        self.nb, self.kappa, self.Z, self.beta, self.c = nb, kappa, Z, beta, c
        self.ell = kappa if kappa > 0 else -kappa - 1
        self.mu = float(mu) if mu is not None else float(self.ell + 1)

        if quad == "cusp":
            # Generalized Gauss-Laguerre rule with weight x^(2mu-1) e^{-x}.
            # Every product of two basis-family functions carries exactly
            # x^(2mu) e^{-x} (or x^(2mu-1) for the 1/r elements), so the
            # overlap, Coulomb and O+- coupling integrands are POLYNOMIAL
            # times this weight: the rule is EXACT for nq >= nb + a few.
            # (1/r^2-type elements, i.e. C2/C3 and the DKB cross block,
            # remain inexact and require mu > 1/2 for convergence at all.)
            nq = nquad or (nb + 30)
            xk, wk = roots_genlaguerre(nq, 2.0 * self.mu - 1.0)
            # instrument envelope: scipy's Golub-Welsch weights break down
            # (NaN nodes/weights) for nq >~ 350-400 depending on alpha.
            # Fail loudly -- a NaN quadrature would otherwise blind the
            # certifications (NaN compares False against any tolerance).
            if not (np.isfinite(xk).all() and np.isfinite(wk).all()):
                raise RuntimeError(
                    f"cusp quadrature breakdown: roots_genlaguerre({nq}, "
                    f"{2 * self.mu - 1:.3f}) returned non-finite nodes/"
                    "weights; reduce nquad (envelope ~350) or use "
                    "quad='legendre'")
            x = xk
            # dx-measure weights for functions that already contain e^{-x/2}
            # each: W_dx = wk * e^{x} * x^{-(2mu-1)}, computed in logs.
            # wk underflows to 0 for far-tail nodes (xk >~ 709); those nodes
            # carry wk*poly contributions <~ 1e-300 and are dropped exactly.
            with np.errstate(divide="ignore"):
                w = np.exp(np.log(wk) + xk - (2.0 * self.mu - 1.0) * np.log(xk))
        else:
            nq = nquad or max(4000, 16 * nb)
            if nq not in _LEGGAUSS_CACHE:
                _LEGGAUSS_CACHE[nq] = leggauss(nq)
            t, w0 = _LEGGAUSS_CACHE[nq]
            X = xmax or (4.0 * nb + 24.0 * (self.ell + 1) + 200.0)
            x = 0.5 * X * (t + 1.0)
            w = 0.5 * X * w0

        q, dq, d2q = laguerre_funcs(nb, self.mu, x)
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
        # cusp rule: C1/C4 integrands are exact (polynomial x weight);
        # C2/C3 involve 1/r^2 and are legitimately inexact -> report only.
        enforced = (c1, c4) if quad == "cusp" else self.certs
        # NaN-proof: 'not (x <= tol)' fails on NaN, unlike 'x > tol'
        if not (max(enforced) <= cert_tol):
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

    def _rayleigh(self, Hs, S, x):
        num = math.fsum((x[:, None] * Hs * x[None, :]).ravel().tolist())
        den = math.fsum((x[:, None] * S * x[None, :]).ravel().tolist())
        return num / den

    def _ingap_shifted(self, Hs, S):
        """Gap eigenvalues of the SHIFTED pencil (H - c^2 S, S): these are
        binding energies directly, so the rest-mass cancellation floor
        (2 c^2 eps_mach ~ 8e-12 Ha) never enters a subtraction.  The dense
        eigensolver still carries an absolute backward error ~ ||H||*eps
        ~ 2c^2*eps, so each gap eigenvalue is refined by a Rayleigh
        quotient of its eigenvector evaluated with exactly-rounded
        summation (math.fsum).  This works because in the shifted pencil
        the only O(c^2) block, -2c^2 S in the lower-lower corner, is
        contracted against the SMALL components of an electronic state:
        every term of the quadratic form is O(|E_binding|) and fsum
        removes the summation cancellation exactly.  Net accuracy
        (adversarially verified): ~1e-14 RELATIVE to the binding energy,
        i.e. 50-100 eps -- the floor is set by the pre-rounding of the
        individual products x_i*Hs_ij*x_j and by assembly rounding of the
        O(c^2) entries, which no refinement of the assembled pencil can
        remove.  Caveats: the second-order eigenvector term
        (~(n eps 2c^2/gap)^2) is negligible for level gaps >~ 1e-4 Ha
        (verified) but untested below; for RKB the ill-conditioned
        small-component Gram (cond(Sss) ~ 1e6-1e7) inflates the
        eigensolver floor to ~1e-9 absolute.
        """
        Hs = 0.5 * (Hs + Hs.T)
        S = 0.5 * (S + S.T)
        ev, U = eigh(Hs, S)
        sel = np.nonzero((ev > -2 * self.c ** 2 + 1e-6) & (ev < -1e-9))[0]
        return np.sort([self._rayleigh(Hs, S, U[:, k]) for k in sel])

    def solve_dkb(self, s: int = -1) -> np.ndarray:
        """DKB, Shabaev sign convention (PRL 93, 130405):

            Psi+_j = (  f_j              , +(1/2c) O- f_j )
            Psi-_j = ( s*(1/2c) O+ f_j   ,  f_j           ),   s = -1.

        The MINUS sign is load-bearing: it flips one component only, which
        changes the span (not a basis-sign congruence). It follows from the
        positronic limit P = c O+ Q / (E - c^2 - V) -> -(1/2c) O+ Q at
        E ~ -c^2, and it is what gives the H-- block the correct
        charge-conjugate nonrelativistic limit (kinetic +(1/2)<O+f|O+f>).
        With s=+1 the positronic limit has the wrong-sign kinetic term and
        spurious states appear deep inside the gap (observed empirically;
        the anti-balanced mechanism is proven for bounded V in FINDINGS
        #12). The s parameter exists ONLY for pollution experiments;
        physical calculations must use the default s = -1.
        """
        assert s in (-1, +1)
        c, V = self.c, self.V
        f, g, h = self.f, self.g, self.h
        ME, MEv = self._me, self._mev
        Sff = ME(f, f)
        Spp = Sff + ME(g, g)
        Spm = s * ME(f, h) + ME(g, f)
        Smm = ME(h, h) + Sff
        Hpp = (MEv(f, V, f) + c * c * Sff
               + 4 * c * c * ME(g, g)                   # <f|cO+|g> + <g|cO-|f>
               + MEv(g, V, g) - c * c * ME(g, g))
        Hmm = (MEv(h, V, h) + c * c * ME(h, h)
               + s * 4 * c * c * ME(h, h)               # <s h|cO+|f>+<f|cO-|s h>
               + MEv(f, V, f) - c * c * Sff)
        Hpm = (s * (MEv(f, V, h) + c * c * ME(f, h))    # <f|(V+c^2)|s h>
               + 2 * c * c * ME(f, h)                   # <f|cO+|f> = 2c^2<f|h>
               + s * 0.5 * ME(g, self.OmOp)             # <g|cO-|s h>
               + MEv(g, V, f) - c * c * ME(g, f))
        H = np.block([[Hpp, Hpm], [Hpm.T, Hmm]])
        S = np.block([[Spp, Spm], [Spm.T, Smm]])
        return self._ingap(H, S)

    def solve_rkb(self, shifted: bool = False) -> np.ndarray:
        """RKB: large in {f_j}, small in {O- f_j} with its Gram metric.

        shifted=True: analytic H - c^2 S pencil + Rayleigh refinement
        (see _ingap_shifted)."""
        c, V = self.c, self.V
        f = self.f
        chi = 2 * c * self.g                          # O- f
        ME, MEv = self._me, self._mev
        Sll = ME(f, f)
        Sss = ME(chi, chi)
        # <f_i | c O+ | chi_j> = c <O- f_i | chi_j> = c (Sss)_ij  (by parts)
        S = np.block([[Sll, 0 * Sll], [0 * Sll, Sss]])
        if shifted:
            Hs = np.block([[MEv(f, V, f), c * Sss],
                           [c * Sss.T, MEv(chi, V, chi) - 2 * c * c * Sss]])
            return self._ingap_shifted(Hs, S)
        Hll = MEv(f, V, f) + c * c * Sll
        Hss = MEv(chi, V, chi) - c * c * Sss
        H = np.block([[Hll, c * Sss], [c * Sss.T, Hss]])
        return self._ingap(H, S)

    def solve_naive(self, shifted: bool = False) -> np.ndarray:
        """Both components in the SAME Laguerre set (pollution showcase).

        shifted=True assembles the pencil H - c^2 S analytically (the mass
        term survives only as -2c^2 in the lower-lower block) and refines
        every gap eigenvalue with the exactly-summed Rayleigh quotient;
        see _ingap_shifted."""
        c, V = self.c, self.V
        f, df = self.f, self.df
        kappa, r = self.kappa, self.r
        ME, MEv = self._me, self._mev
        Opf = -df + kappa / r * f
        Hvv = MEv(f, V, f)
        Hoff = c * ME(f, Opf)
        Sff = ME(f, f)
        S = np.block([[Sff, 0 * Sff], [0 * Sff, Sff]])
        if shifted:
            Hs = np.block([[Hvv, Hoff],
                           [Hoff.T, Hvv - 2 * c * c * Sff]])
            return self._ingap_shifted(Hs, S)
        H = np.block([[Hvv + c * c * Sff, Hoff],
                      [Hoff.T, Hvv - c * c * Sff]])
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
