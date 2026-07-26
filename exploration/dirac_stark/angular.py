"""Relativistic dipole angular coefficients for the Dirac Stark problem.

The field term F z = F r cos(theta) multiplies all four spinor
components.  Between Dirac spin-angular functions Omega_{kappa,m} the
angular factor is <Omega_{kappa',m}|cos theta|Omega_{kappa,m}> (rank-1
tensor C^1_0).  The LARGE component of channel kappa carries
Omega_{kappa}; the SMALL component carries Omega_{-kappa}.  So the
coupled radial blocks are

  large-large:  A_L(kappa',kappa) * <g_{kappa'}| r |g_{kappa}>
  small-small:  A_S(kappa',kappa) * <f_{kappa'}| r |f_{kappa}>

with A_L using (kappa',kappa) and A_S using (-kappa',-kappa).

Reduced matrix element (Grant / Johnson convention):
  <kappa'||C^k||kappa> = (-1)^{j'+1/2} sqrt((2j'+1)(2j+1))
                         * ( j' k j ; 1/2 0 -1/2 ) * Pi(l',k,l)
  Pi = 1 if l'+k+l even else 0.
Full:
  <kappa' m|C^1_0|kappa m> = (-1)^{j'-m} ( j' 1 j ; -m 0 m )
                             * <kappa'||C^1||kappa>
"""

from functools import lru_cache

from sympy import Rational, sqrt as ssqrt
from sympy.physics.wigner import wigner_3j


def kappa_to_jl(kappa):
    """Return (j, l, l_bar) as sympy Rationals/ints. j=|kappa|-1/2;
    l is the LARGE-component orbital, l_bar = l(-kappa) the small one."""
    j = Rational(2 * abs(kappa) - 1, 2)
    l = kappa if kappa > 0 else -kappa - 1
    lb = -kappa if kappa < 0 else kappa - 1     # l(-kappa)
    return j, l, lb


@lru_cache(maxsize=None)
def _reduced(kappa_p, kappa, use_bar):
    """<kappa'||C^1||kappa>, using large (use_bar=False) or small
    (use_bar=True) orbital parity."""
    jp, lp, lpb = kappa_to_jl(kappa_p)
    j, l, lb = kappa_to_jl(kappa)
    L, Lp = (lb, lpb) if use_bar else (l, lp)
    if (Lp + 1 + L) % 2 != 0:                   # parity Pi(l',1,l)
        return Rational(0)
    sign = (-1) ** int(jp + Rational(1, 2))
    w = wigner_3j(jp, 1, j, Rational(1, 2), 0, Rational(-1, 2))
    return sign * ssqrt((2 * jp + 1) * (2 * j + 1)) * w


@lru_cache(maxsize=None)
def dipole_coeff(kappa_p, kappa, m2, use_bar):
    """<Omega_{kappa' m}|cos theta|Omega_{kappa m}>, m = m2/2 (m2 odd int).
    use_bar=False for large-large, True for small-small."""
    jp, _, _ = kappa_to_jl(kappa_p)
    j, _, _ = kappa_to_jl(kappa)
    m = Rational(m2, 2)
    if abs(m) > jp or abs(m) > j:
        return Rational(0)
    red = _reduced(kappa_p, kappa, use_bar)
    if red == 0:
        return Rational(0)
    sign = (-1) ** int(jp - m)
    w = wigner_3j(jp, 1, j, -m, 0, m)
    return sign * w * red


def A_L(kappa_p, kappa, m2):
    return float(dipole_coeff(kappa_p, kappa, m2, use_bar=False))


def A_S(kappa_p, kappa, m2):
    return float(dipole_coeff(kappa_p, kappa, m2, use_bar=True))


def kappa_ladder(kappa0, nchan, m2):
    """Build the coupled kappa-ladder of length nchan starting from
    kappa0, following nonzero dipole links (alternating parity).  For
    m=1/2 (m2=1) starting at s1/2 (kappa0=-1):
    -1, +1, -2, +2, -3, ... (each l appears as both j-partners)."""
    j0 = abs(kappa0) - 0.5
    jmin = abs(m2) / 2.0
    ladder = []
    # order by increasing |kappa| (increasing l), both signs, j>=|m|
    k = 1
    while len(ladder) < nchan:
        for kap in (-k, k):
            j = abs(kap) - 0.5
            if j >= jmin - 1e-9:
                ladder.append(kap)
        k += 1
    return ladder[:nchan]
