"""Paper F, F0b: assumption-free cross-check of the dipole angular
coefficients via Clebsch-Gordan decomposition of the spin-angular
functions (a completely different route from the Wigner-3j reduced
matrix element in angular.py).

Omega_{kappa m} = sum_{ms=+-1/2} CG(l, m-ms, 1/2, ms | j, m) Y_{l,m-ms} chi_ms
<Omega'|cos|Omega> = sum_ms CG'(...) CG(...) G(l', l, m-ms)
with scalar Gaunt G(l',l,mu) = <Y_{l' mu}|cos theta|Y_{l mu}>:
  G(l+1,l,mu) = sqrt(((l+1)^2 - mu^2)/((2l+1)(2l+3)))
  G(l-1,l,mu) = sqrt((l^2 - mu^2)/((2l-1)(2l+1)))
This is textbook and shares no code path with the wigner_3j formula.
"""

import sys
from pathlib import Path

from sympy import Rational, sqrt as ssqrt
from sympy.physics.quantum.cg import CG

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from angular import A_L, kappa_to_jl, kappa_ladder  # noqa: E402


def gaunt_cos(lp, l, mu):
    """<Y_{lp,mu}|cos theta|Y_{l,mu}>, mu integer/half not -- here mu int."""
    if lp == l + 1:
        num = (l + 1) ** 2 - mu ** 2
        den = (2 * l + 1) * (2 * l + 3)
    elif lp == l - 1:
        num = l ** 2 - mu ** 2
        den = (2 * l - 1) * (2 * l + 1)
    else:
        return Rational(0)
    if num < 0:
        return Rational(0)
    return ssqrt(Rational(num, den))


def coeff_cg(kappa_p, kappa, m2):
    """<Omega_{kappa' m}|cos|Omega_{kappa m}> via CG decomposition."""
    jp, lp, _ = kappa_to_jl(kappa_p)
    j, l, _ = kappa_to_jl(kappa)
    m = Rational(m2, 2)
    tot = Rational(0)
    for ms2 in (1, -1):
        ms = Rational(ms2, 2)
        ml = m - ms                      # integer
        if abs(ml) > l or abs(ml) > lp:
            continue
        cg_p = CG(lp, ml, Rational(1, 2), ms, jp, m).doit()
        cg = CG(l, ml, Rational(1, 2), ms, j, m).doit()
        tot += cg_p * cg * gaunt_cos(lp, l, int(ml))
    return tot


print("F0b CG cross-check of A_L (large-component) vs wigner_3j route:\n")
m2 = 1
lad = kappa_ladder(-1, 8, m2)
maxdiff = 0.0
for kp in lad:
    for k in lad:
        a_wig = A_L(kp, k, m2)
        a_cg = float(coeff_cg(kp, k, m2))
        d = abs(a_wig - a_cg)
        maxdiff = max(maxdiff, d)
        if d > 1e-10 and (abs(a_wig) > 1e-10 or abs(a_cg) > 1e-10):
            print(f"   MISMATCH kappa'={kp:+d} kappa={k:+d}: "
                  f"wig={a_wig:+.6f}  cg={a_cg:+.6f}  d={d:.1e}")
print(f"max |wigner - CG| over all pairs: {maxdiff:.2e}  "
      f"({'OK -- formulas agree' if maxdiff < 1e-10 else 'FAIL'})")

# also print the full large-component coupling matrix for the record
print("\nA_L matrix (m=1/2), rows/cols = ladder", lad, ":")
for kp in lad:
    print("  " + "  ".join(f"{A_L(kp,k,m2):+.4f}" for k in lad))
