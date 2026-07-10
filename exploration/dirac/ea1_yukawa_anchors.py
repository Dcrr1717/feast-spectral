"""E-A1 (final form): external validation and adjudication of the
published Dirac-Yukawa anchor points.

RESOLVED CONVENTIONS (deduced from residual patterns, then independently
verified against Baye Table VI and the Krauthauser-Hill 40-digit
benchmarks; see FINDINGS.md section 7):
- c = 137.0359895 (CODATA 1986) underlies ALL published values in the
  comparison chain (never stated in any of the three papers).
- Egs = E - c^2 in Hartree; V(r) = -V0 e^{-lam r}/r.
- CSW Table 4 "Case 2: V0=0.4c, lam=0.07c" is an ERRATUM: the energy
  corresponds to Baye's case V0=0.7c, lam=0.04c (digits transposed,
  consistently in prose, table and figure caption).
- Underlying c=1 benchmarks (Krauthauser-Hill 2002, 40-digit, via Baye
  Table VI): f1 = -0.00408291802884811, f2 = -0.25879891617626010.

RESULTS (this script reproduces them):
  Case 1: our f agrees with KH to ~1e-14; ALLSM (2026) off by 4.4e-11 Ha.
  Case 2: our f agrees with KH to  2e-15; Baye LM good to 3.7e-11 Ha;
          ALLSM off by 7.6e-9 Ha.  Both adjudications favor KH.
  The Legendre quadrature stalls at 1e-4 for gamma=0.714 while the exact
  cusp rule is machine-precise: live proof of the exact-quadrature claim.
"""

import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from dirac import RadialDiracBasis  # noqa: E402

t0 = time.perf_counter()
C86 = 137.0359895                      # CODATA 1986, the chain's convention

CASES = [
    ("Case 1 (V0=0.1c, lam=0.01c)", 0.1, 0.01, 8.0,
     -0.00408291802884811,            # Krauthauser-Hill 40-digit (c=1)
     -76.672555928755170,             # ALLSM 2026
     -76.672555928706560),            # reference row (KH via CSW Table 4)
    ("Case 2 (V0=0.7c, lam=0.04c) [CSW erratum: printed as 0.4c/0.07c]",
     0.7, 0.04, 92.0,
     -0.25879891617626010,
     -4859.949240857070,
     -4859.949240864677),             # Baye Lagrange-mesh value
]

for name, nu0, mu, beta, f_kh, e_allsm, e_ref in CASES:
    gam = np.sqrt(1.0 - nu0 ** 2)
    Z, lam = nu0 * C86, mu * C86
    print(f"== {name}  (gamma={gam:.6f}) ==")
    for quad in ("cusp", "legendre"):
        vals = []
        for nb in (60, 80):
            nq = (nb + 80) if quad == "cusp" else None
            B = RadialDiracBasis(nb, -1, Z, beta, cert_tol=np.inf, mu=gam,
                                 quad=quad, nquad=nq, c=C86)
            B.V = -Z * np.exp(-lam * B.r) / B.r
            vals.append(B.solve_naive()[0])
        print(f"   {quad:>9}: E(nb=80) = {vals[1]:.12f}"
              f"   drift(60->80) = {abs(vals[1]-vals[0]):.0e}")
        if quad == "cusp":
            E, f = vals[1], vals[1] / C86 ** 2
    print(f"   f = {f:.17f}")
    print(f"   f - KH(40d) = {f - f_kh:+.2e}   E - ALLSM = {E - e_allsm:+.2e}"
          f"   E - ref = {E - e_ref:+.2e}")
    print(f"   published ALLSM-ref discrepancy: {e_allsm - e_ref:+.1e}\n")

print(f"total: {time.perf_counter()-t0:.1f} s")
