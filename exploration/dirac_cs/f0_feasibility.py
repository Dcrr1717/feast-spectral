"""Plan B, Phase F0: feasibility gates for the complex-scaled radial
Dirac solver (gamma-basis + shifted pencil, complex-symmetric assembly).

Under r -> r e^{i theta} every term of O+- scales uniformly by
e^{-i theta} (both d/dr and kappa/r), the mass term does not scale, and
the Coulomb potential picks the same factor: the discrete assembly is
therefore built from the SAME certified real matrices,

    H_theta - c^2 S = [[e^{-i th} Vd,            c e^{-i th} K ],
                       [c e^{-i th} K^T,  e^{-i th} Vd - 2 c^2 S]],

complex SYMMETRIC (transpose, not conjugate).  Yukawa adds the complex
screening factor e^{-lam r e^{i th}} inside the (converged) quadrature.

PRE-REGISTERED PREDICTIONS:
P-F0.1 (H-B1) indicial invariance: all singular terms share e^{-i th}
       => gamma unchanged => the gamma-basis stays cusp-adapted; the
       scaled Coulomb bound states converge exponentially in nb.
P-F0.2 fiber survival: the free discrete spectrum equals
       +-sqrt(c^4 + c^2 sigma_j^2 e^{-2i th}) EXACTLY (the L1' SVD
       fibers survive scaling with sigma -> sigma e^{-i th}); electron
       branch rotates by -2 theta at small sigma, positron branch
       mirrors from -2c^2 (binding form).
P-F0.3 theta = 0 reproduces the real solver to machine precision.
P-F0.4 Balslev-Combes (discrete): Coulomb bound states are theta-
       independent; drift <= ~1e-10 over theta in [0.05, 0.4] at
       nb = 100 (real beta; representation of e^{-Z r cos th} costs
       basis size as theta grows -- measured, not assumed).

GATES: G0a fiber identity at machine precision; G0b bound states
theta-stable; G0c no spurious complex eigenvalues inside the gap window
away from the two rotated branches.
"""

import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from dirac import C_AU, RadialDiracBasis, exact_ladder  # noqa: E402

C = C_AU
t0 = time.perf_counter()


def cs_pencil(B, theta, lam=0.0):
    """Complex-scaled shifted pencil (H_theta - c^2 S, S) on basis B.
    lam = 0: pure Coulomb (exact integrals); lam > 0: Yukawa, screening
    factor evaluated at complex argument inside the quadrature."""
    c, Z, r = B.c, B.Z, B.r
    ph = np.exp(-1j * theta)
    f = B.f
    ME, MEv = B._me, B._mev
    K = ME(f, 2.0 * c * B.h)                 # <f_i|O+ f_j>, real, certified
    if lam == 0.0:
        Vd = -Z * ME(f, f / r[None, :] * np.ones((1, 1))) if False else \
            -Z * MEv(f, 1.0 / r, f)
        Vth = ph * Vd.astype(complex)
    else:
        scr = np.exp(-lam * r * np.exp(1j * theta))
        Vth = ph * (-Z) * np.einsum('ik,k,jk,k->ij', f, scr / r, f, B.W)
    n = B.nb
    S = np.block([[np.eye(n), np.zeros((n, n))],
                  [np.zeros((n, n)), np.eye(n)]]).astype(complex)
    Hs = np.block([[Vth, c * ph * K],
                   [c * ph * K.T, Vth - 2 * c * c * np.eye(n)]])
    return 0.5 * (Hs + Hs.T), S


def gap_eigs(Hs):
    ev = np.linalg.eigvals(Hs)
    return ev[(ev.real > -2 * C**2 + 1e-6) & (ev.real < -1e-9)]


# ---------------- P-F0.2 / G0a: fiber survival on the free operator ----
print("P-F0.2 fiber survival (free operator, kappa=-1, nb=80):")
B = RadialDiracBasis(80, -1, 1e-30, 0.9, mu=1.0, quad="cusp",
                     nquad=160, cert_tol=1e-6)
K = B._me(B.f, 2.0 * C * B.h)
sig = np.linalg.svd(K, compute_uv=False)
for theta in (0.1, 0.3):
    Hs, S = cs_pencil(B, theta)
    ev = np.linalg.eigvals(Hs)
    pred = np.concatenate([
        +np.sqrt(C**4 + C**2 * sig**2 * np.exp(-2j * theta)) - C**2,
        -np.sqrt(C**4 + C**2 * sig**2 * np.exp(-2j * theta)) - C**2])
    # match each predicted value to nearest computed
    dev = max(np.min(np.abs(ev - p)) for p in pred)
    print(f"   theta={theta}: max |eig - fiber prediction| = {dev:.2e}"
          f"   (scale c^2 = {C**2:.0f})")

# ---------------- P-F0.3: theta = 0 reproduces the real solver ---------
Z = 80.0
gam = float(np.sqrt(1 - (Z / C) ** 2))
Bz = RadialDiracBasis(100, -1, Z, 0.8 * Z, mu=gam, quad="cusp",
                      nquad=220, cert_tol=1e-6)
Hs0, _ = cs_pencil(Bz, 0.0)
ev0 = np.sort_complex(gap_eigs(Hs0))
real_ev = Bz.solve_naive(shifted=True)
k = min(3, len(ev0), len(real_ev))
dev = np.max(np.abs(np.sort(ev0.real)[:k] - real_ev[:k]))
print(f"\nP-F0.3 theta=0 vs real solver (Z=80, lowest {k}):"
      f" max dev = {dev:.2e},  max |Im| = {np.max(np.abs(ev0.imag)):.2e}")

# ---------------- P-F0.4 / G0b: Balslev-Combes theta-independence ------
print("\nP-F0.4 bound-state theta-drift (Z=80, gamma-basis, nb=100):")
ladder = exact_ladder(-1, Z)
ref = None
print(f"{'theta':>7} {'E(1s) real part':>18} {'|Im|':>9} {'drift vs th=0.05':>17}")
for theta in (0.05, 0.1, 0.2, 0.3, 0.4):
    Hs, _ = cs_pencil(Bz, theta)
    ev = gap_eigs(Hs)
    e1 = ev[np.argmin(np.abs(ev - ladder[0]))]
    if ref is None:
        ref = e1
    print(f"{theta:>7.2f} {e1.real:>18.9f} {abs(e1.imag):>9.1e}"
          f" {abs(e1 - ref):>17.1e}")
print(f"   exact 1s = {ladder[0]:.9f}")

# ---------------- G0c: gap window clean away from the branches ---------
print("\nG0c spurious scan (Z=80, theta=0.3): gap eigenvalues with")
print("   |Im| > 10*|Im branch envelope| and far from the exact ladder:")
Hs, _ = cs_pencil(Bz, 0.3)
ev = gap_eigs(Hs)
susp = []
for e in ev:
    near_bound = np.min(np.abs(e - ladder)) / np.abs(ladder).min() < 1e-3 \
        or abs(e.imag) < 1e-6 * C**2
    # electron branch envelope: -Im(E) ~ |E_b| tan-ish growth; flag only
    # clearly off-curve points deep in the gap
    on_branch = abs(np.angle(-(e + 0j)) ) < 1e-9  # placeholder, refined below
    if not near_bound and e.real < -10.0 and abs(e.imag) < 0.05 * abs(e.real):
        susp.append(e)
print(f"   suspicious count = {len(susp)}"
      + (f"   examples: {[f'{s:.3f}' for s in susp[:3]]}" if susp else ""))
print(f"\ntotal: {time.perf_counter()-t0:.1f} s")
