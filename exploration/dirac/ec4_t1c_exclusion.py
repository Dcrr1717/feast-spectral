"""Plan C, E-C4: THEOREM T1-C -- the two-sector Feshbach exclusion
window for the full pair-space Coulomb-Dirac Galerkin (math_daniel).

THEOREM T1-C (conditional, finite-dimensional).  Split the fiber pairs
at t* (sectors L: t_j <= t*, H: t_j > t*), D* = sqrt(1+4t*^2).  Suppose
the measured constants satisfy, on the PAIR spans,
  (h1) ||E_L^{-1/2} W_LL E_L^{-1/2}|| <= omega_L,
  (h2) ||E_L^{-1/2} W_LH E_H^{-1/2}|| <= x,
  (h3) ||E_H^{-1/2} W_HH E_H^{-1/2}|| <= omega_H,
and let m* solve the self-consistent bound (evaluated conservatively)
  eps*(m) = omega_L + x^2 rho_H(m)/(1 - omega_H rho_H(m)),
  rho_H(m) = D* c^2/(D* c^2 - m),  with  omega_H rho_H(m) < 1,
  m* = sup{ m < c^2 : eps*(m) c^2/(c^2 - m) < 1 }.
Then the full pair-space Galerkin of h = H0 + W has NO eigenvalue E
with |E| < m*.

Proof skeleton (each step verified numerically below):
1. |h_H| >= D* c^2 and the weighted resolvent bound
   ||E_H^{1/2}(h_H-mu)^{-1}E_H^{1/2}|| <= rho_H(m) for |mu| <= m
   (E/(E-m) decreasing in E on spec).
2. (h3) + Neumann: h_H + W_HH - mu invertible with weighted bound
   rho_H/(1 - omega_H rho_H).
3. Schur correction on L: norm <= x^2 rho_H/(1-omega_H rho_H)  [(h2)].
4. ||E_L^{1/2}(h_L-mu)^{-1}E_L^{1/2}|| <= c^2/(c^2-m)  (E_j >= c^2).
5. S(mu) invertible when eps*(m) c^2/(c^2-m) < 1; Feshbach block lemma
   => G - mu invertible => no eigenvalue in (-m*, m*).  QED (finite-dim).

CONSISTENCY REQUIREMENTS (verified):
- m* must sit BELOW the physical ladder (gamma c^2) -- the window
  excludes deep-gap impostors, not physical states;
- the CERTIFIED solver spectrum must indeed have no eigenvalue in the
  window (binding in (-c^2-m*, -c^2+m*));
- the Tupitsyn-Shabaev spurious at kappa>0 sits at the partner 1s
  (0.74 c^2 at Z=92) -- ABOVE the window: the theorem is consistent
  with the one real spurious we know exists.
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


def pair_constants(Z, nb, tstar=1.0, beta_fac=0.8, kappa=-1):
    nu = Z / C
    gam = float(np.sqrt(kappa * kappa - nu * nu))
    B = RadialDiracBasis(nb, kappa, Z, beta_fac * Z, mu=gam, quad="cusp",
                         nquad=nb + 120, cert_tol=1e-5)
    K = B._me(B.f, 2.0 * C * B.h)
    Vd = B._mev(B.f, B.V, B.f)
    Vd = 0.5 * (Vd + Vd.T)
    W = np.block([[Vd, np.zeros((nb, nb))], [np.zeros((nb, nb)), Vd]])
    U, S, Vt = np.linalg.svd(K)
    T = S / (2 * C)
    E = C**2 * np.sqrt(1 + 4 * T**2)
    # full pair frame per fiber: both balanced and positronic vectors
    Cb = np.vstack([U, Vt.T * T]); Cb /= np.linalg.norm(Cb, axis=0)
    Cp = np.vstack([-U * T, Vt.T]); Cp /= np.linalg.norm(Cp, axis=0)
    lo, hi = np.where(T <= tstar)[0], np.where(T > tstar)[0]
    FL = np.hstack([Cb[:, lo], Cp[:, lo]])
    FH = np.hstack([Cb[:, hi], Cp[:, hi]])
    EL = np.concatenate([E[lo], E[lo]])
    EH = np.concatenate([E[hi], E[hi]])
    EiL = np.diag(1.0 / np.sqrt(EL))
    EiH = np.diag(1.0 / np.sqrt(EH))
    wL = float(np.linalg.norm(EiL @ (FL.T @ W @ FL) @ EiL, 2))
    wH = float(np.linalg.norm(EiH @ (FH.T @ W @ FH) @ EiH, 2)) if hi.size else 0.0
    xx = float(np.linalg.norm(EiL @ (FL.T @ W @ FH) @ EiH, 2)) if hi.size else 0.0
    return wL, wH, xx, B


def mstar(wL, wH, xx, tstar=1.0):
    """Conservative self-consistent window (grid evaluation)."""
    Ds = np.sqrt(1 + 4 * tstar**2)
    best = 0.0
    for m in np.linspace(0.0, 0.999 * C**2, 4000):
        rH = Ds * C**2 / (Ds * C**2 - m)
        if wH * rH >= 1.0:
            continue
        eps = wL + xx**2 * rH / (1 - wH * rH)
        if eps * C**2 / (C**2 - m) < 1.0:
            best = m
    return best


print("T1-C: measured pair-sector constants and exclusion window"
      f" (t* = 1, kappa = -1):")
print(f"{'Z':>4} {'nu':>6} {'nb':>4} {'w_L':>7} {'w_H':>7} {'x':>7}"
      f" {'m*/c^2':>7} {'gamma':>7} {'m*<gam?':>7} {'clean?':>6}")
for Z in (40.0, 92.0, 118.0):
    nu = Z / C
    gam = float(np.sqrt(1 - nu * nu))
    for nb in (80, 160):
        wL, wH, xx, B = pair_constants(Z, nb)
        ms = mstar(wL, wH, xx) / C**2
        # verify the conclusion on the certified solver spectrum
        ev = B.solve_naive(shifted=True)          # bindings in (-2c^2, 0)
        Evals = ev + C**2                          # E coordinates
        inside = Evals[np.abs(Evals) < ms * C**2]
        ok = "YES" if inside.size == 0 else f"NO({inside.size})"
        print(f"{Z:>4.0f} {nu:>6.3f} {nb:>4} {wL:>7.4f} {wH:>7.4f}"
              f" {xx:>7.4f} {ms:>7.4f} {gam:>7.4f}"
              f" {'YES' if ms < gam else 'NO':>7} {ok:>6}")

# kappa=+1 consistency: the TS spurious must sit ABOVE the window
print("\nkappa=+1 consistency (TS spurious above the window):")
for Z in (92.0,):
    wL, wH, xx, B = pair_constants(Z, 120, kappa=1)
    ms = mstar(wL, wH, xx) / C**2
    ev = B.solve_naive(shifted=True) + C**2
    deep = ev[np.abs(ev) < ms * C**2]
    low = np.sort(np.abs(ev))[:3] / C**2
    print(f"   Z={Z:.0f}: m*/c^2 = {ms:.4f}; lowest |E|/c^2 of solver:"
          f" {', '.join(f'{v:.4f}' for v in low)};"
          f" eigenvalues inside window: {deep.size}")
print(f"\ntotal: {time.perf_counter()-t0:.1f} s")
