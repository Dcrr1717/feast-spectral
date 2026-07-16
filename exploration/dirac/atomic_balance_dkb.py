"""Atomic (Coulomb-regularized) dual kinetic balance at a point nucleus.

The free DKB dual family Psi- = (-(1/2c) O+ f, f) violates the Coulomb
form domain for kappa = -1 at a point nucleus: O+ f does not vanish at
r = 0, so <O+ f|Z/r|O+ f> is log-divergent (FINDINGS #3a, verified by
the log-linear quadrature drift).  The atomic-balance remedy replaces
the free factor 1/2c by the Coulomb-regularized

    w(r) = c / (2c^2 + Z/r),    w -> c r/Z  (r -> 0),
                                w -> 1/2c   (r -> infinity),

which restores the form domain: w * O+ f ~ r at the origin.  For the
electronic direction this is the standard atomic/ZORA-type balance
(denominator 2c^2 - V never vanishes for attractive V); we use the SAME
regular factor for the dual direction -- a legitimate trial-space
choice whose quality is measured, not assumed.

Tests (pre-registered):
A1  Form domain restored: <w O+ f | Z/r | w O+ f> converges with nq
    (vs the log-linear divergence of the free counterpart).
A2  Spectra: kappa = -1, Z in {80, 120}, INTEGER-prefactor basis
    (the setting where free DKB fails): atomic-DKB gap spectrum is
    clean, variational, and quadrature-stable.
A3  Comparison at fixed nq: free-DKB eigenvalue drifts with nq
    (reproducing #3a) while atomic-DKB is stationary.
"""

import sys
import time
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from dirac import C_AU, RadialDiracBasis, exact_ladder  # noqa: E402

C = C_AU
t0 = time.perf_counter()


def atomic_dkb(B, Z):
    """Assemble the atomic-balance DKB pencil on basis B (any quad)."""
    c, r, V = B.c, B.r, B.V
    f, df = B.f, B.df
    kappa = B.kappa
    ME, MEv = B._me, B._mev
    Omf = 2.0 * c * B.g                      # O- f
    Opf = 2.0 * c * B.h                      # O+ f
    w = c / (2 * c * c + Z / r)
    wp = c * (Z / r**2) / (2 * c * c + Z / r) ** 2
    ga = w * Omf                             # atomic small comp of Psi+
    ha = w * Opf                             # atomic large comp of Psi-
    # operator applications on the atomic components:
    # O+ (w chi) = w O+ chi - w' chi ;  O- (w chi) = w O- chi + w' chi
    Opga = w * B.OpOm + wp * 0 - wp * Omf    # O+(w O- f) = w OpOm - w' O- f
    Omha = w * B.OmOp + wp * Opf             # O-(w O+ f) = w OmOp + w' O+ f
    # Psi+_j = (f_j, ga_j),  Psi-_j = (-ha_j, f_j)
    Spp = ME(f, f) + ME(ga, ga)
    Spm = -ME(f, ha) + ME(ga, f)
    Smm = ME(ha, ha) + ME(f, f)
    Hpp = (MEv(f, V, f) + c * c * ME(f, f) + c * ME(f, Opga)
           + c * ME(ga, Omf) + MEv(ga, V, ga) - c * c * ME(ga, ga))
    Hmm = (MEv(ha, V, ha) + c * c * ME(ha, ha) - c * ME(ha, Opf)
           - c * ME(f, Omha) + MEv(f, V, f) - c * c * ME(f, f))
    Hpm = (-(MEv(f, V, ha) + c * c * ME(f, ha)) + c * ME(f, Opf)
           - c * ME(ga, Omha) + MEv(ga, V, f) - c * c * ME(ga, f))
    H = np.block([[Hpp, Hpm], [Hpm.T, Hmm]])
    S = np.block([[Spp, Spm], [Spm.T, Smm]])
    ev = eigh(0.5 * (H + H.T), 0.5 * (S + S.T), eigvals_only=True)
    ev = np.sort(ev) - c ** 2
    return ev[(ev > -2 * c ** 2 + 1e-6) & (ev < -1e-9)]


# ---------------- A1: form-domain restoration ---------------------------
print("A1 form-domain: <chi|Z/r|chi> vs nq (Z=80, kappa=-1, nb=40,"
      " integer basis):")
Z = 80.0
print(f"{'nq':>7} {'free: <O+f|V|O+f>':>19} {'atomic: <wO+f|V|wO+f>':>22}")
prev_free, prev_at = None, None
for nq in (2000, 4000, 8000):
    B = RadialDiracBasis(40, -1, Z, 0.9 * Z, nquad=nq, cert_tol=1e-4)
    Opf = 2.0 * C * B.h
    w = C / (2 * C * C + Z / B.r)
    e_free = float(np.sum(B.W * B.V * Opf[0] * Opf[0]))
    e_at = float(np.sum(B.W * B.V * (w * Opf[0]) ** 2))
    d1 = "" if prev_free is None else f"  (drift {e_free - prev_free:+.1f})"
    d2 = "" if prev_at is None else f"  (drift {e_at - prev_at:+.2e})"
    print(f"{nq:>7} {e_free:>19.1f}{d1:<18} {e_at:>18.6f}{d2}")
    prev_free, prev_at = e_free, e_at

# ---------------- A2: spectra ------------------------------------------
print("\nA2 atomic-DKB spectra, integer basis, nb=60:")
for Z in (80.0, 120.0):
    ladder = exact_ladder(-1, Z)
    B = RadialDiracBasis(60, -1, Z, 0.9 * Z, nquad=6000, cert_tol=1e-4)
    ev = atomic_dkb(B, Z)
    k = min(3, ev.size)
    rel = (ev[:k] - ladder[:k]) / np.abs(ladder[:k])
    viol = bool(np.any(rel < -1e-8))
    print(f"   Z={Z:>4.0f}: lowest3 rel err = "
          + " ".join(f"{x:+.1e}" for x in rel)
          + f"   variational-violation: {viol}   n_ingap={ev.size}")

# ---------------- A3: quadrature stability ------------------------------
print("\nA3 ground-state stability vs nq (Z=80, nb=60):")
print(f"{'nq':>7} {'free-DKB E0':>16} {'atomic-DKB E0':>16}")
for nq in (2000, 4000, 8000):
    B = RadialDiracBasis(60, -1, 80.0, 72.0, nquad=nq, cert_tol=1e-4)
    e_free = B.solve_dkb()[0]
    e_at = atomic_dkb(B, 80.0)[0]
    print(f"{nq:>7} {e_free:>16.6f} {e_at:>16.6f}")
print(f"\ntotal: {time.perf_counter()-t0:.1f} s")
