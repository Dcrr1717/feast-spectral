"""Plan C, E-C2: the bridge lemmas and the decisive coupling measurement
(math_daniel protocol; symbols: ../../math-context.md; follows E-C1,
FINDINGS #17).

STRUCTURAL LEMMAS (stated, then verified):
(C2a) Component reduction: H_ref^2 = diag(c^4+c^2 KK^T, c^4+c^2 K^TK)
      (SUSY structure), hence |H_ref|, and therefore A2 = (PH^2P)^{1/2}
      and A3 = P|H|P, are COMPONENT-block-diagonal; so is
      A1 = |PHP| (since (PHP)^2 is). The Dirac bridge problem reduces
      to two SCALAR (Schrodinger-like) problems.
(C2b) The E-C1 order chain has classical names:
      A1 = ((PHP)^2)^{1/2} <= (PH^2P)^{1/2} = A2 by Kadison-Schwarz
      ((PHP)^2 <= PH^2P) plus Loewner monotonicity of sqrt;
      A3 = P|H|P <= (PH^2P)^{1/2} = A2 by Hansen's inequality
      (operator concavity of sqrt). Both citable theorems.

DECISIVE MEASUREMENTS:
(C2c) omega_N := || E^{-1/2} W_bal E^{-1/2} ||, the COULOMB FORM ON THE
      BALANCED SPAN relative to the fiber energies E_j -- intrinsic to
      the small space (no reference needed). If omega_N converges, the
      Coulomb coupling is controlled fiber-wise on the balanced family
      and the A3 detour is unnecessary for the form bound.
(C2d) eta_bal := || E^{-1/2} (A3 - A1)_bal E^{-1/2} || on the balanced
      span (needs reference): how much of the E-C1 bridge cost survives
      on the physically relevant span.
(C2e) Schur/locality structure of Delta = A3 - A1 in the fiber basis:
      where do the cross terms live (t_i vs t_j map)?

PRE-REGISTERED PREDICTIONS:
P-C2.1 off-component blocks of A1, A2, A3 vanish to roundoff.
P-C2.2 omega_N converges in nb (Hardy heuristics: limit ~ alpha_3-ish,
       < 1 for nu <= 0.75-ish) and scales ~ nu across Z.
P-C2.3 eta_bal << c_bridge - 1 (the E-C1 cross-fiber cost lives mostly
       OUTSIDE the balanced span).
"""

import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from dirac import C_AU, RadialDiracBasis  # noqa: E402

C = C_AU
t0 = time.perf_counter()


def small_objects(Z, nb, beta_fac=0.8):
    nu = Z / C
    gam = float(np.sqrt(1 - nu * nu))
    B = RadialDiracBasis(nb, -1, Z, beta_fac * Z, mu=gam, quad="cusp",
                         nquad=nb + 120, cert_tol=1e-5)
    K = B._me(B.f, 2.0 * C * B.h)
    Vd = B._mev(B.f, B.V, B.f)
    Vd = 0.5 * (Vd + Vd.T)
    U, S, Vt = np.linalg.svd(K)
    T = S / (2 * C)
    E = C**2 * np.sqrt(1 + 4 * T**2)
    # balanced fiber vectors (columns, orthonormal)
    Cb = np.vstack([U, Vt.T * T])
    Cb = Cb / np.linalg.norm(Cb, axis=0)
    W = np.block([[Vd, np.zeros((nb, nb))], [np.zeros((nb, nb)), Vd]])
    return K, Vd, W, Cb, E, T


def ref_bridge(Z, nb, nb_ref, beta_fac=0.8):
    """A1, A3 on the small space from a nested reference."""
    nu = Z / C
    gam = float(np.sqrt(1 - nu * nu))
    B = RadialDiracBasis(nb_ref, -1, Z, beta_fac * Z, mu=gam, quad="cusp",
                         nquad=nb_ref + 120, cert_tol=1e-5)
    K = B._me(B.f, 2.0 * C * B.h)
    n = nb_ref
    H0 = np.block([[C**2 * np.eye(n), C * K], [C * K.T, -C**2 * np.eye(n)]])
    sel = np.concatenate([np.arange(nb), n + np.arange(nb)])
    PHP = 0.5 * (H0[np.ix_(sel, sel)] + H0[np.ix_(sel, sel)].T)
    ev, Us = np.linalg.eigh(PHP)
    A1 = Us @ np.diag(np.abs(ev)) @ Us.T
    evr, Ur = np.linalg.eigh(H0)
    absH = Ur @ np.diag(np.abs(evr)) @ Ur.T
    A3 = 0.5 * (absH[np.ix_(sel, sel)] + absH[np.ix_(sel, sel)].T)
    H0sq = H0 @ H0
    Msq = 0.5 * (H0sq[np.ix_(sel, sel)] + H0sq[np.ix_(sel, sel)].T)
    e2, U2 = np.linalg.eigh(Msq)
    A2 = U2 @ np.diag(np.sqrt(np.maximum(e2, 0))) @ U2.T
    return A1, A2, A3


# ---------------- P-C2.1: component-diagonality (Lemma C2a) -------------
print("P-C2.1 component-off-block norms (Z=92, nb=80, ref=320):")
A1, A2, A3 = ref_bridge(92.0, 80, 320)
for name, A in (("A1", A1), ("A2", A2), ("A3", A3)):
    off = A[:80, 80:]
    print(f"   {name}: ||upper-lower block|| = {np.linalg.norm(off, 2):.2e}"
          f"   (scale ||A|| = {np.linalg.norm(A, 2):.1e})")

# ---------------- P-C2.2: omega_N, the decisive intrinsic number --------
print("\nP-C2.2 omega_N = ||E^-1/2 W_bal E^-1/2|| (balanced span, Coulomb):")
print(f"{'Z':>4} {'nu':>6} " + " ".join(f"nb={n:<4}" for n in (40, 80, 120, 160))
      + "  omega/nu (last)")
for Z in (40.0, 92.0, 103.0, 118.0):
    nu = Z / C
    row = []
    for nb in (40, 80, 120, 160):
        K, Vd, W, Cb, E, T = small_objects(Z, nb)
        Wb = Cb.T @ W @ Cb
        Ei = np.diag(1.0 / np.sqrt(E))
        row.append(float(np.linalg.norm(Ei @ Wb @ Ei, 2)))
    print(f"{Z:>4.0f} {nu:>6.3f} " + " ".join(f"{w:7.4f}" for w in row)
          + f"   {row[-1] / nu:7.4f}")

# ---------------- P-C2.3: eta_bal vs the raw bridge cost ----------------
print("\nP-C2.3 eta_bal = ||E^-1/2 (A3-A1)_bal E^-1/2|| vs raw c_bridge-1:")
for Z in (40.0, 92.0):
    for nb, ref in ((40, 160), (80, 320)):
        A1, A2, A3 = ref_bridge(Z, nb, ref)
        K, Vd, W, Cb, E, T = small_objects(Z, nb)
        Db = Cb.T @ (A3 - A1) @ Cb
        Ei = np.diag(1.0 / np.sqrt(E))
        eta = float(np.linalg.norm(Ei @ Db @ Ei, 2))
        # raw bridge (E-C1 style): whole-space norm
        ev1, U1 = np.linalg.eigh(A1)
        Mi1 = U1 @ np.diag(1.0 / np.sqrt(np.maximum(ev1, 1e-300))) @ U1.T
        cbr = float(np.linalg.norm(Mi1 @ A3 @ Mi1, 2))
        print(f"   Z={Z:>3.0f} nb={nb:>3} ref={ref:>3}: eta_bal = {eta:.4f}"
              f"   (raw c_bridge - 1 = {cbr - 1:.3f})")

# ---------------- C2e: where do the cross terms live? -------------------
print("\nC2e locality of Delta = A3 - A1 in the fiber-energy scaling")
print("   (Z=92, nb=80, ref=320): largest |Delta_ij|/sqrt(EiEj) entries:")
A1, A2, A3 = ref_bridge(92.0, 80, 320)
K, Vd, W, Cb, E, T = small_objects(92.0, 80)
# full fiber frame: balanced AND positronic vectors span the small space
U, S, Vt = np.linalg.svd(K)
Tt = S / (2 * C)
Cp = np.vstack([-U * Tt, Vt.T])
Cp = Cp / np.linalg.norm(Cp, axis=0)
F = np.hstack([Cb, Cp])                      # 2nb frame (orthonormal-ish)
D = F.T @ (A3 - A1) @ F
EE = np.concatenate([E, E])
Dsc = np.abs(D) / np.sqrt(EE[:, None] * EE[None, :])
idx = np.dstack(np.unravel_index(np.argsort(-Dsc, axis=None)[:5], Dsc.shape))[0]
for i, j in idx:
    ti = Tt[i % 80]
    tj = Tt[j % 80]
    sec = ("bal" if i < 80 else "pos") + "-" + ("bal" if j < 80 else "pos")
    print(f"   ({sec}) t_i={ti:.3f}, t_j={tj:.3f}: {Dsc[i, j]:.4f}")
print(f"\ntotal: {time.perf_counter()-t0:.1f} s")
