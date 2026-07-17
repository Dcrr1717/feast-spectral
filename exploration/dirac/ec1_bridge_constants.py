"""Plan C, E-C1: measure the bridge-lemma constants that decide the
T1-Coulomb route BEFORE investing theory (math_daniel protocol).

Objects on the small space (first nb of a nested gamma-basis reference,
kappa=-1, Coulomb; reference size nb_ref as continuum proxy, now
unconstrained thanks to the log-space cusp rule):

  A1 = |P H0 P|          (fiber-diagonal; the P2 route uses this)
  A2 = (P H0^2 P)^{1/2}  (operator-convexity intermediary:
                          A1 <= A2 ALWAYS, and Pi|D|Pi <= A2)
  A3 = P |H0_ref| P      (compressed absolute value; Hardy-compatible)

Measured quantities (Coulomb W = diag(Vd,Vd)):
  alpha_i(nb; nb_ref) = || A_i^{-1/2} W A_i^{-1/2} ||    i = 1,2,3
  c_bridge = || A1^{-1/2} A3 A1^{-1/2} ||   (cost of moving the form
             from the Hardy-compatible A3 to the fiber-diagonal A1)
  per-fiber diagnostic: <b_j, A3 b_j> / E_j(A1) vs t_j (where the
             bridge hurts: low or high modes?)

PRE-REGISTERED PREDICTIONS:
P-C1.1 alpha_1 diverges (log, known from FINDINGS #13); alpha_3
       converges (E-T2); the NEW datum is alpha_2: if it converges the
       convexity route (a) is alive.
P-C1.2 operator orders hold numerically: A1 <= A2 and A3 <= A2.
P-C1.3 c_bridge grows with nb (A3 sees momentum beyond the sigma_max
       cap); the per-fiber ratio localizes the damage in the HIGH-t
       fibers, suggesting a split-window strategy for the theorem.
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


def build_ref(Z, nb_ref, beta_fac=0.8):
    nu = Z / C
    gam = float(np.sqrt(1 - nu * nu))
    B = RadialDiracBasis(nb_ref, -1, Z, beta_fac * Z, mu=gam, quad="cusp",
                         nquad=nb_ref + 120, cert_tol=1e-5)
    K = B._me(B.f, 2.0 * C * B.h)
    H0 = np.block([[C**2 * np.eye(nb_ref), C * K],
                   [C * K.T, -C**2 * np.eye(nb_ref)]])
    Vd = B._mev(B.f, B.V, B.f)
    Vd = 0.5 * (Vd + Vd.T)
    return H0, Vd


def sqrtm_sym(M):
    ev, U = np.linalg.eigh(M)
    return U @ np.diag(np.sqrt(np.maximum(ev, 0))) @ U.T


def invsqrtm_sym(M):
    ev, U = np.linalg.eigh(M)
    return U @ np.diag(1.0 / np.sqrt(np.maximum(ev, 1e-300))) @ U.T


def analyze(Z, nb, nb_ref):
    H0, Vdf = build_ref(Z, nb_ref, 0.8)
    n = nb_ref
    sel = np.concatenate([np.arange(nb), n + np.arange(nb)])
    PH0P = 0.5 * (H0[np.ix_(sel, sel)] + H0[np.ix_(sel, sel)].T)
    evs, Us = np.linalg.eigh(PH0P)
    A1 = Us @ np.diag(np.abs(evs)) @ Us.T
    H0sq = H0 @ H0
    A2 = sqrtm_sym(0.5 * (H0sq[np.ix_(sel, sel)]
                          + H0sq[np.ix_(sel, sel)].T))
    evr, Ur = np.linalg.eigh(H0)
    absH = Ur @ np.diag(np.abs(evr)) @ Ur.T
    A3 = 0.5 * (absH[np.ix_(sel, sel)] + absH[np.ix_(sel, sel)].T)
    Vd = Vdf[:nb, :nb]
    W = np.block([[Vd, np.zeros((nb, nb))], [np.zeros((nb, nb)), Vd]])
    alphas = []
    for A in (A1, A2, A3):
        Mi = invsqrtm_sym(A)
        alphas.append(float(np.linalg.norm(Mi @ W @ Mi, 2)))
    Mi1 = invsqrtm_sym(A1)
    c_bridge = float(np.linalg.norm(Mi1 @ A3 @ Mi1, 2))
    # operator-order checks (min eig of differences)
    o12 = float(np.min(np.linalg.eigvalsh(A2 - A1)))
    o32 = float(np.min(np.linalg.eigvalsh(A2 - A3)))
    # per-fiber bridge ratio: fibers of the SMALL operator
    Ksm = PH0P[:nb, nb:] / C
    U, S, Vt = np.linalg.svd(Ksm)
    ratios, ts = [], []
    for j in range(nb):
        b = np.concatenate([U[:, j], np.zeros(nb)])   # electronic-ish probe
        Ej = float(np.sqrt(C**4 + C**2 * S[j]**2))
        ratios.append(float(b @ (A3 @ b)) / Ej)
        ts.append(S[j] / (2 * C))
    return alphas, c_bridge, o12, o32, np.array(ts), np.array(ratios)


print("E-C1 bridge constants (kappa=-1, Coulomb, gamma-basis):")
print(f"{'Z':>4} {'nb':>4} {'ref':>4} {'a1(|PHP|)':>10} {'a2(conv)':>9}"
      f" {'a3(P|H|P)':>10} {'c_bridge':>9} {'A2-A1>=0':>9} {'A2-A3>=0':>9}")
for Z in (40.0, 92.0):
    for nb, nb_ref in ((40, 160), (80, 320), (120, 480)):
        al, cb, o12, o32, ts, rat = analyze(Z, nb, nb_ref)
        print(f"{Z:>4.0f} {nb:>4} {nb_ref:>4} {al[0]:>10.4f} {al[1]:>9.4f}"
              f" {al[2]:>10.4f} {cb:>9.3f} {o12:>9.1e} {o32:>9.1e}")

# per-fiber localization at the largest case
al, cb, o12, o32, ts, rat = analyze(92.0, 120, 480)
o = np.argsort(ts)
print("\nper-fiber <b,A3 b>/E_j(A1) (Z=92, nb=120, ref=480):")
for k in (0, len(o) // 4, len(o) // 2, 3 * len(o) // 4, len(o) - 1):
    j = o[k]
    print(f"   t_j = {ts[j]:.4f}:  ratio = {rat[j]:.3f}")
print(f"\ntotal: {time.perf_counter()-t0:.1f} s")
