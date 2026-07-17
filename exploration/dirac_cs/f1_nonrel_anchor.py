"""Plan B, Phase F1: first relativistic resonance + nonrelativistic
anchor, plus adjudication of the F0 G0c flag.

PRE-REGISTERED PREDICTIONS:
P-F1.0 The F0 'suspicious' gap eigenvalue at -33.09 - 0.32j (Z=80,
       theta=0.3) is the n ~ 10 Rydberg BOUND state distorted by basis
       truncation (diffuse state, real beta), not a mid-gap impostor:
       its real part sits within a few Ha of the exact ladder and its
       |Im| shrinks with nb.
P-F1.1 The 2p Yukawa resonance at Z=1, lam=0.25 reproduces our own
       triple-verified nonrelativistic benchmark E_r = 0.00306445,
       Gamma = 4.369e-3 to better than 1e-4 relative, for BOTH kappa=+1
       (2p1/2) and kappa=-2 (2p3/2); their mutual splitting at Z=1 is
       far below that (O(alpha^2)).
GATE G1: P-F1.1 at 1e-4 or better; otherwise stop and re-derive the
scaled-potential convention.
"""

import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "exploration" / "dirac_cs"))

from dirac import C_AU, RadialDiracBasis, exact_dirac_coulomb  # noqa: E402
from f0_feasibility import cs_pencil  # noqa: E402

C = C_AU
t0 = time.perf_counter()

# ---------------- P-F1.0: adjudicate the G0c flag ----------------------
print("P-F1.0 adjudicating the -33.09-0.32j flag (Z=80, theta=0.3):")
Z = 80.0
gam = float(np.sqrt(1 - (Z / C) ** 2))
exact_n10 = [exact_dirac_coulomb(n, -1, Z) for n in (9, 10, 11)]
print(f"   exact ladder around: n=9..11 -> "
      + ", ".join(f"{e:.3f}" for e in exact_n10))
for nb in (100, 160):
    B = RadialDiracBasis(nb, -1, Z, 0.8 * Z, mu=gam, quad="cusp",
                         nquad=nb + 120, cert_tol=1e-6)
    Hs, _ = cs_pencil(B, 0.3)
    ev = np.linalg.eigvals(Hs)
    ev = ev[(ev.real > -40) & (ev.real < -25)]
    ev = ev[np.argsort(np.abs(ev.imag))][:3]
    print(f"   nb={nb}: gap eigenvalues near -33: "
          + ", ".join(f"{e:.4f}" for e in np.sort_complex(ev)))

# ---------------- P-F1.1: the nonrelativistic anchor -------------------
print("\nP-F1.1 Yukawa 2p resonance, Z=1, lam=0.25"
      " (nonrel benchmark: E_r=0.00306445, Gamma=4.369e-3):")
ENR, GNR = 0.00306445, 4.369e-3


def resonance(kappa, nb=160, beta=0.25, thetas=np.arange(0.15, 0.55, 0.05)):
    nu = 1.0 / C
    gam = float(np.sqrt(kappa * kappa - nu * nu))
    B = RadialDiracBasis(nb, kappa, 1.0, beta, mu=gam, quad="cusp",
                         nquad=nb + 120, cert_tol=1e-6)
    traj = []
    for th in thetas:
        Hs, _ = cs_pencil(B, th, lam=0.25)
        ev = np.linalg.eigvals(Hs)
        win = ev[(ev.real > 5e-4) & (ev.real < 0.02)
                 & (ev.imag < -1e-5) & (ev.imag > -0.05)]
        if win.size:
            traj.append(win[np.argmin(np.abs(win - (ENR - 0.5j * GNR)))])
        else:
            traj.append(np.nan + 0j)
    traj = np.array(traj)
    # stationarity: minimal |dE/dtheta|
    d = np.abs(np.diff(traj))
    j = int(np.nanargmin(d))
    return traj[j + 1], d[j], traj


for kappa, name in ((1, "2p1/2 (kappa=+1)"), (-2, "2p3/2 (kappa=-2)")):
    E, stat, traj = resonance(kappa)
    Er, Ga = E.real, -2 * E.imag
    print(f"   {name}: E_r = {Er:.8f}  Gamma = {Ga:.6e}"
          f"   |dE/dth| = {stat:.1e}")
    print(f"      vs nonrel: dE_r/E_r = {(Er - ENR) / ENR:+.2e}"
          f"   dGamma/Gamma = {(Ga - GNR) / GNR:+.2e}")

# fine-structure splitting at Z=1 (prediction: tiny)
E1, _, _ = resonance(1)
E2, _, _ = resonance(-2)
print(f"\n   2p1/2 - 2p3/2 splitting at Z=1:"
      f" dE_r = {E1.real - E2.real:+.2e},"
      f" dGamma = {-2 * (E1.imag - E2.imag):+.2e}")
print(f"\ntotal: {time.perf_counter()-t0:.1f} s")
