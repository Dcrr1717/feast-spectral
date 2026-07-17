"""Plan B, Phase F3: production tables and the spin-symmetry-breaking map.

E-B3a  Reference table for the 2p doublet (kappa = +1: 2p1/2,
       kappa = -2: 2p3/2) of the Dirac-Yukawa problem, Z = 1..130,
       reduced screening s = lambda/Z in {0.24, 0.30}.
       HONEST-DIGITS PROTOCOL:每 state evaluated on a 2x2x2 grid
       (theta pair x nb pair x beta pair); quoted error = max spread;
       Gamow-shooting anchor where Gamma <~ 4 (the double-precision
       envelope of the incoming-coefficient extraction).
E-B3b  The flagship map: fine-structure splitting of widths vs Zalpha.
       Expectation from E-A3/F2: percent-level at low Z, ORDERS OF
       MAGNITUDE at strong coupling; partner may even remain BOUND
       (s < s_c(2p1/2)) while the other resonates -- report that too.
E-B3c  Complex-plane trajectories vs s at Z = 80 for both partners.

Pre-registered: P3.1 both partners converge with spread <= 1e-8 for
narrow/moderate states; P3.2 the width ratio Gamma(3/2)/Gamma(1/2) at
fixed s grows monotonically with Zalpha.
"""

import csv
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "exploration" / "dirac_cs"))

from dirac import C_AU, RadialDiracBasis  # noqa: E402
from cs_core import cs_pencil, resonance_eigs  # noqa: E402
from f2_shooting import refine  # noqa: E402

C = C_AU
t0 = time.perf_counter()


def measure(Z, kappa, lam, re_win=(1e-6, 200.0), im_win=(-80.0, 1e-12)):
    """Honest-digits resonance: 2 theta x 2 nb x 2 beta grid; returns
    (value from the largest basis, max spread, n_found)."""
    nu = Z / C
    gam = float(np.sqrt(kappa * kappa - nu * nu))
    vals = []
    for nb in (200, 240):
        for bf in (0.25, 0.40):
            B = RadialDiracBasis(nb, kappa, Z, bf * max(Z, 1.0), mu=gam,
                                 quad="cusp", nquad=nb + 100, cert_tol=1e-5)
            picks = []
            for th in (0.45, 0.60):
                Hs = cs_pencil(B, th, lam=lam)
                win = resonance_eigs(Hs, re_win, im_win)
                if win.size:
                    picks.append(win)
            if len(picks) == 2:
                d = np.abs(picks[1][None, :] - picks[0][:, None])
                i, j = np.unravel_index(np.argmin(d), d.shape)
                if d[i, j] < 1e-2 * max(abs(picks[1][j]), 1e-3):
                    vals.append(picks[1][j])
    if not vals:
        return None, np.inf, 0
    vals = np.array(vals)
    ref = vals[-1]
    keep = vals[np.abs(vals - ref) < 0.05 * max(abs(ref), 1e-3)]
    return ref, float(np.max(np.abs(keep - ref))) if keep.size > 1 else np.inf, keep.size


rows = []
print("E-B3a/b production (2p doublet), honest-digits protocol:")
print(f"{'Z':>4} {'s':>5} {'kap':>4} {'E_r':>16} {'Gamma':>13}"
      f" {'spread':>9} {'shoot':>9}")
for Z in (1.0, 20.0, 40.0, 60.0, 80.0, 92.0, 110.0, 130.0):
    nu = Z / C
    for s in (0.24, 0.30):
        lam = s * Z
        pair = {}
        for kappa in (1, -2):
            e, spread, n = measure(Z, kappa, lam)
            tag = ""
            sh = np.nan
            if e is not None and -2 * e.imag < 4.0 and e.real > 1e-5:
                try:
                    Rr = max(18.0 / lam, 2.0) if lam > 0 else 20.0
                    es = refine(complex(e) + C**2, Z, lam, kappa, Rr) - C**2
                    sh = abs(e - es)
                except Exception:
                    pass
            if e is None:
                print(f"{Z:>4.0f} {s:>5.2f} {kappa:>4}   (no resonance --"
                      f" bound or dissolved)")
                rows.append(dict(Z=Z, s=s, kappa=kappa, Er=np.nan,
                                 Gamma=np.nan, spread=np.nan, shoot=np.nan))
                continue
            pair[kappa] = e
            print(f"{Z:>4.0f} {s:>5.2f} {kappa:>4} {e.real:>16.9f}"
                  f" {-2 * e.imag:>13.6e} {spread:>9.1e} {sh:>9.1e}")
            rows.append(dict(Z=Z, s=s, kappa=kappa, Er=e.real,
                             Gamma=-2 * e.imag, spread=spread, shoot=sh))
        if len(pair) == 2:
            g1, g2 = -2 * pair[1].imag, -2 * pair[-2].imag
            print(f"{'':>14} doublet: dE_r = {pair[-2].real - pair[1].real:+.4e}"
                  f"   Gamma(3/2)/Gamma(1/2) = {g2 / g1:.3f}")

out = Path(__file__).with_name("f3_table.csv")
with out.open("w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    for r in rows:
        w.writerow(r)

# ---------------- E-B3c: trajectories at Z=80 ---------------------------
print("\nE-B3c trajectories (Z=80, s sweep):")
traj = {1: [], -2: []}
for s in (0.245, 0.26, 0.28, 0.30, 0.33, 0.36):
    for kappa in (1, -2):
        e, spread, _ = measure(80.0, kappa, s * 80.0)
        if e is not None:
            traj[kappa].append((s, e.real, -2 * e.imag))
for kappa, name in ((1, "2p1/2"), (-2, "2p3/2")):
    print(f"   {name}: " + "  ".join(
        f"s={s:.3f}:({er:.3f},{ga:.3f})" for s, er, ga in traj[kappa]))

with Path(__file__).with_name("f3_traj.csv").open("w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["kappa", "s", "Er", "Gamma"])
    for kappa, pts in traj.items():
        for s, er, ga in pts:
            w.writerow([kappa, s, er, ga])
print(f"\ntotal: {time.perf_counter()-t0:.1f} s")
