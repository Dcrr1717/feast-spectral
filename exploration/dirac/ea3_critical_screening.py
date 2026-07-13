"""E-A3: relativistic critical screening map lambda_c(Z, kappa).

Yukawa (exponentially screened Coulomb) V(r) = -Z e^{-lam r}/r.  As the
screening lam grows the ground state of channel kappa loses binding;
lambda_c is the threshold where the binding energy E - c^2 reaches 0.

Physics: the short-range behaviour is still -Z/r, so the cusp exponent
gamma = sqrt(kappa^2 - (Z/c)^2) is unchanged and the gamma-basis + exact
cusp quadrature stay exact at r -> 0.  Near threshold the state is
diffuse (small beta) and the binding -> 0, so the SHIFTED pencil
(rest-mass cancellation removed) is essential.

Threshold extrapolation: for the ell=0 ground state (kappa=-1) the
short-range threshold law gives E_b ~ -(lam_c - lam)^2, i.e. sqrt(-E_b)
is LINEAR in lam near threshold; we fit sqrt(-E_b) vs lam on the
near-threshold window and take the zero crossing.  For each channel we
ALSO fit a free power E_b ~ -A (lam_c - lam)^p to report the exponent and
cross-check.  The reduced critical screening s_c = lam_c / Z is the
dimensionless observable; s_c(Z->0) recovers the nonrelativistic constant
(1s: 1.190612..., the Yukawa 1s value), and its rise with Zalpha is the
relativistic enhancement.

Convergence control: lam_c is computed for two basis sizes and two beta
scales; the spread is the honest error bar.  A run-time consistency check
compares against the known nonrelativistic s_c at the smallest Z.
"""

import csv
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from dirac import C_AU, RadialDiracBasis  # noqa: E402

C = C_AU
t0 = time.perf_counter()

# nonrelativistic reduced critical screening s_c = lam_c/Z (Yukawa), for
# the Z->0 consistency anchor only (NOT used to drive sampling):
# 1s (kappa=-1): 1.190612; 2p_{1/2}/2p_{3/2} (kappa=+1,-2): 0.220216
# (nonrel 2p is degenerate; del Valle-Nader 2018).
S_NONREL = {-1: 1.190612, 1: 0.220216, -2: 0.220216}


def ground_binding(Z, kappa, lam, nb, beta, nq=None):
    """Physical ground-state binding of channel kappa.  For kappa>0 the
    naive equal-basis discretization carries ONE spurious state (degenerate
    with the -kappa 1s; see E-A2/FINDINGS #8), which is the deepest gap
    eigenvalue and screening-inert -> the physical ground state is ev[1].
    For kappa<0 there is no spurious state -> ev[0]."""
    nu = Z / C
    gam = float(np.sqrt(kappa * kappa - nu * nu))
    B = RadialDiracBasis(nb, kappa, Z, beta, mu=gam, quad="cusp",
                         nquad=(nq or nb + 150), cert_tol=1e-6)
    B.V = -Z * np.exp(-lam * B.r) / B.r
    ev = B.solve_naive(shifted=True)
    idx = 1 if kappa > 0 else 0
    return float(ev[idx]) if ev.size > idx else np.nan


def sample_curve(Z, kappa, nb, beta_fac):
    """Self-bracketing: scan the reduced screening s upward until the
    physical ground state is lost, collecting (lam, E_b).  beta adapts to
    the shrinking binding so the diffuse tail stays resolved:
    beta ~ max(beta_fac*Z, sqrt(2|E_b|)).  Fine step near threshold."""
    lams, Ebs = [], []
    lost = 0
    for s in np.arange(0.05, 1.60, 0.01):
        lam = s * Z
        Eb_prev = Ebs[-1] if Ebs else -0.05 * Z * Z / (kappa * kappa)
        beta = max(beta_fac * Z, np.sqrt(2 * abs(Eb_prev)))
        Eb = ground_binding(Z, kappa, lam, nb, beta)
        if not np.isfinite(Eb) or Eb >= -1e-9:
            if lams:                       # only stop after we had a state
                lost += 1
                if lost >= 2:
                    break
            continue
        lost = 0
        lams.append(lam)
        Ebs.append(Eb)
    return np.array(lams), np.array(Ebs)


def _lin_zero(x, y):
    a, b = np.polyfit(x, y, 1)
    return -b / a


def extrapolate_lamc(lams, Ebs):
    """Extrapolate E_b -> 0 with the free-power threshold law
    E_b = -A (lam_c - lam)^p (lam_c, A, p all fitted) on the
    near-threshold window.  The exponent is p ~ 2 for ell = 0 and the
    anomalous p ~ 1.25 for ell >= 1 (screened-Coulomb tail); letting p
    float removes the bias that a fixed ansatz introduces.  A p=2
    sqrt-linear estimate is kept as a robust fallback when curve_fit
    fails."""
    if len(lams) < 5:
        return np.nan, np.nan
    k = min(len(lams), 10)
    lw, ew = lams[-k:], Ebs[-k:]
    fallback = _lin_zero(lw, np.sqrt(-ew))          # p = 2 sqrt-linear
    try:
        def model(lam, lamc, A, p):
            return -A * np.maximum(lamc - lam, 1e-30) ** p
        span = lw[-1] - lw[0]
        popt, _ = curve_fit(
            model, lw, ew, p0=(lw[-1] + 0.1 * span, 1.0, 2.0),
            bounds=([lw[-1], 1e-8, 0.5], [lw[-1] + 2 * span, 1e8, 4.0]),
            maxfev=50000)
        return float(popt[0]), float(popt[2])
    except Exception:
        return fallback, 2.0


ZS = [1.0, 40.0, 70.0, 90.0, 110.0, 120.0, 130.0]
KAPPAS = [-1, 1, -2]

rows = []
print(f"c = {C}, V = -Z e^(-lam r)/r, ground state of each channel")
print(f"{'Z':>5} {'kap':>4} {'Zalpha':>7} {'lam_c':>10} {'s_c=lam_c/Z':>11}"
      f" {'p_thresh':>8} {'spread':>9} {'s_c/s_nr':>8}")
for kappa in KAPPAS:
    for Z in ZS:
        nu = Z / C
        if kappa * kappa - nu * nu <= 0.01:      # channel dived / too deep
            continue
        ests = []
        for nb in (90, 120):
            for bf in (0.20, 0.32):
                lams, Ebs = sample_curve(Z, kappa, nb, bf)
                lc, p = extrapolate_lamc(lams, Ebs)
                if np.isfinite(lc):
                    ests.append((lc, p))
        if not ests:
            print(f"{Z:>5.0f} {kappa:>4}  (no bracket)")
            continue
        lcs = np.array([e[0] for e in ests])
        lamc = float(np.median(lcs))
        spread = float(lcs.max() - lcs.min())
        pmed = float(np.nanmedian([e[1] for e in ests]))
        sc = lamc / Z
        ratio = sc / S_NONREL[kappa]
        rows.append(dict(Z=Z, kappa=kappa, Zalpha=nu, lam_c=lamc, s_c=sc,
                         p_threshold=pmed, spread=spread, ratio=ratio))
        print(f"{Z:>5.0f} {kappa:>4} {nu:>7.4f} {lamc:>10.4f} {sc:>11.5f}"
              f" {pmed:>8.3f} {spread:>9.1e} {ratio:>8.4f}")

# consistency: smallest-Z s_c must match the nonrelativistic constant
for kappa in KAPPAS:
    r = [x for x in rows if x['kappa'] == kappa]
    if r:
        r0 = min(r, key=lambda x: x['Z'])
        dev = abs(r0['s_c'] - S_NONREL[kappa])
        print(f"anchor kappa={kappa}: s_c(Z={r0['Z']:.0f})={r0['s_c']:.5f}"
              f" vs nonrel {S_NONREL[kappa]:.5f}  dev={dev:.1e}")

out = Path(__file__).with_name("ea3_critical_screening.csv")
with out.open("w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    for r in rows:
        w.writerow(r)
print(f"\nCSV -> {out.name}   total: {time.perf_counter()-t0:.1f} s")
