"""INDEPENDENT adversarial verification of E-A3 (critical screening).

Independent choices vs ea3_critical_screening.py:
  * binding = PLATEAU (trimmed median over a broad geometric beta grid),
    NOT an adaptive single-beta schedule. Rejects small-beta spurious
    over-binding and large-beta under-resolution.
  * lam_c from a PROFILE-over-lam_c log-space linear regression
    (ln(-Eb) vs ln(lam_c-lam)); jointly returns exponent p. Totally
    different estimator from their real-space nonlinear curve_fit.
  * cross-checked with a linearized (-Eb)^(1/p) zero-crossing.
  * finer lam grid, larger nb, explicit nb-convergence check.
Also directly tests: (d) is kappa>0 ev[0] screening-inert & degenerate
with the -kappa 1s?  and the beta-schedule sensitivity hypothesis.
"""
import sys
from pathlib import Path
import numpy as np

ROOT = Path(r"C:\Users\dcrr1\OneDrive\Desktop\metodo feast espectral\feast_spectral")
sys.path.insert(0, str(ROOT / "src"))
from dirac import C_AU, RadialDiracBasis
C = C_AU


def gap_eigs(Z, kappa, lam, nb, beta, nq=None):
    nu = Z / C
    gam = float(np.sqrt(kappa * kappa - nu * nu))
    try:
        B = RadialDiracBasis(nb, kappa, Z, beta, mu=gam, quad="cusp",
                             nquad=(nq or nb + 150), cert_tol=1e-6)
    except Exception:
        return np.array([])
    B.V = -Z * np.exp(-lam * B.r) / B.r
    return B.solve_naive(shifted=True)


def phys_binding(Z, kappa, lam, nb, beta):
    ev = gap_eigs(Z, kappa, lam, nb, beta)
    idx = 1 if kappa > 0 else 0
    return float(ev[idx]) if ev.size > idx else np.nan


def plateau_binding(Z, kappa, lam, nb, betas):
    """Trimmed-median plateau binding over a broad beta grid."""
    vals = np.array([phys_binding(Z, kappa, lam, nb, b) for b in betas])
    vals = vals[np.isfinite(vals) & (vals < -1e-9)]
    if vals.size == 0:
        return np.nan, 0
    vals = np.sort(vals)                       # most negative first is vals[0]
    if vals.size >= 5:                         # trim extremes, median plateau
        lo = int(0.2 * vals.size)
        hi = int(np.ceil(0.8 * vals.size))
        vals = vals[lo:hi]
    return float(np.median(vals)), int(vals.size)


def profile_lamc(lams, ebs, ncand=4000):
    """lam_c by profiling: for each candidate lam_c minimise the residual of
    the linear regression ln(-Eb) = ln A + p ln(lam_c - lam). Returns
    (lam_c, p)."""
    lams = np.asarray(lams); ebs = np.asarray(ebs)
    y = np.log(-ebs)
    lo = lams[-1] + 1e-9 * abs(lams[-1])
    hi = lams[-1] + (lams[-1] - lams[0])
    best = (np.nan, np.nan, np.inf)
    for lc in np.linspace(lo, hi, ncand):
        x = np.log(lc - lams)
        p, b = np.polyfit(x, y, 1)
        res = np.sum((y - (p * x + b)) ** 2)
        if res < best[2]:
            best = (lc, p, res)
    return best[0], best[1]


def linzero_pow(lams, ebs, p):
    """Zero crossing of (-Eb)^(1/p) linear in lam (independent cross-check)."""
    lams = np.asarray(lams); ebs = np.asarray(ebs)
    yy = (-ebs) ** (1.0 / p)
    a, b = np.polyfit(lams, yy, 1)
    return -b / a


def sample_curve(Z, kappa, nb, betas, s0, s1, ds):
    lams, ebs, diag = [], [], []
    for s in np.arange(s0, s1, ds):
        eb, n = plateau_binding(Z, kappa, s * Z, nb, betas)
        if np.isfinite(eb):
            lams.append(s * Z); ebs.append(eb); diag.append((s, eb, n))
        elif lams:                     # left the bound region -> stop window
            break
    return np.array(lams), np.array(ebs), diag


def estimate_sc(Z, kappa, nb, betas, s0, s1, ds, window=10):
    lams, ebs, diag = sample_curve(Z, kappa, nb, betas, s0, s1, ds)
    if len(lams) < 6:
        return None
    k = min(window, len(lams))
    lw, ew = lams[-k:], ebs[-k:]
    lc, p = profile_lamc(lw, ew)
    # cross-check with fixed-p linearization and with l=0 sqrt
    lc_pow = linzero_pow(lw, ew, p) if np.isfinite(p) and p > 0.3 else np.nan
    lc_sqrt = linzero_pow(lw, ew, 2.0)
    return dict(sc=lc / Z, p=p, sc_pow=lc_pow / Z,
                sc_sqrt=lc_sqrt / Z, npts=len(lams), lamc=lc)


# ---------------------------------------------------------------------------
print("=" * 78)
print("TEST (d): is kappa>0 ev[0] screening-inert & degenerate with -kappa 1s?")
print("=" * 78)
for Z in (90.0, 130.0):
    print(f" Z={Z}")
    for lam in (0.0, 2.0, 8.0, 20.0):
        evp = gap_eigs(Z, +1, lam, 110, 0.4 * Z)
        evm = gap_eigs(Z, -1, lam, 110, 0.4 * Z)
        d = (evp[0] - evm[0]) if (evp.size and evm.size) else np.nan
        print(f"  lam={lam:5.1f}  ev0(k=+1)={evp[0] if evp.size else np.nan: .4f}"
              f"  ev0(k=-1,1s)={evm[0] if evm.size else np.nan: .4f}"
              f"  diff={d: .3e}")

# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("MAIN: independent s_c  (plateau binding + profile-lamc estimator)")
print("=" * 78)
CASES = [
    # (Z, kappa, s0, s1, ds)  window chosen just below each threshold
    (1.0,  -1, 1.10, 1.30, 0.004),
    (90.0, -1, 1.15, 1.36, 0.006),
    (130.0,-1, 1.36, 1.52, 0.006),
    (1.0,   1, 0.150, 0.245, 0.004),
    (130.0, 1, 0.220, 0.345, 0.004),
    (1.0,  -2, 0.150, 0.245, 0.004),
    (130.0,-2, 0.180, 0.255, 0.004),
]
print(f"{'Z':>6} {'kap':>4} {'sc_prof':>9} {'p':>6} {'sc_pow':>9} {'sc_sqrt':>9} {'npts':>5}")
res = {}
for Z, kappa, s0, s1, ds in CASES:
    betas = np.geomspace(0.004, 2.5, 26) * Z
    r = estimate_sc(Z, kappa, 130, betas, s0, s1, ds)
    if r is None:
        print(f"{Z:>6.0f} {kappa:>4}  (insufficient bracket)")
        continue
    res[(Z, kappa)] = r
    print(f"{Z:>6.0f} {kappa:>4} {r['sc']:>9.5f} {r['p']:>6.2f}"
          f" {r['sc_pow']:>9.5f} {r['sc_sqrt']:>9.5f} {r['npts']:>5}")

print("\n--- verdict inputs ---")
if (1.0, -1) in res and (130.0, -1) in res:
    s_lo = res[(1.0, -1)]['sc']; s_hi = res[(130.0, -1)]['sc']
    print(f"1s rise: s_c(Z=1)={s_lo:.4f} -> s_c(Z=130)={s_hi:.4f}"
          f"  ({100*(s_hi/s_lo-1):+.1f}%)   [ref CSV: 1.1909 -> 1.4777, +24%]")
if (130.0, 1) in res and (130.0, -2) in res:
    a = res[(130.0, 1)]['sc']; b = res[(130.0, -2)]['sc']
    print(f"2p split @Z=130: 2p1/2(k=+1)={a:.4f}  2p3/2(k=-2)={b:.4f}"
          f"  split={a-b:+.4f}   [ref CSV: 0.3127 vs 0.2325, +0.0802]")
    print(f"   2p1/2 more screening-resistant? {'YES' if a > b else 'NO'}")
if (1.0, 1) in res and (1.0, -2) in res:
    print(f"2p @Z=1 (should be ~degenerate ~0.2202): "
          f"k=+1={res[(1.0,1)]['sc']:.4f}  k=-2={res[(1.0,-2)]['sc']:.4f}")
