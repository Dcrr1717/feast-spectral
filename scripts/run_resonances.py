"""Part B: critical screening parameters and resonances of the Yukawa potential.

E5  Critical screening lambda_c for 1s, 2s, 3s (ell=0), 2p, 3p (ell=1), 3d
    (ell=2): the value of lambda where the state stops being bound. Found by
    bisection on the *count* of negative eigenvalues (robust against the
    discretized continuum hovering near E = 0). Uncertainty estimated from
    two basis sizes.

E6  Resonances above critical screening via complex scaling. For lambda >
    lambda_c(2p) the 2p state turns into a resonance E = E_r - i*Gamma/2.
    Procedure: (i) theta-trajectory analysis -- eigenvalues of H(theta) for a
    grid of rotation angles; continuum states rotate, the resonance stays
    stationary; (ii) the stationary point is refined and cross-checked with
    feast_region (contour integration around the complex eigenvalue);
    (iii) stability vs basis size N gives the quoted uncertainty.

Outputs: results/E5_critical_screening.csv, results/E6_resonances.csv,
figures/E6_resonance_trajectory.(pdf|png).
"""

import sys
import time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.linalg import eigvalsh

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from basis import BasisSet            # noqa: E402
from feast import feast_region        # noqa: E402

RES = ROOT / "results"
FIG = ROOT / "figures"
RES.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)

plt.rcParams.update({
    "font.family": "serif", "font.size": 10, "axes.grid": True,
    "grid.alpha": 0.3, "figure.dpi": 150,
})

E_NEG = -1e-7   # "bound" threshold against finite-basis noise


# ----------------------------------------------------------------- E5 --
def count_bound(B: BasisSet, T: np.ndarray, lam: float) -> int:
    ev = eigvalsh(T + B.yukawa(lam=lam))
    return int(np.sum(ev < E_NEG))


def critical_lambda(ell: int, want: int, lo: float, hi: float,
                    N: int, beta: float, bits: int = 26) -> float:
    """Bisect lambda where the number of bound ell-states drops below `want`."""
    B = BasisSet(N, ell, beta)
    T = B.kinetic()
    assert count_bound(B, T, lo) >= want, f"bracket lo fails (ell={ell}, want={want})"
    assert count_bound(B, T, hi) < want, f"bracket hi fails (ell={ell}, want={want})"
    for _ in range(bits):
        mid = 0.5 * (lo + hi)
        if count_bound(B, T, mid) >= want:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def e5_critical_screening():
    print("[E5] Critical screening parameters")
    cases = [  # state, ell, count, bracket
        ("1s", 0, 1, 1.00, 1.40),
        ("2s", 0, 2, 0.25, 0.40),
        ("3s", 0, 3, 0.10, 0.20),
        ("2p", 1, 1, 0.15, 0.30),
        ("3p", 1, 2, 0.08, 0.16),
        ("3d", 2, 1, 0.05, 0.14),
    ]
    beta = 0.05
    rows = ["state,ell,lambda_c_N250,lambda_c_N320,delta"]
    out = {}
    for st, ell, want, lo, hi in cases:
        l1 = critical_lambda(ell, want, lo, hi, N=250, beta=beta)
        l2 = critical_lambda(ell, want, lo, hi, N=320, beta=beta)
        rows.append(f"{st},{ell},{l1:.6f},{l2:.6f},{abs(l1-l2):.1e}")
        out[st] = (l2, abs(l1 - l2))
        print(f"   {st}: lambda_c = {l2:.6f}  (basis-size drift {abs(l1-l2):.1e})")
    (RES / "E5_critical_screening.csv").write_text("\n".join(rows), encoding="utf-8")
    return out


# ----------------------------------------------------------------- E6 --
def theta_trajectory(lam: float, ell: int, N: int, beta: float,
                     thetas: np.ndarray):
    """Eigenvalues of H(theta) for each theta (full complex spectra)."""
    B = BasisSet(N, ell, beta)
    T = B.kinetic().astype(complex)
    spectra = []
    for th in thetas:
        H = np.exp(-2j * th) * T + B.yukawa(lam=lam, theta=th)
        spectra.append(np.linalg.eigvals(H))
    return B, T, spectra


def find_resonance(spectra, thetas, re_win=(1e-4, 0.5), im_win=(-0.3, 1e-4)):
    """Track eigenvalue trajectories through theta; return the most
    theta-stationary one (the resonance) and its trajectory."""
    cands = [z for z in spectra[0]
             if re_win[0] < z.real < re_win[1] and im_win[0] < z.imag < im_win[1]]
    best = None
    for z0 in cands:
        traj = [z0]
        for spec in spectra[1:]:
            traj.append(spec[np.argmin(np.abs(spec - traj[-1]))])
        traj = np.array(traj)
        steps = np.abs(np.diff(traj))
        j = int(np.argmin(steps))            # most stationary segment
        score = steps[j]
        if best is None or score < best[0]:
            best = (score, 0.5 * (traj[j] + traj[j + 1]),
                    0.5 * (thetas[j] + thetas[j + 1]), traj)
    return best  # (stationarity, E_res, theta_opt, trajectory)


def e6_resonances():
    print("[E6] Yukawa resonances above critical screening (ell = 1, 2p)")
    thetas = np.arange(0.10, 0.61, 0.05)
    rows = ["lambda,N,beta,theta_opt,E_r_Ha,Gamma_Ha,stationarity,"
            "feast_E_r,feast_Gamma,feast_resid"]
    plot_payload = None
    for lam in (0.25, 0.28):
        for N, beta in ((150, 0.30), (200, 0.30), (200, 0.45)):
            B, T, spectra = theta_trajectory(lam, 1, N, beta, thetas)
            got = find_resonance(spectra, thetas)
            if got is None:
                print(f"   lam={lam}, N={N}, beta={beta}: no candidate found")
                continue
            stat, Eres, th_opt, traj = got
            # FEAST cross-check on a small circle around the resonance,
            # at a rotation angle safely past the stationary point
            th_f = min(th_opt + 0.10, thetas[-1])
            H = np.exp(-2j * th_f) * T + B.yukawa(lam=lam, theta=th_f)
            rad = max(6.0 * stat, 8e-4)
            lamF, XF, resF, itF, _ = feast_region(H, Eres, rad, m0=6)
            if lamF.size:
                jF = int(np.argmin(np.abs(lamF - Eres)))
                EF, rF = lamF[jF], float(resF[jF])
            else:
                EF, rF = np.nan + 1j * np.nan, np.nan
            rows.append(
                f"{lam},{N},{beta},{th_opt:.3f},{Eres.real:.8f},"
                f"{-2*Eres.imag:.3e},{stat:.2e},"
                f"{EF.real:.8f},{-2*EF.imag:.3e},{rF:.1e}")
            print(f"   lam={lam}, N={N}, beta={beta}: "
                  f"E_r={Eres.real:.6f}, Gamma={-2*Eres.imag:.3e} "
                  f"(theta*={th_opt:.2f}, feast agrees to "
                  f"{abs(EF-Eres):.1e}, resid {rF:.0e})")
            if lam == 0.25 and N == 200 and beta == 0.30:
                plot_payload = (lam, thetas, spectra, traj, Eres)
    (RES / "E6_resonances.csv").write_text("\n".join(rows), encoding="utf-8")

    if plot_payload:
        lam, thetas, spectra, traj, Eres = plot_payload
        fig, ax = plt.subplots(figsize=(5.6, 4.2))
        colors = plt.cm.viridis(np.linspace(0, 1, len(thetas)))
        for th, spec, ccol in zip(thetas, spectra, colors):
            ax.plot(spec.real, spec.imag, ".", ms=3, color=ccol, alpha=0.6)
        ax.plot(traj.real, traj.imag, "r.-", lw=1.2, ms=5,
                label="resonance trajectory")
        ax.plot(Eres.real, Eres.imag, "k*", ms=12, label="stationary point")
        ax.set_xlim(-0.55, 0.25)
        ax.set_ylim(-0.35, 0.02)
        ax.set_xlabel(r"Re $E$  (Ha)")
        ax.set_ylabel(r"Im $E$  (Ha)")
        ax.set_title(rf"Complex-scaled spectrum, Yukawa $\lambda={lam}$, $\ell=1$"
                     "\n(color: $\\theta = 0.10 \\to 0.60$)")
        ax.legend(fontsize=8, loc="lower left")
        fig.tight_layout()
        for ext in ("pdf", "png"):
            fig.savefig(FIG / f"E6_resonance_trajectory.{ext}")
        plt.close(fig)


if __name__ == "__main__":
    t0 = time.perf_counter()
    e5_critical_screening()
    e6_resonances()
    print(f"total wall time: {time.perf_counter() - t0:.1f} s")
