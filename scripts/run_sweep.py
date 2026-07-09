"""E8: resonance sweep E_r(lambda), Gamma(lambda) for the 2p state from just
above critical screening (lambda_c = 0.220216806605, del Valle & Nader 2018)
up to lambda = 0.30, via complex-scaling theta-trajectories.

Basis parameters chosen by regime: near threshold the resonance is diffuse
(small beta, larger N); away from threshold it is compact. The external
anchor at lambda = 0.22131705 (Alhaidari 2008 / Bylicki 2007) is included.

Outputs: results/E8_sweep.csv and figures/fig4_sweep.(pdf|png)
(publication style: no embedded titles, 11 pt fonts, B/W friendly).
"""

import sys
import time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from basis import BasisSet                                   # noqa: E402
from run_resonances import find_resonance                    # noqa: E402

RES, FIG = ROOT / "results", ROOT / "figures"

plt.rcParams.update({
    "font.family": "serif", "font.size": 11,
    "axes.labelsize": 11, "xtick.labelsize": 10, "ytick.labelsize": 10,
    "legend.fontsize": 9.5, "axes.grid": True, "grid.alpha": 0.3,
    "figure.dpi": 150,
})

LAM_ANCHOR = 0.22131705      # published 12-digit benchmark point
LAMS = [0.2210, LAM_ANCHOR, 0.2225, 0.2250, 0.2300, 0.2400,
        0.2500, 0.2600, 0.2800, 0.3000]


def sweep():
    rows = ["lambda,N,beta,theta_opt,E_r_Ha,Gamma_Ha,stationarity"]
    out = []
    # two basis regimes, each built once and reused across lambda
    regimes = {
        "near": dict(N=250, beta=0.08,
                     thetas=np.arange(0.10, 0.46, 0.025),
                     win=((1e-6, 2e-3), (-2e-3, 1e-7))),
        "far":  dict(N=200, beta=0.30,
                     thetas=np.arange(0.10, 0.61, 0.05),
                     win=((1e-4, 5e-1), (-3e-1, 1e-4))),
    }
    built = {}
    for lam in LAMS:
        key = "near" if lam < 0.235 else "far"
        cfg = regimes[key]
        if key not in built:
            B = BasisSet(cfg["N"], 1, cfg["beta"])
            built[key] = (B, B.kinetic().astype(complex))
        B, T = built[key]
        spectra = []
        for th in cfg["thetas"]:
            H = np.exp(-2j * th) * T + B.yukawa(lam=lam, theta=th)
            spectra.append(np.linalg.eigvals(H))
        got = find_resonance(spectra, cfg["thetas"],
                             re_win=cfg["win"][0], im_win=cfg["win"][1])
        if got is None:
            print(f"lam={lam}: no stationary point")
            continue
        stat, Eres, th_opt, _ = got
        rows.append(f"{lam},{cfg['N']},{cfg['beta']},{th_opt:.3f},"
                    f"{Eres.real:.10e},{-2*Eres.imag:.10e},{stat:.1e}")
        out.append((lam, Eres.real, -2 * Eres.imag))
        print(f"lam={lam:.6f}: E_r={Eres.real:.6e}  Gamma={-2*Eres.imag:.6e}"
              f"  (stat {stat:.0e}, th*={th_opt:.2f})")
    (RES / "E8_sweep.csv").write_text("\n".join(rows), encoding="utf-8")
    return np.array(out)


def plot(data):
    lam, Er, Ga = data.T
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.7))
    axes[0].plot(lam, Er, "o-", ms=4.5, lw=1.2, color="k")
    axes[0].set_xlabel(r"screening parameter $\lambda$  (a.u.)")
    axes[0].set_ylabel(r"$E_r$  (Ha)")
    axes[1].semilogy(lam, Ga, "s--", ms=4.5, lw=1.2, color="k")
    axes[1].set_xlabel(r"screening parameter $\lambda$  (a.u.)")
    axes[1].set_ylabel(r"$\Gamma$  (Ha)")
    for ax in axes:
        ax.axvline(0.220216806605, color="0.5", lw=0.9, ls=":")
    j = np.argmin(np.abs(lam - LAM_ANCHOR))
    axes[0].plot(lam[j], Er[j], "*", ms=13, mfc="none", mec="k",
                 label="published benchmark")
    axes[1].plot(lam[j], Ga[j], "*", ms=13, mfc="none", mec="k",
                 label="published benchmark")
    axes[0].legend()
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"fig4_sweep.{ext}")
    plt.close(fig)


if __name__ == "__main__":
    t0 = time.perf_counter()
    data = sweep()
    plot(data)
    print(f"total wall time: {time.perf_counter()-t0:.1f} s")
