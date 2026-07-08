"""Phase-1 validation suite.

E1  Hydrogen, orthonormal Laguerre basis: exponential convergence of the
    first 6 eigenvalues (ell = 0) and 5 eigenvalues (ell = 1) vs exact
    E_n = -1/(2 n^2) Hartree. Basis-scale scan over beta.
E2  FEAST vs LAPACK (scipy.linalg.eigh) on the same Hamiltonian: agreement,
    residuals, contour-size sensitivity.
E3  Screened Coulomb (Debye-Hueckel / Yukawa): Laguerre-basis energies
    cross-validated against an independent finite-difference solver.
E4  Complex scaling H(theta): bound-state energies must stay real and
    theta-independent (internal consistency check of the machinery that
    will later compute resonances).

All outputs land in results/ (CSV, TXT) and figures/ (PDF+PNG, English).
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

from basis import BasisSet                      # noqa: E402
from feast import feast_interval, _self_test    # noqa: E402
import reference_fd as fd                       # noqa: E402

RES = ROOT / "results"
FIG = ROOT / "figures"
RES.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)

plt.rcParams.update({
    "font.family": "serif", "font.size": 10, "axes.grid": True,
    "grid.alpha": 0.3, "figure.dpi": 150,
})


def exact_hydrogen(nmax, ell):
    n = np.arange(ell + 1, nmax + 1)
    return -0.5 / n**2


# ----------------------------------------------------------------- E1 --
def e1_hydrogen_convergence():
    print("[E1] Hydrogen: Laguerre-basis convergence")
    sizes = [20, 30, 40, 60, 80, 100, 120, 150, 180]
    betas = [1.0, 0.5, 1.0 / 3.0, 1.0 / 6.0]
    k = 6
    rows = ["ell,beta,N,overlap_err," + ",".join(f"err_E{i+1}" for i in range(k))]
    best = {}
    for ell in (0, 1):
        kk = k if ell == 0 else k - 1
        ex = exact_hydrogen(6, ell)[:kk]
        for beta in betas:
            errs_by_N = []
            for N in sizes:
                B = BasisSet(N, ell, beta)
                ovl = B.overlap_error()
                lam = np.linalg.eigvalsh(B.hamiltonian_hydrogen())[:kk]
                err = np.abs(lam - ex)
                errs_by_N.append(err)
                rows.append(f"{ell},{beta:.6f},{N},{ovl:.3e}," +
                            ",".join(f"{e:.3e}" for e in err))
            best[(ell, beta)] = np.array(errs_by_N)

    (RES / "E1_hydrogen_convergence.csv").write_text("\n".join(rows), encoding="utf-8")

    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.8))
    for ax, ell in zip(axes, (0, 1)):
        for beta, mark in zip(betas, "osd^"):
            errs = best[(ell, beta)]
            worst = errs.max(axis=1)          # worst state error
            lab = f"$\\beta = {beta:.3g}$"
            ax.semilogy(sizes, np.maximum(worst, 1e-16), marker=mark, ms=4, label=lab)
        ax.set_xlabel("basis size $N$")
        ax.set_title(f"$\\ell = {ell}$")
    axes[0].set_ylabel(r"max$_n$ $|E_n^{\rm num} - E_n^{\rm exact}|$  (Ha)")
    axes[0].legend(fontsize=8)
    fig.suptitle("Hydrogen bound states: convergence of the orthonormal Laguerre basis")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"E1_convergence.{ext}")
    plt.close(fig)

    # headline numbers for the summary
    B = BasisSet(150, 0, 1.0 / 6.0)
    lam = np.linalg.eigvalsh(B.hamiltonian_hydrogen())[:6]
    return lam, np.abs(lam - exact_hydrogen(6, 0))


# ----------------------------------------------------------------- E2 --
def e2_feast_validation():
    print("[E2] FEAST vs LAPACK")
    _self_test()
    B = BasisSet(150, 0, 1.0 / 6.0)
    H = B.hamiltonian_hydrogen()
    lap = np.linalg.eigvalsh(H)
    t0 = time.perf_counter()
    lam, X, res, iters, info = feast_interval(H, -0.6, -0.012, m0=12)
    dt = time.perf_counter() - t0

    inside_lap = lap[(lap > -0.6) & (lap < -0.012)]
    agree = np.max(np.abs(np.sort(lam) - inside_lap))

    lines = [
        "FEAST validation on hydrogen Hamiltonian (Laguerre basis N=150, beta=1/6, ell=0)",
        f"search interval      : [-0.600, -0.012] Ha",
        f"eigenvalues found    : {lam.size} (expected 6: n=1..6)",
        f"FEAST refinement its : {iters}   wall time: {dt*1e3:.1f} ms",
        f"max |FEAST - LAPACK| : {agree:.3e} Ha",
        f"max residual ||Hx-Ex||: {res.max():.3e}",
        "",
        "contour-size sensitivity (n_contour -> max residual, iterations):",
    ]
    for nc in (4, 8, 16, 24):
        l2, _, r2, it2, _ = feast_interval(H, -0.6, -0.012, m0=12, n_contour=nc)
        lines.append(f"  n_contour={nc:2d}: found={l2.size}, res={r2.max():.2e}, its={it2}")
    (RES / "E2_feast_validation.txt").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join("   " + l for l in lines))
    return agree


# ----------------------------------------------------------------- E3 --
def e3_yukawa():
    print("[E3] Yukawa / Debye screening: Laguerre vs finite differences")
    lams = [0.02, 0.05, 0.10, 0.20]
    k = 3
    rows = ["lambda,state,E_laguerre_Ha,E_fd_Ha,abs_diff"]
    table = []
    for lam in lams:
        B = BasisSet(180, 0, 1.0 / 6.0)
        EL = np.linalg.eigvalsh(B.hamiltonian_yukawa(lam=lam))[:k]
        EF = fd.fd_levels(fd.yukawa(lam=lam), ell=0, k=k, dr=0.004, rmax=400.0)
        for i in range(k):
            bound = EL[i] < -1e-6
            rows.append(f"{lam},{i+1},{EL[i]:.10f},{EF[i]:.10f},{abs(EL[i]-EF[i]):.2e}")
            if bound:
                table.append((lam, i + 1, EL[i], EF[i]))
    (RES / "E3_yukawa_crosscheck.csv").write_text("\n".join(rows), encoding="utf-8")

    fig, ax = plt.subplots(figsize=(5.2, 3.8))
    lgrid = np.linspace(0.0, 0.24, 25)
    for state in range(3):
        Es = []
        B = BasisSet(180, 0, 1.0 / 6.0)
        T = B.kinetic()
        for la in lgrid:
            E = np.linalg.eigvalsh(T + B.yukawa(lam=la))[state]
            Es.append(E if E < 0 else np.nan)
        ax.plot(lgrid, Es, marker="", label=f"$n={state+1}$ ($\\ell=0$)")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel(r"screening parameter $\lambda$  (a.u.)")
    ax.set_ylabel(r"$E$  (Ha)")
    ax.set_title("Bound-state energies of the screened Coulomb potential")
    ax.legend(fontsize=8)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"E3_yukawa_levels.{ext}")
    plt.close(fig)
    return table


# ----------------------------------------------------------------- E4 --
def e4_complex_scaling():
    print("[E4] Complex scaling: theta-independence of bound states")
    B = BasisSet(120, 0, 1.0 / 3.0)
    ex = exact_hydrogen(4, 0)
    lines = ["theta,state,Re_E,Im_E,|E - E_exact|"]
    worst = 0.0
    for theta in (0.10, 0.20, 0.30, 0.40):
        H = B.hamiltonian_complex_scaled(theta=theta)
        w = np.linalg.eigvals(H)
        # match the 4 lowest bound states: nearest eigenvalue to each exact E
        for i, Eex in enumerate(ex):
            j = np.argmin(np.abs(w - Eex))
            dev = abs(w[j] - Eex)
            worst = max(worst, dev if i < 3 else 0.0)  # n=4 needs bigger N
            lines.append(f"{theta:.2f},{i+1},{w[j].real:.12f},{w[j].imag:+.3e},{dev:.3e}")
    (RES / "E4_complex_scaling_check.csv").write_text("\n".join(lines), encoding="utf-8")
    return worst


# ---------------------------------------------------------------- main --
if __name__ == "__main__":
    t0 = time.perf_counter()
    lam6, err6 = e1_hydrogen_convergence()
    agree = e2_feast_validation()
    ytable = e3_yukawa()
    worst_cs = e4_complex_scaling()

    print("\n================ SUMMARY ================")
    print("E1 hydrogen ell=0, N=150, beta=1/6 -> errors vs exact (Ha):")
    for i, (E, e) in enumerate(zip(lam6, err6)):
        print(f"   n={i+1}:  E = {E:.12f}   |err| = {e:.2e}")
    print(f"E2 max |FEAST - LAPACK|          : {agree:.2e} Ha")
    print(f"E3 Yukawa bound states validated : {len(ytable)} (see CSV)")
    print(f"E4 complex-scaling worst drift   : {worst_cs:.2e} Ha (n=1..3)")
    print(f"total wall time: {time.perf_counter()-t0:.1f} s")
