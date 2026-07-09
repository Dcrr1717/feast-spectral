"""Regenerate all manuscript figures in publication style.

Journal (Bull CompAMa) requirements and professor's feedback applied:
- no titles embedded in the images (captions carry that information)
- consistent font sizes (11 pt, matching 12 pt body text)
- legible in black and white: distinct markers and line styles, and a
  labeled grayscale-safe colorbar instead of color-coded legends
- all text in English

Reuses stored results (CSV / NPZ) where available, so no expensive
recomputation is needed except the cheap E(lambda) curve.
"""

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from basis import BasisSet  # noqa: E402

RES, FIG = ROOT / "results", ROOT / "figures"

plt.rcParams.update({
    "font.family": "serif", "font.size": 11,
    "axes.labelsize": 11, "xtick.labelsize": 10, "ytick.labelsize": 10,
    "legend.fontsize": 9.5, "axes.grid": True, "grid.alpha": 0.3,
    "figure.dpi": 150,
})
MARKERS = ["o", "s", "d", "^"]
LSTYLES = ["-", "--", "-.", ":"]


# ------------------------------------------------------------- Figure 1 --
def fig_convergence():
    rows = [l.split(",") for l in
            (RES / "E1_hydrogen_convergence.csv").read_text().splitlines()[1:]]
    data = {}
    for r in rows:
        ell, beta, N = int(r[0]), float(r[1]), int(r[2])
        errs = np.array([float(x) for x in r[4:]])
        data.setdefault((ell, beta), []).append((N, errs.max()))

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.7))
    for ax, ell in zip(axes, (0, 1)):
        betas = sorted({b for (l, b) in data if l == ell}, reverse=True)
        for b, mk, ls in zip(betas, MARKERS, LSTYLES):
            pts = sorted(data[(ell, b)])
            ax.semilogy([p[0] for p in pts],
                        [max(p[1], 1e-16) for p in pts],
                        marker=mk, ls=ls, ms=4.5, lw=1.2,
                        label=f"$\\beta = {b:.3g}$")
        ax.set_xlabel("basis size $N$")
        ax.text(0.05, 0.06, f"$\\ell = {ell}$", transform=ax.transAxes,
                fontsize=12)
    axes[0].set_ylabel(r"$\max_n |E_n^{\rm num} - E_n^{\rm exact}|$  (Ha)")
    axes[0].legend()
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"fig1_convergence.{ext}")
    plt.close(fig)
    print("fig1_convergence done")


# ------------------------------------------------------------- Figure 2 --
def fig_yukawa_levels():
    B = BasisSet(180, 0, 1.0 / 6.0)
    T = B.kinetic()
    lgrid = np.linspace(0.0, 0.24, 25)
    fig, ax = plt.subplots(figsize=(5.4, 3.8))
    for state, mk, ls in zip(range(3), MARKERS, LSTYLES):
        Es = []
        for la in lgrid:
            E = np.linalg.eigvalsh(T + B.yukawa(lam=la))[state]
            Es.append(E if E < 0 else np.nan)
        ax.plot(lgrid, Es, marker=mk, ls=ls, ms=3.5, lw=1.2,
                label=f"${state+1}\\mathrm{{s}}$")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel(r"screening parameter $\lambda$  (a.u.)")
    ax.set_ylabel(r"$E$  (Ha)")
    ax.legend()
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"fig2_yukawa_levels.{ext}")
    plt.close(fig)
    print("fig2_yukawa_levels done")


# ------------------------------------------------------------- Figure 3 --
def fig_trajectory():
    d = np.load(RES / "E6_spectra.npz")
    thetas, traj, Eres = d["thetas"], d["traj"], complex(d["Eres"])
    spectra = [d[f"spec_{i}"] for i in range(len(thetas))]

    fig, ax = plt.subplots(figsize=(6.2, 4.3))
    cmap = plt.cm.viridis
    norm = plt.Normalize(thetas[0], thetas[-1])
    for th, spec in zip(thetas, spectra):
        ax.plot(spec.real, spec.imag, ".", ms=3.2, color=cmap(norm(th)),
                alpha=0.65)
    ax.plot(Eres.real, Eres.imag, "k*", ms=13, label="resonance (stationary)")
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cb = fig.colorbar(sm, ax=ax, pad=0.02)
    cb.set_label(r"rotation angle $\theta$")
    ax.set_xlim(-0.25, 0.25)
    ax.set_ylim(-0.35, 0.02)
    ax.set_xlabel(r"$\mathrm{Re}\,E$  (Ha)")
    ax.set_ylabel(r"$\mathrm{Im}\,E$  (Ha)")
    ax.legend(loc="lower left")

    axin = ax.inset_axes([0.56, 0.55, 0.42, 0.42])
    for z, th in zip(traj, thetas):
        axin.plot(z.real, z.imag, "o", ms=4.5, color=cmap(norm(th)))
    axin.plot(traj.real, traj.imag, "-", lw=0.8, color="0.4", zorder=0)
    axin.plot(Eres.real, Eres.imag, "k*", ms=13)
    pad = 6.0 * max(abs(np.diff(traj)).max(), 1e-4)
    axin.set_xlim(Eres.real - pad, Eres.real + pad)
    axin.set_ylim(Eres.imag - pad, Eres.imag + pad)
    axin.tick_params(labelsize=7)
    ax.indicate_inset_zoom(axin, edgecolor="0.3")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"fig3_trajectory.{ext}")
    plt.close(fig)
    print("fig3_trajectory done")


if __name__ == "__main__":
    fig_convergence()
    fig_yukawa_levels()
    fig_trajectory()
