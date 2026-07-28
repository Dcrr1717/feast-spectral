"""Regenerate all manuscript figures at their FINAL physical size.

The Bull CompAMa text block is 390 pt = 5.40 in wide (measured from the
journal's paper class). Each figure is generated at exactly the width used
by its \\includegraphics call, so a 10 pt font in matplotlib is a true
10 pt on the printed page -- no more shrunken labels.

    fig1_convergence  -> 5.40 in wide (included at \\linewidth)
    fig3_trajectory   -> 5.40 in wide (included at \\linewidth)
    fig4_sweep        -> 5.40 in wide (included at \\linewidth)

Revision 27-jul (meeting feedback): every figure at full \\linewidth,
font sizes raised one step for legibility, and the tight-bbox padding
increased so no axis label is clipped at the figure edge.

Also applied: no embedded titles (captions carry them), black-and-white
friendly markers/line styles, all text in English. Data is read from the
stored results (CSV/NPZ); only the cheap E(lambda) curve would need
recomputation and it is no longer part of the manuscript.
"""

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

RES, FIG = ROOT / "results", ROOT / "figures"

TEXTWIDTH_IN = 390.0 / 72.27          # journal text block width
# Each figure is generated at exactly the width of its \includegraphics
# call, so matplotlib font sizes survive to the printed page unscaled.
W_CONV = 0.85 * TEXTWIDTH_IN          # convergence pair, included at 0.85
W_WIDE = 0.90 * TEXTWIDTH_IN          # Figs. 3 and 4, included at 0.90

plt.rcParams.update({
    # dejavuserif keeps math glyphs in the same (serif) family as the
    # text labels; the default mathtext set is sans and mixes families
    "font.family": "serif", "mathtext.fontset": "dejavuserif",
    "font.size": 11,
    "axes.labelsize": 11.5, "xtick.labelsize": 10.5, "ytick.labelsize": 10.5,
    "legend.fontsize": 10, "axes.grid": True, "grid.alpha": 0.3,
    "figure.dpi": 200, "savefig.bbox": "tight", "savefig.pad_inches": 0.06,
})
MARKERS = ["o", "s", "d", "^"]
LSTYLES = ["-", "--", "-.", ":"]


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"{name}.{ext}")
    plt.close(fig)
    print(name, "done")


# ------------------------------------------------------------- Figure 1 --
def fig_convergence():
    rows = [l.split(",") for l in
            (RES / "E1_hydrogen_convergence.csv").read_text().splitlines()[1:]]
    data = {}
    for r in rows:
        ell, beta, N = int(r[0]), float(r[1]), int(r[2])
        errs = np.array([float(x) for x in r[4:]])
        data.setdefault((ell, beta), []).append((N, errs.max()))

    # Two separate figures (one per partial wave), not a joined pair.
    for ell, name in ((0, "fig1a_convergence_l0"),
                      (1, "fig1b_convergence_l1")):
        fig, ax = plt.subplots(figsize=(W_CONV, 1.60),
                               constrained_layout=True)
        betas = sorted({b for (l, b) in data if l == ell}, reverse=True)
        for b, mk, ls in zip(betas, MARKERS, LSTYLES):
            pts = sorted(data[(ell, b)])
            ax.semilogy([p[0] for p in pts],
                        [max(p[1], 1e-16) for p in pts],
                        marker=mk, ls=ls, ms=5, lw=1.3,
                        label=f"$\\beta = {b:.3g}$")
        ax.set_xlabel("basis size $N$")
        ax.set_ylabel(r"$\max_n |\Delta E_n|$ (Ha)")
        ax.legend(handlelength=1.6, labelspacing=0.25)
        save(fig, name)


# ------------------------------------------------------------- Figure 3 --
def fig_trajectory():
    d = np.load(RES / "E6_spectra.npz")
    thetas, traj, Eres = d["thetas"], d["traj"], complex(d["Eres"])
    spectra = [d[f"spec_{i}"] for i in range(len(thetas))]

    fig, ax = plt.subplots(figsize=(W_WIDE, 1.75))
    cmap = plt.cm.viridis
    norm = plt.Normalize(thetas[0], thetas[-1])
    for th, spec in zip(thetas, spectra):
        ax.plot(spec.real, spec.imag, ".", ms=3.4, color=cmap(norm(th)),
                alpha=0.65)
    ax.plot(Eres.real, Eres.imag, "k*", ms=14, label="resonance")
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cb = fig.colorbar(sm, ax=ax, pad=0.02, fraction=0.055)
    cb.set_label(r"$\theta$", fontsize=11.5)
    cb.ax.tick_params(labelsize=10)
    ax.set_xlim(-0.25, 0.25)
    ax.set_ylim(-0.35, 0.02)
    ax.set_xlabel(r"$\mathrm{Re}\,E$ (Ha)")
    ax.set_ylabel(r"$\mathrm{Im}\,E$ (Ha)")
    ax.legend(loc="lower left", handlelength=1.0)

    # upper-left quadrant is free of data (fans rotate down-right, bound
    # states sit exactly on Im E = 0); keeping the inset there stops its
    # frame and outside tick labels from covering trajectory points
    axin = ax.inset_axes([0.16, 0.52, 0.38, 0.42])
    for z, th in zip(traj, thetas):
        axin.plot(z.real, z.imag, "o", ms=4, color=cmap(norm(th)))
    axin.plot(traj.real, traj.imag, "-", lw=0.8, color="0.4", zorder=0)
    axin.plot(Eres.real, Eres.imag, "k*", ms=11)
    pad = 6.0 * max(abs(np.diff(traj)).max(), 1e-4)
    axin.set_xlim(Eres.real - pad, Eres.real + pad)
    axin.set_ylim(Eres.imag - pad, Eres.imag + pad)
    axin.tick_params(labelsize=8.5)
    ax.indicate_inset_zoom(axin, edgecolor="0.3")
    fig.tight_layout()
    save(fig, "fig3_trajectory")


# ------------------------------------------------------------- Figure 4 --
def fig_sweep():
    rows = [l.split(",") for l in
            (RES / "E8_sweep.csv").read_text().splitlines()[1:]]
    lam = np.array([float(r[0]) for r in rows])
    Er = np.array([float(r[4]) for r in rows])
    Ga = np.array([float(r[5]) for r in rows])
    LAM_ANCHOR = 0.22131705
    LC = 0.220216806605

    fig, axes = plt.subplots(1, 2, figsize=(W_WIDE, 1.65),
                             constrained_layout=True)
    axes[0].plot(lam, Er, "o-", ms=5, lw=1.3, color="k")
    axes[0].set_ylabel(r"$E_r$ (Ha)")
    axes[1].semilogy(lam, Ga, "s--", ms=5, lw=1.3, color="k")
    axes[1].set_ylabel(r"$\Gamma$ (Ha)")
    j = int(np.argmin(np.abs(lam - LAM_ANCHOR)))
    for ax, y in ((axes[0], Er), (axes[1], Ga)):
        ax.axvline(LC, color="0.5", lw=0.9, ls=":")
        # white fill + top zorder so the star stays readable where it
        # lands on top of the filled data markers
        ax.plot(lam[j], y[j], "*", ms=14, mfc="w", mec="k", zorder=5,
                label="benchmark")
    # one shared x-label centered under both panels: nothing can clip
    fig.supxlabel(r"screening parameter $\lambda$ ($a_0^{-1}$)", fontsize=11.5)
    # keep the legend box fully inside the axes ("best" let it straddle
    # the left spine); lower right is the only region the curve never
    # crosses
    axes[0].legend(handlelength=1.0, loc="lower right")
    save(fig, "fig4_sweep")


if __name__ == "__main__":
    fig_convergence()
    fig_trajectory()
    fig_sweep()
