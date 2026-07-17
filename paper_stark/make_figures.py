"""Paper D figures from d3_map.csv / d3_manifold.csv."""

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "exploration" / "stark2d"
FIGS = HERE / "figures"
FIGS.mkdir(exist_ok=True)

plt.rcParams.update({
    "font.size": 9, "axes.labelsize": 9, "legend.fontsize": 7.5,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "figure.dpi": 300,
    "lines.linewidth": 1.2, "lines.markersize": 3.5,
})

# ---------------- fig D1: the (F, lambda) ground map --------------------
rows = list(csv.DictReader((DATA / "d3_map.csv").open()))
by_F = {}
for r in rows:
    by_F.setdefault(float(r["F"]), []).append(
        (float(r["lam"]), float(r["Er"]), float(r["Gamma"]),
         float(r["Eyuk"])))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.8, 2.7))
colors = plt.cm.viridis(np.linspace(0.0, 0.85, len(by_F)))
for c, (F, pts) in zip(colors, sorted(by_F.items())):
    pts.sort()
    lam = np.array([p[0] for p in pts])
    er = np.array([p[1] for p in pts])
    ga = np.array([p[2] for p in pts])
    m = ~np.isnan(er)
    ax1.plot(lam[m], er[m], "o-", color=c, label=f"$F={F:.2f}$")
    ax2.semilogy(lam[m], np.maximum(ga[m], 1e-12), "o-", color=c)
# Yukawa reference (F=0) from any row set
pts = sorted(by_F[min(by_F)])
lam = np.array([p[0] for p in pts])
ey = np.array([p[3] for p in pts])
ax1.plot(lam, ey, "k--", lw=0.9, label=r"$F=0$ (bound)")
ax1.set_xlabel(r"$\lambda$ (a.u.)")
ax1.set_ylabel(r"$E_r$ (a.u.)")
ax1.legend(frameon=False, ncol=2, columnspacing=0.8, handlelength=1.4)
ax2.set_xlabel(r"$\lambda$ (a.u.)")
ax2.set_ylabel(r"$\Gamma$ (a.u.)")
for ax in (ax1, ax2):
    ax.tick_params(direction="in", top=True, right=True)
fig.tight_layout(pad=0.4)
fig.savefig(FIGS / "figD_map.pdf")
print("figD_map.pdf", FIGS / "figD_map.pdf")

# ---------------- fig D2: n=2 doublet under screening -------------------
rows = list(csv.DictReader((DATA / "d3_manifold.csv").open()))
comps = {}
for r in rows:
    comps.setdefault(r["comp"], []).append(
        (float(r["lam"]), float(r["Er"]), float(r["Gamma"])))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.8, 2.7))
style = {"red": ("tab:red", "$k=-1$ (red)"), "blue": ("tab:blue", "$k=+1$ (blue)")}
for name, pts in comps.items():
    pts.sort()
    lam = np.array([p[0] for p in pts])
    er = np.array([p[1] for p in pts])
    ga = np.array([p[2] for p in pts])
    c, lab = style.get(name, ("k", name))
    ax1.plot(lam, er, "o-", color=c, label=lab)
    ax1.fill_between(lam, er - ga / 2, er + ga / 2, color=c, alpha=0.18, lw=0)
    ax2.semilogy(lam, np.maximum(ga, 1e-12), "o-", color=c, label=lab)
ax1.axhline(0.0, color="k", lw=0.6, ls=":")
ax1.set_xlabel(r"$\lambda$ (a.u.)")
ax1.set_ylabel(r"$E_r$ (a.u.)")
ax1.legend(frameon=False)
ax2.set_xlabel(r"$\lambda$ (a.u.)")
ax2.set_ylabel(r"$\Gamma$ (a.u.)")
for ax in (ax1, ax2):
    ax.tick_params(direction="in", top=True, right=True)
fig.tight_layout(pad=0.4)
fig.savefig(FIGS / "figD_manifold.pdf")
print("figD_manifold.pdf", FIGS / "figD_manifold.pdf")
