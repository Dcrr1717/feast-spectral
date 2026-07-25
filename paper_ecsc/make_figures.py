"""Paper E figures: figE_map (ground (F,lam) map), figE_2p (field-free
2p resonance birth), figE_manifold (n=2 doublet)."""

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "exploration" / "ecsc"
FIGS = HERE / "figures"
FIGS.mkdir(exist_ok=True)
sys.path.insert(0, str(DATA))

plt.rcParams.update({
    "font.size": 9, "axes.labelsize": 9, "legend.fontsize": 7.5,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "figure.dpi": 300,
    "lines.linewidth": 1.2, "lines.markersize": 3.2,
})

# ---------------- figE_map ----------------------------------------------
rows = list(csv.DictReader((DATA / "e2_map.csv").open()))
by_F = {}
for r in rows:
    by_F.setdefault(float(r["F"]), []).append(
        (float(r["lam"]), float(r["Er"]), float(r["Gamma"]),
         float(r["Eref"]) if r["Eref"] not in ("", "nan") else np.nan,
         float(r["qual"]) if r["qual"] not in ("", "nan") else np.nan))

fig, (ax1, axm, ax2) = plt.subplots(1, 3, figsize=(6.8, 2.4))
colors = plt.cm.viridis(np.linspace(0.0, 0.85, len(by_F)))
for c, (F, pts) in zip(colors, sorted(by_F.items())):
    pts.sort()
    lam = np.array([p[0] for p in pts])
    er = np.array([p[1] for p in pts])
    ga = np.array([p[2] for p in pts])
    eref = np.array([p[3] for p in pts])
    qq = np.array([p[4] for p in pts])
    m = ~np.isnan(er)
    ax1.plot(lam[m], er[m], "o-", color=c, label=f"$F={F:.2f}$")
    sh = er - eref
    mm = m & ~np.isnan(eref)
    axm.plot(lam[mm], 1e3 * sh[mm], "o-", color=c)
    mg = m & (ga > np.maximum(10 * qq, 1e-10))
    ax2.semilogy(lam[mg], ga[mg], "o-", color=c)
pts = sorted(by_F[min(by_F)])
lam = np.array([p[0] for p in pts])
eref = np.array([p[3] for p in pts])
mm = ~np.isnan(eref)
ax1.plot(lam[mm], eref[mm], "k--", lw=0.9, label=r"$F=0$ (bound)")
ax1.axhline(0.0, color="k", lw=0.6, ls=":")
ax1.set_xlabel(r"$\lambda$ (a.u.)")
ax1.set_ylabel(r"$E_r$ (a.u.)")
ax1.legend(frameon=False, ncol=2, columnspacing=0.6, handlelength=1.2,
           fontsize=6.6)
axm.set_xlabel(r"$\lambda$ (a.u.)")
axm.set_ylabel(r"$E_r - E_b^{F=0}$ ($10^{-3}$ a.u.)")
ax2.set_xlabel(r"$\lambda$ (a.u.)")
ax2.set_ylabel(r"$\Gamma$ (a.u.)")
for ax in (ax1, axm, ax2):
    ax.tick_params(direction="in", top=True, right=True)
fig.tight_layout(pad=0.4)
fig.savefig(FIGS / "figE_map.pdf")
print("figE_map.pdf")

# ---------------- figE_2p: field-free 2p resonance birth ----------------
from ecsc_core import StarkChannelsE  # noqa: E402

S1 = StarkChannelsE(200, 1, 0.15)
nb = S1.nb


def bound_l1(lam):
    H = S1.build(1.0, lam, 0.0, 0.0, scr="ecsc")
    H11 = H[nb:, nb:]
    return float(np.linalg.eigvalsh(0.5 * (H11.real + H11.real.T))[0])


def nearest_l1(lam, th, target):
    H = S1.build(1.0, lam, 0.0, th, scr="ecsc")
    ev = np.linalg.eigvals(H[nb:, nb:])
    return ev[np.argmin(np.abs(ev - target))]


lamc = 0.148205
lams_b = np.linspace(0.05, lamc, 12)
eb = [bound_l1(l) for l in lams_b]
prev = complex(eb[-1])
prev_lam = lams_b[-1]
slope = (eb[-1] - eb[-2]) / (lams_b[-1] - lams_b[-2])
traj = []
for lam in np.round(np.arange(lamc + 0.004, 0.2401, 0.004), 4):
    pred = prev + slope * (lam - prev_lam)
    picks = [nearest_l1(lam, th, pred) for th in (0.30, 0.50)]
    e, q = picks[1], abs(picks[1] - picks[0])
    traj.append((lam, e, q))
    slope = (e - prev) / (lam - prev_lam)
    prev, prev_lam = e, lam

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.8, 2.5))
ax1.plot(lams_b, eb, "o-", color="tab:gray", label="bound $2p$")
lamr = np.array([t[0] for t in traj])
er = np.array([t[1].real for t in traj])
ga = np.array([-2 * t[1].imag for t in traj])
qq = np.array([t[2] for t in traj])
ax1.plot(lamr, er, "o-", color="tab:red", label=r"resonance $E_r$")
ax1.fill_between(lamr, er - ga / 2, er + ga / 2, color="tab:red",
                 alpha=0.18, lw=0)
ax1.axhline(0.0, color="k", lw=0.6, ls=":")
ax1.axvline(lamc, color="k", lw=0.6, ls="--")
ax1.set_xlabel(r"$\lambda$ (a.u.)")
ax1.set_ylabel(r"$E$ (a.u.)")
ax1.legend(frameon=False)
mres = ga > np.maximum(10 * qq, 1e-12)
ax2.semilogy(lamr[mres], ga[mres], "o-", color="tab:red")
ax2.set_xlabel(r"$\lambda$ (a.u.)")
ax2.set_ylabel(r"$\Gamma$ (a.u.)")
for ax in (ax1, ax2):
    ax.tick_params(direction="in", top=True, right=True)
fig.tight_layout(pad=0.4)
fig.savefig(FIGS / "figE_2p.pdf")
print("figE_2p.pdf")

# ---------------- figE_manifold -----------------------------------------
rows = list(csv.DictReader((DATA / "e2_manifold.csv").open()))
comps = {}
for r in rows:
    comps.setdefault(r["comp"], []).append(
        (float(r["lam"]), float(r["Er"]), float(r["Gamma"])))
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.8, 2.5))
style = {"red": ("tab:red", "$k=-1$"), "blue": ("tab:blue", "$k=+1$")}
for name, pts in comps.items():
    pts.sort()
    lam = np.array([p[0] for p in pts])
    er = np.array([p[1] for p in pts])
    ga = np.array([p[2] for p in pts])
    c, lab = style.get(name, ("k", name))
    ax1.plot(lam, er, "o-", color=c, label=lab)
    ax1.fill_between(lam, er - ga / 2, er + ga / 2, color=c, alpha=0.18,
                     lw=0)
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
fig.savefig(FIGS / "figE_manifold.pdf")
print("figE_manifold.pdf")
