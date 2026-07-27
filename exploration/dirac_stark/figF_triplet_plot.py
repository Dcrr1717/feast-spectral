"""Regenerate figF_triplet (lost plotting script; rebuilt 26-jul) with
the F3g nonrelativistic baseline added to the right panel.

Left: widths (log) of the three field-dressed n=2 components at Z=40,
f=0.004 (f3d_triplet.csv).  Right: fine-structure ratio Gamma_C/Gamma_A
vs the matched nonrelativistic Gamma_up/Gamma_down (f3g_nr_baseline.csv,
dashed; u<=0.15 only — beyond that the nonrel uphill tracking degrades).
"""

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).parent

tri = {}
with (HERE / "f3d_triplet.csv").open() as fh:
    for r in csv.DictReader(fh):
        tri.setdefault(r["comp"], []).append(
            (float(r["u"]), float(r["Er"]), float(r["Gamma"])))

nr = {}
with (HERE / "f3g_nr_baseline.csv").open() as fh:
    for r in csv.DictReader(fh):
        nr.setdefault(r["comp"], []).append(
            (float(r["lam"]), float(r["Gamma"])))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.0))

STYLE = {"A": ("tab:blue", "o", "A (downhill)"),
         "B": ("tab:gray", "s", "B (central)"),
         "C": ("tab:red", "^", r"C (uphill, $2p_{3/2}$)")}
for comp, pts in sorted(tri.items()):
    c, mk, lab = STYLE[comp]
    u = [p[0] for p in pts]
    g = [p[2] for p in pts]
    ax1.semilogy(u, g, marker=mk, ms=4, color=c, lw=1.2, label=lab)
ax1.set_xlabel(r"$u=\lambda/Z$")
ax1.set_ylabel(r"$\Gamma$ (a.u.)")
ax1.legend(fontsize=7, frameon=False)
ax1.set_title(r"$Z=40$, $f=0.004$, $m=1/2$", fontsize=9)

ga = {p[0]: p[2] for p in tri["A"]}
gc = {p[0]: p[2] for p in tri["C"]}
us = sorted(set(ga) & set(gc))
ax2.plot(us, [gc[u] / ga[u] for u in us], marker="o", ms=4,
         color="tab:purple", lw=1.2,
         label=r"Dirac $\Gamma_{\rm C}/\Gamma_{\rm A}$")

gd = {p[0]: p[1] for p in nr["down"]}
gu = {p[0]: p[1] for p in nr["up"]}
uns = sorted(u for u in (set(gd) & set(gu)) if u <= 0.151)
ax2.plot(uns, [gu[u] / gd[u] for u in uns], ls="--", marker="x", ms=4,
         color="tab:green", lw=1.2,
         label=r"nonrel. $\tilde\Gamma_{\rm up}/\tilde\Gamma_{\rm down}$")
ax2.axhline(1.0, color="0.7", lw=0.6, ls=":")
ax2.set_xlabel(r"$u=\lambda/Z$")
ax2.set_ylabel("width ratio")
ax2.legend(fontsize=7, frameon=False)
ax2.set_title("fine-structure asymmetry", fontsize=9)

fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(HERE / "figs" / f"figF_triplet.{ext}", dpi=200)
print("written figF_triplet.pdf/png")
