"""Methodology pipeline figure (interim matplotlib version, vector output).

Okabe-Ito colorblind-safe palette; English labels; PDF + PNG.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)

# Okabe-Ito
BLUE, ORANGE, GREEN, VERM, PURPLE = "#56B4E9", "#E69F00", "#009E73", "#D55E00", "#CC79A7"

fig, ax = plt.subplots(figsize=(9.6, 5.2))
ax.set_xlim(0, 10)
ax.set_ylim(0, 6)
ax.axis("off")


def box(x, y, w, h, text, fc, fs=8.6):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08",
                                fc=fc, ec="k", lw=1.0, alpha=0.92))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs)


def arrow(x1, y1, x2, y2, text="", style="-|>"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style,
                                 mutation_scale=14, lw=1.2, color="k"))
    if text:
        ax.text((x1 + x2) / 2 + 0.06, (y1 + y2) / 2, text, fontsize=7.6,
                ha="left", va="center", style="italic")


# --- column 1: formulation
box(0.25, 4.6, 2.5, 1.1,
    "Radial Schrödinger problem\n$-u''/2 + [\\frac{\\ell(\\ell+1)}{2r^2}+V(r)]\\,u = E\\,u$",
    "#EEEEEE", fs=8.2)
box(0.25, 2.9, 2.5, 1.2,
    "Orthonormal Laguerre basis\n$\\varphi_j \\propto x^{\\ell+1}e^{-x/2}L_j^{(2\\ell+2)}(x)$\n"
    "$\\langle\\varphi_i|\\varphi_j\\rangle=\\delta_{ij}$ (certified)",
    BLUE)
box(0.25, 1.2, 2.5, 1.1,
    "Dense $H$ matrix\nanalytic/quadrature elements\n(Gauss–Legendre, self-checked)",
    BLUE)
arrow(1.5, 4.6, 1.5, 4.15)
arrow(1.5, 2.9, 1.5, 2.35)

# --- column 2: solver
box(3.6, 2.9, 2.6, 1.9,
    "FEAST eigensolver\ncontour integration\n$P=\\frac{1}{2\\pi i}\\oint_C (z-H)^{-1}dz$\n"
    "+ Rayleigh–Ritz refinement",
    ORANGE)
arrow(2.75, 1.75, 3.9, 2.85)
arrow(2.75, 3.5, 3.6, 3.7)

# --- column 3: validation & applications
box(7.0, 4.6, 2.7, 1.1,
    "Validation I — hydrogen\nexact $E_n=-1/(2n^2)$\nerror $\\sim 10^{-10}$ Ha",
    GREEN)
box(7.0, 3.05, 2.7, 1.1,
    "Validation II — cross-method\nindependent FD solver\n+ LAPACK reference",
    GREEN)
box(7.0, 1.5, 2.7, 1.1,
    "Application — screened Coulomb\n$V=-e^{-\\lambda r}/r$ (Debye plasma)\nbound states $E(\\lambda)$",
    PURPLE)
box(7.0, 0.1, 2.7, 1.05,
    "Resonances — complex scaling\n$r\\to re^{i\\theta}$, FEAST in $\\mathbb{C}$\npositions + widths",
    VERM)
arrow(6.2, 4.1, 7.0, 4.9)
arrow(6.2, 3.8, 7.0, 3.6)
arrow(6.2, 3.4, 7.0, 2.1)
arrow(6.2, 3.1, 7.0, 0.75)

ax.set_title("Computational pipeline: orthonormal Laguerre basis + FEAST "
             "for bound states and resonances", fontsize=11)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(FIG / f"pipeline_methodology.{ext}", dpi=250)
print("saved", FIG / "pipeline_methodology.pdf")
