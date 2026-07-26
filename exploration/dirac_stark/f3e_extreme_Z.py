"""Paper F, F3e: extend the suppression map to Z=54, 82 (extreme
relativistic regime -- Maltsev's field-free suppression reaches ~0.09
at Z=82).  Same impostor-free protocol as f3c; Maltsev values as u=0
seeds.  Appends to f3c_1s_map.csv format in its own file.
"""

import csv
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "stark2d"))

from ds_core import CoupledDiracStark  # noqa: E402
from angular import kappa_ladder  # noqa: E402
from d0_architecture import StarkChannels  # noqa: E402

t0 = time.perf_counter()

FSC = 0.05
US = (0.0, 0.05, 0.10, 0.15, 0.20)
SEEDS = {54: (-1.478389e1, -4.599133e-2),   # Maltsev dE, -Gamma/2
         82: (-2.571202e1, -2.427590e-2)}

print("Nonrel baseline (Z=1, F=%.2f):" % FSC)
NR = StarkChannels(70, 14, 0.35)
nr = {}
prev = None
for u in US:
    pred = (-0.52 + 0.0j) if prev is None else prev
    picks = []
    for th in (0.18, 0.28):
        ev = np.linalg.eigvals(NR.build(1.0, u, FSC, th))
        sel = ev[(np.abs(ev - pred) < 0.08) & (ev.imag < 1e-8)]
        if sel.size:
            picks.append(sel[np.argmin(np.abs(sel - pred))])
    nr[u] = picks[1]
    prev = picks[1]
    print(f"   u={u:.2f}: {nr[u].real:+.7f}{nr[u].imag:+.7f}j")

rows = []
for Z, (dE0, hG0) in SEEDS.items():
    F = FSC * Z ** 3
    kap = kappa_ladder(-1, 12, 1)
    nb = 52
    S = CoupledDiracStark(nb=nb, kappas=kap, Z=Z, beta=0.55 * Z, m2=1,
                          nq=9000)
    print(f"\nZ={Z} F={F:.4g} (cert {S.cert:.0e}, "
          f"gamma1={S.gamma[0]:.4f}):")

    def ref_1s(lam):
        H0 = S.build(0.0, 0.0, lam)
        blk = H0[:2 * nb, :2 * nb]
        ev = np.sort(np.linalg.eigvals(blk).real)
        ev = ev[(ev > -2 * S.c ** 2 + 1e-6) & (ev < -1e-9)]
        return ev[0]

    prevE = prevRef = None
    for u in US:
        lam = u * Z
        E0 = ref_1s(lam)
        if prevE is None:
            pred = E0 + dE0 + 1j * hG0
        else:
            pred = E0 + (prevE - prevRef)
        rad = max(0.10 * abs(dE0), 5 * abs(pred.imag))
        picks = []
        for th in (0.08, 0.13):
            ev = np.linalg.eigvals(S.build(F, th, lam))
            sel = ev[(np.abs(ev - pred) < rad) & (ev.imag < 1e-7)]
            if sel.size:
                picks.append(sel[np.argmin(np.abs(sel - pred))])
        if len(picks) != 2:
            print(f"   u={u:.2f}: LOST")
            continue
        e = picks[1]
        q = abs(picks[1] - picks[0])
        enr = nr[u]
        supp = (-2 * e.imag) / (Z ** 2 * (-2 * enr.imag))
        rows.append(dict(Z=Z, u=u, lam=lam, F=F, E0ref=E0, Er=e.real,
                         Gamma=-2 * e.imag, qual=q, Er_nr=enr.real,
                         Gamma_nr=-2 * enr.imag))
        print(f"   u={u:.2f} lam={lam:6.2f}: dE={e.real-E0:+.5e}  "
              f"Gamma={-2*e.imag:.5e}  supp={supp:.4f}  (q {q:.0e})")
        prevE, prevRef = e, E0

with (HERE / "f3e_extremeZ.csv").open("w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    for r in rows:
        w.writerow(r)

print(f"\ntotal: {time.perf_counter()-t0:.1f} s")
