"""Paper F, F3f: basis/channel convergence ladder at Z=82 (the honest
error bar of the extreme row).  Same impostor-free protocol as f3e,
at the two corners u=0 and u=0.20, f=0.05, theta=0.08/0.13, for
(nchan, nb) = (12, 52) production -> (16, 58) -> (20, 64).

Deliverable: Gamma at each rung + the Maltsev u=0 comparison, so the
manuscript's "1.4% from Maltsev, to be tightened with a convergence
ladder" sentence can be replaced by measured bars.
"""

import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from ds_core import CoupledDiracStark  # noqa: E402
from angular import kappa_ladder  # noqa: E402

t0 = time.perf_counter()

Z = 82
FSC = 0.05
F = FSC * Z ** 3
dE0, hG0 = -2.571202e1, -2.427590e-2   # Maltsev seed (dE, -Gamma/2)
G_MALTSEV = 2 * 2.427590e-2

for nchan, nb in ((12, 52), (16, 58), (20, 64)):
    kap = kappa_ladder(-1, nchan, 1)
    S = CoupledDiracStark(nb=nb, kappas=kap, Z=Z, beta=0.55 * Z, m2=1,
                          nq=9000)
    print(f"\n[nchan={nchan} nb={nb} dim={2*nb*nchan} cert={S.cert:.0e}]")

    def ref_1s(lam):
        H0 = S.build(0.0, 0.0, lam)
        blk = H0[:2 * nb, :2 * nb]
        ev = np.sort(np.linalg.eigvals(blk).real)
        ev = ev[(ev > -2 * S.c ** 2 + 1e-6) & (ev < -1e-9)]
        return ev[0]

    prevE = prevRef = None
    for u in (0.0, 0.05, 0.10, 0.15, 0.20):
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
            prevE = prevRef = None
            continue
        e = picks[1]
        q = abs(picks[1] - picks[0])
        extra = ""
        if u == 0.0:
            extra = f"  vs Maltsev: {(-2*e.imag - G_MALTSEV)/G_MALTSEV:+.2%}"
        print(f"   u={u:.2f}: dE={e.real-E0:+.7e}  Gamma={-2*e.imag:.6e}"
              f"  (q {q:.0e}){extra}")
        prevE, prevRef = e, E0

print(f"\ntotal: {time.perf_counter()-t0:.1f} s")
