"""Paper F, F3c PRODUCTION: relativistic screened Stark maps.

THE FIX from F3a: the field-free reference must come from the kappa=-1
CHANNEL BLOCK -- the full-matrix lowest in-gap eigenvalue is the
kappa=+1 spurious impostor (exactly degenerate with 1s at lam=0 by the
Tupitsyn-Shabaev theorem, detaching DOWNWARD under screening).  Our
Papers A/T pollution theory predicted this; here it is excluded by
construction and the physical anchor is per-channel.

Scaled units (nonrel invariants): u = lam/Z, f = F/Z^3.  Nonrel physics
collapses onto Z-independent curves E/Z^2, Gamma/Z^2; the Z-dependence
of the DEVIATION from that collapse is pure relativity -- the paper's
headline observable.  Nonrel baseline computed live with the Paper D
engine (StarkChannels) at (f, u).

Outputs: f3c_1s_map.csv  (Z, u, lam, F, E0ref, Er, Gamma, qual,
                          Er_nr, Gamma_nr)   [nonrel scaled to Z=1]
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

FSC = 0.05                       # scaled field F/Z^3 (Maltsev anchor)
US = (0.0, 0.05, 0.10, 0.15, 0.20, 0.30)
ZS = (10, 20, 40)

# ---------- nonrel baseline at Z=1: F=FSC, lambda=u ---------------------
print("Nonrel baseline (Paper D engine, Z=1, F=%.2f):" % FSC)
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
    e = picks[1]
    nr[u] = e
    print(f"   u={u:.2f}: E={e.real:+.7f}{e.imag:+.7f}j "
          f"(q {abs(picks[1]-picks[0]):.0e})")
    prev = e

# ---------- relativistic runs -------------------------------------------
rows = []
for Z in ZS:
    F = FSC * Z ** 3
    kap = kappa_ladder(-1, 12, 1)
    nb = 52
    S = CoupledDiracStark(nb=nb, kappas=kap, Z=Z, beta=0.55 * Z, m2=1,
                          nq=8500)
    print(f"\nZ={Z} F={F:.4g} (cert {S.cert:.0e}):")

    def ref_1s(lam):
        """Physical field-free 1s: kappa=-1 block ONLY (impostor-free)."""
        H0 = S.build(0.0, 0.0, lam)
        blk = H0[:2 * nb, :2 * nb]
        ev = np.sort(np.linalg.eigvals(blk).real)
        ev = ev[(ev > -2 * S.c ** 2 + 1e-6) & (ev < -1e-9)]
        return ev[0]

    prevE = None
    prevRef = None
    for u in US:
        lam = u * Z
        E0 = ref_1s(lam)
        if prevE is None:
            pred = E0 + nr[u].real * 0 + (Z ** 2) * (-6.105e-3) \
                + 1j * (Z ** 2) * (-3.86e-5)
        else:
            pred = E0 + (prevE - prevRef)
        rad = max(0.06 * Z ** 2 * 1e-2, 4 * abs(pred.imag), 0.05)
        picks = []
        for th in (0.10, 0.16) if Z >= 20 else (0.14, 0.22):
            ev = np.linalg.eigvals(S.build(F, th, lam))
            sel = ev[(np.abs(ev - pred) < rad) & (ev.imag < 1e-7)]
            if sel.size:
                picks.append(sel[np.argmin(np.abs(sel - pred))])
        if len(picks) != 2:
            print(f"   u={u:.2f} lam={lam:.2f}: LOST")
            continue
        e = picks[1]
        q = abs(picks[1] - picks[0])
        enr = nr[u]
        rows.append(dict(Z=Z, u=u, lam=lam, F=F, E0ref=E0,
                         Er=e.real, Gamma=-2 * e.imag, qual=q,
                         Er_nr=enr.real, Gamma_nr=-2 * enr.imag))
        supp = (-2 * e.imag) / (Z ** 2 * (-2 * enr.imag)) \
            if enr.imag < -1e-15 else float("nan")
        print(f"   u={u:.2f} lam={lam:6.2f}: dE={e.real-E0:+.5e}  "
              f"Gamma={-2*e.imag:.5e}  supp={supp:.4f}  (q {q:.0e})")
        prevE, prevRef = e, E0

with (HERE / "f3c_1s_map.csv").open("w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    for r in rows:
        w.writerow(r)

print(f"\ntotal: {time.perf_counter()-t0:.1f} s")
