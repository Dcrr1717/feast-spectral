"""Paper F, F3d v2: the n=2 (m=1/2) Stark QUARTET under screening at
Z=40 -- fine-structure width splitting via proper continuation.

The m=1/2, n=2 manifold has FOUR states (2s1/2, 2p1/2 from kappa=+-1;
2p3/2 from kappa=-2), which the field mixes into Stark components with
linear shifts of order +-3*f2*Z^2 (~ +-19 a.u. at Z=40, f2=0.004!).
Naive window seeds fail (f3d v1 lesson).  Protocol here:

1. F-RAMP at lam=0: start each component at its field-free level
   (per-channel block references, impostor-free) and ramp the field in
   geometric steps 0 -> F2, tracking by prediction+pairing at each step
   (linear extrapolation in F).  This labels the quartet by adiabatic
   continuation -- the only defensible labeling at strong mixing.
2. LAM-CONTINUATION at full F2: from each ramped component, step u in
   {0, 0.05, 0.10, 0.15, 0.20} using reference-anchored prediction
   E0_channel(lam) + previous field shift.

Output: f3d_quartet.csv (u, lam, comp, Er, Gamma, qual)
"""

import csv
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from ds_core import CoupledDiracStark, channel_sommerfeld  # noqa: E402
from angular import kappa_ladder  # noqa: E402

t0 = time.perf_counter()

Z = 40
F2 = 0.004
F = F2 * Z ** 3
US = (0.0, 0.05, 0.10, 0.15, 0.20)
kap = kappa_ladder(-1, 12, 1)
nb = 52
S = CoupledDiracStark(nb=nb, kappas=kap, Z=Z, beta=0.55 * Z, m2=1,
                      nq=8500)
print(f"Z={Z} F={F:.4g} (f2={F2}, cert {S.cert:.0e})")

COMPS = {
    "2s1/2": (kap.index(-1), channel_sommerfeld(-1, Z)[2]),
    "2p1/2": (kap.index(+1), channel_sommerfeld(+1, Z)[2]),
    "2p3/2": (kap.index(-2), channel_sommerfeld(-2, Z)[2]),
}
print("field-free:", {k: f"{v[1]:.4f}" for k, v in COMPS.items()})


def ref_channel(ich, lam, guess):
    H0 = S.build(0.0, 0.0, lam)
    blk = H0[2 * nb * ich:2 * nb * (ich + 1),
             2 * nb * ich:2 * nb * (ich + 1)]
    ev = np.sort(np.linalg.eigvals(blk).real)
    ev = ev[(ev > -2 * S.c ** 2 + 1e-6) & (ev < -1e-9)]
    return ev[np.argmin(np.abs(ev - guess))]


def track(Fval, lam, pred, rad, tol, thetas=(0.10, 0.16)):
    ev = np.linalg.eigvals(S.build(Fval, thetas[0], lam))
    sel = ev[(np.abs(ev - pred) < rad) & (ev.imag < 1e-7)]
    if not sel.size:
        return None
    e0 = sel[np.argmin(np.abs(sel - pred))]
    ev = np.linalg.eigvals(S.build(Fval, thetas[1], lam))
    sel = ev[(np.abs(ev - e0) < rad) & (ev.imag < 1e-7)]
    if not sel.size:
        return None
    e1 = sel[np.argmin(np.abs(sel - e0))]
    d = abs(e1 - e0)
    return (e1, d) if d <= tol else None


# ---------- 1. F-ramp at lam=0 ------------------------------------------
print("\nF-ramp at lam=0 (adiabatic labeling):")
RAMP = [F * s for s in (0.05, 0.15, 0.3, 0.5, 0.75, 1.0)]
state = {}
for name, (ich, e2) in COMPS.items():
    prevF, prevE = 0.0, complex(e2)
    slope = 0.0 + 0.0j
    ok = True
    for Fv in RAMP:
        pred = prevE + slope * (Fv - prevF)
        rad = max(1.0, 0.15 * abs(Fv - prevF) * 3 * F2 * Z ** 2 / F,
                  6 * abs(pred.imag))
        rad = max(rad, 0.08 * (Fv - prevF) / F * 3 * 0.004 * Z * Z)
        got = track(Fv, 0.0, pred, rad, tol=0.5)
        if got is None:
            print(f"   {name}: ramp LOST at F={Fv:.3g}")
            ok = False
            break
        e, q = got
        slope = (e - prevE) / (Fv - prevF)
        prevF, prevE = Fv, e
    if ok:
        state[name] = prevE
        print(f"   {name}: E(F2)={prevE.real:+.5f}{prevE.imag:+.5f}j "
              f"(dE={prevE.real-e2:+.4f}, q {q:.0e})")

# ---------- 2. lam-continuation at full field ---------------------------
print("\nlam-continuation at F2:")
rows = []
prevRef = {n: COMPS[n][1] for n in state}
prevE = dict(state)
guess = {n: COMPS[n][1] for n in state}
for u in US:
    lam = u * Z
    for name in list(state):
        ich = COMPS[name][0]
        E0 = ref_channel(ich, lam, guess[name])
        guess[name] = E0
        if u == 0.0:
            e = prevE[name]
            q = 0.0
        else:
            pred = E0 + (prevE[name] - prevRef[name])
            rad = max(1.0, 6 * abs(pred.imag))
            got = track(F, lam, pred, rad, tol=0.6)
            if got is None:
                print(f"   {name} u={u:.2f}: LOST")
                continue
            e, q = got
        rows.append(dict(u=u, lam=lam, comp=name, Er=e.real,
                         Gamma=-2 * e.imag, qual=q))
        print(f"   {name} u={u:.2f}: dE={e.real-E0:+.5e}  "
              f"Gamma={-2*e.imag:.5e}  (q {q:.0e})")
        prevE[name], prevRef[name] = e, E0

with (HERE / "f3d_quartet.csv").open("w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    for r in rows:
        w.writerow(r)

print("\nGamma(2p3/2)/Gamma(2p1/2) vs u:")
by = {}
for r in rows:
    by.setdefault(r["u"], {})[r["comp"]] = r["Gamma"]
for u in US:
    d = by.get(u, {})
    if "2p1/2" in d and "2p3/2" in d and d["2p1/2"] > 1e-14:
        print(f"   u={u:.2f}: {d['2p3/2']/d['2p1/2']:.4f}")

print(f"\ntotal: {time.perf_counter()-t0:.1f} s")
