"""Paper F, F3d: fine-structure WIDTH splitting of Stark resonances
under screening -- the relativistic echo of Paper B.

At Z=40 the 2p1/2-2p3/2 splitting is large (~4.3 a.u. field-free).
With a (scaled) field on, each fine-structure component becomes a Stark
resonance with its own width; screening then drives the pair toward
threshold at different rates.  Nonrelativistically the two are
degenerate -- ANY difference in Gamma(2p1/2) vs Gamma(2p3/2) is pure
relativity, and Paper B showed (field-free, near lambda_c) that such
splittings can be threshold-amplified by orders of magnitude.

Scaled field f2 = F/Z^3 chosen mild (0.004) so the n=2 pair is a
resonance but not over-broadened: nonrel n=2 critical field ~ 1/(16*16)
= 0.0039, so f2 ~ 1.0x critical -- widths sizable but trackable.

References per component from the per-channel blocks (impostor-free):
2p1/2 = kappa=+1 block, second in-gap state ABOVE the impostor -- the
impostor sits at the (screened) 1s-position track, so select the
in-gap eigenvalue nearest the exact field-free Dirac 2p1/2; 2p3/2 =
kappa=-2 block lowest.

Output: f3d_2p_map.csv (Z, u, lam, F, comp, E0ref, Er, Gamma, qual)
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

# channel indices in the ladder
I_P12 = kap.index(+1)
I_P32 = kap.index(-2)
E2_p12 = channel_sommerfeld(+1, Z)[2]
E2_p32 = channel_sommerfeld(-2, Z)[2]
print(f"field-free exact: 2p1/2={E2_p12:.6f}  2p3/2={E2_p32:.6f}  "
      f"FS split={E2_p32-E2_p12:.6f}")


def ref_channel(i, lam, guess):
    """Field-free level of channel block i nearest `guess`."""
    H0 = S.build(0.0, 0.0, lam)
    blk = H0[2 * nb * i:2 * nb * (i + 1), 2 * nb * i:2 * nb * (i + 1)]
    ev = np.sort(np.linalg.eigvals(blk).real)
    ev = ev[(ev > -2 * S.c ** 2 + 1e-6) & (ev < -1e-9)]
    return ev[np.argmin(np.abs(ev - guess))]


rows = []
state = {}
for comp, ich, e2 in (("2p1/2", I_P12, E2_p12), ("2p3/2", I_P32, E2_p32)):
    prevE = prevRef = None
    guess = e2
    for u in US:
        lam = u * Z
        E0 = ref_channel(ich, lam, guess)
        guess = E0                     # follow the screened level
        if prevE is None:
            # nonrel n=2 Stark shift is dominated by linear term for the
            # mixed states, but 2p m=1/2 in the coupled basis: seed with
            # a small negative shift and modest width
            pred = E0 - 0.008 * Z ** 2 * F2 / 0.004 - 0.002j * Z ** 2 \
                * F2 / 0.004
            pred = E0 + (-0.13 - 0.033j) * (Z ** 2 * F2 / 0.004) * 1e-2
        else:
            pred = E0 + (prevE - prevRef)
        rad = max(0.35, 5 * abs(pred.imag))
        picks = []
        for th in (0.10, 0.16):
            ev = np.linalg.eigvals(S.build(F, th, lam))
            sel = ev[(np.abs(ev - pred) < rad) & (ev.imag < 1e-7)]
            if sel.size:
                picks.append(sel[np.argmin(np.abs(sel - pred))])
        if len(picks) != 2:
            print(f"   {comp} u={u:.2f}: LOST")
            prevE = prevRef = None
            continue
        e = picks[1]
        q = abs(picks[1] - picks[0])
        rows.append(dict(Z=Z, u=u, lam=lam, F=F, comp=comp, E0ref=E0,
                         Er=e.real, Gamma=-2 * e.imag, qual=q))
        print(f"   {comp} u={u:.2f} lam={lam:5.1f}: dE={e.real-E0:+.5e}  "
              f"Gamma={-2*e.imag:.5e}  (q {q:.0e})")
        prevE, prevRef = e, E0

with (HERE / "f3d_2p_map.csv").open("w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    for r in rows:
        w.writerow(r)

# splitting summary
print("\nGamma ratio Gamma(2p3/2)/Gamma(2p1/2) vs u:")
by = {}
for r in rows:
    by.setdefault(r["u"], {})[r["comp"]] = r["Gamma"]
for u in US:
    d = by.get(u, {})
    if "2p1/2" in d and "2p3/2" in d and d["2p1/2"] > 1e-14:
        print(f"   u={u:.2f}: ratio = {d['2p3/2']/d['2p1/2']:.4f}")

print(f"\ntotal: {time.perf_counter()-t0:.1f} s")
