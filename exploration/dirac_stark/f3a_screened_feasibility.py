"""Paper F, F3a: FIRST relativistic screened Stark resonances.

Take the Maltsev-validated Z=10 point and switch on Debye/Yukawa
screening lambda (a.u.).  Screening lifts the level (less binding) and
should fatten the width -- but does the RELATIVISTIC width suppression
behave differently than the nonrelativistic Paper D catalysis?  This is
the first look at Dirac + screening + field, nobody has computed it.

Field-free screened reference E0(lambda) taken from the engine (F=0).
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

Z = 10
F = 0.05 * Z ** 3          # = 50, the Maltsev field
nchan, nb, beta = 14, 52, 5.0
kap = kappa_ladder(-1, nchan, 1)
S = CoupledDiracStark(nb=nb, kappas=kap, Z=Z, beta=beta, m2=1, nq=8500)
print(f"Z={Z} F={F} (nchan={nchan} nb={nb} beta={beta} cert={S.cert:.0e})")
print(f"First relativistic screened Stark resonances (1s), scan lambda:\n")
print(f"{'lambda':>7} {'E0(F=0)':>14} {'E_r':>14} {'Gamma':>13} "
      f"{'dE_shift':>12} {'spread':>8}")

prev = None
for lam in (0.0, 1.0, 2.0, 4.0, 7.0, 10.0, 14.0):
    # field-free screened reference (1s binding)
    H0 = S.build(0.0, 0.0, lam)
    ev0 = np.sort(np.linalg.eigvals(H0).real)
    ev0 = ev0[(ev0 > -2 * S.c ** 2 + 1e-6) & (ev0 < -1e-9)]
    E0 = ev0[0]
    # predict from previous or from E0 + last shift
    if prev is None:
        pred = E0 - 6.07e-1 - 3.75e-3j          # Maltsev lam=0 seed
    else:
        pred = E0 + (prev[1] - prev[0])          # E0 + previous complex shift
    picks = []
    for th in (0.12, 0.18):
        ev = np.linalg.eigvals(S.build(F, th, lam))
        sel = ev[(np.abs(ev - pred) < max(0.05, 3 * abs(pred.imag)))
                 & (ev.imag < 1e-7)]
        if sel.size:
            picks.append(sel[np.argmin(np.abs(sel - pred))])
    if len(picks) != 2:
        print(f"{lam:7.1f}  lost")
        continue
    e = picks[1]
    spread = abs(picks[1] - picks[0])
    dE = e.real - E0
    print(f"{lam:7.1f} {E0:14.6f} {e.real:14.6f} {-2*e.imag:13.5e} "
          f"{dE:+12.4e} {spread:8.1e}")
    prev = (E0, e)

print(f"\ntotal: {time.perf_counter()-t0:.1f} s")
