"""Paper F, F2 feasibility: coupled Dirac Stark at Z=1, F=0.1.

The 1s resonance should match the nonrelativistic hydrogen Stark value
(Paper D / classical: E_r=-0.527418, Gamma/2=0.0072690) plus the tiny
1s Dirac shift (~ -6.6e-6).  Found by theta-stability of the coupled
complex-scaled pencil.  This validates the field coupling + complex
scaling end-to-end before the Maltsev Z=1-82 ladder.
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

NONREL = -0.5274181749 - 0.0072690567j    # Paper D F=0.10 ground

for nchan in (6, 10, 14):
    kap = kappa_ladder(-1, nchan, 1)
    S = CoupledDiracStark(nb=48, kappas=kap, Z=1, beta=0.35, m2=1,
                          nq=7000)
    picks = []
    for th in (0.20, 0.28):
        H = S.build(F=0.10, theta=th, lam=0.0)
        ev = np.linalg.eigvals(H)
        sel = ev[(np.abs(ev - NONREL) < 0.05) & (ev.imag < 1e-8)]
        if sel.size:
            picks.append(sel[np.argmin(np.abs(sel - NONREL))])
    if len(picks) == 2:
        e = picks[1]
        spread = abs(picks[1] - picks[0])
        print(f"nchan={nchan} (kappas={kap}): E={e.real:.10f}{e.imag:+.10f}j"
              f"  th-spread={spread:.1e}")
        print(f"   vs nonrel {NONREL.real:.7f}{NONREL.imag:+.7f}j: "
              f"dRe={e.real-NONREL.real:+.2e}  dGamma/2={-e.imag-(-NONREL.imag):+.2e}")
    else:
        print(f"nchan={nchan}: not found ({len(picks)} picks)")

print(f"\ntotal: {time.perf_counter()-t0:.1f} s")
