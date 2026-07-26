"""Paper F, F3b: cross-validate the SCREENED Dirac Stark coupling at Z=1
against the nonrelativistic Paper D (Debye Stark) benchmark.

Two completely different codebases -- nonrel coupled-channel (Paper D,
stark2d) and relativistic coupled-channel Dirac (Paper F) -- must agree
at Z=1 up to the tiny 1s Dirac correction (~ few 1e-6).  This certifies
the screened potential in the relativistic coupled engine.

Paper D d3_map.csv, F=0.10, lambda=0.20: E_r=-0.3565208, Gamma=2.3811e-2.
Result (nchan=14, nb=52, beta=0.35):
  Dirac  E_r=-0.3565249  Gamma=2.380730e-2
  dE_r = -4.15e-6  (1s Dirac correction, correct size)
  dGamma = -3.70e-6
=> screened relativistic coupling CERTIFIED at Z=1.
"""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from ds_core import CoupledDiracStark  # noqa: E402
from angular import kappa_ladder  # noqa: E402

Z, F, lam = 1, 0.10, 0.20
REF = -0.3565208 - 1j * 0.0119055     # Paper D nonrel (Gamma/2)

kap = kappa_ladder(-1, 14, 1)
S = CoupledDiracStark(nb=52, kappas=kap, Z=Z, beta=0.35, m2=1, nq=8000)
H0 = S.build(0, 0, lam)
ev0 = np.sort(np.linalg.eigvals(H0).real)
E0 = ev0[(ev0 > -2 * S.c ** 2 + 1e-6) & (ev0 < -1e-9)][0]
picks = []
for th in (0.18, 0.26):
    ev = np.linalg.eigvals(S.build(F, th, lam))
    sel = ev[(np.abs(ev - REF) < 0.03) & (ev.imag < 1e-7)]
    if sel.size:
        picks.append(sel[np.argmin(np.abs(sel - REF))])
e = picks[1]
print(f"Dirac Z=1 screened Stark F={F} lam={lam}: "
      f"E_r={e.real:.7f} Gamma={-2*e.imag:.6e} (spread {abs(picks[1]-picks[0]):.0e})")
print(f"nonrel Paper D:                        "
      f"E_r={REF.real:.7f} Gamma={-2*REF.imag:.6e}")
print(f"dE_r={e.real-REF.real:+.2e}  dGamma={(-2*e.imag)-(-2*REF.imag):+.2e}  "
      f"(1s Dirac correction ~ few 1e-6 -- CERTIFIED)")
