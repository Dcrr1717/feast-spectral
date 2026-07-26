"""Paper F, F2b: validate against Maltsev PRA 108, 042813 (2023).

Their scaled field is F/Z^3 = 0.05 (a.u.); ADE is the Stark SHIFT
relative to the field-free Dirac level, -Gamma/2 the width.  We compute
E(F)-E(F=0) from our engine and compare.

Maltsev 1s table (F/Z^3 = 0.05):
  Z    dE            -Gamma/2
  1   -6.105058e-3  -3.858129e-5
  10  -6.068738e-1  -3.752554e-3
  18  -1.939745      -1.141133e-2
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

MALTSEV = {   # Z: (dE, -Gamma/2)  at F/Z^3 = 0.05
    1:  (-6.105058e-3, -3.858129e-5),
    10: (-6.068738e-1, -3.752554e-3),
    18: (-1.939745,    -1.141133e-2),
}


def find_1s(S, F, lam, pred, thetas, rad, tol):
    picks = []
    for th in thetas:
        ev = np.linalg.eigvals(S.build(F, th, lam))
        sel = ev[(np.abs(ev - pred) < rad) & (ev.imag < 1e-7)]
        if not sel.size:
            return None
        picks.append(sel[np.argmin(np.abs(sel - pred))])
    d = abs(picks[1] - picks[0])
    return (picks[1], d) if d < tol else (picks[1], d)  # return anyway w/ spread


for Z, (dE_ref, hG_ref) in MALTSEV.items():
    F = 0.05 * Z ** 3
    nchan = 12
    kap = kappa_ladder(-1, nchan, 1)
    beta = max(0.35, 0.55 * Z)
    nb = 48
    S = CoupledDiracStark(nb=nb, kappas=kap, Z=Z, beta=beta, m2=1, nq=8000)
    # field-free reference (1s binding)
    H0 = S.build(0.0, 0.0, 0.0)
    ev0 = np.sort(np.linalg.eigvals(H0).real)
    ev0 = ev0[(ev0 > -2 * S.c ** 2 + 1e-6) & (ev0 < -1e-9)]
    E0 = ev0[0]                                  # 1s binding
    pred = E0 + dE_ref + 1j * hG_ref
    thetas = (0.10, 0.16) if Z >= 10 else (0.18, 0.26)
    got = find_1s(S, F, 0.0, pred, thetas, rad=abs(dE_ref) * 3 + 0.02,
                  tol=1.0)
    if got is None:
        print(f"Z={Z} F={F:.4g}: 1s resonance NOT FOUND")
        continue
    e, spread = got
    dE = e.real - E0
    hG = e.imag
    print(f"Z={Z:2d} F={F:.4g} (beta={beta:.2f}): "
          f"E0={E0:.8f}")
    print(f"   dE   = {dE:+.7e}  (ref {dE_ref:+.7e}, "
          f"rel {abs(dE-dE_ref)/abs(dE_ref):.1e})")
    print(f"   -G/2 = {hG:+.7e}  (ref {hG_ref:+.7e}, "
          f"rel {abs(hG-hG_ref)/abs(hG_ref):.1e})   th-spread {spread:.1e}")

print(f"\ntotal: {time.perf_counter()-t0:.1f} s")
