"""Paper F, F1: field-free certification of the coupled engine.

With F=0 the coupled Hamiltonian is block-diagonal; each channel's
shifted pencil must reproduce the exact Dirac-Coulomb (Sommerfeld)
binding energies on the SHARED Legendre grid.  Gate at several Z to
confirm the shared grid is adequate before switching the field on.
"""

import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from ds_core import CoupledDiracStark, channel_sommerfeld  # noqa: E402
from angular import kappa_ladder  # noqa: E402

t0 = time.perf_counter()

for Z in (1, 40, 92):
    kap = kappa_ladder(-1, 4, 1)          # s1/2 p1/2 p3/2 d3/2
    S = CoupledDiracStark(nb=40, kappas=kap, Z=Z, beta=max(0.4, 0.6 * Z),
                          m2=1, nq=9000)
    H = S.build(F=0.0, theta=0.0, lam=0.0)   # real, block-diagonal
    print(f"\nZ={Z}  (cert={S.cert:.0e}, gamma={[f'{g:.4f}' for g in S.gamma]}):")
    # per channel: compare lowest few bound states to Sommerfeld.
    # SHIFTED pencil eigenvalues ARE binding energies (no -c^2 subtraction).
    for i, k in enumerate(kap):
        som = channel_sommerfeld(k, Z)
        blk = H[2 * 40 * i:2 * 40 * (i + 1), 2 * 40 * i:2 * 40 * (i + 1)]
        be = np.sort(np.linalg.eigvals(blk).real)
        be = be[(be > -2 * S.c ** 2 + 1e-6) & (be < -1e-9)]
        # match lowest 2 levels
        msg = []
        for n in sorted(som)[:2]:
            e_ex = som[n]
            e_num = be[np.argmin(np.abs(be - e_ex))]
            msg.append(f"n={n}:{abs(e_num-e_ex)/abs(e_ex):.1e}")
        print(f"   kappa={k:+d}: " + "  ".join(msg))

print(f"\ntotal: {time.perf_counter()-t0:.1f} s")
