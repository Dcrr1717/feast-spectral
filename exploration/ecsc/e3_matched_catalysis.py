"""Paper E, Phase E3: Debye vs ECSC field ionization at MATCHED binding.

The clean question: at equal field-free ground binding E_b (i.e. after
removing the trivial 'screening lifts the level' effect), does the
oscillatory ECSC tail change the field-ionization width?  The
quasistatic exponent depends only on |E_b| and F, so any residual
Gamma difference isolates barrier-shape / tail physics.

Protocol: for anchors E_b in {-0.40, -0.30, -0.20, -0.12}, bisect
lam_Y (Yukawa) and lam_E (ECSC) so the F=0 ground level equals E_b to
1e-9; then compute the ground Stark resonance at F in {0.04, 0.08}
for both potentials with the production basis; report Gamma ratios.
"""

import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from ecsc_core import StarkChannelsE, track_point  # noqa: E402

t0 = time.perf_counter()

Sg = StarkChannelsE(120, 0, 0.25)     # ground-level bisection basis
S = StarkChannelsE(70, 14, 0.35)      # production coupled basis
print(f"[cert: ground {Sg.cert:.0e}, coupled {S.cert:.0e}]")


def ground(lam, scr):
    H = Sg.build(1.0, lam, 0.0, 0.0, scr=scr)
    return float(np.linalg.eigvalsh(0.5 * (H.real + H.real.T)).min())


def match_lam(Eb, scr, lo, hi):
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        if ground(mid, scr) > Eb:      # too shallow -> reduce lam
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


ANCHORS = (-0.40, -0.30, -0.20, -0.12)
FIELDS = (0.04, 0.08)
print(f"\n{'E_b':>6} {'lam_Y':>9} {'lam_E':>9} {'F':>5} "
      f"{'Gamma_Y':>11} {'Gamma_E':>11} {'ratio E/Y':>9}")
for Eb in ANCHORS:
    lamY = match_lam(Eb, "yukawa", 0.0, 1.15)
    lamE = match_lam(Eb, "ecsc", 0.0, 0.72)
    for F in FIELDS:
        out = {}
        for scr, lam in (("yukawa", lamY), ("ecsc", lamE)):
            pred, radius, tol = complex(Eb), 0.08, 5e-3
            # crude two-step continuation in F for robustness
            got = track_point(S, F / 2, lam, pred, radius, tol,
                              (0.18, 0.28), scr=scr)
            if got:
                pred = got[0]
            got = track_point(S, F, lam, pred,
                              max(0.05, 6 * abs(pred.imag)),
                              max(5e-3, 0.3 * abs(pred.imag)),
                              (0.18, 0.28), scr=scr)
            out[scr] = got
        if out["yukawa"] and out["ecsc"]:
            gY = -2 * out["yukawa"][0].imag
            gE = -2 * out["ecsc"][0].imag
            ratio = gE / gY if gY > 1e-13 else float("nan")
            print(f"{Eb:+6.2f} {lamY:9.5f} {lamE:9.5f} {F:5.2f} "
                  f"{gY:11.4e} {gE:11.4e} {ratio:9.3f}")
        else:
            print(f"{Eb:+6.2f} {lamY:9.5f} {lamE:9.5f} {F:5.2f}  LOST "
                  f"Y={bool(out['yukawa'])} E={bool(out['ecsc'])}")

print(f"\ntotal: {time.perf_counter()-t0:.1f} s")
