"""Paper E, Phase E4 (post-audit): (a) re-run the n=2 doublet with the
paired-angle tracker (the lam=0.125 blue point was a single-angle
artifact); (b) the referee's F->0 probe at lam=0.90 and 0.95: does the
field-on ground resonance connect, as F->0, with the field-free
structure (Jiao tail resonance above lam_r+=0.90035 / merged pole
below)?  Report the tan(2 theta2) ratio test per point."""

import csv
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from ecsc_core import StarkChannelsE, track_point  # noqa: E402

t0 = time.perf_counter()
S = StarkChannelsE(70, 14, 0.35)
print(f"[cert {S.cert:.0e}]")

THETAS = (0.18, 0.28)
THETAS_RETRY = (0.30, 0.42)

# ---------------- (a) n=2 doublet re-run --------------------------------
print("\nE4a n=2 doublet, paired-angle tracker:")
man = []
comps = {"red": -0.2067 + 0.0j, "blue": -0.0890 + 0.0j}
state = {k: None for k in comps}
for lam in np.round(np.arange(0.0, 0.301, 0.025), 3):
    for name, seed in comps.items():
        if state[name] == "dead":
            continue
        prev = state[name]
        if prev is None:
            pred, radius, tol = seed, 0.05, 8e-3
        else:
            pred = prev[0] + (lam - prev[1]) * prev[2]
            radius = max(0.035, 3 * abs(prev[0].imag))
            tol = max(8e-3, 0.2 * abs(prev[0].imag))
        got = track_point(S, 0.02, lam, pred, radius, tol, THETAS)
        if got is None:
            got = track_point(S, 0.02, lam, pred, 1.5 * radius, 2 * tol,
                              THETAS_RETRY)
        if got is None:
            print(f"   lam={lam:.3f} {name}: lost -> stop")
            state[name] = "dead"
            continue
        e, q = got
        slope = 0.9 if prev is None else (e - prev[0]) / (lam - prev[1])
        state[name] = (e, lam, slope)
        man.append((name, lam, e.real, -2 * e.imag, q))
        print(f"   lam={lam:.3f} {name}: E_r={e.real:+.6f}"
              f"  Gamma={-2*e.imag:.3e}  (q {q:.0e})")

with (HERE / "e2_manifold.csv").open("w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["comp", "lam", "Er", "Gamma", "qual"])
    for r in man:
        w.writerow(r)

# ---------------- (b) F -> 0 probe --------------------------------------
print("\nE4b F->0 probe (ratio test vs tan(2*0.28)=0.6269):")
for lam in (0.90, 0.95):
    # seed from the map at F=0.02 (or fresh search for lam=0.95)
    prev = None
    for F in (0.02, 0.015, 0.010, 0.007, 0.005, 0.003, 0.002, 0.001):
        if prev is None:
            pred, radius, tol = 0.049 - 0.048j, 0.06, 8e-3
        else:
            pred = prev
            radius = max(0.02, 4 * abs(prev.imag) * 0.5)
            tol = max(5e-3, 0.15 * abs(prev.imag))
        got = track_point(S, F, lam, pred, radius, tol, THETAS)
        if got is None:
            got = track_point(S, F, lam, pred, 1.5 * radius, 2 * tol,
                              THETAS_RETRY)
        if got is None:
            print(f"   lam={lam} F={F:.3f}: lost")
            break
        e, q = got
        ratio = (-e.imag) / e.real if e.real > 0 else float("nan")
        print(f"   lam={lam} F={F:.3f}: E={e.real:+.6f}{e.imag:+.6f}j"
              f"  (G/2)/Er={ratio:.4f}  q={q:.0e}")
        prev = e

print(f"\ntotal: {time.perf_counter()-t0:.1f} s")
