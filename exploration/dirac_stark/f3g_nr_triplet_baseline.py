"""F3g: nonrelativistic n=2 baseline for the triplet ratio claim.

The manuscript (and referee round 1) needs the NONRELATIVISTIC
uphill/downhill width ratio of the n=2, m=0 parabolic doublet at the
scaled field of the relativistic triplet run (f = F/Z^3 = 0.004,
u = lambda/Z = 0..0.20 step 0.025, mirrored at Z=1).  If the nonrel
ratio is monotone in u, the non-monotonicity measured for the Dirac
triplet at Z=40 is genuinely relativistic; if not, the attribution
must be dropped.

Engine and tracking are d3_map_v2.py's (Paper D), unchanged except
F=0.004 and the u grid.  Seeds: E = -1/8 -+ 3F (linear Stark, k=-+1).
"""

import csv
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "exploration" / "stark2d"))

from d0_architecture import StarkChannels  # noqa: E402

t0 = time.perf_counter()
F = 0.004

S = StarkChannels(70, 14, 0.35)
print(f"[basis lmax=14, nb=70, DOF={15*70}, cert={S.cert:.0e}]")

THETAS = (0.18, 0.28)
THETAS_RETRY = (0.30, 0.42)


def track_point(F, lam, pred, radius, tol, thetas):
    es = []
    for th in thetas:
        ev = np.linalg.eigvals(S.build(1.0, lam, F, th))
        ev = ev[(np.abs(ev - pred) < radius) & (ev.imag < 1e-8)]
        if not ev.size:
            return None
        es.append(ev[np.argmin(np.abs(ev - pred))])
    d = abs(es[1] - es[0])
    return (es[1], d) if d <= tol else None


print(f"F3g n=2 doublet (m=0, parabolic k=-1 down / k=+1 up) at F={F}:")
man = []
comps = {"down": -0.125 - 3 * F + 0.0j, "up": -0.125 + 3 * F + 0.0j}
state = {k: None for k in comps}
LAMS2 = np.round(np.arange(0.0, 0.201, 0.025), 3)
for lam in LAMS2:
    for name, seed in comps.items():
        if state[name] == "dead":
            continue
        prev = state[name]
        if prev is None:
            pred, radius, tol = seed, 0.008, 5e-3
        else:
            pred = prev[0] + (lam - prev[1]) * prev[2]
            radius = max(0.02, 3 * abs(prev[0].imag))
            tol = max(5e-3, 0.2 * abs(prev[0].imag))
        got = track_point(F, lam, pred, radius, tol, THETAS)
        if got is None:
            got = track_point(F, lam, pred, 1.5 * radius, 2 * tol,
                              THETAS_RETRY)
        if got is None:
            print(f"   u={lam:.3f} {name}: lost -> stop")
            state[name] = "dead"
            continue
        e, q = got
        slope = 0.9 if prev is None else (e - prev[0]) / (lam - prev[1]) \
            if lam > prev[1] else 0.9
        state[name] = (e, lam, slope)
        man.append((name, lam, e.real, -2 * e.imag, q))
        print(f"   u={lam:.3f} {name}: E_r={e.real:+.7f}"
              f"  Gamma={-2*e.imag:.4e}  (q {q:.0e})")

with (Path(__file__).with_name("f3g_nr_baseline.csv")).open("w",
                                                            newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["comp", "lam", "Er", "Gamma", "qual"])
    for r in man:
        w.writerow(r)

# ratio table
by = {}
for name, lam, er, g, q in man:
    by.setdefault(lam, {})[name] = g
print("\nGamma_up/Gamma_down(u):")
for lam in LAMS2:
    d = by.get(lam, {})
    if "up" in d and "down" in d and d["down"] > 0:
        print(f"   u={lam:.3f}  ratio={d['up']/d['down']:.4f}")

print(f"\ntotal: {time.perf_counter()-t0:.1f} s")
