"""Paper E, Phase E2: the (F, lam) Stark map under ECSC screening.

Same physical-continuation tracker as Paper D v2, with the prediction
anchored to the ECSC F=0 bound curve while it exists (lam < lam_c(1s)
~ 0.7203) and to complex linear extrapolation beyond.  Theta pair
(0.18, 0.28) with retry (0.30, 0.42) -- all < pi/4 (ECSC cap).

E2b: n=2, m=0 doublet at F=0.02 under ECSC (continuation in lam).
"""

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
print(f"[basis lmax=14, nb=70, DOF={15*70}, cert={S.cert:.0e}]")

THETAS = (0.18, 0.28)
THETAS_RETRY = (0.30, 0.42)

# ---------------- ECSC F=0 reference curve ------------------------------
LAMS = np.round(np.arange(0.0, 0.901, 0.05), 3)
print("\nECSC F=0 ground reference:")
eref = {}
for lam in LAMS:
    H = S.build(1.0, lam, 0.0, 0.0, scr="ecsc")
    e = float(np.linalg.eigvalsh(0.5 * (H.real + H.real.T)).min())
    eref[lam] = e if e < -1e-6 else None      # None above lam_c
print("   " + "  ".join(f"{l:.2f}:{eref[l]:+.4f}" if eref[l] is not None
                        else f"{l:.2f}:unbound" for l in LAMS[::3]))

# ---------------- ground-state map --------------------------------------
print("\nE2 ground-state map:")
rows = []
for F in (0.02, 0.04, 0.06, 0.08, 0.10):
    prev, prev_lam, slope, losses = None, None, None, 0
    for lam in LAMS:
        if losses >= 2:
            rows.append(dict(F=F, lam=lam, Er=np.nan, Gamma=np.nan,
                             Eref=eref[lam] if eref[lam] else np.nan,
                             qual=np.nan))
            continue
        if prev is None:
            pred, radius, tol = -0.52 + 0.0j, 0.08, 5e-3
        elif eref[lam] is not None and eref[prev_lam] is not None:
            pred = eref[lam] + (prev - eref[prev_lam])
            radius = max(0.04, 4 * abs(prev.imag),
                         3 * abs(eref[lam] - eref[prev_lam]))
            tol = max(5e-3, 0.2 * abs(prev.imag))
        else:
            pred = prev + (slope or 0.9) * (lam - prev_lam)
            radius = max(0.04, 4 * abs(prev.imag))
            tol = max(8e-3, 0.25 * abs(prev.imag))
        got = track_point(S, F, lam, pred, radius, tol, THETAS)
        if got is None:
            got = track_point(S, F, lam, pred, 1.5 * radius, 2 * tol,
                              THETAS_RETRY)
        if got is None:
            losses += 1
            print(f"   F={F:.2f} lam={lam:.2f}: "
                  + ("dissolved" if losses >= 2 else "lost"))
            rows.append(dict(F=F, lam=lam, Er=np.nan, Gamma=np.nan,
                             Eref=eref[lam] if eref[lam] else np.nan,
                             qual=np.nan))
            continue
        losses = 0
        e, q = got
        rows.append(dict(F=F, lam=lam, Er=e.real, Gamma=-2 * e.imag,
                         Eref=eref[lam] if eref[lam] else np.nan, qual=q))
        print(f"   F={F:.2f} lam={lam:.2f}: E_r={e.real:+.7f}"
              f"  Gamma={-2*e.imag:.4e}  (q {q:.0e})")
        if prev is not None:
            slope = (e - prev) / (lam - prev_lam)
        prev, prev_lam = e, lam
    print()

with (HERE / "e2_map.csv").open("w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    for r in rows:
        w.writerow(r)

# ---------------- E2b: n=2 doublet at F=0.02 ----------------------------
print("E2b n=2 doublet (m=0, k=-1 red / k=+1 blue) at F=0.02, ECSC:")
man = []
comps = {"red": -0.2067 + 0.0j, "blue": -0.0890 + 0.0j}
state = {k: None for k in comps}
LAMS2 = np.round(np.arange(0.0, 0.301, 0.025), 3)
for lam in LAMS2:
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

print(f"\ntotal: {time.perf_counter()-t0:.1f} s")
