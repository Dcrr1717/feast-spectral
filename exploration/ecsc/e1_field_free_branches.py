"""Paper E, Phase E1: field-free ECSC pole structure across lam_c.

Question (the paper's physics headline candidate): what happens to the
1s (l=0) and 2p (l=1) poles when ECSC screening pushes them through
E=0?  The oscillatory tail supports a barrier, so a GENUINE resonance
(not just a virtual state) is possible even for l=0 -- unlike Yukawa.

Method: continuation in lam with fine steps (the Paper D tracker),
prediction = previous pole + slope*dlam; artifact discrimination:
(i) the pole must MOVE with lam (static E = rotated-continuum artifact,
Paper B lesson), (ii) theta-trajectory stationarity over >= 4 angles
< pi/4, (iii) basis-variation stability (two bases).
"""

import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from ecsc_core import StarkChannelsE  # noqa: E402

t0 = time.perf_counter()


def spectrum(S, lam, th, l):
    H = S.build(1.0, lam, 0.0, th, scr="ecsc")
    return np.linalg.eigvals(H)


def nearest(S, lam, th, target, l=0):
    ev = spectrum(S, lam, th, l)
    return ev[np.argmin(np.abs(ev - target))]


def bound_energy(S, lam, l):
    H = S.build(1.0, lam, 0.0, 0.0, scr="ecsc")
    ev = np.linalg.eigvalsh(0.5 * (H.real + H.real.T))
    return float(ev[0]) if l == 0 else float(ev[0])


# ---------------- l=0: track the 1s pole through lam_c ------------------
print("l=0: 1s pole continuation across lam_c ~ 0.7203")
S0 = StarkChannelsE(200, 0, 0.15)
S0b = StarkChannelsE(150, 0, 0.22)
print(f"   [bases: nb=200 beta=0.15 cert={S0.cert:.0e}; "
      f"nb=150 beta=0.22 cert={S0b.cert:.0e}]")

# below threshold: real bound energies, get the slope for prediction
lam_lo = (0.69, 0.70, 0.71)
eb = [bound_energy(S0, l, 0) for l in lam_lo]
slope = (eb[-1] - eb[-2]) / 0.01
print("   below: " + "  ".join(f"{l}:{e:+.7f}" for l, e in zip(lam_lo, eb))
      + f"   slope~{slope:+.3f}")

prev = complex(eb[-1])
prev_lam = lam_lo[-1]
print("   above (theta pair 0.30/0.50, prediction = prev + slope*dlam):")
traj = []
for lam in (0.715, 0.72, 0.725, 0.73, 0.74, 0.75, 0.77, 0.80, 0.85,
            0.90, 1.00, 1.10, 1.25, 1.40):
    pred = prev + slope * (lam - prev_lam)
    picks = []
    for th in (0.30, 0.50):
        picks.append(nearest(S0, lam, th, pred))
    d = abs(picks[1] - picks[0])
    e = picks[1]
    # basis check
    e2 = nearest(S0b, lam, 0.50, pred)
    db = abs(e2 - e)
    moved = abs(e - prev)
    print(f"   lam={lam:5.3f}: E={e.real:+.7f}{e.imag:+.7f}j  "
          f"q(th)={d:.0e}  q(basis)={db:.0e}  dE/dlam~{moved/(lam-prev_lam):+.3f}")
    traj.append((lam, e))
    slope = (e - prev) / (lam - prev_lam)
    prev, prev_lam = e, lam

# ---------------- theta-trajectory at two checkpoints -------------------
print("\n   theta-trajectories (artifact discrimination):")
for lam in (0.75, 1.00):
    seed = [e for l, e in traj if abs(l - lam) < 1e-9][0]
    es = []
    for th in (0.25, 0.35, 0.45, 0.55, 0.65, 0.75):
        e = nearest(S0, lam, th, seed)
        es.append(e)
        print(f"   lam={lam} th={th:.2f}: {e.real:+.8f}{e.imag:+.8f}j")
    core = es[1:5]
    drift = max(abs(a - b) for a in core for b in core)
    print(f"   -> plateau drift th in [0.35,0.65]: {drift:.1e}\n")

# ---------------- l=1: 2p analog ----------------------------------------
print("l=1: 2p pole continuation (lam_c(2p) first)")
S1 = StarkChannelsE(200, 1, 0.15)
# find lam_c(2p): lowest l=1 eigenvalue crosses 0 (use l=1 block only:
# lmax=1 basis has l=0 and l=1; take the lowest eigenvalue with l=1
# character = second block; simpler: build with lmax=1 and look at the
# l=1 sub-block alone)


def bound_l1(S, lam):
    nb = S.nb
    H = S.build(1.0, lam, 0.0, 0.0, scr="ecsc")
    H11 = H[nb:, nb:]
    return float(np.linalg.eigvalsh(0.5 * (H11.real + H11.real.T))[0])


lo, hi = 0.1, 0.6
for _ in range(35):
    mid = 0.5 * (lo + hi)
    if bound_l1(S1, mid) < 0:
        lo = mid
    else:
        hi = mid
lamc2p = 0.5 * (lo + hi)
print(f"   lam_c(2p) ~ {lamc2p:.6f}")

eb = [bound_l1(S1, lamc2p - d) for d in (0.02, 0.01)]
slope = (eb[1] - eb[0]) / 0.01
prev = complex(eb[1])
prev_lam = lamc2p - 0.01


def nearest_l1(S, lam, th, target):
    nb = S.nb
    H = S.build(1.0, lam, 0.0, th, scr="ecsc")
    ev = np.linalg.eigvals(H[nb:, nb:])
    return ev[np.argmin(np.abs(ev - target))]


print("   above lam_c(2p):")
for dl in (0.005, 0.01, 0.02, 0.04, 0.07, 0.11, 0.16, 0.22, 0.30):
    lam = lamc2p + dl
    pred = prev + slope * (lam - prev_lam)
    picks = [nearest_l1(S1, lam, th, pred) for th in (0.30, 0.50)]
    d = abs(picks[1] - picks[0])
    e = picks[1]
    print(f"   lam={lam:5.3f}: E={e.real:+.7f}{e.imag:+.7f}j  q={d:.0e}")
    slope = (e - prev) / (lam - prev_lam)
    prev, prev_lam = e, lam

print(f"\ntotal: {time.perf_counter()-t0:.1f} s")
