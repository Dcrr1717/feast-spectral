"""Paper E, Phase E0: certification gates for the ECSC engine.

P-E0.0  Hydrogen limit (lam=0): lowest levels vs -1/2n^2, dev < 1e-9
        (RAISE on failure, per D4 discipline).
P-E0.1  Yukawa regression: scr='yukawa' build must reproduce
        d0_architecture.StarkChannels.build to machine precision.
P-E0.2  ECSC small-lam physics: E_1s(lam) = -1/2 + lam - lam^3<r^2>
        + O(lam^4) with <r^2>_1s = 3 (NO lam^2 term -- differs from
        Yukawa; exact expansion e^{-x}cos x = 1 - x + x^3/3 - x^4/6).
P-E0.3  Field-free 1s energy curve + critical screening lam_c(1s);
        internal high-precision (two bases), literature cross-check
        recorded separately once digits are confirmed.
P-E0.4  THE PHYSICS GATE: hunt the field-free l=0 resonance above
        lam_c (the oscillatory-tail branch, absent for Yukawa l=0);
        theta-stability with theta < pi/4.
"""

import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "stark2d"))

from ecsc_core import StarkChannelsE, track_point  # noqa: E402
from d0_architecture import StarkChannels  # noqa: E402  (runs P-D0.0)

t0 = time.perf_counter()

# ---------------- P-E0.0 hydrogen limit ---------------------------------
S = StarkChannelsE(60, 2, 0.35)
ev0 = np.sort(np.linalg.eigvals(S.build(1.0, 0.0, 0.0, 0.0)).real)
exact = [-0.5, -0.125, -0.125, -1.0 / 18, -1.0 / 18]
dev0 = max(abs(a - b) for a, b in zip(ev0[:5], exact))
print(f"P-E0.0 hydrogen limit: dev = {dev0:.1e}"
      + ("  OK" if dev0 < 1e-9 else "  FAIL"))
if not (dev0 < 1e-9):
    raise RuntimeError("P-E0.0 FAILED")

# ---------------- P-E0.1 yukawa regression ------------------------------
S_ref = StarkChannels(40, 2, 0.35)
S_new = StarkChannelsE(40, 2, 0.35)
dmax = 0.0
for lam, F, th in ((0.25, 0.04, 0.22), (0.6, 0.1, 0.35), (0.0, 0.08, 0.18)):
    Ha = S_ref.build(1.0, lam, F, th)
    Hb = S_new.build(1.0, lam, F, th, scr="yukawa")
    dmax = max(dmax, float(np.max(np.abs(Ha - Hb))))
print(f"P-E0.1 yukawa regression: max|dH| = {dmax:.1e}"
      + ("  OK" if dmax < 1e-14 else "  FAIL"))
if not (dmax < 1e-14):
    raise RuntimeError("P-E0.1 FAILED")

# ---------------- P-E0.2 small-lam expansion ----------------------------
print("\nP-E0.2 E_1s(lam) vs -1/2 + lam - 3*lam^3/3*<r^2>/3 ... exact "
      "-1/2 + lam - lam^3*<r^2>_pert:")
Sb = StarkChannelsE(80, 0, 0.35)
for lam in (0.01, 0.02, 0.05):
    H = Sb.build(1.0, lam, 0.0, 0.0)
    e = float(np.linalg.eigvalsh(0.5 * (H.real + H.real.T)).min())
    pred = -0.5 + lam - lam**3   # -(lam^3/3)<r^2>_1s = -lam^3 (<r^2>=3)
    print(f"   lam={lam}: E={e:+.9f}  PT={pred:+.9f}  diff={e-pred:+.2e}")

# ---------------- P-E0.3 1s curve + critical screening ------------------
print("\nP-E0.3 field-free ECSC 1s energy (nb=80 vs nb=140, beta tuned):")
Sc = StarkChannelsE(140, 0, 0.25)
lam_grid = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.65, 0.7, 0.72)
E1s = {}
for lam in lam_grid:
    Ha = Sb.build(1.0, lam, 0.0, 0.0)
    Hb2 = Sc.build(1.0, lam, 0.0, 0.0)
    ea = float(np.linalg.eigvalsh(0.5 * (Ha.real + Ha.real.T)).min())
    eb = float(np.linalg.eigvalsh(0.5 * (Hb2.real + Hb2.real.T)).min())
    E1s[lam] = eb
    print(f"   lam={lam}: E={eb:+.9f}   grid-drift={abs(ea-eb):.1e}")

# bisect lam_c(1s): E crosses 0
lo, hi = 0.6, 0.9
Sf = StarkChannelsE(160, 0, 0.20)
for _ in range(40):
    mid = 0.5 * (lo + hi)
    H = Sf.build(1.0, mid, 0.0, 0.0)
    e = float(np.linalg.eigvalsh(0.5 * (H.real + H.real.T)).min())
    if e < 0:
        lo = mid
    else:
        hi = mid
print(f"   lam_c(1s) ~ {0.5*(lo+hi):.6f}   (variational upper-ish; "
      "threshold-diffuse caveat as in the NR paper)")

# ---------------- P-E0.4 the oscillatory-tail resonance -----------------
print("\nP-E0.4 l=0 resonance hunt above lam_c (theta-stable, th<pi/4):")
Sr = StarkChannelsE(120, 0, 0.25)
for lam in (0.75, 0.80, 0.90, 1.00):
    got = None
    for pred in (0.002 + 0.0j, 0.005 + 0.0j, 0.01 + 0.0j, 0.02 + 0.0j,
                 0.05 + 0.0j):
        got = track_point(Sr, 0.0, lam, pred, 0.05, 5e-3,
                          (0.35, 0.55), scr="ecsc")
        if got:
            break
    if got:
        e, q = got
        print(f"   lam={lam}: E_r={e.real:+.7f}  Gamma={-2*e.imag:.3e}"
              f"  (q {q:.0e})")
    else:
        print(f"   lam={lam}: none found in windows")

print(f"\ntotal: {time.perf_counter()-t0:.1f} s")
