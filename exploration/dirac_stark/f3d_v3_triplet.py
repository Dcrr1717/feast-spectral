"""Paper F, F3d v3: the n=2 (m=1/2) manifold = THREE states at Z=40.

Count: 2s1/2 and 2p1/2 (exactly Dirac-degenerate) + 2p3/2 (m=1/2).
Under the field the degenerate pair splits into TWO Stark components;
v2's per-label ramp collapsed both onto one.  v3 disambiguates: at the
first ramp step the TWO theta-stable eigenvalues nearest the degenerate
energy are captured and tracked separately (labels A=lower, B=upper by
continuation); 2p3/2 (nondegenerate) tracks as C.  Finer lambda steps
(du=0.025) + width-scaled radii extend the broad-state tracking.

Output: f3d_triplet.csv (u, lam, comp, Er, Gamma, qual)
"""

import csv
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from ds_core import CoupledDiracStark, channel_sommerfeld  # noqa: E402
from angular import kappa_ladder  # noqa: E402

t0 = time.perf_counter()

Z = 40
F2 = 0.004
F = F2 * Z ** 3
US = tuple(round(0.025 * k, 3) for k in range(9))       # 0..0.20
kap = kappa_ladder(-1, 12, 1)
nb = 52
S = CoupledDiracStark(nb=nb, kappas=kap, Z=Z, beta=0.55 * Z, m2=1,
                      nq=8500)
E_12 = channel_sommerfeld(-1, Z)[2]          # degenerate 2s1/2=2p1/2
E_32 = channel_sommerfeld(-2, Z)[2]
print(f"Z={Z} F={F:.4g}: E(j=1/2)={E_12:.5f} (x2 degenerate), "
      f"E(2p3/2)={E_32:.5f}, cert {S.cert:.0e}")


def eigs(Fval, lam, th):
    return np.linalg.eigvals(S.build(Fval, th, lam))


def stable_near(Fval, lam, pred, rad, tol, kmax=3, thetas=(0.10, 0.16)):
    """Up to kmax theta-stable eigenvalues within rad of pred."""
    A = eigs(Fval, lam, thetas[0])
    A = A[(np.abs(A - pred) < rad) & (A.imag < 1e-7)]
    B = eigs(Fval, lam, thetas[1])
    B = B[(np.abs(B - pred) < rad) & (B.imag < 1e-7)]
    out = []
    if not A.size or not B.size:
        return out
    A = sorted(A, key=lambda z: abs(z - pred))[:kmax + 2]
    used = set()
    for a in A:
        j = int(np.argmin(np.abs(B - a)))
        if j in used:
            continue
        d = abs(B[j] - a)
        if d <= tol:
            out.append((B[j], d))
            used.add(j)
        if len(out) >= kmax:
            break
    return out


# ---- F-ramp with pair disambiguation -----------------------------------
print("\nF-ramp (pair split at first step):")
RAMP = [F * s for s in (0.05, 0.15, 0.3, 0.5, 0.75, 1.0)]
# first step: capture the TWO nearest to the degenerate level + one for C
comp = {}
got = stable_near(RAMP[0], 0.0, complex(E_12), rad=3.0, tol=0.3, kmax=2)
if len(got) >= 2:
    pair = sorted([g[0] for g in got], key=lambda z: z.real)
    comp["A"] = [complex(E_12), pair[0]]
    comp["B"] = [complex(E_12), pair[1]]
    print(f"   pair at F={RAMP[0]:.3g}: A={pair[0]:.5f}  B={pair[1]:.5f}")
else:
    print("   pair capture FAILED", got)
gc = stable_near(RAMP[0], 0.0, complex(E_32), rad=3.0, tol=0.3, kmax=1)
if gc:
    comp["C"] = [complex(E_32), gc[0][0]]
    print(f"   C (2p3/2) at F={RAMP[0]:.3g}: {gc[0][0]:.5f}")

for name in list(comp):
    hist = comp[name]
    ok = True
    for i, Fv in enumerate(RAMP[1:], start=1):
        e_prev, F_prev = hist[-1], RAMP[i - 1]
        slope = (hist[-1] - hist[-2]) / (F_prev - (RAMP[i - 2] if i >= 2
                                                  else 0.0))
        pred = e_prev + slope * (Fv - F_prev)
        rad = max(1.5, 6 * abs(pred.imag))
        got = stable_near(Fv, 0.0, pred, rad, tol=0.5, kmax=1)
        if not got:
            print(f"   {name}: ramp LOST at F={Fv:.3g}")
            ok = False
            break
        hist.append(got[0][0])
    if ok:
        e = hist[-1]
        print(f"   {name}: E(F2)={e.real:+.5f}{e.imag:+.5f}j  "
              f"dE={e.real-(E_12 if name in 'AB' else E_32):+.4f}")
        comp[name] = e
    else:
        del comp[name]

# ---- lambda continuation -----------------------------------------------
print("\nlambda continuation:")
rows = []
prevE = dict(comp)
prevL = {n: 0.0 for n in comp}
slopeL = {n: 0.0 + 0.0j for n in comp}
for u in US:
    if u == 0.0:
        for n, e in comp.items():
            rows.append(dict(u=u, lam=0.0, comp=n, Er=e.real,
                             Gamma=-2 * e.imag, qual=0.0))
        continue
    lam = u * Z
    for n in list(prevE):
        pred = prevE[n] + slopeL[n] * (lam - prevL[n])
        rad = max(0.8, 6 * abs(pred.imag))
        tol = max(0.3, 0.25 * abs(pred.imag))
        got = stable_near(F, lam, pred, rad, tol, kmax=1)
        if not got:
            print(f"   {n} u={u:.3f}: LOST")
            del prevE[n]
            continue
        e, q = got[0]
        rows.append(dict(u=u, lam=lam, comp=n, Er=e.real,
                         Gamma=-2 * e.imag, qual=q))
        print(f"   {n} u={u:.3f}: Er={e.real:+.4f}  "
              f"Gamma={-2*e.imag:.5e}  (q {q:.0e})")
        slopeL[n] = (e - prevE[n]) / (lam - prevL[n])
        prevE[n], prevL[n] = e, lam

with (HERE / "f3d_triplet.csv").open("w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    for r in rows:
        w.writerow(r)

print("\nwidth ratios vs u (C/A ~ FS-partner contrast):")
by = {}
for r in rows:
    by.setdefault(r["u"], {})[r["comp"]] = r["Gamma"]
for u in US:
    d = by.get(u, {})
    msg = []
    for a, b in (("C", "A"), ("B", "A")):
        if a in d and b in d and d[b] > 1e-14:
            msg.append(f"{a}/{b}={d[a]/d[b]:.3f}")
    if msg:
        print(f"   u={u:.3f}: " + "  ".join(msg))

print(f"\ntotal: {time.perf_counter()-t0:.1f} s")
