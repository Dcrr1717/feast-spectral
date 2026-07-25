"""Paper E, E1b: hunt the SHARP above-threshold s resonance for
lam > lam_r+ = 0.90035 (Jiao et al. EPJD 75, 313 (2021)) -- a distinct
pole of the oscillatory-tail well, separate from the broad 1s-descendant
tracked in E1.  Strategy: dense theta-stable scan of the l=0 spectrum
near threshold, several bases, theta-trajectory verification."""

import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from ecsc_core import StarkChannelsE  # noqa: E402

t0 = time.perf_counter()

bases = [StarkChannelsE(200, 0, 0.15), StarkChannelsE(260, 0, 0.10)]
print(f"[certs: {bases[0].cert:.0e}, {bases[1].cert:.0e}]")


def stable_scan(S, lam, th1, th2, emax=0.08, gmax=0.06):
    """All theta-stable eigenvalues near threshold."""
    evs = []
    for th in (th1, th2):
        ev = np.linalg.eigvals(S.build(1.0, lam, 0.0, th, scr="ecsc"))
        ev = ev[(ev.real > -0.01) & (ev.real < emax)
                & (ev.imag < 1e-9) & (ev.imag > -gmax)]
        evs.append(ev)
    out = []
    if not evs[0].size or not evs[1].size:
        return out
    d = np.abs(evs[1][None, :] - evs[0][:, None])
    while True:
        i, j = np.unravel_index(np.argmin(d), d.shape)
        if d[i, j] > 2e-4:
            break
        out.append((evs[1][j], d[i, j]))
        d[i, :] = np.inf
        d[:, j] = np.inf
    return out


for lam in (0.85, 0.90, 0.92, 0.95, 1.00, 1.10, 1.25):
    hits = stable_scan(bases[0], lam, 0.35, 0.55)
    tag = "  ".join(f"({e.real:+.6f},{-2*e.imag:.2e};q{q:.0e})"
                    for e, q in hits) or "none"
    print(f"lam={lam:5.3f}: {tag}")

# theta-trajectory + basis check on the best candidates
print("\ntrajectories for the sharp candidates:")
for lam in (0.95, 1.00, 1.10):
    hits = stable_scan(bases[0], lam, 0.35, 0.55)
    if not hits:
        continue
    # sharpest candidate
    e0 = min(hits, key=lambda t: -t[0].imag / max(t[0].real, 1e-6))[0]
    print(f"  lam={lam}: candidate {e0.real:+.7f}{e0.imag:+.7f}j")
    for S, tagb in ((bases[0], "nb200"), (bases[1], "nb260")):
        es = []
        for th in (0.30, 0.40, 0.50, 0.60, 0.70):
            ev = np.linalg.eigvals(S.build(1.0, lam, 0.0, th, scr="ecsc"))
            es.append(ev[np.argmin(np.abs(ev - e0))])
        drift = max(abs(a - b) for a in es for b in es)
        print(f"    {tagb}: " + "  ".join(f"{e.real:+.6f}{e.imag:+.6f}j"
                                          for e in es[:3])
              + f"   full-traj drift {drift:.1e}")

print(f"\ntotal: {time.perf_counter()-t0:.1f} s")
