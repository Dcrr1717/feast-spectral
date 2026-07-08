"""E7: external validation against published Yukawa benchmarks.

Reference values (all DOIs verified 2026-07-08):
- Bound states near threshold and E(lambda):
  del Valle & Nader, J. Math. Phys. 59, 102103 (2018), 10.1063/1.5050621
  Alhaidari, Bahlouli, Abdelmonem, J. Phys. A 41, 032001 (2008),
  10.1088/1751-8113/41/3/032001  (J-matrix + complex scaling, 12 digits)
  Gebremedhin & Weatherford, Phys. Rev. E 108, 045301 (2023)
- Critical screening:
  Diaz, Fernandez, Castro, J. Phys. A 24, 2061 (1991)
  Jiao et al., J. Phys. B 54, 175002 (2021)
- Resonances (complex rotation, 12 digits):
  Bylicki, Stachow, Karwowski, Mukherjee, Chem. Phys. 331, 346 (2007)
  as tabulated by Alhaidari et al. (2008), Table 2-4.

This script recomputes every benchmarked quantity with the Laguerre+FEAST
pipeline and reports digit-level agreement. Output:
results/E7_benchmark_comparison.txt
"""

import sys
import time
from pathlib import Path

import numpy as np
from scipy.linalg import eigvalsh

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from basis import BasisSet                                  # noqa: E402
from run_resonances import theta_trajectory, find_resonance  # noqa: E402

RES = ROOT / "results"
RES.mkdir(exist_ok=True)
OUT = []


def log(s=""):
    print(s)
    OUT.append(s)


# ------------------------------------------------- bound-state anchors --
def bound_energy(lam, ell, index, N=320, beta=0.05):
    B = BasisSet(N, ell, beta)
    return eigvalsh(B.kinetic() + B.yukawa(lam=lam))[index]


BOUND_REFS = [
    # (label, lam, ell, index, E_published, source)
    ("1s l=0.1",  0.10, 0, 0, -0.407058030613403, "Gebremedhin23/delValle18"),
    ("1s l=0.2",  0.20, 0, 0, -0.3268085113691935, "delValle18"),
    ("2s l=0.1",  0.10, 0, 1, -0.04992827133191,  "delValle18/Alhaidari08"),
    ("3s l=0.05", 0.05, 0, 2, -0.01935255481475,  "delValle18"),
    ("2p l=0.2",  0.20, 1, 0, -4.101646530802e-3, "Alhaidari08 (12 digits)"),
    ("2p l=0.22007434", 0.22007434, 1, 0, -1.867969736495e-5, "Alhaidari08"),
    ("3p l=0.11260485", 0.11260485, 1, 1, -6.6870176e-6,      "Alhaidari08"),
    ("3d l=0.091",      0.091,      2, 0, -7.767498169649e-5, "Alhaidari08"),
]


def e7a_bound():
    log("== E7a: bound-state benchmarks =====================================")
    log(f"{'state':>18} {'E (this work)':>20} {'E (published)':>20} {'diff':>10}")
    for label, lam, ell, idx, Eref, src in BOUND_REFS:
        E = bound_energy(lam, ell, idx)
        log(f"{label:>18} {E:>20.12e} {Eref:>20.12e} {abs(E-Eref):>10.1e}  [{src}]")
    log()


# ------------------------------------------------- refined lambda_c -----
def critical_lambda_refined(ell, want, lo, hi, N=320, beta=0.05,
                            eneg=-1e-10, bits=34):
    B = BasisSet(N, ell, beta)
    T = B.kinetic()

    def nbound(lam):
        return int(np.sum(eigvalsh(T + B.yukawa(lam=lam)) < eneg))

    assert nbound(lo) >= want and nbound(hi) < want
    for _ in range(bits):
        mid = 0.5 * (lo + hi)
        if nbound(mid) >= want:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


LC_REFS = [
    # (state, ell, want, lo, hi, lambda_c published, source)
    ("2p", 1, 1, 0.2199, 0.2206, 0.220216806605, "delValle18/Diaz91"),
    ("3p", 1, 2, 0.1125, 0.1130, 0.112710498359, "delValle18/Diaz91"),
    ("3d", 2, 1, 0.0912, 0.0916, 0.091345120771732, "delValle18=Diaz91"),
]


def e7b_lambdac():
    log("== E7b: critical screening (ell>=1, refined cutoff) ================")
    log(f"{'state':>6} {'lc (this work)':>16} {'lc (published)':>18} {'diff':>10}")
    for st, ell, want, lo, hi, ref, src in LC_REFS:
        lc = critical_lambda_refined(ell, want, lo, hi)
        log(f"{st:>6} {lc:>16.9f} {ref:>18.12f} {abs(lc-ref):>10.1e}  [{src}]")
    log("(s-states: threshold state infinitely diffuse; documented bias ~1e-4,")
    log(" to be handled by E(lambda)->0 extrapolation in a later pass)")
    log()


# ------------------------------------------------- resonance benchmarks --
RES_REFS = [
    # (state, ell, lam, E_r_published, Gamma/2_published, source)
    ("2p", 1, 0.22131705, 1.375843424574e-4, 1.528916191194e-5,
     "Alhaidari08 (Bylicki07: 1.3758648e-4, 1.5287316e-5)"),
    ("3p", 1, 0.11324071, 3.096253670542e-5, 4.640977622365e-6,
     "Alhaidari08 (Bylicki07: 3.0964751e-5, 4.6283598e-6)"),
    ("3d", 2, 0.092, 1.410607617484e-4, 1.170284982249e-6,
     "Alhaidari08 (Bylicki07: 1.410647e-4, 1.1700e-6)"),
]


def e7c_resonances():
    log("== E7c: resonance benchmarks (near-threshold, 12-digit refs) =======")
    thetas = np.arange(0.10, 0.46, 0.025)
    for st, ell, lam, Er_ref, hG_ref, src in RES_REFS:
        log(f"-- {st} at lambda = {lam}  [{src}]")
        for N, beta in ((250, 0.08), (300, 0.12)):
            t0 = time.perf_counter()
            B, T, spectra = theta_trajectory(lam, ell, N, beta, thetas)
            got = find_resonance(spectra, thetas,
                                 re_win=(1e-6, 2e-3), im_win=(-2e-3, 1e-7))
            dt = time.perf_counter() - t0
            if got is None:
                log(f"   N={N}, beta={beta}: NO stationary point found ({dt:.0f}s)")
                continue
            stat, Eres, th_opt, _ = got
            dE = abs(Eres.real - Er_ref)
            dG = abs(-Eres.imag - hG_ref)
            log(f"   N={N}, beta={beta}: E_r={Eres.real:.9e}  G/2={-Eres.imag:.9e}"
                f"  (stat {stat:.1e}, th*={th_opt:.2f}, {dt:.0f}s)")
            log(f"      vs published:   dE_r={dE:.1e}  dG/2={dG:.1e}")
        log()


if __name__ == "__main__":
    t0 = time.perf_counter()
    e7a_bound()
    e7b_lambdac()
    e7c_resonances()
    log(f"total wall time: {time.perf_counter()-t0:.1f} s")
    (RES / "E7_benchmark_comparison.txt").write_text("\n".join(OUT),
                                                     encoding="utf-8")
