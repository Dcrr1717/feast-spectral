"""Independent finite-difference reference solver (cross-validation only).

Uniform radial grid r_i = i*dr, i = 1..N (the choice r_min = dr imposes the
Dirichlet condition u(0) = 0 implicitly, which is what makes the standard
radial FD scheme second-order accurate for the Coulomb problem). Solved with
LAPACK's direct tridiagonal eigensolver -- completely independent from the
Laguerre-basis pipeline, so agreement between the two is a genuine
cross-method validation.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import eigh_tridiagonal


def fd_levels(potential, ell: int = 0, k: int = 6,
              dr: float = 0.005, rmax: float = 300.0) -> np.ndarray:
    """Lowest k eigenvalues (Hartree) of
    -1/2 u'' + [ell(ell+1)/(2 r^2) + potential(r)] u = E u,  u(0)=u(rmax)=0.
    """
    n = int(round(rmax / dr))
    r = dr * np.arange(1, n + 1)
    diag = 1.0 / dr**2 + potential(r)
    if ell:
        diag = diag + 0.5 * ell * (ell + 1) / r**2
    off = np.full(n - 1, -0.5 / dr**2)
    return eigh_tridiagonal(diag, off, select='i',
                            select_range=(0, k - 1), eigvals_only=True)


def coulomb(Z: float = 1.0):
    return lambda r: -Z / r


def yukawa(Z: float = 1.0, lam: float = 0.0):
    return lambda r: -Z * np.exp(-lam * r) / r
