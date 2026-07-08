"""FEAST eigensolver (contour integration + Rayleigh-Ritz), dense reference
implementation.

Algorithm (Polizzi, Phys. Rev. B 79, 115112 (2009)): the spectral projector
onto the eigenspace enclosed by a contour C is

    P = (1 / 2*pi*i) * oint_C (z*I - H)^(-1) dz .

FEAST approximates Q = P @ Y for a random block Y by an n_contour-point
Gauss-Legendre quadrature on a circle enclosing the search interval
[emin, emax], then solves the reduced Rayleigh-Ritz problem in span(Q) and
iterates (subspace iteration with a rational filter).

This implementation is deliberately dense and transparent (the paper's
reference solver); it is validated against scipy.linalg.eigh in the test
scripts. For real-symmetric H the conjugate symmetry of the contour is
exploited implicitly by taking the real part of Q.
"""

from __future__ import annotations

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import lu_factor, lu_solve, eigh, qr


def feast_interval(H: np.ndarray, emin: float, emax: float, m0: int,
                   n_contour: int = 16, max_iter: int = 10,
                   tol: float = 1e-11, seed: int = 12345):
    """Find all eigenpairs of the symmetric matrix H inside [emin, emax].

    Parameters
    ----------
    H         : (n, n) real symmetric matrix
    emin,emax : search interval (must contain the wanted eigenvalues and
                its endpoints must not coincide with an eigenvalue)
    m0        : subspace size, must exceed the number of eigenvalues inside
    n_contour : quadrature nodes on the full circular contour
    tol       : maximum absolute residual ||H x - lambda x||_2

    Returns
    -------
    evals, evecs, residuals, n_iterations, info(dict)
    """
    n = H.shape[0]
    rng = np.random.default_rng(seed)
    Y = rng.standard_normal((n, m0))

    c = 0.5 * (emin + emax)
    R = 0.5 * (emax - emin)
    # Gauss-Legendre nodes in the angle variable, full circle (0, 2*pi).
    # Even n_contour keeps all nodes off the real axis.
    t, wt = leggauss(n_contour)
    ang = np.pi * (t + 1.0)
    zs = c + R * np.exp(1j * ang)
    # (1/2*pi*i) * dz  with  dz = i R e^{i ang} d(ang),  d(ang) = pi dt
    wz = wt * (np.pi * 1j * R * np.exp(1j * ang)) / (2j * np.pi)

    # LU factorizations are reused across refinement iterations
    lus = [lu_factor(z * np.eye(n) - H) for z in zs]

    history = []
    evals = np.empty(0)
    X = np.empty((n, 0))
    res = np.empty(0)
    for it in range(1, max_iter + 1):
        Q = np.zeros((n, m0), dtype=complex)
        for lu, w in zip(lus, wz):
            Q += w * lu_solve(lu, Y)
        Q = np.real(Q)  # H real symmetric: projector is real

        Qo, _ = qr(Q, mode='economic')
        Hr = Qo.T @ H @ Qo
        lam, W = eigh(0.5 * (Hr + Hr.T))
        Xall = Qo @ W

        inside = (lam > emin) & (lam < emax)
        evals, X = lam[inside], Xall[:, inside]
        res = np.linalg.norm(H @ X - X * evals, axis=0)
        history.append(float(res.max()) if res.size else np.nan)
        if res.size and res.max() < tol:
            return evals, X, res, it, {"residual_history": history}
        Y = Xall  # rational-filtered subspace becomes next block

    return evals, X, res, max_iter, {"residual_history": history}


def feast_region(H: np.ndarray, center: complex, radius: float, m0: int,
                 n_contour: int = 24, max_iter: int = 12, tol: float = 1e-9,
                 seed: int = 7):
    """FEAST for a general/complex matrix: eigenvalues inside a circle.

    Finds every eigenvalue of H enclosed by the circle |z - center| < radius,
    using the same contour-integration projector as `feast_interval` but on a
    circle placed anywhere in the complex plane. This is what turns FEAST into
    a *resonance* solver: a small circle around a complex eigenvalue of the
    complex-scaled Hamiltonian isolates the resonance E = E_r - i*Gamma/2 while
    excluding the rotated continuum.

    No real-part projection is taken (H is complex). The reduced Rayleigh-Ritz
    problem is a small non-Hermitian eigenproblem solved with numpy.linalg.eig.

    Returns
    -------
    evals, evecs, residuals, n_iterations, info(dict)
    """
    n = H.shape[0]
    rng = np.random.default_rng(seed)
    Y = rng.standard_normal((n, m0)) + 1j * rng.standard_normal((n, m0))

    t, wt = leggauss(n_contour)
    ang = np.pi * (t + 1.0)
    zs = center + radius * np.exp(1j * ang)
    wz = wt * (np.pi * 1j * radius * np.exp(1j * ang)) / (2j * np.pi)

    eye = np.eye(n)
    lus = [lu_factor(z * eye - H) for z in zs]

    history = []
    evals = np.empty(0, dtype=complex)
    X = np.empty((n, 0), dtype=complex)
    res = np.array([np.inf])
    for it in range(1, max_iter + 1):
        Q = np.zeros((n, m0), dtype=complex)
        for lu, w in zip(lus, wz):
            Q += w * lu_solve(lu, Y)
        W, _ = qr(Q, mode='economic')
        Hr = W.conj().T @ H @ W
        lam, Vr = np.linalg.eig(Hr)
        Xall = W @ Vr
        inside = np.abs(lam - center) < radius
        evals, X = lam[inside], Xall[:, inside]
        if evals.size:
            res = np.linalg.norm(H @ X - X * evals, axis=0)
            history.append(float(res.max()))
            if res.max() < tol:
                return evals, X, res, it, {"residual_history": history}
        else:
            history.append(np.inf)
        Y = Xall
    return evals, X, res, max_iter, {"residual_history": history}


def _self_test():
    """Projector orientation / correctness check on a random symmetric matrix."""
    rng = np.random.default_rng(0)
    A = rng.standard_normal((60, 60))
    A = 0.5 * (A + A.T)
    exact = np.linalg.eigvalsh(A)
    emin, emax = exact[10] - 1e-3, exact[15] + 1e-3
    lam, X, r, it, _ = feast_interval(A, emin, emax, m0=12)
    assert lam.size == 6, f"expected 6 eigenvalues, got {lam.size}"
    assert np.allclose(np.sort(lam), exact[10:16], atol=1e-10)
    assert r.max() < 1e-9
    return True


if __name__ == "__main__":
    _self_test()
    print("feast.py self-test passed.")
