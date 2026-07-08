# FEAST + orthonormal Laguerre basis for radial Schrödinger spectra

Reference implementation for the study *"Bound states and resonances of
screened Coulomb potentials with the FEAST contour-integration eigensolver
in an orthonormal Laguerre basis"* (working title). The hydrogen atom is
used **only as an exactly solvable validation benchmark**; the physics
targets are screened (Debye–Hückel/Yukawa) potentials and, via complex
scaling, their resonances.

## Layout

```
src/basis.py         orthonormal generalized-Laguerre basis (L2-orthonormality
                     certified numerically at build time), matrix elements for
                     kinetic, Coulomb, Yukawa and complex-scaled Hamiltonians
src/feast.py         dense FEAST eigensolver (contour integration +
                     Rayleigh-Ritz), with self-test against LAPACK
src/reference_fd.py  independent finite-difference reference solver
                     (cross-validation only)
scripts/run_validation.py     Phase-1 validation suite (E1-E4)
scripts/make_pipeline_figure.py  methodology figure
results/, figures/   generated tables (CSV/TXT) and figures (PDF+PNG)
```

## Reproduce

```
pip install -r requirements.txt
python scripts/run_validation.py
```

Runtime ≈ 4 min on a laptop. All results are deterministic (fixed seeds).

## Phase-1 validation results (2026-07-08)

| Check | Result |
|---|---|
| E1 Hydrogen ℓ=0, N=150, β=1/6 | max error 2.3×10⁻¹⁰ Ha over n=1–6 (exponential convergence) |
| E2 FEAST vs LAPACK | agreement 5.4×10⁻¹⁵ Ha, residuals < 3×10⁻¹² |
| E3 Yukawa vs independent FD | agreement at the FD discretization error (≈2×10⁻⁶ Ha), Laguerre values more accurate |
| E4 Complex scaling θ ∈ [0.1, 0.4] | bound states real and θ-stable to 1×10⁻¹⁰ Ha |

## Units and conventions

Hartree atomic units throughout (ħ = mₑ = e = 4πε₀ = 1); exact hydrogen
levels Eₙ = −1/(2n²) Ha; 1 Ha = 27.211386245981 eV (CODATA 2022).
