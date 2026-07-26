# Paper F (flagship) — Relativistic screened Stark: Dirac + Debye + field

Gap scan #2 (GAP_SCAN_2026H2.md). The uncomputed TRIPLE: Dirac +
Yukawa/Debye screening + static electric field, fine-structure-resolved
Stark resonance positions AND widths over (Z, lambda, F). Target: PRA.

## Validation ladder (all published, gifts)
1. Ivanov-Ho, PRA 69, 023407 (2004): Dirac Stark Z=1,10, no screening.
2. Maltsev-Tumakov-Popov-Shabaev, PRA 108, 042813 (2023),
   arXiv:2305.07610: Dirac Stark Z=1-82, 1s/2s/2p1/2/2p3/2, positions
   AND widths; relativistic width suppression up to ~1 order. NO
   screening. THE benchmark to reproduce before switching screening on.
3. Our own Paper D (nonrel Stark-Debye) as the c->inf / low-Z limit.
4. Field-free per-channel: exact Dirac-Coulomb (Sommerfeld) + Paper B
   Yukawa resonances.

## New machinery (the ~1 week piece)
- Relativistic dipole coupling between kappa-channels. z = r cos(theta)
  multiplies all 4 spinor components; couples channel kappa to kappa'
  via the angular factor <Omega_{kappa' m}|cos theta|Omega_{kappa m}>
  for the LARGE component and <Omega_{-kappa' m}|cos theta|Omega_{-kappa m}>
  for the SMALL component (small uses opposite-kappa spin-angular fn).
  Computed via Wigner 3j (sympy), validated against the nonrel Gaunt
  limit (l+1)/sqrt((2l+1)(2l+3)).
- Shared radial grid across kappa-channels (cross-channel <r> needs a
  common grid). Each channel keeps its r^gamma prefactor; use a shared
  high-order grid, certify field-free per-channel vs Sommerfeld.
- Negative-continuum guard under complex scaling in a coupled finite
  basis (kinetic balance; our SVD pollution diagnostics = selling pt).

## Phases
- F0  angular coupling module + kappa-ladder + nonrel-limit gate
- F1  shared-grid coupled engine; field-free per-channel cert vs
      Sommerfeld (all kappa in the ladder to machine-ish precision)
- F2  Stark limit (lambda=0): reproduce Ivanov-Ho Z=1,10 then Maltsev
      Z=1-82 (1s, 2s, 2p1/2, 2p3/2) positions AND widths
- F3  switch screening on: (Z, lambda, F) resonance map, fine-structure
      resolved; relativistic width suppression vs screening; 2p1/2-2p3/2
      Stark width splitting under screening; critical (lambda,F) curves
- F4  independent cross-check (shooting or second discretization),
      adversarial verification, manuscript

## Physics headline candidates
- Does Maltsev's ~1-order relativistic width suppression survive /
  strengthen as screening pushes the level toward threshold?
- 2p1/2 vs 2p3/2 Stark-resonance WIDTH splitting under screening (the
  relativistic echo of Paper B's spin-symmetry-breaking result).
- Competition: screening lifts level -> fatter width; relativity binds
  -> narrower. Map the crossover in (Z, lambda, F).

## m sector
Start m=+1/2 (lowest). kappa-ladder for m=1/2:
kappa in {-1(s1/2), +1(p1/2), -2(p3/2), +2(d3/2), -3(d5/2), ...}.
The field couples adjacent-l with dipole selection.
