# Machine-checked lemmas (Lean 4 + Mathlib)

`FiberLemmas.lean` contains machine-verified proofs of the core
algebraic identities of the paper:

| Lean theorem | Paper statement |
|---|---|
| `key_identity` | t − ρ = 4t³/(1+w)² (Thm 3.1, balanced leakage) |
| `rho_lt` | ρ < t for t > 0 (sign resolution in Thm 3.1) |
| `defect_identity` | ‖(I−bbᵀ)Λ₊b‖² = g²(1−g²) (Thm 4.1, per-fiber defect) |
| `anti_in_gap` | R_anti ∈ (−c²,c²) ⇔ 0 < t² < 1/2 (Thm 5.1(a)) |
| `anti_midgap` | mid-gap crossing at t = 1/√5 (Thm 5.1(a)) |

## Reproducing the check

```
# install elan (Lean toolchain manager), then:
lake new feastlemmas math
cd feastlemmas
lake exe cache get          # prebuilt Mathlib (several GB)
cp <this>/FiberLemmas.lean Feastlemmas/
lake build Feastlemmas.FiberLemmas
```

Verified 2026-07-16 on Lean stable + current Mathlib:
`Build completed successfully (3069 jobs)`.
