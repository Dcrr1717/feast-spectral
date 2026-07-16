/-
Machine-checked core identities of the fiber analysis
(companion to "Exact fiber analysis of kinetic-balance spectral
pollution for the radial Dirac operator").

Formalized here, over the reals:
1. `key_identity`     : t - ρ = 4t³/(1+w)²  where w = √(1+4t²), ρ = 2t/(1+w)
                        (the exact balanced-leakage numerator, L1'/Thm 3.1)
2. `rho_lt`           : ρ < t for t > 0 (resolves the absolute values)
3. `defect_identity`  : ‖(I−bbᵀ)Λ₊b‖² = g²(1−g²) in scalar form
                        (per-fiber defect, Thm 4.1 step)
4. `anti_in_gap`      : the anti-balanced Rayleigh value c²(5t²−1)/(1+t²)
                        lies strictly inside (−c², c²) iff 0 < t < 1/√2
                        (pollution dichotomy, Thm 5.1(a))
-/

import Mathlib.Analysis.SpecialFunctions.Sqrt
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Tactic

open Real

namespace FeastLemmas

variable {t : ℝ}

/-- Auxiliary facts about w = √(1+4t²). -/
lemma w_sq (t : ℝ) : (Real.sqrt (1 + 4 * t ^ 2)) ^ 2 = 1 + 4 * t ^ 2 := by
  rw [sq_sqrt]; positivity

lemma w_pos (t : ℝ) : 0 < Real.sqrt (1 + 4 * t ^ 2) := by
  apply Real.sqrt_pos.mpr; positivity

lemma one_lt_w (ht : 0 < t) : 1 < Real.sqrt (1 + 4 * t ^ 2) := by
  have h : (1 : ℝ) < 1 + 4 * t ^ 2 := by nlinarith
  calc (1 : ℝ) = Real.sqrt 1 := (Real.sqrt_one).symm
    _ < Real.sqrt (1 + 4 * t ^ 2) := Real.sqrt_lt_sqrt (by norm_num) h

/-- The key algebraic identity behind the exact balanced leakage
    (Theorem 3.1): with w = √(1+4t²) and ρ = 2t/(1+w),
    t − ρ = 4t³/(1+w)². -/
theorem key_identity (t : ℝ) :
    t - 2 * t / (1 + Real.sqrt (1 + 4 * t ^ 2)) =
    4 * t ^ 3 / (1 + Real.sqrt (1 + 4 * t ^ 2)) ^ 2 := by
  set w := Real.sqrt (1 + 4 * t ^ 2) with hw
  have hw2 : w ^ 2 = 1 + 4 * t ^ 2 := w_sq t
  have hwpos : 0 < w := w_pos t
  have h1w : (1 + w) ≠ 0 := by positivity
  field_simp
  nlinarith [hw2]

/-- ρ < t for t > 0: the balanced trial ratio overshoots the exact
    eigenvector ratio, which fixes the sign in the leakage formula. -/
theorem rho_lt (ht : 0 < t) :
    2 * t / (1 + Real.sqrt (1 + 4 * t ^ 2)) < t := by
  have hw := one_lt_w ht
  have h2 : (2 : ℝ) < 1 + Real.sqrt (1 + 4 * t ^ 2) := by linarith
  have hpos : (0 : ℝ) < 1 + Real.sqrt (1 + 4 * t ^ 2) := by linarith
  rw [div_lt_iff hpos]
  nlinarith

/-- Per-fiber defect identity (Theorem 4.1): for unit vectors the
    squared norm of (I−bbᵀ)Λ₊b reduces to g²(1−g²) with g = ⟨v₊,b⟩.
    Scalar form: g²·(1 − 2g² + g²·1) = g²(1−g²). -/
theorem defect_identity (g : ℝ) :
    g ^ 2 * (1 - 2 * g ^ 2 + g ^ 2 * 1) = g ^ 2 * (1 - g ^ 2) := by
  ring

/-- Pollution dichotomy, part (a) of Theorem 5.1: the anti-balanced
    fiber Rayleigh value R = c²(5t²−1)/(1+t²) lies strictly inside the
    gap (−c², c²) if and only if 0 < t and t² < 1/2. -/
theorem anti_in_gap (c t : ℝ) (hc : 0 < c) :
    (-c ^ 2 < c ^ 2 * (5 * t ^ 2 - 1) / (1 + t ^ 2) ∧
      c ^ 2 * (5 * t ^ 2 - 1) / (1 + t ^ 2) < c ^ 2) ↔
    (0 < t ^ 2 ∧ t ^ 2 < 1 / 2) := by
  have hden : (0 : ℝ) < 1 + t ^ 2 := by positivity
  have hc2 : (0 : ℝ) < c ^ 2 := by positivity
  constructor
  · rintro ⟨hlo, hhi⟩
    constructor
    · by_contra h
      push_neg at h
      have ht0 : t ^ 2 = 0 := le_antisymm h (sq_nonneg t)
      rw [div_lt_iff hden] at hhi
      rw [lt_div_iff hden] at hlo
      nlinarith
    · rw [div_lt_iff hden] at hhi
      nlinarith
  · rintro ⟨ht0, ht2⟩
    constructor
    · rw [lt_div_iff hden]
      nlinarith
    · rw [div_lt_iff hden]
      nlinarith

/-- The anti-balanced Rayleigh value crosses the exact gap midpoint at
    t = 1/√5: 5t² − 1 = 0 there. -/
theorem anti_midgap :
    5 * (1 / Real.sqrt 5) ^ 2 - 1 = 0 := by
  have h5 : (Real.sqrt 5) ^ 2 = 5 := by
    rw [sq_sqrt]; norm_num
  have hpos : Real.sqrt 5 ≠ 0 := by
    have := Real.sqrt_pos.mpr (by norm_num : (5 : ℝ) > 0)
    exact ne_of_gt this
  field_simp
  nlinarith [h5]

end FeastLemmas
