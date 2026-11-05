import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

-- Exact scalar recurrence obtained from the binary probability generating function.
theorem binary_pgf_identity (e q : ℝ) :
    (1 + e)^2 / 2 + (1 - q)^2 / 2 - 1 =
      e + e^2 / 2 - q + q^2 / 2 := by
  ring

-- The one-step barrier estimate; r = 1/n and rnext = 1/(n+1).
theorem barrier_step (e q enext r rnext : ℝ)
    (he : 0 ≤ e) (hr : 0 ≤ r)
    (hb : e ≤ (3 / 2 : ℝ) * r)
    (hf : 3 * r^2 ≤ q - q^2 / 2)
    (hrec : enext = e + e^2 / 2 - q + q^2 / 2)
    (hrnext : r - r^2 ≤ rnext) :
    enext - (3 / 2 : ℝ) * rnext ≤
      e - (3 / 2 : ℝ) * r - (3 / 8 : ℝ) * r^2 := by
  have hprod : 0 ≤ ((3 / 2 : ℝ) * r - e) * ((3 / 2 : ℝ) * r + e) :=
    mul_nonneg (by linarith) (by linarith)
  nlinarith
