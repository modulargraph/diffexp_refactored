import Std
import Init.Data.Rat.Lemmas

/- Small structural contracts. No claim of C++ refinement or analytic convergence. -/
namespace DiffExp

/-- Any triangular recurrence over a natural-number schedule has a solution. -/
def solve {α : Type} (step : (n : Nat) → (Fin n → α) → α) (n : Nat) : α :=
  step n (fun i => solve step i.val)
termination_by n

 theorem solve_equation {α : Type} (step : (n : Nat) → (Fin n → α) → α) (n : Nat) :
    solve step n = step n (fun i => solve step i.val) := by
  rw [solve]

/-- No extra boundary constants can be chosen for an algebraic auxiliary stream. -/
theorem solve_unique {α : Type} (step : (n : Nat) → (Fin n → α) → α)
    (w : Nat → α) (hw : ∀ n, w n = step n (fun i => w i.val)) :
    ∀ n, w n = solve step n := by
  intro n
  induction n using Nat.strongRecOn with
  | ind n ih =>
    rw [hw n, solve_equation]
    congr 1
    funext i
    exact ih i.val i.isLt

/-- Exact rational division discharges one QW=Y coefficient equation. -/
theorem quotient_cell (q y tail : Rat) (hq : q ≠ 0) :
    q * ((y - tail) / q) + tail = y := by
  rw [Rat.mul_comm q, Rat.div_mul_cancel hq, Rat.sub_add_cancel]

/-- The constructed triangular quotient stream satisfies every coefficient equation.
`tail` must separately be connected to the off-diagonal polynomial convolution. -/
theorem quotient_stream_equation (q : Rat) (hq : q ≠ 0) (y : Nat → Rat)
    (tail : (n : Nat) → (Fin n → Rat) → Rat) (n : Nat) :
    let w := solve (fun m prev => (y m - tail m prev) / q)
    q * w n + tail n (fun i => w i.val) = y n := by
  dsimp only
  rw [solve_equation]
  exact quotient_cell q (y n) _ hq

/-- Every nonconstant denominator monomial reads an earlier lexicographic cell. -/
theorem denominator_predecessor (n k a b : Nat) (ha : a ≤ n) (hb : b ≤ k)
    (h : a ≠ 0 ∨ b ≠ 0) :
    n - a < n ∨ (n - a = n ∧ k - b < k) := by
  omega

/-- Flattening a finite epsilon window gives a valid natural-number schedule. -/
theorem flattened_predecessor (n k a b width : Nat)
    (hk : k < width) (ha : a ≤ n) (hb : b ≤ k)
    (h : a ≠ 0 ∨ b ≠ 0) :
    (n - a) * width + (k - b) < n * width + k := by
  rcases denominator_predecessor n k a b ha hb h with hn | ⟨hn, hb'⟩
  · have hrow : (n - a + 1) * width ≤ n * width :=
      Nat.mul_le_mul_right width hn
    rw [Nat.add_mul] at hrow
    omega
  · rw [hn]
    omega

/-- Signed epsilon indices: demanded product terms need source orders at most H-v. -/
theorem epsilon_demand (r s v H : Int) (hr : v ≤ r) (ht : r + s ≤ H) :
    s ≤ H - v := by omega

/-- A finite path in a fixed graph, with additive natural-number edge cost. -/
inductive Path (edge : Nat → Nat → Nat → Prop) : Nat → Nat → Nat → Prop
  | nil (v : Nat) : Path edge v v 0
  | cons {u v w a b : Nat} : edge u v a → Path edge v w b → Path edge u w (a + b)

/-- Local Bellman inequalities provide a lower bound for every allowed path. -/
theorem potential_bound {edge : Nat → Nat → Nat → Prop} (p : Nat → Nat)
    (h : ∀ u v c, edge u v c → p u ≤ c + p v)
    {u v cost : Nat} (path : Path edge u v cost) : p u ≤ cost + p v := by
  induction path with
  | nil v => omega
  | @cons u v w a b e rest ih =>
    have hh := h u v a e
    omega

/-- A witness attaining the certified bound is optimal in the fixed graph. -/
theorem certified_optimal {edge : Nat → Nat → Nat → Prop} (p : Nat → Nat)
    (h : ∀ u v c, edge u v c → p u ≤ c + p v)
    {start finish chosen : Nat} (_witness : Path edge start finish chosen)
    (hfinish : p finish = 0) (hchosen : chosen = p start) :
    ∀ cost, Path edge start finish cost → chosen ≤ cost := by
  intro cost path
  have hb := potential_bound p h path
  omega

/-- Certified local transitions preserve a state invariant across a path. -/
theorem path_invariant {edge : Nat → Nat → Nat → Prop} (safe : Nat → Prop)
    (h : ∀ u v c, edge u v c → safe u → safe v)
    {u v cost : Nat} (path : Path edge u v cost) (hu : safe u) : safe v := by
  induction path with
  | nil => exact hu
  | @cons u v w a b e rest ih => exact ih (h u v a e hu)

end DiffExp

#print axioms DiffExp.solve_equation
#print axioms DiffExp.solve_unique
#print axioms DiffExp.denominator_predecessor
#print axioms DiffExp.flattened_predecessor
#print axioms DiffExp.epsilon_demand
#print axioms DiffExp.potential_bound
#print axioms DiffExp.certified_optimal
#print axioms DiffExp.path_invariant

#print axioms DiffExp.quotient_cell
#print axioms DiffExp.quotient_stream_equation
