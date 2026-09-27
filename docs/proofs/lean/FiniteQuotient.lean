import Contracts

namespace DiffExp.FiniteQuotient

/-- One nonconstant denominator monomial. Duplicate monomials are summed. -/
structure Term where
  a : Nat
  b : Nat
  coefficient : Rat
  nonconstant : a ≠ 0 ∨ b ≠ 0

/-- Shifted epsilon window has positive width; Taylor order is the quotient index. -/
def predecessor (width n : Nat) (t : Term) : Nat :=
  (n / width - t.a) * width + (n % width - t.b)

def Admitted (width n : Nat) (t : Term) : Prop := t.a ≤ n / width ∧ t.b ≤ n % width
instance (width n : Nat) (t : Term) : Decidable (Admitted width n t) :=
  inferInstanceAs (Decidable (_ ∧ _))

theorem predecessor_lt (width n : Nat) (hw : 0 < width) (t : Term)
    (ht : Admitted width n t) : predecessor width n t < n := by
  have h := flattened_predecessor (n / width) (n % width) t.a t.b width
    (Nat.mod_lt n hw) ht.1 ht.2 t.nonconstant
  have hn : n / width * width + n % width = n := by
    rw [Nat.mul_comm, Nat.div_add_mod]
  simpa [predecessor, hn] using h

/-- Literal finite polynomial convolution excluding the separately stored q00 pivot.
Out-of-bounds monomials contribute zero. No expanded rational coefficient AST is used. -/
def offConvolution (width : Nat) (support : List Term) (n : Nat) (w : Nat → Rat) : Rat :=
  match support with
  | [] => 0
  | t :: ts =>
    (if Admitted width n t then t.coefficient * w (predecessor width n t) else 0) +
      offConvolution width ts n w

theorem offConvolution_local (width : Nat) (hw : 0 < width) (support : List Term)
    (n : Nat) (w z : Nat → Rat) (heq : ∀ j, j < n → w j = z j) :
    offConvolution width support n w = offConvolution width support n z := by
  induction support with
  | nil => rfl
  | cons t ts ih =>
    simp only [offConvolution]
    by_cases ht : Admitted width n t
    · rw [ite_eq_left ht, ite_eq_left ht, heq _ (predecessor_lt width n hw t ht), ih]
    · rw [ite_eq_right ht, ite_eq_right ht, ih]

def earlier (n : Nat) (prev : Fin n → Rat) (j : Nat) : Rat :=
  if hj : j < n then prev ⟨j, hj⟩ else 0

def tail (width : Nat) (support : List Term) (n : Nat) (prev : Fin n → Rat) : Rat :=
  offConvolution width support n (earlier n prev)

theorem tail_is_convolution (width : Nat) (hw : 0 < width) (support : List Term)
    (n : Nat) (w : Nat → Rat) :
    tail width support n (fun i => w i.val) = offConvolution width support n w := by
  apply offConvolution_local width hw support n
  intro j hj
  simp [earlier, hj]

def quotient (width : Nat) (support : List Term) (q : Rat) (y : Nat → Rat) : Nat → Rat :=
  solve (fun n prev => (y n - tail width support n prev) / q)

def convolution (width : Nat) (support : List Term) (q : Rat)
    (n : Nat) (w : Nat → Rat) : Rat := q * w n + offConvolution width support n w

/-- Actual finite-support QW=Y, for every cell of the fixed finite epsilon window. -/
theorem quotient_equation (width : Nat) (hw : 0 < width) (support : List Term)
    (q : Rat) (hq : q ≠ 0) (y : Nat → Rat) (n : Nat) :
    convolution width support q n (quotient width support q y) = y n := by
  have h := quotient_stream_equation q hq y (tail width support) n
  change q * quotient width support q y n +
    tail width support n (fun i => quotient width support q y i.val) = y n at h
  rw [tail_is_convolution width hw] at h
  exact h

theorem cell_rearrange (q value off y : Rat) (hq : q ≠ 0)
    (h : q * value + off = y) : value = (y - off) / q := by
  have he : q * value = y - off := by
    rw [← h, Rat.add_sub_cancel]
  calc
    value = (q * value) / q := by rw [Rat.mul_comm q, Rat.mul_div_cancel hq]
    _ = (y - off) / q := by rw [he]

/-- Even if equations are known only below a finite cutoff, their coefficients are unique.
Taking cutoff = TaylorRows * width gives a rectangular retained window. -/
theorem quotient_unique_prefix (width : Nat) (hw : 0 < width) (support : List Term)
    (q : Rat) (hq : q ≠ 0) (y w : Nat → Rat) (cutoff : Nat)
    (heq : ∀ n, n < cutoff → convolution width support q n w = y n) :
    ∀ n, n < cutoff → w n = quotient width support q y n := by
  intro n
  induction n using Nat.strongRecOn with
  | ind n ih =>
    intro hn
    have hwcell := cell_rearrange q (w n) (offConvolution width support n w) (y n)
      hq (heq n hn)
    have hqcell := cell_rearrange q (quotient width support q y n)
      (offConvolution width support n (quotient width support q y)) (y n)
      hq (quotient_equation width hw support q hq y n)
    rw [hwcell, hqcell]
    congr 1
    congr 1
    apply offConvolution_local width hw support n
    intro j hj
    exact ih j hj (Nat.lt_trans hj hn)

/-- Unknown source coefficients above the retained cutoff cannot alter its solution prefix. -/
theorem source_prefix_independent (width : Nat) (hw : 0 < width) (support : List Term)
    (q : Rat) (hq : q ≠ 0) (y z : Nat → Rat) (cutoff : Nat)
    (hyz : ∀ n, n < cutoff → y n = z n) :
    ∀ n, n < cutoff → quotient width support q y n = quotient width support q z n := by
  apply quotient_unique_prefix width hw support q hq z
    (quotient width support q y) cutoff
  intro n hn
  rw [quotient_equation width hw support q hq y n, hyz n hn]

end DiffExp.FiniteQuotient

#print axioms DiffExp.FiniteQuotient.predecessor_lt
#print axioms DiffExp.FiniteQuotient.offConvolution_local
#print axioms DiffExp.FiniteQuotient.quotient_equation
#print axioms DiffExp.FiniteQuotient.quotient_unique_prefix

#print axioms DiffExp.FiniteQuotient.source_prefix_independent
