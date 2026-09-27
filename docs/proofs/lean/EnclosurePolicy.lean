import Contracts

/- Rational interval model for retained exact coefficients only. No analytic-tail,
complex-ball, floating-point implementation, or C++ refinement assertion. -/
namespace DiffExp.EnclosurePolicy

structure Interval where
  lower : Rat
  upper : Rat

def Encloses (a : Interval) (x : Rat) : Prop := a.lower ≤ x ∧ x ≤ a.upper

/-- Fixed absolute-width acceptance policy. Same tolerance applies to every strategy. -/
def Accepted (tolerance : Rat) (a : Interval) : Prop := a.upper - a.lower ≤ tolerance
instance (tolerance : Rat) (a : Interval) : Decidable (Accepted tolerance a) :=
  inferInstanceAs (Decidable (a.upper - a.lower ≤ tolerance))

def intersect (a b : Interval) : Interval where
  lower := if a.lower ≤ b.lower then b.lower else a.lower
  upper := if a.upper ≤ b.upper then a.upper else b.upper

/-- Interval intersection preserves enclosure when both implementations enclose the SAME value. -/
theorem intersection_encloses (a b : Interval) (x : Rat)
    (ha : Encloses a x) (hb : Encloses b x) : Encloses (intersect a b) x := by
  constructor
  · change (if a.lower ≤ b.lower then b.lower else a.lower) ≤ x
    split
    · exact hb.1
    · exact ha.1
  · change x ≤ (if a.upper ≤ b.upper then a.upper else b.upper)
    split
    · exact ha.2
    · exact hb.2

/-- Intersection can only raise the lower endpoint and lower the upper endpoint. -/
theorem intersection_bounds (a b : Interval) :
    a.lower ≤ (intersect a b).lower ∧ (intersect a b).upper ≤ a.upper := by
  constructor
  · change a.lower ≤ (if a.lower ≤ b.lower then b.lower else a.lower)
    split
    · assumption
    · exact Rat.le_refl
  · change (if a.upper ≤ b.upper then a.upper else b.upper) ≤ a.upper
    split
    · exact Rat.le_refl
    · rename_i h
      rcases Rat.le_total (a := a.upper) (b := b.upper) with hab | hba
      · exact False.elim (h hab)
      · exact hba

theorem sub_mono (u' u l l' : Rat) (hu : u' ≤ u) (hl : l ≤ l') :
    u' - l' ≤ u - l := by
  rw [Rat.sub_eq_add_neg, Rat.sub_eq_add_neg]
  exact Rat.le_trans (Rat.add_le_add_right.mpr hu)
    (Rat.add_le_add_left.mpr (Rat.neg_le_neg hl))

/-- A fixed absolute-width policy cannot be invalidated by intersection. -/
theorem intersection_accepted (a b : Interval) (tolerance : Rat)
    (ha : Accepted tolerance a) : Accepted tolerance (intersect a b) := by
  have bounds := intersection_bounds a b
  exact Rat.le_trans (sub_mono _ _ _ _ bounds.2 bounds.1) ha

/-- A lazy fallback receives Unit only in the rejected branch. -/
def select (tolerance : Rat) (candidate : Interval) (fallback : Unit → Interval) : Interval :=
  if Accepted tolerance candidate then candidate else fallback ()

/-- For an accepted candidate the selected interval is independent of every fallback result. -/
theorem accepted_fallback_independent (tolerance : Rat) (candidate : Interval)
    (f g : Unit → Interval) (ha : Accepted tolerance candidate) :
    select tolerance candidate f = candidate ∧ select tolerance candidate g = candidate := by
  simp [select, ha]

/-- Complete conditional policy correctness, including the branch that needs recovery. -/
theorem select_sound (tolerance x : Rat) (candidate : Interval) (fallback : Unit → Interval)
    (hc : Encloses candidate x)
    (hf : ¬ Accepted tolerance candidate →
      Encloses (fallback ()) x ∧ Accepted tolerance (fallback ())) :
    Encloses (select tolerance candidate fallback) x ∧
      Accepted tolerance (select tolerance candidate fallback) := by
  by_cases ha : Accepted tolerance candidate
  · simpa [select, ha] using And.intro hc ha
  · simpa [select, ha] using hf ha

/-- If the candidate already meets the contract, both stopping and intersecting are sound.
This establishes optionality for correctness, not equality or identical error bars. -/
theorem second_enclosure_optional (tolerance x : Rat) (candidate second : Interval)
    (hc : Encloses candidate x) (ha : Accepted tolerance candidate)
    (hs : Encloses second x) :
    (Encloses candidate x ∧ Accepted tolerance candidate) ∧
    (Encloses (intersect candidate second) x ∧ Accepted tolerance (intersect candidate second)) := by
  exact ⟨⟨hc, ha⟩, intersection_encloses candidate second x hc hs,
    intersection_accepted candidate second tolerance ha⟩

/-- Concrete non-identical outputs: both satisfy tolerance 4 and enclose retained value 1. -/
def candidateExample : Interval := ⟨0, 4⟩
def secondExample : Interval := ⟨1, 2⟩

theorem example_candidate_encloses : Encloses candidateExample 1 := by
  unfold Encloses
  decide +kernel

theorem example_candidate_accepted : Accepted 4 candidateExample := by
  decide +kernel

theorem example_second_encloses : Encloses secondExample 1 := by
  unfold Encloses
  decide +kernel

end DiffExp.EnclosurePolicy

#print axioms DiffExp.EnclosurePolicy.intersection_encloses
#print axioms DiffExp.EnclosurePolicy.intersection_bounds
#print axioms DiffExp.EnclosurePolicy.intersection_accepted
#print axioms DiffExp.EnclosurePolicy.accepted_fallback_independent
#print axioms DiffExp.EnclosurePolicy.select_sound
#print axioms DiffExp.EnclosurePolicy.second_enclosure_optional
#print axioms DiffExp.EnclosurePolicy.example_candidate_encloses
#print axioms DiffExp.EnclosurePolicy.example_candidate_accepted
#print axioms DiffExp.EnclosurePolicy.example_second_encloses
