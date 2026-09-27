import Contracts

/- Finite rational coefficient identities. No floating-point, ACB, or C++ refinement claim. -/
namespace DiffExp.UncertaintyWindow

def sumN (n : Nat) (f : Nat → Rat) : Rat :=
  match n with
  | 0 => 0
  | n + 1 => sumN n f + f n

theorem sumN_congr (n : Nat) (f g : Nat → Rat)
    (h : ∀ e, e < n → f e = g e) : sumN n f = sumN n g := by
  induction n with
  | zero => rfl
  | succ n ih =>
    simp only [sumN]
    rw [ih (fun e he => h e (Nat.lt_trans he (Nat.lt_succ_self n))), h n (Nat.lt_succ_self n)]

/-- A finite sum can stop early if all omitted terms are proved zero. -/
theorem sumN_shorten (n bound : Nat) (f : Nat → Rat) (hb : bound ≤ n)
    (hz : ∀ e, bound ≤ e → e < n → f e = 0) : sumN n f = sumN bound f := by
  induction n with
  | zero =>
    have h : bound = 0 := by omega
    subst bound
    rfl
  | succ n ih =>
    by_cases h : bound = n + 1
    · rw [h]
    · have hbn : bound ≤ n := by omega
      rw [sumN, ih hbn (fun e hbe hen => hz e hbe (by omega)), hz n hbn (by omega)]
      exact Rat.add_zero _

/-- Guard before subtraction: an exponent above the output degree never reads delta. -/
def convolutionTerm (T delta : Nat → Rat) (k e : Nat) : Rat :=
  if e ≤ k then T e * delta (k - e) else 0

def convolution (width : Nat) (T delta : Nat → Rat) (k : Nat) : Rat :=
  sumN width (convolutionTerm T delta k)

/-- Reads map coefficients only below width-low, which is zero if low ≥ width. -/
def truncatedConvolution (width low : Nat) (T delta : Nat → Rat) (k : Nat) : Rat :=
  sumN (width - low) (convolutionTerm T delta k)

/-- The key integer dependency implication, including low ≥ width and low = 0. -/
theorem omitted_input_below_low (width low k e : Nat) (hk : k < width)
    (he : e ≤ k) (homit : width - low ≤ e) : k - e < low := by omega

/-- Every term involving an omitted high map coefficient vanishes exactly. -/
theorem omitted_term_zero (width low : Nat) (T delta : Nat → Rat) (k e : Nat)
    (hk : k < width) (homit : width - low ≤ e)
    (hzero : ∀ d, d < low → delta d = 0) : convolutionTerm T delta k e = 0 := by
  unfold convolutionTerm
  split
  · rename_i he
    rw [hzero (k - e) (omitted_input_below_low width low k e hk he homit)]
    exact Rat.mul_zero _
  · rfl

/-- Exact equality through every retained output degree, with no higher map reads. -/
theorem truncate_uncertainty_map (width low : Nat) (T delta : Nat → Rat)
    (hzero : ∀ d, d < low → delta d = 0) (k : Nat) (hk : k < width) :
    convolution width T delta k = truncatedConvolution width low T delta k := by
  apply sumN_shorten width (width - low) (convolutionTerm T delta k) (Nat.sub_le _ _)
  intro e he _
  exact omitted_term_zero width low T delta k e hk he hzero

/-- Completely arbitrary changes to unknown higher map coefficients cannot affect the result. -/
theorem unknown_map_independent (width low : Nat) (T U delta : Nat → Rat)
    (hzero : ∀ d, d < low → delta d = 0)
    (hmap : ∀ e, e < width - low → T e = U e) (k : Nat) (hk : k < width) :
    convolution width T delta k = convolution width U delta k := by
  rw [truncate_uncertainty_map width low T delta hzero k hk,
      truncate_uncertainty_map width low U delta hzero k hk]
  apply sumN_congr
  intro e he
  unfold convolutionTerm
  rw [hmap e he]

/-- One output row of a finite matrix acting on epsilon-series columns. -/
def matrixConvolution (columns width : Nat) (T delta : Nat → Nat → Rat) (k : Nat) : Rat :=
  sumN columns (fun j => convolution width (T j) (delta j) k)

def truncatedMatrixConvolution (columns width : Nat) (low : Nat → Nat)
    (T delta : Nat → Nat → Rat) (k : Nat) : Rat :=
  sumN columns (fun j => truncatedConvolution width (low j) (T j) (delta j) k)

/-- Per-column markers suffice only when EVERY incoming uncertainty below each marker is zero. -/
theorem truncate_each_column (columns width : Nat) (low : Nat → Nat)
    (T delta : Nat → Nat → Rat)
    (hzero : ∀ j, j < columns → ∀ d, d < low j → delta j d = 0)
    (k : Nat) (hk : k < width) :
    matrixConvolution columns width T delta k =
      truncatedMatrixConvolution columns width low T delta k := by
  apply sumN_congr
  intro j hj
  exact truncate_uncertainty_map width (low j) (T j) (delta j) (hzero j hj) k hk

theorem sumN_add (n : Nat) (f g : Nat → Rat) :
    sumN n (fun e => f e + g e) = sumN n f + sumN n g := by
  induction n with
  | zero => exact (Rat.add_zero 0).symm
  | succ n ih =>
    simp only [sumN, ih]
    simp only [Rat.add_assoc, Rat.add_left_comm, Rat.add_comm]

theorem convolution_add (width : Nat) (T m delta : Nat → Rat) (k : Nat) :
    convolution width T (fun d => m d + delta d) k =
      convolution width T m k + convolution width T delta k := by
  unfold convolution
  rw [← sumN_add]
  apply sumN_congr
  intro e _
  unfold convolutionTerm
  split
  · exact Rat.mul_add _ _ _
  · exact (Rat.add_zero 0).symm

theorem matrixConvolution_add (columns width : Nat) (T m delta : Nat → Nat → Rat) (k : Nat) :
    matrixConvolution columns width T (fun j d => m j d + delta j d) k =
      matrixConvolution columns width T m k + matrixConvolution columns width T delta k := by
  unfold matrixConvolution
  rw [← sumN_add]
  apply sumN_congr
  intro j _
  exact convolution_add width (T j) (m j) (delta j) k

/-- Affine forcing belongs to the midpoint solve once; it does not recur in the homogeneous map. -/
theorem affine_decomposition (columns width : Nat) (T m delta : Nat → Nat → Rat)
    (forcing : Rat) (k : Nat) :
    matrixConvolution columns width T (fun j d => m j d + delta j d) k + forcing =
      (matrixConvolution columns width T m k + forcing) +
        matrixConvolution columns width T delta k := by
  rw [matrixConvolution_add]
  simp only [Rat.add_assoc, Rat.add_comm]

/-- The combined practical contract: full midpoint/forcing, shorter homogeneous uncertainty maps. -/
theorem affine_with_truncated_uncertainty (columns width : Nat) (low : Nat → Nat)
    (T m delta : Nat → Nat → Rat) (forcing : Rat)
    (hzero : ∀ j, j < columns → ∀ d, d < low j → delta j d = 0)
    (k : Nat) (hk : k < width) :
    matrixConvolution columns width T (fun j d => m j d + delta j d) k + forcing =
      (matrixConvolution columns width T m k + forcing) +
        truncatedMatrixConvolution columns width low T delta k := by
  rw [affine_decomposition, truncate_each_column columns width low T delta hzero k hk]

end DiffExp.UncertaintyWindow

#print axioms DiffExp.UncertaintyWindow.omitted_input_below_low
#print axioms DiffExp.UncertaintyWindow.truncate_uncertainty_map
#print axioms DiffExp.UncertaintyWindow.unknown_map_independent
#print axioms DiffExp.UncertaintyWindow.truncate_each_column
#print axioms DiffExp.UncertaintyWindow.affine_decomposition
#print axioms DiffExp.UncertaintyWindow.affine_with_truncated_uncertainty
