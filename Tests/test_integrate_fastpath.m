(* Differential tests for Integrate`Private`contractIntegralTerms.
   The reference deliberately uses the old ESNew/ESTimes/ESAdd reduction so
   these tests lock both coefficients and honest-window intersections. *)
repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];
Catch[DiffExp2`Config`LoadConfiguration[{"WorkingPrecision" -> 100}],
  "DiffExp2Error"];

passed = 0; failed = 0;
assert[label_String, cond_] := If[TrueQ[cond],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];

esNew = DiffExp2`EpsSeries`ESNew;
esZero = DiffExp2`EpsSeries`ESZero;
esAdd = DiffExp2`EpsSeries`ESAdd;
esTimes = DiffExp2`EpsSeries`ESTimes;
esScale = DiffExp2`EpsSeries`ESScale;
esShift = DiffExp2`EpsSeries`ESShift;
esInv = DiffExp2`EpsSeries`ESInvert;
esSameQ = DiffExp2`EpsSeries`ESSameQ;
fastContract = DiffExp2`Integrate`Private`contractIntegralTerms;

referenceContract[terms_List, kmin_Integer, kmax_Integer,
    ncomp_Integer] := If[terms === {}, Table[esZero[kmax], {ncomp}],
  Table[Module[{products},
    products = Map[
      esTimes[esNew[kmin, #[[1, All, c]]], #[[2]]] &, terms];
    Fold[esAdd, First[products], Rest[products]]],
    {c, ncomp}]];

sameResultsQ[a_, b_] := Length[a] === Length[b] &&
  And @@ MapThread[
    #1["EpsWindow"] === #2["EpsWindow"] && TrueQ[esSameQ[#1, #2]] &,
    {a, b}];

(* Different primitive minima and depths force both edges of the output
   window.  The exact-zero third term must still constrain CompleteMax. *)
termsExact = {
  {{{1, 2}, {3, 4}, {5, 6}}, esNew[0, {7, 11, 13, 17}]},
  {{{19, 23}, {29, 31}, {37, 41}}, esNew[-1, {43, 47}]},
  {ConstantArray[0, {3, 2}], esNew[2, {53, 59, 61}]}};
fast = fastContract[termsExact, -2, 0, 2];
reference = referenceContract[termsExact, -2, 0, 2];
assert["mixed_windows_exact", sameResultsQ[fast, reference] &&
  And @@ (# ["EpsWindow"] === <|"Min" -> -3, "CompleteMax" -> -2|> & /@ fast)];

(* Arbitrary-precision complex input follows the production banana path. *)
termsNumeric = MapAt[N[#, 100] &, termsExact, {All, 1, All, All}];
termsNumeric = MapAt[N[#, 100] &, termsNumeric, {All, 2, "Coeffs", All}];
fast = fastContract[termsNumeric, -2, 0, 2];
reference = referenceContract[termsNumeric, -2, 0, 2];
assert["mixed_windows_numeric", sameResultsQ[fast, reference]];

(* Symbolic epsilon-free coefficients exercise the final normalization
   seam; equivalent expanded expressions need not have identical trees. *)
termsSymbolic = {
  {{{Global`q + 1}, {2 Global`q - 3}, {Global`q^2}},
    esNew[-1, {1 - Global`q, 2, 3}]},
  {{{-Global`q}, {4}, {1 + Global`q}}, esNew[0, {5, 7}]}};
fast = fastContract[termsSymbolic, 0, 2, 1];
reference = referenceContract[termsSymbolic, 0, 2, 1];
assert["symbolic_normalization", sameResultsQ[fast, reference]];

(* Deterministic property sweep over varying Laurent minima, primitive
   depths, coefficient widths, zero slabs, and three components. *)
SeedRandom[20260710];
propertyOK = And @@ Table[
  Module[{kmin = RandomInteger[{-4, 1}], width = RandomInteger[{2, 7}],
      ncomp = 3, nterms = RandomInteger[{2, 8}], terms, fast0, ref0},
    terms = Table[Module[{slab, bmin, bwidth, bc},
      slab = RandomInteger[{-7, 7}, {width, ncomp}];
      If[RandomInteger[{1, 5}] === 1, slab = 0 slab];
      bmin = RandomInteger[{-2, 2}];
      bwidth = RandomInteger[{1, 7}];
      bc = RandomInteger[{-5, 5}, bwidth];
      {slab, esNew[bmin, bc]}], {nterms}];
    fast0 = fastContract[terms, kmin, kmin + width - 1, ncomp];
    ref0 = referenceContract[terms, kmin, kmin + width - 1, ncomp];
    sameResultsQ[fast0, ref0]],
  {40}];
assert["randomized_differential_40_cases", propertyOK];

empty = fastContract[{}, -3, 5, 3];
assert["empty_terms_canonical_zero", Length[empty] === 3 &&
  And @@ (# ["EpsWindow"] === <|"Min" -> 5, "CompleteMax" -> 5|> & /@ empty)];

(* Differential pin for the direct primitive coefficient formula.  This is
   the former ESInvert/ESTimes construction, retained only in the test. *)
legacyPrimitive[m_, b_, p_, t_, logT_, phase_, kmax_] := Module[
  {width = kmax + p + 4, tbeps, alpha, inv, jsum},
  tbeps = If[b === 0, esNew[0, PadRight[{1}, width + 1]],
    esNew[0, Table[(b logT)^r/r!, {r, 0, width}]]];
  alpha = If[b === 0,
    esNew[0, PadRight[{m + 1}, width + 1]],
    If[m + 1 === 0, esNew[1, PadRight[{b}, width]],
      esNew[0, PadRight[{m + 1, b}, width + 1]]]];
  inv = esInv[alpha];
  jsum = Module[{acc = None, pwr = inv},
    Do[Module[{term = esScale[(-1)^j logT^(p - j)/(p - j)!, pwr]},
      acc = If[acc === None, term, esAdd[acc, term]];
      pwr = esTimes[pwr, inv]], {j, 0, p}];
    acc];
  esScale[phase*t^(m + 1), esTimes[tbeps, esShift[jsum, p]]]];

primitiveCases = Join[
  Flatten[Table[{-1, b, p, 1/4, 1, 6},
    {b, {-3, 2}}, {p, 0, 4}], 1],
  Flatten[Table[{m, b, p, 2/5, Exp[I Pi (m + 1)], 6},
    {m, {-3, -2, 0, 2}}, {b, {0, 2}}, {p, 0, 3}], 2]];
primitiveOK = And @@ Map[Function[case,
  Module[{m = case[[1]], b = case[[2]], p = case[[3]], t = case[[4]],
      phase = case[[5]], km = case[[6]], logT, direct, legacy},
    logT = N[Log[t], 200];
    direct = DiffExp2`Integrate`Private`antiderivativeAtLog[
      m, b, p, t, logT, phase, km];
    legacy = legacyPrimitive[m, b, p, t, logT, phase, km];
    direct["EpsWindow"] === legacy["EpsWindow"] &&
      TrueQ[esSameQ[direct, legacy]]]], primitiveCases];
assert["direct_primitive_differential_42_cases", primitiveOK];

(* Same-sign m=-1,b!=0 intervals must be paired before epsilon-window
   arithmetic.  Separate endpoint primitives would leave a fake eps^-1 row
   and lose one complete order. *)
mkLS[sec_, kmax_] := <|"Center" -> 0,
  "ChartMap" -> <|"Center" -> 0, "Scale" -> 1|>, "Radius" -> 1,
  "Sectors" -> {sec},
  "EpsWindow" -> <|"Min" -> 0, "CompleteMax" -> kmax|>,
  "TWindow" -> <|"CompleteMax" -> 0|>,
  "ErrorEstimate" -> ConstantArray[0, kmax + 1], "Prescriptions" -> {}|>;
sameSignOK = And @@ Table[Module[{rows, ls, got, l1 = Log[1/5],
    l2 = Log[1/4], expected0},
  rows = Table[{{Boole[k === 0]}}, {k, 0, 4}];
  ls = mkLS[<|"a" -> -1, "b" -> 2, "p" -> p,
    "Coeffs" -> rows|>, 4];
  got = DiffExp2`Integrate`IntegrateLocalSolution[ls, {1/5, 1/4}]["Values"][[1]];
  expected0 = (l2^(p + 1) - l1^(p + 1))/(p!*(p + 1));
  got["EpsWindow"] === <|"Min" -> p, "CompleteMax" -> 4 + p|> &&
    TrueQ[PossibleZeroQ[DiffExp2`EpsSeries`ESCoefficient[got, p] - expected0]]],
  {p, 0, 3}];
assert["same_positive_pole_combined_windows_p0_to_p3", sameSignOK];

(* The negative arm first applies the configured +i0 crossing.  Besides
   exercising the same paired primitive, this checks that its honest window
   survives the phase/log-chain convolution used by analytic regularization. *)
sameNegativeOK = And @@ Table[Module[{rows, ls, got,
    l1 = N[Log[1/5], 200], l2 = N[Log[1/4], 200], expected0},
  rows = Table[{{Boole[k === 0]}}, {k, 0, 4}];
  ls = mkLS[<|"a" -> -1, "b" -> 2, "p" -> p,
    "Coeffs" -> rows|>, 4];
  got = DiffExp2`Integrate`IntegrateLocalSolution[ls, {-1/4, -1/5}]["Values"][[1]];
  expected0 = -((l2 + I Pi)^(p + 1) - (l1 + I Pi)^(p + 1))/
    (p!*(p + 1));
  (* ApplyCrossing retains the source's common [0,4] window while moving
     the sector's first nonzero coefficient to order p. *)
  got["EpsWindow"] === <|"Min" -> 0, "CompleteMax" -> 4|> &&
    And @@ Table[TrueQ[PossibleZeroQ[
      DiffExp2`EpsSeries`ESCoefficient[got, k]]], {k, 0, p - 1}] &&
    Abs[N[DiffExp2`EpsSeries`ESCoefficient[got, p] - expected0, 80]] < 10^-70],
  {p, 0, 3}];
assert["same_negative_pole_crossed_windows_p0_to_p3", sameNegativeOK];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Print["Some tests FAILED."]; Quit[1],
  Print["All tests PASSED!"]; Quit[0]];
