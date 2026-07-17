(* Independent Feynman-parameter oracle for the strictly Euclidean massive
   one-loop pentagon fixture.

   Run with:
     wolframscript -file Scripts/pentagon_massive_oracle.m

   PENTAGON_MASSIVE_ORACLE_DIGITS changes the comparison target (default 16,
   maximum 18 because the independent pins carry about 20 digits).
   PENTAGON_MASSIVE_ORACLE_MAX_ORDER selects epsilon orders 0 through n
   (default 2).  DiffExp, FeynmanTrick, and FIRE are not loaded. *)

pmEnvOrDefault[name_, default_] := Module[{value = Environment[name]},
  If[StringQ[value] && StringLength[StringTrim[value]] > 0, value, default]];

(* With Sum[x_i]=1, the Euclidean Symanzik polynomial is bounded below by
   one on this fixture.  Thus the parameter integral and its epsilon
   derivatives are ordinary, absolutely convergent real integrals. *)
pmPentagonSymanzik[x_List] /; Length[x] === 5 :=
  x[[1]] + 3 x[[2]]/2 + 4 x[[3]]/3 + 5 x[[4]]/4 + 6 x[[5]]/5 +
  x[[1]] x[[2]] + x[[2]] x[[3]] + x[[3]] x[[4]] +
  x[[4]] x[[5]] + x[[5]] x[[1]] +
  3 (x[[1]] x[[3]] + x[[1]] x[[4]] + x[[2]] x[[4]] +
      x[[2]] x[[5]] + x[[3]] x[[5]])/2;

(* Stick breaking maps the unit four-cube bijectively to the simplex, up to
   measure-zero faces, without introducing endpoint singularities. *)
pmSimplexPoint[u1_, u2_, u3_, u4_] := {
  u1,
  (1 - u1) u2,
  (1 - u1) (1 - u2) u3,
  (1 - u1) (1 - u2) (1 - u3) u4,
  (1 - u1) (1 - u2) (1 - u3) (1 - u4)};
pmStickJacobian[u1_, u2_, u3_] := (1 - u1)^3 (1 - u2)^2 (1 - u3);

pmGammaCoefficient[k_Integer?NonNegative] :=
  SeriesCoefficient[Gamma[3 + Global`pmEpsilon],
    {Global`pmEpsilon, 0, k}];

(* Coefficient of eps^n in Gamma[3+eps] F^(-3-eps). *)
pmCoefficientDensity[n_Integer?NonNegative, x_List] /; Length[x] === 5 :=
  With[{f = pmPentagonSymanzik[x]},
    Sum[pmGammaCoefficient[k] (-Log[f])^(n - k)/(n - k)!,
      {k, 0, n}]/f^3];

pmCubeCoefficientExpression[n_Integer?NonNegative] :=
  pmStickJacobian[Global`pmu1, Global`pmu2, Global`pmu3] *
    pmCoefficientDensity[n,
      pmSimplexPoint[Global`pmu1, Global`pmu2, Global`pmu3, Global`pmu4]];

pmIntegrateCoefficient[n_Integer?NonNegative, wp_Integer?Positive,
    goal_Integer?Positive] := NIntegrate[
  Evaluate[pmCubeCoefficientExpression[n]],
  {Global`pmu1, 0, 1}, {Global`pmu2, 0, 1},
  {Global`pmu3, 0, 1}, {Global`pmu4, 0, 1},
  WorkingPrecision -> wp, AccuracyGoal -> goal, PrecisionGoal -> goal,
  MaxRecursion -> 30,
  Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0}];

(* Independent high-precision pins obtained from the same convergent
   parameter representation, outside the DiffExp/FeynmanTrick ladder. *)
pmOracleReference[] := <|
  0 -> 1813378668630195764/10^20,
  1 -> 7613115416144053565/10^21,
  2 -> 5214475578477681142/10^21|>;

pmNormalizeCoefficientData[data_Association] := data;
pmNormalizeCoefficientData[data_List] /;
    AllTrue[data, MatchQ[#, {_Integer, _?NumericQ}] &] :=
  Association[Rule @@@ data];
pmNormalizeCoefficientData[data_List] /;
    AllTrue[data, MatchQ[#, Rule[_Integer, _?NumericQ]] &] :=
  Association[data];
pmNormalizeCoefficientData[_] := $Failed;

pmOracleCompare[data_, digits_Integer?Positive,
    requiredPowers_: Automatic] := Module[
  {actual, reference = pmOracleReference[], requested, powers, rows},
  actual = pmNormalizeCoefficientData[data];
  If[actual === $Failed, Return[$Failed]];
  requested = If[requiredPowers === Automatic, Keys[reference], requiredPowers];
  If[!ListQ[requested] || Complement[requested, Keys[reference]] =!= {},
    Return[$Failed]];
  powers = Intersection[requested, Keys[actual]];
  rows = Association@Table[
    power -> With[
      {expected = N[reference[power], digits + 10], value = actual[power]},
      <|"Actual" -> value, "Reference" -> expected,
        "AbsDifference" -> Abs[value - expected],
        "Pass" -> TrueQ[Abs[value - expected] < 10^-digits]|>],
    {power, powers}];
  <|"Pass" -> (powers === requested &&
      AllTrue[Values[rows], TrueQ[#["Pass"]] &]),
    "Digits" -> digits, "Coefficients" -> rows|>];

pmNumberString[value_, digits_Integer?Positive] := ToString[
  NumberForm[value, digits, NumberPadding -> {"", "0"},
    NumberFormat -> (Row[{#1, If[#3 === "", "", "e" <> #3]}] &)],
  OutputForm];

runPentagonMassiveOracle[] := Module[
  {digits, maxOrder, wp, goal, elapsed, values, comparison, records, status},
  digits = Quiet[Check[ToExpression[
    pmEnvOrDefault["PENTAGON_MASSIVE_ORACLE_DIGITS", "16"], InputForm],
    $Failed]];
  maxOrder = Quiet[Check[ToExpression[
    pmEnvOrDefault["PENTAGON_MASSIVE_ORACLE_MAX_ORDER", "2"], InputForm],
    $Failed]];
  If[!IntegerQ[digits] || !Between[digits, {8, 18}],
    Print["PENTAGON_MASSIVE_ORACLE_DIGITS must be an integer from 8 to 18"];
    Exit[2]];
  If[!IntegerQ[maxOrder] || !Between[maxOrder, {0, 2}],
    Print["PENTAGON_MASSIVE_ORACLE_MAX_ORDER must be 0, 1, or 2"];
    Exit[2]];

  wp = Max[60, digits + 30];
  goal = digits + 6;
  elapsed = First@AbsoluteTiming[
    values = Association@Table[
      order -> pmIntegrateCoefficient[order, wp, goal],
      {order, 0, maxOrder}]];
  comparison = pmOracleCompare[values, digits, Range[0, maxOrder]];
  records = Table[With[
    {row = comparison["Coefficients"][order]},
    <|"EpsilonPower" -> order,
      "Value" -> pmNumberString[row["Actual"], digits],
      "Reference" -> pmNumberString[row["Reference"], digits],
      "AbsDifference" -> pmNumberString[row["AbsDifference"], digits],
      "Pass" -> row["Pass"]|>],
    {order, 0, maxOrder}];
  status = If[TrueQ[comparison["Pass"]], "PASS", "FAIL"];
  Print["PENTAGON_MASSIVE_ORACLE ", ExportString[
    <|"Example" -> "pentagon_massive", "Master" -> {1, 1, 1, 1, 1},
      "Dimension" -> "4-2 eps", "TargetDigits" -> digits,
      "WorkingPrecision" -> wp, "AccuracyGoal" -> goal,
      "Coefficients" -> records, "ElapsedSeconds" -> N[elapsed, 6],
      "Status" -> status|>, "RawJSON", "Compact" -> True]];
  If[status === "PASS", 0, 1]
];

$pmPentagonMassiveOracleRan = False;
If[pmEnvOrDefault["PENTAGON_MASSIVE_ORACLE_DEFINITIONS_ONLY", "0"] =!= "1",
  $pmPentagonMassiveOracleRan = True;
  Exit[runPentagonMassiveOracle[]]];
