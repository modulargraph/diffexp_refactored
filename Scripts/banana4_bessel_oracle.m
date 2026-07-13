(* Independent D=2, eps^0 oracle for the four-loop banana.

   Run with:
     wolframscript -file Scripts/banana4_bessel_oracle.m

   BANANA4_BESSEL_DIGITS changes the requested significant digits (default 50).
   BANANA4_BESSEL_EXAMPLE selects banana4 (default) or banana4_unequal.
   The expensive FeynmanTrick/DiffExp ladder is not loaded or run. *)

b4EnvOrDefault[name_, default_] := Module[{value = Environment[name]},
  If[StringQ[value] && StringLength[StringTrim[value]] > 0, value, default]];

b4BesselPrefactor[n_Integer?Positive] := 2^(n - 1);
b4BesselDensityN[n_Integer?Positive, r_] :=
  b4BesselPrefactor[n] r BesselJ[0, r] BesselK[0, r]^n;
b4BesselDensityMasses[masses_List, r_] :=
  b4BesselPrefactor[Length[masses]] r BesselJ[0, r] Times @@
    (BesselK[0, # r] & /@ masses);
b4BesselDensity[r_] := b4BesselDensityN[5, r];
b4SquaredDensityN[n_Integer?Positive, t_?NumericQ] := If[t == 0, 0,
  2 t b4BesselDensityN[n, t^2]];
b4SquaredDensityMasses[masses_List, t_?NumericQ] := If[t == 0, 0,
  2 t b4BesselDensityMasses[masses, t^2]];
b4SquaredDensity[t_?NumericQ] := b4SquaredDensityN[5, t];
b4LogDensity[u_?NumericQ] := Exp[u] b4BesselDensity[Exp[u]];
b4LogDensityMasses[masses_List, u_?NumericQ] :=
  Exp[u] b4BesselDensityMasses[masses, Exp[u]];

b4SquaredMassesForExample["banana4"] := ConstantArray[1, 5];
b4SquaredMassesForExample["banana4_unequal"] := {2, 3/2, 4/3, 5/4, 1};
b4SquaredMassesForExample[_] := $Failed;

(* K0(r) < Sqrt[Pi/(2 r)] Exp[-r] gives this rigorous r>R bound. *)
b4UpperTailBoundMasses[R_, masses_List] := Module[
  {n = Length[masses], massSum = Total[masses],
   massProduct = Times @@ masses},
  b4BesselPrefactor[Length[masses]] (Pi/2)^(Length[masses]/2) *
    massSum^(n/2 - 2)/Sqrt[massProduct] *
    Gamma[2 - n/2, massSum R]
];
b4UpperTailBound[R_] :=
  b4UpperTailBoundMasses[R, ConstantArray[1, 5]];

(* For 0<r<1, K0(r)<Log[2/r].  The same bound applies to every
   K0(m_i r) only when m_i>=1, which is true for both shipped fixtures. *)
b4LowerTailBound[rmin_] := Gamma[6, 2 Log[2/rmin]];
b4LowerTailBoundMasses[rmin_, masses_List] :=
  If[AllTrue[masses, TrueQ[# >= 1] &],
    b4LowerTailBound[rmin], $Failed];

b4OracleAcceptQ[direct_?NumericQ, logValue_?NumericQ,
    upper_?NumericQ, lower_?NumericQ, digits_Integer?Positive] :=
  Abs[direct - logValue] + Abs[upper] + Abs[lower] < 10^-digits;

b4NumberString[value_, digits_Integer] := ToString[
  NumberForm[value, digits, NumberPadding -> {"", "0"},
    NumberFormat -> (Row[{#1, If[#3 === "", "", "e" <> #3]}] &)],
  OutputForm];

runBanana4BesselOracle[] := Module[
  {example, squaredMasses, masses, massSum, digits, wp, goal, cutoff,
    umin, rmin, rCuts, tCuts, uCuts,
    direct, logValue, normalization, normalizationReference,
    normalizationDigits, normalizationDelta, upper, lower, delta, elapsed,
    status, record},
  example = b4EnvOrDefault["BANANA4_BESSEL_EXAMPLE", "banana4"];
  squaredMasses = b4SquaredMassesForExample[example];
  If[squaredMasses === $Failed,
    Print["BANANA4_BESSEL_EXAMPLE must be banana4 or banana4_unequal"];
    Exit[2]];
  masses = Sqrt[squaredMasses];
  massSum = Total[masses];
  If[b4LowerTailBoundMasses[Exp[-1], masses] === $Failed,
    Print["BESSEL_ORACLE lower-tail certificate requires masses >= 1"];
    Exit[2]];

  digits = Quiet[Check[
    ToExpression[b4EnvOrDefault["BANANA4_BESSEL_DIGITS", "50"], InputForm],
    $Failed]];
  If[!IntegerQ[digits] || digits < 10,
    Print["BANANA4_BESSEL_DIGITS must be an integer >= 10"];
    Exit[2]];

  wp = Max[100, digits + 50];
  goal = digits + 5;
  cutoff = Max[32,
    Ceiling[N[(digits + 15) Log[10]/massSum, 30]]];
  umin = -Max[100, 2 digits];
  rmin = Exp[umin];
  rCuts = Join[{0, 1/1000, 1/100, 1/10, 1/2, 1, 2, 4, 8, 16},
    {cutoff}];
  tCuts = Sqrt[rCuts];
  uCuts = Join[{umin},
    Select[{-100, -50, -25, -12, -6, -3, -1, 0, 1, 2},
      TrueQ[N[umin, 30] < # < N[Log[cutoff], 30]] &],
    {Log[cutoff]}];

  elapsed = First@AbsoluteTiming[
    direct = Check[NIntegrate[b4SquaredDensityMasses[masses, t],
      Evaluate[Prepend[tCuts, t]], WorkingPrecision -> wp,
      AccuracyGoal -> goal, PrecisionGoal -> goal, MaxRecursion -> 40,
      Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0}], $Failed];
    logValue = Check[NIntegrate[b4LogDensityMasses[masses, u],
      Evaluate[Prepend[uCuts, u]], WorkingPrecision -> wp,
      AccuracyGoal -> goal, PrecisionGoal -> goal, MaxRecursion -> 40,
      Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0}], $Failed];
    normalization = Check[NIntegrate[b4SquaredDensityN[4, t],
      Evaluate[Prepend[tCuts, t]], WorkingPrecision -> wp,
      AccuracyGoal -> goal, PrecisionGoal -> goal, MaxRecursion -> 40,
      Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0}], $Failed];
  ];
  If[direct === $Failed || logValue === $Failed || normalization === $Failed,
    Print["BESSEL_ORACLE integration failed"];
    Exit[1]];

  normalizationReference = N[
    8268104535868968731543015345479988868728618483845/10^48, wp];
  normalizationDigits = Min[40, Max[8, digits - 5]];
  normalizationDelta = Abs[normalization - normalizationReference];
  upper = N[b4UpperTailBoundMasses[cutoff, masses], wp];
  lower = N[b4LowerTailBoundMasses[rmin, masses], wp];
  delta = Abs[direct - logValue];
  status = If[b4OracleAcceptQ[direct, logValue, upper, lower, digits] &&
      normalizationDelta < 10^-normalizationDigits,
    "PASS", "FAIL"];
  record = <|
    "Example" -> example, "Level" -> 0,
    "Master" -> {1, 1, 1, 1, 1}, "EpsilonPower" -> 0,
    "Dimension" -> 2, "PSquared" -> -1,
    "SquaredMasses" -> (ToString[#, InputForm] & /@ squaredMasses),
    "Formula" -> "16 Integral_0^Infinity r J0(r) Product_i K0(m_i r) dr",
    "TargetDigits" -> digits, "WorkingPrecision" -> wp,
    "AccuracyGoal" -> goal, "PrecisionGoal" -> goal,
    "Certification" -> "cross-chart agreement with analytic tail bounds",
    "Cutoff" -> cutoff, "Value" -> b4NumberString[direct, digits],
    "LogChartValue" -> b4NumberString[logValue, digits],
    "CrosscheckAbs" -> b4NumberString[delta, digits],
    "FourLineValue" -> b4NumberString[normalization, digits],
    "FourLineReference" -> b4NumberString[normalizationReference, digits],
    "FourLineAbsDifference" -> b4NumberString[normalizationDelta, digits],
    "NormalizationDigits" -> normalizationDigits,
    "UpperTailBound" -> b4NumberString[upper, digits],
    "LowerTailBound" -> b4NumberString[lower, digits],
    "ElapsedSeconds" -> N[elapsed, 6], "Status" -> status|>;
  Print["BESSEL_ORACLE ",
    ExportString[record, "RawJSON", "Compact" -> True]];
  If[status === "PASS", 0, 1]
];

$b4BesselOracleRan = False;
If[b4EnvOrDefault["BANANA4_BESSEL_DEFINITIONS_ONLY", "0"] =!= "1",
  $b4BesselOracleRan = True;
  Exit[runBanana4BesselOracle[]]];
