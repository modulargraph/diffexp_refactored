(* Independent D=2, eps^0 oracle for the equal-mass four-loop banana.

   Run with:
     wolframscript -file Scripts/banana4_bessel_oracle.m

   BANANA4_BESSEL_DIGITS changes the requested significant digits (default 50).
   The expensive FeynmanTrick/DiffExp ladder is not loaded or run. *)

b4EnvOrDefault[name_, default_] := Module[{value = Environment[name]},
  If[StringQ[value] && StringLength[StringTrim[value]] > 0, value, default]];

b4BesselPrefactor[n_Integer?Positive] := 2^(n - 1);
b4BesselDensity[r_] := 16 r BesselJ[0, r] BesselK[0, r]^5;
b4SquaredDensity[0] = 0;
b4SquaredDensity[t_?NumericQ] := 2 t b4BesselDensity[t^2];
b4LogDensity[u_?NumericQ] := Exp[u] b4BesselDensity[Exp[u]];

(* K0(r) < Sqrt[Pi/(2 r)] Exp[-r] gives this rigorous r>R bound. *)
b4UpperTailBound[R_] :=
  16 (Pi/2)^(5/2) Sqrt[5] Gamma[-1/2, 5 R];

(* For 0<r<1, K0(r)<Log[2/r]; this bounds the log-chart cutoff. *)
b4LowerTailBound[rmin_] := Gamma[6, 2 Log[2/rmin]];

b4OracleAcceptQ[direct_?NumericQ, logValue_?NumericQ,
    upper_?NumericQ, lower_?NumericQ, digits_Integer?Positive] :=
  Abs[direct - logValue] + Abs[upper] + Abs[lower] < 10^-digits;

b4NumberString[value_, digits_Integer] := ToString[
  NumberForm[value, digits, NumberPadding -> {"", "0"},
    NumberFormat -> (Row[{#1, If[#3 === "", "", "e" <> #3]}] &)],
  OutputForm];

runBanana4BesselOracle[] := Module[
  {digits, wp, goal, cutoff, umin, rmin, rCuts, tCuts, uCuts,
    direct, logValue, upper, lower, delta, elapsed, status, record},
  digits = Quiet[Check[
    ToExpression[b4EnvOrDefault["BANANA4_BESSEL_DIGITS", "50"], InputForm],
    $Failed]];
  If[!IntegerQ[digits] || digits < 10,
    Print["BANANA4_BESSEL_DIGITS must be an integer >= 10"];
    Exit[2]];

  wp = Max[100, digits + 50];
  goal = digits + 5;
  cutoff = Max[32, Ceiling[N[(digits + 15) Log[10]/5, 30]]];
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
    direct = Check[NIntegrate[b4SquaredDensity[t],
      Evaluate[Prepend[tCuts, t]], WorkingPrecision -> wp,
      AccuracyGoal -> goal, PrecisionGoal -> goal, MaxRecursion -> 40,
      Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0}], $Failed];
    logValue = Check[NIntegrate[b4LogDensity[u],
      Evaluate[Prepend[uCuts, u]], WorkingPrecision -> wp,
      AccuracyGoal -> goal, PrecisionGoal -> goal, MaxRecursion -> 40,
      Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0}], $Failed];
  ];
  If[direct === $Failed || logValue === $Failed,
    Print["BESSEL_ORACLE integration failed"];
    Exit[1]];

  upper = N[b4UpperTailBound[cutoff], wp];
  lower = N[b4LowerTailBound[rmin], wp];
  delta = Abs[direct - logValue];
  status = If[b4OracleAcceptQ[direct, logValue, upper, lower, digits],
    "PASS", "FAIL"];
  record = <|
    "Example" -> "banana4", "Level" -> 0,
    "Master" -> {1, 1, 1, 1, 1}, "EpsilonPower" -> 0,
    "Dimension" -> 2, "PSquared" -> -1,
    "Formula" -> "16 Integral_0^Infinity r J0(r) K0(r)^5 dr",
    "TargetDigits" -> digits, "WorkingPrecision" -> wp,
    "AccuracyGoal" -> goal, "PrecisionGoal" -> goal,
    "Cutoff" -> cutoff, "Value" -> b4NumberString[direct, digits],
    "LogChartValue" -> b4NumberString[logValue, digits],
    "CrosscheckAbs" -> b4NumberString[delta, digits],
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
