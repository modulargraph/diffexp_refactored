(* Structural checks for the opt-in banana4 Bessel oracle.  NIntegrate is
   never called: the script is loaded in definitions-only mode. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
SetEnvironment["BANANA4_BESSEL_DEFINITIONS_ONLY" -> "1"];
Get[FileNameJoin[{repoRoot, "Scripts", "banana4_bessel_oracle.m"}]];

passed = 0; failed = 0;
assert[label_, cond_] := If[TrueQ[cond],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];

r = Global`r;
tail32 = N[b4UpperTailBound[32], 40];
lower100 = N[b4LowerTailBound[Exp[-100]], 40];
reference = N[
  396555268342976525299928230469335811564460602187101868181/
    10^55, 70];
unequalReference = N[
  28647167500233570163016705540503359169504240713994/10^48, 70];
unequalSquaredMasses = {2, 3/2, 4/3, 5/4, 1};
unequalMasses = Sqrt[unequalSquaredMasses];
unequalMassSum = Total[unequalMasses];
unequalMassProduct = Times @@ unequalMasses;

assert["banana4_oracle_definitions_only",
  $b4BesselOracleRan === False];
assert["banana4_oracle_exact_bessel_density",
  b4BesselPrefactor[5] === 16 && b4BesselPrefactor[4] === 8 &&
  b4BesselDensity[r] === 16 r BesselJ[0, r] BesselK[0, r]^5];
assert["banana4_unequal_oracle_exact_bessel_density",
  b4SquaredMassesForExample["banana4_unequal"] === unequalSquaredMasses &&
  b4BesselDensityMasses[unequalMasses, r] ===
    16 r BesselJ[0, r] Times @@
      (BesselK[0, # r] & /@ unequalMasses)];
assert["banana4_unequal_oracle_chart_jacobians",
  b4SquaredDensityMasses[unequalMasses, 2] ===
    4 b4BesselDensityMasses[unequalMasses, 4] &&
  b4LogDensityMasses[unequalMasses, 0] ===
    b4BesselDensityMasses[unequalMasses, 1]];
assert["banana4_oracle_tail_bound_at_32",
  10^-72 < tail32 < 2 10^-71 &&
  10^-77 < lower100 < 2 10^-76];
assert["banana4_unequal_oracle_exact_tail_certificates",
  b4UpperTailBoundMasses[32, unequalMasses] ===
    16 (Pi/2)^(5/2) Sqrt[unequalMassSum]/Sqrt[unequalMassProduct] *
      Gamma[-1/2, 32 unequalMassSum] &&
  0 < N[b4UpperTailBoundMasses[32, unequalMasses], 40] < tail32 &&
  b4LowerTailBoundMasses[Exp[-100], unequalMasses] ===
    b4LowerTailBound[Exp[-100]] &&
  b4LowerTailBoundMasses[Exp[-100], {1/2, 1, 1, 1, 1}] === $Failed];
assert["banana4_oracle_value_format_pin",
  b4NumberString[reference, 50] ===
      "39.655526834297652529992823046933581156446060218710" &&
  b4NumberString[unequalReference, 50] ===
      "28.647167500233570163016705540503359169504240713994"];
assert["banana4_oracle_acceptance_is_strict",
  b4OracleAcceptQ[reference, reference + 10^-60, tail32, 10^-75, 50] &&
  !b4OracleAcceptQ[reference, reference + 10^-20, tail32, 10^-75, 50] &&
  !b4OracleAcceptQ[reference, reference, 10^-40, 10^-75, 50] &&
  !b4OracleAcceptQ[reference, reference, tail32, 10^-40, 50]];

SetEnvironment["BANANA4_BESSEL_DEFINITIONS_ONLY" -> None];
Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Quit[1], Quit[0]];
