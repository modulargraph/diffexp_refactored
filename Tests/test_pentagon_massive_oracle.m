(* Fast structural checks for the independent massive-pentagon parameter
   oracle.  The four-dimensional NIntegrate is never called here. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
SetEnvironment["PENTAGON_MASSIVE_ORACLE_DEFINITIONS_ONLY" -> "1"];
Get[FileNameJoin[{repoRoot, "Scripts", "pentagon_massive_oracle.m"}]];

passed = 0; failed = 0;
assert[label_, cond_] := If[TrueQ[cond],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];

x = Array[Global`pmx, 5];
u = Array[Global`pmu, 4];
f = pmPentagonSymanzik[x];
simplexPoint = pmSimplexPoint @@ u;
reference = N[Values[pmOracleReference[]], 30];
referencePairs = MapThread[List, {Range[0, 2], reference}];

assert["pentagon_massive_oracle_definitions_only",
  $pmPentagonMassiveOracleRan === False];
assert["pentagon_massive_oracle_symanzik_polynomial",
  f === x[[1]] + 3 x[[2]]/2 + 4 x[[3]]/3 + 5 x[[4]]/4 +
    6 x[[5]]/5 + x[[1]] x[[2]] + x[[2]] x[[3]] +
    x[[3]] x[[4]] + x[[4]] x[[5]] + x[[5]] x[[1]] +
    3 (x[[1]] x[[3]] + x[[1]] x[[4]] + x[[2]] x[[4]] +
      x[[2]] x[[5]] + x[[3]] x[[5]])/2];
assert["pentagon_massive_oracle_cube_maps_to_simplex",
  Together[Total[simplexPoint]] === 1 &&
  pmStickJacobian @@ Most[u] === (1 - u[[1]])^3 (1 - u[[2]])^2
    (1 - u[[3]])];
assert["pentagon_massive_oracle_epsilon_density",
  pmCoefficientDensity[0, x] === 2/f^3 &&
  Together[pmCoefficientDensity[1, x] -
    (3 - 2 EulerGamma - 2 Log[f])/f^3] === 0];
assert["pentagon_massive_oracle_cube_density",
  Together[(pmCubeCoefficientExpression[0] /.
      Thread[{Global`pmu1, Global`pmu2, Global`pmu3, Global`pmu4} -> u]) -
    pmStickJacobian @@ Most[u] *
      pmCoefficientDensity[0, simplexPoint]] === 0];
assert["pentagon_massive_oracle_reference_pin",
  pmNumberString[reference[[1]], 20] === "0.018133786686301957640" &&
  pmNumberString[reference[[2]], 20] === "0.0076131154161440535650" &&
  pmNumberString[reference[[3]], 20] === "0.0052144755784776811420"];
assert["pentagon_massive_oracle_comparison_contract",
  TrueQ[pmOracleCompare[referencePairs, 18]["Pass"]] &&
  !TrueQ[pmOracleCompare[MapThread[List, {Range[0, 2],
      ReplacePart[reference, 2 -> reference[[2]] + 10^-8]}], 12]["Pass"]] &&
  pmOracleCompare[{reference[[1]], reference[[2]]}, 12] === $Failed];

SetEnvironment["PENTAGON_MASSIVE_ORACLE_DEFINITIONS_ONLY" -> None];
Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Quit[1], Quit[0]];
