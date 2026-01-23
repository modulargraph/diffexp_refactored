(* Benchmark: Rational recurrence speedup *)
(* Test A: Equal mass banana WITH Mobius (4x4 system, fast) *)
(* Test B: Unequal mass banana WITHOUT Mobius (15x15, check it works) *)

SetDirectory[DirectoryName[$InputFileName]];
scriptDir = Directory[];
parentDir = ParentDirectory[scriptDir];
$Path = Prepend[$Path, FileNameJoin[{parentDir, "DiffExp"}]];

Print["==========================================="];
Print["Rational Recurrence Speedup Benchmark"];
Print["===========================================\n"];

Get["DiffExp.m"];

(* ====== TEST A: Equal mass banana WITH Mobius (4x4) ====== *)
Print["=== TEST A: Equal mass (4x4) WITH Mobius ===\n"];

EqualMassBoundaryConditions = {
  "?",
  "?",
  \[Epsilon] (1 + 3 \[Epsilon]) (1 + 4 \[Epsilon]) * (
    -4 E^(3 EulerGamma \[Epsilon]) Gamma[\[Epsilon]]^3/t +
    6 E^(3 EulerGamma \[Epsilon]) (-1/t)^(1 + \[Epsilon]) \[Epsilon] Gamma[-\[Epsilon]]^2 Gamma[\[Epsilon]]^3/Gamma[-2 \[Epsilon]] +
    8 E^(3 EulerGamma \[Epsilon]) (-1/t)^(1 + 2 \[Epsilon]) \[Epsilon] Gamma[-\[Epsilon]]^3 Gamma[\[Epsilon]] Gamma[2 \[Epsilon]]/Gamma[-3 \[Epsilon]] +
    3 E^(3 EulerGamma \[Epsilon]) (-1/t)^(1 + 3 \[Epsilon]) \[Epsilon] Gamma[-\[Epsilon]]^4 Gamma[3 \[Epsilon]]/Gamma[-4 \[Epsilon]]
  ),
  E^(3 EulerGamma \[Epsilon]) \[Epsilon]^3 Gamma[\[Epsilon]]^3
};

(* Run A1: Default *)
Print["A1: Default strategy..."];
LoadConfiguration[{
  DeltaPrescriptions -> {t - 16 + I * \[Delta]},
  MatrixDirectory -> FileNameJoin[{scriptDir, "Banana_EqualMass_Matrices"}] <> "/",
  Verbosity -> 1, UseMobius -> True, UsePade -> True,
  DivisionOrder -> 4, ExpansionOrder -> 70,
  "UseRationalRecurrence" -> False
}];
PreparedBCs = PrepareBoundaryConditions[EqualMassBoundaryConditions, <|t -> -1/x|>];
timeA1 = AbsoluteTime[];
resultA1 = TransportTo[PreparedBCs, <|t -> -1/x|>];
resultA1 = TransportTo[resultA1, <|t -> 1/2|>];
timeA1 = AbsoluteTime[] - timeA1;
Print["  Time: ", timeA1, " s\n"];

(* Run A2: Rational Recurrence *)
Print["A2: Rational recurrence..."];
LoadConfiguration[{
  DeltaPrescriptions -> {t - 16 + I * \[Delta]},
  MatrixDirectory -> FileNameJoin[{scriptDir, "Banana_EqualMass_Matrices"}] <> "/",
  Verbosity -> 1, UseMobius -> True, UsePade -> True,
  DivisionOrder -> 4, ExpansionOrder -> 70,
  "UseRationalRecurrence" -> True
}];
PreparedBCs2 = PrepareBoundaryConditions[EqualMassBoundaryConditions, <|t -> -1/x|>];
timeA2 = AbsoluteTime[];
resultA2 = TransportTo[PreparedBCs2, <|t -> -1/x|>];
resultA2 = TransportTo[resultA2, <|t -> 1/2|>];
timeA2 = AbsoluteTime[] - timeA2;
Print["  Time: ", timeA2, " s\n"];

diffA = Max[Abs[Flatten[N[resultA1["SeriesValues"] - resultA2["SeriesValues"], 100]]]];
Print["  TEST A Results (equal mass, Mobius):"];
Print["    Default:    ", timeA1, " s"];
Print["    Recurrence: ", timeA2, " s"];
Print["    Speedup:    ", N[timeA1/timeA2, 3], "x"];
Print["    Max diff:   ", ScientificForm[diffA, 3]];
Print[""];

(* ====== TEST B: Unequal mass banana WITHOUT Mobius (15x15) ====== *)
Print["=== TEST B: Unequal mass (15x15) WITHOUT Mobius ===\n"];

(* Use already-computed BCs from Test A *)
UnequalMassBoundaryConditions = {
  <|psq -> 1/2, mm1 -> 1, mm2 -> 1, mm3 -> 1, mm4 -> 1|>,
  {
    resultA1["SeriesValues"][[1]], resultA1["SeriesValues"][[1]], resultA1["SeriesValues"][[1]],
    resultA1["SeriesValues"][[1]], resultA1["SeriesValues"][[1]], resultA1["SeriesValues"][[1]],
    resultA1["SeriesValues"][[2]], resultA1["SeriesValues"][[2]], resultA1["SeriesValues"][[2]], resultA1["SeriesValues"][[2]],
    resultA1["SeriesValues"][[3]],
    resultA1["SeriesValues"][[4]], resultA1["SeriesValues"][[4]], resultA1["SeriesValues"][[4]], resultA1["SeriesValues"][[4]]
  }
};
unequalMassLine = <|psq -> 1/2, mm1 -> 1 + x, mm2 -> 1 + x/2, mm3 -> 1 + x/3, mm4 -> 1|>;

(* Run B1: Default, no Mobius *)
Print["B1: Default strategy (no Mobius)..."];
LoadConfiguration[{
  ChopPrecision -> 500, DivisionOrder -> 4, ExpansionOrder -> 70,
  MatrixDirectory -> FileNameJoin[{scriptDir, "Banana_Matrices"}] <> "/",
  RadiusOfConvergence -> 10, UseMobius -> False, UsePade -> True,
  WorkingPrecision -> 1000, Verbosity -> 1,
  "UseRationalRecurrence" -> False
}];

timeB1 = AbsoluteTime[];
resultB1 = TransportTo[UnequalMassBoundaryConditions, unequalMassLine, 1];
timeB1 = AbsoluteTime[] - timeB1;
Print["  Time: ", timeB1, " s\n"];

(* Run B2: Rational recurrence, no Mobius *)
Print["B2: Rational recurrence (no Mobius)..."];
LoadConfiguration[{
  ChopPrecision -> 500, DivisionOrder -> 4, ExpansionOrder -> 70,
  MatrixDirectory -> FileNameJoin[{scriptDir, "Banana_Matrices"}] <> "/",
  RadiusOfConvergence -> 10, UseMobius -> False, UsePade -> True,
  WorkingPrecision -> 1000, Verbosity -> 1,
  "UseRationalRecurrence" -> True
}];

timeB2 = AbsoluteTime[];
resultB2 = TransportTo[UnequalMassBoundaryConditions, unequalMassLine, 1];
timeB2 = AbsoluteTime[] - timeB2;
Print["  Time: ", timeB2, " s\n"];

diffB = Max[Abs[Flatten[N[resultB1["SeriesValues"] - resultB2["SeriesValues"], 200]]]];
Print["  TEST B Results (unequal mass, no Mobius):"];
Print["    Default:    ", timeB1, " s"];
Print["    Recurrence: ", timeB2, " s"];
Print["    Speedup:    ", N[timeB1/timeB2, 3], "x"];
Print["    Max diff:   ", ScientificForm[diffB, 3]];

(* ====== Summary ====== *)
Print["\n==========================================="];
Print["SUMMARY"];
Print["==========================================="];
Print["Test A (4x4, Mobius):   ", N[timeA1/timeA2, 3], "x speedup, diff = ", ScientificForm[diffA, 2]];
Print["Test B (15x15, no Mob): ", N[timeB1/timeB2, 3], "x speedup, diff = ", ScientificForm[diffB, 2]];

testAPass = diffA < 10^-30;
testBPass = diffB < 10^-30;
If[testAPass, Print["  [PASS] Test A"], Print["  [FAIL] Test A: ", ScientificForm[diffA, 3]]];
If[testBPass, Print["  [PASS] Test B"], Print["  [FAIL] Test B: ", ScientificForm[diffB, 3]]];
Print["==========================================="];
