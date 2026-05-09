(* Unequal-mass banana parity and speed regression *)
(* Compares the standard solver with the rational recurrence path at two nontrivial endpoint transports. *)

SetDirectory[DirectoryName[$InputFileName]];
scriptDir = Directory[];
parentDir = ParentDirectory[scriptDir];
$Path = Prepend[$Path, FileNameJoin[{parentDir, "DiffExp"}]];

Print["==========================================="];
Print["Unequal-Mass Banana Parity/Speed Test"];
Print["===========================================\n"];

Get["DiffExp.m"];

testsPassed = 0;
testsTotal = 0;

pass[name_, condition_] := Module[{},
  testsTotal++;
  If[TrueQ[condition],
    testsPassed++;
    Print["  [PASS] ", name],
    Print["  [FAIL] ", name]
  ];
];

maxAbsDiff[a_, b_] := Max[Abs[Flatten[N[a - b, 80]]]];

validTransportQ[result_] :=
  AssociationQ[result] && KeyExistsQ[result, "SeriesValues"] && KeyExistsQ[result, "KinematicPoint"];

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

Print["Generating equal-mass boundary data at t = 1/2..."];
LoadConfiguration[{
  DeltaPrescriptions -> {t - 16 + I * \[Delta]},
  MatrixDirectory -> FileNameJoin[{scriptDir, "Banana_EqualMass_Matrices"}] <> "/",
  Verbosity -> 1,
  UseMobius -> True,
  UsePade -> True,
  UseRationalRecurrence -> False,
  WorkingPrecision -> 300,
  ChopPrecision -> 180,
  DivisionOrder -> 4,
  ExpansionOrder -> 50
}];

preparedBCs = PrepareBoundaryConditions[EqualMassBoundaryConditions, <|t -> -1/x|>];
boundaryAsymptotic = TransportTo[preparedBCs, <|t -> -1/x|>];
If[!validTransportQ[boundaryAsymptotic],
  Print["Failed to generate equal-mass asymptotic boundary data."];
  Quit[1];
];
boundaryAtHalf = TransportTo[boundaryAsymptotic, <|t -> 1/2|>];
If[!validTransportQ[boundaryAtHalf],
  Print["Failed to generate equal-mass boundary data at t = 1/2."];
  Quit[1];
];

UnequalMassBoundaryConditions = {
  <|psq -> 1/2, mm1 -> 1, mm2 -> 1, mm3 -> 1, mm4 -> 1|>,
  {
    boundaryAtHalf["SeriesValues"][[1]], boundaryAtHalf["SeriesValues"][[1]], boundaryAtHalf["SeriesValues"][[1]],
    boundaryAtHalf["SeriesValues"][[1]], boundaryAtHalf["SeriesValues"][[1]], boundaryAtHalf["SeriesValues"][[1]],
    boundaryAtHalf["SeriesValues"][[2]], boundaryAtHalf["SeriesValues"][[2]], boundaryAtHalf["SeriesValues"][[2]], boundaryAtHalf["SeriesValues"][[2]],
    boundaryAtHalf["SeriesValues"][[3]],
    boundaryAtHalf["SeriesValues"][[4]], boundaryAtHalf["SeriesValues"][[4]], boundaryAtHalf["SeriesValues"][[4]], boundaryAtHalf["SeriesValues"][[4]]
  }
};

unequalConfig = {
  ChopPrecision -> 180,
  DivisionOrder -> 4,
  ExpansionOrder -> 50,
  MatrixDirectory -> FileNameJoin[{scriptDir, "Banana_Matrices"}] <> "/",
  RadiusOfConvergence -> 10,
  UseMobius -> False,
  UsePade -> True,
  WorkingPrecision -> 300,
  Verbosity -> 1
};

targets = {
  <|
    "Name" -> "mass-ramp-three-quarter",
    "Line" -> <|psq -> 1/2, mm1 -> 1 + x, mm2 -> 1 + x/2, mm3 -> 1 + x/3, mm4 -> 1|>,
    "At" -> 3/4
  |>,
  <|
    "Name" -> "mass-ramp-endpoint",
    "Line" -> <|psq -> 1/2, mm1 -> 1 + x, mm2 -> 1 + x/2, mm3 -> 1 + x/3, mm4 -> 1|>,
    "At" -> 1
  |>
};

totalDefaultTime = 0;
totalRecurrenceTime = 0;
diffs = {};
speedups = {};

Do[
  Print["\nTarget: ", target["Name"]];

  LoadConfiguration[unequalConfig];
  UpdateConfiguration[{UseRationalRecurrence -> False, IntegrationStrategy -> "Default"}];
  defaultStart = AbsoluteTime[];
  defaultResult = TransportTo[UnequalMassBoundaryConditions, target["Line"], target["At"]];
  defaultTime = AbsoluteTime[] - defaultStart;
  If[!validTransportQ[defaultResult],
    Print["Standard solver failed for target ", target["Name"]];
    Quit[1];
  ];

  LoadConfiguration[unequalConfig];
  UpdateConfiguration[{UseRationalRecurrence -> True, IntegrationStrategy -> "Default"}];
  recurrenceStart = AbsoluteTime[];
  recurrenceResult = TransportTo[UnequalMassBoundaryConditions, target["Line"], target["At"]];
  recurrenceTime = AbsoluteTime[] - recurrenceStart;
  If[!validTransportQ[recurrenceResult],
    Print["Recurrence solver failed for target ", target["Name"]];
    Quit[1];
  ];

  diff = maxAbsDiff[defaultResult["SeriesValues"], recurrenceResult["SeriesValues"]];
  speedup = defaultTime/recurrenceTime;

  AppendTo[diffs, diff];
  AppendTo[speedups, speedup];
  totalDefaultTime += defaultTime;
  totalRecurrenceTime += recurrenceTime;

  Print["  default time:    ", NumberForm[defaultTime, 4], " s"];
  Print["  recurrence time: ", NumberForm[recurrenceTime, 4], " s"];
  Print["  speedup:         ", NumberForm[speedup, 3], "x"];
  Print["  max abs diff:    ", ScientificForm[diff, 3]];

  pass[target["Name"] <> " parity", diff < 10^-9];
  ,
  {target, targets}
];

Print["\n==========================================="];
Print["Summary"];
Print["==========================================="];
Print["  total default time:    ", NumberForm[totalDefaultTime, 4], " s"];
Print["  total recurrence time: ", NumberForm[totalRecurrenceTime, 4], " s"];
Print["  total speedup:         ", NumberForm[totalDefaultTime/totalRecurrenceTime, 3], "x"];
Print["  worst max abs diff:    ", ScientificForm[Max[diffs], 3]];

pass["all target parities", Max[diffs] < 10^-9];
pass["aggregate recurrence speedup", totalDefaultTime/totalRecurrenceTime > 1.15];

Print["\n==========================================="];
Print["Results: ", testsPassed, " / ", testsTotal, " tests passed"];
If[testsPassed === testsTotal,
  Print["All tests PASSED!"],
  Print["Some tests FAILED!"]
];
Print["==========================================="];

If[testsPassed === testsTotal, Quit[0], Quit[1]];
