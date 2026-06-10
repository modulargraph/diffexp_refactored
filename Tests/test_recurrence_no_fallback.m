(* Guardrail: recurrence mode must not silently use the old default path. *)

SetDirectory[DirectoryName[$InputFileName]];
scriptDir = Directory[];
parentDir = ParentDirectory[scriptDir];

Print["==========================================="];
Print["Recurrence No-Fallback Guardrail"];
Print["===========================================\n"];

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

read[path_] := ReadString[FileNameJoin[{parentDir, path}]];
containsQ[text_, needle_] := StringContainsQ[text, needle, IgnoreCase -> False];

dispatch = read["DiffExp/IntegrationStrategies/Dispatch.m"];
singular = read["DiffExp/IntegrationStrategies/Recurrence.m"];
generalSingular = read["DiffExp/IntegrationStrategies/ResonantRecurrence.m"];
state = read["DiffExp/State.m"];
ftDiffExpIntegration = read["FeynmanTrick/DiffExpIntegration.m"];

pass[
  "restricted singular recurrence does not call SolveDefault",
  !containsQ[singular, "SolveDefault["]
];

pass[
  "general singular recurrence does not call SolveDefault",
  !containsQ[generalSingular, "SolveDefault["]
];

pass[
  "dispatch has an explicit recurrence no-fallback error",
  containsQ[dispatch, "UseRationalRecurrence -> True was requested"] &&
    containsQ[dispatch, "Refusing to fall back"]
];

pass[
  "recurrence no-fallback branch precedes default branch",
  StringPosition[dispatch, "UseRationalRecurrence -> True was requested"][[1, 1]] <
    StringPosition[dispatch, "(* Default strategy *)"][[1, 1]]
];

pass[
  "UseRationalRecurrence usage documents exclusive recurrence behavior",
  containsQ[state, "no recursive strategy accepts a block"] &&
    containsQ[state, "instead of silently falling back"]
];

pass[
  "DiffExp defaults to recurrence mode",
  containsQ[state, "UseRationalRecurrence -> True"]
];

pass[
  "FeynmanTrick transport defaults to recurrence mode",
  containsQ[ftDiffExpIntegration, "\"UseRationalRecurrence\" -> True"]
];

Print["\n==========================================="];
Print["Tests passed: ", testsPassed, "/", testsTotal];
Print["==========================================="];

If[testsPassed === testsTotal,
  Quit[0],
  Quit[1]
];
