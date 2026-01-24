(* Test: Equal Mass Banana with General Singular Recurrence *)
(*
   Tests the general singular recurrence strategy on the equal mass banana integral.
   The 4x4 system has a regular singular point at t=infinity (x=0 with t=-1/x).

   Residue matrix M_0 at x=0 has:
   - Eigenvalue 1 with algebraic multiplicity 3 (Jordan block of size 3)
   - Eigenvalue 0 with multiplicity 1
   - Resonance: 0 and 1 differ by integer 1

   This triggers the full resonant machinery:
   - Jordan chains (non-diagonalizable)
   - Resonance between eigenvalue classes
   - Logarithmic terms in the fundamental matrix

   Comparison: results should match the Default strategy (Frobenius/Wronskian).
*)

(* Set directory and paths *)
SetDirectory[DirectoryName[$InputFileName]];
scriptDir = Directory[];
parentDir = ParentDirectory[scriptDir];

(* Add subpackages to path *)
$Path = Prepend[$Path, FileNameJoin[{parentDir, "DiffExp"}]];

Print["==========================================="];
Print["Equal Mass Banana - Resonant Recurrence Test"];
Print["===========================================\n"];

(* Load the refactored DiffExp package *)
Print["Loading DiffExp package..."];
Get["DiffExp.m"];
Print["Package loaded!\n"];

(* First, compute the reference result using the Default strategy *)
Print["=== Phase 1: Computing reference with Default strategy ===\n"];

DefaultStrategyConfig = {
  DeltaPrescriptions -> {t - 16 + I * \[Delta]},
  MatrixDirectory -> FileNameJoin[{scriptDir, "Banana_EqualMass_Matrices"}] <> "/",
  Verbosity -> 1,
  UseMobius -> True,
  UsePade -> True,
  IntegrationStrategy -> "Default",
  UseRationalRecurrence -> False
};

LoadConfiguration[DefaultStrategyConfig];

(* Boundary conditions *)
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

Print["Preparing boundary conditions..."];
PreparedBCsDefault = PrepareBoundaryConditions[EqualMassBoundaryConditions, <|t -> -1/x|>];

Print["Transporting to t = -1 with Default strategy..."];
defaultStart = AbsoluteTime[];
ResultsDefault = TransportTo[PreparedBCsDefault, <|t -> -1|>];
defaultTime = AbsoluteTime[] - defaultStart;
Print["Default transport complete in ", NumberForm[defaultTime, 4], " seconds.\n"];

(* Now compute with the general singular recurrence *)
Print["=== Phase 2: Computing with General Singular Recurrence ===\n"];

RecurrenceConfiguration = {
  DeltaPrescriptions -> {t - 16 + I * \[Delta]},
  MatrixDirectory -> FileNameJoin[{scriptDir, "Banana_EqualMass_Matrices"}] <> "/",
  Verbosity -> 3,
  UseMobius -> True,
  UsePade -> True,
  IntegrationStrategy -> "Default",
  UseRationalRecurrence -> True  (* This enables the recurrence strategies *)
};

LoadConfiguration[RecurrenceConfiguration];

Print["Preparing boundary conditions..."];
PreparedBCsRecurrence = PrepareBoundaryConditions[EqualMassBoundaryConditions, <|t -> -1/x|>];

Print["Transporting to t = -1 with General Singular Recurrence..."];
recurrenceStart = AbsoluteTime[];
ResultsRecurrence = TransportTo[PreparedBCsRecurrence, <|t -> -1|>];
recurrenceTime = AbsoluteTime[] - recurrenceStart;
Print["Recurrence transport complete in ", NumberForm[recurrenceTime, 4], " seconds.\n"];

(* Compare results *)
Print["=== Comparison ==="];
Print["Default result point: ", ResultsDefault["KinematicPoint"]];
Print["Recurrence result point: ", ResultsRecurrence["KinematicPoint"]];

testsPassed = 0;
testsTotal = 0;

(* Check that both results have the right structure *)
testsTotal++;
If[AssociationQ[ResultsRecurrence] && KeyExistsQ[ResultsRecurrence, "SeriesValues"],
  Print["  [PASS] Recurrence transport returned expected structure"];
  testsPassed++;
  ,
  Print["  [FAIL] Recurrence transport returned unexpected structure"];
];

(* Compare numerical values *)
defaultValues = ResultsDefault["SeriesValues"];
recurrenceValues = ResultsRecurrence["SeriesValues"];

Print["\nDefault values at t=-1:"];
Print[TableForm[N[defaultValues, 10]]];
Print["\nRecurrence values at t=-1:"];
Print[TableForm[N[recurrenceValues, 10]]];

(* Calculate maximum relative difference *)
maxRelDiff = 0;
maxAbsDiff = 0;
Do[
  d = defaultValues[[i, j]];
  r = recurrenceValues[[i, j]];
  absDiff = Abs[N[d - r, 50]];
  If[NumericQ[absDiff],
    maxAbsDiff = Max[maxAbsDiff, absDiff];
    If[Abs[d] > 10^-50,
      relDiff = absDiff / Abs[d];
      If[NumericQ[relDiff], maxRelDiff = Max[maxRelDiff, relDiff]];
    ];
  ];
  , {i, Length[defaultValues]}, {j, Length[defaultValues[[1]]]}
];

Print["\nMax absolute difference: ", N[maxAbsDiff, 5]];
Print["Max relative difference: ", N[maxRelDiff, 5]];

testsTotal++;
If[maxAbsDiff < 10^-10 || maxRelDiff < 10^-10,
  Print["  [PASS] Recurrence and Default produce matching results"];
  testsPassed++;
  ,
  Print["  [FAIL] Recurrence and Default results differ significantly"];
];

(* Timing comparison *)
Print["\n=== Timing Comparison ==="];
Print["  Default strategy:    ", NumberForm[defaultTime, 4], " seconds"];
Print["  Recurrence strategy: ", NumberForm[recurrenceTime, 4], " seconds"];
If[recurrenceTime < defaultTime,
  Print["  Speedup: ", NumberForm[defaultTime / recurrenceTime, 3], "x faster with recurrence"];
  ,
  Print["  Ratio: ", NumberForm[recurrenceTime / defaultTime, 3], "x slower with recurrence"];
];

(* Summary *)
Print["\n==========================================="];
Print["Results: ", testsPassed, " / ", testsTotal, " tests passed"];
If[testsPassed === testsTotal,
  Print["All tests PASSED!"];,
  Print["Some tests FAILED!"];
];
Print["==========================================="];
