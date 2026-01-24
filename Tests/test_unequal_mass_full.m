(* Full unequal mass banana test - tests actual transport computation *)
(* Runs with both Default and Recurrence strategies, compares results and timings *)

(* Set directory and paths *)
SetDirectory[DirectoryName[$InputFileName]];
scriptDir = Directory[];
parentDir = ParentDirectory[scriptDir];

(* Add subpackages to path *)
$Path = Prepend[$Path, FileNameJoin[{parentDir, "DiffExp"}]];

Print["==========================================="];
Print["Full Unequal Mass Banana Test"];
Print["(Default vs Recurrence comparison)"];
Print["===========================================\n"];

(* Load the refactored DiffExp package *)
Print["Loading refactored DiffExp package..."];
Get["DiffExp.m"];
Print["Package loaded!\n"];

(* ====== Step 1: Generate boundary conditions from equal mass system ====== *)
Print["=== Step 1: Generate boundary conditions at t = 1/2 using equal mass ==="];
EqualMassConfiguration = {
  DeltaPrescriptions -> {t - 16 + I * \[Delta]},
  MatrixDirectory -> FileNameJoin[{scriptDir, "Banana_EqualMass_Matrices"}] <> "/",
  Verbosity -> 1,
  UseMobius -> True,
  UsePade -> True
};

Print["Loading equal mass configuration..."];
LoadConfiguration[EqualMassConfiguration];

(* Equal mass boundary conditions *)
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
PreparedBCs = PrepareBoundaryConditions[EqualMassBoundaryConditions, <|t -> -1/x|>];

Print["Updating configuration for high precision..."];
UpdateConfiguration[{DivisionOrder -> 4, ExpansionOrder -> 70}];

Print["Transporting to asymptotic limit..."];
TmpAsymptotic = TransportTo[PreparedBCs, <|t -> -1/x|>];

Print["Transporting to t = 1/2..."];
BoundaryConditionsAtHalf = TransportTo[TmpAsymptotic, <|t -> 1/2|>];
Print["Equal mass boundary conditions at t = 1/2 generated.\n"];

(* ====== Step 2: Load unequal mass configuration ====== *)
Print["=== Step 2: Load unequal mass configuration ==="];
UnequalMassConfiguration = {
  ChopPrecision -> 500,
  DivisionOrder -> 4,
  ExpansionOrder -> 70,
  MatrixDirectory -> FileNameJoin[{scriptDir, "Banana_Matrices"}] <> "/",
  RadiusOfConvergence -> 10,
  UseMobius -> True,
  UsePade -> True,
  WorkingPrecision -> 1000,
  Verbosity -> 1
};

(* ====== Step 3: Prepare unequal mass boundary conditions ====== *)
(* The mapping for the unequal mass family:
   6 copies of integral 1, 4 copies of integral 2, 1 copy of integral 3, 4 copies of integral 4 *)
UnequalMassBoundaryConditions = {
  <|psq -> 1/2, mm1 -> 1, mm2 -> 1, mm3 -> 1, mm4 -> 1|>,
  {
    BoundaryConditionsAtHalf["SeriesValues"][[1]],
    BoundaryConditionsAtHalf["SeriesValues"][[1]],
    BoundaryConditionsAtHalf["SeriesValues"][[1]],
    BoundaryConditionsAtHalf["SeriesValues"][[1]],
    BoundaryConditionsAtHalf["SeriesValues"][[1]],
    BoundaryConditionsAtHalf["SeriesValues"][[1]],
    BoundaryConditionsAtHalf["SeriesValues"][[2]],
    BoundaryConditionsAtHalf["SeriesValues"][[2]],
    BoundaryConditionsAtHalf["SeriesValues"][[2]],
    BoundaryConditionsAtHalf["SeriesValues"][[2]],
    BoundaryConditionsAtHalf["SeriesValues"][[3]],
    BoundaryConditionsAtHalf["SeriesValues"][[4]],
    BoundaryConditionsAtHalf["SeriesValues"][[4]],
    BoundaryConditionsAtHalf["SeriesValues"][[4]],
    BoundaryConditionsAtHalf["SeriesValues"][[4]]
  }
};

(* ====== Step 4: Transport with DEFAULT strategy ====== *)
Print["=== Step 4a: Transport with Default strategy ==="];
Print["Target: mm1=2, mm2=3/2, mm3=4/3, mm4=1\n"];

LoadConfiguration[UnequalMassConfiguration];
UpdateConfiguration[{UseRationalRecurrence -> False, IntegrationStrategy -> "Default"}];

defaultStart = AbsoluteTime[];
DefaultResults = TransportTo[
  UnequalMassBoundaryConditions,
  <|psq -> 1/2, mm1 -> 1 + x, mm2 -> 1 + x/2, mm3 -> 1 + x/3, mm4 -> 1|>,
  1
];
defaultTime = AbsoluteTime[] - defaultStart;

Print["Default transport complete in ", NumberForm[defaultTime, 4], " seconds."];
Print["Result point: ", DefaultResults["KinematicPoint"], "\n"];

(* ====== Step 5: Transport with RECURRENCE strategy ====== *)
Print["=== Step 4b: Transport with Recurrence strategy ==="];

LoadConfiguration[UnequalMassConfiguration];
UpdateConfiguration[{UseRationalRecurrence -> True}];

recurrenceStart = AbsoluteTime[];
RecurrenceResults = TransportTo[
  UnequalMassBoundaryConditions,
  <|psq -> 1/2, mm1 -> 1 + x, mm2 -> 1 + x/2, mm3 -> 1 + x/3, mm4 -> 1|>,
  1
];
recurrenceTime = AbsoluteTime[] - recurrenceStart;

Print["Recurrence transport complete in ", NumberForm[recurrenceTime, 4], " seconds.\n"];

(* ====== Step 6: Compare results ====== *)
Print["=== Comparison ==="];

testsPassed = 0;
testsTotal = 0;

(* Check structure *)
testsTotal++;
If[KeyExistsQ[DefaultResults, "SeriesValues"] && KeyExistsQ[RecurrenceResults, "SeriesValues"],
  Print["  [PASS] Both results have SeriesValues"];
  testsPassed++;
  ,
  Print["  [FAIL] Missing SeriesValues in one or both results"];
];

(* Compare numerical values at eps^0 through eps^4 for all integrals *)
maxAbsDiff = 0;
maxRelDiff = 0;
Do[
  Do[
    defVal = DefaultResults["SeriesValues"][[integral, epsOrd + 1]];
    recVal = RecurrenceResults["SeriesValues"][[integral, epsOrd + 1]];
    (* Compare at x=0 (leading coefficient) *)
    If[Head[defVal] === SeriesData && Head[recVal] === SeriesData,
      defCoeffs = defVal[[3]] /. DiffExp`Symbols`Logx -> 0;
      recCoeffs = recVal[[3]] /. DiffExp`Symbols`Logx -> 0;
      If[Length[defCoeffs] > 0 && Length[recCoeffs] > 0,
        nComp = Min[Length[defCoeffs], Length[recCoeffs], 5];
        Do[
          If[NumericQ[defCoeffs[[c]]] && NumericQ[recCoeffs[[c]]],
            absDiff = Abs[defCoeffs[[c]] - recCoeffs[[c]]];
            maxAbsDiff = Max[maxAbsDiff, absDiff];
            If[Abs[defCoeffs[[c]]] > 10^(-100),
              maxRelDiff = Max[maxRelDiff, absDiff / Abs[defCoeffs[[c]]]];
            ];
          ];
          , {c, nComp}
        ];
      ];
    ];
    , {epsOrd, 0, 4}
  ];
  , {integral, Length[DefaultResults["SeriesValues"]]}
];

testsTotal++;
If[maxAbsDiff < 10^(-20),
  Print["  [PASS] Results agree (max abs diff: ", ScientificForm[maxAbsDiff, 3], ")"];
  testsPassed++;
  ,
  Print["  [FAIL] Results disagree (max abs diff: ", ScientificForm[maxAbsDiff, 3], ")"];
];

testsTotal++;
If[maxRelDiff < 10^(-20),
  Print["  [PASS] Relative agreement (max rel diff: ", ScientificForm[maxRelDiff, 3], ")"];
  testsPassed++;
  ,
  Print["  [FAIL] Relative disagreement (max rel diff: ", ScientificForm[maxRelDiff, 3], ")"];
];

(* ====== Timing comparison ====== *)
Print["\n=== Timing Comparison ==="];
Print["  Default strategy:    ", NumberForm[defaultTime, 4], " seconds"];
Print["  Recurrence strategy: ", NumberForm[recurrenceTime, 4], " seconds"];
If[recurrenceTime < defaultTime,
  Print["  Speedup: ", NumberForm[defaultTime / recurrenceTime, 3], "x faster with recurrence"];
  ,
  Print["  Ratio: ", NumberForm[recurrenceTime / defaultTime, 3], "x slower with recurrence"];
];

Print["\n==========================================="];
Print["Results: ", testsPassed, " / ", testsTotal, " tests passed"];
If[testsPassed === testsTotal,
  Print["All tests PASSED!"];
  ,
  Print["Some tests FAILED!"];
];
Print["==========================================="];
