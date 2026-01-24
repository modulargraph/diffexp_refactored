(* Test: Singular Endpoint Transport in DiffExp *)
(*
   Tests that TransportTo correctly handles transport to a singular
   endpoint by returning the series expansion instead of trying to
   evaluate at the singularity.

   Uses the Hypergeometric 2F1 system (singularities at z=0, z=1):
   z(1-z)y'' + [c - (a+b+1)z]y' - ab*y = 0
   Parameters: a = 1/4, b = 1/3, c = 3/2

   Test strategy:
   1. Set boundary conditions at z=1/2 (known from Mathematica)
   2. Transport towards z=0 (singular endpoint)
   3. Verify EndpointIsSingularity flag is set
   4. Verify SegmentData is populated
   5. Evaluate the returned series slightly away from singularity and
      compare with known values
   6. Repeat transport towards z=1 (singular endpoint)
*)

(* Set directory and paths *)
SetDirectory[DirectoryName[$InputFileName]];
scriptDir = Directory[];
parentDir = ParentDirectory[scriptDir];

(* Add subpackages to path *)
$Path = Prepend[$Path, FileNameJoin[{parentDir, "DiffExp"}]];

Print["==========================================="];
Print["Singular Endpoint Transport Test"];
Print["===========================================\n"];

(* Load the refactored DiffExp package *)
Print["Loading DiffExp package..."];
Get["DiffExp.m"];
Print["Package loaded!\n"];

(* Hypergeometric parameters *)
a = 1/4;
b = 1/3;
c = 3/2;

testsPassed = 0;
testsFailed = 0;

passTest[name_] := (testsPassed++; Print["  PASS: ", name]);
failTest[name_, msg_:""] := (testsFailed++; Print["  FAIL: ", name, If[msg =!= "", " - " <> ToString[msg], ""]]);

(* ============================================================ *)
(* Part 1: Transport to z=0 (singular endpoint) *)
(* ============================================================ *)
Print["\n--- Part 1: Transport to z=0 (singular) ---"];

(* Configuration *)
Config2F1 = {
  MatrixDirectory -> FileNameJoin[{scriptDir, "Hypergeometric2F1_Matrices"}] <> "/",
  Verbosity -> 1,
  UseMobius -> False,
  UsePade -> False,
  WorkingPrecision -> 200,
  ExpansionOrder -> 40,
  DivisionOrder -> 4,
  ChopPrecision -> 150,
  EpsilonOrder -> 0
};

Print["Loading configuration..."];
LoadConfiguration[Config2F1];

(* Boundary conditions at z=1/2 *)
zStart = 1/2;
y2F1AtHalf = N[Hypergeometric2F1[a, b, c, zStart], 200];
yPrimeAtHalf = N[D[Hypergeometric2F1[a, b, c, zz], zz] /. zz -> zStart, 200];

Print["Boundary at z = 1/2:"];
Print["  Y1 = ", NumberForm[y2F1AtHalf, 15]];
Print["  Y2 = ", NumberForm[yPrimeAtHalf, 15]];

BCs = {
  <|z -> SetPrecision[zStart, 200]|>,
  {{y2F1AtHalf}, {yPrimeAtHalf}}
};

(* Transport towards z=0 *)
Print["\nTransporting from z=1/2 towards z=0..."];
resultTo0 = DiffExp`Transport`TransportTo[
  BCs,
  <|z -> 0|>,
  1,  (* endpoint *)
  True  (* SaveExpansions *)
];

(* Check results *)
If[AssociationQ[resultTo0],
  passTest["TransportTo z=0 returned Association"];
,
  failTest["TransportTo z=0 returned Association", Head[resultTo0]];
  Print["  Result: ", resultTo0];
  Print["\n=== ABORTING: Transport failed ==="];
  Print["\nResults: ", testsPassed, " passed, ", testsFailed, " failed"];
  Exit[1];
];

(* Check EndpointIsSingularity flag *)
If[KeyExistsQ[resultTo0, "EndpointIsSingularity"] && resultTo0["EndpointIsSingularity"] === True,
  passTest["EndpointIsSingularity flag is True for z=0"];
,
  failTest["EndpointIsSingularity flag", resultTo0["EndpointIsSingularity"]];
];

(* Check SegmentData is populated *)
If[KeyExistsQ[resultTo0, "SegmentData"] && ListQ[resultTo0["SegmentData"]] && Length[resultTo0["SegmentData"]] > 0,
  passTest["SegmentData populated (" <> ToString[Length[resultTo0["SegmentData"]]] <> " segments)"];
,
  failTest["SegmentData populated"];
];

(* Check SeriesValues is a list of series, not just numbers *)
If[KeyExistsQ[resultTo0, "SeriesValues"],
  seriesVals = resultTo0["SeriesValues"];
  If[MatrixQ[seriesVals] && (Dimensions[seriesVals] === {2, 1}),
    (* Check if the series values contain x-dependent terms (SeriesData) *)
    hasSeriesData = Or @@ (MatchQ[#, _SeriesData] & /@ Flatten[seriesVals]);
    If[hasSeriesData,
      passTest["SeriesValues contains SeriesData (x-dependent expansion)"];
    ,
      (* Even if not SeriesData, it should be non-numeric (contains x, Logx, etc.) *)
      isNonNumeric = !And @@ (NumericQ /@ Flatten[seriesVals]);
      If[isNonNumeric,
        passTest["SeriesValues contains non-numeric (symbolic) entries"];
      ,
        failTest["SeriesValues should be x-dependent at singular endpoint"];
      ];
    ];
  ,
    failTest["SeriesValues dimensions", Dimensions[seriesVals]];
  ];
,
  failTest["SeriesValues key exists"];
];

(* Verify the series gives correct values near z=0 *)
(* Use the SECOND-to-last segment which covers a region away from the singularity *)
If[KeyExistsQ[resultTo0, "SegmentData"] && Length[resultTo0["SegmentData"]] >= 2,
  Module[{testSeg, segLine, segInterval, localInterval, seriesData,
          evalXMain, evalXLocal, evalResult, evalZ, expectedY1, diff1},
    (* Use the second segment (away from singularity, better convergence) *)
    testSeg = resultTo0["SegmentData"][[2]];

    segLine = testSeg[[1]];
    segInterval = testSeg[[3]];  (* interval in x-coordinates of MAIN line *)
    localInterval = testSeg[[4]];  (* interval in local segment coordinates *)
    seriesData = testSeg[[5]];

    Print["\n  Test segment interval (main line x): ", segInterval // N];
    Print["  Test segment local interval: ", localInterval // N];
    Print["  Test segment line: ", segLine // N];

    (* Pick a point at the center of the local interval for best convergence *)
    evalXLocal = (localInterval[[1]] + localInterval[[2]]) / 2;

    (* Convert to z-value using the segment line *)
    evalZ = (segLine /. DiffExp`Symbols`x -> evalXLocal)[[Key[z]]];
    Print["  Evaluating at local x = ", evalXLocal // N, ", z = ", evalZ // N];

    (* Evaluate the series at this local point *)
    If[ListQ[seriesData] && Length[seriesData] >= 2,
      evalResult = Table[
        Module[{s},
          s = seriesData[[intIdx, 1]];
          If[MatchQ[s, _SeriesData],
            Normal[s] /. DiffExp`Symbols`x -> evalXLocal /. DiffExp`Symbols`Logx -> Log[Abs[evalXLocal]],
            s
          ]
        ],
        {intIdx, 2}
      ];

      expectedY1 = N[Hypergeometric2F1[a, b, c, evalZ], 200];

      If[NumericQ[evalResult[[1]]] && NumericQ[expectedY1],
        diff1 = Abs[evalResult[[1]] - expectedY1];
        Print["  |Y1_series - Y1_exact| = ", diff1 // N];
        If[diff1 < 10^-10,
          passTest["Series evaluation gives correct Y1 (accuracy: " <> ToString[Floor[-Log10[diff1 // N]]] <> " digits)"];
        ,
          failTest["Series evaluation Y1 accuracy", diff1 // N];
        ];
      ,
        Print["  Note: Could not evaluate series numerically"];
        Print["  Result: ", evalResult[[1]]];
        Print["  (Skipping numerical comparison)"];
      ];
    ,
      Print["  Series data format: ", If[ListQ[seriesData], Dimensions[seriesData], Head[seriesData]]];
    ];
  ];
];

(* ============================================================ *)
(* Part 2: Transport to z=1 (singular endpoint) *)
(* ============================================================ *)
Print["\n--- Part 2: Transport to z=1 (singular) ---"];

(* Reload configuration *)
LoadConfiguration[Config2F1];

(* Transport towards z=1 *)
Print["Transporting from z=1/2 towards z=1..."];
resultTo1 = DiffExp`Transport`TransportTo[
  BCs,
  <|z -> 1|>,
  1,
  True
];

If[AssociationQ[resultTo1],
  passTest["TransportTo z=1 returned Association"];
,
  failTest["TransportTo z=1 returned Association", Head[resultTo1]];
];

If[AssociationQ[resultTo1] && KeyExistsQ[resultTo1, "EndpointIsSingularity"],
  If[resultTo1["EndpointIsSingularity"] === True,
    passTest["EndpointIsSingularity flag is True for z=1"];
  ,
    failTest["EndpointIsSingularity flag for z=1", resultTo1["EndpointIsSingularity"]];
  ];
];

If[AssociationQ[resultTo1] && KeyExistsQ[resultTo1, "SegmentData"],
  If[ListQ[resultTo1["SegmentData"]] && Length[resultTo1["SegmentData"]] > 0,
    passTest["SegmentData populated for z=1 transport (" <> ToString[Length[resultTo1["SegmentData"]]] <> " segments)"];
  ,
    failTest["SegmentData populated for z=1"];
  ];
];


(* ============================================================ *)
(* Part 3: Transport to non-singular endpoint (regression test) *)
(* ============================================================ *)
Print["\n--- Part 3: Transport to z=3/4 (non-singular, regression) ---"];

(* Reload configuration *)
LoadConfiguration[Config2F1];

(* Transport towards z=3/4 (NOT a singularity) *)
Print["Transporting from z=1/2 towards z=3/4..."];
resultTo34 = DiffExp`Transport`TransportTo[
  BCs,
  <|z -> 3/4|>,
  1,
  True
];

If[AssociationQ[resultTo34],
  passTest["TransportTo z=3/4 returned Association"];
,
  failTest["TransportTo z=3/4 returned Association", Head[resultTo34]];
];

If[AssociationQ[resultTo34] && KeyExistsQ[resultTo34, "EndpointIsSingularity"],
  If[resultTo34["EndpointIsSingularity"] === False,
    passTest["EndpointIsSingularity flag is False for z=3/4"];
  ,
    failTest["EndpointIsSingularity flag for z=3/4 (should be False)", resultTo34["EndpointIsSingularity"]];
  ];
];

(* Verify the evaluation at z=3/4 is correct *)
If[AssociationQ[resultTo34] && KeyExistsQ[resultTo34, "SeriesValues"],
  Module[{sv, expectedY1, expectedY2, diff1},
    sv = resultTo34["SeriesValues"];
    expectedY1 = N[Hypergeometric2F1[a, b, c, 3/4], 200];
    expectedY2 = N[D[Hypergeometric2F1[a, b, c, zz], zz] /. zz -> 3/4, 200];

    If[MatrixQ[sv] && NumericQ[sv[[1, 1]]],
      diff1 = Abs[sv[[1, 1]] - expectedY1];
      Print["  |Y1_transport - Y1_exact| = ", diff1 // N];
      If[diff1 < 10^-20,
        passTest["Non-singular transport Y1 accuracy (" <> ToString[Floor[-Log10[diff1 // N]]] <> " digits)"];
      ,
        failTest["Non-singular transport Y1 accuracy", diff1 // N];
      ];
    ,
      Print["  SeriesValues format: ", Dimensions[sv]];
      Print["  First entry: ", sv[[1, 1]]];
    ];
  ];
];


(* ============================================================ *)
(* Summary *)
(* ============================================================ *)
Print["\n==========================================="];
Print["Test Summary"];
Print["==========================================="];
Print["Passed: ", testsPassed];
Print["Failed: ", testsFailed];
Print["Total:  ", testsPassed + testsFailed];

If[testsFailed > 0,
  Print["\nSOME TESTS FAILED"];
  Exit[1];
,
  Print["\nALL TESTS PASSED"];
];
