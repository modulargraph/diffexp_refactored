(* test_feynmantrick_pipeline.m *)
(* End-to-end test of the Feynman trick integration pipeline *)
(* Uses a 1-loop massless box as the simplest non-trivial example *)
(* Combination sequence: {1,2},{1,3},{1,4} -> triangle -> bubble -> tadpole *)

Print["=== FeynmanTrick Pipeline Test ==="];
Print[""];

(* --- Setup --- *)
SetDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$Path = Prepend[$Path, Directory[]];

(* Load the FeynmanTrick package *)
Quiet[Get["FeynmanTrick/FeynmanTrick.m"], {General::shdw, Symbol::shdw}];

(* Also load DiffExp for transport/integration *)
Quiet[Get["DiffExp.m"], {General::shdw, Symbol::shdw}];

(* Configure *)
SetFTOption["Threads", 1];
SetFTOption["FThreads", 1];
SetFTOption["Verbosity", 2];
SetFTOption["WorkingPrecision", 200];

(* --- Test Framework --- *)
passed = 0;
failed = 0;

test[name_String, expr_, expected_] := If[expr === expected,
  passed++; Print["  PASS: ", name],
  failed++; Print["  FAIL: ", name, " (got ", expr, ", expected ", expected, ")"]
];

testTrue[name_String, expr_] := If[TrueQ[expr],
  passed++; Print["  PASS: ", name],
  failed++; Print["  FAIL: ", name, " (got ", expr, ")"]
];

testNumerical[name_String, expr_, expected_, tol_:10^-10] := If[
  NumericQ[expr] && NumericQ[expected] && Abs[expr - expected] < tol,
  passed++; Print["  PASS: ", name, " (", expr, ")"],
  failed++; Print["  FAIL: ", name, " (got ", expr, ", expected ", expected, ")"]
];


(* ============================================================ *)
(* Part 1: Define the 1-loop massless box topology               *)
(* ============================================================ *)
Print["\n--- Part 1: Define Topology ---"];

Module[{topology, ftData, outputDir, epsOrder, matrixEpsOrder, boundaryEpsOrder, workDir},

  (* 1-loop massless box: 4 propagators *)
  (* Convention: D_j = -q_j^2 + m_j^2 (following code convention) *)
  (* For massless box: D_j = -q_j^2 *)
  topology = FeynmanTrick`FIREInterface`DefineTopology[
    "box1L",
    {Global`l1},                                    (* loop momenta *)
    {Global`p1, Global`p2, Global`p3},              (* external momenta *)
    {-Global`l1^2,                                  (* D1 = -l1^2 *)
     -(Global`l1 + Global`p1)^2,                    (* D2 = -(l1+p1)^2 *)
     -(Global`l1 + Global`p1 + Global`p2)^2,        (* D3 = -(l1+p1+p2)^2 *)
     -(Global`l1 + Global`p1 + Global`p2 + Global`p3)^2},  (* D4 = -(l1+p1+p2+p3)^2 *)
    {Global`p1^2 -> 0,                              (* massless legs *)
     Global`p2^2 -> 0,
     Global`p3^2 -> 0,
     Global`p1 Global`p2 -> Global`s/2,             (* Mandelstam variables *)
     Global`p2 Global`p3 -> Global`t/2,
     Global`p1 Global`p3 -> -(Global`s + Global`t)/2}
  ];

  testTrue["Topology defined", AssociationQ[topology]];
  test["Topology has 4 propagators", topology["NumPropagators"], 4];

  (* Define the iteration: combine {1,2}, then {1,3}, then {1,4} *)
  (* This gives: box -> triangle -> bubble -> tadpole *)
  ftData = FeynmanTrick`FeynmanTrickIteration`DefineFTIteration[
    topology,
    {{1, 2}, {1, 3}, {1, 4}},      (* combination sequence *)
    {Global`s -> -1, Global`t -> -1/3}  (* Euclidean kinematic point *)
  ];

  testTrue["FT iteration defined", AssociationQ[ftData]];
  test["3 levels defined", ftData["NumLevels"], 3];

  (* Set up output directory *)
  outputDir = FileNameJoin[{$TemporaryDirectory, "FT_pipeline_test_" <> ToString[$ProcessID]}];
  epsOrder = 2;  (* Keep eps order low for speed *)
  boundaryEpsOrder = epsOrder + 3;
  matrixEpsOrder = boundaryEpsOrder + 1;

  Print["\n  Output directory: ", outputDir];
  Print["  Epsilon order: ", epsOrder];


  (* ============================================================ *)
  (* Part 2: Build all levels and compute matrices                 *)
  (* ============================================================ *)
  Print["\n--- Part 2: Build Levels & Compute Matrices ---"];

  (* Build all levels *)
  ftData = FeynmanTrick`FeynmanTrickIteration`BuildAllLevels[ftData];

  testTrue["Level 1 built", KeyExistsQ[ftData["Levels"], 1]];
  testTrue["Level 2 built", KeyExistsQ[ftData["Levels"], 2]];
  testTrue["Level 3 built", KeyExistsQ[ftData["Levels"], 3]];

  (* Compute each level (bottom-up: 3, 2, 1) *)
  Do[
    Print["\n  Computing level ", level, "..."];
    ftData = FeynmanTrick`FeynmanTrickIteration`ComputeLevelData[ftData, level];
    If[ftData["Levels"][level]["Computed"],
      passed++;
      Print["  PASS: Level ", level, " computed (",
            Length[ftData["Levels"][level]["Masters"]], " masters)"];
    ,
      failed++;
      Print["  FAIL: Level ", level, " computation failed"];
    ];
  , {level, 3, 1, -1}];

  (* Export matrices *)
  Do[
    FeynmanTrick`FeynmanTrickIteration`ExportLevel[ftData, level, outputDir, "diffexp", matrixEpsOrder];
  , {level, 3, 1, -1}];

  (* Check matrix files exist - each level has its own parameter (xx1, xx2, xx3) *)
  Do[
    Module[{matDir, matFile},
      matDir = FileNameJoin[{outputDir, "Level_" <> ToString[level] <> "_Matrices"}];
      matFile = FileNameJoin[{matDir, "dxx" <> ToString[level] <> "_0.m"}];
      testTrue["Level " <> ToString[level] <> " matrix exported",
        FileExistsQ[matFile]];
    ];
  , {level, 1, 3}];


  (* ============================================================ *)
  (* Part 3: Compute deepest level boundary                        *)
  (* ============================================================ *)
  Print["\n--- Part 3: Deepest Level Boundary ---"];

  Module[{deepBoundary},
    deepBoundary = FeynmanTrick`BoundaryConditions`DeepestLevelBoundary[ftData, boundaryEpsOrder];

    If[deepBoundary =!= $Failed && AssociationQ[deepBoundary],
      passed++;
      Print["  PASS: Deepest level boundary computed"];
      Print["    Number of masters: ", Length[deepBoundary["Masters"]]];
      Print["    Eps prefactors: ", deepBoundary["EpsPrefactors"]];
      Print["    U value: ", deepBoundary["Uval"]];
      Print["    F value: ", deepBoundary["Fval"]];

      testTrue["Boundary values are lists",
        AllTrue[deepBoundary["BoundaryValues"], ListQ]];
      testTrue["Boundary values are numeric",
        AllTrue[Flatten[deepBoundary["BoundaryValues"]], NumericQ]];

      (* Store for next phase *)
      $deepBoundary = deepBoundary;
    ,
      failed += 3;
      Print["  FAIL: Deepest level boundary failed"];
      $deepBoundary = $Failed;
    ];
  ];


  (* ============================================================ *)
  (* Part 4: Transport at deepest level                            *)
  (* ============================================================ *)
  Print["\n--- Part 4: DiffExp Transport ---"];

  If[$deepBoundary =!= $Failed,
    Module[{matDir, transportResult, currentBCs},
      currentBCs = $deepBoundary["BoundaryValues"];
      matDir = FileNameJoin[{outputDir, "Level_3_Matrices"}];

      Print["  Transporting level 3 (", Length[currentBCs], " masters)..."];

      transportResult = FeynmanTrick`DiffExpIntegration`TransportLevel[
        matDir, currentBCs, Lookup[$deepBoundary, "WorkingEpsilonOrder", epsOrder],
        "WorkingPrecision" -> 200,
        "ExpansionOrder" -> 30,
        "Verbosity" -> 1,
        "EpsPrefactors" -> $deepBoundary["EpsPrefactors"]
      ];

      If[transportResult =!= $Failed && AssociationQ[transportResult],
        passed++;
        Print["  PASS: Transport completed (",
              Length[transportResult["SegmentData"]], " segments)"];

        testTrue["Has segment data", Length[transportResult["SegmentData"]] > 0];
        test["NumIntegrals matches", transportResult["NumIntegrals"],
             Length[currentBCs]];

        $transportResult3 = transportResult;
      ,
        failed += 3;
        Print["  FAIL: Transport failed at level 3"];
        $transportResult3 = $Failed;
      ];
    ];
  ,
    failed += 4;
    Print["  SKIP: Transport skipped (boundary failed)"];
    $transportResult3 = $Failed;
  ];


  (* ============================================================ *)
  (* Part 5: Integration (level 3 -> level 2 boundary)            *)
  (* ============================================================ *)
  Print["\n--- Part 5: Feynman Trick Integration ---"];

  If[$transportResult3 =!= $Failed,
    Module[{levelBoundary, currentBCs},
      (* Store boundary values in transport result for "direct" case *)
      $transportResult3["BoundaryValuesAbove"] = $deepBoundary["BoundaryValues"];
      $transportResult3["EpsPrefactorsAbove"] = $deepBoundary["EpsPrefactors"];

      Print["  Computing level 2 boundary from level 3 transport..."];
      Print["  Level 2 masters: ", ftData["Levels"][2]["Masters"]];
      Print["  Level 3 masters: ", ftData["Levels"][3]["Masters"]];
      Print["  Combined positions at level 3: ", ftData["Levels"][3]["CombinedPositions"]];

      levelBoundary = FeynmanTrick`DiffExpIntegration`ComputeLevelBoundary[
        ftData, 2, $transportResult3, epsOrder
      ];

      If[AssociationQ[levelBoundary],
        passed++;
        Print["  PASS: Level 2 boundary computed"];
        Print["    Boundary values: ", levelBoundary["BoundaryValues"]];
        Print["    Eps prefactors: ", levelBoundary["EpsPrefactors"]];
        If[KeyExistsQ[levelBoundary, "RawBoundaryValues"],
          Print["    Raw eps range: ",
            levelBoundary["RawMinPower"], " ... ", levelBoundary["RawMaxPower"]];
        ];

        testTrue["Level 2 boundary has values",
          Length[levelBoundary["BoundaryValues"]] > 0];
        testTrue["Level 2 boundary values are numeric",
          AllTrue[Flatten[levelBoundary["BoundaryValues"]], NumericQ]];
        testTrue["Level 2 limit boundary is nonzero",
          Length[levelBoundary["BoundaryValues"]] >= 2 &&
            AnyTrue[levelBoundary["BoundaryValues"][[2]],
              !TrueQ[PossibleZeroQ[#]] &
            ]];

        $level2BCs = levelBoundary["BoundaryValues"];
        $level2Prefactors = levelBoundary["EpsPrefactors"];
      ,
        failed += 3;
        Print["  FAIL: ComputeLevelBoundary failed"];
        Print["    Result: ", levelBoundary];
        $level2BCs = $Failed;
        $level2Prefactors = {};
      ];
    ];
  ,
    failed += 4;
    Print["  SKIP: Integration skipped (transport failed)"];
    $level2BCs = $Failed;
  ];


  (* ============================================================ *)
  (* Part 6: Continue pipeline (level 2 -> level 1 -> level 0)    *)
  (* ============================================================ *)
  Print["\n--- Part 6: Complete Pipeline ---"];

  If[$level2BCs =!= $Failed,
    Module[{matDir, transportResult, levelBoundary, currentBCs, currentPrefactors},
      currentBCs = $level2BCs;
      currentPrefactors = $level2Prefactors;

      (* Transport level 2 *)
      matDir = FileNameJoin[{outputDir, "Level_2_Matrices"}];
      FeynmanTrick`FeynmanTrickIteration`ExportLevel[
        ftData, 2, outputDir, "diffexp", Length[First[currentBCs]] - 1
      ];
      Print["  Transporting level 2 (", Length[currentBCs], " masters)..."];

      transportResult = FeynmanTrick`DiffExpIntegration`TransportLevel[
        matDir, currentBCs, Length[First[currentBCs]] - 1,
        "WorkingPrecision" -> 200,
        "ExpansionOrder" -> 30,
        "Verbosity" -> 1,
        "EpsPrefactors" -> currentPrefactors
      ];

      If[transportResult =!= $Failed,
        passed++;
        Print["  PASS: Level 2 transport complete"];

        (* Integrate to get level 1 boundary *)
        transportResult["BoundaryValuesAbove"] = currentBCs;
        transportResult["EpsPrefactorsAbove"] = currentPrefactors;
        levelBoundary = FeynmanTrick`DiffExpIntegration`ComputeLevelBoundary[
          ftData, 1, transportResult, epsOrder
        ];

        If[AssociationQ[levelBoundary],
          currentBCs = levelBoundary["BoundaryValues"];
          currentPrefactors = Lookup[
            levelBoundary,
            "EpsPrefactors",
            Table[0, {Length[currentBCs]}]
          ];
          passed++;
          Print["  PASS: Level 1 boundary computed (", Length[currentBCs], " masters)"];

          (* Transport level 1 *)
          matDir = FileNameJoin[{outputDir, "Level_1_Matrices"}];
          FeynmanTrick`FeynmanTrickIteration`ExportLevel[
            ftData, 1, outputDir, "diffexp", Length[First[currentBCs]] - 1
          ];
          Print["  Transporting level 1 (", Length[currentBCs], " masters)..."];

          transportResult = FeynmanTrick`DiffExpIntegration`TransportLevel[
            matDir, currentBCs, Length[First[currentBCs]] - 1,
            "WorkingPrecision" -> 200,
            "ExpansionOrder" -> 30,
            "Verbosity" -> 1,
            "EpsPrefactors" -> currentPrefactors
          ];

          If[transportResult =!= $Failed,
            passed++;
            Print["  PASS: Level 1 transport complete"];

            (* Integrate to get level 0 boundary (final result!) *)
            transportResult["BoundaryValuesAbove"] = currentBCs;
            transportResult["EpsPrefactorsAbove"] = currentPrefactors;
            levelBoundary = FeynmanTrick`DiffExpIntegration`ComputeLevelBoundary[
              ftData, 0, transportResult, epsOrder
            ];

            If[AssociationQ[levelBoundary],
              passed++;
              Print["  PASS: Level 0 boundary computed (FINAL RESULT)"];
              Print[""];
              Print["  === FINAL RESULTS ==="];
              Print["  Masters at level 0: ", ftData["Levels"][0]["Masters"]];
              Do[
                Print["  Master ", i, ": ", levelBoundary["BoundaryValues"][[i]]];
              , {i, Length[levelBoundary["BoundaryValues"]]}];
              If[KeyExistsQ[levelBoundary, "RawBoundaryValues"],
                Print["  Raw eps range: ",
                  levelBoundary["RawMinPower"], " ... ", levelBoundary["RawMaxPower"]];
                Print["  Final eps prefactors: ", levelBoundary["EpsPrefactors"]];
                testTrue["Final result has epsilon pole",
                  TrueQ[levelBoundary["RawMinPower"] < 0]];
              ];
              testTrue["Final result is nonzero",
                AnyTrue[Flatten[levelBoundary["BoundaryValues"]],
                  !TrueQ[PossibleZeroQ[#]] &
                ]];
              Print["  =================="];
            ,
              failed++;
              Print["  FAIL: Level 0 boundary failed"];
            ];
          ,
            failed++;
            Print["  FAIL: Level 1 transport failed"];
          ];
        ,
          failed++;
          Print["  FAIL: Level 1 boundary computation failed"];
        ];
      ,
        failed++;
        Print["  FAIL: Level 2 transport failed"];
      ];
    ];
  ,
    failed += 4;
    Print["  SKIP: Pipeline continuation skipped"];
  ];


  (* Cleanup message *)
  Print["\n  Work directory preserved at: ", outputDir];
];


(* --- Summary --- *)
Print["\n=== Test Summary ==="];
Print["Passed: ", passed];
Print["Failed: ", failed];
Print["Total: ", passed + failed];

If[failed > 0,
  Print["\nSome tests FAILED."];
  Exit[1];
,
  Print["\nAll tests PASSED!"];
];
