(* End-to-end test for the Feynman trick pipeline *)
(* Tests the full bottom-up integration on the massless box *)
(* Combination: {1,2}, {1,3}, {1,4} at Euclidean point s=-1, t=-1/3 *)

SetDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$Path = Prepend[$Path, Directory[]];

Print["=== End-to-End FeynmanTrick Pipeline Test ===\n"];

(* Load the package *)
Print["Loading FeynmanTrick package..."];
Get["FeynmanTrick.m"];
Print["Package loaded.\n"];

(* Configuration *)
SetFTOption["Threads", 1];
SetFTOption["FThreads", 1];
SetFTOption["Verbosity", 2];

passed = 0;
failed = 0;

test[name_, expr_, expected_] := Module[{},
  If[expr === expected,
    Print["PASS: ", name]; passed++;,
    Print["FAIL: ", name];
    Print["  Expected: ", expected];
    Print["  Got:      ", expr]; failed++;
  ];
];

testTrue[name_, expr_] := Module[{},
  If[TrueQ[expr],
    Print["PASS: ", name]; passed++;,
    Print["FAIL: ", name, " (got: ", expr, ")"]; failed++;
  ];
];

testNumerical[name_, expr_, expected_, tol_:10^-10] := Module[{},
  If[NumericQ[expr] && NumericQ[expected] && Abs[expr - expected] < tol,
    Print["PASS: ", name, " (error: ", Abs[expr - expected] // ScientificForm, ")"]; passed++;,
    Print["FAIL: ", name];
    Print["  Expected: ", expected];
    Print["  Got:      ", expr]; failed++;
  ];
];


(* ============================================================ *)
(* Part 1: Test BoundaryConditions Module                        *)
(* ============================================================ *)

Print["\n--- Part 1: Boundary Conditions ---\n"];

(* Test rescaled Feynman parameters *)
Print["Testing RescaledFeynmanParameters..."];
Module[{rescaled, params, numVals},
  {rescaled, params} = FeynmanTrick`BoundaryConditions`RescaledFeynmanParameters[4];
  Print["  Rescaled params: ", rescaled];
  Print["  Variables: ", params];

  (* Verify they sum to 1 *)
  testTrue["Rescaled params sum to 1",
    Together[Total[rescaled] - 1] === 0];

  (* Evaluate at x_j = 11/23 *)
  numVals = rescaled /. Thread[params -> Table[11/23, 3]];
  Print["  At x_j = 11/23: ", numVals // N];
  testNumerical["Sum of numerical rescaled params", Total[numVals] // N, 1.0];
];


(* Test Symanzik polynomial computation *)
Print["\n--- Testing Symanzik Polynomials ---"];
Module[{props, loops, replacements, result, U, F, vars},
  (* Simple bubble: 2 propagators, 1 loop *)
  props = {-l1^2, -(l1 + p1)^2};
  loops = {l1};
  replacements = {p1^2 -> -1};  (* Euclidean: p1^2 = -1 *)

  result = FeynmanTrick`BoundaryConditions`ComputeSymanzikPolynomials[
    props, loops, replacements
  ];

  If[result =!= $Failed,
    {U, F, vars} = result;
    Print["  Bubble U = ", U];
    Print["  Bubble F = ", F];
    Print["  Variables: ", vars];
    testTrue["U is non-zero", U =!= 0];
    testTrue["F is non-zero", F =!= 0];
  ,
    Print["  Symanzik computation failed (FIRE6 UF not available)"];
    failed += 2;
  ];
];


(* Test tadpole boundary evaluation *)
Print["\n--- Testing Tadpole Boundary ---"];
Module[{minPow, coeffs, Fval},
  (* One-loop tadpole with M^2 = 1 *)
  (* I_1 = Gamma(1-d/2)/Gamma(1) * 1^(d/2-1) = Gamma(-1+eps) *)
  (* This has a pole at eps=0: Gamma(-1+eps) = -1/eps + (gamma-1) + O(eps) *)
  {minPow, coeffs} = FeynmanTrick`BoundaryConditions`EvaluateTadpoleBoundary[
    1, 1, 1, 1, 4  (* U=1, F=1, v=1, L=1, epsOrder=4 *)
  ];

  Print["  Min eps power: ", minPow];
  Print["  Leading coefficients: ", N[coeffs[[1;;3]]]];

  test["Tadpole with M^2=1 has pole", minPow, -1];
  testNumerical["Pole residue is -1", N[coeffs[[1]]], -1.0];
];

Module[{minPow, coeffs},
  (* One-loop tadpole with M^2 = 2, v = 2 *)
  (* I_2 = Gamma(2-d/2)/Gamma(2) * U^(2-d) / F^(2-d/2) *)
  (* = Gamma(eps)/1 * 1^(2eps) / 2^eps = Gamma(eps) * 2^(-eps) *)
  (* Gamma(eps) = 1/eps - gamma + O(eps) *)
  {minPow, coeffs} = FeynmanTrick`BoundaryConditions`EvaluateTadpoleBoundary[
    1, 2, 2, 1, 4  (* U=1, F=2, v=2, L=1, epsOrder=4 *)
  ];

  Print["  Tadpole v=2, M^2=2: min power = ", minPow];
  test["v=2 tadpole has pole", minPow, -1];
];


(* ============================================================ *)
(* Part 2: Test Level Setup and Matrix Computation               *)
(* ============================================================ *)

Print["\n\n--- Part 2: Multi-Level Iteration Setup ---\n"];

Module[{topo, ftData, workDir, exportDir, nLevels},

  (* Define the massless box topology *)
  topo = FeynmanTrick`FIREInterface`DefineTopology[
    "box_euclidean",
    {l1},
    {p1, p2, p3},
    {-l1^2, -(l1 + p1)^2, -(l1 + p1 + p2)^2, -(l1 + p1 + p2 + p3)^2},
    {p1^2 -> 0, p2^2 -> 0, p3^2 -> 0,
     p1*p2 -> s/2, p2*p3 -> t/2, p1*p3 -> -(s + t)/2}
  ];

  (* Define 3-level iteration: combine {1,2}, then {1,3}, then {1,4} *)
  ftData = FeynmanTrick`FeynmanTrickIteration`DefineFTIteration[
    topo,
    {{1, 2}, {1, 3}, {1, 4}},  (* combination sequence *)
    {s -> -1, t -> -1/3}  (* Euclidean kinematic point *)
  ];

  nLevels = ftData["NumLevels"];
  test["Three levels defined", nLevels, 3];

  (* Build all levels *)
  Print["\nBuilding all levels..."];
  ftData = FeynmanTrick`FeynmanTrickIteration`BuildAllLevels[ftData];

  testTrue["Level 1 built", KeyExistsQ[ftData["Levels"], 1]];
  testTrue["Level 2 built", KeyExistsQ[ftData["Levels"], 2]];
  testTrue["Level 3 built", KeyExistsQ[ftData["Levels"], 3]];

  (* Print propagator info at each level *)
  Do[
    Module[{levelProps},
      levelProps = ftData["Levels"][k]["Propagators"];
      Print["\n  Level ", k, " propagators (", Length[levelProps], " total):"];
      Do[Print["    D", j, " = ", levelProps[[j]]], {j, Length[levelProps]}];
    ],
    {k, nLevels}
  ];

  (* Compute matrices at each level *)
  Print["\n\n--- Computing matrices at each level ---"];
  workDir = FileNameJoin[{$TemporaryDirectory, "FT_endtoend_" <> ToString[$ProcessID]}];
  If[DirectoryQ[workDir], DeleteDirectory[workDir, DeleteContents -> True]];
  exportDir = FileNameJoin[{workDir, "matrices"}];

  Do[
    Print["\n  Computing level ", k, "..."];
    ftData = FeynmanTrick`FeynmanTrickIteration`ComputeLevelData[ftData, k];

    If[ftData["Levels"][k]["Computed"],
      Module[{masters, diffMat},
        masters = ftData["Levels"][k]["Masters"];
        diffMat = ftData["Levels"][k]["DiffMatrix"];
        Print["    Masters: ", masters];
        Print["    Matrix dimensions: ", Dimensions[diffMat]];
        testTrue["Level " <> ToString[k] <> " computed", True];
      ];
    ,
      Print["    FAILED to compute level ", k];
      failed++;
    ];
    ,
    {k, nLevels, 1, -1}  (* bottom-up order *)
  ];

  (* Export matrices *)
  Print["\n  Exporting matrices..."];
  Do[
    If[ftData["Levels"][k]["Computed"],
      FeynmanTrick`FeynmanTrickIteration`ExportLevel[ftData, k, exportDir, "diffexp", 2];
    ];
    ,
    {k, nLevels}
  ];

  (* Verify matrix files exist *)
  Do[
    Module[{levelDir},
      levelDir = FileNameJoin[{exportDir, "Level_" <> ToString[k] <> "_Matrices"}];
      testTrue["Level " <> ToString[k] <> " matrix files exist",
        DirectoryQ[levelDir] && FileExistsQ[FileNameJoin[{levelDir, "dxx_0.m"}]]
      ];
    ];
    ,
    {k, nLevels}
  ];


  (* ============================================================ *)
  (* Part 3: Test Boundary Computation at Deepest Level           *)
  (* ============================================================ *)

  Print["\n\n--- Part 3: Deepest Level Boundary ---\n"];

  Module[{boundary},
    boundary = FeynmanTrick`BoundaryConditions`DeepestLevelBoundary[ftData, 2];

    If[AssociationQ[boundary],
      Print["  Boundary computed successfully."];
      Print["  Eps prefactors: ", boundary["EpsPrefactors"]];
      Print["  Boundary values (first 3 eps orders): ",
        N[boundary["BoundaryValues"][[All, 1;;Min[3, Length[boundary["BoundaryValues"][[1]]]]]]]]];
      Print["  U value: ", boundary["Uval"]];
      Print["  F value: ", boundary["Fval"]];
      testTrue["Boundary is Association", True];
      testTrue["Has boundary values", Length[boundary["BoundaryValues"]] >= 1];
    ,
      Print["  Boundary computation failed: ", boundary];
      failed += 2;
    ];
  ];


  (* ============================================================ *)
  (* Part 4: Test DiffExp Transport at Deepest Level              *)
  (* ============================================================ *)

  Print["\n\n--- Part 4: DiffExp Transport ---\n"];

  Module[{boundary, matDir, transportResult},
    boundary = FeynmanTrick`BoundaryConditions`DeepestLevelBoundary[ftData, 2];

    If[AssociationQ[boundary],
      matDir = FileNameJoin[{exportDir, "Level_" <> ToString[nLevels] <> "_Matrices"}];

      Print["  Transporting level ", nLevels, " from xx=11/23..."];

      transportResult = FeynmanTrick`DiffExpIntegration`TransportLevel[
        matDir,
        boundary["BoundaryValues"],
        2,  (* epsOrder *)
        "WorkingPrecision" -> 100,
        "ExpansionOrder" -> 30,
        "Verbosity" -> 1
      ];

      If[AssociationQ[transportResult],
        Print["  Transport successful!"];
        Print["  Number of segments: ", Length[transportResult["SegmentData"]]];
        testTrue["Transport returned segments", Length[transportResult["SegmentData"]] > 0];
      ,
        Print["  Transport failed: ", transportResult];
        failed++;
      ];
    ,
      Print["  Skipping transport (boundary computation failed)"];
      failed++;
    ];
  ];

  (* Cleanup *)
  Print["\n  Work directory: ", workDir];
];


(* ============================================================ *)
(* Summary *)
(* ============================================================ *)

Print["\n\n============================"];
Print["Results: ", passed, " passed, ", failed, " failed."];
Print["============================"];

If[failed > 0, Exit[1]];
