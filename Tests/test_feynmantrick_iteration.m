(* Test script for FeynmanTrick iteration *)
(* Tests the Feynman trick combination and diff matrix w.r.t. xx *)

SetDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$Path = Prepend[$Path, Directory[]];

Print["=== Testing FeynmanTrick Iteration ===\n"];

(* Load the package *)
Print["Loading FeynmanTrick package..."];
Get["FeynmanTrick.m"];
Print["Package loaded.\n"];

(* Reduce threads for stability *)
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

(* ============================================================ *)
(* Test: Feynman trick on the massless box *)
(*                                                               *)
(* Combine D1 and D2 (propagators 1,2) using the Feynman trick  *)
(* The combined propagator: xx*D1 + (1-xx)*D2                   *)
(* Topology still has 4 propagators:                             *)
(*   {D12_combined, D2, D3, D4}                                 *)
(* Compute diff matrix w.r.t. xx at a numerical kinematic point  *)
(* ============================================================ *)

Print["--- Test: BuildLevel (combine props 1,2 of box) ---\n"];

Module[{topo, ftData, level1Data, level1Topo, level1Props,
        workDir, setupTopo, basisTopo, masters, diffMat,
        decomp, shiftedInts, exportDir},

  (* Define the box topology *)
  topo = FeynmanTrick`FIREInterface`DefineTopology[
    "box",
    {l1},
    {p1, p2, p3},
    {-l1^2, -(l1+p1)^2, -(l1+p1+p2)^2, -(l1+p1+p2+p3)^2},
    {p1^2 -> 0, p2^2 -> 0, p3^2 -> 0,
     p1*p2 -> s/2, p2*p3 -> t/2, p1*p3 -> -(s+t)/2}
  ];

  (* Define iteration: combine propagators 1 and 2 *)
  ftData = FeynmanTrick`FeynmanTrickIteration`DefineFTIteration[
    topo,
    {{1, 2}},  (* combine D1 and D2 *)
    {s -> 1, t -> -1/3}  (* numerical kinematic point *)
  ];

  test["DefineFTIteration created", Head[ftData], Association];
  test["NumLevels", ftData["NumLevels"], 1];

  (* Build level 1 *)
  Print["\nBuilding level 1..."];
  ftData = FeynmanTrick`FeynmanTrickIteration`BuildLevel[ftData, 1];

  level1Data = ftData["Levels"][1];
  level1Topo = level1Data["Topology"];
  level1Props = level1Data["Propagators"];

  Print["  Level 1 propagators:"];
  Do[Print["    D", k, " = ", level1Props[[k]]], {k, Length[level1Props]}];

  testTrue["Level 1 built",
    KeyExistsQ[ftData["Levels"], 1]];
  testTrue["Combined positions recorded",
    level1Data["CombinedPositions"] === {1, 2}];
  (* Each level has a unique parameter: xx1 for level 1, xx2 for level 2, etc. *)
  testTrue["Feynman parameter is xx1",
    level1Data["FeynmanParameter"] === Global`xx1];
  testTrue["Propagator 1 contains xx1",
    !FreeQ[level1Props[[1]], Global`xx1]];

  (* Check the decomposition fast path *)
  Print["\n--- Test: FeynmanTrickDecomposition ---"];
  decomp = FeynmanTrick`PropagatorAlgebra`FeynmanTrickDecomposition[4, 1, 2, FeynmanTrick`xx];
  Print["  coeffMatrix[1,1] = ", decomp[[1]][[1, 1]]];
  Print["  coeffMatrix[1,2] = ", decomp[[1]][[1, 2]]];
  test["Fast decomp coeff[1,1]", decomp[[1]][[1, 1]], 1/FeynmanTrick`xx];
  test["Fast decomp coeff[1,2]", decomp[[1]][[1, 2]], -1/FeynmanTrick`xx];
  testTrue["Other coeffs zero",
    Total[Abs[Flatten[Delete[decomp[[1]], {1}]]]] === 0];

  (* Test DifferentiatedIntegrals with the fast decomposition *)
  Print["\n--- Test: DifferentiatedIntegrals with Feynman decomp ---"];
  shiftedInts = FeynmanTrick`PropagatorAlgebra`DifferentiatedIntegrals[
    {1, 1, 1, 1}, decomp[[1]], decomp[[2]]
  ];
  Print["  Shifted integrals for {1,1,1,1}:"];
  Do[Print["    ", s], {s, shiftedInts}];
  testTrue["Has shifted integrals", Length[shiftedInts] > 0];
  (* d/dxx I_{1,1,1,1}: only prop 1 contributes (coeff 1/xx for itself, -1/xx for D2) *)
  (* -v_1 * [c_{1,1} * I_{unchanged} + c_{1,2} * I_{2,0,1,1}] *)
  (* = -1 * [1/xx * I_{1,1,1,1} + (-1/xx) * I_{2,0,1,1}] *)
  (* = {  {1,1,1,1}, -1/xx  } and { {2,0,1,1}, 1/xx } *)
  testTrue["Contains {1,1,1,1} term",
    MemberQ[shiftedInts[[All, 1]], {1, 1, 1, 1}]];
  testTrue["Contains {2,0,1,1} term",
    MemberQ[shiftedInts[[All, 1]], {2, 0, 1, 1}]];

  (* Now run FIRE to find basis and compute matrix *)
  Print["\n--- Test: SetupFIRE + FindBasis for combined topology ---"];

  workDir = FileNameJoin[{$TemporaryDirectory, "FTtest_iter_" <> ToString[$ProcessID]}];
  If[DirectoryQ[workDir], DeleteDirectory[workDir, DeleteContents -> True]];
  Print["  Work directory: ", workDir];

  setupTopo = FeynmanTrick`FIREInterface`SetupFIRE[level1Topo, workDir];

  If[setupTopo["StartFileReady"],
    Print["  FIRE setup complete."];

    basisTopo = FeynmanTrick`FIREInterface`FindBasis[setupTopo];

    If[basisTopo =!= $Failed,
      masters = basisTopo["Masters"];
      Print["  Masters: ", masters];
      testTrue["Masters found", Length[masters] >= 1];

      (* Compute diff matrix using the fast path *)
      Print["\n--- Test: ComputeDiffMatrix with Feynman fast path ---"];
      diffMat = FeynmanTrick`FIREInterface`ComputeDiffMatrix[basisTopo, FeynmanTrick`xx, decomp];

      If[diffMat =!= $Failed,
        Print["  Matrix dimensions: ", Dimensions[diffMat]];
        Print["  Matrix entries:"];
        Do[
          Do[
            If[diffMat[[ii, jj]] =!= 0,
              Print["    A[", ii, ",", jj, "] = ", diffMat[[ii, jj]]]
            ];
          , {jj, Length[diffMat[[1]]]}];
        , {ii, Length[diffMat]}];
        testTrue["Matrix is square",
          Length[diffMat] == Length[masters]];
        testTrue["Matrix has entries with xx",
          !FreeQ[diffMat, FeynmanTrick`xx]];

        (* Test export *)
        Print["\n--- Test: ExportLevel ---"];
        exportDir = FileNameJoin[{workDir, "export_test"}];

        (* Need to update ftData with computation results *)
        ftData["Levels"][1]["Topology"] = basisTopo;
        ftData["Levels"][1]["Masters"] = masters;
        ftData["Levels"][1]["DiffMatrix"] = diffMat;
        ftData["Levels"][1]["Computed"] = True;

        FeynmanTrick`FeynmanTrickIteration`ExportLevel[ftData, 1, exportDir, "both", 2];
        testTrue["Export directory created", DirectoryQ[exportDir]];
        testTrue["Full matrix exported",
          FileExistsQ[FileNameJoin[{exportDir, "Level_1_Matrices", "dxx1_full.m"}]]];
      ,
        Print["FAIL: ComputeDiffMatrix returned $Failed"];
        failed += 3;
      ];
    ,
      Print["FAIL: FindBasis returned $Failed (FIRE6 crashed)"];
      failed += 4;
    ];
  ,
    Print["FAIL: SetupFIRE failed"];
    failed += 5;
  ];

  Print["\n  (Work directory preserved: ", workDir, ")"];
];

(* ============================================================ *)
(* Summary *)
(* ============================================================ *)

Print["\n============================"];
Print["Results: ", passed, " passed, ", failed, " failed."];
Print["============================"];

If[failed > 0, Exit[1]];
