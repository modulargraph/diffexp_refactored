(* Test script for FeynmanTrick FIRE interface *)
(* Tests DefineTopology, SetupFIRE, FindBasis, ReduceIntegrals, ComputeDiffMatrix *)

SetDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$Path = Prepend[$Path, Directory[]];

Print["=== Testing FeynmanTrick FIRE Interface ===\n"];

(* Load the package *)
Print["Loading FeynmanTrick package..."];
Get["FeynmanTrick.m"];
Print["Package loaded.\n"];

(* Reduce threads for stability on macOS *)
SetFTOption["Threads", 1];
SetFTOption["FThreads", 1];
SetFTOption["Verbosity", 2];

passed = 0;
failed = 0;

test[name_, expr_, expected_] := Module[{result},
  result = (expr === expected);
  If[result,
    Print["PASS: ", name];
    passed++;
  ,
    Print["FAIL: ", name];
    Print["  Expected: ", expected];
    Print["  Got:      ", expr];
    failed++;
  ];
];

testTrue[name_, expr_] := Module[{},
  If[TrueQ[expr],
    Print["PASS: ", name];
    passed++;
  ,
    Print["FAIL: ", name];
    Print["  Expected: True"];
    Print["  Got:      ", expr];
    failed++;
  ];
];

(* ============================================================ *)
(* Test 1: DefineTopology - 1-loop box *)
(* Convention: D_j = -q_j^2 + m_j^2 *)
(* ============================================================ *)

Print["--- Test: DefineTopology (1-loop box) ---"];

Module[{topo, workDir, setupTopo, basisTopo, masters, reductions, diffMat},

  topo = FeynmanTrick`FIREInterface`DefineTopology[
    "box",
    {l1},                     (* loop momenta *)
    {p1, p2, p3},             (* external momenta *)
    {                         (* propagators: -q^2 + m^2 convention *)
      -l1^2,                  (* D1 *)
      -(l1+p1)^2,            (* D2 *)
      -(l1+p1+p2)^2,         (* D3 *)
      -(l1+p1+p2+p3)^2       (* D4 *)
    },
    {p1^2 -> 0, p2^2 -> 0, p3^2 -> 0,
     p1*p2 -> s/2, p2*p3 -> t/2, p1*p3 -> -(s+t)/2}
  ];

  test["Topology created", Head[topo], Association];
  test["Name", topo["Name"], "box"];
  test["NumPropagators", topo["NumPropagators"], 4];
  test["LoopMomenta", topo["LoopMomenta"], {l1}];
  test["ExternalMomenta", topo["ExternalMomenta"], {p1, p2, p3}];

  (* ============================================================ *)
  Print["\n--- Test: SetupFIRE ---"];

  workDir = FileNameJoin[{$TemporaryDirectory, "FTtest_full_" <> ToString[$ProcessID]}];
  If[DirectoryQ[workDir], DeleteDirectory[workDir, DeleteContents -> True]];
  Print["  Work directory: ", workDir];

  setupTopo = FeynmanTrick`FIREInterface`SetupFIRE[topo, workDir];

  test["SetupFIRE returns association", Head[setupTopo], Association];
  test["StartFileReady", setupTopo["StartFileReady"], True];
  test["ProblemNumber assigned", IntegerQ[setupTopo["ProblemNumber"]], True];
  testTrue[".start file exists",
    FileExistsQ[FileNameJoin[{workDir, "box.start"}]]];
  testTrue["setup fingerprint records the exact prepared start",
    AssociationQ[setupTopo["SetupFingerprintRecord"]] &&
    TrueQ[setupTopo["SetupFingerprintRecord", "VerifiedStartFile"]] &&
    StringQ[setupTopo["SetupFingerprintRecord", "StartFileSHA256"]] &&
    StringLength[setupTopo["SetupFingerprintRecord", "StartFileSHA256"]] === 64];
  testTrue["setup fingerprint digest matches its exact record",
    setupTopo["SetupFingerprint"] === IntegerString[
      Hash[setupTopo["SetupFingerprintRecord"], "SHA256"], 16, 64]];
  (* Note: .config file is now written when FIRE is run, not during setup *)

  (* ============================================================ *)
  Print["\n--- Test: FindBasis ---"];

  basisTopo = FeynmanTrick`FIREInterface`FindBasis[setupTopo];

  If[basisTopo =!= $Failed,
    masters = basisTopo["Masters"];
    Print["  Masters found: ", Length[masters]];
    Print["  Master indices: ", masters];

    testTrue["At least 1 master found", Length[masters] >= 1];
    testTrue["Masters are lists of length 4",
      AllTrue[masters, Length[#] == 4 &]];

    (* Known: massless box in d dimensions has 3 masters:
       box {1,1,1,1}, and two bubbles *)

    (* ============================================================ *)
    Print["\n--- Test: ReduceIntegrals ---"];

    (* Reduce a dotted integral *)
    reductions = FeynmanTrick`FIREInterface`ReduceIntegrals[
      basisTopo,
      {{1, 1, 2, 1}}  (* box with D3 raised *)
    ];

    If[reductions =!= $Failed,
      Print["  Reduction of {1,1,2,1}: ", reductions[{1, 1, 2, 1}]];
      testTrue["Reduction succeeded", AssociationQ[reductions]];
      testTrue["Result contains G terms",
        !FreeQ[reductions[{1, 1, 2, 1}], Global`G]];
    ,
      Print["FAIL: ReduceIntegrals returned $Failed"];
      failed++;
    ];

    (* ============================================================ *)
    Print["\n--- Test: ComputeDiffMatrix ---"];

    (* Compute diff matrix w.r.t. s *)
    diffMat = FeynmanTrick`FIREInterface`ComputeDiffMatrix[basisTopo, s];

    If[diffMat =!= $Failed,
      Print["  Matrix dimensions: ", Dimensions[diffMat]];
      Print["  Matrix:\n", MatrixForm[diffMat]];
      testTrue["Matrix is square",
        Length[diffMat] == Length[masters] && Length[diffMat[[1]]] == Length[masters]];
      testTrue["Matrix has nonzero entries",
        !AllTrue[Flatten[diffMat], # === 0 &]];
    ,
      Print["FAIL: ComputeDiffMatrix returned $Failed"];
      failed++;
    ];
  ,
    Print["FAIL: FindBasis returned $Failed"];
    Print["  (FIRE7 may have failed - check logs in ", workDir, ")"];
    failed += 5; (* count missed subtests *)
  ];

  (* Clean up *)
  (* If[DirectoryQ[workDir], DeleteDirectory[workDir, DeleteContents -> True]]; *)
  Print["\n  (Work directory preserved: ", workDir, ")"];
];

(* ============================================================ *)
(* Summary *)
(* ============================================================ *)

Print["\n============================"];
Print["Results: ", passed, " passed, ", failed, " failed."];
Print["============================"];

If[failed > 0, Exit[1]];
