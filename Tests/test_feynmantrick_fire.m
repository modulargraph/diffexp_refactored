(* Test script for FeynmanTrick FIRE interface *)
(* Tests DefineTopology, SetupFIRE, FindBasis with a simple 1-loop box *)

SetDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$Path = Prepend[$Path, Directory[]];

Print["=== Testing FeynmanTrick FIRE Interface ===\n"];

(* Load the package *)
Print["Loading FeynmanTrick package..."];
Get["FeynmanTrick.m"];
Print["Package loaded.\n"];

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

(* ============================================================ *)
(* Test 1: DefineTopology - 1-loop box *)
(* Convention: D_j = -q_j^2 + m_j^2 *)
(* ============================================================ *)

Print["--- Test: DefineTopology (1-loop box) ---"];

(* Note: l1, p1, p2, p3, s, t must be Global symbols for FIRE *)
Module[{topo},
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

  Print["\n--- Test: SetupFIRE ---"];

  Module[{workDir, setupTopo},
    workDir = FileNameJoin[{$TemporaryDirectory, "FTtest_box"}];
    If[DirectoryQ[workDir], DeleteDirectory[workDir, DeleteContents -> True]];

    Print["  Work directory: ", workDir];

    setupTopo = FeynmanTrick`FIREInterface`SetupFIRE[topo, workDir];

    test["SetupFIRE returns association", Head[setupTopo], Association];
    test["StartFileReady", setupTopo["StartFileReady"], True];
    test[".start file exists",
      FileExistsQ[FileNameJoin[{workDir, "box.start"}]], True];
    test[".config file exists",
      FileExistsQ[FileNameJoin[{workDir, "box.config"}]], True];

    (* Check config file contents *)
    Module[{configContent},
      configContent = Import[FileNameJoin[{workDir, "box.config"}], "Text"];
      Print["  Config file:\n", configContent];
      test["Config has #threads", StringContainsQ[configContent, "#threads"], True];
      test["Config has #variables", StringContainsQ[configContent, "#variables"], True];
    ];

    Print["\n--- Test: FindBasis ---"];

    Module[{basisTopo, masters},
      basisTopo = FeynmanTrick`FIREInterface`FindBasis[setupTopo];

      If[basisTopo =!= $Failed,
        masters = basisTopo["Masters"];
        Print["  Masters found: ", Length[masters]];
        Print["  Master indices: ", masters];

        (* 1-loop box should have a small number of masters *)
        (* In d=4-2eps with 4 massless propagators, there are typically *)
        (* a few master integrals depending on how FIRE counts them *)
        test["At least 1 master found", Length[masters] >= 1, True];
        test["Masters are lists of length 4",
          AllTrue[masters, Length[#] == 4 &], True];
      ,
        Print["FAIL: FindBasis returned $Failed"];
        failed++;
      ];
    ];

    (* Clean up *)
    If[DirectoryQ[workDir], DeleteDirectory[workDir, DeleteContents -> True]];
  ];
];

(* ============================================================ *)
(* Summary *)
(* ============================================================ *)

Print["\n============================"];
Print["Results: ", passed, " passed, ", failed, " failed."];
Print["============================"];

If[failed > 0, Exit[1]];
