(* ::Package:: *)
(* FIREInterface - Wrapper for FIRE6 IBP reduction *)
(* Handles topology definition, FIRE setup, basis finding, *)
(* integral reduction, and differential matrix computation *)
(*
   ARCHITECTURE: Multi-topology batch processing

   Instead of clearing and restarting FIRE for each topology, we:
   1. Generate all .start files first (SetupFIRE or SetupFIREBatch)
   2. Load all topologies with unique problem numbers
   3. Run FIRE6 once with all problems in a single config
   4. Parse combined results

   This avoids fermat subprocess race conditions and is how FIRE6 is
   designed to work with multiple related topologies.
*)

BeginPackage["FeynmanTrick`FIREInterface`", {"FeynmanTrick`", "FeynmanTrick`PropagatorAlgebra`"}];

DefineTopology::usage =
  "DefineTopology[name, loopMomenta, externalMomenta, propagators, replacements] \
defines a topology for FIRE reduction. Propagator convention: D_j = -q_j^2 + m_j^2. \
Returns a Topology association.";

SetupFIRE::usage =
  "SetupFIRE[topology, workDir] generates .start file for FIRE6. \
Does NOT run FIRE - call FindBasis or ReduceIntegrals after setup.";

SetupFIREBatch::usage =
  "SetupFIREBatch[{topo1, topo2, ...}, workDir] sets up multiple topologies \
at once, assigning each a unique problem number. Returns list of updated topologies.";

FindBasis::usage =
  "FindBasis[topology] runs FIRE6 to determine master integrals. \
Returns the topology with Masters field populated.";

FindBasisBatch::usage =
  "FindBasisBatch[{topo1, topo2, ...}] finds master integrals for all topologies \
in a single FIRE6 run. Returns list of updated topologies.";

ReduceIntegrals::usage =
  "ReduceIntegrals[topology, integrals] reduces a list of integrals \
(given as index vectors) to the master basis using FIRE6. \
Returns an association: integral -> linear combination of masters.";

ReduceIntegralsBatch::usage =
  "ReduceIntegralsBatch[{{topo1, ints1}, {topo2, ints2}, ...}] reduces \
integrals for multiple topologies in a single FIRE6 run.";

ComputeDiffMatrix::usage =
  "ComputeDiffMatrix[topology, variable] computes the differential equation \
matrix A such that d/dx Masters = A . Masters. \
The variable is typically the Feynman parameter xx.";

GetFIREResult::usage =
  "GetFIREResult[topology, integral] returns the reduction of a single integral.";

TopologyQ::usage = "TopologyQ[expr] returns True if expr is a valid topology association.";

ClearFIREState::usage = "ClearFIREState[] clears all FIRE internal state. Call before setting up new topologies.";

Begin["`Private`"];

(* Global state: track which topologies have been set up *)
$SetupTopologies = <||>;  (* name -> topology association *)
$NextProblemNumber = 1;

(* ============================================================ *)
(* DefineTopology                                                *)
(* ============================================================ *)

DefineTopology[name_String, loopMomenta_List, externalMomenta_List,
               propagators_List, replacements_List:{}] :=
  <|
    "Name" -> name,
    "LoopMomenta" -> loopMomenta,
    "ExternalMomenta" -> externalMomenta,
    "Propagators" -> propagators,
    "Replacements" -> replacements,
    "NumPropagators" -> Length[propagators],
    "WorkDirectory" -> "",
    "ProblemNumber" -> 0,  (* Assigned during SetupFIRE *)
    "Masters" -> {},
    "MasterRules" -> {},
    "StartFileReady" -> False,
    "Variables" -> extractVariables[propagators, replacements, loopMomenta, externalMomenta]
  |>;

TopologyQ[t_Association] := KeyExistsQ[t, "Propagators"] && KeyExistsQ[t, "LoopMomenta"];
TopologyQ[_] := False;


(* Extract symbolic variables from propagators and replacements *)
(* Returns Global` symbols suitable for FIRE/Fermat *)
extractVariables[props_List, repls_List, loops_List, exts_List] :=
Module[{allSyms, momenta, vars},
  allSyms = Cases[{props, repls}, _Symbol, Infinity] // DeleteDuplicates;
  momenta = Join[loops, exts];
  vars = Complement[allSyms, momenta];
  (* Return as Global symbols for FIRE compatibility *)
  Symbol[SymbolName[#]] & /@ vars
];

(* Build substitution rules mapping contexted symbols to Global equivalents *)
(* This ensures Fermat (which FIRE uses internally) gets simple variable names *)
buildFIRESubstitution[topology_Association] :=
Module[{allSyms, momenta, contextedSyms},
  allSyms = Cases[topology["Propagators"], _Symbol, Infinity] // DeleteDuplicates;
  momenta = Join[topology["LoopMomenta"], topology["ExternalMomenta"]];
  (* Find symbols not in Global` context that aren't momenta *)
  contextedSyms = Select[
    Complement[allSyms, momenta],
    (Context[#] =!= "Global`") &
  ];
  (* Map each to its Global equivalent *)
  Rule[#, Symbol[SymbolName[#]]] & /@ contextedSyms
];


(* ============================================================ *)
(* ClearFIREState - reset all FIRE state                        *)
(* ============================================================ *)

ClearFIREState[] := Module[{},
  $SetupTopologies = <||>;
  $NextProblemNumber = 1;

  (* Load FIRE6.m if needed *)
  ensureFIRELoaded[];

  (* Clear FIRE internal state *)
  Quiet[
    Unprotect[FIRE`Internal, FIRE`External, FIRE`Propagators, FIRE`Replacements];
    Clear[FIRE`Internal, FIRE`External, FIRE`Propagators, FIRE`Replacements];
    (* Clear all problem-indexed data *)
    Clear[FIRE`ExampleDimension, FIRE`SBasis0L, FIRE`SBasis0D, FIRE`SBasis0C,
          FIRE`SBasisL, FIRE`SBasisS, FIRE`SBasisR, FIRE`SBasisRL, FIRE`SBasisM, FIRE`HPI];
    FIRE`Burning = False;
  , {Unprotect::ssym, Clear::ssym}];
];


(* ============================================================ *)
(* ensureFIRELoaded - load FIRE6.m if not already               *)
(* ============================================================ *)

ensureFIRELoaded[] := Module[{firePath},
  If[!MemberQ[$Packages, "FIRE`"],
    firePath = FeynmanTrick`Private`$FTConfig["FIREPath"];
    Block[{$Path = Prepend[$Path, firePath]},
      Quiet[Get["FIRE6.m"], {General::shdw}];
    ];
  ];
];


(* ============================================================ *)
(* runFIRE6 - helper to run FIRE6 binary robustly               *)
(* ============================================================ *)

runFIRE6[fireBin_String, dir_String, configName_String] :=
Module[{exitCode, maxRetries = 3, attempt},
  Do[
    exitCode = runFIRE6Once[fireBin, dir, configName];
    If[exitCode === 0, Return[0, Module]];
    If[attempt < maxRetries,
      If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
        Print["FIRE6 failed (exit ", exitCode, "), retrying (attempt ", attempt+1, "/", maxRetries, ")..."];
      ];
      (* Clean up before retry *)
      Quiet[Run["pkill -9 -f fer64 2>/dev/null"]];
      Quiet[Run["pkill -9 -f FIRE6 2>/dev/null"]];
      Quiet[DeleteDirectory[FileNameJoin[{dir, "temp"}], DeleteContents -> True]];
      Pause[2];
    ];
  , {attempt, maxRetries}];
  If[exitCode =!= 0,
    Print["FIRE6 failed after ", maxRetries, " attempts (last exit code: ", exitCode, ")."];
  ];
  exitCode
];

runFIRE6Once[fireBin_String, dir_String, configName_String] :=
Module[{exitCode, result, logFile, oldDir, cmd, stdoutFile, stderrFile},
  logFile = FileNameJoin[{dir, "fire_stdout.log"}];
  stdoutFile = FileNameJoin[{dir, "fire_stdout.tmp"}];
  stderrFile = FileNameJoin[{dir, "fire_stderr.tmp"}];

  (* Clean temp directory *)
  Quiet[DeleteDirectory[FileNameJoin[{dir, "temp"}], DeleteContents -> True]];

  (* Run FIRE6 using shell execution - works better with Rosetta on M1 Macs *)
  (* RunProcess has issues spawning x86_64 child processes (fermat) from arm64 *)
  oldDir = Directory[];
  SetDirectory[dir];

  cmd = StringJoin[
    "cd ", dir, " && ",
    fireBin, " -c ", configName,
    " > ", stdoutFile, " 2> ", stderrFile
  ];

  exitCode = Run[cmd];
  (* Run returns exit code * 256 on Unix *)
  exitCode = If[IntegerQ[exitCode], BitShiftRight[exitCode, 8], exitCode];

  (* Read output *)
  result = <|
    "ExitCode" -> exitCode,
    "StandardOutput" -> If[FileExistsQ[stdoutFile], ReadString[stdoutFile], ""],
    "StandardError" -> If[FileExistsQ[stderrFile], ReadString[stderrFile], ""]
  |>;

  (* Clean up temp files *)
  Quiet[DeleteFile[stdoutFile]];
  Quiet[DeleteFile[stderrFile]];

  SetDirectory[oldDir];

  (* Save output to log file *)
  Export[logFile, result["StandardOutput"] <> result["StandardError"], "Text"];

  exitCode = result["ExitCode"];

  (* Print output on verbose or failure *)
  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 2 || exitCode =!= 0,
    If[StringLength[result["StandardOutput"]] > 0,
      Print["FIRE6 output (tail): "];
      Print[StringTake[result["StandardOutput"], -Min[500, StringLength[result["StandardOutput"]]]]];
    ];
    If[exitCode =!= 0,
      Print["FIRE6 exit code: ", exitCode];
      If[StringLength[result["StandardError"]] > 0,
        Print["FIRE6 stderr: ", StringTake[result["StandardError"], -Min[300, StringLength[result["StandardError"]]]]];
      ];
    ];
  ];

  exitCode
];


(* ============================================================ *)
(* SetupFIRE - Generate .start file for a single topology       *)
(* ============================================================ *)

SetupFIRE[topology_Association, workDir_String:""] :=
Module[{dir, name, result, fireSubst, fireProps, fireRepls, pn},
  name = topology["Name"];

  (* Determine working directory *)
  dir = If[workDir === "",
    FileNameJoin[{FeynmanTrick`Private`$FTConfig["WorkDirectory"], name}],
    workDir
  ];
  If[!DirectoryQ[dir], CreateDirectory[dir, CreateIntermediateDirectories -> True]];

  (* Load FIRE6.m if needed *)
  ensureFIRELoaded[];

  (* Assign problem number *)
  pn = $NextProblemNumber++;

  (* Substitute contexted symbols with Global equivalents for FIRE/Fermat *)
  fireSubst = buildFIRESubstitution[topology];
  fireProps = topology["Propagators"] /. fireSubst;
  fireRepls = topology["Replacements"] /. fireSubst;

  (* CRITICAL: Completely clear all FIRE state before setting up new topology *)
  (* Without this, FIRE reuses stale IBP relations and hangs *)
  Quiet[
    (* Unprotect everything that FIRE protects *)
    Unprotect[FIRE`Internal, FIRE`External, FIRE`Propagators, FIRE`Replacements,
              FIRE`PrepareIBPd, FIRE`BackMatrix, FIRE`Squares, FIRE`startinglist];
    (* Clear the "already prepared" flag and cached IBP data - this is critical! *)
    Clear[FIRE`PrepareIBPd, FIRE`BackMatrix, FIRE`Squares, FIRE`startinglist];
    Clear[FIRE`Internal, FIRE`External, FIRE`Propagators, FIRE`Replacements];
    (* Clear IBP and sector data *)
    Clear[FIRE`ExampleDimension, FIRE`SBasis0L, FIRE`SBasis0D, FIRE`SBasis0C,
          FIRE`SBasisL, FIRE`SBasisD, FIRE`SBasisA, FIRE`SBasisH, FIRE`SBasisO,
          FIRE`SBasisC, FIRE`SBasisS, FIRE`SBasisR, FIRE`SBasisRL, FIRE`SBasisM,
          FIRE`SBasisN, FIRE`HPI, FIRE`LRules];
    FIRE`Burning = False;
    FIRE`ProblemNumber = 0;
  , {Unprotect::ssym, Clear::ssym, Set::write}];

  (* Set FIRE variables *)
  FIRE`Internal = topology["LoopMomenta"];
  FIRE`External = topology["ExternalMomenta"];
  FIRE`Propagators = fireProps;
  If[fireRepls =!= {}, FIRE`Replacements = fireRepls];

  (* Prepare IBP relations *)
  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["Setting up FIRE for topology: ", name, " (problem ", pn, ")"];
  ];
  FIRE`PrepareIBP[];

  (* Prepare sector basis - this should now say "Prepared" not "Already prepared" *)
  FIRE`Prepare[AutoDetectRestrictions -> True];

  (* Save start file *)
  Module[{oldDir = Directory[]},
    SetDirectory[dir];
    FIRE`SaveStart[name];
    SetDirectory[oldDir];
  ];

  (* Return updated topology *)
  result = topology;
  result["WorkDirectory"] = dir;
  result["ProblemNumber"] = pn;
  result["StartFileReady"] = True;

  (* Store in global registry *)
  $SetupTopologies[name] = result;

  result
];


(* ============================================================ *)
(* SetupFIREBatch - Set up multiple topologies at once          *)
(* ============================================================ *)

SetupFIREBatch[topologies_List, workDir_String:""] :=
Module[{results},
  (* Clear state before batch setup *)
  ClearFIREState[];

  (* Set up each topology *)
  results = SetupFIRE[#, workDir] & /@ topologies;

  results
];


(* ============================================================ *)
(* writeMultiProblemConfig - config file for multiple problems  *)
(* ============================================================ *)

writeMultiProblemConfig[topologies_List, dir_String, configName_String,
                        integralsFile_String, outputFile_String] :=
Module[{vars, content, threads, fthreads, problemLines, allVars},
  threads = FeynmanTrick`Private`$FTConfig["Threads"];
  fthreads = FeynmanTrick`Private`$FTConfig["FThreads"];

  (* Collect all variables from all topologies *)
  allVars = Join[
    {FeynmanTrick`Private`$FTConfig["DimensionVariable"]},
    Flatten[#["Variables"] & /@ topologies]
  ] // DeleteDuplicates;

  (* Build #problem lines *)
  problemLines = StringJoin[
    "#problem           ", ToString[#["ProblemNumber"]], " ",
    #["Name"], ".start\n"
  ] & /@ topologies;

  content = StringJoin[
    "#threads           ", ToString[threads], "\n",
    "#fthreads          ", ToString[fthreads], "\n",
    "#variables         ", StringRiffle[SymbolName /@ allVars, ","], "\n",
    "#start\n",
    Sequence @@ problemLines,
    "#integrals         ", integralsFile, "\n",
    "#output            ", outputFile, "\n"
  ];

  Export[FileNameJoin[{dir, configName}], content, "Text"];
];


(* Write single-problem config (for backwards compatibility) *)
writeSingleProblemConfig[topology_Association, dir_String, configName_String,
                         integralsFile_String:"", outputFile_String:""] :=
Module[{vars, content, threads, fthreads, intFile, outFile, name, pn},
  name = topology["Name"];
  pn = topology["ProblemNumber"];
  threads = FeynmanTrick`Private`$FTConfig["Threads"];
  fthreads = FeynmanTrick`Private`$FTConfig["FThreads"];

  vars = Join[
    {FeynmanTrick`Private`$FTConfig["DimensionVariable"]},
    topology["Variables"]
  ] // DeleteDuplicates;

  intFile = If[integralsFile === "", name <> ".m", integralsFile];
  outFile = If[outputFile === "", name <> ".tables", outputFile];

  content = StringJoin[
    "#threads           ", ToString[threads], "\n",
    "#fthreads          ", ToString[fthreads], "\n",
    "#variables         ", StringRiffle[SymbolName /@ vars, ","], "\n",
    "#start\n",
    "#problem           ", ToString[pn], " ", name, ".start\n",
    "#integrals         ", intFile, "\n",
    "#output            ", outFile, "\n"
  ];

  Export[FileNameJoin[{dir, configName}], content, "Text"];
];


(* ============================================================ *)
(* FindBasis - find master integrals for a single topology      *)
(* ============================================================ *)

FindBasis[topology_Association] :=
Module[{dir, name, fireBin, intFile, result, masters, tablesFile, pn},

  If[!topology["StartFileReady"],
    Print["Error: Must call SetupFIRE before FindBasis."];
    Return[$Failed];
  ];

  dir = topology["WorkDirectory"];
  name = topology["Name"];
  pn = topology["ProblemNumber"];
  fireBin = FileNameJoin[{FeynmanTrick`Private`$FTConfig["FIREPath"], "bin", "FIRE6"}];

  (* Write integrals file with all sector corners *)
  Module[{nP = topology["NumPropagators"], allSectors, intContent},
    allSectors = Rest[Tuples[{0, 1}, nP]];
    intContent = StringJoin["{",
      StringRiffle[("{" <> ToString[pn] <> "," <> ToString[#, InputForm] <> "}") & /@ allSectors, ",\n"],
      "}\n"];
    intFile = FileNameJoin[{dir, name <> ".m"}];
    Export[intFile, intContent, "Text"];
  ];

  (* Write config *)
  writeSingleProblemConfig[topology, dir, name <> ".config"];

  (* Run FIRE6 *)
  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["Running FIRE6 to find basis for ", name, "..."];
  ];
  runFIRE6[fireBin, dir, name];

  (* Load tables and extract masters *)
  tablesFile = FileNameJoin[{dir, name <> ".tables"}];
  If[!FileExistsQ[tablesFile],
    Print["Error: FIRE6 did not produce output tables."];
    Return[$Failed];
  ];

  masters = FIRE`Tables2Masters[tablesFile];
  (* Filter for this problem number and extract indices *)
  masters = Cases[masters, {pn, indices_List} :> indices];

  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["Found ", Length[masters], " master integrals."];
  ];

  result = topology;
  result["Masters"] = masters;
  result
];


(* ============================================================ *)
(* FindBasisBatch - find masters for multiple topologies        *)
(* ============================================================ *)

FindBasisBatch[topologies_List] :=
Module[{dir, fireBin, intContent, allSectors, tablesFile, masters, results},

  If[Length[topologies] == 0, Return[{}]];

  (* Use first topology's work directory *)
  dir = topologies[[1]]["WorkDirectory"];
  fireBin = FileNameJoin[{FeynmanTrick`Private`$FTConfig["FIREPath"], "bin", "FIRE6"}];

  (* Build combined integrals file *)
  intContent = StringJoin["{",
    StringRiffle[
      Flatten[
        Table[
          With[{topo = topologies[[i]], pn = topologies[[i]]["ProblemNumber"]},
            allSectors = Rest[Tuples[{0, 1}, topo["NumPropagators"]]];
            ("{" <> ToString[pn] <> "," <> ToString[#, InputForm] <> "}") & /@ allSectors
          ],
          {i, Length[topologies]}
        ],
        1
      ],
      ",\n"
    ],
    "}\n"
  ];
  Export[FileNameJoin[{dir, "batch_basis.m"}], intContent, "Text"];

  (* Write multi-problem config *)
  writeMultiProblemConfig[topologies, dir, "batch_basis.config",
                          "batch_basis.m", "batch_basis.tables"];

  (* Run FIRE6 once *)
  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["Running FIRE6 to find basis for ", Length[topologies], " topologies..."];
  ];
  runFIRE6[fireBin, dir, "batch_basis"];

  (* Load results *)
  tablesFile = FileNameJoin[{dir, "batch_basis.tables"}];
  If[!FileExistsQ[tablesFile],
    Print["Error: FIRE6 did not produce output tables."];
    Return[$Failed];
  ];

  masters = FIRE`Tables2Masters[tablesFile];

  (* Update each topology with its masters *)
  results = Table[
    With[{topo = topologies[[i]], pn = topologies[[i]]["ProblemNumber"]},
      Module[{topoMasters, result},
        topoMasters = Cases[masters, {pn, indices_List} :> indices];
        If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
          Print["  ", topo["Name"], ": ", Length[topoMasters], " masters"];
        ];
        result = topo;
        result["Masters"] = topoMasters;
        result
      ]
    ],
    {i, Length[topologies]}
  ];

  results
];


(* ============================================================ *)
(* ReduceIntegrals - reduce integrals for a single topology     *)
(* ============================================================ *)

ReduceIntegrals[topology_Association, integrals_List] :=
Module[{dir, name, fireBin, intFile, tablesFile, rules, intContent, result, pn},

  If[!topology["StartFileReady"],
    Print["Error: Must call SetupFIRE before ReduceIntegrals."];
    Return[$Failed];
  ];

  dir = topology["WorkDirectory"];
  name = topology["Name"];
  pn = topology["ProblemNumber"];
  fireBin = FileNameJoin[{FeynmanTrick`Private`$FTConfig["FIREPath"], "bin", "FIRE6"}];

  (* Write integrals file with problem number *)
  intContent = StringJoin[
    "{",
    StringRiffle[
      ("{" <> ToString[pn] <> "," <> ToString[#, InputForm] <> "}") & /@ integrals,
      ",\n"
    ],
    "}\n"
  ];
  intFile = name <> "_reduce.m";
  Export[FileNameJoin[{dir, intFile}], intContent, "Text"];

  (* Write config *)
  writeSingleProblemConfig[topology, dir, name <> "_reduce.config",
                           intFile, name <> "_reduce.tables"];

  (* Run FIRE6 *)
  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 2,
    Print["Reducing ", Length[integrals], " integrals..."];
  ];
  runFIRE6[fireBin, dir, name <> "_reduce"];

  (* Load results *)
  tablesFile = FileNameJoin[{dir, name <> "_reduce.tables"}];
  If[!FileExistsQ[tablesFile],
    Print["Error: FIRE6 reduction did not produce output."];
    Return[$Failed];
  ];

  rules = FIRE`Tables2Rules[tablesFile];

  (* Convert to our format *)
  (* FIRE uses G[problemNumber, indices] *)
  result = Association[
    (Global`G[pn, #] /. rules) & /@ integrals //
    (* Replace G[pn, ...] with G[1, ...] for consistency *)
    (# /. Global`G[pn, idx_] :> Global`G[1, idx]) & //
    MapThread[Rule, {integrals, #}] &
  ];

  result
];


(* ============================================================ *)
(* ReduceIntegralsBatch - reduce integrals for multiple topos   *)
(* ============================================================ *)

ReduceIntegralsBatch[topoIntegralPairs_List] :=
Module[{dir, fireBin, intContent, tablesFile, rules, results, allIntegrals},

  If[Length[topoIntegralPairs] == 0, Return[{}]];

  (* Use first topology's work directory *)
  dir = topoIntegralPairs[[1, 1]]["WorkDirectory"];
  fireBin = FileNameJoin[{FeynmanTrick`Private`$FTConfig["FIREPath"], "bin", "FIRE6"}];

  (* Build combined integrals file *)
  allIntegrals = Flatten[
    Table[
      With[{topo = topoIntegralPairs[[i, 1]],
            ints = topoIntegralPairs[[i, 2]],
            pn = topoIntegralPairs[[i, 1]]["ProblemNumber"]},
        {pn, #} & /@ ints
      ],
      {i, Length[topoIntegralPairs]}
    ],
    1
  ];

  intContent = StringJoin[
    "{",
    StringRiffle[
      ("{" <> ToString[#[[1]]] <> "," <> ToString[#[[2]], InputForm] <> "}") & /@ allIntegrals,
      ",\n"
    ],
    "}\n"
  ];
  Export[FileNameJoin[{dir, "batch_reduce.m"}], intContent, "Text"];

  (* Write multi-problem config *)
  writeMultiProblemConfig[topoIntegralPairs[[All, 1]], dir, "batch_reduce.config",
                          "batch_reduce.m", "batch_reduce.tables"];

  (* Run FIRE6 *)
  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["Reducing integrals for ", Length[topoIntegralPairs], " topologies..."];
  ];
  runFIRE6[fireBin, dir, "batch_reduce"];

  (* Load results *)
  tablesFile = FileNameJoin[{dir, "batch_reduce.tables"}];
  If[!FileExistsQ[tablesFile],
    Print["Error: FIRE6 reduction did not produce output."];
    Return[$Failed];
  ];

  rules = FIRE`Tables2Rules[tablesFile];

  (* Parse results for each topology *)
  results = Table[
    With[{topo = topoIntegralPairs[[i, 1]],
          ints = topoIntegralPairs[[i, 2]],
          pn = topoIntegralPairs[[i, 1]]["ProblemNumber"]},
      Association[
        (Global`G[pn, #] /. rules) & /@ ints //
        (* Normalize G[pn, ...] -> G[1, ...] *)
        (# /. Global`G[pn, idx_] :> Global`G[1, idx]) & //
        MapThread[Rule, {ints, #}] &
      ]
    ],
    {i, Length[topoIntegralPairs]}
  ];

  results
];


(* ============================================================ *)
(* ComputeDiffMatrix                                             *)
(* ============================================================ *)

ComputeDiffMatrix[topology_Association, variable_, precomputedDecomp_:{Automatic}] :=
Module[{masters, nMasters, coeffMat, constVec, shiftedPerMaster,
        uniqueIntegrals, reductions, diffMatrix, i, j, expr},

  masters = topology["Masters"];
  If[masters === {},
    Print["Error: No masters found. Call FindBasis first."];
    Return[$Failed];
  ];
  nMasters = Length[masters];

  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["Computing differential matrix in ", variable,
          " for ", nMasters, " masters..."];
  ];

  (* Step 1: Decompose propagator derivatives *)
  If[precomputedDecomp =!= {Automatic} && Length[precomputedDecomp] == 2,
    {coeffMat, constVec} = precomputedDecomp;
  ,
    {coeffMat, constVec} = DecomposePropagatorDerivative[
      topology["Propagators"],
      topology["LoopMomenta"],
      variable,
      topology["Replacements"]
    ];
  ];

  (* Step 2: For each master, find shifted integrals in d/dx M_i *)
  shiftedPerMaster = Table[
    DifferentiatedIntegrals[masters[[i]], coeffMat, constVec],
    {i, nMasters}
  ];

  (* Step 3: Collect all unique integrals needing reduction *)
  uniqueIntegrals = DeleteDuplicates[
    Flatten[shiftedPerMaster[[All, All, 1]], 1]
  ];

  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["Need to reduce ", Length[uniqueIntegrals], " integrals."];
  ];

  (* Step 4: Reduce all integrals *)
  reductions = ReduceIntegrals[topology, uniqueIntegrals];

  If[reductions === $Failed,
    Print["Error: FIRE reduction failed."];
    Return[$Failed];
  ];

  (* Map FIRE's Global variables back to the requested variable symbol *)
  If[Context[variable] =!= "Global`",
    Module[{globalVar = Symbol[SymbolName[variable]]},
      reductions = Map[(# /. globalVar -> variable) &, reductions];
    ];
  ];

  (* Step 5: Assemble the matrix *)
  diffMatrix = Table[0, {nMasters}, {nMasters}];

  Do[
    expr = Total[
      (#[[2]] * (reductions[#[[1]]])) & /@ shiftedPerMaster[[i]]
    ];
    Do[
      diffMatrix[[i, j]] = Coefficient[expr, Global`G[1, masters[[j]]]];
    , {j, nMasters}];
  , {i, nMasters}];

  diffMatrix = Map[Together, diffMatrix, {2}];

  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["Differential matrix computed (", nMasters, "x", nMasters, ")."];
  ];

  diffMatrix
];


(* ============================================================ *)
(* GetFIREResult                                                 *)
(* ============================================================ *)

GetFIREResult[topology_Association, integral_List] :=
Module[{result},
  result = ReduceIntegrals[topology, {integral}];
  If[result === $Failed, Return[$Failed]];
  result[integral]
];


End[];
EndPackage[];
