(* ::Package:: *)
(* FIREInterface - Wrapper for FIRE6 IBP reduction *)
(* Handles topology definition, FIRE setup, basis finding, *)
(* integral reduction, and differential matrix computation *)

BeginPackage["FeynmanTrick`FIREInterface`", {"FeynmanTrick`", "FeynmanTrick`PropagatorAlgebra`"}];

DefineTopology::usage =
  "DefineTopology[name, loopMomenta, externalMomenta, propagators, replacements] \
defines a topology for FIRE reduction. Propagator convention: D_j = -q_j^2 + m_j^2. \
Returns a Topology association.";

SetupFIRE::usage =
  "SetupFIRE[topology, workDir] generates .start and .config files for FIRE6. \
Loads FIRE6.m, sets up IBP relations, and saves the start file.";

FindBasis::usage =
  "FindBasis[topology] runs FIRE6 to determine master integrals. \
Returns the list of master integral index vectors.";

ReduceIntegrals::usage =
  "ReduceIntegrals[topology, integrals] reduces a list of integrals \
(given as index vectors) to the master basis using FIRE6. \
Returns an association: integral -> {coeff, masterIndex} pairs.";

ComputeDiffMatrix::usage =
  "ComputeDiffMatrix[topology, variable] computes the differential equation \
matrix A such that d/dx Masters = A . Masters. \
The variable is typically the Feynman parameter xx. \
Optional: ComputeDiffMatrix[topology, variable, {coeffMatrix, constVector}] \
uses a pre-computed decomposition (e.g. from FeynmanTrickDecomposition).";

GetFIREResult::usage =
  "GetFIREResult[topology, integral] returns the reduction of a single integral.";

TopologyQ::usage = "TopologyQ[expr] returns True if expr is a valid topology association.";

Begin["`Private`"];

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
    "Masters" -> {},
    "MasterRules" -> {},
    "StartFileReady" -> False,
    "Variables" -> extractVariables[propagators, replacements, loopMomenta, externalMomenta]
  |>;

TopologyQ[t_Association] := KeyExistsQ[t, "Propagators"] && KeyExistsQ[t, "LoopMomenta"];
TopologyQ[_] := False;


(* Extract symbolic variables from propagators and replacements *)
extractVariables[props_List, repls_List, loops_List, exts_List] :=
Module[{allSyms, momenta, vars},
  allSyms = Cases[{props, repls}, _Symbol, Infinity] // DeleteDuplicates;
  momenta = Join[loops, exts];
  vars = Complement[allSyms, momenta];
  vars
];


(* ============================================================ *)
(* SetupFIRE                                                     *)
(* Generates .start file via FIRE6.m API                         *)
(* ============================================================ *)

SetupFIRE[topology_Association, workDir_String:""] :=
Module[{dir, firePath, name, config, result},
  name = topology["Name"];
  firePath = FeynmanTrick`Private`$FTConfig["FIREPath"];

  (* Determine working directory *)
  dir = If[workDir === "",
    FileNameJoin[{FeynmanTrick`Private`$FTConfig["WorkDirectory"], name}],
    workDir
  ];
  If[!DirectoryQ[dir], CreateDirectory[dir, CreateIntermediateDirectories -> True]];

  (* Load FIRE6.m if not already loaded *)
  If[!MemberQ[$Packages, "FIRE`"],
    Block[{$Path = Prepend[$Path, firePath]},
      Quiet[Get["FIRE6.m"], {General::shdw}];
    ];
  ];

  (* Set FIRE variables (all in FIRE` context) *)
  FIRE`Internal = topology["LoopMomenta"];
  FIRE`External = topology["ExternalMomenta"];
  FIRE`Propagators = topology["Propagators"];

  (* Replacements: convert our rule format to FIRE's *)
  If[topology["Replacements"] =!= {},
    FIRE`Replacements = topology["Replacements"];
  ];

  (* Prepare IBP relations *)
  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["Setting up FIRE for topology: ", name];
  ];
  FIRE`PrepareIBP[];

  (* Prepare sector basis *)
  FIRE`Prepare[AutoDetectRestrictions -> True];

  (* Save start file *)
  Module[{oldDir = Directory[]},
    SetDirectory[dir];
    FIRE`SaveStart[name];
    SetDirectory[oldDir];
  ];

  (* Write config file *)
  writeConfigFile[topology, dir, name <> ".config"];

  (* Return updated topology *)
  result = topology;
  result["WorkDirectory"] = dir;
  result["StartFileReady"] = True;
  result
];


(* Write FIRE6 config file *)
writeConfigFile[topology_Association, dir_String, configName_String,
                integralsFile_String:"", outputFile_String:""] :=
Module[{vars, content, threads, fthreads, intFile, outFile, name},
  name = topology["Name"];
  threads = FeynmanTrick`Private`$FTConfig["Threads"];
  fthreads = FeynmanTrick`Private`$FTConfig["FThreads"];

  (* Variables: d plus any kinematic/Feynman parameter variables *)
  vars = Join[
    {FeynmanTrick`Private`$FTConfig["DimensionVariable"]},
    topology["Variables"]
  ] // DeleteDuplicates;

  intFile = If[integralsFile === "", name <> ".m", integralsFile];
  outFile = If[outputFile === "", name <> ".tables", outputFile];

  content = StringJoin[
    "#threads           ", ToString[threads], "\n",
    "#fthreads          ", ToString[fthreads], "\n",
    "#variables         ", StringRiffle[ToString /@ vars, ","], "\n",
    "#start\n",
    "#problem           1 ", name, ".start\n",
    "#integrals         ", intFile, "\n",
    "#output            ", outFile, "\n"
  ];

  Export[FileNameJoin[{dir, configName}], content, "Text"];
];


(* ============================================================ *)
(* FindBasis                                                     *)
(* Runs FIRE6 to determine master integrals                      *)
(* ============================================================ *)

FindBasis[topology_Association] :=
Module[{dir, name, fireBin, intFile, configFile, result, masters,
        topSector, tablesFile, rules, allIntegrals},

  If[!topology["StartFileReady"],
    Print["Error: Must call SetupFIRE before FindBasis."];
    Return[$Failed];
  ];

  dir = topology["WorkDirectory"];
  name = topology["Name"];
  fireBin = FileNameJoin[{FeynmanTrick`Private`$FTConfig["FIREPath"], "bin", "FIRE6"}];

  (* Write integrals file with the top-sector integral *)
  topSector = Table[1, {topology["NumPropagators"]}];
  intFile = FileNameJoin[{dir, name <> ".m"}];
  Export[intFile, "{{1," <> ToString[topSector, InputForm] <> "}}\n", "Text"];

  (* Write config for this reduction *)
  writeConfigFile[topology, dir, name <> ".config"];

  (* Run FIRE6 *)
  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["Running FIRE6 to find basis..."];
  ];
  Module[{cmd, exitCode},
    cmd = "cd " <> StringReplace[dir, " " -> "\\ "] <>
          " && " <> StringReplace[fireBin, " " -> "\\ "] <> " -c " <> name <>
          " > fire_stdout.log 2> fire_stderr.log";
    exitCode = Run[cmd];
    If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 2,
      Print["FIRE6 exit code: ", exitCode];
    ];
    If[exitCode =!= 0,
      Print["FIRE6 failed with exit code ", exitCode];
      If[FileExistsQ[FileNameJoin[{dir, "fire_stderr.log"}]],
        Print["stderr: ", Import[FileNameJoin[{dir, "fire_stderr.log"}], "Text"]];
      ];
    ];
  ];

  (* Load tables and extract masters *)
  tablesFile = FileNameJoin[{dir, name <> ".tables"}];
  If[!FileExistsQ[tablesFile],
    Print["Error: FIRE6 did not produce output tables."];
    Print["Binary: ", fireBin, " exists: ", FileExistsQ[fireBin]];
    Return[$Failed];
  ];

  (* Extract master integrals using FIRE's Tables2Masters *)
  masters = FIRE`Tables2Masters[tablesFile];
  (* Tables2Masters returns {{problemNum, {indices}}, ...} - extract just the indices *)
  masters = Cases[masters, {_, indices_List} :> indices];

  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["Found ", Length[masters], " master integrals."];
  ];

  (* Return updated topology *)
  result = topology;
  result["Masters"] = masters;
  result
];


(* Extract master integral indices from reduction rules *)
extractMastersFromRules[rules_List, nProps_Integer] :=
Module[{allG, masters},
  (* rules are of the form: G[1, {indices}] -> expression *)
  (* Masters are the G[1, {...}] that appear on the RHS *)
  allG = Cases[rules[[All, 2]], Global`G[1, _List], Infinity] // DeleteDuplicates;
  masters = Cases[allG, Global`G[1, indices_List] :> indices];
  (* Sort by total index sum (simpler integrals first) *)
  SortBy[masters, Total]
];


(* ============================================================ *)
(* ReduceIntegrals                                               *)
(* Reduces a list of integrals to the master basis               *)
(* ============================================================ *)

ReduceIntegrals[topology_Association, integrals_List] :=
Module[{dir, name, fireBin, intFile, configFile, tablesFile, rules,
        intContent, result},

  If[!topology["StartFileReady"],
    Print["Error: Must call SetupFIRE before ReduceIntegrals."];
    Return[$Failed];
  ];

  dir = topology["WorkDirectory"];
  name = topology["Name"];
  fireBin = FileNameJoin[{FeynmanTrick`Private`$FTConfig["FIREPath"], "bin", "FIRE6"}];

  (* Write integrals file *)
  intContent = StringJoin[
    "{",
    StringRiffle[
      ("{1," <> ToString[#, InputForm] <> "}") & /@ integrals,
      ",\n"
    ],
    "}\n"
  ];
  intFile = name <> "_reduce.m";
  Export[FileNameJoin[{dir, intFile}], intContent, "Text"];

  (* Write config *)
  writeConfigFile[topology, dir, name <> "_reduce.config",
                  intFile, name <> "_reduce.tables"];

  (* Run FIRE6 *)
  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 2,
    Print["Reducing ", Length[integrals], " integrals..."];
  ];
  Module[{cmd, exitCode},
    cmd = "cd " <> StringReplace[dir, " " -> "\\ "] <>
          " && " <> StringReplace[fireBin, " " -> "\\ "] <>
          " -c " <> name <> "_reduce" <>
          " > fire_stdout.log 2> fire_stderr.log";
    exitCode = Run[cmd];
    If[exitCode =!= 0,
      Print["FIRE6 reduction failed with exit code ", exitCode];
    ];
  ];

  (* Load results *)
  tablesFile = FileNameJoin[{dir, name <> "_reduce.tables"}];
  If[!FileExistsQ[tablesFile],
    Print["Error: FIRE6 reduction did not produce output."];
    Return[$Failed];
  ];

  rules = FIRE`Tables2Rules[tablesFile];

  (* Convert to our format: integral index -> linear combination of masters *)
  result = Association[
    (Global`G[1, #] /. rules) & /@ integrals //
    MapThread[Rule, {integrals, #}] &
  ];

  result
];


(* ============================================================ *)
(* ComputeDiffMatrix                                             *)
(* Computes the differential equation matrix                     *)
(* d/dx vec{M} = A(x) . vec{M}                                  *)
(* ============================================================ *)

ComputeDiffMatrix[topology_Association, variable_, precomputedDecomp_:{Automatic}] :=
Module[{masters, nMasters, coeffMat, constVec, allShifted, shiftedPerMaster,
        uniqueIntegrals, reductions, diffMatrix, i, j, expr, updatedTopology},

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
    (* Use pre-computed decomposition (e.g. from FeynmanTrickDecomposition) *)
    {coeffMat, constVec} = precomputedDecomp;
  ,
    {coeffMat, constVec} = DecomposePropagatorDerivative[
      topology["Propagators"],
      topology["LoopMomenta"],
      variable,
      topology["Replacements"]
    ];
  ];

  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 2,
    Print["Derivative decomposition computed."];
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

  (* Remove integrals that are already masters (no reduction needed) *)
  (* But keep them in the list for uniform processing *)

  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["Need to reduce ", Length[uniqueIntegrals], " integrals."];
  ];

  (* Step 4: Reduce all integrals *)
  reductions = ReduceIntegrals[topology, uniqueIntegrals];

  If[reductions === $Failed,
    Print["Error: FIRE reduction failed."];
    Return[$Failed];
  ];

  (* Step 5: Assemble the matrix *)
  diffMatrix = Table[0, {nMasters}, {nMasters}];

  Do[
    (* d/dx M_i = Sum over shifted integrals *)
    expr = Total[
      (#[[2]] * (reductions[#[[1]]])) & /@ shiftedPerMaster[[i]]
    ];

    (* Extract coefficient of each master G[1, masters[[j]]] *)
    Do[
      diffMatrix[[i, j]] = Coefficient[expr, Global`G[1, masters[[j]]]];
    , {j, nMasters}];
  , {i, nMasters}];

  (* Simplify entries *)
  diffMatrix = Together /@ diffMatrix // Map[#, {2}] &;

  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["Differential matrix computed (", nMasters, "x", nMasters, ")."];
  ];

  diffMatrix
];


(* ============================================================ *)
(* GetFIREResult                                                 *)
(* Returns reduction of a single integral (for debugging)        *)
(* ============================================================ *)

GetFIREResult[topology_Association, integral_List] :=
Module[{result},
  result = ReduceIntegrals[topology, {integral}];
  If[result === $Failed, Return[$Failed]];
  result[integral]
];


End[];
EndPackage[];
