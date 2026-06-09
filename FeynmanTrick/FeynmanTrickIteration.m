(* ::Package:: *)
(* FeynmanTrickIteration - Multi-level Feynman trick orchestrator *)
(* Manages the iterative propagator combination and matrix computation *)

BeginPackage["FeynmanTrick`FeynmanTrickIteration`", {
  "FeynmanTrick`",
  "FeynmanTrick`PropagatorAlgebra`",
  "FeynmanTrick`FIREInterface`",
  "FeynmanTrick`MatrixExport`",
  "FeynmanTrick`EpsPrefactors`"
}];

DefineFTIteration::usage =
  "DefineFTIteration[topology, combinationSequence, numericalPoint] defines a \
multi-level Feynman trick iteration. \
combinationSequence: list of {i,j} pairs (propagator positions to combine). \
numericalPoint: replacement rules for kinematic variables.";

BuildLevel::usage =
  "BuildLevel[ftData, level] constructs the topology at a given level by \
combining the specified propagators. Returns updated ftData.";

ComputeLevelData::usage =
  "ComputeLevelData[ftData, level] runs the full computation for one level: \
SetupFIRE + FindBasis + ComputeDiffMatrix. Returns updated ftData.";

ExportLevel::usage =
  "ExportLevel[ftData, level, directory, format, epsOrder] exports level matrices. \
format: \"diffexp\" (order-by-order) or \"general\" (full eps).";

IdentifyNeededIntegrals::usage =
  "IdentifyNeededIntegrals[ftData, level] identifies which integrals from level+1 \
are needed for the Feynman trick recursion formula at the given level.";

RunFullIteration::usage =
  "RunFullIteration[ftData, outputDir] runs all levels bottom-up, computes \
matrices, and exports them. Returns the fully populated ftData.";

BuildAllLevels::usage =
  "BuildAllLevels[ftData] builds all levels of the iteration without computing. \
Useful for inspecting the propagator structure.";

CloseLevelMasterBasis::usage =
  "CloseLevelMasterBasis[topology, parameter, combinedPositions, seedMasters] \
builds a FIRE master basis closed under the Feynman trick seed integrals and \
their parameter derivatives.";

GetLevelTopology::usage =
  "GetLevelTopology[ftData, level] returns the topology association for a level.";

GetLevelMatrix::usage =
  "GetLevelMatrix[ftData, level] returns the differential matrix for a level.";

Begin["`Private`"];

(* ============================================================ *)
(* DefineFTIteration                                             *)
(* Sets up the multi-level iteration data structure              *)
(* ============================================================ *)

DefineFTIteration[topology_Association, combinationSeq_List, numericalPoint_List:{}] :=
Module[{ftData, nLevels, eliminatedPositions, topMaster},
  nLevels = Length[combinationSeq];
  eliminatedPositions = Lookup[topology, "EliminatedPositions", {}];
  topMaster = ReplacePart[
    ConstantArray[1, topology["NumPropagators"]],
    Thread[eliminatedPositions -> 0]
  ];

  ftData = <|
    "TopTopology" -> topology,
    "CombinationSequence" -> combinationSeq,
    "NumericalPoint" -> numericalPoint,
    "NumLevels" -> nLevels,
    "FixedParamValue" -> FeynmanTrick`Private`$FTConfig["FixedParameterValue"],
    "Levels" -> <|
      0 -> <|
        "Topology" -> topology,
        "Propagators" -> topology["Propagators"],
        "FeynmanParameter" -> None,
        "FixedParams" -> {},
        "EliminatedPositions" -> eliminatedPositions,
        "Masters" -> {topMaster},
        "DiffMatrix" -> {},
        "Computed" -> False
      |>
    |>
  |>;

  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["Defined Feynman trick iteration with ", nLevels, " levels."];
    Print["Top topology has ", topology["NumPropagators"], " propagators."];
  ];

  ftData
];


(* ============================================================ *)
(* BuildLevel                                                    *)
(* Constructs the topology at level k by combining propagators   *)
(*                                                               *)
(* At level k (with combination {i, j}):                         *)
(* - Position i gets replaced by xx*D_i + (1-xx)*D_j            *)
(* - Position j stays unchanged (sub-sector)                     *)
(* - Previous Feynman parameters fixed to 11/23                  *)
(* - Kinematics fixed to numerical point                         *)
(* ============================================================ *)

BuildLevel[ftData_Association, level_Integer] :=
Module[{prevLevel, combo, prevProps, newProps, newProp, i, j,
        param, fixedVal, prevParam, numericalPoint, newTopology,
        prevFixedParams, fixRules, name, result, newReplacements,
        prevEliminated, eliminatedPositions},

  If[level < 1 || level > ftData["NumLevels"],
    Print["Error: level must be between 1 and ", ftData["NumLevels"]];
    Return[ftData];
  ];

  (* Ensure previous level exists *)
  If[!KeyExistsQ[ftData["Levels"], level - 1],
    Print["Error: level ", level - 1, " not yet built."];
    Return[ftData];
  ];

  prevLevel = ftData["Levels"][level - 1];
  combo = ftData["CombinationSequence"][[level]];
  {i, j} = combo;
  prevProps = prevLevel["Propagators"];
  prevEliminated = Lookup[prevLevel, "EliminatedPositions", {}];
  eliminatedPositions = Sort[DeleteDuplicates[Append[prevEliminated, j]]];

  (* Each level gets a unique Feynman parameter symbol: xx1, xx2, xx3, ... *)
  (* This is essential: the paper uses x_1, x_2, ..., x_{n-1} for each combination step *)
  param = Symbol["xx" <> ToString[level]];
  fixedVal = ftData["FixedParamValue"];
  numericalPoint = ftData["NumericalPoint"];

  (* Build the combined propagator *)
  newProp = CombinePropagatorSymbolic[prevProps[[i]], prevProps[[j]], param];

  (* New propagator list: replace position i with combined, keep j unchanged *)
  newProps = prevProps;
  newProps[[i]] = newProp;

  (* Fix previous level's Feynman parameter to numerical value *)
  prevFixedParams = prevLevel["FixedParams"];
  If[prevLevel["FeynmanParameter"] =!= None,
    fixRules = {prevLevel["FeynmanParameter"] -> fixedVal};
    newProps = newProps /. fixRules;
    prevFixedParams = Append[prevFixedParams,
      prevLevel["FeynmanParameter"] -> fixedVal];
  ];

  (* Propagator replacements: original scalar product rules with numerical kinematic values *)
  newReplacements = ftData["TopTopology"]["Replacements"] /. numericalPoint;

  (* Create new topology *)
  name = ftData["TopTopology"]["Name"] <> "_L" <> ToString[level];
  newTopology = DefineTopology[
    name,
    ftData["TopTopology"]["LoopMomenta"],
    ftData["TopTopology"]["ExternalMomenta"],
    newProps,
    newReplacements
  ];
  newTopology["EliminatedPositions"] = eliminatedPositions;

  (* Update ftData *)
  result = ftData;
  result["Levels"][level] = <|
    "Topology" -> newTopology,
    "Propagators" -> newProps,
    "FeynmanParameter" -> param,
    "FixedParams" -> prevFixedParams,
    "CombinedPositions" -> {i, j},
    "EliminatedPositions" -> eliminatedPositions,
    "Masters" -> {},
    "DiffMatrix" -> {},
    "Computed" -> False
  |>;

  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["Built level ", level, ": combined positions {", i, ",", j, "}"];
  ];

  result
];


(* ============================================================ *)
(* BuildAllLevels                                                *)
(* Convenience: builds all levels                                *)
(* ============================================================ *)

BuildAllLevels[ftData_Association] :=
Module[{result = ftData, k},
  Do[
    result = BuildLevel[result, k];
  , {k, 1, ftData["NumLevels"]}];
  result
];


(* ============================================================ *)
(* Master-basis closure for Feynman trick levels                 *)
(* ============================================================ *)

FeynmanTrickNeededIntegral[masterVec_List, combinedPositions_List] :=
Module[{posI, posJ, vi, vj, neededVec, case},
  {posI, posJ} = combinedPositions;
  vi = masterVec[[posI]];
  vj = masterVec[[posJ]];

  case = Which[
    vi > 0 && vj > 0, "integrate",
    vi > 0 && vj == 0, "limitUpper",
    vi == 0 && vj > 0, "limitLower",
    True, "direct"
  ];

  neededVec = masterVec;
  Switch[case,
    "integrate",
      neededVec[[posI]] = vi + vj;
      neededVec[[posJ]] = 0;,
    "limitUpper",
      neededVec[[posI]] = vi;
      neededVec[[posJ]] = 0;,
    "limitLower",
      neededVec[[posI]] = vj;
      neededVec[[posJ]] = 0;,
    "direct",
      neededVec[[posJ]] = 0
  ];

  neededVec
];

FeynmanTrickSeedMasters[mastersAtPreviousLevel_List,
    combinedPositions_List] :=
  DeleteDuplicates[
    FeynmanTrickNeededIntegral[#, combinedPositions] & /@
      mastersAtPreviousLevel
  ];

derivativeShiftedIntegrals[masters_List, decomp_] :=
Module[{shifted},
  shifted = Flatten[
    (FeynmanTrick`PropagatorAlgebra`DifferentiatedIntegrals[
        #, decomp[[1]], decomp[[2]]
      ][[All, 1]] &) /@ masters,
    1
  ];
  DeleteDuplicates[shifted]
];

CloseLevelMasterBasis[topology_Association, parameter_,
    combinedPositions_List, seedMasters_List] :=
Module[{setupTopo, masters = {}, detailed, decomp, shiftedIntegrals,
        newMasters, iter = 0, maxIter = 12, resultTopo, paddedSeedMasters},

  If[seedMasters === {},
    Print["Error: cannot close Feynman trick master basis without seed masters."];
    Return[$Failed];
  ];

  setupTopo = FeynmanTrick`FIREInterface`SetupFIRE[topology];
  If[setupTopo === $Failed, Return[$Failed]];

  paddedSeedMasters = PadRight[#, setupTopo["NumPropagators"], 0] & /@
    DeleteDuplicates[seedMasters];

  detailed = FeynmanTrick`FIREInterface`ReduceIntegralsDetailed[
    setupTopo, paddedSeedMasters
  ];
  If[detailed === $Failed, Return[$Failed]];

  masters = DeleteDuplicates[detailed["Masters"]];
  If[masters === {},
    masters = paddedSeedMasters;
  ];

  decomp = FeynmanTrick`PropagatorAlgebra`FeynmanTrickDecomposition[
    setupTopo["NumPropagators"], combinedPositions[[1]],
    combinedPositions[[2]], parameter
  ];

  While[iter < maxIter,
    shiftedIntegrals = derivativeShiftedIntegrals[masters, decomp];
    If[shiftedIntegrals === {}, Break[]];

    detailed = FeynmanTrick`FIREInterface`ReduceIntegralsDetailed[
      setupTopo, shiftedIntegrals
    ];
    If[detailed === $Failed, Return[$Failed]];

    newMasters = Complement[detailed["Masters"], masters];
    If[newMasters === {}, Break[]];

    masters = DeleteDuplicates[Join[masters, newMasters]];
    iter++;
  ];

  If[iter >= maxIter,
    Print["Warning: Feynman trick master-basis closure reached iteration limit for ",
      topology["Name"]];
  ];

  resultTopo = setupTopo;
  resultTopo["Masters"] = masters;

  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["Closed Feynman trick basis for ", topology["Name"], ": ",
      Length[masters], " masters"];
    If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 2,
      Print["  Seed masters: ", seedMasters];
      Print["  Closed masters: ", masters];
    ];
  ];

  <|"Topology" -> resultTopo, "Masters" -> masters|>
];


(* ============================================================ *)
(* ComputeLevelData                                              *)
(* Full computation for one level:                               *)
(* SetupFIRE -> FindBasis -> ComputeDiffMatrix                   *)
(* ============================================================ *)

ComputeLevelData[ftData_Association, level_Integer] :=
Module[{levelData, topology, masters, diffMat, param, result, updatedTopo,
        decomp, combo, combinedPos, otherPos, seedMasters, closedBasis},

  If[!KeyExistsQ[ftData["Levels"], level],
    Print["Error: level ", level, " not built. Call BuildLevel first."];
    Return[ftData];
  ];

  levelData = ftData["Levels"][level];
  topology = levelData["Topology"];
  param = levelData["FeynmanParameter"];
  seedMasters = Lookup[levelData, "SeedMasters", Automatic];

  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["\n=== Computing level ", level, " ==="];
  ];

  combo = levelData["CombinedPositions"];
  If[ListQ[seedMasters] && seedMasters =!= {} && ListQ[combo],
    closedBasis = CloseLevelMasterBasis[topology, param, combo, seedMasters];
    If[closedBasis === $Failed,
      Print["Error: CloseLevelMasterBasis failed for level ", level];
      result = ftData;
      result["Levels"][level]["Computed"] = False;
      result["Levels"][level]["Error"] = "CloseLevelMasterBasis failed";
      Return[result];
    ];
    updatedTopo = closedBasis["Topology"];
    masters = closedBasis["Masters"];
  ,
    (* Step 1: Setup FIRE *)
    updatedTopo = SetupFIRE[topology];
    If[updatedTopo === $Failed,
      Print["Error: SetupFIRE failed for level ", level];
      result = ftData;
      result["Levels"][level]["Computed"] = False;
      result["Levels"][level]["Error"] = "SetupFIRE failed";
      Return[result];
    ];

    (* Step 2: Find master basis *)
    updatedTopo = FindBasis[updatedTopo];
    If[updatedTopo === $Failed,
      Print["Error: FindBasis failed for level ", level];
      result = ftData;
      result["Levels"][level]["Computed"] = False;
      result["Levels"][level]["Error"] = "FindBasis failed";
      Return[result];
    ];
    masters = updatedTopo["Masters"];
  ];

  If[!ListQ[masters] || Length[masters] == 0,
    Print["Error: No masters found for level ", level];
    result = ftData;
    result["Levels"][level]["Computed"] = False;
    result["Levels"][level]["Error"] = "No masters found";
    Return[result];
  ];

  (* Step 3: Compute differential matrix *)
  (* Use the Feynman trick fast path for the decomposition *)
  If[combo =!= Missing["KeyAbsent", "CombinedPositions"] && ListQ[combo],
    {combinedPos, otherPos} = combo;
    decomp = FeynmanTrickDecomposition[
      updatedTopo["NumPropagators"], combinedPos, otherPos, param
    ];
    diffMat = ComputeDiffMatrix[updatedTopo, param, decomp];
  ,
    (* General case: no combination info, use full decomposition *)
    diffMat = ComputeDiffMatrix[updatedTopo, param];
  ];

  If[diffMat === $Failed,
    Print["Error: ComputeDiffMatrix failed for level ", level];
    result = ftData;
    result["Levels"][level]["Computed"] = False;
    result["Levels"][level]["Error"] = "ComputeDiffMatrix failed";
    Return[result];
  ];

  (* Update ftData *)
  result = ftData;
  result["Levels"][level]["Topology"] = updatedTopo;
  result["Levels"][level]["Masters"] = masters;
  result["Levels"][level]["DiffMatrix"] = diffMat;
  result["Levels"][level]["Computed"] = True;

  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["Level ", level, " complete: ", Length[masters], " masters, ",
          Length[masters], "x", Length[masters], " matrix."];
  ];

  result
];


(* ============================================================ *)
(* ExportLevel                                                   *)
(* Export matrices for a level                                   *)
(* ============================================================ *)

ExportLevel[ftData_Association, level_Integer, directory_String,
            format_String:"diffexp", epsOrder_Integer:4] :=
Module[{levelData, diffMat, param, levelDir, prefactors, matToExport},

  If[!KeyExistsQ[ftData["Levels"], level] ||
     !ftData["Levels"][level]["Computed"],
    Print["Error: level ", level, " not computed yet."];
    Return[];
  ];

  levelData = ftData["Levels"][level];
  diffMat = levelData["DiffMatrix"];
  param = levelData["FeynmanParameter"];

  (* Create level-specific subdirectory *)
  levelDir = FileNameJoin[{directory, "Level_" <> ToString[level] <> "_Matrices"}];

  (* Check for eps poles and apply prefactors if needed *)
  Module[{eps = FeynmanTrick`Private`$FTConfig["EpsilonSymbol"]},
    If[CheckEpsPoles[diffMat, eps],
      prefactors = FindEpsPrefactors[diffMat, eps];
      matToExport = ApplyEpsPrefactors[diffMat, prefactors, eps];
      (* Also export prefactors for reference *)
      If[!DirectoryQ[levelDir], CreateDirectory[levelDir, CreateIntermediateDirectories -> True]];
      Put[prefactors, FileNameJoin[{levelDir, "eps_prefactors.m"}]];
      If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
        Print["Applied eps prefactors to remove poles."];
      ];
    ,
      matToExport = diffMat;
    ];
  ]; (* End Module *)

  Switch[format,
    "diffexp",
      ExportDiffExpMatrix[matToExport, param, levelDir, epsOrder],
    "general",
      ExportGeneralMatrix[matToExport, param, levelDir],
    "both",
      ExportDiffExpMatrix[matToExport, param, levelDir, epsOrder];
      ExportGeneralMatrix[matToExport, param, levelDir],
    _,
      Print["Unknown format: ", format, ". Use \"diffexp\", \"general\", or \"both\"."]
  ];
];


(* ============================================================ *)
(* IdentifyNeededIntegrals                                       *)
(* Identifies integrals from level+1 needed for the Feynman      *)
(* trick recursion at the given level                            *)
(*                                                               *)
(* At level l, the recursion involves:                           *)
(* I^(l-1)_{v1,...,v_{n-l+1}} =                                 *)
(*   Gamma(v1+v2)/(Gamma(v1)*Gamma(v2)) *                       *)
(*   Int_0^1 dx x^{v1-1} (1-x)^{v2-1} *                        *)
(*   I^(l)_{v1+v2, v3, ..., v_{n-l}}                            *)
(*                                                               *)
(* So we need integrals from level l where the combined          *)
(* propagator has index v1+v2 and other indices are v3,...        *)
(* ============================================================ *)

IdentifyNeededIntegrals[ftData_Association, level_Integer] :=
Module[{nextLevel, combo, nextMasters, neededIntegrals, i, j, nProps},

  If[level >= ftData["NumLevels"],
    Print["Error: no level ", level + 1, " exists."];
    Return[{}];
  ];

  (* The combination at level+1 tells us which propagators were merged *)
  combo = ftData["CombinationSequence"][[level + 1]];
  {i, j} = combo;

  (* Masters at level+1 *)
  If[!KeyExistsQ[ftData["Levels"], level + 1] ||
     ftData["Levels"][level + 1]["Masters"] === {},
    Print["Warning: level ", level + 1, " masters not computed yet."];
    Return[{}];
  ];

  nextMasters = ftData["Levels"][level + 1]["Masters"];
  nProps = Length[nextMasters[[1]]];

  (* For each master at level+1, determine what integrals at level *)
  (* are involved in the Feynman trick formula *)
  (* The combined index at position i becomes separate indices v_i, v_j *)
  (* at level, where v_i + v_j = combined index *)
  neededIntegrals = Flatten[Table[
    Module[{v = master, combinedIdx, possibleSplits},
      combinedIdx = v[[i]]; (* This was v_1 + v_2 in the recursion *)

      (* All ways to split combinedIdx into v_i + v_j with v_i, v_j >= 1 *)
      possibleSplits = Table[
        Module[{newV = v},
          newV[[i]] = vi;
          (* Insert v_j back at position j *)
          (* At level l, position j has the uncombined propagator *)
          newV[[j]] = combinedIdx - vi;
          newV
        ],
        {vi, 1, combinedIdx - 1}
      ];

      possibleSplits
    ],
    {master, nextMasters}
  ], 1] // DeleteDuplicates;

  neededIntegrals
];


(* ============================================================ *)
(* RunFullIteration                                              *)
(* Runs the complete bottom-up iteration                         *)
(* ============================================================ *)

RunFullIteration[ftData_Association, outputDir_String:""] :=
Module[{result = ftData, k, dir, nLevels},
  nLevels = ftData["NumLevels"];
  dir = If[outputDir === "",
    FileNameJoin[{FeynmanTrick`Private`$FTConfig["WorkDirectory"],
                  ftData["TopTopology"]["Name"] <> "_iteration"}],
    outputDir
  ];

  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["\n===================================="];
    Print["Starting full Feynman trick iteration"];
    Print["Levels: ", nLevels];
    Print["Output: ", dir];
    Print["====================================\n"];
  ];

  (* Build all levels *)
  result = BuildAllLevels[result];

  (* Compute top-to-bottom so each closed basis supplies the seeds for the
     next Feynman-trick level.  Transport still runs bottom-up later. *)
  Do[
    result["Levels"][k]["SeedMasters"] = FeynmanTrickSeedMasters[
      result["Levels"][k - 1]["Masters"],
      result["Levels"][k]["CombinedPositions"]
    ];
    If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
      Print["Level ", k, " seed masters: ",
        Length[result["Levels"][k]["SeedMasters"]]];
    ];
    result = ComputeLevelData[result, k];
    ExportLevel[result, k, dir, "both"];
  , {k, 1, nLevels}];

  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["\n===================================="];
    Print["Iteration complete."];
    Print["===================================="];
  ];

  result
];


(* ============================================================ *)
(* Accessor functions                                            *)
(* ============================================================ *)

GetLevelTopology[ftData_Association, level_Integer] :=
  ftData["Levels"][level]["Topology"];

GetLevelMatrix[ftData_Association, level_Integer] :=
  ftData["Levels"][level]["DiffMatrix"];


End[];
EndPackage[];
