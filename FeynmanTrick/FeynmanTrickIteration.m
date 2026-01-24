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
Module[{ftData, nLevels},
  nLevels = Length[combinationSeq];

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
        "Masters" -> {},
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
        prevFixedParams, fixRules, name, result},

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

  param = FeynmanTrick`xx;
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

  (* Fix kinematics to numerical point *)
  If[numericalPoint =!= {},
    newProps = newProps /. numericalPoint;
  ];

  (* Create new topology *)
  name = ftData["TopTopology"]["Name"] <> "_L" <> ToString[level];
  newTopology = DefineTopology[
    name,
    ftData["TopTopology"]["LoopMomenta"],
    ftData["TopTopology"]["ExternalMomenta"],
    newProps,
    {} (* Replacements already applied via numericalPoint *)
  ];

  (* Update ftData *)
  result = ftData;
  result["Levels"][level] = <|
    "Topology" -> newTopology,
    "Propagators" -> newProps,
    "FeynmanParameter" -> param,
    "FixedParams" -> prevFixedParams,
    "CombinedPositions" -> {i, j},
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
(* ComputeLevelData                                              *)
(* Full computation for one level:                               *)
(* SetupFIRE -> FindBasis -> ComputeDiffMatrix                   *)
(* ============================================================ *)

ComputeLevelData[ftData_Association, level_Integer] :=
Module[{levelData, topology, masters, diffMat, param, result, updatedTopo,
        decomp, combo, combinedPos, otherPos},

  If[!KeyExistsQ[ftData["Levels"], level],
    Print["Error: level ", level, " not built. Call BuildLevel first."];
    Return[ftData];
  ];

  levelData = ftData["Levels"][level];
  topology = levelData["Topology"];
  param = levelData["FeynmanParameter"];

  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["\n=== Computing level ", level, " ==="];
  ];

  (* Step 1: Setup FIRE *)
  updatedTopo = SetupFIRE[topology];

  (* Step 2: Find master basis *)
  updatedTopo = FindBasis[updatedTopo];
  masters = updatedTopo["Masters"];

  (* Step 3: Compute differential matrix *)
  (* Use the Feynman trick fast path for the decomposition *)
  combo = levelData["CombinedPositions"];
  If[combo =!= Missing["KeyAbsent", "CombinedPositions"] && ListQ[combo],
    {combinedPos, otherPos} = combo;
    decomp = FeynmanTrickDecomposition[
      topology["NumPropagators"], combinedPos, otherPos, param
    ];
    diffMat = ComputeDiffMatrix[updatedTopo, param, decomp];
  ,
    (* General case: no combination info, use full decomposition *)
    diffMat = ComputeDiffMatrix[updatedTopo, param];
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

  (* Compute from bottom (highest level) to top *)
  Do[
    result = ComputeLevelData[result, k];
    ExportLevel[result, k, dir, "both"];
  , {k, nLevels, 1, -1}];

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
