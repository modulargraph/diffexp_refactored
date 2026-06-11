(* Master-by-master FeynmanTrick boundary output for comparison with pySecDec. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
$Path = DeleteDuplicates[Prepend[$Path, repoRoot]];

Quiet[Get["FeynmanTrick/FeynmanTrick.m"], {General::shdw, Symbol::shdw}];
Quiet[Get["DiffExp.m"], {General::shdw, Symbol::shdw}];
Get[FileNameJoin[{repoRoot, "Scripts", "FTExamples.m"}]];

ClearAll[envOrDefault, laurentCoefficient,
  finiteCoefficient, cleanNumber, cleanExpr, coeffTable, printBoundaryRows,
  runExample];

envOrDefault[name_, default_] := Module[{value = Environment[name]},
  If[StringQ[value] && StringLength[StringTrim[value]] > 0, value, default]
];

wp = ToExpression[envOrDefault["FT_WORKING_PRECISION", "500"]];
epsOrder = ToExpression[envOrDefault["FT_EPS_ORDER", "0"]];
expansionOrder = ToExpression[envOrDefault["FT_EXPANSION_ORDER", "50"]];
divisionOrder = ToExpression[envOrDefault["FT_DIVISION_ORDER", "4"]];
boundaryExtraOrder = ToExpression[envOrDefault["FT_BOUNDARY_EXTRA_ORDER", "4"]];
stopAfterBoundaryLevel = envOrDefault["FT_STOP_AFTER_BOUNDARY_LEVEL", ""];
stopAfterBoundaryLevel = If[StringLength[StringTrim[stopAfterBoundaryLevel]] > 0,
  ToExpression[stopAfterBoundaryLevel],
  Missing["NotSet"]
];

FeynmanTrick`SetFTOption["Threads", 1];
FeynmanTrick`SetFTOption["FThreads", 1];
FeynmanTrick`SetFTOption["Verbosity", 0];
FeynmanTrick`SetFTOption["WorkingPrecision", wp];
FeynmanTrick`SetFTOption["ReductionCache", False];
FeynmanTrick`SetFTOption["FIRETimeoutSeconds", 1800];
FeynmanTrick`SetFTOption["IntegrationPoleAllowance",
  ToExpression[envOrDefault["FT_POLE_ALLOWANCE", "4"]]];

(* Solver self-check diagnostics (slow; banana-era instrumentation):
   GENCHECK verifies every solve piece against the block ODE, RESIDSER
   prints the escaping residual series. *)
If[Environment["DEBUG_FUCHS_CHECK"] === "1",
  DiffExp`State`$DebugFuchsianizedCheck = True];
If[Environment["DEBUG_BLOCK_RESID"] === "1",
  DiffExp`State`$DebugBlockResidualSeries = True];

laurentCoefficient[laur_Association, power_Integer] := Module[
  {idx = power - laur["MinPower"] + 1},
  If[idx >= 1 && idx <= Length[laur["Coefficients"]],
    laur["Coefficients"][[idx]],
    0
  ]
];

finiteCoefficient[result_Association] := Module[{raw},
  If[KeyExistsQ[result, "BoundaryValues"] &&
      MatrixQ[result["BoundaryValues"]] &&
      Length[result["BoundaryValues"]] >= 1 &&
      Length[result["BoundaryValues"][[1]]] >= 1,
    result["BoundaryValues"][[1, 1]],
    raw = result["RawBoundaryValues"][[1]];
    laurentCoefficient[raw, 0]
  ]
];

cleanNumber[value_] := Module[{n = N[value, 50], re, im},
  If[!NumericQ[n], Return[value]];
  re = Re[n];
  im = Im[n];
  If[Abs[im] < 10^-80,
    re,
    <|"Re" -> re, "Im" -> im|>
  ]
];

cleanExpr[expr_Association] := AssociationThread[Keys[expr], cleanExpr /@ Values[expr]];
cleanExpr[expr_List] := cleanExpr /@ expr;
cleanExpr[expr_?NumericQ] := cleanNumber[expr];
cleanExpr[expr_] := expr;

coeffTable[raw_Association, minPower_Integer, maxPower_Integer] :=
  Table[{p, cleanNumber[laurentCoefficient[raw, p]]}, {p, minPower, maxPower}];

printBoundaryRows[example_String, level_Integer, masters_List, boundary_Association] := Module[
  {rawValues, prefactors, rawMin, rows, rowMin},
  rawValues = boundary["RawBoundaryValues"];
  prefactors = boundary["EpsPrefactors"];
  rawMin = boundary["RawMinPower"];
  rows = Table[
    rowMin = Lookup[rawValues[[i]], "MinPower", rawMin];
    <|
      "Example" -> example,
      "Level" -> level,
      "Master" -> masters[[i]],
      "EpsPrefactor" -> prefactors[[i]],
      "RawMinPower" -> rowMin,
      "Coefficients" -> coeffTable[rawValues[[i]], rowMin, 0]
    |>,
    {i, Length[masters]}
  ];
  Do[
    Print[
      "STEPWISE ",
      ExportString[cleanExpr[row], "RawJSON", "Compact" -> True]
    ],
    {row, rows}
  ];
];

runExample[name_String] := Module[
  {
    topology, ftData, outputDir, nLevels, boundaryOrder, deepBoundary,
    currentBCs, currentPrefactors, transportOrder, matrixDir,
    extraSingularFactors, transportResult, levelBoundary, finalBoundary
  },
  Print["EXAMPLE ", name];
  FeynmanTrick`SetFTOption["DimensionExpression", FTExampleDimension[name]];
  topology = FTExampleTopology[name, "step"];
  If[topology === $Failed, Return[$Failed, Module]];
  ftData = FeynmanTrick`FeynmanTrickIteration`DefineFTIteration[
    topology, FTExampleSequence[name], {}
  ];
  outputDir = FileNameJoin[{$TemporaryDirectory,
    "FT_stepwise_" <> name <> "_" <> ToString[$ProcessID]}];
  If[DirectoryQ[outputDir], DeleteDirectory[outputDir, DeleteContents -> True]];
  CreateDirectory[outputDir, CreateIntermediateDirectories -> True];
  ftData = FeynmanTrick`FeynmanTrickIteration`RunFullIteration[ftData, outputDir];
  If[ftData === $Failed, Return[$Failed]];

  nLevels = ftData["NumLevels"];
  boundaryOrder = epsOrder + nLevels + boundaryExtraOrder;
  deepBoundary = FeynmanTrick`BoundaryConditions`DeepestLevelBoundary[
    ftData, boundaryOrder
  ];
  If[!AssociationQ[deepBoundary], Return[$Failed]];

  printBoundaryRows[
    name, nLevels, ftData["Levels"][nLevels]["Masters"],
    <|
      "RawBoundaryValues" -> Map[
        <|
          "MinPower" -> -deepBoundary["EpsPrefactors"][[#]],
          "Coefficients" -> deepBoundary["BoundaryValues"][[#]]
        |> &,
        Range[Length[deepBoundary["BoundaryValues"]]]
      ],
      "EpsPrefactors" -> deepBoundary["EpsPrefactors"],
      "RawMinPower" -> Min[-deepBoundary["EpsPrefactors"]]
    |>
  ];

  currentBCs = deepBoundary["BoundaryValues"];
  currentPrefactors = deepBoundary["EpsPrefactors"];
  Do[
    (* Carry ALL available boundary depth through the chain: the
       single-level RequiredTransportEpsilonOrder underestimates deep
       towers (each level's IBP eps-poles and prefactor shifts consume
       orders cumulatively), and cutting here starves the levels below
       (observed: box_bubble L2 capped at eps order 2 while needing 9
       despite a generous deepest-level budget).  requiredOrder is kept
       as a lower-bound diagnostic only. *)
    Module[{requiredOrder =
        FeynmanTrick`DiffExpIntegration`Private`RequiredTransportEpsilonOrder[
          ftData, level, epsOrder, currentPrefactors
        ]},
      transportOrder = Length[First[currentBCs]] - 1;
      If[transportOrder < requiredOrder,
        Print["WARNING level ", level, " transport depth ", transportOrder,
          " is below the single-level requirement ", requiredOrder,
          "; increase FT_BOUNDARY_EXTRA_ORDER."];
      ];
    ];
    FeynmanTrick`FeynmanTrickIteration`ExportLevel[
      ftData, level, outputDir, "diffexp", transportOrder
    ];
    matrixDir = FileNameJoin[{outputDir, "Level_" <> ToString[level] <> "_Matrices"}];
    extraSingularFactors =
      FeynmanTrick`DiffExpIntegration`CollectLevelIBPSingularFactors[
        ftData, level
      ];
    Print["EXTRA_SINGULAR_FACTORS level ", level, ": ",
      InputForm[extraSingularFactors]];
    transportResult = FeynmanTrick`DiffExpIntegration`TransportLevel[
      matrixDir, currentBCs, transportOrder,
      "WorkingPrecision" -> wp,
      "ExpansionOrder" -> expansionOrder,
      "DivisionOrder" -> divisionOrder,
      "Verbosity" -> ToExpression[envOrDefault["FT_TRANSPORT_VERBOSITY", "0"]],
      "EpsPrefactors" -> currentPrefactors,
      "ExtraSingularFactors" -> extraSingularFactors,
      "UseRationalRecurrence" -> True
    ];
    If[transportResult === $Failed, Return[$Failed]];
    transportResult["BoundaryValuesAbove"] = currentBCs;
    transportResult["EpsPrefactorsAbove"] = currentPrefactors;
    Module[{saveDir = Environment["FT_SAVE_TRANSPORT_DIR"]},
      If[StringQ[saveDir] && StringLength[StringTrim[saveDir]] > 0,
        If[!DirectoryQ[saveDir],
          CreateDirectory[saveDir, CreateIntermediateDirectories -> True]
        ];
        Put[
          <|
            "Level" -> level,
            "TransportResult" -> transportResult,
            "TransportOrder" -> transportOrder,
            "MatrixDir" -> matrixDir
          |>,
          FileNameJoin[{saveDir,
            "transport_level_" <> ToString[level] <> ".m"}]
        ];
        Print["SAVED_TRANSPORT level ", level, " -> ", saveDir];
      ];
    ];
    levelBoundary = FeynmanTrick`DiffExpIntegration`ComputeLevelBoundary[
      ftData, level - 1, transportResult, epsOrder
    ];
    If[!AssociationQ[levelBoundary], Return[$Failed]];
    printBoundaryRows[
      name, level - 1, ftData["Levels"][level - 1]["Masters"],
      levelBoundary
    ];
    If[IntegerQ[stopAfterBoundaryLevel] && level - 1 === stopAfterBoundaryLevel,
      Print["STOPPED_AFTER_BOUNDARY_LEVEL ", stopAfterBoundaryLevel];
      Return[<|"FtData" -> ftData, "Boundary" -> levelBoundary|>]
    ];
    currentBCs = levelBoundary["BoundaryValues"];
    currentPrefactors = levelBoundary["EpsPrefactors"];
    finalBoundary = levelBoundary;
    ,
    {level, nLevels, 1, -1}
  ];
  If[AssociationQ[finalBoundary],
    Print[
      "FINAL ",
      ExportString[
        cleanExpr[<|
          "Example" -> name,
          "Finite" -> N[finiteCoefficient[finalBoundary], 50],
          "RawMinPower" -> finalBoundary["RawMinPower"]
        |>],
        "RawJSON", "Compact" -> True
      ]
    ];
  ];
];

requested = StringTrim /@ StringSplit[envOrDefault["FT_EXAMPLES", "bubble,sunrise"], ","];
Do[
  If[runExample[name] === $Failed,
    Print["FAILED ", name];
    Quit[1]
  ],
  {name, requested}
];
Quit[0];
