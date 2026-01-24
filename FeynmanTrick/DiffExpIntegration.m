(* ::Package:: *)
(* DiffExpIntegration - Bridge between FeynmanTrick and DiffExp *)
(* Handles transport, integration with Feynman trick prefactors, *)
(* and the full bottom-up integration pipeline.                   *)

BeginPackage["FeynmanTrick`DiffExpIntegration`", {"FeynmanTrick`"}];

TransportLevel::usage =
  "TransportLevel[matrixDir, boundaryValues, epsOrder, opts] loads DiffExp with \
matrices from matrixDir, sets boundary conditions at the fixed parameter point, \
and transports to cover [0,1]. Returns the TransportTo result with SegmentData.";

IntegrateLevelMaster::usage =
  "IntegrateLevelMaster[transportResult, masterIdx, v1, v2, ibpCoeff, epsOrder] \
integrates the Feynman trick recursion for a single master: \
Gamma(v1+v2)/(Gamma(v1)*Gamma(v2)) * Integral[x^(v1-1)*(1-x)^(v2-1) * c(x)*f(x), {x,0,1}]. \
Returns the integrated boundary value as a list of eps-order coefficients.";

EvaluateLimitFromTransport::usage =
  "EvaluateLimitFromTransport[transportResult, ibpCoeffs, boundary, epsOrder] \
evaluates lim_{x->boundary} of the linear combination sum_j ibpCoeffs[[j]] * f_j(x). \
boundary = 0 or 1. Uses segment decomposition, keeping only pure Taylor terms \
(setting x^{a+b*eps} terms with b!=0 to zero). Returns eps-order coefficient list.";

ComputeLevelBoundary::usage =
  "ComputeLevelBoundary[ftData, level, transportResult, epsOrder] computes boundary \
conditions for all masters at 'level' using the transport results from level+1. \
Handles IBP reductions and the Feynman trick recursion formula. \
Returns <|\"BoundaryValues\" -> {...}, \"Masters\" -> {...}, \"Level\" -> ...|>.";

RunIntegrationPipeline::usage =
  "RunIntegrationPipeline[ftData, outputDir, epsOrder, opts] runs the full bottom-up \
integration pipeline: computes matrices at all levels, evaluates boundary at deepest \
level, transports and integrates level by level. Returns final boundary conditions for level 0.";

Begin["`Private`"];

(* ============================================================ *)
(* DiffExp Loading and Transport                                 *)
(* ============================================================ *)

(*
  TransportLevel loads DiffExp fresh for a given level's matrices,
  sets boundary conditions, and transports to cover [0,1].

  The transport goes from fixedParamValue (default 11/23) outward to
  cover the full interval. DiffExp automatically handles singularities.
*)
TransportLevel[matrixDir_String, boundaryValues_List, epsOrder_Integer,
    opts:OptionsPattern[{
      "FixedParamValue" -> 11/23,
      "WorkingPrecision" -> 500,
      "ExpansionOrder" -> 50,
      "Verbosity" -> 1
    }]] :=
Module[{fixedVal, precision, expOrder, verbosity,
        diffExpConfig, bcs, startPoint, result,
        diffExpPath, resultToLower, resultToUpper,
        savedDataLower, savedDataUpper, combinedSegments},

  fixedVal = OptionValue["FixedParamValue"];
  precision = OptionValue["WorkingPrecision"];
  expOrder = OptionValue["ExpansionOrder"];
  verbosity = OptionValue["Verbosity"];

  (* Determine DiffExp path *)
  diffExpPath = FileNameJoin[{
    ParentDirectory[DirectoryName[$InputFileName]],
    "DiffExp.m"
  }];

  (* Load DiffExp if not already loaded *)
  If[!ValueQ[DiffExp`State`FEC],
    If[FileExistsQ[diffExpPath],
      Block[{$ContextPath},
        Quiet[Get[diffExpPath], {General::shdw, Symbol::shdw}];
      ];
    ,
      Print["Error: DiffExp.m not found at ", diffExpPath];
      Return[$Failed];
    ];
  ];

  (* Configure DiffExp for this level *)
  diffExpConfig = {
    MatrixDirectory -> matrixDir,
    System`WorkingPrecision -> precision,
    ExpansionOrder -> expOrder,
    EpsilonOrder -> epsOrder,
    UseMobius -> False,  (* Required for integration! *)
    UsePade -> False,
    Verbosity -> verbosity,
    SegmentationStrategy -> "Predivision"
  };

  If[verbosity >= 1,
    Print["  Loading DiffExp matrices from: ", matrixDir];
  ];

  DiffExp`LoadConfiguration[diffExpConfig];

  (* Prepare boundary conditions at xx = fixedVal *)
  (* boundaryValues is a list of lists: {{bc_eps0, bc_eps1, ...}, ...} *)
  (* One entry per master integral *)
  bcs = boundaryValues;

  (* Transport from fixedVal towards 0 (lower bound) *)
  If[verbosity >= 1,
    Print["  Transporting from xx=", fixedVal, " towards 0..."];
  ];

  startPoint = <|Global`xx -> SetPrecision[fixedVal, precision]|>;

  (* Transport to xx = 0 (or close to it) *)
  resultToLower = Quiet[
    DiffExp`TransportTo[
      {startPoint, bcs},
      <|Global`xx -> 0|>,
      1,  (* endpoint *)
      True  (* SaveExpansions *)
    ],
    {General::shdw, Symbol::shdw}
  ];

  If[verbosity >= 1,
    Print["  Transporting from xx=", fixedVal, " towards 1..."];
  ];

  (* Reload config for transport in other direction *)
  DiffExp`LoadConfiguration[diffExpConfig];

  (* Transport to xx = 1 *)
  resultToUpper = Quiet[
    DiffExp`TransportTo[
      {startPoint, bcs},
      <|Global`xx -> 1|>,
      1,
      True  (* SaveExpansions *)
    ],
    {General::shdw, Symbol::shdw}
  ];

  (* Combine segment data from both transports *)
  If[AssociationQ[resultToLower] && AssociationQ[resultToUpper],
    savedDataLower = resultToLower["SegmentData"];
    savedDataUpper = resultToUpper["SegmentData"];

    (* Reverse the lower segments (they go from fixedVal to 0) *)
    (* and concatenate with upper segments (fixedVal to 1) *)
    combinedSegments = Join[Reverse[savedDataLower], savedDataUpper];

    If[verbosity >= 1,
      Print["  Transport complete: ", Length[combinedSegments], " total segments."];
    ];

    <|
      "SegmentData" -> combinedSegments,
      "LowerResult" -> resultToLower,
      "UpperResult" -> resultToUpper,
      "NumIntegrals" -> Length[bcs],
      "EpsilonOrder" -> epsOrder,
      "FixedParamValue" -> fixedVal
    |>
  ,
    Print["Error: Transport failed."];
    If[verbosity >= 2,
      Print["  Lower result: ", Head[resultToLower]];
      Print["  Upper result: ", Head[resultToUpper]];
    ];
    $Failed
  ]
];


(* ============================================================ *)
(* Integration with Feynman Trick Prefactors                     *)
(* ============================================================ *)

(*
  Integrates the Feynman trick recursion for a single contribution:
  Gamma(v1+v2)/(Gamma(v1)*Gamma(v2)) * Integral[x^(v1-1)*(1-x)^(v2-1) * c(x)*f_j(x), {x,0,1}]

  Where:
  - v1, v2 are the propagator exponents being combined
  - c(x) is the IBP coefficient (rational function in xx)
  - f_j(x) is the j-th master from the transport result
*)
IntegrateLevelMaster[transportResult_Association, masterIdx_Integer,
    v1_Integer, v2_Integer, ibpCoeff_, epsOrder_Integer] :=
Module[{prefactorSpec, result, gammaPrefactor, segData, modifiedSegments,
        singleMasterData},

  (* Gamma prefactor: Gamma(v1+v2)/(Gamma(v1)*Gamma(v2)) *)
  gammaPrefactor = Gamma[v1 + v2] / (Gamma[v1] * Gamma[v2]);

  (* Set up the prefactor specification for DefiniteIntegralWithPrefactor *)
  prefactorSpec = <|
    "PowerAtLower" -> v1 - 1,      (* x^(v1-1) *)
    "PowerAtUpper" -> v2 - 1,      (* (1-x)^(v2-1) *)
    "RationalFactor" -> ibpCoeff,   (* IBP coefficient c(x) *)
    "Variable" -> Global`xx         (* variable in rational prefactor *)
  |>;

  (* Extract segment data for just this master *)
  segData = transportResult["SegmentData"];

  (* Modify each segment to contain only the requested master *)
  modifiedSegments = Table[
    Module[{seg, seriesRaw, uncompressed, singleSeries},
      seg = segData[[segIdx]];
      seriesRaw = seg[[5]];

      (* Uncompress if needed *)
      If[StringQ[seriesRaw],
        uncompressed = Uncompress[Import[seriesRaw]];,
        If[Head[seriesRaw] === String && StringLength[seriesRaw] > 100,
          uncompressed = Uncompress[seriesRaw];,
          uncompressed = seriesRaw;
        ]
      ];

      (* Extract just this master's series *)
      singleSeries = {uncompressed[[masterIdx]]};

      (* Rebuild segment with single master *)
      {seg[[1]], seg[[2]], seg[[3]], seg[[4]], singleSeries}
    ],
    {segIdx, Length[segData]}
  ];

  singleMasterData = <|
    "SegmentData" -> modifiedSegments,
    "NumIntegrals" -> 1,
    "EpsilonOrder" -> epsOrder
  |>;

  (* Call the extended DefiniteIntegral *)
  result = DiffExp`RegularizedIntegration`DefiniteIntegralWithPrefactor[
    singleMasterData,
    {0, 1},
    prefactorSpec
  ];

  (* Result is indexed as result[[intIdx, epsOrd+1]] for 1 integral *)
  (* Multiply by Gamma prefactor *)
  If[ListQ[result] && Length[result] >= 1,
    gammaPrefactor * result[[1]],
    Print["Warning: Integration returned unexpected format: ", Head[result]];
    Table[0, {epsOrder + 1}]
  ]
];


(* ============================================================ *)
(* Evaluate Limit from Transport Result                          *)
(* ============================================================ *)

(*
  Evaluates lim_{x->boundary} sum_j ibpCoeffs[[j]] * f_j(x)
  where f_j are the masters from the transport result.

  Per the paper (section 3.2): "Take the segment centered at x=0
  (or x'=1-x=0), and filter out the finite coefficient of the Taylor
  series g_0(x,eps). Put any contributions of the form x^{a_i+b_i*eps}
  with b_i != 0 to zero (even when a_i < 0)."
*)
EvaluateLimitFromTransport[transportResult_Association, ibpCoeffs_List,
    boundary_Integer, epsOrder_Integer] :=
Module[{segData, targetSeg, uncompressed, numMasters,
        limitValues, seriesAtMaster, decomposition,
        taylorPart, limitVal},

  segData = transportResult["SegmentData"];
  numMasters = transportResult["NumIntegrals"];

  (* Select the boundary segment *)
  If[boundary === 0,
    (* First segment: nearest to x=0 *)
    targetSeg = First[SortBy[segData, Min[#[[3]]] &]];
  ,
    (* Last segment: nearest to x=1 *)
    targetSeg = First[SortBy[segData, -Max[#[[3]]] &]];
  ];

  (* Uncompress series data *)
  If[StringQ[targetSeg[[5]]],
    uncompressed = Uncompress[Import[targetSeg[[5]]]];,
    If[Head[targetSeg[[5]]] === String && StringLength[targetSeg[[5]]] > 100,
      uncompressed = Uncompress[targetSeg[[5]]];,
      uncompressed = targetSeg[[5]];
    ]
  ];

  (* For each master, evaluate the limit *)
  limitValues = Table[0, {epsOrder + 1}];

  Do[
    If[ibpCoeffs[[mIdx]] =!= 0,
      (* Get series for this master (list of eps orders) *)
      seriesAtMaster = uncompressed[[mIdx]];

      (* Decompose into x^{a+b*eps} * g(x,eps) terms *)
      decomposition = DiffExp`SingularityDecomposition`DecomposeSingularity[seriesAtMaster];

      (* Extract the pure Taylor contribution: a >= 0 and b == 0 *)
      taylorPart = Select[decomposition, (#["a"] >= 0 && #["b"] == 0) &];

      (* The limit at x=0 is the constant term of the Taylor series g *)
      (* For a=0, b=0: g(0,eps) is the constant coefficient of the series *)
      limitVal = Table[0, {epsOrder + 1}];

      Do[
        Module[{gSeries, constCoeff},
          gSeries = term["g"];
          (* g is a list of SeriesData objects, one per eps order *)
          Do[
            If[epsIdx <= Length[gSeries] && MatchQ[gSeries[[epsIdx]], _SeriesData],
              constCoeff = SeriesCoefficient[gSeries[[epsIdx]], 0];
              If[NumericQ[constCoeff],
                limitVal[[epsIdx]] += constCoeff;
              ];
            ];
          , {epsIdx, Min[epsOrder + 1, Length[gSeries]]}];
        ];
      , {term, taylorPart}];

      (* Add weighted contribution *)
      limitValues += ibpCoeffs[[mIdx]] * limitVal;
    ];
  , {mIdx, numMasters}];

  limitValues
];


(* ============================================================ *)
(* Compute Level Boundary                                        *)
(* ============================================================ *)

(*
  Computes boundary conditions for all masters at a given level,
  using the transport results from the level above.

  The Feynman trick recursion (eq. 2.8 of the paper):
  I^(k-1)_{v1,...} = Gamma(v1+v2)/(Gamma(v1)*Gamma(v2)) *
    Integral[x^(v1-1)*(1-x)^(v2-1) * I^(k)_{v1+v2,...}, {x,0,1}]

  Where the combination at level k merges positions {posI, posJ}.
  Given a master M at level k-1:
  - vi = M[[posI]], vj = M[[posJ]]
  - The needed integral at level k has:
    - Position posI: vi+vj (for integration) or vi/vj (for limits)
    - Position posJ: 0 (absorbed into combined propagator)

  Special cases (eq. 2.10):
  I_{0,0,...}^(k-1) = I_{0,...}^(k)  (direct: no integration)
  I_{v1,0,...}^(k-1) = lim_{x->1} I_{v1,...}^(k)
  I_{0,v2,...}^(k-1) = lim_{x->0} I_{v2,...}^(k)
*)
ComputeLevelBoundary[ftData_Association, level_Integer,
    transportResult_Association, epsOrder_Integer] :=
Module[{levelData, levelAbove, mastersAtLevel, mastersAbove,
        combinedPositions, posI, posJ, topologyAbove, bcValues},

  levelData = ftData["Levels"][level];
  levelAbove = ftData["Levels"][level + 1];
  mastersAtLevel = levelData["Masters"];
  mastersAbove = levelAbove["Masters"];
  combinedPositions = levelAbove["CombinedPositions"];
  {posI, posJ} = combinedPositions;
  topologyAbove = levelAbove["Topology"];

  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["Computing boundary for level ", level, " from level ", level + 1];
    Print["  Masters at level ", level, ": ", Length[mastersAtLevel]];
    Print["  Masters at level ", level + 1, ": ", Length[mastersAbove]];
    Print["  Combined positions: {", posI, ", ", posJ, "}"];
  ];

  (* For each master at the current level *)
  bcValues = Table[
    Module[{masterVec, vi, vj, neededVec, case, reduction, expr,
            ibpCoeffs, totalBC},

      masterVec = mastersAtLevel[[masterIdx]];
      vi = masterVec[[posI]];
      vj = masterVec[[posJ]];

      (* Determine which case of the recursion applies *)
      case = Which[
        vi > 0 && vj > 0, "integrate",
        vi > 0 && vj == 0, "limitUpper",
        vi == 0 && vj > 0, "limitLower",
        True, "direct"
      ];

      If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 2,
        Print["  Master ", masterIdx, " (", masterVec, "): case=", case,
              " vi=", vi, " vj=", vj];
      ];

      (* Construct the needed integral at level+1 *)
      neededVec = masterVec;
      Switch[case,
        "integrate",
          (* Combined propagator gets power vi+vj, position j gets 0 *)
          neededVec[[posI]] = vi + vj;
          neededVec[[posJ]] = 0;,
        "limitUpper",
          (* At x=1: D_combined = D_i. Position i has vi, j has 0 *)
          neededVec[[posI]] = vi;
          neededVec[[posJ]] = 0;,
        "limitLower",
          (* At x=0: D_combined = D_j. Position i has vj, j has 0 *)
          neededVec[[posI]] = vj;
          neededVec[[posJ]] = 0;,
        "direct",
          (* Both absent: same integral with j set to 0 *)
          neededVec[[posJ]] = 0;
      ];

      (* Reduce the needed integral to masters at level+1 via IBP *)
      reduction = FeynmanTrick`FIREInterface`ReduceIntegrals[
        topologyAbove,
        {neededVec}
      ];

      If[reduction === $Failed,
        If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
          Print["  Warning: IBP reduction failed for ", neededVec];
        ];
        Return[Table[0, {epsOrder + 1}], Module];
      ];

      (* Extract the reduction expression *)
      expr = reduction[neededVec];

      (* Extract coefficient of each master G[1, masters_j] *)
      ibpCoeffs = Table[
        Coefficient[expr, Global`G[1, mastersAbove[[j]]]],
        {j, Length[mastersAbove]}
      ];

      If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 2,
        Print["    IBP coefficients: ", ibpCoeffs];
      ];

      (* Compute the boundary value based on the case *)
      Switch[case,
        "integrate",
          (* Full integration: sum_j ibpCoeffs[[j]] * IntegrateLevelMaster[...] *)
          totalBC = Table[0, {epsOrder + 1}];
          Do[
            If[ibpCoeffs[[j]] =!= 0,
              Module[{contribution},
                contribution = IntegrateLevelMaster[
                  transportResult, j, vi, vj, ibpCoeffs[[j]], epsOrder
                ];
                totalBC = totalBC + contribution;
              ];
            ];
          , {j, Length[mastersAbove]}];
          totalBC,

        "limitUpper",
          (* lim_{x->1} sum_j ibpCoeffs[[j]] * f_j(x) *)
          EvaluateLimitFromTransport[transportResult, ibpCoeffs, 1, epsOrder],

        "limitLower",
          (* lim_{x->0} sum_j ibpCoeffs[[j]] * f_j(x) *)
          EvaluateLimitFromTransport[transportResult, ibpCoeffs, 0, epsOrder],

        "direct",
          (* Evaluate at the fixed parameter value *)
          (* The "direct" case means both propagators are absent, *)
          (* so this integral at level+1 equals the one at level directly *)
          (* Evaluate at the fixed point from the boundary values *)
          Module[{directVal},
            directVal = Table[0, {epsOrder + 1}];
            Do[
              If[ibpCoeffs[[j]] =!= 0,
                (* Use the boundary values from level+1 *)
                (* These were already computed for the transport *)
                directVal += ibpCoeffs[[j]] * transportResult["BoundaryValuesAbove"][[j]];
              ];
            , {j, Length[mastersAbove]}];
            directVal
          ]
      ]
    ],
    {masterIdx, Length[mastersAtLevel]}
  ];

  <|
    "BoundaryValues" -> bcValues,
    "Masters" -> mastersAtLevel,
    "Level" -> level
  |>
];


(* ============================================================ *)
(* Full Integration Pipeline                                     *)
(* ============================================================ *)

RunIntegrationPipeline[ftData_Association, outputDir_String, epsOrder_Integer:4,
    opts:OptionsPattern[{
      "WorkingPrecision" -> 500,
      "ExpansionOrder" -> 50
    }]] :=
Module[{nLevels, currentBCs, currentPrefactors, matrixDir,
        transportResult, levelBoundary, precision, expOrder,
        updatedFtData},

  precision = OptionValue["WorkingPrecision"];
  expOrder = OptionValue["ExpansionOrder"];
  nLevels = ftData["NumLevels"];

  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["\n=== Running Full Integration Pipeline ==="];
    Print["  Levels: ", nLevels];
    Print["  Epsilon order: ", epsOrder];
    Print["  Working precision: ", precision];
  ];

  (* Step 1: Run the iteration to compute all matrices *)
  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["\n--- Phase 1: Computing differential matrices ---"];
  ];

  updatedFtData = FeynmanTrick`FeynmanTrickIteration`RunFullIteration[ftData, outputDir];
  If[updatedFtData === $Failed,
    Print["Error: RunFullIteration failed."];
    Return[$Failed];
  ];

  (* Step 2: Compute boundary at deepest level *)
  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["\n--- Phase 2: Computing boundary at deepest level ---"];
  ];

  Module[{deepBoundary},
    deepBoundary = FeynmanTrick`BoundaryConditions`DeepestLevelBoundary[
      updatedFtData, epsOrder
    ];

    If[deepBoundary === $Failed,
      Print["Error: DeepestLevelBoundary failed."];
      Return[$Failed];
    ];

    currentBCs = deepBoundary["BoundaryValues"];
    currentPrefactors = deepBoundary["EpsPrefactors"];

    If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
      Print["  Deepest level boundary computed successfully."];
      Print["  Eps prefactors: ", currentPrefactors];
    ];
  ];

  (* Step 3: Transport and integrate level by level (bottom-up) *)
  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["\n--- Phase 3: Transport and integration ---"];
  ];

  Do[
    If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
      Print["\n  Processing level ", level, "..."];
    ];

    (* Get matrix directory for this level *)
    matrixDir = FileNameJoin[{outputDir, "Level_" <> ToString[level] <> "_Matrices"}];

    If[!DirectoryQ[matrixDir],
      Print["Error: Matrix directory not found: ", matrixDir];
      Return[$Failed];
    ];

    (* Transport from fixed point to cover [0,1] *)
    transportResult = TransportLevel[
      matrixDir, currentBCs, epsOrder,
      "WorkingPrecision" -> precision,
      "ExpansionOrder" -> expOrder,
      "Verbosity" -> FeynmanTrick`Private`$FTConfig["Verbosity"]
    ];

    If[transportResult === $Failed,
      Print["Error: Transport failed at level ", level];
      Return[$Failed];
    ];

    (* Store boundary values used for transport (needed for "direct" case) *)
    transportResult["BoundaryValuesAbove"] = currentBCs;

    (* Compute boundary for level-1 by integration *)
    levelBoundary = ComputeLevelBoundary[
      updatedFtData, level - 1, transportResult, epsOrder
    ];

    If[AssociationQ[levelBoundary],
      currentBCs = levelBoundary["BoundaryValues"];
      If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
        Print["  Level ", level - 1, " boundary computed: ",
              Length[currentBCs], " masters"];
      ];
    ,
      Print["Error: ComputeLevelBoundary failed at level ", level - 1];
      Return[$Failed];
    ];
    ,
    {level, nLevels, 1, -1}
  ];

  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["\n=== Pipeline Complete ==="];
    Print["  Final boundary conditions for level 0: ", Length[currentBCs], " masters"];
  ];

  <|
    "BoundaryValues" -> currentBCs,
    "Level" -> 0,
    "FtData" -> updatedFtData
  |>
];


End[];
EndPackage[];
