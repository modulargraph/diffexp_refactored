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
  "IntegrateLevelMaster[transportResult, masterIdx, v1, v2, ibpCoefficients, epsOrder] \
integrates the Feynman trick recursion for a single master: \
Gamma(v1+v2)/(Gamma(v1)*Gamma(v2)) * Integral[x^(v1-1)*(1-x)^(v2-1) * sum_j c_j(x)*f_j(x), {x,0,1}]. \
Returns the integrated boundary value as a list of eps-order coefficients.";

ComputeLevelBoundary::usage =
  "ComputeLevelBoundary[ftData, level, transportResult, epsOrder] computes boundary \
conditions for all masters at 'level' using the transport results from level+1. \
Handles IBP reductions and the Feynman trick recursion formula. \
Returns <|\"BoundaryValues\" -> {...}, \"EpsPrefactors\" -> {...}|>.";

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
Module[{prefactorSpec, result, gammaPrefactor, eps},

  eps = FeynmanTrick`FTeps;

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
  (* We need to modify the transport result to point to a single integral *)
  Module[{singleMasterData, segData, modifiedSegments},
    segData = transportResult["SegmentData"];

    (* Modify each segment to contain only the requested master *)
    modifiedSegments = Table[
      Module[{seg, seriesData, uncompressed, singleSeries},
        seg = segData[[segIdx]];
        seriesData = seg[[5]];

        (* Uncompress *)
        If[StringQ[seriesData],
          uncompressed = Uncompress[Import[seriesData]];,
          If[Head[seriesData] === String && StringLength[seriesData] > 100,
            uncompressed = Uncompress[seriesData];,
            uncompressed = seriesData;
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

    (* Result is indexed as {{eps0, eps1, ...}} for 1 integral *)
    (* Multiply by Gamma prefactor *)
    If[ListQ[result] && Length[result] >= 1,
      gammaPrefactor * result[[1]],
      Print["Warning: Integration returned unexpected format."];
      Table[0, {epsOrder + 1}]
    ]
  ]
];


(* ============================================================ *)
(* Compute Level Boundary                                        *)
(* ============================================================ *)

(*
  Computes boundary conditions for all masters at a given level,
  using the transport results from the level above.

  The Feynman trick recursion:
  I^(k-1)_{v1,...} = Gamma(v1+v2)/(Gamma(v1)*Gamma(v2)) *
    Integral[x^(v1-1)*(1-x)^(v2-1) * I^(k)_{v1+v2,...}, {x,0,1}]

  Special cases:
  I^(k-1)_{0,0,...} = I^(k)_{0,...}  (no integration)
  I^(k-1)_{v1,0,...} = lim_{x->1} I^(k)_{v1,...}
  I^(k-1)_{0,v2,...} = lim_{x->0} I^(k)_{v2,...}
*)
ComputeLevelBoundary[ftData_Association, level_Integer,
    transportResult_Association, epsOrder_Integer] :=
Module[{levelData, levelAbove, mastersAtLevel, mastersAbove,
        combinedPositions, neededIntegrals, bcValues,
        topologyAbove, reducedIntegrals, v1, v2},

  levelData = ftData["Levels"][level];
  levelAbove = ftData["Levels"][level + 1];
  mastersAtLevel = levelData["Masters"];
  mastersAbove = levelAbove["Masters"];
  combinedPositions = levelAbove["CombinedPositions"];

  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["Computing boundary for level ", level, " from level ", level + 1];
    Print["  Masters at level ", level, ": ", Length[mastersAtLevel]];
    Print["  Masters at level ", level + 1, ": ", Length[mastersAbove]];
    Print["  Combined positions: ", combinedPositions];
  ];

  (* For each master at the current level, identify what integrals
     at the level above are needed, and how they reduce to masters above *)
  topologyAbove = levelAbove["Topology"];

  bcValues = Table[
    Module[{masterVec, integralsNeeded, totalBC, integral, reduced,
            ibpCoeffs, masterContribIdx, masterContrib},

      masterVec = mastersAtLevel[[masterIdx]];

      (* Identify needed integrals at level+1 via the Feynman trick recursion *)
      integralsNeeded = FeynmanTrick`FeynmanTrickIteration`IdentifyNeededIntegrals[
        ftData, level, masterVec
      ];

      If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 2,
        Print["  Master ", masterIdx, " (", masterVec, ") needs: ", Length[integralsNeeded], " integrals"];
      ];

      (* Sum contributions from all needed integrals *)
      totalBC = Table[0, {epsOrder + 1}];

      Do[
        Module[{neededVec, neededV1, neededV2, isSpecialCase},
          neededVec = integralsNeeded[[intIdx, 1]];
          neededV1 = integralsNeeded[[intIdx, 2]];
          neededV2 = integralsNeeded[[intIdx, 3]];

          (* Check for special cases *)
          isSpecialCase = Which[
            neededV1 == 0 && neededV2 == 0, "zero",
            neededV1 > 0 && neededV2 == 0, "limitUpper",
            neededV1 == 0 && neededV2 > 0, "limitLower",
            True, "integrate"
          ];

          (* Reduce the needed integral to masters at level+1 *)
          reduced = FeynmanTrick`FIREInterface`ReduceIntegrals[
            topologyAbove,
            {neededVec}
          ];

          If[reduced =!= $Failed && Length[reduced] > 0,
            Module[{reduction, ibpPairs},
              reduction = reduced[[1]];  (* IBP reduction of the first (only) integral *)

              (* reduction is a list of {coefficient, masterIndex} pairs *)
              (* or similar format from FIRE *)

              Switch[isSpecialCase,
                "zero",
                  (* Nothing to add *)
                  Null,

                "limitLower",
                  (* lim_{x->0} I^(k)_{v2,...} *)
                  (* Use EvaluateLimitAtSingularity on the transport result *)
                  Module[{limitVal},
                    Do[
                      Module[{coeff, mIdx, segData, firstSeg, seriesAtMaster, decomp, limitResult},
                        {coeff, mIdx} = ibpPairs[[pairIdx]];
                        (* Get the first segment (near x=0) *)
                        segData = transportResult["SegmentData"];
                        firstSeg = segData[[1]];
                        (* Decompose and evaluate limit *)
                        (* ... simplified for now *)
                      ],
                      {pairIdx, Length[ibpPairs]}
                    ];
                  ],

                "limitUpper",
                  (* lim_{x->1} I^(k)_{v1,...} *)
                  (* Similar to limitLower but at x=1 *)
                  Null,

                "integrate",
                  (* Full integration: Gamma(v1+v2)/(Gamma(v1)*Gamma(v2)) *
                     Integral[x^(v1-1)*(1-x)^(v2-1) * sum_j c_j(x) * f_j(x)] *)

                  (* For each master j at level+1 with IBP coefficient c_j *)
                  Do[
                    Module[{coeff, mIdx, contribution},
                      (* coeff is the IBP coefficient (may be a function of xx) *)
                      (* mIdx is the index into mastersAbove *)
                      {coeff, mIdx} = reduction[[pairIdx]];

                      contribution = IntegrateLevelMaster[
                        transportResult, mIdx,
                        neededV1, neededV2,
                        coeff, epsOrder
                      ];

                      totalBC = totalBC + contribution;
                    ],
                    {pairIdx, Length[reduction]}
                  ];
              ];
            ],
            If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
              Print["  Warning: IBP reduction failed for ", neededVec];
            ];
          ];
        ],
        {intIdx, Length[integralsNeeded]}
      ];

      totalBC
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

  (* Step 3: Transport and integrate level by level *)
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

    (* Compute boundary for level-1 by integration *)
    If[level > 1,
      levelBoundary = ComputeLevelBoundary[
        updatedFtData, level - 1, transportResult, epsOrder
      ];

      If[AssociationQ[levelBoundary],
        currentBCs = levelBoundary["BoundaryValues"];
        If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
          Print["  Level ", level - 1, " boundary computed: ", Length[currentBCs], " masters"];
        ];
      ,
        Print["Error: ComputeLevelBoundary failed at level ", level - 1];
        Return[$Failed];
      ];
    ,
      (* At level 1: the integration gives boundary for level 0 (the original topology) *)
      levelBoundary = ComputeLevelBoundary[
        updatedFtData, 0, transportResult, epsOrder
      ];
      currentBCs = If[AssociationQ[levelBoundary],
        levelBoundary["BoundaryValues"],
        $Failed
      ];
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
