(*::Package::*)(*DiffExpIntegration - Bridge between FeynmanTrick and
                DiffExp *)(*Handles transport,
                           integration with Feynman trick prefactors,
                               *)(*and the full bottom -
                                  up integration pipeline.*)

    BeginPackage["FeynmanTrick`DiffExpIntegration`", {"FeynmanTrick`"}];

TransportLevel::usage =
    "TransportLevel[matrixDir, boundaryValues, epsOrder, opts] loads DiffExp with \
matrices from matrixDir, sets boundary conditions at the fixed parameter point, \
and transports to cover [0,1]. Returns the TransportTo result with SegmentData. \
Use option \"EpsPrefactors\" -> {k1,...} when boundaryValues are for \
J_i = eps^k_i I_i rather than the raw master basis.";

IntegrateCombinedMasters::usage =
    "IntegrateCombinedMasters[transportResult, ibpCoeffs, v1, v2, epsOrder, epsPrefactors, feynmanParam] \
integrates the combined Feynman trick integrand: \
Gamma(v1+v2)/(Gamma(v1)*Gamma(v2)) * Integral[x^(v1-1)*(1-x)^(v2-1) * Sum_j c_j(x)*f_j(x), {x,0,1}]. \
Combines the integrand before integration to handle cancellation of poles in IBP coefficients. \
Returns the integrated boundary value as a list of eps-order coefficients.";

EvaluateLimitFromTransport::usage =
    "EvaluateLimitFromTransport[transportResult, ibpCoeffs, boundary, epsOrder, epsPrefactors] \
evaluates lim_{x->boundary} of the linear combination sum_j ibpCoeffs[[j]] * f_j(x). \
boundary = 0 or 1. Uses segment decomposition, keeping only pure Taylor terms \
(setting x^{a+b*eps} terms with b!=0 to zero). Returns eps-order coefficient list.";

ComputeLevelBoundary::usage =
    "ComputeLevelBoundary[ftData, level, transportResult, epsOrder] computes boundary \
conditions for all masters at 'level' using the transport results from level+1. \
Handles IBP reductions and the Feynman trick recursion formula. \
Returns <|\"BoundaryValues\" -> {...}, \"Masters\" -> {...}, \"Level\" -> ...|>.";

CollectLevelIBPSingularFactors::usage =
    "CollectLevelIBPSingularFactors[ftData, level] collects Feynman-parameter \
denominator factors from the IBP reductions needed to integrate from level to level-1.";

RunIntegrationPipeline::usage =
    "RunIntegrationPipeline[ftData, outputDir, epsOrder, opts] runs the full bottom-up \
integration pipeline: computes matrices at all levels, evaluates boundary at deepest \
level, transports and integrates level by level. Returns final boundary conditions for level 0.";

Begin["`Private`"];

$DiffExpIntegrationDirectory = DirectoryName[$InputFileName];

(* ============================================================ *)
(* Resolve DiffExp analytic continuation symbols                 *)
(* ============================================================ *)

(*
  DiffExp transport series contain theta symbols (θp, θm) for analytic
  continuation. Resolve to numerical values for +i*delta prescription:
    θp -> 1, θm -> 0  (approaching from upper half plane)

  IMPORTANT: Do NOT substitute Logx -> Log[x]! The Logx symbol must
  remain symbolic for DecomposeSingularity, IntegrateWithLogPower, and
  EvaluateIntegralAtPoint to work correctly. Logx is resolved to Log[x]
  only at the final boundary evaluation step (where Logx -> 0 at x=0
  avoids computing 0 * Log[0] = Indeterminate).
*)
$thetaPlusRules = {
  DiffExp`Symbols`\[Theta]p -> 1,
  DiffExp`Symbols`\[Theta]m -> 0
};

$thetaMinusRules = {
  DiffExp`Symbols`\[Theta]p -> 0,
  DiffExp`Symbols`\[Theta]m -> 1
};

$thetaRules = $thetaPlusRules;

activeNumericPrecision[] := Module[{precision},
  precision = Quiet[
    Check[
      DiffExp`State`FEWorkingPrecision,
      FeynmanTrick`Private`$FTConfig["WorkingPrecision"]
    ]
  ];
  If[IntegerQ[precision] && precision > 0,
    precision,
    500
  ]
];

numericAtActivePrecision[expr_, precision_:Automatic] := Module[
  {p = If[precision === Automatic, activeNumericPrecision[], precision], val},
  val = Quiet[Check[N[expr, p], $Failed]];
  If[val === $Failed,
    val,
    SetPrecision[val, p]
  ]
];

realNumericAtActivePrecision[expr_, tol_:Automatic] := Module[
  {eps = If[tol === Automatic, DiffExp`State`FEC[RationalizationTolerance], tol],
   val},
  val = numericAtActivePrecision[expr];
  If[val === $Failed || !NumericQ[val],
    val,
    If[TrueQ[Abs[Im[val]] < eps], Re[val], val]
  ]
];

thetaRulesForSegment[seg_List] := Module[
  {tol, localBounds, nonzeroBounds, sample},
  tol = DiffExp`State`FEC[RationalizationTolerance];
  localBounds = numericAtActivePrecision[seg[[4]]];
  nonzeroBounds = Select[localBounds,
    !TrueQ[PossibleZeroQ[#]] &&
      !TrueQ[NumericQ[#] && Abs[realNumericAtActivePrecision[#, tol]] < tol] &
  ];
  sample = If[Length[nonzeroBounds] > 0,
    Mean[nonzeroBounds],
    1
  ];
  If[TrueQ[realNumericAtActivePrecision[sample, tol] < 0],
    $thetaMinusRules,
    $thetaPlusRules
  ]
];

thetaRulesAtLocalPoint[pt_, direction_:Automatic] := Module[
  {tol, sign},
  tol = DiffExp`State`FEC[RationalizationTolerance];
  sign = Which[
    direction =!= Automatic && NumericQ[direction],
      Sign[realNumericAtActivePrecision[direction, tol]],
    TrueQ[PossibleZeroQ[pt]] ||
      TrueQ[NumericQ[pt] && Abs[realNumericAtActivePrecision[pt, tol]] < tol], 1,
    TrueQ[realNumericAtActivePrecision[pt, tol] < 0], -1,
    True, 1
  ];

  If[sign < 0,
    $thetaMinusRules,
    $thetaPlusRules
  ]
];

localEndpointDirection[seg_List, localEndpoint_] := Module[
  {tol, localBounds, otherBounds},
  tol = DiffExp`State`FEC[RationalizationTolerance];
  localBounds = seg[[4]];
  otherBounds = Select[localBounds,
    !TrueQ[PossibleZeroQ[# - localEndpoint]] &&
      !TrueQ[NumericQ[# - localEndpoint] &&
        Abs[realNumericAtActivePrecision[# - localEndpoint, tol]] < tol] &
  ];
  If[Length[otherBounds] > 0,
    First[otherBounds] - localEndpoint,
    Automatic
  ]
];

evaluateLocalExpressionAtPoint[expr_, pt_, direction_:Automatic,
    precision_:Automatic] := Module[{xLocal, val, p},
  p = If[precision === Automatic, activeNumericPrecision[], precision];
  xLocal = DiffExp`Symbols`x;
  val = expr /.
      thetaRulesAtLocalPoint[pt, direction] /.
      DiffExp`Symbols`Logx -> Log[xLocal] /.
      xLocal -> SetPrecision[pt, p];
  numericAtActivePrecision[val, p]
];

snapMainExpression[expr_, localBounds_List, snapTargets_List:{0, 1}] := Module[
  {xLocal, tol, snapped, val, target, samplePoints, targets},
  xLocal = DiffExp`Symbols`x;
  tol = Quiet[Check[DiffExp`State`FEC[RationalizationTolerance], 10^-40]];
  If[!NumericQ[tol] || tol <= 0, tol = 10^-40];
  snapped = expr;
  samplePoints = DeleteDuplicates[Join[localBounds, {0}],
    TrueQ[PossibleZeroQ[#1 - #2]] &
  ];
  targets = DeleteDuplicates[snapTargets,
    TrueQ[PossibleZeroQ[#1 - #2]] &
  ];

  Do[
    val = numericAtActivePrecision[snapped /. xLocal -> lb];
    If[val =!= $Failed && NumericQ[val],
      target = SelectFirst[targets,
        TrueQ[Abs[numericAtActivePrecision[val - #]] < tol] &,
        Missing["NoSnapTarget"]
      ];
      If[target =!= Missing["NoSnapTarget"],
        snapped = Expand[snapped + (target - (snapped /. xLocal -> lb))];
      ];
    ];
  , {lb, samplePoints}];

  DiffExp`Utilities`PChop[Expand[snapped]]
];

snapValuesFromFactors[factors_List, variable_] := Module[
  {roots, numericRoots, tol},
  tol = Quiet[Check[DiffExp`State`FEC[RationalizationTolerance], 10^-40]];
  If[!NumericQ[tol] || tol <= 0, tol = 10^-40];
  roots = Flatten[
    Quiet[
      variable /. Solve[# == 0, variable] & /@ factors,
      {Solve::ratnz, Solve::svars}
    ]
  ];
  numericRoots = Select[RootReduce /@ roots,
    With[{n = numericAtActivePrecision[#]},
      n =!= $Failed &&
        NumericQ[n] &&
        TrueQ[Abs[Im[n]] < tol] &&
        TrueQ[Re[n] >= -tol] &&
        TrueQ[Re[n] <= 1 + tol]
    ] &
  ];
  DeleteDuplicates[Chop[Re /@ numericRoots, tol],
    TrueQ[Abs[numericAtActivePrecision[#1 - #2]] < tol] &
  ]
];

appendMatrixFactors[factors_List, verbosity_Integer:0] := Module[
  {cleanFactors},
  cleanFactors = DeleteCases[
    DeleteDuplicates[Factor /@ Flatten[{factors}]],
    0 | 1 | -1
  ];
  If[Length[cleanFactors] > 0,
    DiffExp`State`MatricesIrreducibleFactors =
      DeleteDuplicates[
        Join[Flatten[{DiffExp`State`MatricesIrreducibleFactors}], cleanFactors],
        TrueQ[PossibleZeroQ[Expand[#1 - #2]]] ||
          TrueQ[PossibleZeroQ[Expand[#1 + #2]]] &
      ];
    If[verbosity >= 2,
      Print["  Added extra segmentation factors: ", cleanFactors];
    ];
  ];
  cleanFactors
];

deltaPrescriptionsForFactors[detectedVar_, factors_List:{}, sign_:1] := Module[
  {prescriptions, cleanFactors},
  cleanFactors = DeleteCases[Factor /@ Flatten[{factors}], 0 | 1 | -1];
  prescriptions = Join[
    {
      {detectedVar, sign},
      {1 - detectedVar, sign}
    },
    ({#, sign} & /@ cleanFactors)
  ];
  DeleteDuplicates[prescriptions,
    TrueQ[PossibleZeroQ[Expand[#1[[1]] - #2[[1]]]]] ||
      TrueQ[PossibleZeroQ[Expand[#1[[1]] + #2[[1]]]]] &
  ]
];

(* Resolve theta functions in a single series expression *)
resolveThetas[expr_] := expr /. $thetaRules;

(* Resolve theta functions in segment data (list of segments).
   Each segment is {line, transform, mainBounds, localBounds, seriesData}.
   The seriesData (element 5) may contain theta functions in coefficients. *)
resolveSegmentThetas[segData_List] := Table[
  Module[{seg, rawSeries, resolved},
    seg = segData[[segIdx]];
    rawSeries = seg[[5]];

    (* Uncompress if needed, resolve thetas, recompress *)
    If[StringQ[rawSeries],
      If[FileExistsQ[rawSeries],
        resolved = Uncompress[Import[rawSeries]];,
        resolved = Uncompress[rawSeries];
      ];
      resolved = resolved /. $thetaRules;
      {seg[[1]], seg[[2]], seg[[3]], seg[[4]], Compress[resolved]}
    ,
      resolved = rawSeries /. $thetaRules;
      {seg[[1]], seg[[2]], seg[[3]], seg[[4]], resolved}
    ]
  ],
  {segIdx, Length[segData]}
];

uncompressSeriesData[rawSeries_] := Module[{},
  If[StringQ[rawSeries],
    If[FileExistsQ[rawSeries],
      Uncompress[Import[rawSeries]],
      Uncompress[rawSeries]
    ],
    rawSeries
  ]
];

uncompressSegmentSeries[seg_List] :=
  ReplacePart[seg, 5 -> uncompressSeriesData[seg[[5]]]];

prepareTransportSegments[transportResult_Association] :=
  uncompressSegmentSeries /@ transportResult["SegmentData"];

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
      "DivisionOrder" -> 4,
      "Verbosity" -> 1,
      "EpsPrefactors" -> Automatic,
      "ExtraSingularFactors" -> {},
      "HomogeneousSolve" -> "DontExpand",
      "UseRationalRecurrence" -> True,
      "IntegrationStrategy" -> "Default",
      "EstimateError" -> "Fast",
      "DeltaPrescriptionSign" -> 1,
      "PrescribeMatrixFactors" -> False,
      "LowerEndpoint" -> 0,
      "UpperEndpoint" -> 1
    }]] :=
Module[{fixedVal, precision, expOrder, verbosity,
        diffExpConfig, bcs, startPoint, result,
        diffExpPath, resultToLower, resultToUpper,
        savedDataLower, savedDataUpper, combinedSegments,
        epsPrefactors, divisionOrder, extraSingularFactors,
        remappedExtraFactors = {}, snapValues = {0, 1},
        deltaPrescriptions, homogeneousSolve, useRationalRecurrence,
        integrationStrategy, estimateError, deltaPrescriptionSign,
        prescribeMatrixFactors, lowerEndpoint, upperEndpoint},

  fixedVal = OptionValue["FixedParamValue"];
  precision = OptionValue["WorkingPrecision"];
  expOrder = OptionValue["ExpansionOrder"];
  divisionOrder = OptionValue["DivisionOrder"];
  verbosity = OptionValue["Verbosity"];
  extraSingularFactors = OptionValue["ExtraSingularFactors"];
  homogeneousSolve = OptionValue["HomogeneousSolve"];
  useRationalRecurrence = OptionValue["UseRationalRecurrence"];
  integrationStrategy = OptionValue["IntegrationStrategy"];
  estimateError = OptionValue["EstimateError"];
  deltaPrescriptionSign = OptionValue["DeltaPrescriptionSign"];
  prescribeMatrixFactors = OptionValue["PrescribeMatrixFactors"];
  lowerEndpoint = OptionValue["LowerEndpoint"];
  upperEndpoint = OptionValue["UpperEndpoint"];
  snapValues = DeleteDuplicates[
    Join[
      snapValues,
      Select[{lowerEndpoint, upperEndpoint}, NumericQ]
    ],
    TrueQ[PossibleZeroQ[#1 - #2]] &
  ];
  epsPrefactors = OptionValue["EpsPrefactors"];
  If[!(ListQ[epsPrefactors] && Length[epsPrefactors] == Length[boundaryValues]),
    epsPrefactors = Table[0, {Length[boundaryValues]}];
  ];

  (* Determine DiffExp path *)
  diffExpPath = FileNameJoin[{
    ParentDirectory[$DiffExpIntegrationDirectory],
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

  (* Configure DiffExp for this level.
     All config keys must be fully qualified to match DiffExp's internal symbols,
     since this package's BeginPackage restricts the context path. *)
  (* Delta prescriptions for branch cuts at x=0 and x=1 *)
  (* Format: {polynomial, sign} where sign=1 means +I*delta prescription *)
  (* The variable will be detected after loading matrices *)

  diffExpConfig = {
    DiffExp`State`MatrixDirectory -> matrixDir,
    System`WorkingPrecision -> precision,
    DiffExp`State`ChopPrecision -> precision - 50,
    DiffExp`State`ExpansionOrder -> expOrder,
    DiffExp`State`EpsilonOrder -> epsOrder,
    DiffExp`State`UseMobius -> False,  (* Required for integration! *)
    DiffExp`State`UsePade -> False,
    DiffExp`State`DivisionOrder -> divisionOrder,
    DiffExp`State`Verbosity -> verbosity,
    DiffExp`State`SegmentationStrategy -> "Predivision",
    DiffExp`State`UseRationalRecurrence -> useRationalRecurrence,
    DiffExp`State`IntegrationStrategy -> integrationStrategy,
    "EstimateError" -> estimateError,
    "HomogeneousSolve" -> homogeneousSolve
  };

  If[verbosity >= 1,
    Print["  Loading DiffExp matrices from: ", matrixDir];
  ];

  DiffExp`LoadConfiguration[diffExpConfig];

  (* After loading, DiffExp auto-detects the variable from filenames (dxx_*.m).
     Extract the detected variable symbol so our points match DiffExp's internal state. *)
  Module[{detectedVar},
    detectedVar = First[DiffExp`State`FEC[System`Variables]];

    If[ListQ[extraSingularFactors] && Length[extraSingularFactors] > 0,
      Module[{varName},
        varName = SymbolName[detectedVar];
        remappedExtraFactors = extraSingularFactors /.
          s_Symbol /; SymbolName[s] === varName :> detectedVar;
        remappedExtraFactors = appendMatrixFactors[
          remappedExtraFactors, verbosity
        ];
        snapValues = DeleteDuplicates[
          Join[
            snapValues,
            snapValuesFromFactors[remappedExtraFactors, detectedVar]
          ],
          TrueQ[PossibleZeroQ[#1 - #2]] &
        ];
        If[verbosity >= 2 && Length[snapValues] > 2,
          Print["  Snap values from IBP factors: ", snapValues];
        ];
      ];
    ];

    (* Add delta prescriptions for branch cuts at endpoints and for
       IBP-induced singular factors that were added to the segmentation
       alphabet. The sign is chosen consistently as +i delta. *)
    deltaPrescriptions = deltaPrescriptionsForFactors[
      detectedVar,
      If[TrueQ[prescribeMatrixFactors],
        Join[
          Flatten[{DiffExp`State`MatricesIrreducibleFactors}],
          remappedExtraFactors
        ],
        remappedExtraFactors
      ],
      deltaPrescriptionSign
    ];
    (* Also disable abort on analytic continuation failure - the pipeline
       handles incomplete results gracefully *)
    DiffExp`UpdateConfiguration[{
      DiffExp`State`DeltaPrescriptions -> deltaPrescriptions,
      "AbortOnAnalyticContinuationFail" -> False
    }];

    If[verbosity >= 1,
      Print["  DiffExp detected variable: ", detectedVar, " (context: ", Context[detectedVar], ")"];
      Print["  ExternalScalesVal: ", DiffExp`State`ExternalScalesVal];
      Print["  Boundary values dimensions: ", Dimensions[boundaryValues]];
    ];

    (* Prepare boundary conditions at xx = fixedVal *)
    bcs = boundaryValues;

    (* Transport from fixedVal towards 0 (lower bound) *)
    If[verbosity >= 1,
      If[lowerEndpoint =!= None,
        Print["  Transporting from xx=", fixedVal, " towards ", lowerEndpoint, "..."];
      ];
    ];

    (* Use the detected variable symbol for start/end points.
       DiffExp now handles singular endpoints by returning the series expansion
       instead of trying to evaluate at the singularity. *)
    startPoint = Association[detectedVar -> SetPrecision[fixedVal, precision]];

    (* Transport towards lower endpoint (singular endpoint handled by DiffExp) *)
    resultToLower = If[lowerEndpoint === None,
      Missing["Skipped"],
      DiffExp`Transport`TransportTo[
        {startPoint, bcs},
        Association[detectedVar -> lowerEndpoint],
        1,  (* endpoint *)
        True  (* SaveExpansions *)
      ]
    ];

    If[verbosity >= 2,
      Print["  Lower transport result head: ", Head[resultToLower]];
      If[AssociationQ[resultToLower], Print["  Lower keys: ", Keys[resultToLower]]];
    ];

    If[verbosity >= 1,
      If[upperEndpoint =!= None,
        Print["  Transporting from xx=", fixedVal, " towards ", upperEndpoint, "..."];
      ];
    ];

    (* Reload config for transport in other direction *)
    DiffExp`LoadConfiguration[diffExpConfig];
    If[Length[remappedExtraFactors] > 0,
      appendMatrixFactors[remappedExtraFactors, verbosity];
    ];

    (* CRITICAL: Re-add delta prescriptions after LoadConfiguration reset *)
    DiffExp`UpdateConfiguration[{
      DiffExp`State`DeltaPrescriptions -> deltaPrescriptions,
      "AbortOnAnalyticContinuationFail" -> False
    }];

    (* Transport towards upper endpoint (singular endpoint handled by DiffExp) *)
    resultToUpper = If[upperEndpoint === None,
      Missing["Skipped"],
      DiffExp`Transport`TransportTo[
        {startPoint, bcs},
        Association[detectedVar -> upperEndpoint],
        1,
        True  (* SaveExpansions *)
      ]
    ];

    If[verbosity >= 2,
      Print["  Upper transport result head: ", Head[resultToUpper]];
      If[AssociationQ[resultToUpper], Print["  Upper keys: ", Keys[resultToUpper]]];
    ];
  ];  (* End Module with detectedVar *)

  (* Combine segment data from both transports *)
  If[AssociationQ[resultToLower] || AssociationQ[resultToUpper],
    savedDataLower = If[AssociationQ[resultToLower],
      resultToLower["SegmentData"],
      {}
    ];
    savedDataUpper = If[AssociationQ[resultToUpper],
      resultToUpper["SegmentData"],
      {}
    ];

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
      "FixedParamValue" -> fixedVal,
      "EpsPrefactors" -> epsPrefactors,
      "SnapValues" -> snapValues
    |>
  ,
    Print["Error: Transport failed."];
    If[verbosity >= 2,
      Print["  Lower result head: ", Head[resultToLower]];
      Print["  Upper result head: ", Head[resultToUpper]];
    ];
    $Failed
  ]
];


(* ============================================================ *)
(* Integration with Feynman Trick Prefactors                     *)
(* ============================================================ *)

(*
  IntegrateCombinedMasters: Integrates the Feynman trick recursion.

  Gamma(v1+v2)/(Gamma(v1)*Gamma(v2)) *
    Integral[x^(v1-1)*(1-x)^(v2-1) * Sum_j c_j(x)*f_j(x), {x,0,1}]

  CRITICAL: The integrand Sum_j c_j(x)*f_j(x) must be combined BEFORE
  integration. Individual c_j(x) may have poles at x=0 and x=1, but the
  sum cancels these poles. Integrating term by term would diverge.

  The IBP coefficients c_j(x,d) depend on d. We expand d using
  FTConfiguration["DimensionExpression"]:
  c_j(x,eps) = Sum_k eps^k * c_j^{(k)}(x)

  The transport gives J_j = eps^{k_j} * I_j (prefactored masters), while
  IBP coefficients relate to I_j. So the full integrand at eps order n is:
  Sum_j Sum_{k=0}^{n+k_j} c_j^{(k)}(x) * J_j^{(n+k_j-k)}(x)

  where k_j is the eps-prefactor for master j.

  Arguments:
  - transportResult: transport data with SegmentData
  - ibpCoeffs: list of IBP coefficients (one per master, may contain d)
  - v1, v2: propagator exponents being combined
  - epsOrder: number of eps orders
  - epsPrefactors: list of eps-prefactor powers {k_1,...,k_n} (0 = no prefactor)
  - feynmanParam: the Feynman parameter symbol for this level

  Returns: list of eps-order coefficients for the integrated result.
*)
IntegrateCombinedMasters[transportResult_Association, ibpCoeffs_List,
    v1_Integer, v2_Integer, epsOrder_Integer,
    epsPrefactors_List:{}, feynmanParam_:Automatic, returnLaurent_:False] :=
Module[{gammaPrefactor, dimVar, epsSymbol, ibpCoeffOrders,
        segData, actualVar, combinedSegments, combinedData,
        prefactorSpec, result, prefacs, numMasters, workingMaxPower,
        combinedMinPower, combinedMaxPower, activeMasters},

  (* Gamma prefactor: Gamma(v1+v2)/(Gamma(v1)*Gamma(v2)) *)
  gammaPrefactor = Gamma[v1 + v2] / (Gamma[v1] * Gamma[v2]);

  numMasters = Length[ibpCoeffs];

  (* Default: no eps-prefactors (all zero) *)
  prefacs = If[Length[epsPrefactors] == numMasters, epsPrefactors,
    Table[0, {numMasters}]];
  workingMaxPower = Lookup[transportResult, "EpsilonOrder", epsOrder + Max[prefacs]];

  (* Expand each IBP coefficient in eps using the configured dimension. *)
  dimVar = FeynmanTrick`Private`$FTConfig["DimensionVariable"];
  epsSymbol = FeynmanTrick`Private`$FTConfig["EpsilonSymbol"];

  (* ibpCoeffOrders[[j]] is a Laurent coefficient association for master j.
     Each coefficient is a rational function of x (no eps). *)
  ibpCoeffOrders = Table[
    ExpandIBPCoeffLaurent[ibpCoeffs[[j]], workingMaxPower + Max[prefacs]],
    {j, numMasters}
  ];

  activeMasters = Select[Range[numMasters], ibpCoeffs[[#]] =!= 0 &];
  combinedMinPower = If[Length[activeMasters] == 0,
    0,
    Min[Table[ibpCoeffOrders[[j]]["MinPower"] - prefacs[[j]], {j, activeMasters}]]
  ];
  combinedMaxPower = workingMaxPower;

  (* Determine the Feynman parameter variable.
     DiffExp segments have format: {CurrLine, lineRelation, mainBounds, localBounds, seriesData}
     - seg[[1]] (CurrLine) is an Association <|xx3 -> a + b*x|> mapping local x to main var
     - seg[[2]] (lineRelation) is a Rule x -> expr for re-parameterization
     The IBP coefficients contain the Feynman parameter (e.g. Global`xx3).
     We use feynmanParam (from caller) as the substitution target, and
     seg[[1]] to get the expression for the main variable in local coords. *)
  actualVar = If[feynmanParam =!= Automatic,
    feynmanParam,
    (* Fallback: detect from CurrLine (segment element 1) *)
    Module[{seg1 = transportResult["SegmentData"][[1]]},
      If[AssociationQ[seg1[[1]]],
        First[Keys[seg1[[1]]]],
        Global`xx  (* ultimate fallback *)
      ]
    ]
  ];

  (* Also check if IBP coefficients use a different context for the same variable name.
     E.g., IBP might use Global`xx3 but segments use FeynmanTrick`...`xx3 *)
  Module[{ibpVar, varName, allSymbols, ibpSymbols},
    varName = SymbolName[actualVar];
    allSymbols = Cases[ibpCoeffOrders, _Symbol, Infinity] // DeleteDuplicates;
    ibpSymbols = Select[allSymbols, SymbolName[#] === varName &];
    If[Length[ibpSymbols] > 0 && !MemberQ[ibpSymbols, actualVar],
      (* IBP uses a different context - remap to match our actualVar *)
      ibpVar = First[ibpSymbols];
      If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 2,
        Print["    Variable context mismatch: IBP uses ", Context[ibpVar],
              varName, ", remapping to ", Context[actualVar], varName];
      ];
      ibpCoeffOrders = ibpCoeffOrders /. ibpVar -> actualVar;
    ];
  ];

  segData = transportResult["SegmentData"];

  (* A combined output order n is only complete when every active master
     still has transport data at the shifted index n + k_j - k for every
     nonzero IBP coefficient order k.  Above that bound contributions would
     be silently dropped, so trim those incomplete top orders away instead
     of handing partially combined epsilon orders downstream (they would
     poison the residual endpoint-sector resummation in the regularized
     integration). *)
  Module[{numEpsOrdersAvailable, completeMaxPower},
    numEpsOrdersAvailable = Length[uncompressSeriesData[segData[[1, 5]]][[1]]];
    completeMaxPower = If[Length[activeMasters] == 0,
      combinedMaxPower,
      Min[Table[
        numEpsOrdersAvailable - 1 - prefacs[[j]] +
          ibpCoeffOrders[[j]]["MinPower"],
        {j, activeMasters}
      ]]
    ];
    If[completeMaxPower < combinedMaxPower,
      If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 2,
        Print["    [Combine] trimming incomplete top eps orders ",
          completeMaxPower + 1, " .. ", combinedMaxPower,
          " (transport provides ", numEpsOrdersAvailable, " orders)"];
      ];
      (* Keep at least one output order so callers receive a well-formed
         (if explicitly warned-about) result. *)
      combinedMaxPower = Max[completeMaxPower, combinedMinPower];
    ];
    If[combinedMaxPower < epsOrder,
      Print["  Warning: combined series only reaches eps^", combinedMaxPower,
        " but eps^", epsOrder, " was requested. Increase the transport ",
        "epsilon order (IntegrationPoleAllowance) for trustworthy results."];
    ];
  ];

  (* For each segment, combine all masters weighted by IBP coefficients.
     The combined series at eps order n is:
       combined^{(n)}(x) = Sum_j Sum_{k=0}^{n+k_j} c_j^{(k)}(x) * J_j^{(n+k_j-k)}(x)

     where k_j is the eps-prefactor for master j (J_j = eps^{k_j} * I_j).
     The prefactor shift accounts for the fact that J_j starts at eps^0
     while I_j starts at eps^{-k_j}.

     The IBP coefficient multiplied by eps^{-k_j} shifts the eps order:
     c_j(x) * I_j(x) = c_j(x) * eps^{-k_j} * J_j(x)
     At eps order n: Sum_{k=0}^{n+k_j} c_j^{(k)}(x) * J_j^{(n+k_j-k)}(x)
  *)
	  combinedSegments = Table[
	    Module[{seg, seriesRaw, uncompressed, xLocal, xMainExpr,
	            combinedEpsOrders, numEpsOrders, snappedLine},
      seg = segData[[segIdx]];
      seriesRaw = seg[[5]];

	      (* Uncompress the saved local series. Keep DiffExp theta symbols
	         until point evaluation/integration so segments that cross local
	         x=0 can use the correct branch on each side. *)
      uncompressed = uncompressSeriesData[seriesRaw];

	      (* uncompressed[[j]] = list of eps orders for master j
	         uncompressed[[j]][[n+1]] = SeriesData for J_j at eps order n *)
      numEpsOrders = Length[uncompressed[[1]]];

      (* Get the coordinate transformation: x_main = f(x_local)
         CurrLine (seg[[1]]) is <|xx3 -> a + b*x|> giving the main variable
         as a function of the local DiffExp variable x. This is what we need
         to substitute into the IBP coefficients c_j(xx3). *)
      xLocal = DiffExp`Symbols`x;
      xMainExpr = If[AssociationQ[seg[[1]]],
        First[Values[seg[[1]]]],
        If[Head[seg[[1]]] === Rule, seg[[1, 2]], seg[[1]]]
	      ];
	      xMainExpr = DiffExp`Utilities`PChop[Expand[xMainExpr]];
	      xMainExpr = snapMainExpression[
	        xMainExpr, seg[[4]], Lookup[transportResult, "SnapValues", {0, 1}]
	      ];
	      snappedLine = If[AssociationQ[seg[[1]]],
	        Association[First[Keys[seg[[1]]]] -> xMainExpr],
	        seg[[1]]
	      ];

      (* Build combined series for each output eps order *)
      combinedEpsOrders = Table[
        Module[{total = 0},
          Do[
            If[ibpCoeffs[[j]] =!= 0,
              Module[{kj, cjSeries},
                kj = prefacs[[j]];

                (* For this master, the contribution at output order n is:
                   Sum_{k=0}^{n+kj} c_j^{(k)}(x) * J_j^{(n+kj-k)}(x)
                   where the J index is in transport's eps numbering *)
                Do[
                  If[LaurentCoeff[ibpCoeffOrders[[j]], k] =!= 0,
                    Module[{jIdx, cLocal, cSeries, prod},
                      jIdx = n + kj - k + 1;  (* 1-based index into eps orders *)

                      If[jIdx >= 1 && jIdx <= numEpsOrders,
                        (* Convert IBP coefficient to local coordinates *)
                        cLocal = Together[
                          DiffExp`Utilities`PChop[
                            Expand[LaurentCoeff[ibpCoeffOrders[[j]], k] /. actualVar -> xMainExpr]
                          ]
                        ];

                        (* Series-expand the coefficient in local coords *)
                        Module[{transportVal = uncompressed[[j, jIdx]]},
                          If[MatchQ[transportVal, _SeriesData],
                            (* Normal case: transport gives a series *)
                            Module[{expOrd},
                              expOrd = transportVal[[5]] - transportVal[[4]];
                              cSeries = Quiet[
                                Series[cLocal, {xLocal, 0, expOrd}],
                                {Power::infy, Infinity::indet, General::indet}
                              ];
                              prod = cSeries * transportVal;
                              total = total + prod;
                            ];
                          , If[!TrueQ[transportVal === 0 || PossibleZeroQ[transportVal]],
                              (* Transport can also give symbolic local expressions,
                                 such as const + Logx.  These terms carry the
                                 epsilon-dependent endpoint exponents and must be
                                 combined with the IBP coefficient rather than
                                 discarded as non-numeric. *)
                              Module[{expOrd = 30},
                                cSeries = Quiet[
                                  Series[cLocal * transportVal, {xLocal, 0, expOrd}],
                                  {Power::infy, Infinity::indet, General::indet}
                                ];
                                total = total + cSeries;
                              ];
                            ];
                          ];
                        ];
                      ];
                    ];
                  ];
                , {k, ibpCoeffOrders[[j]]["MinPower"], LaurentMaxPower[ibpCoeffOrders[[j]]]}];
              ];
            ];
          , {j, numMasters}];
          If[MatchQ[total, _SeriesData],
            Module[{seriesOrder, normal},
              seriesOrder = Max[0, Ceiling[total[[5]] / total[[6]]]];
              normal = DiffExp`Utilities`PChop[Expand[Normal[total]]];
              If[normal === 0,
                0,
                Quiet[Series[normal, {xLocal, 0, seriesOrder}]]
              ]
            ],
            DiffExp`Utilities`PChop[Expand[total]]
          ]
        ],
        {n, combinedMinPower, combinedMaxPower}
      ];

      (* Package as single-master segment:
         {line, transform, mainBounds, localBounds, {{epsOrder0, epsOrder1, ...}}} *)
	      {snappedLine, seg[[2]], seg[[3]], seg[[4]], {combinedEpsOrders}}
	    ],
    {segIdx, Length[segData]}
  ];

  combinedData = <|
    "SegmentData" -> combinedSegments,
    "NumIntegrals" -> 1,
    "EpsilonOrder" -> combinedMaxPower,
    "EpsilonMinPower" -> combinedMinPower
  |>;

  If[returnLaurent === "CombinedData",
    Return[combinedData, Module]
  ];

	  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 3,
    Do[
      Module[{orders = combinedSegments[[s, 5, 1]], powers, bad},
        powers = Table[
          If[MatchQ[orders[[n]], _SeriesData],
            orders[[n]][[4]]/orders[[n]][[6]],
            If[orders[[n]] === 0, Infinity, 0]
          ],
          {n, Length[orders]}
        ];
        bad = !FreeQ[orders, Indeterminate | ComplexInfinity | DirectedInfinity];
        Print["    [Combine] segment ", s, " bounds=", combinedSegments[[s, 3]],
              " local=", combinedSegments[[s, 4]],
              " map=", combinedSegments[[s, 1]],
              " leading powers=", powers, " bad=", bad];
        Module[{decomp = DiffExp`SingularityDecomposition`DecomposeSingularity[orders]},
          Print["      sectors=",
            ({#["a"], #["b"]} & /@ decomp)
          ];
        ];
        If[MemberQ[powers, _?(NumericQ[#] && # < 0 &)],
          Module[{minPow = Min[Select[powers, NumericQ[#] && # < 0 &]], coeffsAtMin},
            coeffsAtMin = Table[
              If[MatchQ[orders[[n]], _SeriesData],
                Quiet[SeriesCoefficient[orders[[n]], {DiffExp`Symbols`x, 0, minPow}]],
                0
              ],
              {n, 1, Min[3, Length[orders]]}
            ];
            Print["      coeffs at x^", minPow, " first orders=", coeffsAtMin];
          ];
        ];
      ],
      {s, Length[combinedSegments]}
    ];
  ];

  (* Check if all combined series are zero (no contributions from any master).
     This happens when all non-zero IBP coefficients multiply zero transport data.
     In this case, the boundary value is zero - skip DefiniteIntegralWithPrefactor. *)
  Module[{allZero},
    allZero = And @@ Table[
      And @@ Table[
        combinedSegments[[s, 5, 1, n]] === 0,
        {n, Length[combinedSegments[[s, 5, 1]]]}
      ],
      {s, Length[combinedSegments]}
    ];
    If[allZero,
      If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 2,
        Print["    [Combine] All combined series are zero -> returning zero boundary"];
      ];
      If[TrueQ[returnLaurent],
        Return[LaurentZero[0, epsOrder], Module],
        Return[Table[0, {epsOrder + 1}], Module]
      ];
    ];
  ];

  (* Now integrate the combined series with just the power-law prefactors.
     No rational factor needed - it's already been folded into the series. *)
  prefactorSpec = <|
    "PowerAtLower" -> v1 - 1,
    "PowerAtUpper" -> v2 - 1,
    "RationalFactor" -> 1,
    "Variable" -> actualVar
  |>;

  result = DiffExp`RegularizedIntegration`DefiniteIntegralWithPrefactorLaurent[
    combinedData, {0, 1}, prefactorSpec, combinedMinPower
  ];

  If[ListQ[result] && Length[result] >= 1 && AssociationQ[result[[1]]],
    Module[{scaled = LaurentScale[gammaPrefactor, result[[1]]]},
      If[TrueQ[returnLaurent],
        scaled,
        LaurentToNonNegativeList[scaled, epsOrder]
      ]
    ]
  ,
    If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
      Print["Warning: Combined integration returned unexpected format: ", Head[result]];
    ];
    If[TrueQ[returnLaurent], LaurentZero[0, epsOrder], Table[0, {epsOrder + 1}]]
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
    boundary_Integer, epsOrder_Integer, epsPrefactors_List:{}, returnLaurent_:False] :=
Module[{segData, targetSeg, uncompressed, numMasters,
        limitValues, seriesAtMaster, decomposition,
        taylorPart, limitVal, dimVar, epsSymbol, ibpCoeffsExpanded,
        prefacs, maxPrefactor, workingMaxPower, tol, zeroQ, nonNegativeQ,
        xLocal, targetInfo, localEndpoint, actualVar, coeffBoundaryRules,
        boundaryPrecision},

  segData = transportResult["SegmentData"];
  numMasters = transportResult["NumIntegrals"];
  prefacs = If[Length[epsPrefactors] == numMasters,
    epsPrefactors,
    Table[0, {numMasters}]
  ];
  maxPrefactor = Max[prefacs];
  workingMaxPower = Lookup[transportResult, "EpsilonOrder", epsOrder + maxPrefactor];
  tol = DiffExp`State`FEC[RationalizationTolerance];
  zeroQ[z_] := TrueQ[PossibleZeroQ[z]] ||
    TrueQ[NumericQ[z] && Abs[numericAtActivePrecision[z]] < tol];
  nonNegativeQ[z_] := TrueQ[z >= 0] ||
    zeroQ[z] ||
    TrueQ[NumericQ[z] && realNumericAtActivePrecision[z, tol] > -tol];

  xLocal = DiffExp`Symbols`x;

  (* Select the segment whose actual Feynman parameter reaches the endpoint.
     Segment element 3 stores the path parameter bounds, not necessarily the
     Feynman parameter bounds; using it can confuse the upper endpoint with the
     reversed lower transport. *)
  targetInfo = First[SortBy[
    Table[
      Module[{seg = segData[[segIdx]], xMainExpr, localBounds, values,
              distances, best},
        xMainExpr = If[AssociationQ[seg[[1]]],
          First[Values[seg[[1]]]],
          If[Head[seg[[1]]] === Rule, seg[[1, 2]], seg[[1]]]
        ];
        localBounds = seg[[4]];
        values = numericAtActivePrecision[(xMainExpr /. xLocal -> #) & /@ localBounds];
        distances = Abs[values - boundary];
        best = First[Ordering[distances, 1]];
        <|
          "Segment" -> seg,
          "LocalEndpoint" -> localBounds[[best]],
          "Distance" -> distances[[best]],
          "EndpointValue" -> values[[best]]
        |>
      ],
      {segIdx, Length[segData]}
    ],
    #["Distance"] &
  ]];

  targetSeg = targetInfo["Segment"];
  localEndpoint = targetInfo["LocalEndpoint"];
  actualVar = If[AssociationQ[targetSeg[[1]]],
    First[Keys[targetSeg[[1]]]],
    Automatic
  ];
  boundaryPrecision = If[NumericQ[targetInfo["EndpointValue"]],
    Precision[targetInfo["EndpointValue"]],
    activeNumericPrecision[]
  ];

  (* Expand IBP coefficients in eps using the configured dimension. *)
  dimVar = FeynmanTrick`Private`$FTConfig["DimensionVariable"];
  epsSymbol = FeynmanTrick`Private`$FTConfig["EpsilonSymbol"];
  ibpCoeffsExpanded = Table[
    ExpandIBPCoeffLaurent[ibpCoeffs[[i]], workingMaxPower + maxPrefactor],
    {i, Length[ibpCoeffs]}
  ];

  (* FIRE reductions and DiffExp segments may use symbols with the same name in
     different contexts. Align IBP coefficients with the segment variable before
     taking endpoint limits. *)
  coeffBoundaryRules = {};
  If[actualVar =!= Automatic,
    Module[{varName, allSymbols, ibpVars},
      varName = SymbolName[actualVar];
      allSymbols = Cases[ibpCoeffsExpanded, _Symbol, Infinity] // DeleteDuplicates;
      ibpVars = Select[allSymbols, SymbolName[#] === varName && # =!= actualVar &];
      coeffBoundaryRules = Thread[ibpVars -> actualVar];
      If[Length[coeffBoundaryRules] > 0,
        ibpCoeffsExpanded = ibpCoeffsExpanded /. coeffBoundaryRules;
      ];
    ];
  ];

  (* Uncompress series data. Theta symbols must be resolved at the actual
     local endpoint, not once for the whole segment. *)
  uncompressed = uncompressSeriesData[targetSeg[[5]]];

  (* For each master, evaluate the limit *)
  limitValues = LaurentZero[0, workingMaxPower];

  Do[
    If[ibpCoeffs[[mIdx]] =!= 0,
      (* Get series for this master (list of eps orders) *)
      seriesAtMaster = uncompressed[[mIdx]];

      limitVal = Table[0, {workingMaxPower + 1}];

      If[zeroQ[localEndpoint],
        (* Endpoint is the local expansion center. Decompose into
           x^{a+b eps} sectors and keep only finite pure Taylor terms. *)
        decomposition = DiffExp`SingularityDecomposition`DecomposeSingularity[seriesAtMaster];

        taylorPart = Select[
          decomposition,
          (nonNegativeQ[#["a"]] && zeroQ[#["b"]]) &
        ];

        Do[
          Module[{gSeries, constCoeff},
            gSeries = term["g"];
            (* g is a list of SeriesData objects, one per eps order *)
            Do[
              If[epsIdx <= Length[gSeries],
                Which[
                  MatchQ[gSeries[[epsIdx]], _SeriesData],
                    constCoeff = evaluateLocalExpressionAtPoint[
                      SeriesCoefficient[gSeries[[epsIdx]], 0] /.
                        DiffExp`Symbols`Logx -> 0,
                      localEndpoint,
                      localEndpointDirection[targetSeg, localEndpoint],
                      boundaryPrecision
                    ];
                    If[NumericQ[constCoeff], limitVal[[epsIdx]] += constCoeff],
                  NumericQ[gSeries[[epsIdx]]],
                    limitVal[[epsIdx]] += gSeries[[epsIdx]]
                ];
              ];
            , {epsIdx, Min[Length[limitVal], Length[gSeries]]}];
          ],
        {term, taylorPart}];
      ,
        (* Nonsingular endpoint inside the segment. Evaluate the saved local
           series at the endpoint's local coordinate instead of taking the
           constant term at the segment center. *)
	        Do[
	          If[epsIdx <= Length[seriesAtMaster],
	            Module[{seriesTerm = seriesAtMaster[[epsIdx]], endpointValue},
	              endpointValue = Which[
	                MatchQ[seriesTerm, _SeriesData],
	                  Quiet[Check[
	                    evaluateLocalExpressionAtPoint[
	                      Normal[seriesTerm],
	                      localEndpoint,
	                      localEndpointDirection[targetSeg, localEndpoint],
	                      boundaryPrecision
	                    ],
	                    0
	                  ]],
	                NumericQ[seriesTerm],
	                  seriesTerm,
	                seriesTerm === 0,
	                  0,
	                True,
	                  Quiet[Check[
	                    evaluateLocalExpressionAtPoint[
	                      seriesTerm,
	                      localEndpoint,
	                      localEndpointDirection[targetSeg, localEndpoint],
	                      boundaryPrecision
	                    ],
	                    0
	                  ]]
	              ];
              If[NumericQ[endpointValue],
                limitVal[[epsIdx]] += endpointValue;
              ];
            ];
          ];
        , {epsIdx, Min[Length[limitVal], Length[seriesAtMaster]]}];
      ];

      (* Add weighted contribution using eps-expanded coefficients *)
      (* ibpCoeffsExpanded[[mIdx]] is a Laurent association in eps. Evaluate
         its x-dependent coefficients at the same boundary. *)
      Module[{coeffsAtM, coeffsAtBoundary, contribution},
        coeffsAtM = ibpCoeffsExpanded[[mIdx]];
        coeffsAtBoundary = If[actualVar === Automatic,
          coeffsAtM,
          LaurentTrim[<|
            "MinPower" -> coeffsAtM["MinPower"],
            "Coefficients" -> Table[
              Module[{c = coeffsAtM["Coefficients"][[cIdx]], endpointCoeff},
                endpointCoeff = Quiet[
                  Check[
                    Limit[c, actualVar -> boundary],
                    c /. actualVar -> SetPrecision[boundary, boundaryPrecision]
                  ],
                  {Power::infy, Infinity::indet, General::indet}
                ];
                If[FreeQ[endpointCoeff, Indeterminate | ComplexInfinity | DirectedInfinity],
                  endpointCoeff,
                  Quiet[Check[
                    c /. actualVar -> SetPrecision[boundary, boundaryPrecision],
                    0
                  ]]
                ]
              ],
              {cIdx, Length[coeffsAtM["Coefficients"]]}
            ]
          |>]
        ];
        contribution = MultiplyLaurentShifted[
          coeffsAtBoundary,
          <|"MinPower" -> 0, "Coefficients" -> limitVal|>,
          prefacs[[mIdx]],
          coeffsAtBoundary["MinPower"] - prefacs[[mIdx]],
          workingMaxPower
        ];
        limitValues = LaurentAdd[limitValues, contribution];
      ];
    ];
  , {mIdx, numMasters}];

  If[TrueQ[returnLaurent],
    LaurentTrim[limitValues],
    LaurentToNonNegativeList[limitValues, epsOrder]
  ]
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

(* Helper: expand IBP coefficient in eps, returning list of eps-order coefficients *)
ExpandIBPCoeffInEps[coeff_, epsOrder_Integer] :=
Module[{dimVar, epsSymbol, expanded},
  dimVar = FeynmanTrick`Private`$FTConfig["DimensionVariable"];
  epsSymbol = FeynmanTrick`Private`$FTConfig["EpsilonSymbol"];
  expanded = coeff /. dimVar -> FeynmanTrick`Private`DimensionExpression[];
  Table[SeriesCoefficient[expanded + O[epsSymbol]^(epsOrder + 1), k], {k, 0, epsOrder}]
];

LaurentMaxPower[laur_Association] :=
  laur["MinPower"] + Length[laur["Coefficients"]] - 1;

LaurentCoeff[laur_Association, power_Integer] := Module[
  {idx = power - laur["MinPower"] + 1},
  If[idx >= 1 && idx <= Length[laur["Coefficients"]],
    laur["Coefficients"][[idx]],
    0
  ]
];

LaurentZero[minPower_Integer, maxPower_Integer] := <|
  "MinPower" -> minPower,
  "Coefficients" -> Table[0, {Max[0, maxPower - minPower + 1]}]
|>;

LaurentAdd[a_Association, b_Association] := Module[
  {minPower, maxPower},
  minPower = Min[a["MinPower"], b["MinPower"]];
  maxPower = Max[LaurentMaxPower[a], LaurentMaxPower[b]];
  <|
    "MinPower" -> minPower,
    "Coefficients" -> Table[
      Together[Expand[LaurentCoeff[a, p] + LaurentCoeff[b, p]]],
      {p, minPower, maxPower}
    ]
  |>
];

LaurentScale[c_, laur_Association] := <|
  "MinPower" -> laur["MinPower"],
  "Coefficients" -> (Together[Expand[c * #]] & /@ laur["Coefficients"])
|>;

zeroCoeffQ[c_] := TrueQ[PossibleZeroQ[c]] ||
  TrueQ[NumericQ[c] && Abs[numericAtActivePrecision[c]] < 10^-40];

LaurentTrim[laur_Association] := Module[
  {minPower = laur["MinPower"], coeffs = laur["Coefficients"]},
  While[Length[coeffs] > 0 && zeroCoeffQ[First[coeffs]],
    coeffs = Rest[coeffs];
    minPower++;
  ];
  If[Length[coeffs] == 0,
    <|"MinPower" -> 0, "Coefficients" -> {0}|>,
    <|"MinPower" -> minPower, "Coefficients" -> coeffs|>
  ]
];

LaurentToRange[laur_Association, minPower_Integer, maxPower_Integer] :=
  Table[LaurentCoeff[laur, p], {p, minPower, maxPower}];

LaurentToNonNegativeList[laur_Association, epsOrder_Integer] :=
  LaurentToRange[laur, 0, epsOrder];

ExpandIBPCoeffLaurent[coeff_, maxPower_Integer] :=
Module[{dimVar, epsSymbol, expanded, minPower, coeffs},
  dimVar = FeynmanTrick`Private`$FTConfig["DimensionVariable"];
  epsSymbol = FeynmanTrick`Private`$FTConfig["EpsilonSymbol"];

  If[coeff === 0,
    Return[LaurentZero[0, maxPower]]
  ];

  expanded = Quiet[
    Check[
      Series[coeff /. dimVar -> FeynmanTrick`Private`DimensionExpression[],
        {epsSymbol, 0, maxPower}],
      $Failed
    ]
  ];

  If[expanded === $Failed,
    expanded = coeff /. dimVar -> FeynmanTrick`Private`DimensionExpression[];
  ];

  minPower = If[Head[expanded] === SeriesData,
    expanded[[4]] / expanded[[6]],
    0
  ];

  If[!IntegerQ[minPower], minPower = Floor[minPower]];

  coeffs = Table[
    Together[
      Normal[
        If[Head[expanded] === SeriesData,
          SeriesCoefficient[expanded, p],
          SeriesCoefficient[expanded + O[epsSymbol]^(maxPower + 1), {epsSymbol, 0, p}]
        ]
      ]
    ],
    {p, minPower, maxPower}
  ];

  LaurentTrim[<|"MinPower" -> minPower, "Coefficients" -> coeffs|>]
];

MultiplyLaurentShifted[coeffs_Association, values_Association,
    shift_Integer, minOut_Integer, maxOut_Integer] :=
Module[{coeffMin, coeffMax},
  coeffMin = coeffs["MinPower"];
  coeffMax = LaurentMaxPower[coeffs];
  <|
    "MinPower" -> minOut,
    "Coefficients" -> Table[
      Sum[
        LaurentCoeff[coeffs, k] * LaurentCoeff[values, n + shift - k],
        {k, coeffMin, coeffMax}
      ] // Together // Expand,
      {n, minOut, maxOut}
    ]
  |> // LaurentTrim
];

ShiftRawBoundariesToFinite[rawBCs_List, epsOrder_Integer] := Module[
  {trimmed, minPower, maxPower, shift, finiteBCs},
  If[TrueQ[DiffExp`RegularizedIntegration`Private`$DebugBadRegularizedIntegration],
    Print["DEBUG_SHIFT_INPUT_BAD=", Position[
      rawBCs,
      Indeterminate | ComplexInfinity | DirectedInfinity[_],
      Infinity
    ]];
  ];
  trimmed = LaurentTrim /@ rawBCs;
  If[TrueQ[DiffExp`RegularizedIntegration`Private`$DebugBadRegularizedIntegration],
    Print["DEBUG_SHIFT_TRIMMED_BAD=", Position[
      trimmed,
      Indeterminate | ComplexInfinity | DirectedInfinity[_],
      Infinity
    ]];
  ];
  minPower = Min[trimmed[[All, "MinPower"]]];
  maxPower = Max[LaurentMaxPower /@ trimmed];
  shift = Max[0, -minPower];

  finiteBCs = Table[
    Table[
      LaurentCoeff[trimmed[[i]], n - shift],
      {n, 0, Max[epsOrder + shift, maxPower + shift]}
    ],
    {i, Length[trimmed]}
  ];
  If[TrueQ[DiffExp`RegularizedIntegration`Private`$DebugBadRegularizedIntegration],
    Print["DEBUG_SHIFT_FINITE_BAD=", Position[
      finiteBCs,
      Indeterminate | ComplexInfinity | DirectedInfinity[_],
      Infinity
    ]];
  ];

  <|
    "BoundaryValues" -> finiteBCs,
    "EpsPrefactors" -> Table[shift, {Length[trimmed]}],
    "RawBoundaryValues" -> trimmed,
    "RawMinPower" -> minPower,
    "RawMaxPower" -> maxPower
  |>
];

BoundaryRequestRecords[mastersAtLevel_List, combinedPositions_List] :=
Module[{posI, posJ},
  {posI, posJ} = combinedPositions;
  Table[
    Module[{masterVec, vi, vj, neededVec, case},
      masterVec = mastersAtLevel[[masterIdx]];
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
          neededVec[[posJ]] = 0;
      ];

      <|
        "MasterIndex" -> masterIdx,
        "MasterVec" -> masterVec,
        "Vi" -> vi,
        "Vj" -> vj,
        "Case" -> case,
        "NeededVec" -> neededVec
      |>
    ],
    {masterIdx, Length[mastersAtLevel]}
  ]
];

CollectLevelIBPSingularFactors[ftData_Association, level_Integer] :=
Module[{levelData, levelAbove, mastersAtLevel, mastersAbove,
        combinedPositions, posI, posJ, topologyAbove, feynmanParamAbove,
        varName, factors = {}, boundaryRequests, neededVecs, reductions},

  If[level <= 0 || !KeyExistsQ[ftData["Levels"], level] ||
     !KeyExistsQ[ftData["Levels"], level - 1],
    Return[{}]
  ];

  levelData = ftData["Levels"][level - 1];
  levelAbove = ftData["Levels"][level];
  mastersAtLevel = levelData["Masters"];
  mastersAbove = levelAbove["Masters"];
  combinedPositions = levelAbove["CombinedPositions"];
  {posI, posJ} = combinedPositions;
  topologyAbove = levelAbove["Topology"];
  feynmanParamAbove = levelAbove["FeynmanParameter"];
  varName = SymbolName[feynmanParamAbove];
  boundaryRequests = BoundaryRequestRecords[mastersAtLevel, combinedPositions];
  neededVecs = DeleteDuplicates[#["NeededVec"] & /@ boundaryRequests];
  reductions = If[neededVecs === {},
    <||>,
    FeynmanTrick`FIREInterface`ReduceIntegrals[topologyAbove, neededVecs]
  ];
  If[reductions === $Failed, Return[{}]];

  Do[
    Module[{request, neededVec, expr, ibpCoeffs},
      request = boundaryRequests[[masterIdx]];
      neededVec = request["NeededVec"];
      If[!KeyExistsQ[reductions, neededVec], Continue[]];
      expr = reductions[neededVec];
      ibpCoeffs = Table[
        Coefficient[expr, Global`G[1, mastersAbove[[j]]]],
        {j, Length[mastersAbove]}
      ];

      factors = Join[
        factors,
        Flatten[
          Table[
            Module[{den, factorList},
              den = Denominator[Together[coeff]];
              factorList = If[den === 1,
                {},
                FactorList[Factor[den]][[All, 1]]
              ];
              Select[factorList,
                !FreeQ[#, s_Symbol /; SymbolName[s] === varName] &
              ]
            ],
            {coeff, ibpCoeffs}
          ]
        ]
      ];
    ],
    {masterIdx, Length[mastersAtLevel]}
  ];

  DeleteDuplicates[
    DeleteCases[Factor /@ factors, 0 | 1 | -1],
    TrueQ[PossibleZeroQ[Expand[#1 - #2]]] ||
      TrueQ[PossibleZeroQ[Expand[#1 + #2]]] &
  ]
];

RequiredTransportEpsilonOrder[ftData_Association, level_Integer,
    epsOrder_Integer, epsPrefactors_List:{}] :=
Module[{levelData, levelBelow, mastersBelow, mastersAbove,
        combinedPositions, posI, posJ, topologyAbove, prefacs,
        required = epsOrder, maxProbe, boundaryRequests, neededVecs,
        reductions, integrationPoleAllowance},

  If[level <= 0 || !KeyExistsQ[ftData["Levels"], level] ||
     !KeyExistsQ[ftData["Levels"], level - 1],
    Return[epsOrder + If[ListQ[epsPrefactors] && Length[epsPrefactors] > 0,
      Max[epsPrefactors],
      0
    ]]
  ];

  levelBelow = ftData["Levels"][level - 1];
  levelData = ftData["Levels"][level];
  mastersBelow = levelBelow["Masters"];
  mastersAbove = levelData["Masters"];
  combinedPositions = levelData["CombinedPositions"];
  {posI, posJ} = combinedPositions;
  topologyAbove = levelData["Topology"];
  prefacs = If[ListQ[epsPrefactors] && Length[epsPrefactors] == Length[mastersAbove],
    epsPrefactors,
    Table[0, {Length[mastersAbove]}]
  ];
  maxProbe = epsOrder + Max[prefacs] + 20;
  boundaryRequests = BoundaryRequestRecords[mastersBelow, combinedPositions];
  integrationPoleAllowance = If[
    AnyTrue[boundaryRequests, #["Case"] === "integrate" &],
    Module[{raw = Environment["FT_INTEGRATION_POLE_ALLOWANCE"], parsed, configured},
      parsed = If[StringQ[raw], Quiet[Check[ToExpression[raw], None]], None];
      configured = Lookup[
        FeynmanTrick`Private`$FTConfig, "IntegrationPoleAllowance", 4
      ];
      Which[
        IntegerQ[parsed] && parsed >= 0, parsed,
        IntegerQ[configured] && configured >= 0, configured,
        True, 4
      ]
    ],
    0
  ];
  required = Max[required, epsOrder + integrationPoleAllowance];
  neededVecs = DeleteDuplicates[#["NeededVec"] & /@ boundaryRequests];
  reductions = If[neededVecs === {},
    <||>,
    FeynmanTrick`FIREInterface`ReduceIntegrals[topologyAbove, neededVecs]
  ];
  If[reductions === $Failed, Return[Max[0, Ceiling[required]]]];

  Do[
    Module[{request, neededVec, expr, ibpCoeffs, coeffLaurent},
      request = boundaryRequests[[masterIdx]];
      neededVec = request["NeededVec"];
      If[!KeyExistsQ[reductions, neededVec], Continue[]];
      expr = reductions[neededVec];
      ibpCoeffs = Table[
        Coefficient[expr, Global`G[1, mastersAbove[[j]]]],
        {j, Length[mastersAbove]}
      ];

      Do[
        If[ibpCoeffs[[j]] =!= 0,
          coeffLaurent = ExpandIBPCoeffLaurent[ibpCoeffs[[j]], maxProbe];
          required = Max[
            required,
            epsOrder + prefacs[[j]] - coeffLaurent["MinPower"] +
              integrationPoleAllowance
          ];
        ],
        {j, Length[mastersAbove]}
      ];
    ],
    {masterIdx, Length[mastersBelow]}
  ];

  Max[0, Ceiling[required]]
];

(* Helper: multiply two eps-expanded coefficient lists (convolution) *)
  MultiplyEpsCoeffs[coeffs1_List, coeffs2_List, epsOrder_Integer] :=
Table[
  Sum[
    If[k >= 0 && k < Length[coeffs1] && (n - k) >= 0 && (n - k) < Length[coeffs2],
      coeffs1[[k + 1]] * coeffs2[[n - k + 1]],
      0
    ],
    {k, 0, n}
  ],
  {n, 0, epsOrder}
];

(* Multiply eps-expanded coeffs by a prefactored master series.
   If J = eps^shift I, then I = eps^-shift J.  The coefficient at
   output order n receives coeff[k] * J[n + shift - k]. *)
MultiplyEpsCoeffsShifted[coeffs1_List, coeffs2_List, shift_Integer, epsOrder_Integer] :=
Table[
  Sum[
    Module[{jIdx = n + shift - k},
      If[k >= 0 && k < Length[coeffs1] &&
         jIdx >= 0 && jIdx < Length[coeffs2],
        coeffs1[[k + 1]] * coeffs2[[jIdx + 1]],
        0
      ]
    ],
    {k, 0, n + shift}
  ],
  {n, 0, epsOrder}
];

ComputeLevelBoundary[ftData_Association, level_Integer,
    transportResult_Association, epsOrder_Integer] :=
Module[{levelData, levelAbove, mastersAtLevel, mastersAbove,
        combinedPositions, posI, posJ, topologyAbove, bcValues,
        feynmanParamAbove, epsPrefactorsAbove, shiftedBoundary,
        workingMaxPower, boundaryRequests, neededVecs, reductions,
        preparedTransportResult},

  levelData = ftData["Levels"][level];
  levelAbove = ftData["Levels"][level + 1];
  mastersAtLevel = levelData["Masters"];
  mastersAbove = levelAbove["Masters"];
  combinedPositions = levelAbove["CombinedPositions"];
  {posI, posJ} = combinedPositions;
  topologyAbove = levelAbove["Topology"];

  (* Get the Feynman parameter for the level above - needed for IBP coefficient expansion *)
  feynmanParamAbove = levelAbove["FeynmanParameter"];

  (* Transported series may be in a prefactored basis J_i = eps^k_i I_i.
     Prefer explicit metadata from TransportLevel/caller, with the old matrix-pole
     inference as a fallback for existing callers. *)
  epsPrefactorsAbove = Lookup[
    transportResult,
    "EpsPrefactorsAbove",
    Lookup[transportResult, "EpsPrefactors", Automatic]
  ];
  If[!(ListQ[epsPrefactorsAbove] && Length[epsPrefactorsAbove] == Length[mastersAbove]),
    Module[{diffMatAbove, epsSymbol},
      epsSymbol = FeynmanTrick`Private`$FTConfig["EpsilonSymbol"];
      diffMatAbove = levelAbove["DiffMatrix"];
      If[MatrixQ[diffMatAbove] &&
         FeynmanTrick`EpsPrefactors`CheckEpsPoles[diffMatAbove, epsSymbol],
        epsPrefactorsAbove = FeynmanTrick`EpsPrefactors`FindEpsPrefactors[diffMatAbove, epsSymbol],
        epsPrefactorsAbove = Table[0, {Length[mastersAbove]}]
      ];
    ];
  ];
  workingMaxPower = Lookup[
    transportResult,
    "EpsilonOrder",
    epsOrder + Max[epsPrefactorsAbove]
  ];
  boundaryRequests = BoundaryRequestRecords[mastersAtLevel, combinedPositions];
  neededVecs = DeleteDuplicates[#["NeededVec"] & /@ boundaryRequests];
  reductions = If[neededVecs === {},
    <||>,
    FeynmanTrick`FIREInterface`ReduceIntegrals[topologyAbove, neededVecs]
  ];
  If[reductions === $Failed,
    If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
      Print["  Warning: IBP batch reduction failed for level ", level];
    ];
    Return[$Failed];
  ];

  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["Computing boundary for level ", level, " from level ", level + 1];
    Print["  Masters at level ", level, ": ", Length[mastersAtLevel]];
    Print["  Masters at level ", level + 1, ": ", Length[mastersAbove]];
    Print["  Combined positions: {", posI, ", ", posJ, "}"];
  ];

  preparedTransportResult = Join[
    transportResult,
    <|"SegmentData" -> prepareTransportSegments[transportResult]|>
  ];

  (* For each master at the current level *)
  bcValues = Table[
    Module[{masterVec, vi, vj, neededVec, case, reduction, expr,
            ibpCoeffs, totalBC},

      Module[{request = boundaryRequests[[masterIdx]]},
        masterVec = request["MasterVec"];
        vi = request["Vi"];
        vj = request["Vj"];
        case = request["Case"];
        neededVec = request["NeededVec"];
      ];

      If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 2,
        Print["  Master ", masterIdx, " (", masterVec, "): case=", case,
              " vi=", vi, " vj=", vj];
      ];

      If[!KeyExistsQ[reductions, neededVec],
        If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
          Print["  Warning: IBP reduction missing for ", neededVec];
        ];
        Return[LaurentZero[0, workingMaxPower], Module];
      ];

      (* Extract the reduction expression *)
      expr = reductions[neededVec];

      (* Extract coefficient of each master G[1, masters_j] *)
      ibpCoeffs = Table[
        Coefficient[expr, Global`G[1, mastersAbove[[j]]]],
        {j, Length[mastersAbove]}
      ];

      If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 2,
        Print["    IBP coefficients: ", ibpCoeffs];
      ];

	      (* Compute the boundary value based on the case *)
	      totalBC = Switch[case,
	        "integrate",
          (* Full integration: combine Sum_j c_j(x)*f_j(x) first, then integrate.
             Individual c_j(x) may have poles at x=0,1 that cancel in the sum.
             Also accounts for eps-prefactors (J vs I basis). *)
          Module[{},
            If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 2,
              Print["    Eps prefactors for level above: ", epsPrefactorsAbove];
            ];

            IntegrateCombinedMasters[
              preparedTransportResult, ibpCoeffs, vi, vj, epsOrder,
              epsPrefactorsAbove, feynmanParamAbove, True
            ]
          ],

        "limitUpper",
          (* lim_{x->1} sum_j ibpCoeffs[[j]] * f_j(x) *)
          EvaluateLimitFromTransport[
            preparedTransportResult, ibpCoeffs, 1, epsOrder, epsPrefactorsAbove, True
          ],

        "limitLower",
          (* lim_{x->0} sum_j ibpCoeffs[[j]] * f_j(x) *)
          EvaluateLimitFromTransport[
            preparedTransportResult, ibpCoeffs, 0, epsOrder, epsPrefactorsAbove, True
          ],

        "direct",
          (* Evaluate at the fixed parameter value *)
          (* The "direct" case means both propagators are absent, *)
          (* so this integral at level+1 equals the one at level directly *)
          (* Evaluate at the fixed point from the boundary values *)
          Module[{directVal},
            directVal = LaurentZero[0, workingMaxPower];
            Do[
              If[ibpCoeffs[[j]] =!= 0,
                (* Use the boundary values from level+1 *)
                (* Expand IBP coeff in eps and convolve with boundary values *)
                Module[{coeffExpanded, bcAbove, shift, contribution},
                  shift = epsPrefactorsAbove[[j]];
                  coeffExpanded = ExpandIBPCoeffLaurent[
                    ibpCoeffs[[j]], workingMaxPower + shift
                  ];
                  bcAbove = transportResult["BoundaryValuesAbove"][[j]];
                  contribution = MultiplyLaurentShifted[
                    coeffExpanded,
                    <|"MinPower" -> 0, "Coefficients" -> bcAbove|>,
                    shift,
                    coeffExpanded["MinPower"] - shift,
                    workingMaxPower
                  ];
                  directVal = LaurentAdd[directVal, contribution];
                ];
              ];
            , {j, Length[mastersAbove]}];
            directVal
	          ]
		      ];
	        totalBC
		    ],
    {masterIdx, Length[mastersAtLevel]}
  ];

  shiftedBoundary = ShiftRawBoundariesToFinite[bcValues, epsOrder];

  <|
    "BoundaryValues" -> shiftedBoundary["BoundaryValues"],
    "Masters" -> mastersAtLevel,
    "Level" -> level,
    "EpsPrefactors" -> shiftedBoundary["EpsPrefactors"],
    "EpsPrefactorsAbove" -> epsPrefactorsAbove,
    "RawBoundaryValues" -> shiftedBoundary["RawBoundaryValues"],
    "RawMinPower" -> shiftedBoundary["RawMinPower"],
    "RawMaxPower" -> shiftedBoundary["RawMaxPower"]
  |>
];


(* ============================================================ *)
(* Full Integration Pipeline                                     *)
(* ============================================================ *)

RunIntegrationPipeline[ftData_Association, outputDir_String, epsOrder_Integer:4,
    opts:OptionsPattern[{
      "WorkingPrecision" -> 500,
      "ExpansionOrder" -> 50,
      "DivisionOrder" -> 4,
      "HomogeneousSolve" -> "DontExpand",
      "UseRationalRecurrence" -> True,
      "IntegrationStrategy" -> "Default",
      "EstimateError" -> "Fast",
      "CheckpointDirectory" -> None,
      "StopAfterBoundaryLevel" -> None
    }]] :=
  Module[{nLevels, currentBCs, currentPrefactors, matrixDir,
        transportResult, levelBoundary, precision, expOrder,
        updatedFtData, transportEpsOrder, boundaryWorkingOrder,
        finalBoundary, extraSingularFactors, homogeneousSolve,
        useRationalRecurrence, integrationStrategy, estimateError, checkpointDir,
        stopAfterBoundaryLevel, divisionOrder},

  precision = OptionValue["WorkingPrecision"];
  expOrder = OptionValue["ExpansionOrder"];
  divisionOrder = OptionValue["DivisionOrder"];
  homogeneousSolve = OptionValue["HomogeneousSolve"];
  useRationalRecurrence = OptionValue["UseRationalRecurrence"];
  integrationStrategy = OptionValue["IntegrationStrategy"];
  estimateError = OptionValue["EstimateError"];
  checkpointDir = OptionValue["CheckpointDirectory"];
  stopAfterBoundaryLevel = OptionValue["StopAfterBoundaryLevel"];
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
    boundaryWorkingOrder = epsOrder + nLevels;
    deepBoundary = FeynmanTrick`BoundaryConditions`DeepestLevelBoundary[
      updatedFtData, boundaryWorkingOrder
    ];

    If[deepBoundary === $Failed,
      Print["Error: DeepestLevelBoundary failed."];
      Return[$Failed];
    ];

    currentBCs = deepBoundary["BoundaryValues"];
    currentPrefactors = deepBoundary["EpsPrefactors"];
    transportEpsOrder = Length[First[currentBCs]] - 1;

    If[transportEpsOrder > epsOrder,
      Do[
        FeynmanTrick`FeynmanTrickIteration`ExportLevel[
          updatedFtData, level, outputDir, "diffexp", transportEpsOrder
        ],
        {level, nLevels, 1, -1}
      ];
    ];

    If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
      Print["  Deepest level boundary computed successfully."];
      Print["  Eps prefactors: ", currentPrefactors];
      Print["  Transport epsilon order: ", transportEpsOrder];
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
    Module[{requiredOrder = RequiredTransportEpsilonOrder[
        updatedFtData, level, epsOrder, currentPrefactors
      ]},
      transportEpsOrder = Min[Length[First[currentBCs]] - 1, requiredOrder];
      If[transportEpsOrder < requiredOrder,
        Print["  Warning: level ", level, " transport needs eps order ",
          requiredOrder, " but only ", transportEpsOrder,
          " is available from the boundary above. Top Laurent orders of ",
          "this level's boundary may be incomplete; increase the deepest ",
          "boundary order to recover them."];
      ];
    ];
    If[transportEpsOrder < Length[First[currentBCs]] - 1,
      currentBCs = currentBCs[[All, 1 ;; transportEpsOrder + 1]];
      If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
        Print["  Trimming transport epsilon order to ", transportEpsOrder,
              " for level ", level, "."];
      ];
    ];

    FeynmanTrick`FeynmanTrickIteration`ExportLevel[
      updatedFtData, level, outputDir, "diffexp", transportEpsOrder
    ];

    If[!DirectoryQ[matrixDir],
      Print["Error: Matrix directory not found: ", matrixDir];
      Return[$Failed];
    ];

    extraSingularFactors = CollectLevelIBPSingularFactors[updatedFtData, level];
    If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 2 &&
       Length[extraSingularFactors] > 0,
      Print["  Extra IBP segmentation factors: ", extraSingularFactors];
    ];

    (* Transport from fixed point to cover [0,1] *)
    transportResult = TransportLevel[
      matrixDir, currentBCs, transportEpsOrder,
      "WorkingPrecision" -> precision,
      "ExpansionOrder" -> expOrder,
      "DivisionOrder" -> divisionOrder,
      "Verbosity" -> FeynmanTrick`Private`$FTConfig["Verbosity"],
      "EpsPrefactors" -> currentPrefactors,
      "ExtraSingularFactors" -> extraSingularFactors,
      "HomogeneousSolve" -> homogeneousSolve,
      "UseRationalRecurrence" -> useRationalRecurrence,
      "IntegrationStrategy" -> integrationStrategy,
      "EstimateError" -> estimateError
    ];

    If[transportResult === $Failed,
      Print["Error: Transport failed at level ", level];
      Return[$Failed];
    ];

    (* Store boundary values used for transport (needed for "direct" case) *)
    transportResult["BoundaryValuesAbove"] = currentBCs;
    transportResult["EpsPrefactorsAbove"] = currentPrefactors;

    (* Compute boundary for level-1 by integration *)
    levelBoundary = ComputeLevelBoundary[
      updatedFtData, level - 1, transportResult, epsOrder
    ];

    If[AssociationQ[levelBoundary],
      currentBCs = levelBoundary["BoundaryValues"];
      currentPrefactors = Lookup[
        levelBoundary,
        "EpsPrefactors",
        Table[0, {Length[currentBCs]}]
      ];
      If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
        Print["  Level ", level - 1, " boundary computed: ",
              Length[currentBCs], " masters"];
      ];
      If[StringQ[checkpointDir],
        If[!DirectoryQ[checkpointDir],
          CreateDirectory[checkpointDir, CreateIntermediateDirectories -> True]
        ];
        Put[
          <|
            "BoundaryLevel" -> level - 1,
            "NextTransportLevel" -> level - 1,
            "BoundaryValues" -> currentBCs,
            "EpsPrefactors" -> currentPrefactors,
            "TransportEpsilonOrder" -> Length[First[currentBCs]] - 1,
            "RequestedEpsilonOrder" -> epsOrder,
            "OutputDirectory" -> outputDir,
            "FtData" -> updatedFtData,
            "DimensionExpression" ->
              FeynmanTrick`Private`DimensionExpression[]
          |>,
          FileNameJoin[{checkpointDir,
            "boundary_level_" <> ToString[level - 1] <> ".m"}]
        ];
      ];
      finalBoundary = levelBoundary;
      If[stopAfterBoundaryLevel === level - 1,
        If[level - 1 >= 1,
          FeynmanTrick`FeynmanTrickIteration`ExportLevel[
            updatedFtData, level - 1, outputDir, "diffexp",
            Length[First[currentBCs]] - 1
          ];
        ];
        Return[
          <|
            "BoundaryValues" -> currentBCs,
            "Level" -> level - 1,
            "EpsPrefactors" -> currentPrefactors,
            "TransportEpsilonOrder" -> Length[First[currentBCs]] - 1,
            "FtData" -> updatedFtData,
            "OutputDirectory" -> outputDir,
            "DimensionExpression" ->
              FeynmanTrick`Private`DimensionExpression[],
            "StoppedAfterBoundaryLevel" -> level - 1
          |>,
          Module
        ];
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
    "EpsPrefactors" -> currentPrefactors,
    "RawBoundaryValues" -> If[AssociationQ[finalBoundary],
      Lookup[finalBoundary, "RawBoundaryValues", Missing["NotAvailable"]],
      Missing["NotAvailable"]
    ],
    "RawMinPower" -> If[AssociationQ[finalBoundary],
      Lookup[finalBoundary, "RawMinPower", Missing["NotAvailable"]],
      Missing["NotAvailable"]
    ],
    "FtData" -> updatedFtData
  |>
];


End[];
EndPackage[];
