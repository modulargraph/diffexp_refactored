(*::Package::*)(*DiffExpIntegration - Bridge between FeynmanTrick and
                DiffExp *)(*Handles transport,
                           integration with Feynman trick prefactors,
                               *)(*and the full bottom -
                                  up integration pipeline.*)

    BeginPackage["FeynmanTrick`DiffExpIntegration`", {"FeynmanTrick`"}];

TransportLevel::usage =
    "TransportLevel[matrixDir, boundaryValues, epsOrder, opts] loads DiffExp with \
matrices from matrixDir, sets boundary conditions at the fixed parameter point, \
and transports to cover [0,1]. Returns the TransportTo result with SegmentData.";

IntegrateCombinedMasters::usage =
    "IntegrateCombinedMasters[transportResult, ibpCoeffs, v1, v2, epsOrder, epsPrefactors, feynmanParam] \
integrates the combined Feynman trick integrand: \
Gamma(v1+v2)/(Gamma(v1)*Gamma(v2)) * Integral[x^(v1-1)*(1-x)^(v2-1) * Sum_j c_j(x)*f_j(x), {x,0,1}]. \
Combines the integrand before integration to handle cancellation of poles in IBP coefficients. \
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
(* Resolve DiffExp analytic continuation symbols                 *)
(* ============================================================ *)

(*
  DiffExp's transport produces series containing symbolic theta functions
  (θp, θm) from analytic continuation, and Logx = Log[x].
  For integration over [0,1] with +i*delta prescription, we need to
  resolve these to numerical values:
    θp -> 1, θm -> 0  (approaching from upper half plane)
    Logx -> Log[x]     (actual logarithm)
*)
$thetaRules = {
  DiffExp`Symbols`\[Theta]p -> 1,
  DiffExp`Symbols`\[Theta]m -> 0
};

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
    DiffExp`State`Verbosity -> verbosity,
    DiffExp`State`SegmentationStrategy -> "Predivision"
  };

  If[verbosity >= 1,
    Print["  Loading DiffExp matrices from: ", matrixDir];
  ];

  DiffExp`LoadConfiguration[diffExpConfig];

  (* After loading, DiffExp auto-detects the variable from filenames (dxx_*.m).
     Extract the detected variable symbol so our points match DiffExp's internal state. *)
  Module[{detectedVar},
    detectedVar = First[DiffExp`State`FEC[System`Variables]];

    (* Add delta prescriptions for branch cuts at endpoints *)
    (* {polynomial, sign}: sign=1 means approach from above (Im > 0) *)
    (* Also disable abort on analytic continuation failure - the pipeline
       handles incomplete results gracefully *)
    DiffExp`UpdateConfiguration[{
      DiffExp`State`DeltaPrescriptions -> {
        {detectedVar, 1},        (* x + I*delta at x=0 *)
        {1 - detectedVar, 1}     (* (1-x) + I*delta at x=1 *)
      },
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
      Print["  Transporting from xx=", fixedVal, " towards 0..."];
    ];

    (* Use the detected variable symbol for start/end points.
       DiffExp now handles singular endpoints by returning the series expansion
       instead of trying to evaluate at the singularity. *)
    startPoint = Association[detectedVar -> SetPrecision[fixedVal, precision]];

    (* Transport towards xx = 0 (singular endpoint handled by DiffExp) *)
    resultToLower = DiffExp`Transport`TransportTo[
      {startPoint, bcs},
      Association[detectedVar -> 0],
      1,  (* endpoint *)
      True  (* SaveExpansions *)
    ];

    If[verbosity >= 2,
      Print["  Lower transport result head: ", Head[resultToLower]];
      If[AssociationQ[resultToLower], Print["  Lower keys: ", Keys[resultToLower]]];
    ];

    If[verbosity >= 1,
      Print["  Transporting from xx=", fixedVal, " towards 1..."];
    ];

    (* Reload config for transport in other direction *)
    DiffExp`LoadConfiguration[diffExpConfig];

    (* CRITICAL: Re-add delta prescriptions after LoadConfiguration reset *)
    DiffExp`UpdateConfiguration[{
      DiffExp`State`DeltaPrescriptions -> {
        {detectedVar, 1},        (* x + I*delta at x=0 *)
        {1 - detectedVar, 1}     (* (1-x) + I*delta at x=1 *)
      },
      "AbortOnAnalyticContinuationFail" -> False
    }];

    (* Transport towards xx = 1 (singular endpoint handled by DiffExp) *)
    resultToUpper = DiffExp`Transport`TransportTo[
      {startPoint, bcs},
      Association[detectedVar -> 1],
      1,
      True  (* SaveExpansions *)
    ];

    If[verbosity >= 2,
      Print["  Upper transport result head: ", Head[resultToUpper]];
      If[AssociationQ[resultToUpper], Print["  Upper keys: ", Keys[resultToUpper]]];
    ];
  ];  (* End Module with detectedVar *)

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

  The IBP coefficients c_j(x,d) depend on d = 4-2*eps. We expand in eps:
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
    epsPrefactors_List:{}, feynmanParam_:Automatic] :=
Module[{gammaPrefactor, dimVar, epsSymbol, ibpCoeffOrders,
        segData, actualVar, combinedSegments, combinedData,
        prefactorSpec, result, prefacs, numMasters},

  (* Gamma prefactor: Gamma(v1+v2)/(Gamma(v1)*Gamma(v2)) *)
  gammaPrefactor = Gamma[v1 + v2] / (Gamma[v1] * Gamma[v2]);

  numMasters = Length[ibpCoeffs];

  (* Default: no eps-prefactors (all zero) *)
  prefacs = If[Length[epsPrefactors] == numMasters, epsPrefactors,
    Table[0, {numMasters}]];

  (* Expand each IBP coefficient in eps: d -> 4-2*eps *)
  dimVar = FeynmanTrick`Private`$FTConfig["DimensionVariable"];
  epsSymbol = FeynmanTrick`Private`$FTConfig["EpsilonSymbol"];

  (* ibpCoeffOrders[[j]] is a list of eps-order coefficients for master j
     Each element is a rational function of x (no eps) *)
  ibpCoeffOrders = Table[
    Module[{expanded},
      If[ibpCoeffs[[j]] === 0,
        Table[0, {epsOrder + Max[prefacs] + 1}],
        expanded = ibpCoeffs[[j]] /. dimVar -> (4 - 2*epsSymbol);
        Table[
          SeriesCoefficient[expanded, {epsSymbol, 0, k}] // Normal // Together,
          {k, 0, epsOrder + Max[prefacs]}
        ]
      ]
    ],
    {j, numMasters}
  ];

  (* Determine the actual variable from the segment transform.
     CRITICAL: The transport segments use a variable created by DiffExp's ToExpression
     (e.g., FeynmanTrick`DiffExpIntegration`Private`xx3), while the IBP coefficients
     from FIRE6 use a different symbol (e.g., Global`xx3). We must detect BOTH
     and normalize the IBP coefficients to use the segment variable. *)
  Module[{segVar, ibpVar, varName},
    (* Get the variable from the segment transform *)
    Module[{seg1, transform},
      seg1 = transportResult["SegmentData"][[1]];
      transform = seg1[[2]];
      segVar = If[AssociationQ[transform],
        First[Keys[transform]], Global`xx];
    ];

    (* The IBP variable may be in a different context.
       Find it by looking for symbols with the same name in the IBP coefficients. *)
    varName = SymbolName[segVar];
    ibpVar = segVar;  (* Default: assume same *)

    (* Check if the IBP coefficients contain a different symbol with the same name *)
    Module[{allSymbols, ibpSymbols},
      allSymbols = Cases[ibpCoeffOrders, _Symbol, Infinity] // DeleteDuplicates;
      ibpSymbols = Select[allSymbols,
        SymbolName[#] === varName && # =!= segVar &];
      If[Length[ibpSymbols] > 0,
        ibpVar = First[ibpSymbols];
        If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 2,
          Print["    Variable context mismatch: IBP uses ", ibpVar, " (",
                Context[ibpVar], "), segments use ", segVar, " (", Context[segVar], ")"];
          Print["    Remapping IBP coefficients..."];
        ];
        (* Remap IBP coefficients to use the segment variable *)
        ibpCoeffOrders = ibpCoeffOrders /. ibpVar -> segVar;
      ];
    ];

    actualVar = segVar;
  ];

  segData = transportResult["SegmentData"];

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
            combinedEpsOrders, numEpsOrders},
      seg = segData[[segIdx]];
      seriesRaw = seg[[5]];

      (* Uncompress and resolve theta functions *)
      If[StringQ[seriesRaw],
        If[FileExistsQ[seriesRaw],
           uncompressed = Uncompress[Import[seriesRaw]];,
           uncompressed = Uncompress[seriesRaw];
        ];,
        uncompressed = seriesRaw;
      ];
      uncompressed = uncompressed /. $thetaRules;

      (* uncompressed[[j]] = list of eps orders for master j
         uncompressed[[j]][[n+1]] = SeriesData for J_j at eps order n *)
      numEpsOrders = Length[uncompressed[[1]]];

      (* Get the coordinate transformation: x_main = f(x_local) *)
      xLocal = DiffExp`Symbols`x;
      xMainExpr = If[AssociationQ[seg[[2]]],
        First[Values[seg[[2]]]],
        seg[[2, 2]]  (* Rule format: variable -> expression *)
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
                  If[ibpCoeffOrders[[j, k + 1]] =!= 0,
                    Module[{jIdx, cLocal, cSeries, prod},
                      jIdx = n + kj - k + 1;  (* 1-based index into eps orders *)

                      If[jIdx >= 1 && jIdx <= numEpsOrders,
                        (* Convert IBP coefficient to local coordinates *)
                        cLocal = ibpCoeffOrders[[j, k + 1]] /. actualVar -> xMainExpr;

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
                          , If[NumericQ[transportVal] && Chop[transportVal] =!= 0,
                              (* Transport is a non-zero constant (e.g., from a master
                                 with constant boundary value). Expand c_j(x)*const
                                 as a series in local coords. *)
                              Module[{expOrd = 30},
                                cSeries = Quiet[
                                  Series[cLocal * transportVal, {xLocal, 0, expOrd}],
                                  {Power::infy, Infinity::indet, General::indet}
                                ];
                                If[MatchQ[cSeries, _SeriesData],
                                  total = total + cSeries;
                                ];
                              ];
                            ];
                            (* Zero or non-numeric: skip *)
                          ];
                        ];
                      ];
                    ];
                  ];
                , {k, 0, Min[n + kj, Length[ibpCoeffOrders[[j]]] - 1]}];
              ];
            ];
          , {j, numMasters}];
          total
        ],
        {n, 0, epsOrder}  (* eps orders 0 through epsOrder *)
      ];

      (* Package as single-master segment:
         {line, transform, mainBounds, localBounds, {{epsOrder0, epsOrder1, ...}}} *)
      {seg[[1]], seg[[2]], seg[[3]], seg[[4]], {combinedEpsOrders}}
    ],
    {segIdx, Length[segData]}
  ];

  combinedData = <|
    "SegmentData" -> combinedSegments,
    "NumIntegrals" -> 1,
    "EpsilonOrder" -> epsOrder
  |>;

  (* Check if all combined series are zero (no contributions from any master).
     This happens when all non-zero IBP coefficients multiply zero transport data.
     In this case, the boundary value is zero - skip DefiniteIntegralWithPrefactor. *)
  Module[{allZero},
    allZero = And @@ Table[
      And @@ Table[
        combinedSegments[[s, 5, 1, n]] === 0,
        {n, epsOrder + 1}
      ],
      {s, Length[combinedSegments]}
    ];
    If[allZero,
      If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 2,
        Print["    [Combine] All combined series are zero -> returning zero boundary"];
      ];
      Return[Table[0, {epsOrder + 1}], Module];
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

  result = DiffExp`RegularizedIntegration`DefiniteIntegralWithPrefactor[
    combinedData, {0, 1}, prefactorSpec
  ];

  If[ListQ[result] && Length[result] >= 1 && ListQ[result[[1]]],
    gammaPrefactor * result[[1]]
  ,
    If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
      Print["Warning: Combined integration returned unexpected format: ", Head[result]];
    ];
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
        taylorPart, limitVal, dimVar, epsSymbol, ibpCoeffsExpanded},

  segData = transportResult["SegmentData"];
  numMasters = transportResult["NumIntegrals"];

  (* Expand IBP coefficients in eps: d -> 4 - 2*eps *)
  dimVar = FeynmanTrick`Private`$FTConfig["DimensionVariable"];
  epsSymbol = FeynmanTrick`Private`$FTConfig["EpsilonSymbol"];
  ibpCoeffsExpanded = Table[
    Module[{coeff, expanded},
      coeff = ibpCoeffs[[i]] /. dimVar -> (4 - 2*epsSymbol);
      (* Series expand in eps to required order *)
      expanded = coeff + O[epsSymbol]^(epsOrder + 1);
      (* Extract coefficient at each eps order *)
      Table[SeriesCoefficient[expanded, k], {k, 0, epsOrder}]
    ],
    {i, Length[ibpCoeffs]}
  ];

  (* Select the boundary segment *)
  If[boundary === 0,
    (* First segment: nearest to x=0 *)
    targetSeg = First[SortBy[segData, Min[#[[3]]] &]];
  ,
    (* Last segment: nearest to x=1 *)
    targetSeg = First[SortBy[segData, -Max[#[[3]]] &]];
  ];

  (* Uncompress series data and resolve theta functions *)
  If[StringQ[targetSeg[[5]]],
    If[FileExistsQ[targetSeg[[5]]],
       uncompressed = Uncompress[Import[targetSeg[[5]]]];,
       uncompressed = Uncompress[targetSeg[[5]]];
    ];,
    uncompressed = targetSeg[[5]];
  ];
  uncompressed = uncompressed /. $thetaRules;

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

      (* Add weighted contribution using eps-expanded coefficients *)
      (* ibpCoeffsExpanded[[mIdx]] is a list of eps-order coefficients *)
      (* Convolve with limitVal to get the product in eps *)
      Module[{coeffsAtM, contribution},
        coeffsAtM = ibpCoeffsExpanded[[mIdx]];
        contribution = Table[
          Sum[
            If[k >= 0 && k <= epsOrder && (n - k) >= 0 && (n - k) <= epsOrder,
              coeffsAtM[[k + 1]] * limitVal[[n - k + 1]],
              0
            ],
            {k, 0, n}
          ],
          {n, 0, epsOrder}
        ];
        limitValues += contribution;
      ];
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

(* Helper: expand IBP coefficient in eps, returning list of eps-order coefficients *)
ExpandIBPCoeffInEps[coeff_, epsOrder_Integer] :=
Module[{dimVar, epsSymbol, expanded},
  dimVar = FeynmanTrick`Private`$FTConfig["DimensionVariable"];
  epsSymbol = FeynmanTrick`Private`$FTConfig["EpsilonSymbol"];
  expanded = coeff /. dimVar -> (4 - 2*epsSymbol);
  Table[SeriesCoefficient[expanded + O[epsSymbol]^(epsOrder + 1), k], {k, 0, epsOrder}]
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

ComputeLevelBoundary[ftData_Association, level_Integer,
    transportResult_Association, epsOrder_Integer] :=
Module[{levelData, levelAbove, mastersAtLevel, mastersAbove,
        combinedPositions, posI, posJ, topologyAbove, bcValues, feynmanParamAbove},

  levelData = ftData["Levels"][level];
  levelAbove = ftData["Levels"][level + 1];
  mastersAtLevel = levelData["Masters"];
  mastersAbove = levelAbove["Masters"];
  combinedPositions = levelAbove["CombinedPositions"];
  {posI, posJ} = combinedPositions;
  topologyAbove = levelAbove["Topology"];

  (* Get the Feynman parameter for the level above - needed for IBP coefficient expansion *)
  feynmanParamAbove = levelAbove["FeynmanParameter"];

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
          (* Full integration: combine Sum_j c_j(x)*f_j(x) first, then integrate.
             Individual c_j(x) may have poles at x=0,1 that cancel in the sum.
             Also accounts for eps-prefactors (J vs I basis). *)
          Module[{epsPrefacs, diffMatAbove, epsSymbol},
            (* Compute eps-prefactors from the diff matrix at level+1.
               These tell us J_j = eps^{k_j} * I_j, where k_j accounts for
               the transformation applied during matrix export. *)
            epsSymbol = FeynmanTrick`Private`$FTConfig["EpsilonSymbol"];
            diffMatAbove = levelAbove["DiffMatrix"];
            If[MatrixQ[diffMatAbove] &&
               FeynmanTrick`EpsPrefactors`CheckEpsPoles[diffMatAbove, epsSymbol],
              epsPrefacs = FeynmanTrick`EpsPrefactors`FindEpsPrefactors[diffMatAbove, epsSymbol];
            ,
              epsPrefacs = Table[0, {Length[mastersAbove]}];
            ];

            If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 2,
              Print["    Eps prefactors for level above: ", epsPrefacs];
            ];

            IntegrateCombinedMasters[
              transportResult, ibpCoeffs, vi, vj, epsOrder,
              epsPrefacs, feynmanParamAbove
            ]
          ],

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
                (* Expand IBP coeff in eps and convolve with boundary values *)
                Module[{coeffExpanded, bcAbove},
                  coeffExpanded = ExpandIBPCoeffInEps[ibpCoeffs[[j]], epsOrder];
                  bcAbove = transportResult["BoundaryValuesAbove"][[j]];
                  directVal += MultiplyEpsCoeffs[coeffExpanded, bcAbove, epsOrder];
                ];
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
