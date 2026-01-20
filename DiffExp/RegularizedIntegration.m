(* ::Package:: *)

(* DiffExp Regularized Integration Subpackage *)
(*
   Implements integration of DiffExp series expansions using the regularization
   prescription for non-integrable singularities at boundaries.

   IMPORTANT: Use UseMobius -> False when computing integrals!
   Mobius transforms complicate the coordinate transformations. With UseMobius -> False,
   segments use linear transforms which are handled correctly.

   COORDINATE SYSTEM:
   - DefiniteIntegral and IntegratePiecewiseSaved take bounds in the main line parameter x_main
   - If main line is: kinematic_var -> a + b*x_main, then the bounds are in x_main, not kinematic_var
   - To convert ∫ f(z) dz to ∫ f(x) dx, multiply the result by dz/dx = b

   Example: If main line is z -> 0.01 + 0.49*x, and you want ∫_{z=0.01}^{z=0.3} f(z) dz:
     x_start = (0.01 - 0.01) / 0.49 = 0
     x_end = (0.3 - 0.01) / 0.49 ≈ 0.592
     result_x = DefiniteIntegral[savedData, {0, 0.592}]
     result_z = 0.49 * result_x  (* multiply by dz/dx *)

   The regularization formula for integrals of the form:
     Integrate[x^(a + b*eps) * g(x), {x, 0, c}]
   when a <= -1 and b != 0 is:
     Integrate[x^(a+b*eps+1)/(1+a+b*eps) * ((2+a+b*eps)/c * g(x) - (1-x/c) * g'(x)), {x, 0, c}]

   This increases the power by 1. Repeated application resolves non-integrable singularities.
   When a = -1 and b != 0, the prefactor 1/(1+a+b*eps) = 1/(b*eps) generates an eps^{-1} pole.
   This is tracked via an epsMinPower offset in the returned results.

   For limits at singularities:
   - Terms with b != 0: set to zero (they vanish under analytic regularization)
   - Taylor terms (b = 0): evaluate normally
*)

BeginPackage["DiffExp`RegularizedIntegration`", {
  "DiffExp`Symbols`",
  "DiffExp`State`",
  "DiffExp`Utilities`",
  "DiffExp`SeriesOps`",
  "DiffExp`Integration`",
  "DiffExp`Pade`",
  "DiffExp`SingularityDecomposition`",
  "DiffExp`LineSegmentation`",
  "DiffExp`AnalyticContinuation`"
}];

(* Main exported functions *)
ApplyRegularizationStep::usage = "ApplyRegularizationStep[a, b, gList, c] applies one step of the regularization formula to the integrand x^(a+b*eps)*g(x,eps), increasing the power by 1. Returns {a+1, b, newGList}.";

RegularizeIntegrand::usage = "RegularizeIntegrand[a, b, gList, c] repeatedly applies the regularization formula until a >= 0. Returns {newA, b, newGList}.";

IntegrateSingularTerm::usage = "IntegrateSingularTerm[a, b, gList, {xmin, xmax}] integrates x^(a+b*eps)*g(x,eps) from xmin to xmax using regularization if needed. Returns list of results per epsilon order.";

IntegrateDecomposition::usage = "IntegrateDecomposition[decomposition, {xmin, xmax}] integrates a decomposed series (from DecomposeSingularity) over the given interval.";

EvaluateLimitAtSingularity::usage = "EvaluateLimitAtSingularity[decomposition, direction] evaluates the limit of a decomposed series at x=0. direction is +1 for x->0+ or -1 for x->0-. Terms with b != 0 are set to zero per regularization prescription.";

IntegrateSegmentData::usage = "IntegrateSegmentData[segmentData, {a, b}] integrates a single segment's data over the portion [a,b] (in main line coordinates).";

IntegratePiecewiseSaved::usage = "IntegratePiecewiseSaved[savedData, {a, b}] integrates the piecewise function (from TransportTo with save=True) from a to b in main line x coordinates. Returns a list indexed by {integralIndex, epsilonOrder}. To get the integral in kinematic coordinates, multiply by dz/dx.";

DefiniteIntegral::usage = "DefiniteIntegral[savedData, {a, b}] computes definite integrals from a to b in main line x coordinates. Returns results indexed by {integralIndex, epsilonOrder}. Note: Use UseMobius -> False when transporting. To convert to kinematic variable z, multiply by dz/dx.";

IndefiniteIntegral::usage = "IndefiniteIntegral[savedData] returns a piecewise function representing the indefinite integral. The integration constant is fixed such that the integral is 0 at the start of the first segment.";

Begin["`Private`"];

(* ============================================================ *)
(* Regularization Formula Implementation *)
(* ============================================================ *)

(*
   We track epsilon expansion with an offset to handle negative powers.
   Format: {epsMinPower, gList} where gList[[k]] is coefficient of eps^(epsMinPower + k - 1)

   When applying regularization with a = -1, we get 1/(b*eps) which shifts epsMinPower by -1.
*)

(* Apply one step of the regularization formula:
   Integrate[x^(a+b*eps)*g(x), {x,0,c}] =
     Integrate[x^(a+1+b*eps)/(1+a+b*eps) * ((2+a+b*eps)/c*g(x) - (1-x/c)*g'(x)), {x,0,c}]

   Input: a (rational), b (rational), epsMinPower (integer), gList (list of series per eps order), c (upper bound)
   Output: {a+1, b, newEpsMinPower, newGList}

   The gList represents coefficients of eps^(epsMinPower), eps^(epsMinPower+1), ...
*)
ApplyRegularizationStep[a_, b_, epsMinPower_Integer, gList_List, c_] := Module[
  {epsOrder, newGList, newEpsMin, transformedG, n},

  epsOrder = Length[gList] - 1;

  (* First compute the transformed g: (2+a+b*eps)/c * g - (1-x/c) * g'
     This doesn't change the eps structure, just transforms each coefficient *)
  transformedG = Table[
    Module[{gAtOrd, gAtOrdMinus1, gpAtOrd},
      gAtOrd = gList[[ord + 1]];
      gAtOrdMinus1 = If[ord > 0, gList[[ord]], 0];
      gpAtOrd = DiffExp`SeriesOps`SD[gAtOrd, DiffExp`Symbols`x];

      (* (2+a)/c * g_ord + b/c * g_{ord-1} - (1-x/c) * g'_ord *)
      DiffExp`SeriesOps`SExpand[
        (2 + a) / c * gAtOrd +
        (If[ord > 0, b / c * gAtOrdMinus1, 0]) -
        (1 - DiffExp`Symbols`x / c) * gpAtOrd
      ]
    ],
    {ord, 0, epsOrder}
  ];

  (* Now handle the prefactor 1/(1+a+b*eps) *)
  If[1 + a == 0,
    (* Special case: a = -1, prefactor is 1/(b*eps) = (1/b) * eps^{-1} *)
    (* This shifts eps powers down by 1 and multiplies by 1/b *)
    newEpsMin = epsMinPower - 1;
    newGList = (1/b) * # & /@ transformedG;
    ,
    (* General case: 1/(1+a+b*eps) expands as geometric series *)
    (* 1/(1+a+b*eps) = 1/(1+a) * Sum[(-b*eps/(1+a))^k, {k,0,inf}] *)
    newEpsMin = epsMinPower;
    newGList = Table[
      Sum[
        Module[{prefactorCoeff},
          prefactorCoeff = (-b / (1 + a))^k / (1 + a);
          If[n - k >= 0 && n - k <= epsOrder,
            prefactorCoeff * transformedG[[n - k + 1]],
            0
          ]
        ],
        {k, 0, n}
      ] // DiffExp`SeriesOps`SExpand,
      {n, 0, epsOrder}
    ];
  ];

  {a + 1, b, newEpsMin, newGList}
];

(* Backward compatible version without explicit epsMinPower (assumes 0) *)
ApplyRegularizationStep[a_, b_, gList_List, c_] :=
  Module[{result},
    result = ApplyRegularizationStep[a, b, 0, gList, c];
    (* Return {a+1, b, newGList} and separately track epsMin if needed *)
    {result[[1]], result[[2]], result[[3]], result[[4]]}
  ];

(* Old interface for compatibility - returns {a+1, b, newGList} *)
ApplyRegularizationStepOld[a_, b_, gList_List, c_] := Module[
  {epsOrder, newGList, transformedG, n},

  epsOrder = Length[gList] - 1;

  transformedG = Table[
    Module[{gAtOrd, gAtOrdMinus1, gpAtOrd},
      gAtOrd = gList[[ord + 1]];
      gAtOrdMinus1 = If[ord > 0, gList[[ord]], 0];
      gpAtOrd = DiffExp`SeriesOps`SD[gAtOrd, DiffExp`Symbols`x];
      DiffExp`SeriesOps`SExpand[
        (2 + a) / c * gAtOrd +
        (If[ord > 0, b / c * gAtOrdMinus1, 0]) -
        (1 - DiffExp`Symbols`x / c) * gpAtOrd
      ]
    ],
    {ord, 0, epsOrder}
  ];

  If[1 + a == 0,
    DiffExp`Utilities`ReportError["ApplyRegularizationStepOld: a = -1 encountered. Use ApplyRegularizationStep with epsMinPower."];
  ];

  newGList = Table[
    Sum[
      Module[{prefactorCoeff},
        prefactorCoeff = (-b / (1 + a))^k / (1 + a);
        If[n - k >= 0 && n - k <= epsOrder,
          prefactorCoeff * transformedG[[n - k + 1]],
          0
        ]
      ],
      {k, 0, n}
    ] // DiffExp`SeriesOps`SExpand,
    {n, 0, epsOrder}
  ];

  {a + 1, b, newGList}
];

(* Repeatedly apply regularization until a >= 0, tracking eps offset *)
(* Returns {newA, b, epsMinPower, newGList} *)
RegularizeIntegrand[a_, b_, gList_List, c_] := Module[
  {currentA = a, currentB = b, currentEpsMin = 0, currentG = gList, maxIter = 50, iter = 0, result},

  While[currentA < 0 && iter < maxIter,
    iter++;
    result = ApplyRegularizationStep[currentA, currentB, currentEpsMin, currentG, c];
    currentA = result[[1]];
    currentB = result[[2]];
    currentEpsMin = result[[3]];
    currentG = result[[4]];
  ];

  If[iter >= maxIter,
    DiffExp`Utilities`PrintWarning["RegularizeIntegrand: reached maximum iterations"];
  ];

  {currentA, currentB, currentEpsMin, currentG}
];

(* ============================================================ *)
(* Integration of Individual Terms *)
(* ============================================================ *)

(* Integrate x^(a+b*eps) * g(x,eps) from xmin to xmax.
   Handles:
   - a >= 0: direct integration
   - a < 0, b != 0: regularization first
   - a < 0, b = 0: only if a > -1 (otherwise non-integrable)

   For definite integrals with singularity at boundary:
   - At x=0: singular terms with b != 0 contribute 0
   - At x=c: evaluate normally
*)
IntegrateSingularTerm[a_, b_, gList_List, {xmin_, xmax_}] := Module[
  {epsOrder, regA, regB, regEpsMin, regG, integratedG, upperEval, lowerEval, result,
   atLowerSingularity, atUpperSingularity, regResult},

  epsOrder = Length[gList] - 1;

  (* Check if boundaries are at singularities *)
  atLowerSingularity = (Abs[xmin] < DiffExp`State`FEC[RationalizationTolerance]);
  atUpperSingularity = (Abs[xmax] < DiffExp`State`FEC[RationalizationTolerance]);

  (* If a < 0 and b != 0, we need regularization *)
  If[a < 0 && b != 0,
    (* Apply regularization to make a >= 0 *)
    regResult = RegularizeIntegrand[a, b, gList, xmax];
    regA = regResult[[1]];
    regB = regResult[[2]];
    regEpsMin = regResult[[3]];
    regG = regResult[[4]];
    ,
    (* No regularization needed *)
    regA = a; regB = b; regEpsMin = 0; regG = gList;
  ];

  (* Now integrate x^(regA + regB*eps) * g(x,eps) *)
  (* The integral is: x^(regA+1+regB*eps)/(regA+1+regB*eps) * G(x,eps)
     where G is the integrated g, plus correction terms from the power *)

  (* Actually, we need to be more careful. The full integrand is:
     x^(regA + regB*eps) * g(x,eps)

     Let's integrate term by term in x. For a series g = sum_k c_k x^k:
     Integrate[x^(regA+regB*eps) * sum_k c_k x^k] = sum_k c_k * x^(regA+regB*eps+k+1) / (regA+regB*eps+k+1)

     The denominator 1/(regA+regB*eps+k+1) needs to be expanded in eps.
  *)

  result = Table[
    Module[{gAtOrd, norms, integrated, upperVal, lowerVal},
      gAtOrd = regG[[ord + 1]];

      (* Integrate the series: for each x^k term in g, we get x^(regA+regB*eps+k+1)/(regA+regB*eps+k+1) *)
      (* This is complex because regB*eps affects the power and the denominator *)

      (* For now, let's use a simpler approach:
         If regA >= 0 and regB = 0 (Taylor term), just integrate directly.
         If regB != 0, we need the full eps-dependent integration. *)

      If[regB == 0,
        (* Simple Taylor integration *)
        integrated = DiffExp`Integration`DiffExpIntegrate[gAtOrd * DiffExp`Symbols`x^regA, DiffExp`Symbols`x];
        upperVal = DiffExp`Pade`SEval[integrated, xmax];
        lowerVal = If[atLowerSingularity && regA < 0, 0, DiffExp`Pade`SEval[integrated, xmin]];
        upperVal - lowerVal
        ,
        (* eps-dependent power: x^(regA+regB*eps) *)
        (* Need to handle this more carefully *)
        IntegratePowerTimesSeriesAtOrder[regA, regB, regG, ord, {xmin, xmax}, atLowerSingularity, atUpperSingularity]
      ]
    ],
    {ord, 0, epsOrder}
  ];

  result
];

(* Helper: integrate x^(a+b*eps) * g(x,eps) at a specific epsilon order *)
(* This handles the eps-expansion of the power and denominator *)
IntegratePowerTimesSeriesAtOrder[a_, b_, gList_List, targetOrd_, {xmin_, xmax_}, atLower_, atUpper_] := Module[
  {epsOrder, result, k, m, contrib, gCoeffs, intCoeff, upperVal, lowerVal},

  epsOrder = Length[gList] - 1;

  (* x^(a+b*eps) = x^a * x^(b*eps) = x^a * exp(b*eps*Logx) = x^a * sum_k (b*Logx)^k/k! * eps^k *)
  (* So at eps^n, we get: sum_{k=0}^{n} (b*Logx)^k/k! * g_{n-k}(x) * x^a *)

  (* The integral of x^a * Logx^k * (series in x) is computed by:
     - Expand the series
     - For each term x^m * Logx^k, integrate to get:
       x^(a+m+1)/(a+m+1) * (Logx^k - k*Logx^(k-1)/(a+m+1) + k*(k-1)*Logx^(k-2)/(a+m+1)^2 - ...)
  *)

  result = Sum[
    Module[{logPow = k, logCoeff, gAtOrder, integ, uVal, lVal},
      logCoeff = b^logPow / logPow!;
      gAtOrder = If[targetOrd - k >= 0 && targetOrd - k <= epsOrder,
                    gList[[targetOrd - k + 1]],
                    0];

      If[gAtOrder === 0,
        0,
        (* Integrate x^a * Logx^logPow * gAtOrder *)
        integ = IntegrateWithLogPower[a, logPow, gAtOrder];

        (* Evaluate at bounds *)
        uVal = EvaluateIntegralAtPoint[integ, xmax, a, atUpper];
        lVal = EvaluateIntegralAtPoint[integ, xmin, a, atLower];

        logCoeff * (uVal - lVal)
      ]
    ],
    {k, 0, targetOrd}
  ];

  result // DiffExp`Utilities`PChop // Expand
];

(* Integrate x^a * Logx^n * (series in x) *)
(* Returns a function/expression that can be evaluated at bounds *)
IntegrateWithLogPower[a_, n_, ser_SeriesData] := Module[
  {nmin, nmax, den, coeffs, result, m, intPow, logTerms},

  nmin = ser[[4]];
  nmax = ser[[5]];
  den = ser[[6]];
  coeffs = ser[[3]];

  (* For each term c_m * x^(m/den) in the series:
     Integrate x^a * Logx^n * c_m * x^(m/den) = c_m * Integrate x^(a+m/den) * Logx^n

     Using the formula:
     Integrate x^p * Logx^n = x^(p+1)/(p+1) * sum_{j=0}^{n} (-1)^j * n!/(n-j)! * Logx^(n-j) / (p+1)^j
  *)

  result = Sum[
    Module[{coeff, powIdx, p, intResult},
      powIdx = nmin + idx - 1;
      coeff = If[idx <= Length[coeffs], coeffs[[idx]], 0];
      p = a + powIdx / den;

      If[coeff === 0 || (NumericQ[coeff] && Abs[coeff] < 10^-DiffExp`State`FEC[ChopPrecision]),
        0,
        (* Integrate x^p * Logx^n *)
        If[p == -1,
          (* Special case: gives Logx^(n+1)/(n+1) *)
          coeff * DiffExp`Symbols`Logx^(n + 1) / (n + 1)
          ,
          (* General case *)
          coeff * DiffExp`Symbols`x^(p + 1) / (p + 1) *
            Sum[(-1)^j * Factorial[n] / Factorial[n - j] * DiffExp`Symbols`Logx^(n - j) / (p + 1)^j, {j, 0, n}]
        ]
      ]
    ],
    {idx, 1, Length[coeffs]}
  ];

  result // Expand
];

IntegrateWithLogPower[a_, n_, 0] := 0;
IntegrateWithLogPower[a_, n_, c_?NumericQ] := IntegrateWithLogPower[a, n,
  SeriesData[DiffExp`Symbols`x, 0, {c}, 0, 1, 1]];

(* Evaluate the integral expression at a point, handling singularity limits *)
EvaluateIntegralAtPoint[expr_, pt_, a_, atSingularity_] := Module[{val},
  If[atSingularity && pt == 0,
    (* At x=0 singularity: terms with Logx or negative powers vanish *)
    (* Only finite constant terms survive *)
    (expr /. DiffExp`Symbols`Logx -> 0 /. DiffExp`Symbols`x -> 0) //
      Quiet[Limit[#, DiffExp`Symbols`x -> 0]] & //
      If[NumericQ[#], #, 0] &
    ,
    (* Normal evaluation *)
    expr /. DiffExp`Symbols`Logx -> Log[DiffExp`Symbols`x] /. DiffExp`Symbols`x -> pt
  ]
];

(* ============================================================ *)
(* Integration of Decomposed Series *)
(* ============================================================ *)

(* Integrate a full decomposition from DecomposeSingularity *)
IntegrateDecomposition[decomposition_List, {xmin_, xmax_}] := Module[
  {epsOrder, result, term, a, b, g},

  If[Length[decomposition] == 0, Return[{}]];

  epsOrder = Length[decomposition[[1]]["g"]] - 1;

  (* Sum contributions from all terms *)
  result = Table[0, {epsOrder + 1}];

  Do[
    a = decomposition[[i]]["a"];
    b = decomposition[[i]]["b"];
    g = decomposition[[i]]["g"];

    result = result + IntegrateSingularTerm[a, b, g, {xmin, xmax}];
    ,
    {i, Length[decomposition]}
  ];

  result
];

(* ============================================================ *)
(* Limit Evaluation at Singularities *)
(* ============================================================ *)

(* Evaluate limit of decomposed series at x=0
   Per regularization prescription:
   - Terms with b != 0: set to 0
   - Taylor terms (a >= 0, b = 0): evaluate g at x=0
*)
EvaluateLimitAtSingularity[decomposition_List, direction_: 1] := Module[
  {epsOrder, result, term, a, b, g, gVal},

  If[Length[decomposition] == 0, Return[{}]];

  epsOrder = Length[decomposition[[1]]["g"]] - 1;
  result = Table[0, {epsOrder + 1}];

  Do[
    a = decomposition[[i]]["a"];
    b = decomposition[[i]]["b"];
    g = decomposition[[i]]["g"];

    If[b == 0 && a >= 0,
      (* Taylor term: evaluate g at x=0 *)
      gVal = Table[
        If[MatchQ[g[[ord + 1]], _SeriesData],
          (* Get constant term of series *)
          SeriesCoefficient[g[[ord + 1]], {DiffExp`Symbols`x, 0, 0}],
          g[[ord + 1]]
        ],
        {ord, 0, epsOrder}
      ];
      result = result + gVal;
    ];
    (* Terms with b != 0 contribute 0 *)
    ,
    {i, Length[decomposition]}
  ];

  result
];

(* ============================================================ *)
(* Segment Integration *)
(* ============================================================ *)

(* Integrate a single segment's data over [a, b] in main line coordinates.
   segmentData format: {CurrLine, lineRelation, {xmin_main, xmax_main}, {xmin_local, xmax_local}, seriesData}
*)
IntegrateSegmentData[segmentData_List, {a_, b_}, intIndex_Integer, epsOrder_Integer] := Module[
  {currLine, lineRelation, mainBounds, localBounds, seriesData,
   jacobian, localA, localB, decomposition, integralResult,
   seriesAtIndex, uncompressed, atLowerSing, atUpperSing,
   mainMin, mainMax, localMin, localMax, slope},

  currLine = segmentData[[1]];
  lineRelation = segmentData[[2]];
  mainBounds = segmentData[[3]];
  localBounds = segmentData[[4]];
  seriesData = segmentData[[5]];

  (* Uncompress if needed *)
  If[StringQ[seriesData],
    uncompressed = Uncompress[Import[seriesData]];,
    If[Head[seriesData] === String && StringLength[seriesData] > 100,
      uncompressed = Uncompress[seriesData];,
      uncompressed = seriesData;
    ]
  ];

  (* Get the series for this integral (all eps orders) *)
  seriesAtIndex = uncompressed[[intIndex]];

  (* Compute Jacobian: dx_main/dx_local *)
  (* lineRelation is x -> expression, so jacobian = D[expression, x] *)
  jacobian = D[lineRelation[[2]], DiffExp`Symbols`x];

  (* Convert main line bounds [a, b] to local coordinates *)
  (* Use linear interpolation based on the stored bounds - this works for linear transforms *)
  (* For Mobius transforms, use the existing RelateLines infrastructure *)
  mainMin = mainBounds[[1]];
  mainMax = mainBounds[[2]];
  localMin = localBounds[[1]];
  localMax = localBounds[[2]];

  (* For linear segments, we can interpolate directly *)
  (* localCoord = localMin + (mainCoord - mainMin) * (localMax - localMin) / (mainMax - mainMin) *)
  If[mainMax - mainMin != 0,
    slope = (localMax - localMin) / (mainMax - mainMin);
    localA = localMin + (a - mainMin) * slope;
    localB = localMin + (b - mainMin) * slope;
    (* For linear transforms, jacobian is constant = 1/slope (dx_main/dx_local) *)
    (* But we compute it from lineRelation to be safe *)
    jacobian = jacobian /. DiffExp`Symbols`x -> (localA + localB)/2;  (* Evaluate at midpoint *)
    ,
    (* Degenerate case: single point *)
    localA = localMin;
    localB = localMin;
    jacobian = 1;
  ];

  (* Check if we're at singularities *)
  atLowerSing = Abs[localA] < DiffExp`State`FEC[RationalizationTolerance];
  atUpperSing = Abs[localB] < DiffExp`State`FEC[RationalizationTolerance];

  (* Decompose the series near singularity *)
  decomposition = DiffExp`SingularityDecomposition`DecomposeSingularity[seriesAtIndex];

  (* Integrate the decomposition *)
  integralResult = IntegrateDecomposition[decomposition, {localA, localB}];

  (* Apply Jacobian - lineRelation maps x_main -> x_local, so jacobian = dx_local/dx_main *)
  (* For the integral: ∫ f(x_main) dx_main = ∫ g(x_local) * |dx_main/dx_local| dx_local = (1/|jacobian|) * ∫ g dx_local *)
  integralResult[[epsOrder + 1]] / Abs[jacobian]
];

(* ============================================================ *)
(* Main Entry Points *)
(* ============================================================ *)

(* Integrate piecewise saved data from TransportTo *)
IntegratePiecewiseSaved[savedData2_, {a_, b_}] := Module[
  {savedData, numIntegrals, numEpsOrders, relevantSegments,
   segmentContributions, result, seg, segBounds, overlapBounds,
   uncompressed},

  (* Handle wrapper format from TransportTo *)
  If[MatchQ[savedData2, {{_Association, __}, {__}}],
    savedData = savedData2[[2]],
    savedData = savedData2
  ];

  (* Validate format *)
  If[!ListQ[savedData] || Length[savedData] == 0,
    DiffExp`Utilities`ReportError["IntegratePiecewiseSaved: Invalid saved data format"];
  ];

  (* Get dimensions from first segment *)
  If[StringQ[savedData[[1, 5]]],
    uncompressed = Uncompress[Import[savedData[[1, 5]]]];,
    If[Head[savedData[[1, 5]]] === String && StringLength[savedData[[1, 5]]] > 100,
      uncompressed = Uncompress[savedData[[1, 5]]];,
      uncompressed = savedData[[1, 5]];
    ]
  ];

  numIntegrals = Length[uncompressed];
  numEpsOrders = Length[uncompressed[[1]]];

  (* Find segments that overlap with [a, b] *)
  relevantSegments = Select[savedData,
    Module[{segMin, segMax},
      segMin = Min[#[[3]]];
      segMax = Max[#[[3]]];
      Not[segMax < a || segMin > b]
    ] &
  ];

  (* Integrate each relevant segment *)
  result = Table[
    Sum[
      seg = relevantSegments[[segIdx]];
      segBounds = seg[[3]];

      (* Find overlap with [a, b] *)
      overlapBounds = {Max[Min[segBounds], a], Min[Max[segBounds], b]};

      If[overlapBounds[[1]] >= overlapBounds[[2]],
        0,
        IntegrateSegmentData[seg, overlapBounds, intIdx, epsOrd]
      ]
      ,
      {segIdx, Length[relevantSegments]}
    ],
    {intIdx, numIntegrals},
    {epsOrd, 0, numEpsOrders - 1}
  ];

  result
];

(* Main entry point for definite integrals *)
DefiniteIntegral[savedData_, {a_, b_}] := Module[{},
  DiffExp`Utilities`PrintInfo["Computing definite integral from ", a // N, " to ", b // N][1];
  IntegratePiecewiseSaved[savedData, {a, b}]
];

(* Indefinite integral - returns piecewise function *)
IndefiniteIntegral[savedData2_] := Module[
  {savedData, numIntegrals, numEpsOrders, uncompressed,
   cumulativeIntegrals, currentValue, seg, segBounds,
   result, segIdx, intIdx, epsOrd},

  (* Handle wrapper format *)
  If[MatchQ[savedData2, {{_Association, __}, {__}}],
    savedData = savedData2[[2]],
    savedData = savedData2
  ];

  (* Get dimensions *)
  If[StringQ[savedData[[1, 5]]],
    uncompressed = Uncompress[Import[savedData[[1, 5]]]];,
    If[Head[savedData[[1, 5]]] === String && StringLength[savedData[[1, 5]]] > 100,
      uncompressed = Uncompress[savedData[[1, 5]]];,
      uncompressed = savedData[[1, 5]];
    ]
  ];

  numIntegrals = Length[uncompressed];
  numEpsOrders = Length[uncompressed[[1]]];

  DiffExp`Utilities`PrintInfo["Computing indefinite integral with ", Length[savedData], " segments"][1];

  (* Build piecewise functions for each integral and eps order *)
  result = Table[
    Module[{pieces, cumulative = 0, prevEnd},
      prevEnd = Min[savedData[[1, 3]]];

      pieces = Table[
        seg = savedData[[segIdx]];
        segBounds = seg[[3]];

        (* Integrate this segment *)
        Module[{localContrib, segStart, segEnd, piece},
          segStart = Min[segBounds];
          segEnd = Max[segBounds];

          (* Add contribution from gap (if any) - shouldn't happen normally *)

          (* Integrate from segment start to variable x *)
          (* This returns a function of x *)
          localContrib = IntegrateSegmentIndefinite[seg, intIdx, epsOrd];

          (* The value at segment start should match cumulative *)
          (* piece = cumulative + (localContrib - localContrib[segStart]) *)
          piece = {
            cumulative + localContrib - (localContrib /. DiffExp`Symbols`x -> segStart),
            segStart <= DiffExp`Symbols`x <= segEnd
          };

          (* Update cumulative for next segment *)
          cumulative = cumulative + (localContrib /. DiffExp`Symbols`x -> segEnd) -
                       (localContrib /. DiffExp`Symbols`x -> segStart);

          piece
        ],
        {segIdx, Length[savedData]}
      ];

      Piecewise[pieces] /. DiffExp`Symbols`x -> # &
    ],
    {intIdx, numIntegrals},
    {epsOrd, 0, numEpsOrders - 1}
  ];

  result
];

(* Helper: integrate a segment indefinitely, returning a function of x *)
IntegrateSegmentIndefinite[segmentData_List, intIndex_Integer, epsOrder_Integer] := Module[
  {currLine, lineRelation, seriesData, jacobian,
   seriesAtIndex, uncompressed, decomposition,
   localX, indefiniteLocal, indefiniteMain},

  currLine = segmentData[[1]];
  lineRelation = segmentData[[2]];
  seriesData = segmentData[[5]];

  (* Uncompress *)
  If[StringQ[seriesData],
    uncompressed = Uncompress[Import[seriesData]];,
    If[Head[seriesData] === String && StringLength[seriesData] > 100,
      uncompressed = Uncompress[seriesData];,
      uncompressed = seriesData;
    ]
  ];

  seriesAtIndex = uncompressed[[intIndex]];

  (* Jacobian *)
  jacobian = D[lineRelation[[2]], DiffExp`Symbols`x] // Simplify;

  (* Get the series at the requested eps order *)
  (* For indefinite integral, we just integrate the series directly *)
  (* This is simpler - we don't need full decomposition *)

  (* Integrate the series *)
  indefiniteLocal = DiffExp`Integration`DiffExpIntegrate[
    seriesAtIndex[[epsOrder + 1]],
    DiffExp`Symbols`x
  ];

  (* Convert to main line coordinates *)
  (* x_local appears in indefiniteLocal, replace with inverse of lineRelation *)
  (* lineRelation: x_main = f(x_local), so x_local = f^{-1}(x_main) *)

  (* For linear lines, x_main = a*x_local + b, so x_local = (x_main - b)/a *)
  Module[{localToMain},
    localToMain = DiffExp`Symbols`x /.
      Solve[lineRelation[[2]] == Global`xmain, DiffExp`Symbols`x][[1]] /.
      Global`xmain -> DiffExp`Symbols`x;

    (* Apply substitution and multiply by Jacobian *)
    jacobian * (Normal[indefiniteLocal] /. DiffExp`Symbols`x -> localToMain) /.
      DiffExp`Symbols`Logx -> Log[DiffExp`Symbols`x]
  ]
];

End[];

EndPackage[];
