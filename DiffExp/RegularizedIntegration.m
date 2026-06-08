(* ::Package:: *)

(* DiffExp Regularized Integration Subpackage *)
(*
   Implements integration of DiffExp series expansions using the regularization
   prescription for non-integrable singularities at boundaries.

   IMPORTANT: Use UseMobius -> False when computing integrals!
   Mobius transforms complicate the coordinate transformations. With UseMobius -> False,
   segments use linear transforms which are handled correctly.

   COORDINATE SYSTEM:
   - All integrals are computed in the main line parameter x: ∫ f(x) dx
   - Bounds {a, b} are in x coordinates
   - The kinematic invariants are functions of x along the main line
   - If you need ∫ f(z) dz where z = z(x), multiply result by dz/dx

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

IntegrateSingularTermLaurent::usage = "IntegrateSingularTermLaurent[a, b, epsMinPower, gList, {xmin, xmax}] integrates x^(a+b*eps)*g(x,eps), where gList starts at eps^epsMinPower. Returns <|\"MinPower\" -> m, \"Coefficients\" -> {...}|>.";

IntegrateDecomposition::usage = "IntegrateDecomposition[decomposition, {xmin, xmax}] integrates a decomposed series (from DecomposeSingularity) over the given interval.";

IntegrateDecompositionLaurent::usage = "IntegrateDecompositionLaurent[decomposition, {xmin, xmax}, epsMinPower] integrates a decomposed series whose epsilon coefficient list starts at eps^epsMinPower. Returns a Laurent coefficient association.";

EvaluateLimitAtSingularity::usage = "EvaluateLimitAtSingularity[decomposition, direction] evaluates the limit of a decomposed series at x=0. direction is +1 for x->0+ or -1 for x->0-. Terms with b != 0 are set to zero per regularization prescription.";

IntegrateSegmentData::usage = "IntegrateSegmentData[segmentData, {a, b}] integrates a single segment's data over the portion [a,b] (in main line coordinates).";

IntegratePiecewiseSaved::usage = "IntegratePiecewiseSaved[savedData, {a, b}] computes ∫ f(x) dx from a to b where x is the main line parameter. Returns a list indexed by {integralIndex, epsilonOrder}.";

DefiniteIntegral::usage = "DefiniteIntegral[savedData, {a, b}] computes ∫ f(x) dx from a to b where x is the main line parameter. Returns results indexed by {integralIndex, epsilonOrder}. Use UseMobius -> False when transporting.";

IndefiniteIntegral::usage = "IndefiniteIntegral[savedData] returns a piecewise function representing the indefinite integral. The integration constant is fixed such that the integral is 0 at the start of the first segment.";

DefiniteIntegralWithPrefactor::usage = "DefiniteIntegralWithPrefactor[savedData, {lower, upper}, prefactorSpec] computes ∫_lower^upper x^alpha * (upper-x)^beta * r(x) * f(x) dx. \
prefactorSpec is <|\"PowerAtLower\" -> alpha, \"PowerAtUpper\" -> beta, \"RationalFactor\" -> r(x), \"Variable\" -> x|>. \
f(x) is the piecewise series from saved transport data. Powers alpha, beta >= 0 for convergence (regularization applied otherwise). \
RationalFactor is a rational function that gets series-expanded at each segment.";

DefiniteIntegralWithPrefactorLaurent::usage = "DefiniteIntegralWithPrefactorLaurent[savedData, {lower, upper}, prefactorSpec, epsMinPower] is the Laurent-aware version of DefiniteIntegralWithPrefactor. It returns one association per integral, each with \"MinPower\" and \"Coefficients\".";

Begin["`Private`"];

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

LaurentScale[c_, laur_Association] := <|
  "MinPower" -> laur["MinPower"],
  "Coefficients" -> (DiffExp`Utilities`PChop /@ (c * laur["Coefficients"]))
|>;

LaurentAdd[a_Association, b_Association] := Module[
  {minPower, maxPower},
  minPower = Min[a["MinPower"], b["MinPower"]];
  maxPower = Max[LaurentMaxPower[a], LaurentMaxPower[b]];
  <|
    "MinPower" -> minPower,
    "Coefficients" -> Table[
      DiffExp`Utilities`PChop[Expand[LaurentCoeff[a, p] + LaurentCoeff[b, p]]],
      {p, minPower, maxPower}
    ]
  |>
];

LaurentToNonNegativeList[laur_Association, epsOrder_Integer] :=
  Table[LaurentCoeff[laur, p], {p, 0, epsOrder}];

KnownSeriesDerivative[ser_SeriesData, var_] := Module[
  {order},
  order = Max[0, Ceiling[ser[[5]] / ser[[6]]]];
  Quiet[
    Series[
      DiffExp`SeriesOps`SD[Normal[ser], var],
      {var, ser[[2]], order}
    ]
  ]
];

KnownSeriesDerivative[expr_, var_] :=
  DiffExp`SeriesOps`SD[expr, var];

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
  {epsOrder, transformedMaxOrder, newGList, newEpsMin, transformedG, n, tol, zeroQ},

  epsOrder = Length[gList] - 1;
  tol = DiffExp`State`FEC[RationalizationTolerance];
  zeroQ[z_] := TrueQ[PossibleZeroQ[z]] || TrueQ[NumericQ[z] && Abs[N[z, 50]] < tol];
  transformedMaxOrder = epsOrder + If[zeroQ[1 + a], 1, 0];

  (* First compute the transformed g: (2+a+b*eps)/c * g - (1-x/c) * g'
     This doesn't change the eps structure, just transforms each coefficient *)
  transformedG = Table[
    Module[{gAtOrd, gAtOrdMinus1, gpAtOrd},
      gAtOrd = If[ord <= epsOrder, gList[[ord + 1]], 0];
      gAtOrdMinus1 = If[ord > 0 && ord - 1 <= epsOrder, gList[[ord]], 0];
      gpAtOrd = KnownSeriesDerivative[gAtOrd, DiffExp`Symbols`x];

      (* (2+a)/c * g_ord + b/c * g_{ord-1} - (1-x/c) * g'_ord *)
      DiffExp`SeriesOps`SExpand[
        (2 + a) / c * gAtOrd +
        (If[ord > 0, b / c * gAtOrdMinus1, 0]) -
        (1 - DiffExp`Symbols`x / c) * gpAtOrd
      ]
    ],
    {ord, 0, transformedMaxOrder}
  ];

  (* Now handle the prefactor 1/(1+a+b*eps) *)
  If[zeroQ[1 + a],
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
          If[n - k >= 0 && n - k <= transformedMaxOrder,
            prefactorCoeff * transformedG[[n - k + 1]],
            0
          ]
        ],
        {k, 0, n}
      ] // DiffExp`SeriesOps`SExpand,
      {n, 0, transformedMaxOrder}
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
      gpAtOrd = KnownSeriesDerivative[gAtOrd, DiffExp`Symbols`x];
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
RegularizeIntegrand[a_, b_, gList_List, c_] :=
  RegularizeIntegrand[a, b, 0, gList, c];

RegularizeIntegrand[a_, b_, epsMinPower_Integer, gList_List, c_] := Module[
  {currentA = a, currentB = b, currentEpsMin = epsMinPower, currentG = gList, maxIter = 50, iter = 0, result},

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
  {epsOrder, laur},
  epsOrder = Length[gList] - 1;
  laur = IntegrateSingularTermLaurent[a, b, 0, gList, {xmin, xmax}];
  LaurentToNonNegativeList[laur, epsOrder]
];

IntegrateSingularTermLaurent[a_, b_, epsMinPower_Integer, gList_List, {xmin_, xmax_}] := Module[
  {epsOrder, regA, regB, regEpsMin, regG, integratedG, upperEval, lowerEval, result,
   atLowerSingularity, atUpperSingularity, regResult, tol, zeroQ, epsMaxPower},

  epsOrder = Length[gList] - 1;
  epsMaxPower = epsMinPower + epsOrder;
  tol = DiffExp`State`FEC[RationalizationTolerance];
  zeroQ[z_] := TrueQ[PossibleZeroQ[z]] || TrueQ[Abs[N[z, 50]] < tol];

  (* The regularization formula is written with the singular point x=0
     as the lower endpoint.  If zero is the upper endpoint, flip the
     interval; if zero lies inside the interval, split into two regulated
     endpoint integrals. *)
  If[zeroQ[xmax] && !zeroQ[xmin],
    Return[LaurentScale[-1, IntegrateSingularTermLaurent[a, b, epsMinPower, gList, {0, xmin}]]]
  ];
  If[!zeroQ[xmin] && !zeroQ[xmax] &&
     TrueQ[N[xmin, 50] < 0 < N[xmax, 50]],
    Return[
      LaurentAdd[
        IntegrateSingularTermLaurent[a, b, epsMinPower, gList, {xmin, 0}],
        IntegrateSingularTermLaurent[a, b, epsMinPower, gList, {0, xmax}]
      ]
    ]
  ];

  (* Check if boundaries are at singularities *)
  atLowerSingularity = If[zeroQ[xmin],
    If[zeroQ[xmax], 1, Sign[N[xmax, 50]]],
    False
  ];
  atUpperSingularity = If[zeroQ[xmax],
    If[zeroQ[xmin], 1, Sign[N[xmin, 50]]],
    False
  ];

  (* If a < 0 and b != 0 at the endpoint x=0, we need regularization.
     Away from x=0 the ordinary termwise antiderivative is finite; applying
     the endpoint regularization formula there creates a spurious boundary term. *)
  If[a < 0 && b != 0 && atLowerSingularity =!= False,
    (* Apply regularization to make a >= 0 *)
    regResult = RegularizeIntegrand[a, b, epsMinPower, gList, xmax];
    regA = regResult[[1]];
    regB = regResult[[2]];
    regEpsMin = regResult[[3]];
    regG = regResult[[4]];
    ,
    (* No regularization needed *)
    regA = a; regB = b; regEpsMin = epsMinPower; regG = gList;
  ];

  epsMaxPower = regEpsMin + Length[regG] - 1;

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
      gAtOrd = LaurentCoeff[<|"MinPower" -> regEpsMin, "Coefficients" -> regG|>, ord];

      (* Integrate the series: for each x^k term in g, we get x^(regA+regB*eps+k+1)/(regA+regB*eps+k+1) *)
      (* This is complex because regB*eps affects the power and the denominator *)

      (* For now, let's use a simpler approach:
         If regA >= 0 and regB = 0 (Taylor term), just integrate directly.
         If regB != 0, we need the full eps-dependent integration. *)

      (* Always use the detailed integration to handle fractional terms correctly *)
      IntegratePowerTimesSeriesAtPower[regA, regB, regEpsMin, regG, ord,
        {xmin, xmax}, atLowerSingularity, atUpperSingularity]
    ],
    {ord, regEpsMin, epsMaxPower}
  ];

  <|"MinPower" -> regEpsMin, "Coefficients" -> result|>
];

(* Helper: integrate x^(a+b*eps) * g(x,eps) at a specific epsilon order *)
(* This handles the eps-expansion of the power and denominator *)
IntegratePowerTimesSeriesAtOrder[a_, b_, gList_List, targetOrd_, {xmin_, xmax_}, atLower_, atUpper_] :=
  IntegratePowerTimesSeriesAtPower[a, b, 0, gList, targetOrd, {xmin, xmax}, atLower, atUpper];

IntegratePowerTimesSeriesAtPower[a_, b_, epsMinPower_Integer, gList_List, targetPower_Integer,
    {xmin_, xmax_}, atLower_, atUpper_] := Module[
  {epsMaxPower, result, k, m, contrib, gCoeffs, intCoeff, upperVal, lowerVal, maxLogPower},

  epsMaxPower = epsMinPower + Length[gList] - 1;

  (* x^(a+b*eps) = x^a * x^(b*eps) = x^a * exp(b*eps*Logx) = x^a * sum_k (b*Logx)^k/k! * eps^k *)
  (* So at eps^n, we get: sum_{k=0}^{n} (b*Logx)^k/k! * g_{n-k}(x) * x^a *)

  (* The integral of x^a * Logx^k * (series in x) is computed by:
     - Expand the series
     - For each term x^m * Logx^k, integrate to get:
       x^(a+m+1)/(a+m+1) * (Logx^k - k*Logx^(k-1)/(a+m+1) + k*(k-1)*Logx^(k-2)/(a+m+1)^2 - ...)
  *)

  maxLogPower = targetPower - epsMinPower;
  If[maxLogPower < 0, Return[0]];

  result = Sum[
    Module[{logPow = k, logCoeff, gAtOrder, integ, uVal, lVal},
      logCoeff = If[b == 0 && logPow == 0, 1, If[b == 0, 0, b^logPow / logPow!]];
      gAtOrder = If[targetPower - k >= epsMinPower && targetPower - k <= epsMaxPower,
                    gList[[targetPower - k - epsMinPower + 1]],
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
    {k, 0, maxLogPower}
  ];

  result // DiffExp`Utilities`PChop // Expand
];

(* Integrate x^a * Logx^n * (series in x) *)
(* Returns a function/expression that can be evaluated at bounds *)
IntegratePowerLogMonomial[p_, n_, coeff_, tol_] := Module[
  {zeroQ},
  zeroQ[z_] := TrueQ[PossibleZeroQ[z]] || TrueQ[NumericQ[z] && Abs[N[z, 50]] < tol];

  If[coeff === 0 || TrueQ[PossibleZeroQ[coeff]] ||
      (NumericQ[coeff] && Abs[coeff] < 10^-DiffExp`State`FEC[ChopPrecision]),
    0,
    If[zeroQ[p + 1],
      coeff * DiffExp`Symbols`Logx^(n + 1) / (n + 1),
      coeff * DiffExp`Symbols`x^(p + 1) / (p + 1) *
        Sum[(-1)^j * Factorial[n] / Factorial[n - j] *
          DiffExp`Symbols`Logx^(n - j) / (p + 1)^j, {j, 0, n}]
    ]
  ]
];

IntegrateWithLogPower[a_, n_, ser_SeriesData] := Module[
  {nmin, nmax, den, coeffs, result, m, intPow, logTerms, tol},

  nmin = ser[[4]];
  nmax = ser[[5]];
  den = ser[[6]];
  coeffs = ser[[3]];
  tol = DiffExp`State`FEC[RationalizationTolerance];
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

      Sum[
        IntegratePowerLogMonomial[
          p, n + logPow,
          DiffExp`SeriesOps`LogxCoeffNS[coeff, logPow],
          tol
        ],
        {logPow, 0, DiffExp`SeriesOps`MaxLogxPower[coeff]}
      ]
    ],
    {idx, 1, Length[coeffs]}
  ];

  result // Expand
];

IntegrateWithLogPower[a_, n_, 0] := 0;
IntegrateWithLogPower[a_, n_, c_?NumericQ] := IntegrateWithLogPower[a, n,
  SeriesData[DiffExp`Symbols`x, 0, {c}, 0, 1, 1]];

thetaRulesAtPoint[pt_, direction_:Automatic] := Module[
  {tol, sign},
  tol = DiffExp`State`FEC[RationalizationTolerance];
  sign = Which[
    direction =!= Automatic && NumericQ[direction], Sign[N[direction, 50]],
    TrueQ[PossibleZeroQ[pt]] || TrueQ[NumericQ[pt] && Abs[N[pt, 50]] < tol], 1,
    TrueQ[N[pt, 50] < 0], -1,
    True, 1
  ];

  If[sign < 0,
    {DiffExp`Symbols`\[Theta]p -> 0, DiffExp`Symbols`\[Theta]m -> 1},
    {DiffExp`Symbols`\[Theta]p -> 1, DiffExp`Symbols`\[Theta]m -> 0}
  ]
];

branchDirection[pt_, atSingularity_] := Module[
  {},
  Which[
    NumericQ[atSingularity], atSingularity,
    TrueQ[atSingularity === False], Automatic,
    True, Automatic
  ]
];

localSidePhase[side_, power_] := Module[
  {tol = DiffExp`State`FEC[RationalizationTolerance], numericSide},
  numericSide = Quiet[Check[N[side, 50], $Failed]];
  If[numericSide =!= $Failed && NumericQ[numericSide] &&
      TrueQ[numericSide < -tol],
    Exp[-I Pi power],
    1
  ]
];

(* Evaluate the integral expression at a point, handling singularity limits *)
EvaluateIntegralAtPoint[expr_, pt_, a_, atSingularity_] := Module[{val, zeroQ, tol},
  tol = DiffExp`State`FEC[RationalizationTolerance];
  zeroQ = TrueQ[PossibleZeroQ[pt]] || TrueQ[Abs[N[pt, 50]] < tol];

  If[zeroQ,
    (* At x=0 singularity: terms with Logx or negative powers vanish *)
    (* Only finite constant terms survive *)
    Quiet[Limit[
      (expr /. thetaRulesAtPoint[pt, branchDirection[pt, atSingularity]]) /.
        DiffExp`Symbols`Logx -> 0,
      DiffExp`Symbols`x -> 0
    ]] //
      If[NumericQ[#], #, 0] &
    ,
    (* Normal evaluation *)
    val = expr /. thetaRulesAtPoint[pt, branchDirection[pt, atSingularity]] /.
      DiffExp`Symbols`Logx -> Log[DiffExp`Symbols`x] /. DiffExp`Symbols`x -> pt;
    If[NumericQ[val] && TrueQ[Abs[N[val, 50]] > 10^100],
      Module[{ser, padeVal},
        ser = Quiet[
          Check[
            Series[expr, {DiffExp`Symbols`x, 0,
              DiffExp`State`FEC[ExpansionOrder] + 2}],
            $Failed
          ]
        ];
        If[MatchQ[ser, _SeriesData],
          padeVal = Quiet[
            Check[
              DiffExp`Pade`SEval2[DiffExp`Pade`GetPade[ser], pt],
              $Failed
            ]
          ];
          If[NumericQ[padeVal] && TrueQ[Abs[N[padeVal, 50]] < Abs[N[val, 50]]],
            val = padeVal
          ];
        ];
      ];
    ];
    val
  ]
];

(* ============================================================ *)
(* Integration of Decomposed Series *)
(* ============================================================ *)

(* Integrate a full decomposition from DecomposeSingularity *)
IntegrateDecomposition[decomposition_List, {xmin_, xmax_}] := Module[
  {epsOrder, laur},

  If[Length[decomposition] == 0, Return[{}]];

  epsOrder = Length[decomposition[[1]]["g"]] - 1;
  laur = IntegrateDecompositionLaurent[decomposition, {xmin, xmax}, 0];
  LaurentToNonNegativeList[laur, epsOrder]
];

IntegrateDecompositionLaurent[decomposition_List, {xmin_, xmax_}, epsMinPower_Integer:0] := Module[
  {epsMaxPower, result, term, a, b, g},

  If[Length[decomposition] == 0,
    Return[LaurentZero[epsMinPower, epsMinPower - 1]]
  ];

  If[TrueQ[N[xmin, 50] > N[xmax, 50]],
    Return[LaurentScale[-1, IntegrateDecompositionLaurent[decomposition, {xmax, xmin}, epsMinPower]]]
  ];

  epsMaxPower = epsMinPower + Length[decomposition[[1]]["g"]] - 1;
  result = LaurentZero[epsMinPower, epsMaxPower];

  Do[
    a = decomposition[[i]]["a"];
    b = decomposition[[i]]["b"];
    g = decomposition[[i]]["g"];

    result = LaurentAdd[
      result,
      IntegrateSingularTermLaurent[a, b, epsMinPower, g, {xmin, xmax}]
    ];
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
    If[FileExistsQ[seriesData],
       uncompressed = Uncompress[Import[seriesData]];,
       uncompressed = Uncompress[seriesData];
    ];,
    uncompressed = seriesData;
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

  (* Apply Jacobian - handle empty decomposition gracefully *)
  If[ListQ[integralResult] && Length[integralResult] > epsOrder,
    integralResult[[epsOrder + 1]] / Abs[jacobian],
    0  (* Empty decomposition means zero integrand on this segment *)
  ]
];

(* ============================================================ *)
(* Main Entry Points *)
(* ============================================================ *)

(* Integrate piecewise saved data from TransportTo *)
IntegratePiecewiseSaved[savedData2_, {a_, b_}] := Module[
  {savedData, numIntegrals, numEpsOrders, relevantSegments,
   segmentContributions, result, seg, segBounds, overlapBounds,
   uncompressed},

  (* Handle Association output from TransportTo *)
  If[AssociationQ[savedData2],
    If[MissingQ[savedData2["SegmentData"]],
      DiffExp`Utilities`ReportError["IntegratePiecewiseSaved: No segment data available."];
    ];
    savedData = savedData2["SegmentData"],
    If[MatchQ[savedData2, {{_Association, __}, {__}}],
      savedData = savedData2[[2]],
      savedData = savedData2
    ]
  ];

  (* Validate format *)
  If[!ListQ[savedData] || Length[savedData] == 0,
    DiffExp`Utilities`ReportError["IntegratePiecewiseSaved: Invalid saved data format"];
  ];

  (* Get dimensions from first segment *)
  If[StringQ[savedData[[1, 5]]],
    If[FileExistsQ[savedData[[1, 5]]],
       uncompressed = Uncompress[Import[savedData[[1, 5]]]];,
       uncompressed = Uncompress[savedData[[1, 5]]];
    ];,
    uncompressed = savedData[[1, 5]];
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

  (* Handle Association output from TransportTo *)
  If[AssociationQ[savedData2],
    If[MissingQ[savedData2["SegmentData"]],
      DiffExp`Utilities`ReportError["IndefiniteIntegral: No segment data available."];
    ];
    savedData = savedData2["SegmentData"],
    If[MatchQ[savedData2, {{_Association, __}, {__}}],
      savedData = savedData2[[2]],
      savedData = savedData2
    ]
  ];

  (* Get dimensions *)
  If[StringQ[savedData[[1, 5]]],
    If[FileExistsQ[savedData[[1, 5]]],
       uncompressed = Uncompress[Import[savedData[[1, 5]]]];,
       uncompressed = Uncompress[savedData[[1, 5]]];
    ];,
    uncompressed = savedData[[1, 5]];
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

segmentMainExpression[segmentData_List] := Module[
  {currLine = segmentData[[1]], lineRelation = segmentData[[2]]},
  If[AssociationQ[currLine],
    First[Values[currLine]],
    If[Head[currLine] === Rule,
      currLine[[2]],
      lineRelation[[2]]
    ]
  ]
];

segmentActualBounds[segmentData_List] := Module[
  {xLocal = DiffExp`Symbols`x, xMainExpr, localBounds, values},
  xMainExpr = segmentMainExpression[segmentData];
  localBounds = segmentData[[4]];
  values = Quiet[N[(xMainExpr /. xLocal -> #), 80] & /@ localBounds];
  {Min[values], Max[values]}
];

segmentLocalCoordinateForValue[segmentData_List, value_, xMainExpr_:Automatic] := Module[
  {xLocal = DiffExp`Symbols`x, expr, localBounds, values, denom},
  expr = If[xMainExpr === Automatic, segmentMainExpression[segmentData], xMainExpr];
  localBounds = segmentData[[4]];
  values = (expr /. xLocal -> #) & /@ localBounds;
  denom = values[[2]] - values[[1]];
  If[TrueQ[PossibleZeroQ[denom]],
    localBounds[[1]],
    localBounds[[1]] + (value - values[[1]]) *
      (localBounds[[2]] - localBounds[[1]]) / denom
  ]
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
    If[FileExistsQ[seriesData],
       uncompressed = Uncompress[Import[seriesData]];,
       uncompressed = Uncompress[seriesData];
    ];,
    uncompressed = seriesData;
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

(* ============================================================ *)
(* Definite Integral with Prefactors                            *)
(* ============================================================ *)

(*
   Computes ∫_lower^upper x^alpha * (upper-x)^beta * r(x) * f_i(x) dx
   for each integral i and each epsilon order.

   The prefactor x^alpha * (upper-x)^beta * r(x) is handled by:
   1. Converting to local segment coordinates
   2. Separating singular parts (power-law at boundaries) from smooth parts
   3. Series-expanding the smooth parts
   4. Multiplying with the DiffExp series
   5. Adjusting the singular power in the decomposition
   6. Integrating using the standard machinery
*)

DefiniteIntegralWithPrefactor[savedData2_, {lower_, upper_}, prefactorSpec_Association] := Module[
  {savedData, numIntegrals, numEpsOrders, uncompressed,
   alpha, beta, rationalFactor, variable, relevantSegments,
   result, seg, segBounds, overlapBounds},

  (* Extract prefactor specification *)
  alpha = Lookup[prefactorSpec, "PowerAtLower", 0];
  beta = Lookup[prefactorSpec, "PowerAtUpper", 0];
  rationalFactor = Lookup[prefactorSpec, "RationalFactor", 1];
  variable = Lookup[prefactorSpec, "Variable", DiffExp`Symbols`x];

  (* Handle Association output from TransportTo *)
  If[AssociationQ[savedData2],
    If[MissingQ[savedData2["SegmentData"]],
      DiffExp`Utilities`ReportError["DefiniteIntegralWithPrefactor: No segment data available."];
    ];
    savedData = savedData2["SegmentData"],
    savedData = savedData2
  ];

  (* Get dimensions from first segment *)
  If[StringQ[savedData[[1, 5]]],
    If[FileExistsQ[savedData[[1, 5]]],
       uncompressed = Uncompress[Import[savedData[[1, 5]]]];,
       uncompressed = Uncompress[savedData[[1, 5]]];
    ];,
    uncompressed = savedData[[1, 5]];
  ];

  numIntegrals = Length[uncompressed];
  numEpsOrders = Length[uncompressed[[1]]];

  DiffExp`Utilities`PrintInfo["Computing definite integral with prefactor from ", lower // N, " to ", upper // N][1];

  (* Find segments that overlap with [lower, upper] *)
  relevantSegments = Select[savedData,
    Module[{segMin, segMax, actualBounds},
      actualBounds = segmentActualBounds[#];
      segMin = actualBounds[[1]];
      segMax = actualBounds[[2]];
      Not[segMax <= lower || segMin >= upper]
    ] &
  ];

  (* Integrate each relevant segment with prefactors *)
  result = Table[
    Sum[
      seg = relevantSegments[[segIdx]];
      segBounds = segmentActualBounds[seg];

      (* Find overlap with [lower, upper] *)
      overlapBounds = {Max[Min[segBounds], lower], Min[Max[segBounds], upper]};

      If[overlapBounds[[1]] >= overlapBounds[[2]],
        0,
        IntegrateSegmentWithPrefactor[seg, overlapBounds, intIdx, epsOrd,
          alpha, beta, rationalFactor, variable, lower, upper]
      ]
      ,
      {segIdx, Length[relevantSegments]}
    ],
    {intIdx, numIntegrals},
    {epsOrd, 0, numEpsOrders - 1}
  ];

  result
];

DefiniteIntegralWithPrefactorLaurent[savedData2_, {lower_, upper_}, prefactorSpec_Association,
    epsMinPower_Integer:0] := Module[
  {savedData, numIntegrals, numEpsOrders, uncompressed,
   alpha, beta, rationalFactor, variable, relevantSegments,
   result, seg, segBounds, overlapBounds, epsMaxPower},

  alpha = Lookup[prefactorSpec, "PowerAtLower", 0];
  beta = Lookup[prefactorSpec, "PowerAtUpper", 0];
  rationalFactor = Lookup[prefactorSpec, "RationalFactor", 1];
  variable = Lookup[prefactorSpec, "Variable", DiffExp`Symbols`x];

  If[AssociationQ[savedData2],
    If[MissingQ[savedData2["SegmentData"]],
      DiffExp`Utilities`ReportError["DefiniteIntegralWithPrefactorLaurent: No segment data available."];
    ];
    savedData = savedData2["SegmentData"],
    savedData = savedData2
  ];

  If[StringQ[savedData[[1, 5]]],
    If[FileExistsQ[savedData[[1, 5]]],
       uncompressed = Uncompress[Import[savedData[[1, 5]]]];,
       uncompressed = Uncompress[savedData[[1, 5]]];
    ];,
    uncompressed = savedData[[1, 5]];
  ];

  numIntegrals = Length[uncompressed];
  numEpsOrders = Length[uncompressed[[1]]];
  epsMaxPower = epsMinPower + numEpsOrders - 1;

  DiffExp`Utilities`PrintInfo["Computing Laurent definite integral with prefactor from ", lower // N, " to ", upper // N][1];

  relevantSegments = Select[savedData,
    Module[{segMin, segMax, actualBounds},
      actualBounds = segmentActualBounds[#];
      segMin = actualBounds[[1]];
      segMax = actualBounds[[2]];
      Not[segMax <= lower || segMin >= upper]
    ] &
  ];

  result = Table[
    Module[{acc = LaurentZero[epsMinPower, epsMaxPower]},
      Do[
        seg = relevantSegments[[segIdx]];
        segBounds = segmentActualBounds[seg];
        overlapBounds = {Max[Min[segBounds], lower], Min[Max[segBounds], upper]};

        If[overlapBounds[[1]] < overlapBounds[[2]],
          acc = LaurentAdd[
            acc,
            IntegrateSegmentWithPrefactorLaurent[seg, overlapBounds, intIdx,
              epsMinPower, alpha, beta, rationalFactor, variable, lower, upper]
          ];
        ];
        ,
        {segIdx, Length[relevantSegments]}
      ];
      acc
    ],
    {intIdx, numIntegrals}
  ];

  result
];


(* Integrate a single segment with prefactors.
   The full integrand is: (x-lower)^alpha * (upper-x)^beta * r(x) * f(x)
   where f(x) is the DiffExp series at this segment.
*)
IntegrateSegmentWithPrefactor[segmentData_List, {a_, b_}, intIndex_Integer, epsOrder_Integer,
    alpha_, beta_, rationalFactor_, variable_, lower_, upper_] := Module[
  {currLine, lineRelation, mainBounds, localBounds, seriesData,
   jacobian, localA, localB, seriesAtIndex, uncompressed,
   mainMin, mainMax, localMin, localMax, slope,
   atLowerBound, atUpperBound, expansionOrder,
   prefactorSeries, combinedDecomposition, integralResult,
   localSeries, decomposition, xLocal, xMainExpr},

  currLine = segmentData[[1]];
  lineRelation = segmentData[[2]];
  mainBounds = segmentData[[3]];
  localBounds = segmentData[[4]];
  seriesData = segmentData[[5]];

  (* Uncompress if needed *)
  If[StringQ[seriesData],
    If[FileExistsQ[seriesData],
       uncompressed = Uncompress[Import[seriesData]];,
       uncompressed = Uncompress[seriesData];
    ];,
    uncompressed = seriesData;
  ];

  (* Get the series for this integral (all eps orders) *)
  seriesAtIndex = uncompressed[[intIndex]];

  (* Compute coordinate transforms *)
  mainMin = mainBounds[[1]];
  mainMax = mainBounds[[2]];
  localMin = localBounds[[1]];
  localMax = localBounds[[2]];
  xLocal = DiffExp`Symbols`x;
  xMainExpr = DiffExp`Utilities`PChop[Expand[segmentMainExpression[segmentData]]];

  localA = segmentLocalCoordinateForValue[segmentData, a, xMainExpr];
  localB = segmentLocalCoordinateForValue[segmentData, b, xMainExpr];
  jacobian = D[xMainExpr, xLocal];
  jacobian = jacobian /. xLocal -> (localA + localB)/2;

  (* Determine expansion order from the series *)
  expansionOrder = If[MatchQ[seriesAtIndex[[1]], _SeriesData],
    seriesAtIndex[[1]][[5]] - seriesAtIndex[[1]][[4]],
    50
  ];

  (* Check if this segment is at the integration boundaries *)
  atLowerBound = Abs[a - lower] < DiffExp`State`FEC[RationalizationTolerance];
  atUpperBound = Abs[b - upper] < DiffExp`State`FEC[RationalizationTolerance];

  (* Compute the prefactor in local coordinates.
     lineRelation gives: x_main = f(x_local)
     So x_main - lower = f(x_local) - lower, and upper - x_main = upper - f(x_local)
  *)
  Module[{prefactorLower, prefactorUpper, prefactorRational,
          smoothPrefactor, singularPowerLower, singularPowerUpper,
          totalSingularPower, smoothSeries, modifiedDecomp},

    (* Handle power-law prefactors:
       - (x_main - lower)^alpha: if at lower boundary, this is singular
       - (upper - x_main)^beta: if at upper boundary, this is singular
    *)

    (* Compute smooth prefactor: the product of all smooth parts *)
    (* Start with the rational factor r(x_main) *)
    prefactorRational = rationalFactor /. variable -> xMainExpr;

    (* For boundary segments, separate singular and smooth parts *)
    singularPowerLower = 0;
    singularPowerUpper = 0;

    If[atLowerBound && Abs[localA] < DiffExp`State`FEC[RationalizationTolerance],
      (* This segment starts at the lower integration bound *)
      (* (x_main - lower) ~ jacobian * (x_local - localA) ~ jacobian * x_local when localA ~ 0 *)
      (* So (x_main - lower)^alpha is |jacobian|^alpha times the
         branch phase needed to express the physical positive distance
         as x_local^alpha. *)
      singularPowerLower = alpha;
      prefactorLower = Abs[jacobian]^alpha *
        localSidePhase[localB - localA, alpha];
      ,
      (* Not at lower boundary: (x_main - lower)^alpha is smooth, series-expand *)
      prefactorLower = (xMainExpr - lower)^alpha;
    ];

    If[atUpperBound && Abs[localB] < DiffExp`State`FEC[RationalizationTolerance],
      (* This segment ends at the upper integration bound *)
      (* (upper - x_main) ~ |jacobian| * (localB - x_local) ~ |jacobian| * (-x_local) when localB ~ 0 *)
      (* But typically near upper bound, local x goes to 0 at x_main = upper *)
      (* So (upper - x_main) ~ |jacobian| * x_local *)
      singularPowerUpper = beta;
      prefactorUpper = Abs[jacobian]^beta *
        localSidePhase[localA - localB, beta];
      ,
      (* Not at upper boundary: smooth *)
      prefactorUpper = (upper - xMainExpr)^beta;
    ];

    totalSingularPower = singularPowerLower + singularPowerUpper;

    (* Combine all smooth prefactor parts *)
    smoothPrefactor = prefactorLower * prefactorUpper * prefactorRational;

    (* Series-expand the smooth prefactor in local coordinates *)
    smoothSeries = Quiet[
      Series[smoothPrefactor, {xLocal, 0, expansionOrder}] // Normal // Expand
    ];

    (* If Series fails (e.g., non-analytic), try direct expansion *)
    If[!FreeQ[smoothSeries, Series] || !FreeQ[smoothSeries, SeriesData],
      smoothSeries = Normal[Series[smoothPrefactor, {xLocal, 0, expansionOrder}]];
    ];

    (* Decompose the DiffExp series *)
    decomposition = DiffExp`SingularityDecomposition`DecomposeSingularity[seriesAtIndex];

    (* Modify each term in the decomposition by:
       1. Adding the singular power from boundary prefactors
       2. Multiplying g-series by the smooth prefactor series
    *)
    modifiedDecomp = Table[
      Module[{termA, termB, termG, newA, newG, gSeries, smoothSer},
        termA = decomposition[[termIdx]]["a"];
        termB = decomposition[[termIdx]]["b"];
        termG = decomposition[[termIdx]]["g"];

        (* Add singular power from boundary prefactors *)
        newA = termA + totalSingularPower;

        (* Multiply each epsilon-order g-series by the smooth prefactor *)
        newG = Table[
          If[ord <= Length[termG],
            gSeries = termG[[ord]];
            (* Multiply series by smooth prefactor series *)
            If[MatchQ[gSeries, _SeriesData],
              Module[{result},
                smoothSer = Series[smoothPrefactor, {xLocal, gSeries[[2]], gSeries[[5]] - gSeries[[4]]}];
                result = gSeries * smoothSer // Normal;
                (* Convert back to SeriesData *)
                Series[result, {xLocal, gSeries[[2]], Min[gSeries[[5]] - gSeries[[4]], expansionOrder]}]
              ],
              gSeries * (smoothPrefactor /. xLocal -> 0)
            ],
            0
          ],
          {ord, Length[termG]}
        ];

        (* If the rational factor introduced poles (negative starting powers in g),
           absorb them into the exponent a. This ensures IntegrateSingularTerm
           sees a g-series starting at x^0. *)
        Module[{gNmins, minNmin},
          gNmins = Table[
            If[MatchQ[newG[[ord]], _SeriesData],
              newG[[ord]][[4]],  (* nmin of this SeriesData *)
              0
            ], {ord, Length[newG]}
          ];
          minNmin = Min[gNmins];
          If[minNmin < 0,
            newA = newA + minNmin;
            (* Shift each g-series so it starts at x^0 *)
            newG = Table[
              If[MatchQ[newG[[ord]], _SeriesData],
                Module[{s = newG[[ord]]},
                  SeriesData[s[[1]], s[[2]], s[[3]],
                    s[[4]] - minNmin, s[[5]] - minNmin, s[[6]]]
                ],
                newG[[ord]]
              ], {ord, Length[newG]}
            ];
          ];
        ];

        <|"a" -> newA, "b" -> termB, "g" -> newG|>
      ],
      {termIdx, Length[decomposition]}
    ];

    (* Integrate the modified decomposition *)
    integralResult = IntegrateDecomposition[modifiedDecomp, {localA, localB}];

    (* Apply dy/dx_local for integration in the Feynman parameter. *)
    If[ListQ[integralResult] && Length[integralResult] > epsOrder,
      jacobian * integralResult[[epsOrder + 1]],
      0  (* Empty decomposition means zero integrand on this segment *)
    ]
  ]
];

IntegrateSegmentWithPrefactorLaurent[segmentData_List, {a_, b_}, intIndex_Integer,
    epsMinPower_Integer, alpha_, beta_, rationalFactor_, variable_, lower_, upper_] := Module[
  {currLine, lineRelation, mainBounds, localBounds, seriesData,
   jacobian, localA, localB, seriesAtIndex, uncompressed,
   mainMin, mainMax, localMin, localMax, slope,
   atLowerBound, atUpperBound, expansionOrder,
   prefactorSeries, combinedDecomposition, integralResult,
   localSeries, decomposition, xLocal, epsMaxPower, xMainExpr},

  currLine = segmentData[[1]];
  lineRelation = segmentData[[2]];
  mainBounds = segmentData[[3]];
  localBounds = segmentData[[4]];
  seriesData = segmentData[[5]];

  If[StringQ[seriesData],
    If[FileExistsQ[seriesData],
       uncompressed = Uncompress[Import[seriesData]];,
       uncompressed = Uncompress[seriesData];
    ];,
    uncompressed = seriesData;
  ];

  seriesAtIndex = uncompressed[[intIndex]];
  epsMaxPower = epsMinPower + Length[seriesAtIndex] - 1;

  mainMin = mainBounds[[1]];
  mainMax = mainBounds[[2]];
  localMin = localBounds[[1]];
  localMax = localBounds[[2]];
  xLocal = DiffExp`Symbols`x;
  xMainExpr = DiffExp`Utilities`PChop[Expand[segmentMainExpression[segmentData]]];

  localA = segmentLocalCoordinateForValue[segmentData, a, xMainExpr];
  localB = segmentLocalCoordinateForValue[segmentData, b, xMainExpr];
  jacobian = D[xMainExpr, xLocal];
  jacobian = jacobian /. xLocal -> (localA + localB)/2;

  expansionOrder = If[MatchQ[seriesAtIndex[[1]], _SeriesData],
    seriesAtIndex[[1]][[5]] - seriesAtIndex[[1]][[4]],
    50
  ];

  atLowerBound = Abs[a - lower] < DiffExp`State`FEC[RationalizationTolerance];
  atUpperBound = Abs[b - upper] < DiffExp`State`FEC[RationalizationTolerance];

  Module[{prefactorLower, prefactorUpper, prefactorRational,
          smoothPrefactor, singularPowerLower, singularPowerUpper,
          totalSingularPower, smoothSeries, modifiedDecomp},

    prefactorRational = rationalFactor /. variable -> xMainExpr;

    singularPowerLower = 0;
    singularPowerUpper = 0;

    If[atLowerBound && Abs[localA] < DiffExp`State`FEC[RationalizationTolerance],
      singularPowerLower = alpha;
      prefactorLower = Abs[jacobian]^alpha *
        localSidePhase[localB - localA, alpha];
      ,
      prefactorLower = (xMainExpr - lower)^alpha;
    ];

    If[atUpperBound && Abs[localB] < DiffExp`State`FEC[RationalizationTolerance],
      singularPowerUpper = beta;
      prefactorUpper = Abs[jacobian]^beta *
        localSidePhase[localA - localB, beta];
      ,
      prefactorUpper = (upper - xMainExpr)^beta;
    ];

    totalSingularPower = singularPowerLower + singularPowerUpper;
    smoothPrefactor = prefactorLower * prefactorUpper * prefactorRational;

    smoothSeries = Quiet[
      Series[smoothPrefactor, {xLocal, 0, expansionOrder}] // Normal // Expand
    ];

    If[!FreeQ[smoothSeries, Series] || !FreeQ[smoothSeries, SeriesData],
      smoothSeries = Normal[Series[smoothPrefactor, {xLocal, 0, expansionOrder}]];
    ];

    decomposition = DiffExp`SingularityDecomposition`DecomposeSingularity[seriesAtIndex];

    modifiedDecomp = Table[
      Module[{termA, termB, termG, newA, newG, gSeries, smoothSer},
        termA = decomposition[[termIdx]]["a"];
        termB = decomposition[[termIdx]]["b"];
        termG = decomposition[[termIdx]]["g"];

        newA = termA + totalSingularPower;

        newG = Table[
          If[ord <= Length[termG],
            gSeries = termG[[ord]];
            If[MatchQ[gSeries, _SeriesData],
              Module[{result},
                smoothSer = Series[smoothPrefactor, {xLocal, gSeries[[2]], gSeries[[5]] - gSeries[[4]]}];
                result = gSeries * smoothSer // Normal;
                Series[result, {xLocal, gSeries[[2]], Min[gSeries[[5]] - gSeries[[4]], expansionOrder]}]
              ],
              gSeries * (smoothPrefactor /. xLocal -> 0)
            ],
            0
          ],
          {ord, Length[termG]}
        ];

        Module[{gNmins, minNmin},
          gNmins = Table[
            If[MatchQ[newG[[ord]], _SeriesData],
              newG[[ord]][[4]],
              0
            ], {ord, Length[newG]}
          ];
          minNmin = Min[gNmins];
          If[minNmin < 0,
            newA = newA + minNmin;
            newG = Table[
              If[MatchQ[newG[[ord]], _SeriesData],
                Module[{s = newG[[ord]]},
                  SeriesData[s[[1]], s[[2]], s[[3]],
                    s[[4]] - minNmin, s[[5]] - minNmin, s[[6]]]
                ],
                newG[[ord]]
              ], {ord, Length[newG]}
            ];
          ];
        ];

        <|"a" -> newA, "b" -> termB, "g" -> newG|>
      ],
      {termIdx, Length[decomposition]}
    ];

    integralResult = IntegrateDecompositionLaurent[
      modifiedDecomp, {localA, localB}, epsMinPower
    ];

    If[AssociationQ[integralResult],
      LaurentScale[jacobian, integralResult],
      LaurentZero[epsMinPower, epsMaxPower]
    ]
  ]
];


End[];

EndPackage[];
