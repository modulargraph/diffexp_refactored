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

EvaluateEndpointLimitSectors::usage = "EvaluateEndpointLimitSectors[seriesList, direction] evaluates lim_{x->0} of a transport series tower (one entry per epsilon order), resolving residual x^(a + b_i*eps) endpoint sectors so the analytic-regularization prescription (drop sectors with b != 0, even when a < 0) is applied per sector instead of to DecomposeSingularity's single collapsed exponent. Returns the per-epsilon-order list of limit values.";

FitResidualEndpointSectors::usage = "FitResidualEndpointSectors[coeffList, branchRules] resolves the epsilon tower of a fixed local power (one polynomial in Logx per epsilon offset) into residual x^(r*eps) sectors. Returns <|\"Sectors\" -> {<|\"ResidualB\", \"Coefficients\"|>..}, \"SalvageOffsets\", \"SalvageExact\", \"Resolved\"|> or $Failed.";

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

$laurentDumpCounter = 0;

activeNumericPrecision[] := Module[{precision},
  precision = Quiet[Check[DiffExp`State`FEWorkingPrecision, 500]];
  If[IntegerQ[precision] && precision > 0, precision, 500]
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

branchReplaceAtActivePrecision[expr_, rules_] :=
  SetPrecision[expr, activeNumericPrecision[]] /. rules;

seriesAtActivePrecision[expr_] := DiffExp`SeriesOps`SApply[
  SetPrecision[#, activeNumericPrecision[]] &,
  expr
];

seriesMultiplyByX[ser_SeriesData] :=
  SeriesData[ser[[1]], ser[[2]], ser[[3]], ser[[4]] + ser[[6]],
    ser[[5]] + ser[[6]], ser[[6]]];
seriesMultiplyByX[expr_] := DiffExp`Symbols`x * expr;

linearSeriesCombine[terms_List, chop_:True] := Module[
  {activeTerms, seriesTerms, first, var, center, den, minPow, maxPow, coeffs,
   outCoeffs},

  activeTerms = Select[terms, #[[1]] =!= 0 && #[[2]] =!= 0 &];
  If[activeTerms === {}, Return[0, Module]];
  seriesTerms = Select[activeTerms, MatchQ[#[[2]], _SeriesData] &];
  If[seriesTerms === {},
    Return[DiffExp`Utilities`PChop[Total[#[[1]] #[[2]] & /@ activeTerms]], Module]
  ];

  If[!AllTrue[activeTerms, MatchQ[#[[2]], _SeriesData] || TrueQ[FreeQ[#[[2]], DiffExp`Symbols`x]] &],
    Return[DiffExp`SeriesOps`SExpand[Total[#[[1]] #[[2]] & /@ activeTerms]], Module]
  ];

  first = seriesTerms[[1, 2]];
  {var, center, den} = {first[[1]], first[[2]], first[[6]]};
  If[!AllTrue[
      seriesTerms[[All, 2]],
      #[[1]] === var && #[[2]] === center && #[[6]] === den &
    ],
    Return[DiffExp`SeriesOps`SExpand[Total[#[[1]] #[[2]] & /@ activeTerms]], Module]
  ];

  minPow = Min[
    Join[
      seriesTerms[[All, 2, 4]],
      If[Length[seriesTerms] < Length[activeTerms], {0}, {}]
    ]
  ];
  maxPow = Max[
    Join[
      seriesTerms[[All, 2, 5]],
      If[Length[seriesTerms] < Length[activeTerms], {1}, {}]
    ]
  ];
  coeffs = ConstantArray[0, Max[0, maxPow - minPow]];

  Do[
    Module[{coef = term[[1]], ser = term[[2]], start, serCoeffs},
      If[MatchQ[ser, _SeriesData],
        start = ser[[4]] - minPow + 1;
        serCoeffs = ser[[3]];
        Do[
          coeffs[[start + idx - 1]] += coef * serCoeffs[[idx]],
          {idx, Length[serCoeffs]}
        ],
        coeffs[[1 - minPow]] += coef * ser
      ];
    ],
    {term, activeTerms}
  ];

  outCoeffs = If[TrueQ[chop], DiffExp`Utilities`PChop /@ coeffs, coeffs];
  SeriesData[var, center, outCoeffs, minPow, maxPow, den]
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

LaurentTrim[laur_Association] := Module[
  {minPower = laur["MinPower"], coeffs = laur["Coefficients"], nonzero, first, last},
  nonzero = Select[Range[Length[coeffs]], !EffectiveZeroExprQ[coeffs[[#]]] &];
  If[nonzero === {}, Return[LaurentZero[0, 0], Module]];
  first = First[nonzero];
  last = Last[nonzero];
  <|
    "MinPower" -> minPower + first - 1,
    "Coefficients" -> coeffs[[first ;; last]]
  |>
];

LaurentToNonNegativeList[laur_Association, epsOrder_Integer] :=
  Table[LaurentCoeff[laur, p], {p, 0, epsOrder}];

maybeDumpLaurentDefiniteIntegral[payload_Association] := Module[
  {dir, targetRaw, target, file, dumpPayload},

  dir = Environment["DIFFEXP_DUMP_LAURENT_DIR"];
  If[!StringQ[dir] || dir === "", Return[Null, Module]];

  $laurentDumpCounter++;
  targetRaw = Environment["DIFFEXP_DUMP_LAURENT_ABORT_AFTER"];
  target = Quiet[Check[ToExpression[targetRaw], None]];
  If[!DirectoryQ[dir],
    CreateDirectory[dir, CreateIntermediateDirectories -> True]
  ];

  dumpPayload = Join[
    <|"CallIndex" -> $laurentDumpCounter|>,
    payload
  ];
  file = FileNameJoin[{
    dir,
    "laurent_integral_" <>
      StringPadLeft[ToString[$laurentDumpCounter], 4, "0"] <> ".m"
  }];
  laurentIntegralDump = dumpPayload;
  Save[file, laurentIntegralDump];
  Clear[laurentIntegralDump];
  Print["DiffExp: dumped Laurent definite integral input to ", file];

  If[IntegerQ[target] && target === $laurentDumpCounter,
    Print["DiffExp: aborting after requested Laurent dump ", target];
    Abort[];
  ];
];

NumericNegativeQ[z_, tol_:Automatic] := Module[
  {eps = If[tol === Automatic, DiffExp`State`FEC[RationalizationTolerance], tol],
   nz},
  nz = realNumericAtActivePrecision[z, eps];
  TrueQ[nz < -eps]
];

NumericZeroQ[z_, tol_:Automatic] := Module[
  {eps = If[tol === Automatic, DiffExp`State`FEC[RationalizationTolerance], tol],
   nz},
  TrueQ[PossibleZeroQ[z]] ||
    (nz = realNumericAtActivePrecision[z, eps];
     nz =!= $Failed && NumericQ[nz] && TrueQ[Abs[nz] < eps])
];

IntegerPower[p_] := Module[
  {tol = DiffExp`State`FEC[RationalizationTolerance], np, rounded, exact},
  If[IntegerQ[p], Return[p, Module]];
  np = realNumericAtActivePrecision[p, tol];
  If[np =!= $Failed && NumericQ[np] &&
      TrueQ[Abs[np - Round[np]] < tol],
    rounded = Round[np];
    exact = Rationalize[rounded, 0];
    If[IntegerQ[exact], exact, rounded],
    p
  ]
];

EffectiveZeroExprQ[expr_, tol_:Automatic] := Module[
  {eps = If[tol === Automatic, DiffExp`State`FEC[RationalizationTolerance], tol],
   expanded, numeric, thetaBranches, branchValues},
  expanded = DiffExp`Utilities`PChop[Expand[expr]];
  If[TrueQ[PossibleZeroQ[expanded]], Return[True, Module]];
  thetaBranches = {
    {DiffExp`Symbols`\[Theta]p -> 1, DiffExp`Symbols`\[Theta]m -> 0},
    {DiffExp`Symbols`\[Theta]p -> 0, DiffExp`Symbols`\[Theta]m -> 1}
  };
  If[!FreeQ[expanded, DiffExp`Symbols`\[Theta]p | DiffExp`Symbols`\[Theta]m],
    branchValues = numericAtActivePrecision[
      branchReplaceAtActivePrecision[expanded, #]
    ] & /@ thetaBranches;
    If[AllTrue[
        branchValues,
        (# =!= $Failed && NumericQ[#] && Abs[#] < eps) &
      ],
      Return[True, Module]
    ];
  ];
  numeric = numericAtActivePrecision[expanded];
  NumericQ[numeric] && Abs[numeric] < eps
];

logPolynomialDerivativeCoefficient[coeff_] :=
  D[coeff, DiffExp`Symbols`Logx];

KnownSeriesDerivative[ser_SeriesData, var_] := Module[
  {order, nmin, nmax, den, coeffs, newCoeffs},
  If[ser[[1]] === var && NumericZeroQ[ser[[2]]],
    nmin = ser[[4]];
    nmax = ser[[5]];
    den = ser[[6]];
    coeffs = ser[[3]];
    newCoeffs = Table[
      Module[{power = (nmin + idx - 1)/den, coeff = coeffs[[idx]]},
        power * coeff + logPolynomialDerivativeCoefficient[coeff]
      ],
      {idx, Length[coeffs]}
    ];
    SeriesData[ser[[1]], ser[[2]], newCoeffs, nmin - den, nmax - den, den],
    order = Max[0, Ceiling[ser[[5]] / ser[[6]]]];
    Quiet[
      Series[
        SetPrecision[
          DiffExp`SeriesOps`SD[Normal[ser], var],
          activeNumericPrecision[]
        ],
        {var, ser[[2]], order}
      ]
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
  {epsOrder, transformedMaxOrder, newGList, newEpsMin, transformedG, n, tol,
   zeroQ, onePlusA, twoPlusA, bPrecise, cPrecise},

  epsOrder = Length[gList] - 1;
  tol = DiffExp`State`FEC[RationalizationTolerance];
  zeroQ[z_] := TrueQ[PossibleZeroQ[z]] ||
    Module[{nz = realNumericAtActivePrecision[z, tol]},
      nz =!= $Failed && NumericQ[nz] && TrueQ[Abs[nz] < tol]
    ];
  onePlusA = DiffExp`Utilities`PChop[
    SetPrecision[1 + a, activeNumericPrecision[]]
  ];
  twoPlusA = DiffExp`Utilities`PChop[
    SetPrecision[2 + a, activeNumericPrecision[]]
  ];
  bPrecise = SetPrecision[b, activeNumericPrecision[]];
  cPrecise = SetPrecision[c, activeNumericPrecision[]];
  transformedMaxOrder = epsOrder + If[zeroQ[onePlusA], 1, 0];

  (* First compute the transformed g: (2+a+b*eps)/c * g - (1-x/c) * g'
     This doesn't change the eps structure, just transforms each coefficient *)
  transformedG = Table[
    Module[{gAtOrd, gAtOrdMinus1, gpAtOrd},
      gAtOrd = seriesAtActivePrecision[
        If[ord <= epsOrder, gList[[ord + 1]], 0]
      ];
      gAtOrdMinus1 = seriesAtActivePrecision[
        If[ord > 0 && ord - 1 <= epsOrder, gList[[ord]], 0]
      ];
      gpAtOrd = seriesAtActivePrecision[
        KnownSeriesDerivative[gAtOrd, DiffExp`Symbols`x]
      ];

      (* (2+a)/c * g_ord + b/c * g_{ord-1} - (1-x/c) * g'_ord *)
      linearSeriesCombine[
        {
          {twoPlusA / cPrecise, gAtOrd},
          {If[ord > 0, bPrecise / cPrecise, 0], gAtOrdMinus1},
          {-1, gpAtOrd},
          {1 / cPrecise, seriesMultiplyByX[gpAtOrd]}
        },
        False
      ]
    ],
    {ord, 0, transformedMaxOrder}
  ];

  (* Now handle the prefactor 1/(1+a+b*eps) *)
  If[zeroQ[onePlusA],
    (* Special case: a = -1, prefactor is 1/(b*eps) = (1/b) * eps^{-1} *)
    (* This shifts eps powers down by 1 and multiplies by 1/b *)
    newEpsMin = epsMinPower - 1;
    newGList = (1/bPrecise) * # & /@ transformedG;
    ,
    (* General case: 1/(1+a+b*eps) expands as a geometric series.
       Combine the finite epsilon convolution coefficient-wise, avoiding both
       global SExpand and recursively nested partial sums. *)
    newEpsMin = epsMinPower;
    newGList = Table[
      linearSeriesCombine[
        Table[
          {(-bPrecise / onePlusA)^k / onePlusA,
            transformedG[[n - k + 1]]},
          {k, 0, n}
        ],
        False
      ],
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

RegularizeIntegrand[a_, b_, gList_List, c_] :=
  RegularizeIntegrand[a, b, 0, gList, c];

RegularizeIntegrand[a_, b_, epsMinPower_Integer, gList_List, c_] := Module[
  {currentA = a, currentB = b, currentEpsMin = epsMinPower, currentG = gList, maxIter = 50, iter = 0, result},

  While[NumericNegativeQ[currentA] && iter < maxIter,
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
  {epsOrder, regA, regB, regEpsMin, regG, result,
   atLowerSingularity, atUpperSingularity, tol, zeroQ, epsMaxPower,
   branchDir, branchPoint, effectiveB, regulatorB},

  epsOrder = Length[gList] - 1;
  epsMaxPower = epsMinPower + epsOrder;
  tol = DiffExp`State`FEC[RationalizationTolerance];
  zeroQ[z_] := TrueQ[PossibleZeroQ[z]] ||
    TrueQ[Abs[numericAtActivePrecision[z]] < tol];

  (* The regularization formula is written with the singular point x=0
     as the lower endpoint.  If zero is the upper endpoint, flip the
     interval; if zero lies inside the interval, split into two regulated
     endpoint integrals. *)
  If[zeroQ[xmax] && !zeroQ[xmin],
    Return[LaurentScale[-1, IntegrateSingularTermLaurent[a, b, epsMinPower, gList, {0, xmin}]]]
  ];
  If[!zeroQ[xmin] && !zeroQ[xmax] &&
     TrueQ[realNumericAtActivePrecision[xmin, tol] < 0 <
       realNumericAtActivePrecision[xmax, tol]],
    Return[
      LaurentAdd[
        IntegrateSingularTermLaurent[a, b, epsMinPower, gList, {xmin, 0}],
        IntegrateSingularTermLaurent[a, b, epsMinPower, gList, {0, xmax}]
      ]
    ]
  ];

  (* Check if boundaries are at singularities *)
  atLowerSingularity = If[zeroQ[xmin],
    If[zeroQ[xmax], 1, Sign[realNumericAtActivePrecision[xmax, tol]]],
    False
  ];
  atUpperSingularity = If[zeroQ[xmax],
    If[zeroQ[xmin], 1, Sign[realNumericAtActivePrecision[xmin, tol]]],
    False
  ];

  branchDir = Which[
    NumericQ[atLowerSingularity], atLowerSingularity,
    NumericQ[atUpperSingularity], atUpperSingularity,
    !zeroQ[xmin], Sign[realNumericAtActivePrecision[xmin, tol]],
    !zeroQ[xmax], Sign[realNumericAtActivePrecision[xmax, tol]],
    True, 1
  ];
  branchPoint = If[zeroQ[xmin], xmax, xmin];
  effectiveB = DiffExp`Utilities`PChop[
    Expand[
      branchReplaceAtActivePrecision[
        b,
        thetaRulesAtPoint[branchPoint, branchDir]
      ]
    ]
  ];
  (* A vanishing extracted exponent with a divergent power needs a formal
     unit regulator, but only when the coefficients carry no residual Logx
     towers: towers mean the true sector exponents survive inside gList
     (e.g. perfectly cancelling leading weights), and shifting the basis
     exponent by one would displace every recovered x^(b eps) sector. *)
  regulatorB = If[
    zeroQ[effectiveB] && atLowerSingularity =!= False &&
      TrueQ[realNumericAtActivePrecision[a, tol] <= -1 + tol] &&
      FreeQ[gList, DiffExp`Symbols`Logx],
    1,
    effectiveB
  ];

  (* If a < 0 and b != 0 at the endpoint x=0, we need regularization.
     Away from x=0 the ordinary termwise antiderivative is finite; applying
     the endpoint regularization formula there creates a spurious boundary term.
     Residual Logx towers also take this path even when the extracted b
     vanishes: the towers carry the true x^(b_i eps) sector exponents, which
     the endpoint-sector recovery inside the subtraction formula restores. *)
  If[NumericNegativeQ[a, tol] && atLowerSingularity =!= False &&
      (!zeroQ[regulatorB] || !FreeQ[gList, DiffExp`Symbols`Logx]),
    Return[
      IntegrateAnalyticRegularizedBySubtractionLaurent[
        a, regulatorB, epsMinPower, gList, {xmin, xmax},
        atLowerSingularity, atUpperSingularity
      ]
    ];
    ,
    regA = a; regB = regulatorB; regEpsMin = epsMinPower; regG = gList;
  ];

  regEpsMin = IntegerPower[regEpsMin];
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
    IntegratePowerTimesSeriesAtPower[regA, regB, regEpsMin, regG, ord,
      {xmin, xmax}, atLowerSingularity, atUpperSingularity],
    {ord, regEpsMin, epsMaxPower}
  ];

  <|"MinPower" -> regEpsMin, "Coefficients" -> result|>
];

IntegratePowerTimesSeriesAtPower[a_, b_, epsMinPower_, gList_List, targetPower_,
    {xmin_, xmax_}, atLower_, atUpper_] := Module[
  {epsMaxPower, result, maxLogPower, epsMin, target, contribs, tol, zeroQ,
   branchDir, branchPoint, branchRules},

  epsMin = IntegerPower[epsMinPower];
  target = IntegerPower[targetPower];
  If[!IntegerQ[epsMin] || !IntegerQ[target],
    DiffExp`Utilities`PrintWarning[
      "IntegratePowerTimesSeriesAtPower: non-integer epsilon power bounds ",
      {epsMinPower, targetPower}, "; dropping contribution."
    ];
    Return[0];
  ];

  epsMaxPower = epsMin + Length[gList] - 1;
  tol = DiffExp`State`FEC[RationalizationTolerance];
  zeroQ[z_] := TrueQ[PossibleZeroQ[z]] ||
    TrueQ[Abs[numericAtActivePrecision[z]] < tol];
  branchDir = Which[
    NumericQ[atLower], atLower,
    NumericQ[atUpper], atUpper,
    !zeroQ[xmin], Sign[realNumericAtActivePrecision[xmin, tol]],
    !zeroQ[xmax], Sign[realNumericAtActivePrecision[xmax, tol]],
    True, 1
  ];
  branchPoint = If[zeroQ[xmin], xmax, xmin];
  branchRules = thetaRulesAtPoint[branchPoint, branchDir];
  maxLogPower = IntegerPower[target - epsMin];
  If[!IntegerQ[maxLogPower],
    DiffExp`Utilities`PrintWarning[
      "IntegratePowerTimesSeriesAtPower: non-integer log bound ", maxLogPower,
      "; dropping contribution."
    ];
    Return[0];
  ];
  If[maxLogPower < 0, Return[0]];

  contribs = Table[
    With[{logPow = k},
      Module[{logCoeff, gAtOrder, integ, uVal, lVal},
        logCoeff = If[logPow == 0,
          1,
          DiffExp`Utilities`PChop[
            Expand[
              branchReplaceAtActivePrecision[b, branchRules]^logPow / logPow!
            ]
          ]
        ];
        If[EffectiveZeroExprQ[logCoeff, tol], Return[0, Module]];

        gAtOrder = If[target - k >= epsMin && target - k <= epsMaxPower,
                      gList[[target - k - epsMin + 1]],
                      0];
        gAtOrder = DiffExp`SeriesOps`SApply[
          Chop[
            DiffExp`Utilities`PChop[
              Expand[branchReplaceAtActivePrecision[#, branchRules]]
            ],
            tol
          ] &,
          gAtOrder
        ];

        If[gAtOrder === 0,
          0,
          (* Integrate x^a * Logx^logPow * gAtOrder *)
          integ = IntegrateWithLogPower[a, logPow, gAtOrder];

          (* Evaluate at bounds *)
          uVal = EvaluateIntegralAtPoint[integ, xmax, a, atUpper];
          lVal = EvaluateIntegralAtPoint[integ, xmin, a, atLower];

          logCoeff * (uVal - lVal)
        ]
      ]
    ],
    {k, 0, maxLogPower}
  ];
  result = Total[contribs];

  result // DiffExp`Utilities`PChop // Expand
];

zeroPowerSafe[base_, power_Integer] := If[power === 0, 1, base^power];

IntegrateAnalyticRegularizedByIBPLaurent[a_, b_, epsMinPower_Integer, gList_List,
    {xmin_, xmax_}, atLower_, atUpper_] := Module[
  {tol, zeroQ, regularizationPoleShiftCount, extraRegOrders, regInputG,
   regResult, regA, regB, regEpsMin, regG, requestedEpsMaxPower, result},

  tol = DiffExp`State`FEC[RationalizationTolerance];
  zeroQ[z_] := TrueQ[PossibleZeroQ[z]] ||
    TrueQ[Abs[numericAtActivePrecision[z]] < tol];
  regularizationPoleShiftCount[startA_] := Module[
    {currentA = startA, iter = 0, maxIter = 50, count = 0},
    While[NumericNegativeQ[currentA, tol] && iter < maxIter,
      If[zeroQ[1 + currentA], count++];
      currentA = currentA + 1;
      iter++;
    ];
    count
  ];

  requestedEpsMaxPower = epsMinPower + Length[gList] - 1;
  extraRegOrders = regularizationPoleShiftCount[a];
  regInputG = Join[gList, ConstantArray[0, extraRegOrders]];
  regResult = RegularizeIntegrand[a, b, epsMinPower, regInputG, xmax];
  regA = regResult[[1]];
  regB = regResult[[2]];
  regEpsMin = IntegerPower[regResult[[3]]];
  regG = regResult[[4]];

  result = Table[
    IntegratePowerTimesSeriesAtPower[
      regA, regB, regEpsMin, regG, ord,
      {xmin, xmax}, atLower, atUpper
    ],
    {ord, regEpsMin, requestedEpsMaxPower}
  ];

  <|"MinPower" -> regEpsMin, "Coefficients" -> result|>
];

(* Recover residual endpoint sectors sum_i w_i(eps) x^(r_i eps) from the
   epsilon tower of a fixed local power.  After DecomposeSingularity has
   factored out x^(a + branchB*eps), genuinely multi-sector data leaves a
   residual tower whose Logx structure encodes the remaining exponents:

     T(n, k) := [Logx^k] coeff(epsilon offset n)
             ==  Sum_i c_{i, n-k} r_i^k / k!

   with absolute offsets n counted from the first gList entry, so each
   (n, k) pair probes a single weight order q = n - k.  The shifted
   diagonals m_k = k! T(q0 + k, k) = Sum_i c_{i, q0} r_i^k form a Prony
   system for the roots r_i for any reference order q0 whose weight slice
   is nonzero; q0 is scanned because leading weight slices can cancel
   (e.g. the s0 = 0 banana {1,0,2,1} tower).  Fits are only trusted on
   the leading run of offsets where the full reconstructed Logx content
   matches the data.  Beyond that run (truncated or partially combined
   top epsilon orders) the non-log content is salvaged against the plain
   branch exponent under the usual truncation-boundary convention, and
   dropped Logx content is reported. *)
FitResidualEndpointSectors[coeffList_List, branchRules_List] := Module[
  {resolved, maxOffset, firstVisible, maxLogAtFirst, logCoeffAt,
   relTol, relZeroQ, q0Candidates, candidateRootSets, ladder,
   m0, m1, m2, m3, det, e1, e2, disc, twoRoots,
   evaluateFit, best, fitResult, sectorCount, roots, sectorCoeffs,
   usableQMax, salvageOffsets, salvageExact, tol},

  tol = DiffExp`State`FEC[RationalizationTolerance];
  resolved = DiffExp`Utilities`PChop[
    Expand[branchReplaceAtActivePrecision[#, branchRules]]
  ] & /@ coeffList;
  maxOffset = Length[resolved] - 1;
  firstVisible = SelectFirst[
    Range[0, maxOffset],
    !EffectiveZeroExprQ[resolved[[# + 1]], tol] &,
    Missing["None"]
  ];
  If[MissingQ[firstVisible],
    Return[<|"Sectors" -> {}, "SalvageOffsets" -> {},
      "SalvageExact" -> True, "Resolved" -> resolved|>, Module]
  ];

  logCoeffAt[offset_Integer, logPower_Integer] := If[
    offset < 0 || offset > maxOffset,
    0,
    DiffExp`SeriesOps`LogxCoeffNS[resolved[[offset + 1]], logPower]
  ];

  (* Relative consistency test.  Tower data carries the active working
     precision, so genuine sector-model violations (truncated or
     partially combined data, resonant logs) sit many orders of magnitude
     above the noise floor of a valid fit. *)
  relTol = Max[tol^(1/3), 10^-12];
  relZeroQ[diff_, scale_] := Module[
    {d = numericAtActivePrecision[diff], s = numericAtActivePrecision[scale]},
    If[!NumericQ[d] || !NumericQ[s], Return[False, Module]];
    Abs[d] <= relTol * (1 + Abs[s])
  ];

  (* The weight slices may start before the first visible offset when the
     leading weights cancel: content at (firstVisible, log power k) probes
     weight order firstVisible - k. *)
  maxLogAtFirst = Max[0, Max[Join[{0},
    Select[
      DiffExp`SeriesOps`LogxPowerRange[resolved[[firstVisible + 1]]],
      !EffectiveZeroExprQ[logCoeffAt[firstVisible, #], tol] &
    ]
  ]]];
  q0Candidates = Range[Max[0, firstVisible - maxLogAtFirst], firstVisible];

  (* Candidate root sets from the shifted moment diagonals, smallest
     sector count first for each usable reference order. *)
  candidateRootSets = {};
  Do[
    ladder = Table[
      DiffExp`Utilities`PChop[
        Expand[Factorial[k] * logCoeffAt[q0 + k, k]]
      ],
      {k, 0, Min[3, maxOffset - q0]}
    ];
    If[Length[ladder] < 2 ||
        AllTrue[ladder, EffectiveZeroExprQ[#, tol] &],
      Continue[]
    ];
    {m0, m1} = ladder[[1 ;; 2]];
    m2 = If[Length[ladder] >= 3, ladder[[3]], 0];
    m3 = If[Length[ladder] >= 4, ladder[[4]], 0];
    If[!EffectiveZeroExprQ[m0, tol],
      AppendTo[candidateRootSets,
        <|"Roots" -> {If[EffectiveZeroExprQ[m1, tol], 0,
            DiffExp`Utilities`PChop[Expand[m1 / m0]]
          ]},
          "ReferenceOrder" -> q0|>
      ];
    ];
    If[Length[ladder] >= 4,
      det = DiffExp`Utilities`PChop[Expand[m0 * m2 - m1^2]];
      If[!EffectiveZeroExprQ[det, tol],
        e1 = DiffExp`Utilities`PChop[Expand[(m0 * m3 - m1 * m2) / det]];
        e2 = DiffExp`Utilities`PChop[Expand[(m1 * m3 - m2^2) / det]];
        disc = DiffExp`Utilities`PChop[Expand[e1^2 - 4 * e2]];
        twoRoots = (DiffExp`Utilities`PChop[Expand[#]] &) /@
          {(e1 + Sqrt[disc])/2, (e1 - Sqrt[disc])/2};
        If[!EffectiveZeroExprQ[twoRoots[[1]] - twoRoots[[2]], tol],
          AppendTo[candidateRootSets,
            <|"Roots" -> twoRoots, "ReferenceOrder" -> q0|>
          ];
        ];
      ];
    ],
    {q0, q0Candidates}
  ];
  If[candidateRootSets === {}, Return[$Failed, Module]];

  (* Solve the weight tower for a candidate root set and count how many
     leading offsets the full Logx reconstruction explains. *)
  evaluateFit[testRoots_List] := Module[
    {count = Length[testRoots], coeffs, qSolveMax, w0, w1, denom,
     validOffsets, q, predicted, dataVal, scaleVal, failed},

    qSolveMax = maxOffset - (count - 1);
    coeffs = ConstantArray[0, {count, maxOffset + 1}];
    Do[
      If[count === 1,
        coeffs[[1, q + 1]] = DiffExp`Utilities`PChop[
          Expand[logCoeffAt[q, 0]]
        ],
        w0 = logCoeffAt[q, 0];
        w1 = logCoeffAt[q + 1, 1];
        denom = testRoots[[1]] - testRoots[[2]];
        coeffs[[1, q + 1]] = DiffExp`Utilities`PChop[
          Expand[(w1 - w0 * testRoots[[2]]) / denom]
        ];
        coeffs[[2, q + 1]] = DiffExp`Utilities`PChop[
          Expand[w0 - coeffs[[1, q + 1]]]
        ];
      ],
      {q, 0, qSolveMax}
    ];

    (* Validate every available (offset, log power) pair whose weight
       order has been solved; stop at the first failing offset. *)
    validOffsets = maxOffset + 1;
    Do[
      failed = False;
      Do[
        Module[{qq = n - k},
          If[qq >= 0 && qq <= qSolveMax,
            dataVal = Factorial[k] * logCoeffAt[n, k];
            predicted = Total[Table[
              coeffs[[i, qq + 1]] * zeroPowerSafe[testRoots[[i]], k],
              {i, count}
            ]];
            scaleVal = Max[
              Abs[numericAtActivePrecision[dataVal]],
              Abs[numericAtActivePrecision[predicted]]
            ];
            If[!relZeroQ[dataVal - predicted, scaleVal],
              failed = True;
            ];
          ];
        ];
        If[failed, Break[]],
        {k, 0, n}
      ];
      If[failed,
        validOffsets = n;
        Break[];
      ],
      {n, 0, maxOffset}
    ];

    (* When validation fails at some offset, the weight orders solved
       from the by-construction rows touching the failing region may have
       absorbed pollution that only became visible one offset later (for
       example missing low-log homogeneous content in truncated upstream
       data), so retreat one extra order beyond the validated run. *)
    <|
      "Roots" -> testRoots,
      "SectorCount" -> count,
      "Coefficients" -> coeffs,
      "ValidOffsets" -> validOffsets,
      "UsableQMax" -> Min[
        qSolveMax,
        validOffsets - count - If[validOffsets <= maxOffset, 1, 0]
      ]
    |>
  ];

  best = Missing["None"];
  Do[
    fitResult = evaluateFit[candidate["Roots"]];
    (* The moment ladder at the candidate's reference order must lie
       inside the validated run, otherwise the roots themselves are not
       trustworthy. *)
    If[fitResult["ValidOffsets"] >=
        candidate["ReferenceOrder"] + 2 * fitResult["SectorCount"],
      If[MissingQ[best] ||
          fitResult["ValidOffsets"] > best["ValidOffsets"],
        best = fitResult;
      ];
      If[fitResult["ValidOffsets"] >= maxOffset + 1, Break[]];
    ],
    {candidate, candidateRootSets}
  ];
  If[MissingQ[best] || best["UsableQMax"] < 0, Return[$Failed, Module]];

  sectorCount = best["SectorCount"];
  roots = best["Roots"];
  usableQMax = best["UsableQMax"];
  salvageOffsets = Select[
    Range[usableQMax + 1, maxOffset],
    !EffectiveZeroExprQ[resolved[[# + 1]], tol] &
  ];
  salvageExact = sectorCount === 1 &&
    EffectiveZeroExprQ[roots[[1]], tol];


  sectorCoeffs = ConstantArray[0, {sectorCount, Length[resolved]}];
  Do[
    sectorCoeffs[[i, q + 1]] = best["Coefficients"][[i, q + 1]],
    {i, sectorCount},
    {q, 0, usableQMax}
  ];

  <|
    "Sectors" -> Table[
      <|
        "ResidualB" -> roots[[i]],
        "Coefficients" -> sectorCoeffs[[i]]
      |>,
      {i, sectorCount}
    ],
    "SalvageOffsets" -> salvageOffsets,
    "SalvageExact" -> salvageExact,
    "Resolved" -> resolved
  |>
];

IntegrateAnalyticRegularizedBySubtractionLaurent[a_, b_, epsMinPower_Integer, gList_List,
    {xmin_, xmax_}, atLower_, atUpper_] := Module[
  {epsMin, epsMaxPower, baseCoeffs, result, tol, zeroQ, branchDir, branchPoint,
   branchRules, branchB, logValue, xLocal, subtractPowerQ, subtractSeries,
   remainderGList, monomialBasis, addMonomial, addMonomialWithB,
   scanSubtractedCoeff, coefficientAtLocalPower, collectSubtractedPowers,
   addSubtractedPower,
   gPower, ser, nmin, den, coeffs, localPower, coeff, logPowers, logCoeff,
   localP, resolvedCoeff, maxRelPower, relPower, epsLocal, s, basisExpr,
   seriesExpr, basisLaurent, droppedTopOrderLogTerms = 0,
   droppedUnresolvedSectorTerms = 0,
   unresolvedSectorOrdersInWindow = Infinity},

  epsMin = IntegerPower[epsMinPower];
  epsMaxPower = epsMin + Length[gList] - 1;
  tol = DiffExp`State`FEC[RationalizationTolerance];
  xLocal = DiffExp`Symbols`x;
  zeroQ[z_] := TrueQ[PossibleZeroQ[z]] ||
    TrueQ[Abs[numericAtActivePrecision[z]] < tol];

  branchDir = Which[
    NumericQ[atLower], atLower,
    NumericQ[atUpper], atUpper,
    !zeroQ[xmin], Sign[realNumericAtActivePrecision[xmin, tol]],
    !zeroQ[xmax], Sign[realNumericAtActivePrecision[xmax, tol]],
    True, 1
  ];
  branchPoint = If[zeroQ[xmin], xmax, xmin];
  branchRules = thetaRulesAtPoint[branchPoint, branchDir];
  branchB = DiffExp`Utilities`PChop[
    Expand[branchReplaceAtActivePrecision[b, branchRules]]
  ];

  subtractPowerQ[localPower_] :=
    NumericNegativeQ[a + localPower + 1, tol] ||
      NumericZeroQ[a + localPower + 1, tol];

  subtractSeries[ser_SeriesData] := Module[
    {snmin = ser[[4]], sden = ser[[6]], scoeffs},
    scoeffs = Table[
      If[subtractPowerQ[(snmin + idx - 1)/sden], 0, ser[[3, idx]]],
      {idx, Length[ser[[3]]]}
    ];
    SeriesData[ser[[1]], ser[[2]], scoeffs, ser[[4]], ser[[5]], ser[[6]]]
  ];
  subtractSeries[expr_] := If[subtractPowerQ[0], 0, expr];

  remainderGList = subtractSeries /@ gList;

  baseCoeffs = Table[
    IntegratePowerTimesSeriesAtPower[
      a, branchB, epsMin, remainderGList, targetPower,
      {xmin, xmax}, atLower, atUpper
    ],
    {targetPower, epsMin, epsMaxPower}
  ];
  result = <|"MinPower" -> epsMin, "Coefficients" -> baseCoeffs|>;

  logValue = EvaluateIntegralAtPoint[
    DiffExp`Symbols`Logx, xmax, a, atUpper
  ];

  monomialBasis[p_, logPower_Integer, maxPower_Integer, basisB_] :=
    monomialBasis[p, logPower, maxPower, basisB] = Module[
      {basisMinPower, coeffList},
      If[maxPower < 0, Return[LaurentZero[0, -1], Module]];

      If[NumericZeroQ[p, tol],
        basisMinPower = -logPower - 1;
        coeffList = Table[
          DiffExp`Utilities`PChop[
            Expand[
              Which[
                relPower === basisMinPower,
                  (-1)^logPower * Factorial[logPower] *
                    basisB^(-logPower - 1),
                relPower >= 0,
                  basisB^relPower *
                    zeroPowerSafe[logValue, relPower + logPower + 1] /
                    ((relPower + logPower + 1) * Factorial[relPower]),
                True,
                  0
              ]
            ]
          ],
          {relPower, basisMinPower, maxPower}
        ];
        Return[
          LaurentTrim[
            <|"MinPower" -> basisMinPower, "Coefficients" -> coeffList|>
          ],
          Module
        ];
      ];

      epsLocal = Unique["epsReg"];
      s = p + basisB * epsLocal;
      basisExpr = Exp[s * logValue] *
        Sum[
          Binomial[logPower, j] *
            zeroPowerSafe[logValue, logPower - j] *
            (-1)^j * Factorial[j] / s^(j + 1),
          {j, 0, logPower}
        ];
      seriesExpr = Quiet[
        Check[
          Normal[Series[basisExpr, {epsLocal, 0, maxPower}]],
          $Failed
        ]
      ];
      If[seriesExpr === $Failed,
        coeffList = Table[
          Quiet[
            Check[
              SeriesCoefficient[basisExpr, {epsLocal, 0, relPower}],
              0
            ]
          ],
          {relPower, 0, maxPower}
        ],
        coeffList = Table[
          Coefficient[seriesExpr, epsLocal, relPower],
          {relPower, 0, maxPower}
        ];
      ];
      LaurentTrim[
        <|
          "MinPower" -> 0,
          "Coefficients" ->
            DiffExp`Utilities`PChop /@ (Expand /@ coeffList)
        |>
      ]
    ];

  addMonomialWithB[basisB_, gPower_Integer, localPower_, logPower_Integer, coeff_] := Module[
    {basis},
    If[EffectiveZeroExprQ[coeff, tol], Return[Null, Module]];
    If[NumericZeroQ[basisB, tol] &&
        NumericZeroQ[a + localPower + 1, tol],
      DiffExp`Utilities`PrintWarning[
        "IntegrateSingularTermLaurent: encountered a resonant endpoint ",
        "coefficient with zero epsilon regulator; dropping contribution."
      ];
      Return[Null, Module];
    ];
    maxRelPower = epsMaxPower - gPower;
    If[maxRelPower < 0, Return[Null, Module]];
    localP = DiffExp`Utilities`PChop[
      SetPrecision[a + localPower + 1, activeNumericPrecision[]]
    ];
    basis = monomialBasis[localP, logPower, maxRelPower, basisB];
    If[Length[basis["Coefficients"]] === 0, Return[Null, Module]];
    result = LaurentAdd[
      result,
      <|
        "MinPower" -> gPower + basis["MinPower"],
        "Coefficients" ->
          DiffExp`Utilities`PChop /@
            (Expand /@ (coeff * basis["Coefficients"]))
      |>
    ];
  ];

  addMonomial[gPower_Integer, localPower_, logPower_Integer, coeff_] :=
    addMonomialWithB[branchB, gPower, localPower, logPower, coeff];

  coefficientAtLocalPower[ser_SeriesData, lp_] := Module[
    {idx = lp * ser[[6]] - ser[[4]] + 1},
    If[IntegerQ[idx] && idx >= 1 && idx <= Length[ser[[3]]],
      ser[[3, idx]],
      0
    ]
  ];
  coefficientAtLocalPower[expr_, lp_] := If[lp === 0, expr, 0];

  collectSubtractedPowers[] := DeleteDuplicates[
    Flatten[
      Table[
        ser = gList[[gIdx]];
        Which[
          MatchQ[ser, _SeriesData],
            nmin = ser[[4]];
            den = ser[[6]];
            Table[
              localPower = (nmin + idx - 1)/den;
              If[subtractPowerQ[localPower], localPower, Nothing],
              {idx, Length[ser[[3]]]}
            ],
          !EffectiveZeroExprQ[ser, tol] && subtractPowerQ[0],
            {0},
          True,
            {}
        ],
        {gIdx, Length[gList]}
      ]
    ]
  ];

  addSubtractedPower[lp_] := Module[
    {coeffList, fit},
    coeffList = Table[
      coefficientAtLocalPower[gList[[gIdx]], lp],
      {gIdx, Length[gList]}
    ];
    fit = FitResidualEndpointSectors[coeffList, branchRules];
    If[AssociationQ[fit],
      Do[
        Module[{basisB = DiffExp`Utilities`PChop[
            Expand[branchB + sector["ResidualB"]]
          ]},
          Do[
            addMonomialWithB[
              basisB,
              epsMin + gIdx - 1,
              lp,
              0,
              sector["Coefficients"][[gIdx]]
            ],
            {gIdx, Length[sector["Coefficients"]]}
          ];
        ],
        {sector, fit["Sectors"]}
      ];
      (* Offsets beyond the validated run: keep their non-log content on
         the plain branch exponent and drop residual Logx content as
         truncation-boundary noise.  Salvage is exact only for a single
         zero-root sector with no Logx content; flag everything else for
         the integrator's trust warnings. *)
      Do[
        Module[{resolvedCoeff2, logPowers2, logCoeff2, salvLogPs},
          resolvedCoeff2 = fit["Resolved"][[o + 1]];
          logPowers2 = DiffExp`SeriesOps`LogxPowerRange[resolvedCoeff2];
          Do[
            logCoeff2 = DiffExp`SeriesOps`LogxCoeffNS[resolvedCoeff2, logPower2];
            If[!EffectiveZeroExprQ[logCoeff2, tol] && logPower2 === 0,
              addMonomial[epsMin + o, lp, 0, logCoeff2]
            ],
            {logPower2, logPowers2}
          ];
          salvLogPs = Select[
            logPowers2,
            # > 0 && !EffectiveZeroExprQ[
              DiffExp`SeriesOps`LogxCoeffNS[resolvedCoeff2, #], tol
            ] &
          ];
          If[Length[salvLogPs] > 0 || !fit["SalvageExact"],
            droppedUnresolvedSectorTerms++;
            Module[{outputOrder = epsMin + o - 1},
              If[outputOrder <= epsMaxPower,
                unresolvedSectorOrdersInWindow = Min[
                  unresolvedSectorOrdersInWindow,
                  outputOrder
                ];
              ];
            ];
          ];
        ],
        {o, fit["SalvageOffsets"]}
      ];
      Return[Null, Module];
    ];

    (* The sector fit refused this tower.  The explicit-log fallback below
       integrates Logx monomials against the single extracted exponent
       branchB, which is only correct when no residual x^(r eps) sector
       structure is present.  Surface that loudly when logs remain. *)
    If[fit === $Failed &&
        AnyTrue[
          coeffList,
          Function[c, Module[{rc},
            rc = DiffExp`Utilities`PChop[
              Expand[branchReplaceAtActivePrecision[c, branchRules]]
            ];
            AnyTrue[
              DiffExp`SeriesOps`LogxPowerRange[rc],
              # > 0 && !EffectiveZeroExprQ[
                DiffExp`SeriesOps`LogxCoeffNS[rc, #], tol
              ] &
            ]
          ]]
        ],
      DiffExp`Utilities`PrintWarning[
        "IntegrateSingularTermLaurent: residual endpoint sector recovery ",
        "failed for local power ", lp, " but the subtracted coefficients ",
        "still carry Logx towers. Falling back to explicit-log integration ",
        "with the averaged exponent; the result may be WRONG if this ",
        "endpoint mixes several x^(a + b eps) sectors."
      ];
    ];

    Do[
      scanSubtractedCoeff[gIdx, epsMin + gIdx - 1, lp, coeffList[[gIdx]]],
      {gIdx, Length[coeffList]}
    ];
  ];

  scanSubtractedCoeff[gIdx_Integer, gPower_Integer, localPower_, coeff_] := Module[{},
    If[EffectiveZeroExprQ[coeff, tol] ||
        !subtractPowerQ[localPower],
      Return[Null, Module]
    ];
    resolvedCoeff = DiffExp`Utilities`PChop[
      Expand[branchReplaceAtActivePrecision[coeff, branchRules]]
    ];
    logPowers = DiffExp`SeriesOps`LogxPowerRange[resolvedCoeff];
    Do[
      logCoeff = DiffExp`SeriesOps`LogxCoeffNS[resolvedCoeff, logPower];
      If[!EffectiveZeroExprQ[logCoeff, tol],
        If[logPower > 0 && gIdx === Length[gList] && epsMaxPower > epsMin,
          droppedTopOrderLogTerms++,
          addMonomial[gPower, localPower, logPower, logCoeff]
        ]
      ],
      {logPower, logPowers}
    ];
  ];

  Do[
    addSubtractedPower[localPower],
    {localPower, collectSubtractedPowers[]}
  ];

  If[droppedTopOrderLogTerms > 0,
    DiffExp`Utilities`PrintWarning[
      "IntegrateSingularTermLaurent: dropped ",
      droppedTopOrderLogTerms,
      " explicit Logx term(s) from the highest available endpoint ",
      "subtraction coefficient. These are treated as truncation-boundary ",
      "terms and are not allowed to generate lower Laurent poles."
    ];
  ];
  If[droppedUnresolvedSectorTerms > 0,
    If[unresolvedSectorOrdersInWindow <= epsMaxPower,
      DiffExp`Utilities`PrintWarning[
        "IntegrateSingularTermLaurent: omitted ",
        droppedUnresolvedSectorTerms,
        " endpoint coefficient(s) while resumming residual Logx towers ",
        "into x^(b eps) sectors, and the omitted data affects the reported ",
        "Laurent window starting at epsilon order ",
        unresolvedSectorOrdersInWindow,
        ". Results at and above that order are NOT trustworthy; rerun with ",
        "more epsilon lookahead (e.g. a larger integration pole allowance)."
      ];
      ,
      DiffExp`Utilities`PrintWarning[
        "IntegrateSingularTermLaurent: omitted ",
        droppedUnresolvedSectorTerms,
        " highest-order endpoint coefficient(s) while resumming residual ",
        "Logx towers into x^(b eps) sectors. These terms require more ",
        "epsilon lookahead and only affect orders beyond the requested window."
      ];
    ];
  ];

  LaurentTrim[result]
];

FiniteEndpointConstant[expr_] := Module[
  {x = DiffExp`Symbols`x, logx = DiffExp`Symbols`Logx, noLog, terms},
  noLog = DiffExp`Utilities`PChop[Expand[expr /. logx -> 0]];
  terms = If[Head[noLog] === Plus, List @@ noLog, {noLog}];
  DiffExp`Utilities`PChop[Expand[Total[Select[terms, FreeQ[#, x] &]]]]
];

(* Integrate x^a * Logx^n * (series in x) *)
(* Returns a function/expression that can be evaluated at bounds *)
IntegratePowerLogMonomial[p_, n_, coeff_, tol_] := Module[
  {zeroQ, nn, pPlusOne},
  zeroQ[z_] := TrueQ[PossibleZeroQ[z]] ||
    Module[{nz = realNumericAtActivePrecision[z, tol]},
      nz =!= $Failed && NumericQ[nz] && TrueQ[Abs[nz] < tol]
    ];
  nn = IntegerPower[n];
  pPlusOne = DiffExp`Utilities`PChop[
    SetPrecision[p + 1, activeNumericPrecision[]]
  ];

  If[EffectiveZeroExprQ[coeff, tol],
    0,
    If[!IntegerQ[nn] || nn < 0, Return[0, Module]];
    If[zeroQ[pPlusOne],
      coeff * DiffExp`Symbols`Logx^(nn + 1) / (nn + 1),
      coeff * DiffExp`Symbols`x^pPlusOne / pPlusOne *
        Sum[(-1)^j * Factorial[nn] / Factorial[nn - j] *
          DiffExp`Symbols`Logx^(nn - j) / pPlusOne^j, {j, 0, nn}]
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

      Module[{maxCoeffLogPower =
          IntegerPower[DiffExp`SeriesOps`MaxLogxPower[coeff]]},
        If[IntegerQ[maxCoeffLogPower],
          Sum[
            IntegratePowerLogMonomial[
              p, n + logPow,
              DiffExp`SeriesOps`LogxCoeffNS[coeff, logPow],
              tol
            ],
            {logPow, 0, maxCoeffLogPower}
          ],
          0
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

thetaRulesAtPoint[pt_, direction_:Automatic] := Module[
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
  numericSide = realNumericAtActivePrecision[side, tol];
  If[numericSide =!= $Failed && NumericQ[numericSide] &&
      TrueQ[numericSide < -tol],
    Exp[-I Pi power],
    1
  ]
];

(* Evaluate the integral expression at a point, handling singularity limits *)
EvaluateIntegralAtPoint[expr_, pt_, a_, atSingularity_] := Module[
  {val, zeroQ, tol, branchExpr, precision, ptPrecise},
  precision = activeNumericPrecision[];
  ptPrecise = SetPrecision[pt, precision];
  tol = DiffExp`State`FEC[RationalizationTolerance];
  zeroQ = TrueQ[PossibleZeroQ[ptPrecise]] ||
    TrueQ[Abs[realNumericAtActivePrecision[ptPrecise, tol]] < tol];
  branchExpr = Chop[
    DiffExp`Utilities`PChop[
      Expand[
        branchReplaceAtActivePrecision[
          expr,
          thetaRulesAtPoint[ptPrecise, branchDirection[ptPrecise, atSingularity]]
        ]
      ]
    ],
    tol
  ];

  If[zeroQ,
    (* At x=0 singularity: terms with Logx or negative powers vanish *)
    (* Only finite constant terms survive *)
    FiniteEndpointConstant[branchExpr]
    ,
    (* Normal evaluation *)
    val = branchExpr /. DiffExp`Symbols`Logx -> Log[DiffExp`Symbols`x] /.
      DiffExp`Symbols`x -> ptPrecise;
    val = numericAtActivePrecision[val, precision];
    If[NumericQ[val] && TrueQ[Abs[numericAtActivePrecision[val, precision]] > 10^100],
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
              DiffExp`Pade`SEval2[DiffExp`Pade`GetPade[ser], ptPrecise],
              $Failed
            ]
          ];
          padeVal = If[padeVal === $Failed, $Failed,
            numericAtActivePrecision[padeVal, precision]
          ];
          If[padeVal =!= $Failed && NumericQ[padeVal] &&
              TrueQ[Abs[padeVal] < Abs[val]],
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

  If[TrueQ[realNumericAtActivePrecision[xmin] >
      realNumericAtActivePrecision[xmax]],
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

  LaurentTrim[result]
];

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

(* Sector-aware endpoint limit.  EvaluateLimitAtSingularity (and the paper
   prescription it implements) assumes DecomposeSingularity's single
   extracted exponent describes the whole series.  For multi-sector
   endpoints - several local solutions x^(a_i + b_i eps) sharing integer
   powers - that collapses distinct b_i into one averaged (or, with the
   most negative power leading, simply the wrong) exponent, and the b != 0
   drop rule then discards or keeps entire towers wholesale.  Here each
   integer power's epsilon tower is resolved into sectors with
   FitResidualEndpointSectors first, and the drop rule is applied to the
   ABSOLUTE exponent (extracted b + residual root) per sector: only the
   b = 0 sector at absolute power x^0 survives the limit. *)
EvaluateEndpointLimitSectors[seriesList_List, direction_:1] := Module[
  {tol, branchRules, decomposition, nOrders, result, coeffAtLP,
   droppedLogOffsets = {}, unresolvedTower = False},

  nOrders = Length[seriesList];
  result = Table[0, {nOrders}];
  tol = DiffExp`State`FEC[RationalizationTolerance];
  branchRules = thetaRulesAtPoint[0, direction];

  coeffAtLP[ser_SeriesData, lp_] := Module[
    {idx = lp * ser[[6]] - ser[[4]] + 1},
    If[IntegerQ[idx] && idx >= 1 && idx <= Length[ser[[3]]],
      ser[[3, idx]],
      0
    ]
  ];
  coeffAtLP[expr_, lp_] := If[lp === 0, expr, 0];

  decomposition =
    DiffExp`SingularityDecomposition`DecomposeSingularity[seriesList];

  Do[
    Module[{a = term["a"], b = term["b"], g = term["g"], branchB,
            localPowers},
      branchB = DiffExp`Utilities`PChop[
        Expand[branchReplaceAtActivePrecision[b, branchRules]]
      ];

      (* Local powers lp with a + lp <= 0: absolute x^0 contributes the
         limit; negative absolute powers are scanned for genuine (b = 0)
         divergences.  Powers above x^0 vanish at the endpoint. *)
      localPowers = DeleteDuplicates[Flatten[Table[
        Module[{ser = g[[n]]},
          Which[
            MatchQ[ser, _SeriesData],
              Module[{nmin = ser[[4]], den = ser[[6]]},
                Table[
                  Module[{lp = (nmin + idx - 1)/den},
                    If[NumericNegativeQ[a + lp, tol] ||
                        NumericZeroQ[a + lp, tol],
                      lp,
                      Nothing
                    ]
                  ],
                  {idx, Length[ser[[3]]]}
                ]
              ],
            !EffectiveZeroExprQ[ser, tol] &&
              (NumericNegativeQ[a, tol] || NumericZeroQ[a, tol]),
              {0},
            True,
              {}
          ]
        ],
        {n, Length[g]}
      ]]];

      Do[
        Module[{coeffList, fit, atZeroPower},
          atZeroPower = NumericZeroQ[a + lp, tol];
          coeffList = Table[coeffAtLP[g[[n]], lp], {n, Length[g]}];
          If[AllTrue[coeffList, EffectiveZeroExprQ[#, tol] &],
            Continue[]
          ];
          fit = FitResidualEndpointSectors[coeffList, branchRules];
          If[AssociationQ[fit],
            Do[
              Module[{absB = DiffExp`Utilities`PChop[
                  Expand[branchB + sector["ResidualB"]]
                ]},
                Which[
                  !NumericZeroQ[absB, tol],
                    (* b != 0 sector: put to zero, even at a + lp < 0 *)
                    Null,
                  atZeroPower,
                    Do[
                      result[[q]] += sector["Coefficients"][[q]],
                      {q, Min[nOrders, Length[sector["Coefficients"]]]}
                    ],
                  AnyTrue[sector["Coefficients"],
                      !EffectiveZeroExprQ[#, tol] &],
                    DiffExp`Utilities`PrintWarning[
                      "EvaluateEndpointLimitSectors: a b = 0 sector at ",
                      "negative power x^", a + lp, " does not vanish at ",
                      "the endpoint; the limit diverges. Dropping this ",
                      "contribution."
                    ]
                ];
              ],
              {sector, fit["Sectors"]}
            ];
            (* Offsets beyond the validated run: keep their non-log content
               on the plain extracted exponent (truncation-boundary
               convention), mirroring the definite-integral path. *)
            Do[
              Module[{rc = fit["Resolved"][[o + 1]], logPs, c0},
                logPs = Select[
                  DiffExp`SeriesOps`LogxPowerRange[rc],
                  # > 0 && !EffectiveZeroExprQ[
                    DiffExp`SeriesOps`LogxCoeffNS[rc, #], tol
                  ] &
                ];
                If[Length[logPs] > 0 || !fit["SalvageExact"],
                  AppendTo[droppedLogOffsets, o];
                ];
                If[atZeroPower && NumericZeroQ[branchB, tol],
                  c0 = DiffExp`SeriesOps`LogxCoeffNS[rc, 0];
                  If[!EffectiveZeroExprQ[c0, tol] && o + 1 <= nOrders,
                    result[[o + 1]] += c0;
                  ];
                ];
              ],
              {o, fit["SalvageOffsets"]}
            ];
            ,
            (* Fit refused the tower: fall back to the collapsed-exponent
               rule, loudly when Logx towers indicate unresolved sectors. *)
            Module[{resolved},
              resolved = DiffExp`Utilities`PChop[
                Expand[branchReplaceAtActivePrecision[#, branchRules]]
              ] & /@ coeffList;
              If[AnyTrue[resolved, Function[rc, AnyTrue[
                    DiffExp`SeriesOps`LogxPowerRange[rc],
                    # > 0 && !EffectiveZeroExprQ[
                      DiffExp`SeriesOps`LogxCoeffNS[rc, #], tol
                    ] &
                  ]]],
                unresolvedTower = True;
              ];
              If[atZeroPower && NumericZeroQ[branchB, tol],
                Do[
                  Module[{c0 = DiffExp`SeriesOps`LogxCoeffNS[resolved[[n]], 0]},
                    If[!EffectiveZeroExprQ[c0, tol],
                      result[[n]] += c0
                    ];
                  ],
                  {n, nOrders}
                ];
              ];
            ];
          ];
        ],
        {lp, localPowers}
      ];
    ],
    {term, decomposition}
  ];

  If[unresolvedTower,
    DiffExp`Utilities`PrintWarning[
      "EvaluateEndpointLimitSectors: endpoint tower with unresolved ",
      "x^(r eps) sector structure (sector fit failed); the limit used the ",
      "single collapsed exponent and is NOT trustworthy."
    ];
  ];
  If[Length[droppedLogOffsets] > 0,
    DiffExp`Utilities`PrintWarning[
      "EvaluateEndpointLimitSectors: dropped unresolved sector content at ",
      Length[DeleteDuplicates[droppedLogOffsets]],
      " epsilon offset(s) beyond the validated run (first affected ",
      "offset ", Min[droppedLogOffsets],
      "); limit orders at and above that offset may be incomplete."
    ];
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
  uncompressed = uncompressSeriesData[seriesData];

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
  uncompressed = uncompressSeriesData[savedData[[1, 5]]];

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
  uncompressed = uncompressSeriesData[savedData[[1, 5]]];

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
  values = numericAtActivePrecision[(xMainExpr /. xLocal -> #) & /@ localBounds];
  {Min[values], Max[values]}
];

segmentLocalCoordinateForValue[segmentData_List, value_, xMainExpr_:Automatic] := Module[
  {xLocal = DiffExp`Symbols`x, expr, localBounds, values, denom, coord},
  expr = If[xMainExpr === Automatic, segmentMainExpression[segmentData], xMainExpr];
  localBounds = segmentData[[4]];
  values = (expr /. xLocal -> #) & /@ localBounds;
  denom = values[[2]] - values[[1]];
  coord = If[TrueQ[PossibleZeroQ[denom]],
    localBounds[[1]],
    localBounds[[1]] + (value - values[[1]]) *
      (localBounds[[2]] - localBounds[[1]]) / denom
  ];
  snapLocalCoordinate[coord, localBounds]
];

snapLocalCoordinate[coord_, localBounds_List] := Module[
  {tol = DiffExp`State`FEC[RationalizationTolerance], match},
  match = SelectFirst[
    localBounds,
    TrueQ[PossibleZeroQ[coord - #]] ||
      TrueQ[NumericQ[numericAtActivePrecision[coord - #]] &&
        Abs[numericAtActivePrecision[coord - #]] < tol] &,
    Missing["NoSnap"]
  ];
  If[MissingQ[match], coord, match]
];

localSeriesData[expr_, var_, center_, order_] := Module[{ser},
  ser = Quiet[
    Series[
      SetPrecision[expr, activeNumericPrecision[]],
      {var, center, order}
    ]
  ];
  If[MatchQ[ser, _SeriesData],
    ser,
    SeriesData[var, center, {ser}, 0, 1, 1]
  ]
];

uncompressSeriesData[seriesData_] := Module[{},
  If[StringQ[seriesData],
    If[FileExistsQ[seriesData],
      Uncompress[Import[seriesData]],
      Uncompress[seriesData]
    ],
    seriesData
  ]
];

segmentWithUncompressedSeries[segmentData_List] :=
  ReplacePart[segmentData, 5 -> uncompressSeriesData[segmentData[[5]]]];

(* Helper: integrate a segment indefinitely, returning a function of x *)
IntegrateSegmentIndefinite[segmentData_List, intIndex_Integer, epsOrder_Integer] := Module[
  {currLine, lineRelation, seriesData, jacobian,
   seriesAtIndex, uncompressed, decomposition,
   localX, indefiniteLocal, indefiniteMain},

  currLine = segmentData[[1]];
  lineRelation = segmentData[[2]];
  seriesData = segmentData[[5]];

  (* Uncompress *)
  uncompressed = uncompressSeriesData[seriesData];

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
  uncompressed = uncompressSeriesData[savedData[[1, 5]]];

  numIntegrals = Length[uncompressed];
  numEpsOrders = Length[uncompressed[[1]]];

  DiffExp`Utilities`PrintInfo["Computing definite integral with prefactor from ", lower // N, " to ", upper // N][1];

  (* Find segments that overlap with [lower, upper] *)
  relevantSegments = segmentWithUncompressedSeries /@ Select[savedData,
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

  uncompressed = uncompressSeriesData[savedData[[1, 5]]];

  numIntegrals = Length[uncompressed];
  numEpsOrders = Length[uncompressed[[1]]];
  epsMaxPower = epsMinPower + numEpsOrders - 1;

  DiffExp`Utilities`PrintInfo["Computing Laurent definite integral with prefactor from ", lower // N, " to ", upper // N][1];

  relevantSegments = segmentWithUncompressedSeries /@ Select[savedData,
    Module[{segMin, segMax, actualBounds},
      actualBounds = segmentActualBounds[#];
      segMin = actualBounds[[1]];
      segMax = actualBounds[[2]];
      Not[segMax <= lower || segMin >= upper]
    ] &
  ];

  maybeDumpLaurentDefiniteIntegral[
    <|
      "Bounds" -> {lower, upper},
      "PrefactorSpec" -> prefactorSpec,
      "EpsMinPower" -> epsMinPower,
      "SavedData" -> <|
        "SegmentData" -> relevantSegments,
        "NumIntegrals" -> numIntegrals,
        "EpsilonOrder" -> epsMaxPower,
        "EpsilonMinPower" -> epsMinPower
      |>,
      "RelevantSegmentCount" -> Length[relevantSegments],
      "NumIntegrals" -> numIntegrals,
      "NumEpsOrders" -> numEpsOrders
    |>
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
  uncompressed = uncompressSeriesData[seriesData];

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
  atLowerBound = NumericZeroQ[a - lower];
  atUpperBound = NumericZeroQ[b - upper];

  (* Compute the prefactor in local coordinates.
     lineRelation gives: x_main = f(x_local)
     So x_main - lower = f(x_local) - lower, and upper - x_main = upper - f(x_local)
  *)
    Module[{prefactorLower, prefactorUpper, prefactorRational,
            smoothPrefactor, singularPowerLower, singularPowerUpper,
            totalSingularPower, modifiedDecomp,
            smoothPrefactorIsOne},

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

    If[atLowerBound && NumericZeroQ[localA],
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

    If[atUpperBound && NumericZeroQ[localB],
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
    smoothPrefactorIsOne = TrueQ[
      PossibleZeroQ[DiffExp`Utilities`PChop[Expand[smoothPrefactor - 1]]]
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
              If[smoothPrefactorIsOne,
                gSeries,
                Module[{result},
                  smoothSer = Series[
                    SetPrecision[smoothPrefactor, activeNumericPrecision[]],
                    {xLocal, gSeries[[2]], gSeries[[5]] - gSeries[[4]]}
                  ];
                  result = DiffExp`SeriesOps`SExpand[gSeries * smoothSer];
                  If[MatchQ[result, _SeriesData],
                    result,
                    localSeriesData[
                      result,
                      xLocal,
                      gSeries[[2]],
                      Min[gSeries[[5]] - gSeries[[4]], expansionOrder]
                    ]
                  ]
                ]
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
        Module[{gMinPowers, minPower},
          gMinPowers = Table[
            If[MatchQ[newG[[ord]], _SeriesData],
              newG[[ord]][[4]] / newG[[ord]][[6]],
              0
            ], {ord, Length[newG]}
          ];
          minPower = Min[gMinPowers];
          If[TrueQ[minPower < 0],
            newA = newA + minPower;
            (* Shift each g-series so it starts at x^0 *)
            newG = Table[
              If[MatchQ[newG[[ord]], _SeriesData],
                Module[{s = newG[[ord]], shiftIndex},
                  shiftIndex = IntegerPower[minPower * newG[[ord]][[6]]];
                  If[!IntegerQ[shiftIndex],
                    Return[
                      localSeriesData[
                        Normal[newG[[ord]]] * xLocal^(-minPower),
                        xLocal,
                        newG[[ord]][[2]],
                        expansionOrder
                      ],
                      Module
                    ]
                  ];
                  SeriesData[s[[1]], s[[2]], s[[3]],
                    s[[4]] - shiftIndex, s[[5]] - shiftIndex, s[[6]]]
                ],
                localSeriesData[
                  newG[[ord]] * xLocal^(-minPower),
                  xLocal,
                  0,
                  expansionOrder
                ]
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

  uncompressed = uncompressSeriesData[seriesData];

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

  atLowerBound = NumericZeroQ[a - lower];
  atUpperBound = NumericZeroQ[b - upper];

    Module[{prefactorLower, prefactorUpper, prefactorRational,
            smoothPrefactor, singularPowerLower, singularPowerUpper,
            totalSingularPower, modifiedDecomp,
            smoothPrefactorIsOne},

    prefactorRational = rationalFactor /. variable -> xMainExpr;

    singularPowerLower = 0;
    singularPowerUpper = 0;

    If[atLowerBound && NumericZeroQ[localA],
      singularPowerLower = alpha;
      prefactorLower = Abs[jacobian]^alpha *
        localSidePhase[localB - localA, alpha];
      ,
      prefactorLower = (xMainExpr - lower)^alpha;
    ];

    If[atUpperBound && NumericZeroQ[localB],
      singularPowerUpper = beta;
      prefactorUpper = Abs[jacobian]^beta *
        localSidePhase[localA - localB, beta];
      ,
      prefactorUpper = (upper - xMainExpr)^beta;
    ];

    totalSingularPower = singularPowerLower + singularPowerUpper;
    smoothPrefactor = prefactorLower * prefactorUpper * prefactorRational;
    smoothPrefactorIsOne = TrueQ[
      PossibleZeroQ[DiffExp`Utilities`PChop[Expand[smoothPrefactor - 1]]]
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
              If[smoothPrefactorIsOne,
                gSeries,
                Module[{result},
                  smoothSer = Series[
                    SetPrecision[smoothPrefactor, activeNumericPrecision[]],
                    {xLocal, gSeries[[2]], gSeries[[5]] - gSeries[[4]]}
                  ];
                  result = DiffExp`SeriesOps`SExpand[gSeries * smoothSer];
                  If[MatchQ[result, _SeriesData],
                    result,
                    localSeriesData[
                      result,
                      xLocal,
                      gSeries[[2]],
                      Min[gSeries[[5]] - gSeries[[4]], expansionOrder]
                    ]
                  ]
                ]
              ],
                  gSeries * (smoothPrefactor /. xLocal -> 0)
                ],
            0
          ],
          {ord, Length[termG]}
        ];

        Module[{gMinPowers, minPower},
          gMinPowers = Table[
            If[MatchQ[newG[[ord]], _SeriesData],
              newG[[ord]][[4]] / newG[[ord]][[6]],
              0
            ], {ord, Length[newG]}
          ];
          minPower = Min[gMinPowers];
          If[TrueQ[minPower < 0],
            newA = newA + minPower;
            newG = Table[
              If[MatchQ[newG[[ord]], _SeriesData],
                Module[{s = newG[[ord]], shiftIndex},
                  shiftIndex = IntegerPower[minPower * newG[[ord]][[6]]];
                  If[!IntegerQ[shiftIndex],
                    Return[
                      localSeriesData[
                        Normal[newG[[ord]]] * xLocal^(-minPower),
                        xLocal,
                        newG[[ord]][[2]],
                        expansionOrder
                      ],
                      Module
                    ]
                  ];
                  SeriesData[s[[1]], s[[2]], s[[3]],
                    s[[4]] - shiftIndex, s[[5]] - shiftIndex, s[[6]]]
                ],
                localSeriesData[
                  newG[[ord]] * xLocal^(-minPower),
                  xLocal,
                  0,
                  expansionOrder
                ]
              ], {ord, Length[newG]}
            ];
          ];
        ];

        <|"a" -> newA, "b" -> termB, "g" -> newG|>
      ],
      {termIdx, Length[decomposition]}
    ];

        (* Residual Logx terms after extracting x^(a+b eps) may be genuine
           resonant/local-sector data.  The subtraction formula below handles
           log monomials directly; only the top epsilon order is treated as a
           truncation boundary inside IntegrateAnalyticRegularizedBySubtractionLaurent. *)

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
