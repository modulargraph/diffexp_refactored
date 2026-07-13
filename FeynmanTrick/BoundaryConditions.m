(* ::Package:: *)
(* BoundaryConditions - Compute boundary conditions at the deepest level *)
(* Uses the generalized tadpole formula:                                  *)
(*   I_v^(n-1) = Gamma(v - L*d/2) / Gamma(v) * U^(v-(L+1)*d/2) / F^(v-L*d/2) *)
(* where U and F are Symanzik polynomials, L is the loop number.         *)
(* Multiloop-ready: follows FIRE's UF construction for Symanzik data.   *)

BeginPackage["FeynmanTrick`BoundaryConditions`", {"FeynmanTrick`"}];

ComputeSymanzikPolynomials::usage =
  "ComputeSymanzikPolynomials[propagators, loopMomenta, replacements] computes the \
Symanzik polynomials U and F using FIRE's UF construction. Returns {U, F, feynmanVars}.";

EvaluateTadpoleBoundary::usage =
  "EvaluateTadpoleBoundary[U, F, v, numLoops, epsOrder] evaluates the generalized \
tadpole I_v = Gamma(v-L*d/2)/Gamma(v) * U^(v-(L+1)*d/2) / F^(v-L*d/2) and expands \
in eps using FTConfiguration[\"DimensionExpression\"]. Returns {epsMinPower, coefficients}.";

EvaluateTadpoleNumerical::usage =
  "EvaluateTadpoleNumerical[U, F, v, numLoops, epsValue, precision] evaluates the \
generalized tadpole at a fixed numerical epsilon value. Returns a numerical result.";

DeepestLevelBoundary::usage =
  "DeepestLevelBoundary[ftData, epsOrder] computes boundary conditions for all masters \
at the deepest level of the Feynman trick iteration. Returns \
<|\"BoundaryValues\" -> {bc1, bc2, ...}, \"EpsPrefactors\" -> {k1, k2, ...}, \
\"EpsMinPower\" -> minPow|>.";

RescaledFeynmanParameters::usage =
  "RescaledFeynmanParameters[numProps] returns the rescaled Feynman parameters x_j' \
as functions of the Feynman parameters x_1,...,x_{n-1}, following eq. (2.11) of the paper \
(always combine leftmost two).";

RescaledFeynmanParametersFromSequence::usage =
  "RescaledFeynmanParametersFromSequence[numProps, combinationSequence] returns the \
rescaled Feynman parameters x_j' for a general combination sequence by tracing \
D_combined = x'_1*D_1 + ... + x'_n*D_n symbolically.";

Begin["`Private`"];

(* ============================================================ *)
(* Symanzik Polynomial Computation                              *)
(* Computes U and F directly by completing the square over      *)
(* loop momenta. Same algorithm as FIRE's UF construction.      *)
(* ============================================================ *)

ComputeSymanzikPolynomials[propagators_List, loopMomenta_List, replacements_List] :=
Module[{nProps, vs, degree, coeff, i, t2, t1, t0, k, U, F, cz},
  nProps = Length[propagators];

  (* Feynman parameter variables: x[1], ..., x[n] *)
  vs = Table[Global`x[j], {j, nProps}];

  (* Rationalize the scalar product replacement rules *)
  cz = Map[Rationalize[#, 0] &, replacements, {0, Infinity}];

  (* Combined denominator with Feynman parameters:
     degree = -Sum[D_i * x_i, {i, 1, n}]
     Note: propagators are D_j = -q_j^2 + m_j^2, so -D_j = q_j^2 - m_j^2 *)
  degree = -Sum[propagators[[i]] * vs[[i]], {i, nProps}];

  (* Complete the square for each loop momentum *)
  coeff = 1;
  Do[
    k = loopMomenta[[i]];
    t2 = Coefficient[degree, k, 2];
    t1 = Coefficient[degree, k, 1];
    t0 = Coefficient[degree, k, 0];

    (* Apply scalar product rules *)
    t2 = t2 //. cz;
    t1 = t1 //. cz;
    t0 = t0 //. cz;

    If[t2 === 0,
      Print["Error: Coefficient of ", k, "^2 is zero. Check propagator definitions."];
      Return[$Failed];
    ];

    coeff = coeff * t2;
    degree = Together[t0 - t1^2 / (4 * t2)];
  , {i, Length[loopMomenta]}];

  (* U = product of quadratic coefficients (first Symanzik polynomial) *)
  U = Together[coeff] //. cz;

  (* F = -coeff * degree (second Symanzik polynomial, includes masses) *)
  F = ExpandAll[Together[-coeff * degree] //. cz] //. cz;

  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 2,
    Print["  Symanzik U = ", U];
    Print["  Symanzik F = ", F];
    Print["  Feynman parameters: ", vs];
  ];

  {U, F, vs}
];


(* ============================================================ *)
(* Rescaled Feynman Parameters                                   *)
(* Following eq. (2.11) of the paper:                           *)
(*   x_1' = prod_{i=1}^{n-1} x_i                              *)
(*   x_j' = (1-x_{j-1}) * prod_{i=j}^{n-1} x_i  for j=2,...,n-1 *)
(*   x_n' = (1-x_{n-1})                                        *)
(* ============================================================ *)

RescaledFeynmanParameters[numProps_Integer] :=
Module[{n, params, rescaled},
  n = numProps;  (* Number of original propagators *)

  (* Feynman parameters: x[1], ..., x[n-1] *)
  params = Table[Global`xFT[j], {j, n - 1}];

  rescaled = Table[
    Which[
      j == 1,
        Product[params[[i]], {i, 1, n - 1}],
      j <= n - 1,
        (1 - params[[j - 1]]) * Product[params[[i]], {i, j, n - 1}],
      j == n,
        (1 - params[[n - 1]])
    ],
    {j, n}
  ];

  {rescaled, params}
];


(* General rescaling from arbitrary combination sequence.
   Traces D_combined = x'_1*D_1 + ... + x'_n*D_n by symbolically
   applying each combination step. Works for any sequence. *)
RescaledFeynmanParametersFromSequence[numProps_Integer, combinationSequence_List] :=
Module[{nLevels, params, coeffs, k, ci, cj, finalPos, rescaled},
  nLevels = Length[combinationSequence];
  params = Table[Global`xFT[k], {k, nLevels}];

  (* coeffs[[p]] = list of coefficients of original props in the
     propagator currently at position p. Start as identity. *)
  coeffs = IdentityMatrix[numProps];

  Do[
    {ci, cj} = combinationSequence[[k]];
    (* Position ci gets: param*old_ci + (1-param)*old_cj *)
    coeffs[[ci]] = params[[k]] * coeffs[[ci]] + (1 - params[[k]]) * coeffs[[cj]];
  , {k, nLevels}];

  (* The fully combined propagator is at the last combination's first position *)
  finalPos = combinationSequence[[-1, 1]];
  rescaled = coeffs[[finalPos]];  (* length = numProps *)

  {Expand /@ rescaled, params}
];


(* Evaluate rescaled params at numerical Feynman parameter values *)
EvaluateRescaledParams[numProps_Integer, paramValues_List] :=
Module[{rescaled, params, rules},
  {rescaled, params} = RescaledFeynmanParameters[numProps];

  (* Substitute numerical values *)
  rules = Thread[params -> paramValues];
  rescaled /. rules
];


(* ============================================================ *)
(* Generalized Tadpole Evaluation                               *)
(* ============================================================ *)

finitePositiveRealNumericalQ[value_] := Module[{numeric},
  If[
    !NumericQ[value] ||
      !FreeQ[value, Indeterminate | ComplexInfinity | _DirectedInfinity],
    Return[False]
  ];
  numeric = Quiet[Check[N[value], $Failed]];
  numeric =!= $Failed && TrueQ[Im[numeric] == 0] && TrueQ[Re[numeric] > 0]
];


exactAffineDimensionData[dimension_, eps_] := Module[
  {polynomial, constant, slope},
  polynomial = Quiet[Check[Expand[Cancel[Together[dimension]]], $Failed]];
  If[
    polynomial === $Failed || !FreeQ[polynomial, _Real] ||
      !PolynomialQ[polynomial, eps] || !TrueQ[Exponent[polynomial, eps] <= 1],
    Return[$Failed]
  ];
  constant = Together[Coefficient[polynomial, eps, 0]];
  slope = Together[Coefficient[polynomial, eps, 1]];
  If[
    !FreeQ[{constant, slope}, eps] ||
      !And @@ (NumericQ /@ {constant, slope}) ||
      !FreeQ[{constant, slope}, _Real | _DirectedInfinity] ||
      !And @@ (TrueQ[Im[N[#, 20]] == 0] & /@ {constant, slope}),
    Return[$Failed]
  ];
  {constant, slope}
];


multiplyLinearCoefficientVector[
    coefficients_List, constant_, slope_, maxDegree_Integer] := Module[
  {output, inputDegree, outputDegree, degree},
  inputDegree = Length[coefficients] - 1;
  outputDegree = Min[maxDegree, inputDegree + 1];
  output = ConstantArray[0, outputDegree + 1];
  Do[
    output[[degree + 1]] += constant * coefficients[[degree + 1]];
    If[degree + 1 <= outputDegree,
      output[[degree + 2]] += slope * coefficients[[degree + 1]]
    ],
    {degree, 0, Min[inputDegree, maxDegree]}
  ];
  output
];


truncatedCoefficientConvolution[left_List, right_List, maxDegree_Integer] := Module[
  {leftPadded, rightPadded},
  leftPadded = PadRight[Take[left, UpTo[maxDegree + 1]], maxDegree + 1, 0];
  rightPadded = PadRight[Take[right, UpTo[maxDegree + 1]], maxDegree + 1, 0];
  Table[
    leftPadded[[1 ;; degree + 1]] .
      Reverse[rightPadded[[1 ;; degree + 1]]],
    {degree, 0, maxDegree}
  ]
];


(* Fast numerical Laurent expansion for the branch-insensitive Euclidean case.
   Write d=d0+d1 eps and

     Gamma[a+b eps] U^(c+u eps) / F^(a+b eps)

   as a finite integer Gamma shift times Gamma[1+b eps] Exp[lambda eps].
   The latter exponential is generated from its log coefficients in O(N^2),
   avoiding Mathematica's much more expensive symbolic Gamma Series. *)
fastNumericalTadpoleBoundary[
    Uval_, Fval_, v_Integer, numLoops_Integer, epsOrder_Integer,
    eps_, dimension_, precision_Integer] := Module[
  {dimensionData, d0, d1, a, b, c, uSlope, guardDigits, workPrecision,
   uWorking, fWorking, lambda, minPow, maxAnalyticDegree,
   logCoefficients, exponentialCoefficients, shiftCoefficients,
   denominatorCoefficients, reciprocalCoefficients, analyticCoefficients,
   constant, m, degree, k},

  If[
    v <= 0 || numLoops < 0 || epsOrder < 0 ||
      !finitePositiveRealNumericalQ[Uval] ||
      !finitePositiveRealNumericalQ[Fval],
    Return[$Failed]
  ];
  dimensionData = exactAffineDimensionData[dimension, eps];
  If[dimensionData === $Failed, Return[$Failed]];
  {d0, d1} = dimensionData;

  a = Together[v - numLoops*d0/2];
  b = Together[-numLoops*d1/2];
  c = Together[v - (numLoops + 1)*d0/2];
  uSlope = Together[-(numLoops + 1)*d1/2];
  If[
    !IntegerQ[a] || !And @@ (NumericQ /@ {b, c, uSlope}) ||
      !And @@ (TrueQ[Im[N[#, 20]] == 0] & /@ {b, c, uSlope}) ||
      (a <= 0 && TrueQ[PossibleZeroQ[b]]),
    Return[$Failed]
  ];

  guardDigits = Max[20, 10 + Ceiling[Log[10, Max[2, epsOrder + 2]]]];
  workPrecision = precision + guardDigits;
  uWorking = N[SetPrecision[Uval, workPrecision], workPrecision];
  fWorking = N[SetPrecision[Fval, workPrecision], workPrecision];
  lambda = N[uSlope*Log[uWorking] - b*Log[fWorking], workPrecision];

  minPow = If[a <= 0, -1, 0];
  maxAnalyticDegree = epsOrder - minPow;
  logCoefficients = N[
    Table[
      If[degree == 1,
        lambda - b*EulerGamma,
        (-1)^degree*Zeta[degree]*b^degree/degree
      ],
      {degree, 1, maxAnalyticDegree}
    ],
    workPrecision
  ];
  exponentialCoefficients = ConstantArray[0, maxAnalyticDegree + 1];
  exponentialCoefficients[[1]] = SetPrecision[1, workPrecision];
  Do[
    exponentialCoefficients[[degree + 1]] =
      ((Range[degree]*logCoefficients[[1 ;; degree]]) .
        Reverse[exponentialCoefficients[[1 ;; degree]]])/degree,
    {degree, 1, maxAnalyticDegree}
  ];

  If[a >= 1,
    shiftCoefficients = {1};
    Do[
      shiftCoefficients = multiplyLinearCoefficientVector[
        shiftCoefficients, k, b, maxAnalyticDegree],
      {k, 1, a - 1}
    ];
    shiftCoefficients = N[shiftCoefficients, workPrecision],

    (* For a=-m, factor the simple pole explicitly:
       Gamma[-m+b eps] = eps^-1 Gamma[1+b eps] /
         (b Product[b eps-j,{j,1,m}]). *)
    m = -a;
    denominatorCoefficients = {b};
    Do[
      denominatorCoefficients = multiplyLinearCoefficientVector[
        denominatorCoefficients, -k, b, maxAnalyticDegree],
      {k, 1, m}
    ];
    denominatorCoefficients = N[denominatorCoefficients, workPrecision];
    reciprocalCoefficients = ConstantArray[0, maxAnalyticDegree + 1];
    reciprocalCoefficients[[1]] = 1/denominatorCoefficients[[1]];
    Do[
      reciprocalCoefficients[[degree + 1]] =
        -Total[Table[
          denominatorCoefficients[[k + 1]] *
            reciprocalCoefficients[[degree - k + 1]],
          {k, 1, Min[degree, Length[denominatorCoefficients] - 1]}
        ]]/denominatorCoefficients[[1]],
      {degree, 1, maxAnalyticDegree}
    ];
    shiftCoefficients = reciprocalCoefficients
  ];

  analyticCoefficients = truncatedCoefficientConvolution[
    exponentialCoefficients, shiftCoefficients, maxAnalyticDegree];
  constant = N[
    uWorking^c/(Gamma[v]*fWorking^a),
    workPrecision
  ];
  {minPow, SetPrecision[constant*analyticCoefficients, precision]}
];


seriesTadpoleBoundary[
    Uval_, Fval_, v_Integer, numLoops_Integer, epsOrder_Integer,
    eps_, dExpr_, precision_Integer] :=
Module[{gammaArg, UPow, FPow, fullExpr, series, coeffs, minPow},

  gammaArg = v - numLoops * dExpr / 2;
  UPow = v - (numLoops + 1) * dExpr / 2;
  FPow = v - numLoops * dExpr / 2;

  (* Keep this expression and Series path as the exact fallback for complex,
     nonpositive, non-affine, or otherwise branch-sensitive inputs. *)
  Module[{uTerm, fTerm},
    uTerm = If[Chop[Uval - 1, 10^(-precision/2)] === 0,
      1,
      SetPrecision[Uval, precision]^UPow
    ];
    fTerm = If[Chop[Fval - 1, 10^(-precision/2)] === 0,
      1,
      SetPrecision[Fval, precision]^FPow
    ];
    fullExpr = Gamma[gammaArg] / Gamma[v] * uTerm / fTerm;
  ];

  series = Series[fullExpr, {eps, 0, epsOrder}];
  If[Head[series] === SeriesData,
    minPow = series[[4]];
    coeffs = PadRight[series[[3]], epsOrder - minPow + 1, 0],
    minPow = 0;
    coeffs = {series}
  ];
  {minPow, coeffs}
];


EvaluateTadpoleBoundary[Uval_?NumericQ, Fval_?NumericQ, v_Integer, numLoops_Integer, epsOrder_Integer] :=
Module[{eps, coeffs, minPow, precision, dExpr, result},

  precision = FeynmanTrick`Private`$FTConfig["WorkingPrecision"];
  If[!IntegerQ[precision] || precision < 50, precision = 200];

  eps = FeynmanTrick`Private`$FTConfig["EpsilonSymbol"];
  dExpr = FeynmanTrick`Private`DimensionExpression[];

  result = fastNumericalTadpoleBoundary[
    Uval, Fval, v, numLoops, epsOrder, eps, dExpr, precision];
  If[result === $Failed,
    result = seriesTadpoleBoundary[
      Uval, Fval, v, numLoops, epsOrder, eps, dExpr, precision]
  ];
  {minPow, coeffs} = result;

  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 2,
    Print["  Tadpole boundary: min eps power = ", minPow];
    Print["  Leading coefficient = ", coeffs[[1]]];
  ];

  result
];


(* Numerical evaluation at fixed eps (for eps-sampling approach) *)
EvaluateTadpoleNumerical[Uval_?NumericQ, Fval_?NumericQ, v_Integer, numLoops_Integer, epsValue_?NumericQ, precision_Integer:200] :=
Module[{d, gammaArg, UPow, FPow, result},
  d = FeynmanTrick`Private`DimensionExpression[] /.
    FeynmanTrick`Private`$FTConfig["EpsilonSymbol"] -> epsValue;
  gammaArg = v - numLoops * d / 2;
  UPow = v - (numLoops + 1) * d / 2;
  FPow = v - numLoops * d / 2;

  result = N[
    Gamma[gammaArg] / Gamma[v] *
    SetPrecision[Uval, precision]^UPow /
    SetPrecision[Fval, precision]^FPow,
    precision
  ];

  result
];


(* ============================================================ *)
(* DeepestLevelBoundary - Main Entry Point                      *)
(* ============================================================ *)

DeepestLevelBoundary[ftData_Association, epsOrder_Integer:4] :=
Module[{deepestLevel, levelData, originalTopology, propagators, loopMomenta,
        externalMomenta, replacements, numLoops, masters,
        ufResult, U, F, feynVars, numOriginalProps, fixedValue, kinPoint,
        Uval, Fval, bcValues, epsPrefactors, epsMinPower,
        combinationSequence, nLevels, precision},

  nLevels = ftData["NumLevels"];
  deepestLevel = nLevels;
  precision = FeynmanTrick`Private`$FTConfig["WorkingPrecision"];
  If[!IntegerQ[precision] || precision < 50, precision = 200];

  (* Check that the deepest level is computed *)
  If[!KeyExistsQ[ftData["Levels"], deepestLevel],
    Print["Error: Level ", deepestLevel, " not built yet."];
    Return[$Failed];
  ];

  levelData = ftData["Levels"][deepestLevel];
  masters = levelData["Masters"];

  (* Validate masters *)
  If[!ListQ[masters] || Length[masters] == 0,
    Print["Error: No valid masters found at deepest level ", deepestLevel];
    Print["  masters = ", masters];
    Return[$Failed];
  ];
  If[!AllTrue[masters, ListQ],
    Print["Error: Masters should be a list of index lists, got: ", masters];
    Return[$Failed];
  ];

  (* Use the ORIGINAL topology (level 0) for Symanzik polynomial computation.
     Per eq. 2.16 of the paper: U_tilde and F_tilde are obtained from the
     original Symanzik polynomials by substituting rescaled Feynman parameters. *)
  originalTopology = ftData["TopTopology"];
  propagators = originalTopology["Propagators"];
  loopMomenta = originalTopology["LoopMomenta"];
  externalMomenta = originalTopology["ExternalMomenta"];
  replacements = originalTopology["Replacements"];
  numLoops = Length[loopMomenta];
  numOriginalProps = Length[propagators];

  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["Computing deepest level boundary at level ", deepestLevel];
    Print["  Number of masters: ", Length[masters]];
    Print["  Number of loops: ", numLoops];
    Print["  Original propagators: ", numOriginalProps];
  ];

  (* Get fixed parameter value and kinematic point *)
  fixedValue = ftData["FixedParamValue"];
  kinPoint = ftData["NumericalPoint"];
  combinationSequence = ftData["CombinationSequence"];

  (* Compute Symanzik polynomials for the ORIGINAL topology *)
  ufResult = ComputeSymanzikPolynomials[propagators, loopMomenta, replacements];
  If[ufResult === $Failed, Return[$Failed]];
  {U, F, feynVars} = ufResult;

  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 2,
    Print["  Symanzik U (symbolic) = ", U];
    Print["  Symanzik F (symbolic) = ", F];
    Print["  Feynman variables from UF: ", feynVars];
  ];

  (* Compute rescaled Feynman parameters from the combination sequence.
     This traces D_combined = x'_1*D_1 + ... + x'_n*D_n symbolically. *)
  Module[{rescaled, params, rescaledNumerical, paramRules, subRules},
    {rescaled, params} = RescaledFeynmanParametersFromSequence[
      numOriginalProps, combinationSequence
    ];

    If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 2,
      Print["  Rescaled Feynman parameters (symbolic): ", rescaled];
    ];

    (* Evaluate at fixed parameter value (all xx_k = fixedValue) *)
    paramRules = Thread[params -> Table[SetPrecision[fixedValue, precision], {nLevels}]];
    rescaledNumerical = rescaled /. paramRules;

    If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 2,
      Print["  Rescaled Feynman parameters (numerical): ", rescaledNumerical];
      Print["  Sum of rescaled params: ", Total[rescaledNumerical]];
    ];

    (* Substitute Feynman parameters once, then freeze the exact kinematic
       point to a fixed point.  This keeps boundary kinematics identical to
       the FIRE and level-topology paths for chained assignments. *)
    subRules = Thread[feynVars -> rescaledNumerical];
    Uval = FeynmanTrick`FIREInterface`Private`applyFIRENumericalPoint[
      Expand[U /. subRules], kinPoint];
    Fval = FeynmanTrick`FIREInterface`Private`applyFIRENumericalPoint[
      Expand[F /. subRules], kinPoint];
    If[Uval =!= $Failed, Uval = N[Uval, precision]];
    If[Fval =!= $Failed, Fval = N[Fval, precision]];
  ];

  If[Uval === $Failed || Fval === $Failed,
    Print["Error: NumericalPoint substitutions failed in the deepest boundary."];
    Return[$Failed]];

  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["  U_tilde(numerical) = ", Uval];
    Print["  F_tilde(numerical) = ", Fval];
  ];

  (* For each master integral, compute the boundary condition *)
  (* The master at the deepest level is typically I_{v} where v is determined by *)
  (* the number of propagators with non-zero indices *)
  bcValues = Table[
    Module[{masterIndices, v, result, minPow, coeffs},
      masterIndices = masters[[masterIdx]];
      v = Total[Select[masterIndices, # > 0 &]];

      If[v == 0,
        (* Zero integral *)
        {0, Table[0, {epsOrder + 1}]},
        (* Evaluate tadpole formula *)
        {minPow, coeffs} = EvaluateTadpoleBoundary[Uval, Fval, v, numLoops, epsOrder];
        {minPow, coeffs}
      ]
    ],
    {masterIdx, Length[masters]}
  ];

  (* Determine a global eps prefactor to make all boundaries finite.
     A per-master boundary prefactor would require transforming the differential
     equation matrix. A single global factor leaves the matrix unchanged. *)
  epsMinPower = Min[bcValues[[All, 1]]];
  epsPrefactors = Table[Max[0, -epsMinPower], {Length[masters]}];

  (* Apply prefactors to get finite boundary values *)
  (* After global prefactor eps^k, the boundary for J_i = eps^k * I_i starts at eps^0 *)
  Module[{shiftedBCs, shiftedOrder},
    shiftedOrder = epsOrder + Max[epsPrefactors];
    shiftedBCs = Table[
      Module[{minPow, coeffs, shift, shiftedCoeffs},
        {minPow, coeffs} = bcValues[[masterIdx]];
        shift = epsPrefactors[[masterIdx]];

        (* Coefficients are indexed by original powers minPow, minPow+1, ...
           The shifted eps^n coefficient uses original eps^(n-shift). *)
        shiftedCoeffs = Table[
          Module[{origPower = n - shift, idx},
            idx = origPower - minPow + 1;
            If[IntegerQ[idx] && idx >= 1 && idx <= Length[coeffs],
              coeffs[[idx]],
              0
            ]
          ],
          {n, 0, shiftedOrder}
        ];
        shiftedCoeffs
      ],
      {masterIdx, Length[masters]}
    ];

    If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
      Print["  Eps prefactors: ", epsPrefactors];
      Print["  Boundary values (finite): ", shiftedBCs[[All, 1]]];
    ];

    <|
      "BoundaryValues" -> shiftedBCs,
      "EpsPrefactors" -> epsPrefactors,
      "EpsMinPower" -> epsMinPower,
      "WorkingEpsilonOrder" -> shiftedOrder,
      "RequestedEpsilonOrder" -> epsOrder,
      "Masters" -> masters,
      "NumLoops" -> numLoops,
      "Uval" -> Uval,
      "Fval" -> Fval
    |>
  ]
];


End[];
EndPackage[];
