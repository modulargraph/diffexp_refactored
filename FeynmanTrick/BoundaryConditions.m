(* ::Package:: *)
(* BoundaryConditions - Compute boundary conditions at the deepest level *)
(* Uses the generalized tadpole formula:                                  *)
(*   I_v^(n-1) = Gamma(v - L*d/2) / Gamma(v) * U^(v-(L+1)*d/2) / F^(v-L*d/2) *)
(* where U and F are Symanzik polynomials, L is the loop number.         *)
(* Multiloop-ready: uses FIRE6's UF function for Symanzik computation.   *)

BeginPackage["FeynmanTrick`BoundaryConditions`", {"FeynmanTrick`"}];

ComputeSymanzikPolynomials::usage =
  "ComputeSymanzikPolynomials[propagators, loopMomenta, replacements] computes the \
Symanzik polynomials U and F using FIRE6's UF function. Returns {U, F, feynmanVars}.";

EvaluateTadpoleBoundary::usage =
  "EvaluateTadpoleBoundary[U, F, v, numLoops, epsOrder] evaluates the generalized \
tadpole I_v = Gamma(v-L*d/2)/Gamma(v) * U^(v-(L+1)*d/2) / F^(v-L*d/2) and expands \
in eps (d=4-2*eps) to the given order. Returns {epsMinPower, coefficients}.";

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
(* loop momenta. Same algorithm as FIRE6's UF function.         *)
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

(*
  I_v^(n-1) = Gamma(v - L*d/2) / Gamma(v) * U^(v-(L+1)*d/2) / F^(v-L*d/2)

  With d = 4 - 2*eps:
  - v - L*d/2 = v - L*(2-eps) = v - 2L + L*eps
  - v - (L+1)*d/2 = v - (L+1)*(2-eps) = v - 2(L+1) + (L+1)*eps
  - v - L*d/2 = v - 2L + L*eps (same as first)
*)

EvaluateTadpoleBoundary[Uval_?NumericQ, Fval_?NumericQ, v_Integer, numLoops_Integer, epsOrder_Integer] :=
Module[{eps, gammaArg, gammaPrefactor, UPow, FPow, fullExpr, series, coeffs, minPow,
        precision},

  precision = FeynmanTrick`Private`$FTConfig["WorkingPrecision"];
  If[!IntegerQ[precision] || precision < 50, precision = 200];

  eps = FeynmanTrick`FTeps;

  (* Arguments *)
  gammaArg = v - numLoops * (2 - eps);  (* v - L*d/2 *)
  UPow = v - (numLoops + 1) * (2 - eps);  (* v - (L+1)*d/2 *)
  FPow = v - numLoops * (2 - eps);  (* v - L*d/2 *)

  (* Build the expression symbolically.
     Handle U=1 and F=1 cases explicitly: 1.0^(symbolic) doesn't simplify in Mathematica,
     which prevents Series from expanding properly. *)
  Module[{uTerm, fTerm},
    uTerm = If[Chop[Uval - 1, 10^(-precision/2)] === 0,
      1,  (* U=1 for 1-loop; avoid 1.^(symbolic) which doesn't simplify *)
      SetPrecision[Uval, precision]^UPow
    ];
    fTerm = If[Chop[Fval - 1, 10^(-precision/2)] === 0,
      1,  (* F=1; same issue as U=1 *)
      SetPrecision[Fval, precision]^FPow
    ];
    fullExpr = Gamma[gammaArg] / Gamma[v] * uTerm / fTerm;
  ];

  (* Series expand in eps around 0 *)
  series = Series[fullExpr, {eps, 0, epsOrder}];

  (* Extract the minimum power and coefficients *)
  If[Head[series] === SeriesData,
    minPow = series[[4]];  (* Minimum power of eps *)
    coeffs = series[[3]];  (* Coefficients *)
    (* Pad with zeros if needed *)
    coeffs = PadRight[coeffs, epsOrder - minPow + 1, 0];
    ,
    minPow = 0;
    coeffs = {series};
  ];

  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 2,
    Print["  Tadpole boundary: min eps power = ", minPow];
    Print["  Leading coefficient = ", coeffs[[1]]];
  ];

  {minPow, coeffs}
];


(* Numerical evaluation at fixed eps (for eps-sampling approach) *)
EvaluateTadpoleNumerical[Uval_?NumericQ, Fval_?NumericQ, v_Integer, numLoops_Integer, epsValue_?NumericQ, precision_Integer:200] :=
Module[{d, gammaArg, UPow, FPow, result},
  d = 4 - 2*epsValue;
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
        combinationSequence, nLevels},

  nLevels = ftData["NumLevels"];
  deepestLevel = nLevels;

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
    paramRules = Thread[params -> Table[fixedValue, {nLevels}]];
    rescaledNumerical = rescaled /. paramRules;

    If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 2,
      Print["  Rescaled Feynman parameters (numerical): ", rescaledNumerical];
      Print["  Sum of rescaled params: ", Total[rescaledNumerical]];
    ];

    (* Substitute rescaled params into U and F, plus kinematic values *)
    subRules = Join[
      Thread[feynVars -> rescaledNumerical],
      kinPoint
    ];

    Uval = U /. subRules // Expand // N;
    Fval = F /. subRules // Expand // N;
  ];

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
