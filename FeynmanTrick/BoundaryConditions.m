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
  "RescaledFeynmanParameters[combinationSequence, numOriginalProps] returns the \
rescaled Feynman parameters x_j' as functions of the Feynman parameters x_1,...,x_{n-1}, \
following eq. (2.11) of the paper.";

Begin["`Private`"];

(* ============================================================ *)
(* Symanzik Polynomial Computation via FIRE6's UF               *)
(* ============================================================ *)

(* Load FIRE6.m if not already loaded *)
loadFIRE6[] := Module[{firePath},
  If[!ValueQ[Global`UF] || Head[Global`UF] =!= Symbol,
    firePath = FileNameJoin[{
      ParentDirectory[DirectoryName[$InputFileName]],
      "Dependencies", "fire", "FIRE6", "FIRE6.m"
    }];
    If[FileExistsQ[firePath],
      Block[{$ContextPath},
        Quiet[Get[firePath], {General::shdw, Symbol::shdw}];
      ];
      If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 2,
        Print["Loaded FIRE6.m for Symanzik polynomial computation."];
      ];
    ,
      (* Try alternative path *)
      firePath = FeynmanTrick`Private`$FTConfig["FIREPath"];
      If[StringQ[firePath],
        firePath = FileNameJoin[{DirectoryName[firePath], "..", "FIRE6.m"}];
        If[FileExistsQ[firePath],
          Block[{$ContextPath}, Quiet[Get[firePath], {General::shdw}]];
        ];
      ];
    ];
  ];
];


ComputeSymanzikPolynomials[propagators_List, loopMomenta_List, replacements_List] :=
Module[{result, U, F, vars},
  loadFIRE6[];

  (* Call FIRE6's UF function:
     UF[loopMomenta, propagators, scalarProductRules]
     Returns {U, F, feynmanVars}
     where U is the first Symanzik polynomial and F is the second *)
  result = UF[loopMomenta, propagators, replacements];

  If[result === {0, 0, {}},
    Print["Error: UF returned degenerate result. Check propagator definitions."];
    Return[$Failed];
  ];

  U = result[[1]];
  F = result[[2]];
  vars = result[[3]];

  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 2,
    Print["  Symanzik U = ", U];
    Print["  Symanzik F = ", F];
    Print["  Feynman parameters: ", vars];
  ];

  {U, F, vars}
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

  (* Build the expression symbolically *)
  fullExpr = Gamma[gammaArg] / Gamma[v] *
             SetPrecision[Uval, precision]^UPow /
             SetPrecision[Fval, precision]^FPow;

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
Module[{deepestLevel, levelData, topology, propagators, loopMomenta,
        externalMomenta, replacements, numLoops, masters,
        UF, U, F, feynVars, rescaledParams, paramValues,
        numOriginalProps, fixedValue, kinPoint,
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
  topology = levelData["Topology"];
  propagators = topology["Propagators"];
  loopMomenta = topology["LoopMomenta"];
  externalMomenta = topology["ExternalMomenta"];
  replacements = topology["Replacements"];
  numLoops = Length[loopMomenta];
  masters = levelData["Masters"];

  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["Computing deepest level boundary at level ", deepestLevel];
    Print["  Number of masters: ", Length[masters]];
    Print["  Number of loops: ", numLoops];
  ];

  (* Get fixed parameter value and kinematic point *)
  fixedValue = ftData["FixedParamValue"];
  kinPoint = ftData["NumericalPoint"];

  (* Compute Symanzik polynomials for the deepest-level topology *)
  UF = ComputeSymanzikPolynomials[propagators, loopMomenta, replacements];
  If[UF === $Failed, Return[$Failed]];
  {U, F, feynVars} = UF;

  (* Substitute the Feynman parameter (xx) at its fixed value *)
  (* The deepest level's propagators contain FeynmanTrick`xx as a variable *)
  (* Also substitute any remaining kinematic values *)
  Module[{allRules, xxSym},
    xxSym = FeynmanTrick`xx;
    allRules = Join[
      {xxSym -> fixedValue},
      kinPoint,
      replacements
    ];

    Uval = U //. allRules // N;
    Fval = F //. allRules // N;
  ];

  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 2,
    Print["  U(numerical) = ", Uval];
    Print["  F(numerical) = ", Fval];
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

  (* Determine eps prefactors to make all boundaries finite *)
  (* If minPow < 0 for master i, we need prefactor eps^(-minPow) *)
  epsMinPower = Min[bcValues[[All, 1]]];
  epsPrefactors = Table[-bcValues[[masterIdx, 1]], {masterIdx, Length[masters]}];

  (* Apply prefactors to get finite boundary values *)
  (* After prefactor eps^k_i, the boundary for J_i = eps^k_i * I_i starts at eps^0 *)
  Module[{shiftedBCs},
    shiftedBCs = Table[
      Module[{minPow, coeffs, shift, shiftedCoeffs},
        {minPow, coeffs} = bcValues[[masterIdx]];
        shift = epsPrefactors[[masterIdx]];  (* = -minPow *)

        (* The boundary of J_i = eps^shift * I_i is:
           coeffs shifted by 'shift' positions *)
        If[shift > 0,
          (* Drop the first 'shift' coefficients (they were at negative eps powers) *)
          shiftedCoeffs = Drop[coeffs, shift];
          (* Pad to epsOrder+1 length *)
          shiftedCoeffs = PadRight[shiftedCoeffs, epsOrder + 1, 0];
          ,
          shiftedCoeffs = coeffs;
          shiftedCoeffs = Take[shiftedCoeffs, Min[Length[shiftedCoeffs], epsOrder + 1]];
          shiftedCoeffs = PadRight[shiftedCoeffs, epsOrder + 1, 0];
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
      "Masters" -> masters,
      "NumLoops" -> numLoops,
      "Uval" -> Uval,
      "Fval" -> Fval
    |>
  ]
];


End[];
EndPackage[];
