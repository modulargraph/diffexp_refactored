(* ::Package:: *)
(* PropagatorAlgebra - Propagator manipulation utilities *)
(* Handles derivative decomposition and shifted integral construction *)

BeginPackage["FeynmanTrick`PropagatorAlgebra`", {"FeynmanTrick`"}];

DecomposePropagatorDerivative::usage =
  "DecomposePropagatorDerivative[propagators, loopMomenta, variable, replacements] \
returns {coeffMatrix, constVector} such that \
d(D_k)/d(variable) = Sum_j coeffMatrix[[k,j]] * D_j + constVector[[k]]. \
The coefficients may depend on variable and external parameters.";

DifferentiatedIntegrals::usage =
  "DifferentiatedIntegrals[indexVector, coeffMatrix, constVector] returns a list of \
{shiftedIndex, coefficient} pairs appearing in d/dx of an integral with the given indices.";

CombinePropagatorSymbolic::usage =
  "CombinePropagatorSymbolic[Di, Dj, param] returns param*Di + (1-param)*Dj expanded.";

FeynmanTrickDecomposition::usage =
  "FeynmanTrickDecomposition[nProps, combinedPos, otherPos, param] returns \
{coeffMatrix, constVector} for the special case of a Feynman trick combined \
propagator. This is the fast path: only the combined propagator depends on param, \
and its derivative decomposition is {1/param, -1/param} at positions {combinedPos, otherPos}.";

Begin["`Private`"];

(* ============================================================ *)
(* DecomposePropagatorDerivative                                 *)
(* Given propagators D_1,...,D_n and a variable, computes:       *)
(*   dD_k/dvar = Sum_j c_{kj} D_j + c_{k0}                     *)
(*                                                               *)
(* Method: Use CoefficientRules to extract the polynomial        *)
(* structure in loop momenta, then solve a linear system.        *)
(* ============================================================ *)

DecomposePropagatorDerivative[propagators_List, loopMomenta_List, variable_,
                              replacements_List:{}] :=
Module[{nProps, derivs, expandedProps, loopVectors, derivVectors,
        propMatrix, coeffMatrix, constVector, k, sol, propConsts, derivConsts,
        monomialKeys, allKeys, nMon},

  nProps = Length[propagators];

  (* Step 1: Expand propagators and compute derivatives *)
  expandedProps = Expand[# /. replacements] & /@ propagators;
  derivs = Expand[D[#, variable] /. replacements] & /@ propagators;

  (* Step 2: Extract polynomial structure in loop momenta *)
  (* Using CoefficientRules: each expression is a polynomial in loop momenta *)
  (* with coefficients that may depend on external momenta, masses, xx, etc. *)

  (* Get the set of all monomials (exponent vectors) appearing *)
  allKeys = DeleteDuplicates[Join @@ (
    Keys[CoefficientRules[#, loopMomenta]] & /@ Join[expandedProps, derivs]
  )];

  (* Separate: loop-dependent monomials (non-zero exponent vector) vs constant *)
  monomialKeys = DeleteCases[allKeys, Table[0, {Length[loopMomenta]}]];
  nMon = Length[monomialKeys];

  If[nMon == 0,
    (* No loop-momentum dependence in any propagator - unusual but handle it *)
    coeffMatrix = Table[0, {nProps}, {nProps}];
    constVector = derivs;
    Return[{coeffMatrix, constVector}];
  ];

  (* Step 3: Build vectors for propagators and derivatives *)
  (* loopVectors[[k]] = vector of coefficients for loop-dependent monomials *)
  loopVectors = Table[
    Table[
      lookupCoeff[CoefficientRules[expandedProps[[k]], loopMomenta], mon],
      {mon, monomialKeys}
    ],
    {k, nProps}
  ];

  derivVectors = Table[
    Table[
      lookupCoeff[CoefficientRules[derivs[[k]], loopMomenta], mon],
      {mon, monomialKeys}
    ],
    {k, nProps}
  ];

  (* Constants: coefficient of the all-zero monomial *)
  propConsts = Table[
    lookupCoeff[CoefficientRules[expandedProps[[k]], loopMomenta],
                Table[0, {Length[loopMomenta]}]],
    {k, nProps}
  ];

  derivConsts = Table[
    lookupCoeff[CoefficientRules[derivs[[k]], loopMomenta],
                Table[0, {Length[loopMomenta]}]],
    {k, nProps}
  ];

  (* Step 4: Solve the linear system *)
  (* For each derivative k: find c_{k1},...,c_{kn} such that *)
  (* derivVector_k = Sum_j c_kj * loopVector_j *)
  (* Then: constVector_k = derivConst_k - Sum_j c_kj * propConst_j *)

  propMatrix = Transpose[loopVectors]; (* nMon x nProps matrix *)

  coeffMatrix = Table[0, {nProps}, {nProps}];
  constVector = Table[0, {nProps}];

  Do[
    If[derivVectors[[k]] === Table[0, {nMon}],
      (* No loop-momentum dependence in this derivative *)
      coeffMatrix[[k]] = Table[0, {nProps}];
      constVector[[k]] = derivConsts[[k]];
    ,
      (* Solve: propMatrix . c = derivVector *)
      sol = Quiet[LinearSolve[propMatrix, derivVectors[[k]]]];
      If[Head[sol] === LinearSolve,
        (* LinearSolve failed - the system is inconsistent *)
        Print["Warning: decomposition failed for propagator ", k,
              ". The propagator set may be incomplete."];
        coeffMatrix[[k]] = Table[0, {nProps}];
        constVector[[k]] = derivs[[k]]; (* leave undecomposed *)
      ,
        coeffMatrix[[k]] = Simplify /@ sol;
        constVector[[k]] = Simplify[derivConsts[[k]] - sol . propConsts];
      ];
    ];
  , {k, nProps}];

  {coeffMatrix, constVector}
];


(* Helper: look up a coefficient in CoefficientRules output *)
lookupCoeff[rules_List, monomial_List] :=
Module[{pos},
  pos = Position[rules[[All, 1]], monomial];
  If[pos === {}, 0, rules[[pos[[1, 1]], 2]]]
];


(* ============================================================ *)
(* FeynmanTrickDecomposition                                     *)
(* Fast path for Feynman trick: only one propagator depends on   *)
(* the Feynman parameter xx, and its derivative is D_i - D_j     *)
(* which decomposes as (1/xx)*D_combined + (-1/xx)*D_j           *)
(* ============================================================ *)

FeynmanTrickDecomposition[nProps_Integer, combinedPos_Integer, otherPos_Integer,
                           param_] :=
Module[{coeffMatrix, constVector},
  coeffMatrix = Table[0, {nProps}, {nProps}];
  constVector = Table[0, {nProps}];

  (* Only the combined propagator has a nonzero derivative *)
  (* d(D_combined)/d(param) = D_i_orig - D_j_orig *)
  (* = (D_combined - (1-param)*D_j)/param - D_j *)
  (* = D_combined/param - D_j/param *)
  coeffMatrix[[combinedPos, combinedPos]] = 1/param;
  coeffMatrix[[combinedPos, otherPos]] = -1/param;

  {coeffMatrix, constVector}
];


(* ============================================================ *)
(* DifferentiatedIntegrals                                       *)
(* Given master index vector v, coefficient matrix and const     *)
(* vector from DecomposePropagatorDerivative, returns all         *)
(* {shifted_index, coefficient} pairs in d/dx I_v                *)
(*                                                               *)
(* d/dx I_v = Sum_k (-v_k) * [                                  *)
(*   Sum_j c_{kj} * I_{v with k+1, j-1}  (from D_j terms)       *)
(*   + c_{k0} * I_{v with k+1}           (from constant)        *)
(* ]                                                             *)
(* ============================================================ *)

DifferentiatedIntegrals[indexVector_List, coeffMatrix_, constVector_] :=
Module[{nProps = Length[indexVector], result = {}, v = indexVector,
        shifted, coeff, k, j},

  Do[
    If[v[[k]] =!= 0,
      (* Terms from D_j contributions: index k goes up, index j goes down *)
      Do[
        coeff = -v[[k]] * coeffMatrix[[k, j]];
        If[coeff =!= 0,
          If[k === j,
            (* k=j: raise and lower same position → net zero change *)
            shifted = v;
          ,
            shifted = v;
            shifted = ReplacePart[shifted, k -> v[[k]] + 1];
            shifted = ReplacePart[shifted, j -> v[[j]] - 1];
          ];
          AppendTo[result, {shifted, Expand[coeff]}];
        ];
      , {j, nProps}];

      (* Term from constant contribution: only index k goes up *)
      coeff = -v[[k]] * constVector[[k]];
      If[coeff =!= 0,
        shifted = ReplacePart[v, k -> v[[k]] + 1];
        AppendTo[result, {shifted, Expand[coeff]}];
      ];
    ];
  , {k, nProps}];

  (* Combine terms with the same shifted index *)
  Module[{grouped = GroupBy[result, First], combined = {}},
    Do[
      AppendTo[combined, {key, Expand[Total[grouped[key][[All, 2]]]]}];
    , {key, Keys[grouped]}];
    DeleteCases[combined, {_, 0}]
  ]
];


(* ============================================================ *)
(* CombinePropagatorSymbolic                                     *)
(* Returns the combined propagator: param*Di + (1-param)*Dj      *)
(* ============================================================ *)

CombinePropagatorSymbolic[Di_, Dj_, param_] := Expand[param * Di + (1 - param) * Dj];


End[];
EndPackage[];
