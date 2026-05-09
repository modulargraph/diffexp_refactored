(* ::Package:: *)

(* DiffExp LocalSeries Subpackage *)
(* Local Fuchsianization and finite-width recursive local solutions. *)

BeginPackage["DiffExp`LocalSeries`", {
  "DiffExp`Symbols`",
  "DiffExp`State`",
  "DiffExp`Utilities`"
}];

LocalZValuation::usage = "LocalZValuation[expr, z] gives the generic z-adic valuation of a rational expression.";
LocalLeadingCoefficient::usage = "LocalLeadingCoefficient[expr, z] gives the leading coefficient at z=0 of a rational expression.";
LocalConnectionMatrix::usage = "LocalConnectionMatrix[M, T, z] transforms theta f = M f by f = T g.";
FuchsianizeLocal::usage = "FuchsianizeLocal[M, z] finds a local meromorphic gauge transform that makes theta g = B g holomorphic in z, when the saturation terminates.";
TrimFuchsianLattice::usage = "TrimFuchsianLattice[M, T, z] removes unnecessary column poles from a Fuchsianizing lattice basis when possible.";
LaurentCoefficientsRational::usage = "LaurentCoefficientsRational[expr, var, {kmin, kmax}] returns rational Laurent coefficients without using asymptotic objects.";
ClearZDenominators::usage = "ClearZDenominators[B, z] returns {q, C} with q theta g = C g and q, C polynomial in z.";
PrepareFiniteWidthData::usage = "PrepareFiniteWidthData[B, z, eps, {epsMin, epsMax}] prepares q/C coefficient dictionaries for RecursiveFiniteWidthSolve.";
RecursiveFiniteWidthSolve::usage = "RecursiveFiniteWidthSolve[qAssoc, CAssoc, dim, alpha, {nMax, epsMin, epsMax, logMax}] solves q theta g = C g in one Frobenius sector.";
SolveLocalFuchsianSeries::usage = "SolveLocalFuchsianSeries[M, z, eps, alpha, {nMax, epsMin, epsMax, logMax}] runs local Fuchsianization, finite-width preparation, and one-sector recursive solving.";
FiniteWidthCoefficient::usage = "FiniteWidthCoefficient[sol, k, n, ell] returns the vector coefficient c[k,n,ell].";
FiniteWidthEvaluate::usage = "FiniteWidthEvaluate[sol, zval, epsval] evaluates a truncated local solution.";
RationalMatrixZEpsLaurentAssoc::usage = "RationalMatrixZEpsLaurentAssoc[M, z, {zMin,zMax}, eps, {epsMin,epsMax}] expands a rational matrix into coefficient matrices.";
ApplyGaugeToSolution::usage = "ApplyGaugeToSolution[gSol, TAssoc, {nMin,nMax}, {epsMin,epsMax}] multiplies a local solution by a finite gauge expansion.";
ValidateFiniteWidthSolution::usage = "ValidateFiniteWidthSolution[sol, qAssoc, CAssoc] checks q theta g - C g coefficient residuals inside the solved window.";

Begin["`Private`"];

infVal = 10^12;

ZeroQGeneric[e_] :=
  TrueQ[Cancel[Together[e]] === 0] || TrueQ[PossibleZeroQ[e]];

RatCancel[e_] := Cancel[Together[e]];

RatCancelVector[v_List] := Map[RatCancel, v];

RatCancelMatrix[m_?MatrixQ] := Map[RatCancel, m, {2}];

AssocGet[assoc_Association, key_, default_] :=
  If[KeyExistsQ[assoc, key], assoc[key], default];

SetAttributes[AssocSet, HoldFirst];
AssocSet[assoc_, key_, value_] := AssociateTo[assoc, key -> value];

UnwrapAssocKey[Key[key_]] := key;
UnwrapAssocKey[key_] := key;
AssocKeys[assoc_Association] := UnwrapAssocKey /@ Keys[assoc];

ThetaExpr[e_, z_] := RatCancel[z D[e, z]];

ThetaMatrix[m_?MatrixQ, z_] := Map[ThetaExpr[#, z] &, m, {2}];

ToDimVector[v_, dim_Integer] := Module[{w},
  w = If[ListQ[v], Flatten[v], {v}];

  If[Length[w] =!= dim,
    Throw[
      Failure[
        "BadVectorDimension",
        <|"Expected" -> dim, "Got" -> Length[w], "Vector" -> v|>
      ]
    ]
  ];

  w
];

PolyMinDegree[p_, x_] := Module[{q, rules},
  q = Expand[p];

  If[ZeroQGeneric[q], Return[Infinity]];

  rules = CoefficientRules[q, {x}];

  If[rules === {}, Infinity, Min[rules[[All, 1, 1]]]]
];

PolyMaxDegree[p_, x_] := Module[{q, rules},
  q = Expand[p];

  If[ZeroQGeneric[q], Return[-Infinity]];

  rules = CoefficientRules[q, {x}];

  If[rules === {}, -Infinity, Max[rules[[All, 1, 1]]]]
];

LocalZValuation[e_, x_] := Module[{r, num, den},
  r = RatCancel[e];

  If[ZeroQGeneric[r], Return[Infinity]];

  num = Numerator[r];
  den = Denominator[r];

  PolyMinDegree[num, x] - PolyMinDegree[den, x]
];

LocalLeadingCoefficient[e_, x_] := Module[{r, num, den, vn, vd},
  r = RatCancel[e];

  If[ZeroQGeneric[r], Return[0]];

  num = Numerator[r];
  den = Denominator[r];

  vn = PolyMinDegree[num, x];
  vd = PolyMinDegree[den, x];

  RatCancel[
    Coefficient[Expand[num], x, vn]/
      Coefficient[Expand[den], x, vd]
  ]
];

MatrixMinValuation[M_?MatrixQ, x_] := Module[{vals},
  vals = DeleteCases[Flatten[Map[LocalZValuation[#, x] &, M, {2}]], Infinity];

  If[vals === {}, Infinity, Min[vals]]
];

ColumnMinValuation[M_?MatrixQ, j_Integer, x_] := Module[{vals},
  vals = DeleteCases[LocalZValuation[#, x] & /@ M[[All, j]], Infinity];

  If[vals === {}, Infinity, Min[vals]]
];

LocalConnectionMatrix[M_?MatrixQ, T_?MatrixQ, z_] :=
  RatCancelMatrix[Inverse[T].(M.T - ThetaMatrix[T, z])];

AdjoinVectorToLattice[T_?MatrixQ, uIn_, z_] := Module[
  {d, u, vals, m, c, pivot, cp, Emat, uNew, TChanged, U},

  d = Length[T];
  u = Flatten[uIn];

  vals = LocalZValuation[#, z] & /@ u;
  m = Min[vals];

  If[m >= 0, Return[T]];

  c = Table[
    If[vals[[i]] === m, LocalLeadingCoefficient[u[[i]], z], 0],
    {i, d}
  ];

  pivot = Select[Range[d], !ZeroQGeneric[c[[#]]] &, 1];

  If[pivot === {},
    Throw[
      Failure[
        "NoPivot",
        <|"Vector" -> u, "Valuations" -> vals|>
      ]
    ]
  ];

  pivot = First[pivot];
  cp = c[[pivot]];

  Emat = IdentityMatrix[d];
  Emat[[pivot, All]] = ConstantArray[0, d];
  Emat[[pivot, pivot]] = RatCancel[1/cp];

  Do[
    If[i =!= pivot,
      Emat[[i, pivot]] = RatCancel[Emat[[i, pivot]] - c[[i]]/cp]
    ],
    {i, d}
  ];

  uNew = RatCancelVector[Emat.u];
  TChanged = RatCancelMatrix[T.Inverse[Emat]];

  U = IdentityMatrix[d];
  U[[All, pivot]] = uNew;

  RatCancelMatrix[TChanged.U]
];

Options[FuchsianizeLocal] = {
  "MaxSteps" -> 200,
  "Verbose" -> False
};

FuchsianizeLocal[M_?MatrixQ, z_, OptionsPattern[]] := Catch@Module[
  {d, T, B, minVal, maxSteps, verbose, step, colVals, bestCol},

  d = Length[M];
  T = IdentityMatrix[d];
  maxSteps = OptionValue["MaxSteps"];
  verbose = TrueQ[OptionValue["Verbose"]];

  For[step = 0, step < maxSteps, step++,
    B = LocalConnectionMatrix[M, T, z];
    minVal = MatrixMinValuation[B, z];

    If[verbose, Print["step = ", step, ", min valuation = ", minVal]];

    If[minVal >= 0,
      Return[{RatCancelMatrix[T], RatCancelMatrix[B]}]
    ];

    colVals = Table[ColumnMinValuation[B, j, z], {j, d}];
    bestCol = First[Ordering[colVals, 1]];

    T = AdjoinVectorToLattice[T, B[[All, bestCol]], z];
  ];

  Failure[
    "FuchsianizationDidNotTerminate",
    <|"MaxSteps" -> maxSteps|>
  ]
];

Options[TrimFuchsianLattice] = {
  "MaxPasses" -> 50,
  "Verbose" -> False
};

TrimFuchsianLattice[M_?MatrixQ, T0_?MatrixQ, z_, OptionsPattern[]] := Module[
  {T, d, maxPasses, verbose, pass, changed, p, Dmat, TTrial, BTrial},

  T = RatCancelMatrix[T0];
  d = Length[T];
  maxPasses = OptionValue["MaxPasses"];
  verbose = TrueQ[OptionValue["Verbose"]];

  For[pass = 1, pass <= maxPasses, pass++,
    changed = False;

    Do[
      Dmat = IdentityMatrix[d];
      Dmat[[p, p]] = z;

      TTrial = RatCancelMatrix[T.Dmat];
      BTrial = LocalConnectionMatrix[M, TTrial, z];

      If[MatrixMinValuation[BTrial, z] >= 0,
        If[verbose, Print["trimmed column ", p]];
        T = TTrial;
        changed = True;
      ],
      {p, d}
    ];

    If[!changed, Break[]];
  ];

  {RatCancelMatrix[T], LocalConnectionMatrix[M, T, z]}
];

LaurentCoefficientsRational[expr_, var_, {kMin_Integer, kMax_Integer}] := Module[
  {r, num, den, vn, vd, v, nShift, dShift, d0, degD,
   coeffs, mMax, m, nm, sum, k},

  If[kMin > kMax, Return[<||>]];

  r = RatCancel[expr];

  If[ZeroQGeneric[r],
    Return[Association@Table[k -> 0, {k, kMin, kMax}]]
  ];

  num = Numerator[r];
  den = Denominator[r];

  vn = PolyMinDegree[num, var];
  vd = PolyMinDegree[den, var];
  v = vn - vd;

  If[kMax < v,
    Return[Association@Table[k -> 0, {k, kMin, kMax}]]
  ];

  nShift = Expand[num/var^vn];
  dShift = Expand[den/var^vd];

  d0 = Coefficient[dShift, var, 0];
  degD = PolyMaxDegree[dShift, var];
  coeffs = <||>;
  mMax = kMax - v;

  For[m = 0, m <= mMax, m++,
    nm = Coefficient[nShift, var, m];
    sum = Sum[
      Coefficient[dShift, var, b] * Lookup[coeffs, v + m - b, 0],
      {b, 1, Min[degD, m]}
    ];
    coeffs[v + m] = RatCancel[(nm - sum)/d0];
  ];

  Association@Table[
    k -> If[k < v, 0, Lookup[coeffs, k, 0]],
    {k, kMin, kMax}
  ]
];

ScalarZEpsLaurentAssoc[expr_, z_, eps_, {rMin_Integer, rMax_Integer}] := Module[
  {e, rules, out, j, coeff, epsCoeffs},

  e = RatCancel[expr];
  out = <||>;

  If[ZeroQGeneric[e], Return[out]];

  If[!PolynomialQ[e, z],
    Throw[Failure["NotPolynomialInZ", <|"Expression" -> e|>]]
  ];

  rules = CoefficientRules[Expand[e], {z}];

  Do[
    j = rule[[1, 1]];
    coeff = rule[[2]];
    epsCoeffs = LaurentCoefficientsRational[coeff, eps, {rMin, rMax}];

    KeyValueMap[
      Function[{epsKey, val},
        Module[{ek = UnwrapAssocKey[epsKey]},
          If[!ZeroQGeneric[val], AssocSet[out, {j, ek}, RatCancel[val]]]
        ]
      ],
      epsCoeffs
    ],
    {rule, rules}
  ];

  out
];

MatrixZEpsLaurentAssoc[M_?MatrixQ, z_, eps_, {rMin_Integer, rMax_Integer}] := Module[
  {rows, cols, out, rules, i, j, zpow, coeff, epsCoeffs, key, mat},

  {rows, cols} = Dimensions[M];
  out = <||>;

  Do[
    If[!ZeroQGeneric[M[[i, j]]],
      If[!PolynomialQ[RatCancel[M[[i, j]]], z],
        Throw[
          Failure[
            "MatrixEntryNotPolynomialInZ",
            <|"Entry" -> {i, j}, "Expression" -> M[[i, j]]|>
          ]
        ]
      ];

      rules = CoefficientRules[Expand[RatCancel[M[[i, j]]]], {z}];

      Do[
        zpow = rule[[1, 1]];
        coeff = rule[[2]];
        epsCoeffs = LaurentCoefficientsRational[coeff, eps, {rMin, rMax}];

        KeyValueMap[
          Function[{epsKey, val},
            Module[{ek = UnwrapAssocKey[epsKey]},
            If[!ZeroQGeneric[val],
              key = {zpow, ek};
              mat = AssocGet[out, key, ConstantArray[0, {rows, cols}]];
              mat[[i, j]] = RatCancel[mat[[i, j]] + val];
              AssocSet[out, key, mat];
            ]
            ]
          ],
          epsCoeffs
        ],
        {rule, rules}
      ];
    ],
    {i, rows}, {j, cols}
  ];

  out
];

AllZCoefficients[expr_, z_] := Module[{e, rules},
  e = RatCancel[expr];

  If[ZeroQGeneric[e], Return[{}]];

  If[!PolynomialQ[e, z],
    Throw[Failure["NotPolynomialInZ", <|"Expression" -> e|>]]
  ];

  rules = CoefficientRules[Expand[e], {z}];
  rules[[All, 2]]
];

AllZCoefficientsMatrix[M_?MatrixQ, z_] :=
  Flatten[AllZCoefficients[#, z] & /@ Flatten[M]];

ClearZDenominators[B_?MatrixQ, z_] := Catch@Module[
  {B0, dens, q, q0, qNorm, C},

  B0 = RatCancelMatrix[B];
  dens = Denominator /@ Flatten[B0];
  q = Fold[PolynomialLCM, 1, dens];
  q = RatCancel[q];
  q0 = Coefficient[Expand[q], z, 0];

  If[ZeroQGeneric[q0],
    Throw[
      Failure[
        "BadZDenominator",
        <|"q" -> q, "Message" -> "The common denominator vanishes at z=0."|>
      ]
    ]
  ];

  qNorm = RatCancel[q/q0];
  C = RatCancelMatrix[qNorm B0];

  If[!AllTrue[Flatten[C], PolynomialQ[#, z] &],
    Throw[Failure["ClearingFailed", <|"q" -> qNorm, "C" -> C|>]]
  ];

  {qNorm, C}
];

Options[PrepareFiniteWidthData] = {
  "EpsilonWorkMin" -> Automatic,
  "EpsilonWorkMax" -> Automatic,
  "EquationEpsilonMin" -> Automatic,
  "EquationEpsilonMax" -> Automatic
};

PrepareFiniteWidthData[
  B_?MatrixQ,
  z_,
  eps_,
  {epsMin_Integer, epsMax_Integer},
  OptionsPattern[]
] := Catch@Module[
  {q, C, coeffs, epsVals, epsPoleMin, epsWorkMin, epsWorkMax,
   eqEpsMin, eqEpsMax, coeffRange, qAssoc, CAssoc},

  {q, C} = ClearZDenominators[B, z];

  coeffs = Join[AllZCoefficients[q, z], AllZCoefficientsMatrix[C, z]];
  epsVals = DeleteCases[LocalZValuation[#, eps] & /@ coeffs, Infinity];
  epsPoleMin = If[epsVals === {}, 0, Min[epsVals]];

  epsWorkMin = OptionValue["EpsilonWorkMin"];
  epsWorkMax = OptionValue["EpsilonWorkMax"];
  eqEpsMin = OptionValue["EquationEpsilonMin"];
  eqEpsMax = OptionValue["EquationEpsilonMax"];

  If[epsWorkMin === Automatic, epsWorkMin = epsMin];
  If[epsWorkMax === Automatic, epsWorkMax = epsMax + Max[0, -epsPoleMin]];
  If[eqEpsMin === Automatic, eqEpsMin = epsMin + Min[0, epsPoleMin]];
  If[eqEpsMax === Automatic, eqEpsMax = epsMax];

  coeffRange = {eqEpsMin - epsWorkMax, eqEpsMax - epsWorkMin};

  qAssoc = ScalarZEpsLaurentAssoc[q, z, eps, coeffRange];
  CAssoc = MatrixZEpsLaurentAssoc[C, z, eps, coeffRange];

  <|
    "q" -> q,
    "C" -> C,
    "qAssoc" -> qAssoc,
    "CAssoc" -> CAssoc,
    "EpsilonWorkRange" -> {epsWorkMin, epsWorkMax},
    "EquationEpsilonRange" -> {eqEpsMin, eqEpsMax},
    "CoefficientEpsilonRange" -> coeffRange,
    "EpsilonValuationMin" -> epsPoleMin
  |>
];

Options[RecursiveFiniteWidthSolve] = {
  "EpsilonWorkMin" -> Automatic,
  "EpsilonWorkMax" -> Automatic,
  "EquationEpsilonMin" -> Automatic,
  "EquationEpsilonMax" -> Automatic,
  "Prescribed" -> <||>,
  "Tolerance" -> 10^-10,
  "CheckResidual" -> True
};

Options[SolveLocalFuchsianSeries] = {
  "MaxSteps" -> 200,
  "Verbose" -> False,
  "EpsilonWorkMin" -> Automatic,
  "EpsilonWorkMax" -> Automatic,
  "EquationEpsilonMin" -> Automatic,
  "EquationEpsilonMax" -> Automatic,
  "Prescribed" -> <||>,
  "Tolerance" -> 10^-10,
  "CheckResidual" -> True,
  "ReturnOriginalSolution" -> True,
  "GaugeZMin" -> Automatic,
  "GaugeZMax" -> Automatic,
  "GaugeEpsilonRange" -> Automatic
};

RecursiveFiniteWidthSolve[
  qAssoc_Association,
  CAssoc_Association,
  dim_Integer,
  alpha_,
  {nMax_Integer, epsMin_Integer, epsMax_Integer, logMax_Integer},
  OptionsPattern[]
] := Catch@Module[
  {epsDegrees, rMin, epsWorkMin, epsWorkMax, eqEpsMin, eqEpsMax,
   epsOrders, eqOrders, epsIndex, blockSize, id, prescribed, coeff,
   idx, getCoeff, solveBlock, n, numEq, b, A, eqPos, eqK,
   ell, row0, addContribution, j, r, qjr, Cjr, ns, kSrc, mat,
   fixed, x, tol, checkResidual},

  If[Length[qAssoc] == 0,
    Throw[Failure["MissingQ", <|"Message" -> "qAssoc must contain at least one coefficient."|>]]
  ];

  epsDegrees = Join[
    AssocKeys[qAssoc][[All, 2]],
    If[Length[CAssoc] == 0, {}, AssocKeys[CAssoc][[All, 2]]]
  ];
  rMin = Min[epsDegrees];

  epsWorkMin = OptionValue["EpsilonWorkMin"];
  epsWorkMax = OptionValue["EpsilonWorkMax"];
  eqEpsMin = OptionValue["EquationEpsilonMin"];
  eqEpsMax = OptionValue["EquationEpsilonMax"];

  If[epsWorkMin === Automatic, epsWorkMin = epsMin];
  If[epsWorkMax === Automatic, epsWorkMax = epsMax + Max[0, -rMin]];
  If[eqEpsMin === Automatic, eqEpsMin = epsMin + Min[0, rMin]];
  If[eqEpsMax === Automatic, eqEpsMax = epsMax];

  epsOrders = Range[epsWorkMin, epsWorkMax];
  eqOrders = Range[eqEpsMin, eqEpsMax];
  epsIndex = AssociationThread[epsOrders -> Range[Length[epsOrders]]];
  blockSize = Length[epsOrders] (logMax + 1) dim;
  id = IdentityMatrix[dim];
  prescribed = Association[OptionValue["Prescribed"]];
  tol = OptionValue["Tolerance"];
  checkResidual = TrueQ[OptionValue["CheckResidual"]];
  coeff = <||>;

  idx[k_, l_, comp_] :=
    1 + ((epsIndex[k] - 1) (logMax + 1) + l) dim + (comp - 1);

  getCoeff[k_, n0_, l_] := Module[{key = {k, n0, l}},
    If[n0 < 0 || l < 0 || l > logMax,
      ConstantArray[0, dim],
      Which[
        KeyExistsQ[coeff, key], AssocGet[coeff, key, ConstantArray[0, dim]],
        KeyExistsQ[prescribed, key], ToDimVector[AssocGet[prescribed, key, ConstantArray[0, dim]], dim],
        True, ConstantArray[0, dim]
      ]
    ]
  ];

  solveBlock[matA_, rhs_, fixed_Association, nNow_Integer] := Module[
    {fixedIdx, fixedVals, freeIdx, x0, rhs2, Afree, sol, residual, resList, rel},

    fixedIdx = Sort[Keys[fixed]];
    fixedVals = Lookup[fixed, fixedIdx];
    freeIdx = Complement[Range[blockSize], fixedIdx];
    x0 = ConstantArray[0, blockSize];

    If[fixedIdx =!= {},
      x0[[fixedIdx]] = fixedVals;
      rhs2 = rhs - matA[[All, fixedIdx]].fixedVals,
      rhs2 = rhs
    ];

    If[freeIdx =!= {},
      Afree = matA[[All, freeIdx]];
      sol = Quiet@Check[
        If[Dimensions[Afree][[1]] == Dimensions[Afree][[2]],
          LinearSolve[Afree, rhs2],
          LeastSquares[Afree, rhs2]
        ],
        LeastSquares[Afree, rhs2]
      ];
      x0[[freeIdx]] = sol;
    ];

    If[checkResidual,
      residual = matA.x0 - rhs;
      resList = Normal[residual];

      If[VectorQ[N[resList], NumericQ] && VectorQ[N[rhs], NumericQ],
        rel = Norm[N[resList]]/(1 + Norm[N[rhs]]);
        If[rel > tol,
          Throw[
            Failure[
              "BlockResidualTooLarge",
              <|"n" -> nNow, "RelativeResidual" -> rel, "Tolerance" -> tol|>
            ]
          ]
        ],
        If[!AllTrue[resList, ZeroQGeneric],
          Throw[Failure["BlockResidualNonzero", <|"n" -> nNow, "Residual" -> resList|>]]
        ]
      ];
    ];

    x0
  ];

  Do[
    numEq = Length[eqOrders] (logMax + 1) dim;
    A = ConstantArray[0, {numEq, blockSize}];
    b = ConstantArray[0, numEq];

    addContribution[matIn_, kSource_, nSource_, lSource_, rowStart_] := Module[
      {known, rowRange, a, c, pos, val},

      If[nSource < 0 || lSource < 0 || lSource > logMax, Return[]];

      If[nSource == n && KeyExistsQ[epsIndex, kSource],
        Do[
          val = matIn[[a, c]];
          If[!ZeroQGeneric[val],
            A[[rowStart + a - 1, idx[kSource, lSource, c]]] =
              RatCancel[A[[rowStart + a - 1, idx[kSource, lSource, c]]] + val];
          ],
          {a, dim}, {c, dim}
        ],
        known = getCoeff[kSource, nSource, lSource];
        rowRange = rowStart ;; rowStart + dim - 1;
        b[[rowRange]] = RatCancelVector[b[[rowRange]] - matIn.known];
      ];
    ];

    Do[
      eqK = eqOrders[[eqPos]];

      Do[
        row0 = 1 + ((eqPos - 1) (logMax + 1) + ell) dim;

        KeyValueMap[
          Function[{key0, val0},
            Module[{kk = UnwrapAssocKey[key0]},
            j = kk[[1]];
            r = kk[[2]];
            qjr = val0;
            ns = n - j;

            If[ns >= 0,
              kSrc = eqK - r;
              mat = qjr (alpha + ns) id;
              addContribution[mat, kSrc, ns, ell, row0];

              If[ell + 1 <= logMax,
                mat = qjr (ell + 1) id;
                addContribution[mat, kSrc, ns, ell + 1, row0];
              ];
            ];
            ]
          ],
          qAssoc
        ];

        KeyValueMap[
          Function[{key0, val0},
            Module[{kk = UnwrapAssocKey[key0]},
            j = kk[[1]];
            r = kk[[2]];
            Cjr = val0;
            ns = n - j;

            If[ns >= 0,
              kSrc = eqK - r;
              addContribution[-Cjr, kSrc, ns, ell, row0];
            ];
            ]
          ],
          CAssoc
        ],
        {ell, 0, logMax}
      ],
      {eqPos, Length[eqOrders]}
    ];

    fixed = <||>;

    KeyValueMap[
      Function[{key0, val0},
        Module[{kk = UnwrapAssocKey[key0], k0, n0, l0, vec},
          If[Length[kk] == 3,
            k0 = kk[[1]];
            n0 = kk[[2]];
            l0 = kk[[3]];

            If[n0 == n && KeyExistsQ[epsIndex, k0] && TrueQ[0 <= l0 <= logMax],
              vec = ToDimVector[val0, dim];
              Do[fixed[idx[k0, l0, comp]] = vec[[comp]], {comp, dim}];
            ];
          ];
        ]
      ],
      prescribed
    ];

    x = solveBlock[A, b, fixed, n];

    Do[
      AssocSet[
        coeff,
        {k, n, ell},
        RatCancelVector[Table[x[[idx[k, ell, comp]]], {comp, dim}]]
      ],
      {k, epsOrders}, {ell, 0, logMax}
    ],
    {n, 0, nMax}
  ];

  <|
    "Alpha" -> alpha,
    "Coeff" -> coeff,
    "Dim" -> dim,
    "N" -> nMax,
    "ZRange" -> {0, nMax},
    "LogMax" -> logMax,
    "EpsRange" -> {epsMin, epsMax},
    "EpsilonWorkRange" -> {epsWorkMin, epsWorkMax},
    "EquationEpsilonRange" -> {eqEpsMin, eqEpsMax}
  |>
];

SolveLocalFuchsianSeries[
  M_?MatrixQ,
  z_,
  eps_,
  alpha_,
  {nMax_Integer, epsMin_Integer, epsMax_Integer, logMax_Integer},
  OptionsPattern[]
] := Catch@Module[
  {fb, T, B, prep, gSol, tMin, tMax, tEpsRange, TAssoc, fSol,
   tVals},

  fb = FuchsianizeLocal[
    M,
    z,
    "MaxSteps" -> OptionValue["MaxSteps"],
    "Verbose" -> OptionValue["Verbose"]
  ];

  If[MatchQ[fb, _Failure], Return[fb]];

  {T, B} = fb;

  prep = PrepareFiniteWidthData[
    B,
    z,
    eps,
    {epsMin, epsMax},
    "EpsilonWorkMin" -> OptionValue["EpsilonWorkMin"],
    "EpsilonWorkMax" -> OptionValue["EpsilonWorkMax"],
    "EquationEpsilonMin" -> OptionValue["EquationEpsilonMin"],
    "EquationEpsilonMax" -> OptionValue["EquationEpsilonMax"]
  ];

  If[MatchQ[prep, _Failure], Return[prep]];

  gSol = RecursiveFiniteWidthSolve[
    prep["qAssoc"],
    prep["CAssoc"],
    Length[M],
    alpha,
    {nMax, epsMin, epsMax, logMax},
    "EpsilonWorkMin" -> prep["EpsilonWorkRange"][[1]],
    "EpsilonWorkMax" -> prep["EpsilonWorkRange"][[2]],
    "EquationEpsilonMin" -> prep["EquationEpsilonRange"][[1]],
    "EquationEpsilonMax" -> prep["EquationEpsilonRange"][[2]],
    "Prescribed" -> OptionValue["Prescribed"],
    "Tolerance" -> OptionValue["Tolerance"],
    "CheckResidual" -> OptionValue["CheckResidual"]
  ];

  If[MatchQ[gSol, _Failure], Return[gSol]];

  fSol = Missing["NotRequested"];
  TAssoc = <||>;

  If[TrueQ[OptionValue["ReturnOriginalSolution"]],
    tMin = OptionValue["GaugeZMin"];
    tMax = OptionValue["GaugeZMax"];
    tEpsRange = OptionValue["GaugeEpsilonRange"];

    If[tMin === Automatic,
      tVals = DeleteCases[LocalZValuation[#, z] & /@ Flatten[T], Infinity];
      tMin = If[tVals === {}, 0, Min[tVals]];
    ];
    If[tMax === Automatic, tMax = nMax];
    If[tEpsRange === Automatic,
      tEpsRange = {
        epsMin - gSol["EpsilonWorkRange"][[2]],
        epsMax - gSol["EpsilonWorkRange"][[1]]
      };
    ];

    TAssoc = RationalMatrixZEpsLaurentAssoc[T, z, {tMin, tMax}, eps, tEpsRange];
    fSol = ApplyGaugeToSolution[gSol, TAssoc, {tMin, nMax}, {epsMin, epsMax}];

    If[MatchQ[fSol, _Failure], Return[fSol]];
  ];

  <|
    "Gauge" -> T,
    "FuchsianMatrix" -> B,
    "Preparation" -> prep,
    "RegularSolution" -> gSol,
    "GaugeAssoc" -> TAssoc,
    "OriginalSolution" -> fSol
  |>
];

FiniteWidthCoefficient[sol_Association, k_Integer, n_Integer, ell_Integer] :=
  AssocGet[sol["Coeff"], {k, n, ell}, ConstantArray[0, sol["Dim"]]];

Options[FiniteWidthEvaluate] = {
  "UseWorkRange" -> False
};

FiniteWidthEvaluate[sol_Association, zval_, epsval_, OptionsPattern[]] := Module[
  {alpha, dim, logMax, epsRange, zRange, out},

  alpha = sol["Alpha"];
  dim = sol["Dim"];
  logMax = sol["LogMax"];
  epsRange = If[TrueQ[OptionValue["UseWorkRange"]], sol["EpsilonWorkRange"], sol["EpsRange"]];
  zRange = Lookup[sol, "ZRange", {0, sol["N"]}];
  out = ConstantArray[0, dim];

  Do[
    out = out + epsval^k zval^(alpha + n) Log[zval]^ell FiniteWidthCoefficient[sol, k, n, ell],
    {k, epsRange[[1]], epsRange[[2]]},
    {n, zRange[[1]], zRange[[2]]},
    {ell, 0, logMax}
  ];

  out
];

RationalMatrixZEpsLaurentAssoc[
  M_?MatrixQ,
  z_,
  {zMin_Integer, zMax_Integer},
  eps_,
  {epsMin_Integer, epsMax_Integer}
] := Module[
  {rows, cols, out, i, j, zCoeffs, epsCoeffs, key, mat},

  {rows, cols} = Dimensions[M];
  out = <||>;

  Do[
    zCoeffs = LaurentCoefficientsRational[M[[i, j]], z, {zMin, zMax}];

    KeyValueMap[
      Function[{zp, zCoeff},
        Module[{zpk = UnwrapAssocKey[zp]},
        epsCoeffs = LaurentCoefficientsRational[zCoeff, eps, {epsMin, epsMax}];

        KeyValueMap[
          Function[{ep, val},
            Module[{epk = UnwrapAssocKey[ep]},
            If[!ZeroQGeneric[val],
              key = {zpk, epk};
              mat = AssocGet[out, key, ConstantArray[0, {rows, cols}]];
              mat[[i, j]] = RatCancel[mat[[i, j]] + val];
              AssocSet[out, key, mat];
            ]
            ]
          ],
          epsCoeffs
        ];
        ]
      ],
      zCoeffs
    ],
    {i, rows}, {j, cols}
  ];

  out
];

ApplyGaugeToSolution[
  gSol_Association,
  TAssoc_Association,
  {nMin_Integer, nMax_Integer},
  {epsMin_Integer, epsMax_Integer}
] := Catch@Module[
  {firstMat, dimOut, dimIn, alpha, logMax, coeff, vec, a, b, mat, gvec},

  If[Length[TAssoc] == 0, Return[gSol]];

  firstMat = First[Values[TAssoc]];
  {dimOut, dimIn} = Dimensions[firstMat];

  If[dimIn =!= gSol["Dim"],
    Throw[
      Failure[
        "GaugeDimensionMismatch",
        <|"GaugeInputDim" -> dimIn, "SolutionDim" -> gSol["Dim"]|>
      ]
    ]
  ];

  alpha = gSol["Alpha"];
  logMax = gSol["LogMax"];
  coeff = <||>;

  Do[
    vec = ConstantArray[0, dimOut];

    KeyValueMap[
      Function[{key0, mat0},
        Module[{kk = UnwrapAssocKey[key0]},
        a = kk[[1]];
        b = kk[[2]];
        mat = mat0;
        gvec = FiniteWidthCoefficient[gSol, k - b, n - a, ell];
        vec = vec + mat.gvec;
        ]
      ],
      TAssoc
    ];

    If[!AllTrue[vec, ZeroQGeneric],
      AssocSet[coeff, {k, n, ell}, RatCancelVector[vec]]
    ],
    {k, epsMin, epsMax}, {n, nMin, nMax}, {ell, 0, logMax}
  ];

  Join[
    gSol,
    <|
      "Coeff" -> coeff,
      "Dim" -> dimOut,
      "ZRange" -> {nMin, nMax},
      "EpsRange" -> {epsMin, epsMax},
      "EpsilonWorkRange" -> {epsMin, epsMax}
    |>
  ]
];

Options[ValidateFiniteWidthSolution] = {
  "EpsilonRange" -> Automatic,
  "ZRange" -> Automatic,
  "LogMax" -> Automatic
};

ValidateFiniteWidthSolution[sol_Association, qAssoc_Association, CAssoc_Association, OptionsPattern[]] := Module[
  {dim, alpha, epsRange, zRange, logMax, residuals, maxResidual,
   coeff, n, eqK, ell, res, j, r, qjr, Cjr, ns, kSrc, val},

  dim = sol["Dim"];
  alpha = sol["Alpha"];
  epsRange = OptionValue["EpsilonRange"];
  zRange = OptionValue["ZRange"];
  logMax = OptionValue["LogMax"];

  If[epsRange === Automatic, epsRange = sol["EquationEpsilonRange"]];
  If[zRange === Automatic, zRange = sol["ZRange"]];
  If[logMax === Automatic, logMax = sol["LogMax"]];

  coeff[k_, n0_, l_] := FiniteWidthCoefficient[sol, k, n0, l];
  residuals = <||>;
  maxResidual = 0;

  Do[
    res = ConstantArray[0, dim];

    KeyValueMap[
      Function[{key0, val0},
        Module[{kk = UnwrapAssocKey[key0]},
        j = kk[[1]];
        r = kk[[2]];
        qjr = val0;
        ns = n - j;

        If[ns >= 0,
          kSrc = eqK - r;
          res = res + qjr (alpha + ns) coeff[kSrc, ns, ell];
          If[ell + 1 <= sol["LogMax"],
            res = res + qjr (ell + 1) coeff[kSrc, ns, ell + 1];
          ];
        ];
        ]
      ],
      qAssoc
    ];

    KeyValueMap[
      Function[{key0, val0},
        Module[{kk = UnwrapAssocKey[key0]},
        j = kk[[1]];
        r = kk[[2]];
        Cjr = val0;
        ns = n - j;

        If[ns >= 0,
          kSrc = eqK - r;
          res = res - Cjr.coeff[kSrc, ns, ell];
        ];
        ]
      ],
      CAssoc
    ];

    val = RatCancelVector[res];
    AssocSet[residuals, {eqK, n, ell}, val];

    If[VectorQ[N[val], NumericQ],
      maxResidual = Max[maxResidual, Norm[N[val]]]
    ],
    {n, zRange[[1]], zRange[[2]]},
    {eqK, epsRange[[1]], epsRange[[2]]},
    {ell, 0, logMax}
  ];

  <|"Residuals" -> residuals, "MaxResidual" -> maxResidual|>
];

End[];

EndPackage[];
