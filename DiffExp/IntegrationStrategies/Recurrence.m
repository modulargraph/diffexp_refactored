(* IntegrationStrategies/Recurrence.m *)
(* Rational recurrence and singular recurrence strategies *)

(* ============================================================================ *)
(* SHARED HELPERS                                                              *)
(* ============================================================================ *)

(* Core rationalization of the factored A matrix.
   Computes the common denominator and numerator matrix.
   Returns {True, AMatTogether, Dpoly, NAMat, maxPolyDeg} on success,
   or {False} if rationalization is not possible or not beneficial.
   The caller specifies whether a pole (x factor) is expected via requirePole. *)
RationalizeAMatrixCore[ctx_Association, requirePole:(True|False)] := Module[
  {AMat, maxOrd, result},

  maxOrd = ctx["ExpansionOrder"];

  (* Only attempt if factored matrix is available *)
  If[MissingQ[ctx["AMatFactored"]],
    Return[{False}]
  ];

  AMat = ctx["AMatFactored"];

  result = TimeConstrained[
    Module[{at, fe, dn, dp, naMat, degA},
      at = Map[Together, AMat, {2}];
      fe = Flatten[at];

      (* Check all entries are rational in x *)
      If[!And @@ (PolynomialQ[Numerator[#], DiffExp`Symbols`x] &&
                  PolynomialQ[Denominator[#], DiffExp`Symbols`x] & /@ fe),
        Return[{False}]
      ];

      (* Compute common denominator *)
      dn = DeleteDuplicates[Denominator /@ fe];
      dn = DeleteCases[dn, a_ /; FreeQ[a, DiffExp`Symbols`x]];

      If[requirePole,
        (* Singular case: denominator must have x-dependence *)
        If[dn === {}, Return[{False}]];
        dp = PolynomialLCM @@ dn;
        (* Must have x as a factor (pole) *)
        If[Exponent[dp, DiffExp`Symbols`x, Min] < 1, Return[{False}]];
        (* Factor out x: use the full x*D(x) for the numerator *)
        naMat = Expand[dp * at];
        (* Now dp becomes D(x) without the x factor *)
        dp = Cancel[dp / DiffExp`Symbols`x];
        ,
        (* Non-singular case *)
        If[dn === {},
          dp = 1;,
          dp = PolynomialLCM @@ dn;
        ];
        (* D(0) must be nonzero *)
        If[(dp /. DiffExp`Symbols`x -> 0) === 0, Return[{False}]];
        naMat = Expand[dp * at];
      ];

      (* For the singular case, check D(0) != 0 after factoring out x *)
      If[requirePole,
        If[(dp /. DiffExp`Symbols`x -> 0) === 0, Return[{False}]];
      ];

      (* Check that polynomial degree is beneficial: degA < maxOrd/2 *)
      degA = Max[0, Max[Exponent[#, DiffExp`Symbols`x] & /@ Flatten[naMat]]];
      If[degA > maxOrd/2, Return[{False}]];

      {True, at, dp, naMat, degA}
    ],
    5.0, (* 5 second timeout *)
    {False}
  ];

  result
];

(* Extract matrix coefficients from the expanded A matrix at specified order range.
   Returns a list of matrices: result[[j+1]] = coefficient of x^(minOrder+j) in AMatExpanded.
   For non-singular: minOrder=0, maxOrder=maxOrd-1 gives A_0, A_1, ..., A_{maxOrd-1}
   For singular (M coefficients): minOrder=-1, maxOrder=maxOrd-1 gives A_{-1}, A_0, ..., A_{maxOrd-1} *)
ExtractAMatCoefficients[ctx_Association, minOrder_Integer, maxOrder_Integer] := Module[
  {AMatExpanded, systemSize, coeffs},

  systemSize = ctx["SystemSize"];
  AMatExpanded = ctx["AMatExpanded"];

  coeffs = Table[
    Table[
      SeriesCoefficient[AMatExpanded[[i, k]], {DiffExp`Symbols`x, 0, j}],
      {i, systemSize}, {k, systemSize}
    ],
    {j, minOrder, maxOrder}
  ];

  coeffs
];

(* Prepare eigenvalue data for the singular recurrence method.
   Combines the applicability check with eigenvalue computation.
   Returns {eigenvalues, P, PInv} on success, or $Failed if not applicable. *)
PrepareSingularRecurrence[ctx_Association] := Quiet[Check[
  Module[{AMatExpanded, minOrders, minOrder, residueMat, systemSize,
          eigenvalues, eigenvectors, diffs, P, PInv},

    systemSize = ctx["SystemSize"];
    AMatExpanded = ctx["AMatExpanded"];

    (* Determine minimum order of the matrix entries *)
    minOrders = Flatten[Table[
      If[Head[AMatExpanded[[i, j]]] === SeriesData,
        AMatExpanded[[i, j]][[4]] / AMatExpanded[[i, j]][[6]],
        0 (* constant entries have order 0 *)
      ],
      {i, systemSize}, {j, systemSize}
    ]];
    minOrder = Min[minOrders];

    (* Must be exactly -1 for a simple pole *)
    If[minOrder =!= -1, Return[$Failed]];

    (* Extract residue matrix A_{-1} *)
    residueMat = Table[
      SeriesCoefficient[AMatExpanded[[i, j]], {DiffExp`Symbols`x, 0, -1}],
      {i, systemSize}, {j, systemSize}
    ];

    (* Residue must not be identically zero *)
    If[DiffExp`Utilities`PChop[residueMat] === ConstantArray[0, {systemSize, systemSize}],
      Return[$Failed]
    ];

    (* Compute eigenvalues *)
    eigenvalues = Eigenvalues[N[residueMat, ctx["WorkingPrecision"]]];
    eigenvalues = Rationalize[eigenvalues, 10^(-ctx["ChopPrecision"]/2)];

    (* Degenerate (snapped-equal) eigenvalues are Jordan/resonant territory:
       float-level noise splits them by roughly the square root of the
       noise, so the eigenvector matrix looks full-rank to any reasonable
       tolerance while being catastrophically ill-conditioned.  Reject and
       let the general singular recurrence handle the block with proper
       generalized-eigenvector chains. *)
    If[Length[DeleteDuplicates[eigenvalues]] < systemSize,
      Return[$Failed]
    ];

    (* Check diagonalizability: rank of eigenvector matrix must equal system size *)
    eigenvectors = Eigenvectors[N[residueMat, ctx["WorkingPrecision"]]];
    If[MatrixRank[eigenvectors, Tolerance -> 10^(-ctx["ChopPrecision"]/2)] < systemSize,
      Return[$Failed]
    ];

    (* Check non-resonance: no difference lambda_j - lambda_i is a positive integer *)
    diffs = Flatten[Table[
      eigenvalues[[j]] - eigenvalues[[i]],
      {i, systemSize}, {j, systemSize}
    ]];
    If[AnyTrue[diffs, (IntegerQ[#] && # > 0) &],
      Return[$Failed]
    ];

    (* Diagonalize: P = eigenvectors as columns *)
    P = Transpose[eigenvectors];
    PInv = Inverse[P];
    P = N[P, ctx["WorkingPrecision"]];
    PInv = N[PInv, ctx["WorkingPrecision"]];

    {eigenvalues, P, PInv}
  ],
  $Failed (* Return $Failed on any error *)
]];

(* ============================================================================ *)
(* RATIONAL RECURRENCE METHOD                                                  *)
(* ============================================================================ *)

(* Check if the rational recurrence method is applicable:
   The expansion point must be non-singular (series starts at x^0).
   Checks the factored matrix in ctx for the given integral block. *)
RationalRecurrenceApplicableQ[ctx_Association] := Quiet[Check[
  Module[{AMat, atZero},

    (* Check that factored matrix is available *)
    If[MissingQ[ctx["AMatFactored"]],
      Return[False]
    ];

    (* Non-singular iff the factored matrix is finite at x=0 *)
    AMat = ctx["AMatFactored"];
    atZero = AMat /. DiffExp`Symbols`x -> 0;

    TrueQ[And @@ (NumericQ /@ Flatten[atZero])]
  ],
  False (* Return False on any error *)
]];

(* SingularRecurrenceApplicableQ is now a thin wrapper around PrepareSingularRecurrence.
   Returns True if the singular recurrence method is applicable. *)
SingularRecurrenceApplicableQ[ctx_Association] :=
  PrepareSingularRecurrence[ctx] =!= $Failed;

(* Try to rationalize the factored A matrix and extract polynomial coefficients.
   Returns {True, dCoeffs, aCoeffs, dD, dA, d0, d0Inv} on success,
   or {False} if rationalization is too expensive or not possible. *)
TryRationalizeMatrix[ctx_Association] := Module[
  {result, Dpoly, NAMat, dA, dD, dCoeffs, aCoeffs, d0, d0Inv},

  result = RationalizeAMatrixCore[ctx, False];
  If[!result[[1]], Return[{False}]];

  {Dpoly, NAMat, dA} = result[[{3, 4, 5}]];

  (* Extract polynomial coefficients *)
  dD = Exponent[Dpoly, DiffExp`Symbols`x];
  dCoeffs = Table[
    N[Coefficient[Dpoly, DiffExp`Symbols`x, i], ctx["WorkingPrecision"]],
    {i, 0, dD}
  ];

  aCoeffs = Table[
    Map[N[Coefficient[#, DiffExp`Symbols`x, j], ctx["WorkingPrecision"]] &, NAMat, {2}],
    {j, 0, dA}
  ];

  d0 = dCoeffs[[1]];
  d0Inv = 1/d0;

  {True, dCoeffs, aCoeffs, dD, dA, d0, d0Inv}
];

(* Extract series coefficients of the A matrix from the expanded form.
   Returns {aCoeffs, dA} where dA = maxOrd-1 (all orders used). *)
ExtractSeriesCoefficients[ctx_Association] := Module[{maxOrd, aCoeffs},
  maxOrd = ctx["ExpansionOrder"];
  aCoeffs = ExtractAMatCoefficients[ctx, 0, maxOrd - 1];
  {aCoeffs, maxOrd - 1}
];

(* Core recurrence computation: computes the fundamental matrix (homogeneous solutions). *)
ComputeFundamentalMatrix[ctx_Association, aCoeffs_, dCoeffs_, dD_, dA_, d0Inv_] := Module[
  {systemSize, maxOrd, fCoeffs, FMat},

  systemSize = ctx["SystemSize"];
  maxOrd = ctx["ExpansionOrder"];

  FMat = Table[
    fCoeffs = ConstantArray[
      ConstantArray[N[0, ctx["WorkingPrecision"]], systemSize],
      maxOrd + 1
    ];
    fCoeffs[[1]] = N[UnitVector[systemSize, col], ctx["WorkingPrecision"]];

    Do[
      fCoeffs[[m + 2]] = d0Inv/(m + 1) * (
        Sum[aCoeffs[[j + 1]] . fCoeffs[[m - j + 1]], {j, 0, Min[m, dA]}] -
        If[dD >= 1,
          Sum[dCoeffs[[i + 1]] (m - i + 1) fCoeffs[[m - i + 2]], {i, 1, Min[m, dD]}],
          0
        ]
      );
      , {m, 0, maxOrd - 1}
    ];

    Table[
      SeriesData[DiffExp`Symbols`x, 0, fCoeffs[[All, k]], 0, maxOrd + 1, 1],
      {k, systemSize}
    ]
    , {col, systemSize}
  ] // Transpose;

  DiffExp`Utilities`PChop[FMat]
];

(* Rational recurrence solver for f'(x) = A(x)f(x) + B(x) at non-singular points.
   Two modes:
   1. Denominator clearing: if A(x) is cheaply rationalizable, clears denominators
      to get D(x)f' = N_A(x)f + N_B(x). Recurrence cost is O(N*dA) per solution.
   2. Direct series: uses the expanded series coefficients A_j directly.
      Recurrence cost is O(N^2) per solution but avoids rationalization overhead.
   Both modes bypass the Frobenius/Wronskian machinery. *)
SolveRationalRecurrence[ctx_Association, bVec_, epsord_, cacheIn_Association] := Module[
	  {systemSize, maxOrd, fCoeffs, bSeriesCoeffs, bCoeffAtN, nbCoeff,
	   FMat, fParticular, fGeneral, cIndices, c,
	   dCoeffs, aCoeffs, dD, dA, d0, d0Inv, rationalResult,
	   bMaxLogK, zeroVec, logDerivativeCoeff,
	   cache = cacheIn},

  systemSize = ctx["SystemSize"];
  maxOrd = ctx["ExpansionOrder"];

  If[epsord === 0 && !KeyExistsQ[cache, "RR"],
    (* Try the denominator-clearing approach first *)
    rationalResult = TryRationalizeMatrix[ctx];

    If[rationalResult[[1]],
      (* Denominator-clearing mode: polynomial recurrence *)
      {dCoeffs, aCoeffs, dD, dA, d0, d0Inv} = rationalResult[[{2, 3, 4, 5, 6, 7}]];
      DiffExp`Utilities`PrintInfo["Using rational recurrence (poly degree ", dA, ") for integrals ", ctx["Label"], "."][3];
      ,
      (* Direct series mode: use expanded series coefficients *)
      {aCoeffs, dA} = ExtractSeriesCoefficients[ctx];
      dD = 0;
      dCoeffs = {N[1, ctx["WorkingPrecision"]]};
      d0Inv = N[1, ctx["WorkingPrecision"]];
      DiffExp`Utilities`PrintInfo["Using series recurrence for integrals ", ctx["Label"], "."][3];
    ];

    (* Compute fundamental matrix *)
    FMat = ComputeFundamentalMatrix[ctx, aCoeffs, dCoeffs, dD, dA, d0Inv];

    cache["RR"] = {FMat, dCoeffs, aCoeffs, dD, dA, d0Inv};
  ];

  {FMat, dCoeffs, aCoeffs, dD, dA, d0Inv} = cache["RR"];

  (* Compute particular solution using bVec *)
  If[DiffExp`Utilities`PChop[bVec] === ConstantArray[0, systemSize],
    fParticular = ConstantArray[
      SeriesData[DiffExp`Symbols`x, 0, {}, 0, maxOrd + 1, 1],
      systemSize
    ];
    ,
	    bMaxLogK = Max[0, Max[Table[
	      DiffExp`SeriesOps`MaxLogxPower[bVec[[k]]],
	      {k, systemSize}
	    ]]];
	    zeroVec = ConstantArray[N[0, ctx["WorkingPrecision"]], systemSize];

	    If[bMaxLogK > DiffExp`State`IMaxLogOrder,
	      DiffExp`Integration`UpdateIntReps[bMaxLogK];
	    ];

	    (* Extract series coefficients from each Logx sector of bVec. *)
	    bSeriesCoeffs = Table[
	      bCoeffAtN = SeriesCoefficient[
	        DiffExp`SeriesOps`LogxCoeff[bVec[[k]], logk],
	        {DiffExp`Symbols`x, 0, n}
	      ];
	      If[NumericQ[bCoeffAtN], bCoeffAtN, DiffExp`Utilities`PChop[bCoeffAtN]],
	      {n, 0, maxOrd - 1}, {logk, 0, bMaxLogK}, {k, systemSize}
	    ];

	    (* Run recurrence with f_0 = 0 for the particular solution.
	       Coefficients are indexed as fCoeffs[[n+1, logk+1, comp]].
	       Solve log powers top-down because d/dx Logx^(k+1) contributes to
	       the Logx^k equation at the same power of x. *)
	    fCoeffs = ConstantArray[
	      ConstantArray[
	        ConstantArray[N[0, ctx["WorkingPrecision"]], systemSize],
	        bMaxLogK + 1
	      ],
	      maxOrd + 1
	    ];

	    Do[
	      Do[
	        (* Compute N_B coefficient at order m and log power logk:
	           nb_{m,k} = sum_i d_i * B_{m-i,k}. *)
	        nbCoeff = Sum[
	          If[m - i <= maxOrd - 1,
	            dCoeffs[[i + 1]] * bSeriesCoeffs[[m - i + 1, logk + 1]],
	            zeroVec
	          ],
	          {i, 0, Min[m, dD]}
	        ];

	        logDerivativeCoeff = If[logk < bMaxLogK,
	          Sum[
	            If[0 <= m - i + 1 <= maxOrd,
	              dCoeffs[[i + 1]] * (logk + 1) *
	                fCoeffs[[m - i + 2, logk + 2]],
	              zeroVec
	            ],
	            {i, 0, Min[m + 1, dD]}
	          ],
	          zeroVec
	        ];

	        (* D f' = N_A f + D B.  Isolate d0*(m+1)*f_{m+1,k};
	           all lower denominator-derivative terms and Logx-derivative
	           mixing move to the right-hand side with a minus sign. *)
	        fCoeffs[[m + 2, logk + 1]] = d0Inv/(m + 1) * (
	          nbCoeff +
	          Sum[aCoeffs[[j + 1]] . fCoeffs[[m - j + 1, logk + 1]], {j, 0, Min[m, dA]}] -
	          If[dD >= 1,
	            Sum[dCoeffs[[i + 1]] (m - i + 1) fCoeffs[[m - i + 2, logk + 1]], {i, 1, Min[m, dD]}],
	            zeroVec
	          ] -
	          logDerivativeCoeff
	        );
	        ,
	        {logk, bMaxLogK, 0, -1}
	      ];
	      , {m, 0, maxOrd - 1}
	    ];

	    fParticular = Table[
	      Sum[
	        DiffExp`Symbols`Logx^logk *
	          SeriesData[DiffExp`Symbols`x, 0,
	            DiffExp`Utilities`PChop[fCoeffs[[All, logk + 1, k]]],
	            0, maxOrd + 1, 1],
	        {logk, 0, bMaxLogK}
	      ],
	      {k, systemSize}
	    ];
  ];

  (* General solution: f = particular + sum_i c_i * F_i *)
  cIndices = Table[Subscript[c, i], {i, systemSize}];
  fGeneral = fParticular + Sum[cIndices[[i]] * FMat[[All, i]], {i, systemSize}];
  fGeneral = DiffExp`SeriesOps`SExpand[fGeneral];

  {cIndices, fGeneral, cache}
];

(* ============================================================================ *)
(* SINGULAR RECURRENCE METHOD (regular singular point / simple pole)           *)
(* For f'(x) = A(x)f + B(x) where A(x) = A_{-1}/x + A_0 + A_1*x + ...      *)
(* Requires: diagonalizable A_{-1}, non-resonant eigenvalues                  *)
(* ============================================================================ *)

(* Try to rationalize the singular matrix A(x) = R(x)/(x*D(x))
   and extract polynomial coefficients in the eigenbasis.
   Returns {True, rHatCoeffs, dCoeffs, dR, dD, d0} on success, or {False} on failure. *)
TryRationalizeSingularMatrix[ctx_Association, PInv_, P_] := Module[
  {result, Dpoly, NAMat, dR, dD, dCoeffs, rCoeffs, rHatCoeffs, d0},

  result = RationalizeAMatrixCore[ctx, True];
  If[!result[[1]], Return[{False}]];

  {Dpoly, NAMat, dR} = result[[{3, 4, 5}]];

  (* Extract D(x) polynomial coefficients *)
  dD = Exponent[Dpoly, DiffExp`Symbols`x];
  dCoeffs = Table[
    N[Coefficient[Dpoly, DiffExp`Symbols`x, i], ctx["WorkingPrecision"]],
    {i, 0, dD}
  ];
  d0 = dCoeffs[[1]];

  (* Extract R(x) = N_A(x) polynomial matrix coefficients *)
  rCoeffs = Table[
    Map[N[Coefficient[#, DiffExp`Symbols`x, j], ctx["WorkingPrecision"]] &, NAMat, {2}],
    {j, 0, dR}
  ];

  (* Transform to eigenbasis: Rhat_j = PInv . R_j . P *)
  rHatCoeffs = Table[PInv . rCoeffs[[j]] . P, {j, Length[rCoeffs]}];

  {True, rHatCoeffs, dCoeffs, dR, dD, d0}
];

(* Extract series coefficients A_0, A_1, ... from the expanded matrix (excluding the pole)
   and transform to eigenbasis. Returns {bHatCoeffs, numCoeffs}. *)
ExtractSingularSeriesCoefficients[ctx_Association, PInv_, P_] := Module[
  {maxOrd, aCoeffs, bHatCoeffs},

  maxOrd = ctx["ExpansionOrder"];

  (* Extract A_k for k = 0, 1, ..., maxOrd-1 (the regular part) *)
  aCoeffs = ExtractAMatCoefficients[ctx, 0, maxOrd - 1];

  (* Transform to eigenbasis: Bhat_k = PInv . A_k . P *)
  bHatCoeffs = Table[PInv . aCoeffs[[k]] . P, {k, Length[aCoeffs]}];

  {bHatCoeffs, maxOrd - 1}
];

(* Compute the fundamental matrix for the singular recurrence.
   Each column i corresponds to the Frobenius solution with exponent lambda_i.
   Mode "rational": uses denominator-cleared polynomial recurrence.
   Mode "series": uses direct series coefficient recurrence. *)
ComputeSingularFundamentalMatrix[ctx_Association, eigenvalues_, P_, PInv_,
    rHatCoeffs_, dCoeffs_, dR_, dD_, d0_, bHatCoeffs_, numBCoeffs_, mode_] := Module[
  {systemSize, maxOrd, gCoeffs, FMat, rhs, divisor, col, m, j},

  systemSize = ctx["SystemSize"];
  maxOrd = ctx["ExpansionOrder"];

  FMat = Table[
    (* For eigenvalue lambda_col, initial condition g_0 = e_col *)
    gCoeffs = ConstantArray[
      ConstantArray[N[0, ctx["WorkingPrecision"]], systemSize],
      maxOrd + 1
    ];
    gCoeffs[[1]] = N[UnitVector[systemSize, col], ctx["WorkingPrecision"]];

    If[mode === "rational",
      (* Denominator-clearing recurrence:
         d0*(m+lambda_col - lambda_j) * g_m^(j) =
           sum_{k=1}^{min(m,dR)} Rhat_k . g_{m-k}
           - sum_{i=1}^{min(m,dD)} d_i*(m-i+lambda_col) * g_{m-i}  *)
      Do[
        rhs = Sum[rHatCoeffs[[k + 1]] . gCoeffs[[m - k + 1]], {k, 1, Min[m, dR]}] -
          If[dD >= 1,
            Sum[dCoeffs[[i + 1]] * (m - i + eigenvalues[[col]]) * gCoeffs[[m - i + 1]],
              {i, 1, Min[m, dD]}],
            ConstantArray[0, systemSize]
          ];
        (* Component-wise division by d0*(m + lambda_col - lambda_j) *)
        gCoeffs[[m + 1]] = Table[
          divisor = d0 * (m + eigenvalues[[col]] - eigenvalues[[j]]);
          If[divisor == 0,
            (* Should not happen due to non-resonance check, but be safe *)
            N[0, ctx["WorkingPrecision"]],
            rhs[[j]] / divisor
          ],
          {j, systemSize}
        ];
        , {m, 1, maxOrd}
      ];
      ,
      (* Series mode recurrence:
         (m + lambda_col - lambda_j) * g_m^(j) =
           [sum_{k=0}^{m-1} Bhat_k . g_{m-1-k}]^(j)  *)
      Do[
        rhs = Sum[bHatCoeffs[[k + 1]] . gCoeffs[[m - 1 - k + 1]], {k, 0, Min[m - 1, numBCoeffs]}];
        gCoeffs[[m + 1]] = Table[
          divisor = m + eigenvalues[[col]] - eigenvalues[[j]];
          If[divisor == 0,
            N[0, ctx["WorkingPrecision"]],
            rhs[[j]] / divisor
          ],
          {j, systemSize}
        ];
        , {m, 1, maxOrd}
      ];
    ];

    (* Transform back to original basis and build SeriesData.
       Column = P . (x^{lambda_col} * series in eigenbasis) *)
    Table[
      (* The k-th component in original basis is sum_j P[[k,j]] * g_coeffs_j *)
      SeriesData[DiffExp`Symbols`x, 0,
        DiffExp`Utilities`PChop[Table[Sum[P[[k, j]] * gCoeffs[[n + 1, j]], {j, systemSize}], {n, 0, maxOrd}]],
        0, maxOrd + 1, 1] * DiffExp`Symbols`x^eigenvalues[[col]],
      {k, systemSize}
    ]
    , {col, systemSize}
  ] // Transpose; (* FMat[[component, solution]] *)

  DiffExp`Utilities`PChop[FMat]
];

(* Compute the particular solution for the singular recurrence.
   The particular solution starts at power s determined by the leading power of bVec.
   Returns the particular solution vector, or $Failed if resonance is detected. *)
ComputeSingularParticular[bVec_, eigenvalues_, P_, PInv_,
    rHatCoeffs_, dCoeffs_, dR_, dD_, d0_, bHatCoeffs_, numBCoeffs_, mode_, ctx_Association] := Module[
  {systemSize, maxOrd, bLeadPow, bLeadDen, s, gCoeffs, rhs, divisor, nbCoeff,
   bVecInEigen, bSeriesCoeffs, numBVecCoeffs, fParticular, m, j},

  systemSize = Length[bVec];
  maxOrd = ctx["ExpansionOrder"];

  (* Determine leading power of bVec *)
  bLeadPow = Min[Table[
    If[Head[bVec[[k]]] === SeriesData,
      bVec[[k]][[4]] / bVec[[k]][[6]],
      0
    ],
    {k, systemSize}
  ]];

  s = bLeadPow + 1;

  (* Check non-resonance for particular solution:
     need (m + s - lambda_j) != 0 for all m >= 0 and all j *)
  If[AnyTrue[eigenvalues, (IntegerQ[# - s] && # - s >= 0) &],
    Return[$Failed]
  ];

  (* Extract bVec coefficients in eigenbasis. *)
  bSeriesCoeffs = Table[
    Table[
      SeriesCoefficient[bVec[[k]], {DiffExp`Symbols`x, 0, s - 1 + n}],
      {k, systemSize}
    ],
    {n, 0, maxOrd}
  ];
  (* Transform to eigenbasis *)
  bVecInEigen = Table[PInv . bSeriesCoeffs[[n]], {n, Length[bSeriesCoeffs]}];

  (* Run recurrence for particular solution *)
  gCoeffs = ConstantArray[
    ConstantArray[N[0, ctx["WorkingPrecision"]], systemSize],
    maxOrd + 1
  ];

  If[mode === "rational",
    Do[
      nbCoeff = Sum[
        If[m - i >= 0 && m - i < Length[bVecInEigen],
          dCoeffs[[i + 1]] * bVecInEigen[[m - i + 1]],
          ConstantArray[0, systemSize]
        ],
        {i, 0, Min[m, dD]}
      ];

      rhs = nbCoeff +
        Sum[rHatCoeffs[[k + 1]] . gCoeffs[[m - k + 1]], {k, 1, Min[m, dR]}] -
        If[dD >= 1,
          Sum[dCoeffs[[i + 1]] * (m - i + s) * gCoeffs[[m - i + 1]], {i, 1, Min[m, dD]}],
          ConstantArray[0, systemSize]
        ];

      gCoeffs[[m + 1]] = Table[
        divisor = d0 * (m + s - eigenvalues[[j]]);
        rhs[[j]] / divisor,
        {j, systemSize}
      ];
      , {m, 0, maxOrd}
    ];
    ,
    Do[
      rhs = bVecInEigen[[m + 1]] +
        If[m >= 1,
          Sum[bHatCoeffs[[k + 1]] . gCoeffs[[m - 1 - k + 1]], {k, 0, Min[m - 1, numBCoeffs]}],
          ConstantArray[0, systemSize]
        ];

      gCoeffs[[m + 1]] = Table[
        divisor = m + s - eigenvalues[[j]];
        rhs[[j]] / divisor,
        {j, systemSize}
      ];
      , {m, 0, maxOrd}
    ];
  ];

  (* Transform back to original basis and build SeriesData at power s *)
  fParticular = Table[
    SeriesData[DiffExp`Symbols`x, 0,
      DiffExp`Utilities`PChop[Table[Sum[P[[k, j]] * gCoeffs[[n + 1, j]], {j, systemSize}], {n, 0, maxOrd}]],
      0, maxOrd + 1, 1] * DiffExp`Symbols`x^s,
    {k, systemSize}
  ];

  DiffExp`Utilities`PChop[fParticular]
];

(* Main entry point for the singular recurrence solver.
   Handles regular singular points with diagonalizable, non-resonant residue. *)
SolveSingularRecurrence[ctx_Association, bVec_, epsord_, cacheIn_Association] := Module[
  {systemSize, maxOrd, FMat, fParticular, fGeneral, cIndices, c,
   eigenvalues, P, PInv,
   rHatCoeffs, dCoeffs, dR, dD, d0, bHatCoeffs, numBCoeffs,
   mode, rationalResult,
   cache = cacheIn},

  systemSize = ctx["SystemSize"];
  maxOrd = ctx["ExpansionOrder"];

  If[epsord === 0 && !KeyExistsQ[cache, "SingRR"],
    (* Get eigenvalue data from cache (stored by dispatch) or compute fresh *)
    If[KeyExistsQ[cache, "SingularEigenData"],
      {eigenvalues, P, PInv} = cache["SingularEigenData"];
      ,
      {eigenvalues, P, PInv} = PrepareSingularRecurrence[ctx];
    ];
    DiffExp`Utilities`PrintInfo["Using singular recurrence for integrals ", ctx["Label"],
      " (eigenvalues: ", eigenvalues, ")."][3];

    (* Try denominator-clearing mode first *)
    rationalResult = TryRationalizeSingularMatrix[ctx, PInv, P];

    If[rationalResult[[1]],
      mode = "rational";
      {rHatCoeffs, dCoeffs, dR, dD, d0} = rationalResult[[{2, 3, 4, 5, 6}]];
      bHatCoeffs = {};
      numBCoeffs = 0;
      DiffExp`Utilities`PrintInfo["  Singular recurrence: denominator-clearing mode (poly degree ", dR, ")."][3];
      ,
      mode = "series";
      {bHatCoeffs, numBCoeffs} = ExtractSingularSeriesCoefficients[ctx, PInv, P];
      rHatCoeffs = {};
      dCoeffs = {};
      dR = 0;
      dD = 0;
      d0 = N[1, ctx["WorkingPrecision"]];
      DiffExp`Utilities`PrintInfo["  Singular recurrence: series mode."][3];
    ];

    (* Compute fundamental matrix *)
    FMat = ComputeSingularFundamentalMatrix[ctx, eigenvalues, P, PInv,
      rHatCoeffs, dCoeffs, dR, dD, d0, bHatCoeffs, numBCoeffs, mode];

    cache["SingRR"] = {FMat, eigenvalues, P, PInv,
      rHatCoeffs, dCoeffs, dR, dD, d0, bHatCoeffs, numBCoeffs, mode};
  ];

  {FMat, eigenvalues, P, PInv,
   rHatCoeffs, dCoeffs, dR, dD, d0, bHatCoeffs, numBCoeffs, mode} = cache["SingRR"];

  (* Compute particular solution *)
  If[DiffExp`Utilities`PChop[bVec] === ConstantArray[0, systemSize],
    fParticular = ConstantArray[
      SeriesData[DiffExp`Symbols`x, 0, {}, 0, maxOrd + 1, 1],
      systemSize
    ];
    ,
    fParticular = ComputeSingularParticular[bVec, eigenvalues, P, PInv,
      rHatCoeffs, dCoeffs, dR, dD, d0, bHatCoeffs, numBCoeffs, mode, ctx];

    (* If this restricted solver hits a resonance, stay on the recursive path. *)
    If[fParticular === $Failed,
      DiffExp`Utilities`PrintInfo[
        "Singular recurrence: particular solution resonance detected; switching to general singular recurrence."
      ][3];
      Return[SolveGeneralSingularRecurrence[ctx, bVec, epsord, cache]]
    ];
  ];

  (* General solution: f = particular + sum_i c_i * FMat_column_i *)
  cIndices = Table[Subscript[c, i], {i, systemSize}];
  fGeneral = fParticular + Sum[cIndices[[i]] * FMat[[All, i]], {i, systemSize}];
  fGeneral = DiffExp`SeriesOps`SExpand[fGeneral];

  {cIndices, fGeneral, cache}
];
