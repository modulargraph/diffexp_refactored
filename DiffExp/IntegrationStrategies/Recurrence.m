(* IntegrationStrategies/Recurrence.m *)
(* Rational recurrence and singular recurrence strategies *)

(* ============================================================================ *)
(* RATIONAL RECURRENCE METHOD                                                  *)
(* ============================================================================ *)

(* Check if the rational recurrence method is applicable:
   The expansion point must be non-singular (series starts at x^0).
   Checks the expanded matrix series for the given integral block. *)
RationalRecurrenceApplicableQ[intind_, line_] := Quiet[Check[
  Module[{AMat, atZero},

    (* Check that factored matrices exist *)
    If[!KeyExistsQ[DiffExp`State`DEqnMatricesFactored, line],
      Return[False]
    ];

    (* Non-singular iff the factored matrix is finite at x=0 *)
    AMat = DiffExp`State`DEqnMatricesFactored[line][0][[intind, intind]];
    atZero = AMat /. DiffExp`Symbols`x -> 0;

    TrueQ[And @@ (NumericQ /@ Flatten[atZero])]
  ],
  False (* Return False on any error *)
]];

(* Try to rationalize the factored A matrix and extract polynomial coefficients.
   Returns {True, dCoeffs, aCoeffs, dD, dA, d0, d0Inv} on success,
   or {False} if rationalization is too expensive or not possible. *)
TryRationalizeMatrix[intind_, line_] := Module[
  {AMat, AMatTogether, flatEntries, denoms, Dpoly, NAMat,
   dCoeffs, aCoeffs, dD, dA, d0, d0Inv, maxOrd, result},

  maxOrd = DiffExp`State`ExpansionOrderVal;

  (* Only attempt rationalization if factored matrices exist *)
  If[!KeyExistsQ[DiffExp`State`DEqnMatricesFactored, line],
    Return[{False}]
  ];

  (* Get the factored A matrix *)
  AMat = DiffExp`State`DEqnMatricesFactored[line][0][[intind, intind]];

  (* Try to put in Together form - use TimeConstrained to avoid hanging *)
  result = TimeConstrained[
    Module[{at, fe, dn, dp},
      at = Map[Together, AMat, {2}];
      fe = Flatten[at];

      (* Check rationality *)
      If[!And @@ (PolynomialQ[Numerator[#], DiffExp`Symbols`x] &&
                  PolynomialQ[Denominator[#], DiffExp`Symbols`x] & /@ fe),
        Return[{False}]
      ];

      (* Compute common denominator *)
      dn = DeleteDuplicates[Denominator /@ fe];
      dn = DeleteCases[dn, a_ /; FreeQ[a, DiffExp`Symbols`x]];
      If[dn === {},
        dp = 1;,
        dp = PolynomialLCM @@ dn;
      ];

      (* Check non-singularity *)
      If[(dp /. DiffExp`Symbols`x -> 0) === 0, Return[{False}]];

      (* Check that polynomial degree is beneficial: dA < maxOrd/2 *)
      NAMat = Expand[dp * at];
      dA = Max[0, Max[Exponent[#, DiffExp`Symbols`x] & /@ Flatten[NAMat]]];
      If[dA > maxOrd/2, Return[{False}]]; (* Not beneficial over direct series *)

      {True, at, dp, NAMat, dA}
    ],
    5.0, (* 5 second timeout *)
    {False}
  ];

  If[!result[[1]], Return[{False}]];

  {AMatTogether, Dpoly, NAMat, dA} = result[[{2, 3, 4, 5}]];

  (* Extract polynomial coefficients *)
  dD = Exponent[Dpoly, DiffExp`Symbols`x];
  dCoeffs = Table[
    N[Coefficient[Dpoly, DiffExp`Symbols`x, i], DiffExp`State`FEWorkingPrecision],
    {i, 0, dD}
  ];

  aCoeffs = Table[
    Map[N[Coefficient[#, DiffExp`Symbols`x, j], DiffExp`State`FEWorkingPrecision] &, NAMat, {2}],
    {j, 0, dA}
  ];

  d0 = dCoeffs[[1]];
  d0Inv = 1/d0;

  {True, dCoeffs, aCoeffs, dD, dA, d0, d0Inv}
];

(* Extract series coefficients of the A matrix from the expanded form.
   Returns {aCoeffs, dA} where dA = maxOrd (all orders used). *)
ExtractSeriesCoefficients[intind_, line_] := Module[
  {AMatExpanded, maxOrd, aCoeffs, systemSize},

  systemSize = Length[intind];
  maxOrd = DiffExp`State`ExpansionOrderVal;
  AMatExpanded = DiffExp`State`DEqnMatricesExpanded[line][0][[intind, intind]];

  (* Extract the matrix coefficient at each order from the SeriesData *)
  aCoeffs = Table[
    Table[
      SeriesCoefficient[AMatExpanded[[i, k]], {DiffExp`Symbols`x, 0, j}],
      {i, systemSize}, {k, systemSize}
    ],
    {j, 0, maxOrd - 1}
  ];

  {aCoeffs, maxOrd - 1}
];

(* Core recurrence computation: computes the fundamental matrix (homogeneous solutions)
   and stores preprocessed data in BufferedData. *)
ComputeFundamentalMatrix[intind_, line_, aCoeffs_, dCoeffs_, dD_, dA_, d0Inv_, BufferedData_] := Module[
  {systemSize, maxOrd, fCoeffs, FMat},

  systemSize = Length[intind];
  maxOrd = DiffExp`State`ExpansionOrderVal;

  FMat = Table[
    fCoeffs = ConstantArray[
      ConstantArray[N[0, DiffExp`State`FEWorkingPrecision], systemSize],
      maxOrd + 1
    ];
    fCoeffs[[1]] = N[UnitVector[systemSize, col], DiffExp`State`FEWorkingPrecision];

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
SolveRationalRecurrence[intind_, bVec_, line_, epsord_, BufferedDataIn_] := Module[
  {systemSize, maxOrd, fCoeffs, bSeriesCoeffs, bCoeffAtN, nbCoeff,
   FMat, fParticular, fGeneral, cIndices, c,
   dCoeffs, aCoeffs, dD, dA, d0, d0Inv, rationalResult,
   BufferedData = BufferedDataIn},

  systemSize = Length[intind];
  maxOrd = DiffExp`State`ExpansionOrderVal;

  If[epsord === 0 && !KeyExistsQ[BufferedData, {"RR", intind}],
    DiffExp`State`BenchmarkData["Segments"][line // N]["HomogeneousSolveAllPreprocessing"]["Integrals"][intind] = AbsoluteTime[];

    (* Try the denominator-clearing approach first *)
    rationalResult = TryRationalizeMatrix[intind, line];

    If[rationalResult[[1]],
      (* Denominator-clearing mode: polynomial recurrence *)
      {dCoeffs, aCoeffs, dD, dA, d0, d0Inv} = rationalResult[[{2, 3, 4, 5, 6, 7}]];
      DiffExp`Utilities`PrintInfo["Using rational recurrence (poly degree ", dA, ") for integrals ", intind, "."][3];
      ,
      (* Direct series mode: use expanded series coefficients *)
      {aCoeffs, dA} = ExtractSeriesCoefficients[intind, line];
      dD = 0;
      dCoeffs = {N[1, DiffExp`State`FEWorkingPrecision]};
      d0Inv = N[1, DiffExp`State`FEWorkingPrecision];
      DiffExp`Utilities`PrintInfo["Using series recurrence for integrals ", intind, "."][3];
    ];

    (* Compute fundamental matrix *)
    FMat = ComputeFundamentalMatrix[intind, line, aCoeffs, dCoeffs, dD, dA, d0Inv, BufferedData];

    BufferedData[{"RR", intind}] = {FMat, dCoeffs, aCoeffs, dD, dA, d0Inv};

    DiffExp`State`BenchmarkData["Segments"][line // N]["HomogeneousSolveAllPreprocessing"]["Integrals"][intind] = AbsoluteTime[] - DiffExp`State`BenchmarkData["Segments"][line // N]["HomogeneousSolveAllPreprocessing"]["Integrals"][intind];
  ];

  {FMat, dCoeffs, aCoeffs, dD, dA, d0Inv} = BufferedData[{"RR", intind}];

  (* Compute particular solution using bVec *)
  If[DiffExp`Utilities`PChop[bVec] === ConstantArray[0, systemSize],
    fParticular = ConstantArray[
      SeriesData[DiffExp`Symbols`x, 0, {}, 0, maxOrd + 1, 1],
      systemSize
    ];
    ,
    (* Extract series coefficients from bVec *)
    bSeriesCoeffs = Table[
      bCoeffAtN = SeriesCoefficient[bVec[[k]], {DiffExp`Symbols`x, 0, n}];
      If[NumericQ[bCoeffAtN], bCoeffAtN, DiffExp`Utilities`PChop[bCoeffAtN]],
      {n, 0, maxOrd - 1}, {k, systemSize}
    ];

    (* Run recurrence with f_0 = 0 for particular solution *)
    fCoeffs = ConstantArray[
      ConstantArray[N[0, DiffExp`State`FEWorkingPrecision], systemSize],
      maxOrd + 1
    ];

    Do[
      (* Compute N_B coefficient at order m: nb_m = sum_{i=0}^{min(m,dD)} d_i * B_{m-i} *)
      nbCoeff = Sum[
        If[m - i <= maxOrd - 1,
          dCoeffs[[i + 1]] * bSeriesCoeffs[[m - i + 1]],
          ConstantArray[0, systemSize]
        ],
        {i, 0, Min[m, dD]}
      ];

      (* f_{m+1} = d0Inv/(m+1) * [nb_m + sum_j a_j.f_{m-j} - sum_{i>=1} d_i*(m-i+1)*f_{m-i+1}] *)
      fCoeffs[[m + 2]] = d0Inv/(m + 1) * (
        nbCoeff +
        Sum[aCoeffs[[j + 1]] . fCoeffs[[m - j + 1]], {j, 0, Min[m, dA]}] -
        If[dD >= 1,
          Sum[dCoeffs[[i + 1]] (m - i + 1) fCoeffs[[m - i + 2]], {i, 1, Min[m, dD]}],
          0
        ]
      );
      , {m, 0, maxOrd - 1}
    ];

    fParticular = Table[
      SeriesData[DiffExp`Symbols`x, 0,
        DiffExp`Utilities`PChop[fCoeffs[[All, k]]],
        0, maxOrd + 1, 1],
      {k, systemSize}
    ];
  ];

  (* General solution: f = particular + sum_i c_i * F_i *)
  cIndices = Table[Subscript[c, i], {i, systemSize}];
  fGeneral = fParticular + Sum[cIndices[[i]] * FMat[[All, i]], {i, systemSize}];
  fGeneral = DiffExp`SeriesOps`SExpand[fGeneral];

  {cIndices, fGeneral, BufferedData}
];

(* ============================================================================ *)
(* SINGULAR RECURRENCE METHOD (regular singular point / simple pole)           *)
(* For f'(x) = A(x)f + B(x) where A(x) = A_{-1}/x + A_0 + A_1*x + ...      *)
(* Requires: diagonalizable A_{-1}, non-resonant eigenvalues                  *)
(* ============================================================================ *)

(* Check if the singular recurrence method is applicable:
   1. The expanded matrix has minimum order -1 (simple pole at x=0)
   2. The residue matrix A_{-1} is diagonalizable
   3. No two eigenvalues differ by a positive integer (non-resonance) *)
SingularRecurrenceApplicableQ[intind_, line_] := Quiet[Check[
  Module[{AMatExpanded, minOrders, minOrder, residueMat, systemSize,
          eigenvalues, eigenvectors, diffs},

    systemSize = Length[intind];

    (* Check that expanded matrices exist *)
    If[!KeyExistsQ[DiffExp`State`DEqnMatricesExpanded, line],
      Return[False]
    ];

    AMatExpanded = DiffExp`State`DEqnMatricesExpanded[line][0][[intind, intind]];

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
    If[minOrder =!= -1, Return[False]];

    (* Extract residue matrix A_{-1} *)
    residueMat = Table[
      SeriesCoefficient[AMatExpanded[[i, j]], {DiffExp`Symbols`x, 0, -1}],
      {i, systemSize}, {j, systemSize}
    ];

    (* Residue must not be identically zero *)
    If[DiffExp`Utilities`PChop[residueMat] === ConstantArray[0, {systemSize, systemSize}],
      Return[False]
    ];

    (* Compute eigenvalues *)
    eigenvalues = Eigenvalues[N[residueMat, DiffExp`State`FEWorkingPrecision]];
    eigenvalues = Rationalize[eigenvalues, 10^(-DiffExp`State`ChopPrecisionVal/2)];

    (* Check diagonalizability: rank of eigenvector matrix must equal system size *)
    eigenvectors = Eigenvectors[N[residueMat, DiffExp`State`FEWorkingPrecision]];
    If[MatrixRank[eigenvectors, Tolerance -> 10^(-DiffExp`State`ChopPrecisionVal/2)] < systemSize,
      Return[False]
    ];

    (* Check non-resonance: no difference lambda_j - lambda_i is a positive integer *)
    diffs = Flatten[Table[
      eigenvalues[[j]] - eigenvalues[[i]],
      {i, systemSize}, {j, systemSize}
    ]];
    If[AnyTrue[diffs, (IntegerQ[#] && # > 0) &],
      Return[False]
    ];

    True
  ],
  False (* Return False on any error *)
]];

(* Diagonalize the residue matrix.
   Returns {eigenvalues, P, PInv} where P is the matrix of eigenvectors (columns)
   and PInv = P^{-1}. Eigenvalues are rationalized for exact arithmetic. *)
DiagonalizeResidue[residueMat_, systemSize_] := Module[
  {eigenvalues, eigenvectors, P, PInv, sorted, perm},

  {eigenvalues, eigenvectors} = Eigensystem[N[residueMat, DiffExp`State`FEWorkingPrecision]];

  (* Rationalize eigenvalues for exact recurrence denominators *)
  eigenvalues = Rationalize[eigenvalues, 10^(-DiffExp`State`ChopPrecisionVal/2)];

  (* P = matrix of eigenvectors as columns, so P = Transpose[eigenvectors] *)
  P = Transpose[eigenvectors];
  PInv = Inverse[P];

  (* Normalize to working precision *)
  P = N[P, DiffExp`State`FEWorkingPrecision];
  PInv = N[PInv, DiffExp`State`FEWorkingPrecision];

  {eigenvalues, P, PInv}
];

(* Try to rationalize the singular matrix A(x) = R(x)/(x*D(x))
   and extract polynomial coefficients in the eigenbasis.
   Returns {True, rHatCoeffs, dCoeffs, dR, dD, d0} on success, or {False} on failure. *)
TryRationalizeSingularMatrix[intind_, line_, PInv_, P_] := Module[
  {AMat, systemSize, result, AMatTogether, flatEntries, denoms, Dpoly, xDpoly,
   NAMat, dD, dR, dCoeffs, rCoeffs, rHatCoeffs, d0, maxOrd},

  systemSize = Length[intind];
  maxOrd = DiffExp`State`ExpansionOrderVal;

  (* Only attempt if factored matrices exist *)
  If[!KeyExistsQ[DiffExp`State`DEqnMatricesFactored, line],
    Return[{False}]
  ];

  AMat = DiffExp`State`DEqnMatricesFactored[line][0][[intind, intind]];

  result = TimeConstrained[
    Module[{at, fe, dn, dp, xdp, naMat, degA, degD},
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
      If[dn === {},
        Return[{False}] (* No x-dependence in denominator means no pole *)
      ];
      dp = PolynomialLCM @@ dn;

      (* The denominator must have exactly one factor of x (simple pole) *)
      If[Exponent[dp, DiffExp`Symbols`x, Min] < 1,
        (* dp doesn't have x as factor - not a pole *)
        Return[{False}]
      ];

      (* Factor out x: dp = x * D(x) *)
      xdp = dp;
      dp = Cancel[dp / DiffExp`Symbols`x];

      (* D(0) must be nonzero *)
      If[(dp /. DiffExp`Symbols`x -> 0) === 0, Return[{False}]];

      (* N_A = x*D(x) * A(x) = polynomial matrix *)
      naMat = Expand[xdp * at];

      (* Check polynomial degrees *)
      degA = Max[0, Max[Exponent[#, DiffExp`Symbols`x] & /@ Flatten[naMat]]];
      degD = Exponent[dp, DiffExp`Symbols`x];

      (* Only beneficial if polynomial degree is small *)
      If[degA > maxOrd/2, Return[{False}]];

      {True, at, dp, naMat, degA, degD}
    ],
    5.0, (* 5 second timeout *)
    {False}
  ];

  If[!result[[1]], Return[{False}]];

  {AMatTogether, Dpoly, NAMat, dR, dD} = result[[{2, 3, 4, 5, 6}]];

  (* Extract D(x) polynomial coefficients *)
  dCoeffs = Table[
    N[Coefficient[Dpoly, DiffExp`Symbols`x, i], DiffExp`State`FEWorkingPrecision],
    {i, 0, dD}
  ];
  d0 = dCoeffs[[1]];

  (* Extract R(x) = N_A(x) polynomial matrix coefficients *)
  rCoeffs = Table[
    Map[N[Coefficient[#, DiffExp`Symbols`x, j], DiffExp`State`FEWorkingPrecision] &, NAMat, {2}],
    {j, 0, dR}
  ];

  (* Transform to eigenbasis: Rhat_j = PInv . R_j . P *)
  rHatCoeffs = Table[PInv . rCoeffs[[j]] . P, {j, Length[rCoeffs]}];

  {True, rHatCoeffs, dCoeffs, dR, dD, d0}
];

(* Extract series coefficients A_0, A_1, ... from the expanded matrix (excluding the pole)
   and transform to eigenbasis. Returns {bHatCoeffs, numCoeffs}. *)
ExtractSingularSeriesCoefficients[intind_, line_, PInv_, P_] := Module[
  {AMatExpanded, maxOrd, systemSize, aCoeffs, bHatCoeffs},

  systemSize = Length[intind];
  maxOrd = DiffExp`State`ExpansionOrderVal;
  AMatExpanded = DiffExp`State`DEqnMatricesExpanded[line][0][[intind, intind]];

  (* Extract A_k for k = 0, 1, ..., maxOrd-1 (the regular part) *)
  aCoeffs = Table[
    Table[
      SeriesCoefficient[AMatExpanded[[i, k]], {DiffExp`Symbols`x, 0, j}],
      {i, systemSize}, {k, systemSize}
    ],
    {j, 0, maxOrd - 1}
  ];

  (* Transform to eigenbasis: Bhat_k = PInv . A_k . P *)
  bHatCoeffs = Table[PInv . aCoeffs[[k]] . P, {k, Length[aCoeffs]}];

  {bHatCoeffs, maxOrd - 1}
];

(* Compute the fundamental matrix for the singular recurrence.
   Each column i corresponds to the Frobenius solution with exponent lambda_i.
   Mode "rational": uses denominator-cleared polynomial recurrence.
   Mode "series": uses direct series coefficient recurrence. *)
ComputeSingularFundamentalMatrix[intind_, eigenvalues_, P_, PInv_,
    rHatCoeffs_, dCoeffs_, dR_, dD_, d0_, bHatCoeffs_, numBCoeffs_, mode_] := Module[
  {systemSize, maxOrd, gCoeffs, FMat, rhs, divisor, col, m, j},

  systemSize = Length[intind];
  maxOrd = DiffExp`State`ExpansionOrderVal;

  FMat = Table[
    (* For eigenvalue lambda_col, initial condition g_0 = e_col *)
    gCoeffs = ConstantArray[
      ConstantArray[N[0, DiffExp`State`FEWorkingPrecision], systemSize],
      maxOrd + 1
    ];
    gCoeffs[[1]] = N[UnitVector[systemSize, col], DiffExp`State`FEWorkingPrecision];

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
            N[0, DiffExp`State`FEWorkingPrecision],
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
            N[0, DiffExp`State`FEWorkingPrecision],
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
    rHatCoeffs_, dCoeffs_, dR_, dD_, d0_, bHatCoeffs_, numBCoeffs_, mode_] := Module[
  {systemSize, maxOrd, bLeadPow, bLeadDen, s, gCoeffs, rhs, divisor, nbCoeff,
   bVecInEigen, bSeriesCoeffs, numBVecCoeffs, fParticular, m, j},

  systemSize = Length[bVec];
  maxOrd = DiffExp`State`ExpansionOrderVal;

  (* Determine leading power of bVec *)
  bLeadPow = Min[Table[
    If[Head[bVec[[k]]] === SeriesData,
      bVec[[k]][[4]] / bVec[[k]][[6]],
      0
    ],
    {k, systemSize}
  ]];

  (* For the ODE x*f' = A_{-1}*f + (regular)*f + x*B,
     the particular solution starts at power s = bLeadPow + 1
     (since x*B shifts B up by one power, and the recurrence preserves order) *)
  (* Actually, f' = A(x)*f + B means:
     x*f' = A_{-1}*f + x*A_0*f + ... + x*B
     If B starts at x^bLeadPow, then x*B starts at x^{bLeadPow+1}.
     The particular solution starts at x^s where s = bLeadPow + 1 is the first
     power where the inhomogeneous term contributes to the recurrence. *)
  (* However, if bLeadPow < min(eigenvalues), f_p could start at bLeadPow itself
     from the derivative balance. The correct s comes from the full ODE:
     f' = (A_{-1}/x + ...)*f + B
     At leading order: s*f_s*x^{s-1} = A_{-1}*f_s*x^{s-1} + B_bLeadPow*x^{bLeadPow}
     If s-1 = bLeadPow: (s*I - A_{-1})*f_s = B_bLeadPow *)
  s = bLeadPow + 1;

  (* Check non-resonance for particular solution:
     need (m + s - lambda_j) != 0 for all m >= 0 and all j *)
  If[AnyTrue[eigenvalues, (IntegerQ[# - s] && # - s >= 0) &],
    Return[$Failed]
  ];

  (* Extract bVec coefficients in eigenbasis.
     Transform: b_eigen = PInv . bVec componentwise *)
  (* We need coefficients of x*B (= xB) at powers s, s+1, ..., s+maxOrd *)
  (* (xB)_{s+m} = B_{s+m-1} = B_{bLeadPow+m} *)
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
    ConstantArray[N[0, DiffExp`State`FEWorkingPrecision], systemSize],
    maxOrd + 1
  ];

  If[mode === "rational",
    (* Denominator-clearing:
       d0*(m+s-lambda_j)*g_m^(j) =
         sum_{k=1}^{min(m,dR)} Rhat_k . g_{m-k}
         - sum_{i=1}^{min(m,dD)} d_i*(m-i+s)*g_{m-i}
         + [D(x)*x*B in eigenbasis]_{m+s}  *)
    (* For the inhomogeneous term with denominator clearing:
       x*D(x)*B has coefficients that combine d_i and B coefficients.
       nb_m = sum_{i=0}^{min(m,dD)} d_i * bVecInEigen_{m-i} *)
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
    (* Series mode:
       (m+s-lambda_j)*g_m^(j) =
         [sum_{k=0}^{m-1} Bhat_k . g_{m-1-k}]^(j) + bVecInEigen_m  *)
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
SolveSingularRecurrence[intind_, bVec_, line_, epsord_, BufferedDataIn_] := Module[
  {systemSize, maxOrd, FMat, fParticular, fGeneral, cIndices, c,
   eigenvalues, P, PInv, residueMat,
   rHatCoeffs, dCoeffs, dR, dD, d0, bHatCoeffs, numBCoeffs,
   mode, rationalResult, AMatExpanded,
   BufferedData = BufferedDataIn},

  systemSize = Length[intind];
  maxOrd = DiffExp`State`ExpansionOrderVal;

  If[epsord === 0 && !KeyExistsQ[BufferedData, {"SingRR", intind}],
    DiffExp`State`BenchmarkData["Segments"][line // N]["HomogeneousSolveAllPreprocessing"]["Integrals"][intind] = AbsoluteTime[];

    (* Extract residue matrix A_{-1} *)
    AMatExpanded = DiffExp`State`DEqnMatricesExpanded[line][0][[intind, intind]];
    residueMat = Table[
      SeriesCoefficient[AMatExpanded[[i, j]], {DiffExp`Symbols`x, 0, -1}],
      {i, systemSize}, {j, systemSize}
    ];

    (* Diagonalize the residue *)
    {eigenvalues, P, PInv} = DiagonalizeResidue[residueMat, systemSize];
    DiffExp`Utilities`PrintInfo["Using singular recurrence for integrals ", intind,
      " (eigenvalues: ", eigenvalues, ")."][3];

    (* Try denominator-clearing mode first *)
    rationalResult = TryRationalizeSingularMatrix[intind, line, PInv, P];

    If[rationalResult[[1]],
      mode = "rational";
      {rHatCoeffs, dCoeffs, dR, dD, d0} = rationalResult[[{2, 3, 4, 5, 6}]];
      bHatCoeffs = {};
      numBCoeffs = 0;
      DiffExp`Utilities`PrintInfo["  Singular recurrence: denominator-clearing mode (poly degree ", dR, ")."][3];
      ,
      mode = "series";
      {bHatCoeffs, numBCoeffs} = ExtractSingularSeriesCoefficients[intind, line, PInv, P];
      rHatCoeffs = {};
      dCoeffs = {};
      dR = 0;
      dD = 0;
      d0 = N[1, DiffExp`State`FEWorkingPrecision];
      DiffExp`Utilities`PrintInfo["  Singular recurrence: series mode."][3];
    ];

    (* Compute fundamental matrix *)
    FMat = ComputeSingularFundamentalMatrix[intind, eigenvalues, P, PInv,
      rHatCoeffs, dCoeffs, dR, dD, d0, bHatCoeffs, numBCoeffs, mode];

    BufferedData[{"SingRR", intind}] = {FMat, eigenvalues, P, PInv,
      rHatCoeffs, dCoeffs, dR, dD, d0, bHatCoeffs, numBCoeffs, mode};

    DiffExp`State`BenchmarkData["Segments"][line // N]["HomogeneousSolveAllPreprocessing"]["Integrals"][intind] = AbsoluteTime[] - DiffExp`State`BenchmarkData["Segments"][line // N]["HomogeneousSolveAllPreprocessing"]["Integrals"][intind];
  ];

  {FMat, eigenvalues, P, PInv,
   rHatCoeffs, dCoeffs, dR, dD, d0, bHatCoeffs, numBCoeffs, mode} = BufferedData[{"SingRR", intind}];

  (* Compute particular solution *)
  If[DiffExp`Utilities`PChop[bVec] === ConstantArray[0, systemSize],
    fParticular = ConstantArray[
      SeriesData[DiffExp`Symbols`x, 0, {}, 0, maxOrd + 1, 1],
      systemSize
    ];
    ,
    fParticular = ComputeSingularParticular[bVec, eigenvalues, P, PInv,
      rHatCoeffs, dCoeffs, dR, dD, d0, bHatCoeffs, numBCoeffs, mode];

    (* If particular solution fails (resonance), fall back to Default *)
    If[fParticular === $Failed,
      DiffExp`Utilities`PrintInfo["Singular recurrence: particular solution resonance detected, falling back to Default."][3];
      (* Ensure SolveDefault's BufferedData is initialized *)
      If[!KeyExistsQ[BufferedData, intind],
        BufferedData = SolveDefault[intind,
          ConstantArray[SeriesData[DiffExp`Symbols`x, 0, {}, 0, maxOrd + 1, 1], systemSize],
          line, 0, BufferedData][[3]];
      ];
      Return[SolveDefault[intind, bVec, line, epsord, BufferedData]]
    ];
  ];

  (* General solution: f = particular + sum_i c_i * FMat_column_i *)
  cIndices = Table[Subscript[c, i], {i, systemSize}];
  fGeneral = fParticular + Sum[cIndices[[i]] * FMat[[All, i]], {i, systemSize}];
  fGeneral = DiffExp`SeriesOps`SExpand[fGeneral];

  {cIndices, fGeneral, BufferedData}
];
