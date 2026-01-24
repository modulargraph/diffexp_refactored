(* IntegrationStrategies/ResonantRecurrence.m *)
(* General singular recurrence strategy for resonant and non-diagonalizable cases *)
(* Handles: eigenvalue multiplicities, eigenvalues differing by integers, log terms *)

(* ============================================================================ *)
(* GENERAL SINGULAR RECURRENCE METHOD                                          *)
(* For f'(x) = A(x)f + B(x) where A(x) has a regular singular point (1/x pole)*)
(* Handles resonant eigenvalues (differing by integers) and multiplicities.    *)
(* The ansatz is: f(x) = x^lambda * sum_n x^n sum_k f_{n,k} (ln x)^k          *)
(* ============================================================================ *)

(* Check if the general singular recurrence method is applicable:
   1. The expanded matrix has minimum order -1 (simple pole at x=0)
   This is a more permissive check than SingularRecurrenceApplicableQ -
   it does NOT require diagonalizability or non-resonance. *)
GeneralSingularRecurrenceApplicableQ[ctx_Association] := Quiet[Check[
  Module[{AMatExpanded, minOrders, minOrder, residueMat, systemSize},

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
    If[minOrder =!= -1, Return[False]];

    (* Extract residue matrix M_0 = coefficient of 1/x in A(x), i.e., x*A(x)|_{x=0} *)
    residueMat = Table[
      SeriesCoefficient[AMatExpanded[[i, j]], {DiffExp`Symbols`x, 0, -1}],
      {i, systemSize}, {j, systemSize}
    ];

    (* Residue must not be identically zero *)
    If[DiffExp`Utilities`PChop[residueMat] === ConstantArray[0, {systemSize, systemSize}],
      Return[False]
    ];

    True
  ],
  False
]];

(* Compute the resonance structure of the eigenvalues.
   Returns an Association with:
   - "Eigenvalues": rationalized eigenvalues (list)
   - "ResonanceClasses": groups of indices where eigenvalues differ by integers
   - "MaxLogPower": for each solution index, the max log power needed
   - "SolutionExponents": for each solution, the leading exponent
   - "InitialLogPower": for each solution, the log power of the initial condition
   - "InitialVectors": for each solution, the initial f_{0, initK} vector
*)
ComputeResonanceStructure[residueMat_, systemSize_, ctx_Association] := Module[
  {eigenvalues, rawEigenvalues, jordanS, jordanJ, jordanSInv,
   blockStarts, blockSizes, blockEigenvalues,
   resonanceClasses, visited, class, classIdx,
   solutionExponents, maxLogPowers, initialLogPowers, initialVectors,
   solIdx, blockIdx, posInBlock, currentEV,
   nullVecs, genEigVecs, tolerance},

  tolerance = 10^(-ctx["ChopPrecision"]/2);

  (* Compute Jordan decomposition *)
  {jordanS, jordanJ} = JordanDecomposition[N[residueMat, ctx["WorkingPrecision"]]];
  jordanSInv = Inverse[jordanS];

  (* Rationalize eigenvalues on the diagonal of J *)
  rawEigenvalues = Diagonal[N[jordanJ, ctx["WorkingPrecision"]]];
  eigenvalues = Rationalize[rawEigenvalues, tolerance];

  (* Identify Jordan block structure from J *)
  blockStarts = {1};
  Do[
    If[i > 1 && (eigenvalues[[i]] =!= eigenvalues[[i-1]] ||
        Abs[jordanJ[[i-1, i]]] < tolerance),
      AppendTo[blockStarts, i]
    ],
    {i, 2, systemSize}
  ];
  AppendTo[blockStarts, systemSize + 1]; (* sentinel *)
  blockSizes = Differences[blockStarts];
  blockEigenvalues = eigenvalues[[blockStarts[[;; -2]]]];

  (* Group eigenvalues into resonance classes:
     two eigenvalues are in the same class if they differ by an integer *)
  visited = ConstantArray[False, Length[blockEigenvalues]];
  resonanceClasses = {};
  Do[
    If[!visited[[i]],
      class = {i};
      visited[[i]] = True;
      Do[
        If[!visited[[j]] && IntegerQ[blockEigenvalues[[i]] - blockEigenvalues[[j]]],
          AppendTo[class, j];
          visited[[j]] = True;
        ],
        {j, i + 1, Length[blockEigenvalues]}
      ];
      (* Sort class by eigenvalue (smallest first) *)
      class = SortBy[class, blockEigenvalues[[#]] &];
      AppendTo[resonanceClasses, class];
    ],
    {i, Length[blockEigenvalues]}
  ];

  (* For each solution, determine:
     - Leading exponent (eigenvalue)
     - Max log power
     - Initial condition vector *)
  solutionExponents = {};
  maxLogPowers = {};
  initialLogPowers = {};
  initialVectors = {};

  solIdx = 0;
  Do[
    Do[
      blockIdx = resonanceClasses[[classIdx, posInClass]];
      currentEV = blockEigenvalues[[blockIdx]];

      Do[
        solIdx++;

        (* For position posInBlock within a Jordan block of size blockSizes[[blockIdx]]:
           - The initial log power is posInBlock - 1
           - The max log power = (posInBlock - 1) + resonance contribution from higher eigenvalues *)

        AppendTo[solutionExponents, currentEV];

        (* Initial log power for this solution within its Jordan block *)
        AppendTo[initialLogPowers, posInBlock - 1];

        (* Max log power: (posInBlock - 1) from Jordan structure +
           sum of block sizes of eigenvalues STRICTLY greater than currentEV in the class *)
        With[{resonanceK = Total[blockSizes[[#]] & /@ Select[resonanceClasses[[classIdx]],
              blockEigenvalues[[#]] > currentEV &]]},
          AppendTo[maxLogPowers, posInBlock - 1 + resonanceK];
        ];

        (* Initial vector: always the EIGENVECTOR (first column of the Jordan block).
           For all solutions in a block, the initial condition at the highest log power
           is the eigenvector; the lower log terms are determined by the recurrence. *)
        With[{eigIdx = blockStarts[[blockIdx]]},
          AppendTo[initialVectors, jordanS[[All, eigIdx]]];
        ];

        , {posInBlock, blockSizes[[blockIdx]]}
      ];
      , {posInClass, Length[resonanceClasses[[classIdx]]]}
    ];
    , {classIdx, Length[resonanceClasses]}
  ];

  (* Normalize initial vectors to working precision *)
  initialVectors = N[initialVectors, ctx["WorkingPrecision"]];

  <|
    "Eigenvalues" -> eigenvalues,
    "JordanS" -> N[jordanS, ctx["WorkingPrecision"]],
    "JordanSInv" -> N[jordanSInv, ctx["WorkingPrecision"]],
    "JordanJ" -> jordanJ,
    "BlockStarts" -> blockStarts,
    "BlockSizes" -> blockSizes,
    "BlockEigenvalues" -> blockEigenvalues,
    "ResonanceClasses" -> resonanceClasses,
    "SolutionExponents" -> solutionExponents,
    "MaxLogPowers" -> maxLogPowers,
    "InitialLogPowers" -> initialLogPowers,
    "InitialVectors" -> initialVectors
  |>
];

(* Solve one step of the recurrence when the matrix L_n may be singular.
   L_n = (lambda + n)I - M_0
   Solves: L_n f = rhs
   If the system is singular:
     - Checks solvability (rhs must be in column space)
     - Returns {solution, nullSpaceVecs, isSingular}
     - The null-space component is left as zero (to be determined later by solvability)
   If non-singular:
     - Returns {solution, {}, False}
*)
SolveRecurrenceStep[Ln_, rhs_, systemSize_, ctx_Association] := Module[
  {det, tolerance, solution, nullVecs, leftNullVecs, residual, projRhs},

  tolerance = 10^(-ctx["ChopPrecision"]/3);

  (* Check if L_n is singular *)
  det = Det[Ln];

  If[Abs[det] > tolerance,
    (* Non-singular: direct solve *)
    solution = LinearSolve[Ln, rhs];
    Return[{solution, {}, False}]
  ];

  (* Singular case *)
  nullVecs = NullSpace[Ln, Tolerance -> tolerance];
  leftNullVecs = NullSpace[Transpose[Ln], Tolerance -> tolerance];

  If[Length[nullVecs] == 0,
    (* Numerically singular but no null space found - treat as non-singular *)
    solution = LinearSolve[Ln, rhs];
    Return[{solution, {}, False}]
  ];

  (* Check solvability: rhs must be orthogonal to left null space *)
  (* Project rhs onto column space *)
  projRhs = rhs;
  Do[
    With[{wNorm = leftNullVecs[[j]] . Conjugate[leftNullVecs[[j]]]},
      If[wNorm > tolerance,
        projRhs = projRhs - (leftNullVecs[[j]] . projRhs / wNorm) * leftNullVecs[[j]];
      ];
    ],
    {j, Length[leftNullVecs]}
  ];

  (* Compute particular solution using pseudo-inverse *)
  solution = PseudoInverse[Ln, Tolerance -> tolerance] . projRhs;

  (* Return solution with null space info *)
  {solution, nullVecs, True}
];

(* Compute the resonant fundamental matrix.
   Returns FMat where FMat[[component, solution]] is the series with log terms.
   Each entry is: x^lambda * (sum_k Logx^k * SeriesData[...])
*)
ComputeResonantFundamentalMatrix[ctx_Association, resInfo_Association,
    mCoeffs_, numMCoeffs_] := Module[
  {systemSize, maxOrd, FMat, lambda, maxK, initK, initVec,
   fCoeffs, Ln, LnData, rhs, stepResult, nullVecs, isSingular,
   freeParams, prevFreeParams, solvCondRhs, overlap,
   residueMat, M0, sol, solSeries, k, n,
   tolerance, wp},

  systemSize = ctx["SystemSize"];
  maxOrd = ctx["ExpansionOrder"];
  wp = ctx["WorkingPrecision"];
  tolerance = 10^(-ctx["ChopPrecision"]/3);

  (* Extract M_0 from the coefficient list *)
  residueMat = Table[
    SeriesCoefficient[ctx["AMatExpanded"][[i, j]], {DiffExp`Symbols`x, 0, -1}],
    {i, systemSize}, {j, systemSize}
  ];
  M0 = N[residueMat, wp];

  (* Pre-compute which values of n produce singular L_n *)
  (* L_n = (lambda + n)I - M_0 *)

  FMat = Table[
    lambda = resInfo["SolutionExponents"][[sol]];
    maxK = resInfo["MaxLogPowers"][[sol]];
    initK = resInfo["InitialLogPowers"][[sol]];
    initVec = N[resInfo["InitialVectors"][[sol]], wp];

    (* Allocate coefficient array: fCoeffs[[n+1, k+1]] = vector of length systemSize *)
    fCoeffs = Table[
      ConstantArray[N[0, wp], systemSize],
      {nn, maxOrd + 1}, {kk, maxK + 1}
    ];

    (* Set initial condition: f_{0, initK} = initVec *)
    fCoeffs[[1, initK + 1]] = initVec;

    (* Process n = 0: solve the chain of initial conditions *)
    Ln = N[(lambda) * IdentityMatrix[systemSize] - M0, wp];
    LnData = {};

    (* For n = 0, iterate k from initK - 1 down to 0 *)
    Do[
      rhs = -(k + 1) * fCoeffs[[1, k + 2]];
      {fCoeffs[[1, k + 1]], nullVecs, isSingular} =
        SolveRecurrenceStep[Ln, rhs, systemSize, ctx];
      ,
      {k, initK - 1, 0, -1}
    ];

    (* Process n = 1, 2, ..., maxOrd *)
    prevFreeParams = {}; (* Track free parameters from previous k within same n *)

    Do[
      Ln = N[(lambda + n) * IdentityMatrix[systemSize] - M0, wp];
      prevFreeParams = {};

      (* Iterate k from maxK down to 0 *)
      Do[
        (* Compute RHS: -(k+1) f_{n,k+1} + sum_{i=1}^n M_i f_{n-i,k} *)
        rhs = If[k < maxK,
          -(k + 1) * fCoeffs[[n + 1, k + 2]],
          ConstantArray[N[0, wp], systemSize]
        ];

        (* Add matrix contribution: sum_{i=1}^{min(n, numMCoeffs)} M_i . f_{n-i, k} *)
        Do[
          rhs += mCoeffs[[i + 1]] . fCoeffs[[n - i + 1, k + 1]];
          , {i, 1, Min[n, numMCoeffs]}
        ];

        (* Solve the system *)
        {fCoeffs[[n + 1, k + 1]], nullVecs, isSingular} =
          SolveRecurrenceStep[Ln, rhs, systemSize, ctx];

        (* If singular and we have a previous free parameter, use solvability to determine it *)
        If[isSingular && Length[prevFreeParams] > 0,
          (* The previous step (k+1) left a free parameter.
             Solvability at this step determines it.
             Recompute: adjust f_{n, k+1} by adding null-space components
             and check if this makes the current step solvable. *)
          Module[{leftNull, fullRhs, residualProj, correction, alpha,
                  prevNull, overlapMat, rhsVec},
            leftNull = NullSpace[Transpose[Ln], Tolerance -> tolerance];
            prevNull = prevFreeParams;

            If[Length[leftNull] > 0 && Length[prevNull] > 0,
              (* The correction to f_{n,k+1} is: sum_j alpha_j * prevNull[[j]] *)
              (* This changes rhs by: -(k+1) * sum_j alpha_j * prevNull[[j]] *)
              (* For solvability: leftNull . (rhs - (k+1) sum_j alpha_j prevNull[[j]]) = 0 *)
              (* Already have: leftNull . rhs = residual (should be ~0 from PseudoInverse) *)
              (* But with the correction: leftNull . correction_rhs = 0 *)
              (* leftNull . rhs_original + leftNull . (-(k+1) prevNull . alpha) = 0 *)

              (* Recompute with the original RHS (before projection) *)
              fullRhs = If[k < maxK,
                -(k + 1) * fCoeffs[[n + 1, k + 2]],
                ConstantArray[N[0, wp], systemSize]
              ];
              Do[
                fullRhs += mCoeffs[[ii + 1]] . fCoeffs[[n - ii + 1, k + 1]];
                , {ii, 1, Min[n, numMCoeffs]}
              ];

              (* Compute overlap: leftNull . prevNull^T and leftNull . fullRhs *)
              overlapMat = Table[
                leftNull[[a]] . prevNull[[b]],
                {a, Length[leftNull]}, {b, Length[prevNull]}
              ] * (-(k + 1));
              rhsVec = Table[leftNull[[a]] . fullRhs, {a, Length[leftNull]}];

              If[Max[Abs[Flatten[overlapMat]]] > tolerance,
                (* Solve for alpha *)
                alpha = Quiet[Check[
                  LinearSolve[overlapMat, -rhsVec],
                  ConstantArray[N[0, wp], Length[prevNull]]
                ]];

                (* Update f_{n, k+1} with the correction *)
                Do[
                  fCoeffs[[n + 1, k + 2]] += alpha[[j]] * prevNull[[j]];
                  , {j, Length[prevNull]}
                ];

                (* Recompute f_{n, k} with corrected RHS *)
                rhs = If[k < maxK,
                  -(k + 1) * fCoeffs[[n + 1, k + 2]],
                  ConstantArray[N[0, wp], systemSize]
                ];
                Do[
                  rhs += mCoeffs[[ii + 1]] . fCoeffs[[n - ii + 1, k + 1]];
                  , {ii, 1, Min[n, numMCoeffs]}
                ];
                {fCoeffs[[n + 1, k + 1]], nullVecs, isSingular} =
                  SolveRecurrenceStep[Ln, rhs, systemSize, ctx];
              ];
            ];
          ];
        ];

        (* Store the null-space vectors for use by the next k level *)
        prevFreeParams = If[isSingular, nullVecs, {}];

        , {k, maxK, 0, -1}
      ];

      , {n, 1, maxOrd}
    ];

    (* Build SeriesData output for this solution *)
    (* FMat entry = x^lambda * sum_k Logx^k * SeriesData[x, 0, coeffs_k, 0, maxOrd+1, 1] *)
    Table[
      DiffExp`Utilities`PChop[
        DiffExp`Symbols`x^lambda * Sum[
          DiffExp`Symbols`Logx^k *
            SeriesData[DiffExp`Symbols`x, 0,
              DiffExp`Utilities`PChop[Table[fCoeffs[[n + 1, k + 1, comp]], {n, 0, maxOrd}]],
              0, maxOrd + 1, 1],
          {k, 0, maxK}
        ]
      ],
      {comp, systemSize}
    ]

    , {sol, systemSize}
  ] // Transpose; (* FMat[[component, solution]] *)

  DiffExp`Utilities`PChop[FMat]
];

(* ============================================================================ *)
(* UNIFIED PARTICULAR SOLUTION SOLVER                                          *)
(* Handles resonant eigenvalues, fractional base exponents, and log terms      *)
(* in the inhomogeneous source bVec. Uses SVD pseudoinverse at resonant steps  *)
(* with dynamic K_max growth and residual verification.                        *)
(* ============================================================================ *)

(* Compute the RHS of the recurrence at step (n, k):
   rhs = beta_{n,k} - (k+1)*f_{n,k+1} + sum_{i=1}^{min(n,numM)} M_i . f_{n-i,k} *)
ComputeRecurrenceRHS[fCoeffs_, bCoeffs_, n_, k_, kMax_, bMaxLogK_,
    mCoeffs_, numMCoeffs_, systemSize_, maxOrd_, wp_] := Module[{rhs},
  (* Source term beta_{n,k} *)
  rhs = If[k <= bMaxLogK && n + 1 <= Length[bCoeffs] && k + 1 <= Length[bCoeffs[[1]]],
    N[bCoeffs[[n + 1, k + 1]], wp],
    ConstantArray[N[0, wp], systemSize]
  ];
  (* Log descent: -(k+1)*f_{n,k+1} *)
  If[k < kMax,
    rhs -= (k + 1) * fCoeffs[[n + 1, k + 2]];
  ];
  (* Matrix convolution: sum M_i . f_{n-i,k} *)
  Do[
    rhs += mCoeffs[[i + 1]] . fCoeffs[[n - i + 1, k + 1]];
    , {i, 1, Min[n, numMCoeffs]}
  ];
  rhs
];

(* Solve the recurrence at a singular step n using the full block system.
   When L_n has a Jordan block structure, the k-by-k approach with PseudoInverse
   fails because the null-space/left-null-space overlap is zero for Jordan blocks
   of size > 1. Instead, we assemble the coupled system across ALL k values:
     L_n f_{n,k} + (k+1) f_{n,k+1} = sourceRHS_k
   and solve it simultaneously.
   Returns a list of vectors: result[[k+1]] = f_{n,k}. *)
SolveSingularNBlock[Ln_, n_, kMax_, bCoeffs_, bMaxLogK_,
    fCoeffs_, mCoeffs_, numMCoeffs_,
    systemSize_, maxOrd_, wp_, tolerance_] := Module[
  {bigSize, bigMatrix, bigRHS, bigSol, result, sourceRHS, k, i},

  bigSize = (kMax + 1) * systemSize;

  (* Build block matrix: L_n on diagonal, (k+1)*I on super-diagonal *)
  bigMatrix = ConstantArray[N[0, wp], {bigSize, bigSize}];
  bigRHS = ConstantArray[N[0, wp], bigSize];

  Do[
    (* Diagonal block at position k: L_n *)
    bigMatrix[[k * systemSize + 1 ;; (k + 1) * systemSize,
               k * systemSize + 1 ;; (k + 1) * systemSize]] = Ln;
    (* Super-diagonal block: (k+1)*I *)
    If[k < kMax,
      bigMatrix[[k * systemSize + 1 ;; (k + 1) * systemSize,
                 (k + 1) * systemSize + 1 ;; (k + 2) * systemSize]] =
        N[(k + 1) * IdentityMatrix[systemSize], wp];
    ];
    (* Source RHS: β_{n,k} + Σ M_i f_{n-i,k} (without the -(k+1)f_{n,k+1} coupling) *)
    sourceRHS = If[k <= bMaxLogK && n + 1 <= Length[bCoeffs] && k + 1 <= Length[bCoeffs[[1]]],
      N[bCoeffs[[n + 1, k + 1]], wp],
      ConstantArray[N[0, wp], systemSize]
    ];
    Do[
      sourceRHS += mCoeffs[[i + 1]] . fCoeffs[[n - i + 1, k + 1]];
      , {i, 1, Min[n, numMCoeffs]}
    ];
    bigRHS[[k * systemSize + 1 ;; (k + 1) * systemSize]] = sourceRHS;
    , {k, 0, kMax}
  ];

  (* Solve using PseudoInverse (handles the singular block system) *)
  bigSol = PseudoInverse[N[bigMatrix, wp], Tolerance -> tolerance] . bigRHS;

  (* Extract results: result[[k+1]] = f_{n,k} *)
  result = Table[
    bigSol[[k * systemSize + 1 ;; (k + 1) * systemSize]],
    {k, 0, kMax}
  ];

  result
];

(* Run the core particular recurrence for a given base exponent sigma.
   At non-singular steps (L_n invertible): direct back-substitution k-by-k.
   At singular steps (L_n has zero eigenvalues): full block system solve.
   Returns fCoeffs[[n+1, k+1]] = vector of length systemSize. *)
RunParticularRecurrence[sigma_, bCoeffs_, bMaxLogK_, kMax_,
    M0_, mCoeffs_, numMCoeffs_,
    systemSize_, maxOrd_, wp_, tolerance_, ctx_] := Module[
  {fCoeffs, Ln, rhs, det, blockResult, n, k},

  fCoeffs = Table[
    ConstantArray[N[0, wp], systemSize],
    {maxOrd + 1}, {kMax + 1}
  ];

  Do[
    Ln = N[(sigma + n) * IdentityMatrix[systemSize] - M0, wp];
    det = Det[Ln];

    If[Abs[det] < tolerance,
      (* Singular step: solve the full block system across all k *)
      blockResult = SolveSingularNBlock[Ln, n, kMax, bCoeffs, bMaxLogK,
        fCoeffs, mCoeffs, numMCoeffs,
        systemSize, maxOrd, wp, tolerance];
      Do[
        fCoeffs[[n + 1, k + 1]] = blockResult[[k + 1]];
        , {k, 0, kMax}
      ];
      ,
      (* Non-singular step: direct back-substitution *)
      Do[
        rhs = ComputeRecurrenceRHS[fCoeffs, bCoeffs, n, k, kMax,
          bMaxLogK, mCoeffs, numMCoeffs, systemSize, maxOrd, wp];
        fCoeffs[[n + 1, k + 1]] = LinearSolve[Ln, rhs];
        , {k, kMax, 0, -1}
      ];
    ];

    , {n, 0, maxOrd}
  ];

  fCoeffs
];

(* Check residual at all resonant steps to verify solution quality.
   Returns the maximum residual norm across all singular steps. *)
CheckParticularResidual[fCoeffs_, sigma_, bCoeffs_, bMaxLogK_, kMax_,
    M0_, mCoeffs_, numMCoeffs_,
    systemSize_, maxOrd_, wp_, tolerance_] := Module[
  {maxResidual = 0, Ln, rhs, residual, n, k},

  Do[
    Ln = N[(sigma + n) * IdentityMatrix[systemSize] - M0, wp];
    (* Only check at potentially singular steps *)
    If[Abs[Det[Ln]] < tolerance,
      Do[
        rhs = ComputeRecurrenceRHS[fCoeffs, bCoeffs, n, k, kMax,
          bMaxLogK, mCoeffs, numMCoeffs, systemSize, maxOrd, wp];
        residual = Ln . fCoeffs[[n + 1, k + 1]] - rhs;
        maxResidual = Max[maxResidual, Max[Abs[residual]]];
        , {k, kMax, 0, -1}
      ];
    ];
    , {n, 0, maxOrd}
  ];

  maxResidual
];

(* Build the SeriesData output from the coefficient array for a given sigma. *)
BuildParticularFromCoeffs[fCoeffs_, sigma_, kMax_,
    systemSize_, maxOrd_, wp_] := Module[{fParticular},

  fParticular = Table[
    DiffExp`Utilities`PChop[
      DiffExp`Symbols`x^sigma * Sum[
        DiffExp`Symbols`Logx^k *
          SeriesData[DiffExp`Symbols`x, 0,
            DiffExp`Utilities`PChop[Table[fCoeffs[[n + 1, k + 1, comp]], {n, 0, maxOrd}]],
            0, maxOrd + 1, 1],
        {k, 0, kMax}
      ]
    ],
    {comp, systemSize}
  ];

  DiffExp`Utilities`PChop[fParticular]
];

(* Unified particular solution solver.
   Handles resonant eigenvalues, fractional base exponents, Logx in bVec.
   Uses SVD pseudoinverse at resonant steps with residual-driven K_max growth.

   Arguments:
     bVec_        - inhomogeneous source vector (may contain Logx, x^sigma)
     resInfo_     - resonance structure (from ComputeResonanceStructure)
     mCoeffs_     - matrix coefficients M_i (M_0 = residue, M_1 = A_0, etc.)
     numMCoeffs_  - number of M coefficients available
     ctx_         - SegmentContext association

   Returns: particular solution vector (list of series expressions), or $Failed *)
ComputeUnifiedParticular[bVec_, resInfo_Association,
    mCoeffs_, numMCoeffs_, ctx_Association] := Module[
  {systemSize, maxOrd, wp, tolerance,
   M0, eigenvalues,
   bLeadPow, bMaxLogK, bCoeffs,
   sigma, kMax, kMaxInitial, kMaxGrowthAttempts,
   fCoeffs, maxResidual, fParticular},

  systemSize = Length[bVec];
  maxOrd = ctx["ExpansionOrder"];
  wp = ctx["WorkingPrecision"];
  tolerance = 10^(-ctx["ChopPrecision"]/3);

  (* Extract M_0 (residue matrix) *)
  M0 = N[mCoeffs[[1]], wp];
  eigenvalues = resInfo["SolutionExponents"];

  (* Extract log structure of bVec *)
  bMaxLogK = Max[0, Max[Table[
    DiffExp`SeriesOps`MaxLogxPower[bVec[[i]]],
    {i, systemSize}
  ]]];

  (* Determine leading power of bVec across all components and log orders *)
  bLeadPow = Min[Table[
    Module[{ser, minPow = maxOrd},
      Do[
        ser = DiffExp`SeriesOps`LogxCoeff[bVec[[comp]], logk];
        If[Head[ser] === SeriesData && ser =!= 0,
          minPow = Min[minPow, ser[[4]] / ser[[6]]];
        ];
        , {logk, 0, bMaxLogK}
      ];
      minPow
    ],
    {comp, systemSize}
  ]];

  (* The particular solution base exponent: sigma = bLeadPow + 1
     (from xf' = M(x)f + xB(x): the x*B shifts power by +1) *)
  sigma = bLeadPow + 1;

  (* Initial K_max: log depth of source + max Jordan block size among resonating eigenvalues.
     MaxLogPowers = blockSize - 1, so we add +1 to get blockSize.
     This ensures the block system at singular n is fully solvable. *)
  kMaxInitial = bMaxLogK + Max[resInfo["MaxLogPowers"]] + 1;
  kMax = kMaxInitial;

  (* Extract beta coefficients: beta_{n,k,comp} = coeff of x^{sigma-1+n} (Logx)^k in bVec *)
  bCoeffs = Table[
    Module[{logCoeff},
      logCoeff = DiffExp`SeriesOps`LogxCoeff[bVec[[comp]], logk];
      If[Head[logCoeff] === SeriesData,
        SeriesCoefficient[logCoeff, {DiffExp`Symbols`x, 0, sigma - 1 + n}],
        If[n == 0 && logk == 0 && NumericQ[logCoeff], logCoeff, 0]
      ]
    ],
    {n, 0, maxOrd}, {logk, 0, bMaxLogK}, {comp, systemSize}
  ];

  DiffExp`Utilities`PrintInfo["  Unified particular: sigma = ", sigma,
    ", K_b = ", bMaxLogK, ", K_max = ", kMax, "."][3];

  (* Run recurrence with residual checking and dynamic K_max growth *)
  kMaxGrowthAttempts = 0;
  While[True,
    fCoeffs = RunParticularRecurrence[sigma, bCoeffs, bMaxLogK, kMax,
      M0, mCoeffs, numMCoeffs,
      systemSize, maxOrd, wp, tolerance, ctx];

    (* Check residual at resonant steps *)
    maxResidual = CheckParticularResidual[fCoeffs, sigma, bCoeffs, bMaxLogK, kMax,
      M0, mCoeffs, numMCoeffs,
      systemSize, maxOrd, wp, tolerance];

    If[maxResidual < tolerance,
      (* Solution is good *)
      Break[];
    ];

    (* Residual too large: try increasing K_max *)
    kMaxGrowthAttempts++;
    If[kMaxGrowthAttempts > 3,
      DiffExp`Utilities`PrintInfo[
        "  Unified particular: residual ", maxResidual,
        " still too large after 3 K_max increases. Returning $Failed."][3];
      Return[$Failed]
    ];

    kMax += 1;
    DiffExp`Utilities`PrintInfo[
      "  Unified particular: residual ", maxResidual,
      " > tolerance. Increasing K_max to ", kMax, "."][3];
  ];

  (* Update IMaxLogOrder if our particular solution has higher log powers *)
  If[kMax > DiffExp`State`IMaxLogOrder,
    DiffExp`Integration`UpdateIntReps[kMax];
  ];

  (* Build and return the particular solution *)
  BuildParticularFromCoeffs[fCoeffs, sigma, kMax, systemSize, maxOrd, wp]
];

(* Extract the M_i coefficients (M(x) = x*A(x) = sum M_i x^i) from the expanded matrix.
   M_0 is the residue, M_1 = A_0, M_2 = A_1, etc.
   Returns {mCoeffs, numMCoeffs} where mCoeffs[[i+1]] = M_i. *)
ExtractMCoefficients[ctx_Association] := Module[
  {AMatExpanded, systemSize, maxOrd, mCoeffs},

  systemSize = ctx["SystemSize"];
  maxOrd = ctx["ExpansionOrder"];
  AMatExpanded = ctx["AMatExpanded"];

  (* M_i = coefficient of x^i in x*A(x) = coefficient of x^{i-1} in A(x) *)
  mCoeffs = Table[
    Table[
      SeriesCoefficient[AMatExpanded[[row, col]], {DiffExp`Symbols`x, 0, i - 1}],
      {row, systemSize}, {col, systemSize}
    ],
    {i, 0, maxOrd}
  ];

  (* Convert to working precision *)
  mCoeffs = N[mCoeffs, ctx["WorkingPrecision"]];

  {mCoeffs, maxOrd}
];

(* Main entry point for the general singular recurrence solver.
   Handles regular singular points with any eigenvalue structure
   (resonant, non-diagonalizable, etc.) *)
SolveGeneralSingularRecurrence[ctx_Association, bVec_, epsord_, cacheIn_Association] := Module[
  {systemSize, maxOrd, FMat, fParticular, fGeneral, cIndices, c,
   resInfo, mCoeffs, numMCoeffs,
   cache = cacheIn},

  systemSize = ctx["SystemSize"];
  maxOrd = ctx["ExpansionOrder"];

  If[epsord === 0 && !KeyExistsQ[cache, "GenSingRR"],
    (* Extract matrix coefficients M_i *)
    {mCoeffs, numMCoeffs} = ExtractMCoefficients[ctx];

    (* Compute resonance structure *)
    resInfo = ComputeResonanceStructure[
      Table[
        SeriesCoefficient[ctx["AMatExpanded"][[i, j]], {DiffExp`Symbols`x, 0, -1}],
        {i, systemSize}, {j, systemSize}
      ],
      systemSize, ctx
    ];

    DiffExp`Utilities`PrintInfo["Using general singular recurrence for integrals ", ctx["Label"],
      " (eigenvalues: ", resInfo["SolutionExponents"], ", max log powers: ",
      resInfo["MaxLogPowers"], ")."][3];

    (* Update IMaxLogOrder if our solutions will have higher log powers *)
    With[{maxLogK = Max[resInfo["MaxLogPowers"]]},
      If[maxLogK > DiffExp`State`IMaxLogOrder,
        DiffExp`Integration`UpdateIntReps[maxLogK];
      ];
    ];

    (* Compute fundamental matrix *)
    FMat = ComputeResonantFundamentalMatrix[ctx, resInfo, mCoeffs, numMCoeffs];

    cache["GenSingRR"] = {FMat, resInfo, mCoeffs, numMCoeffs};
  ];

  {FMat, resInfo, mCoeffs, numMCoeffs} = cache["GenSingRR"];

  (* Compute particular solution using the unified solver *)
  If[DiffExp`Utilities`PChop[bVec] === ConstantArray[0, systemSize],
    fParticular = ConstantArray[
      SeriesData[DiffExp`Symbols`x, 0, {}, 0, maxOrd + 1, 1],
      systemSize
    ];
    ,
    fParticular = ComputeUnifiedParticular[bVec, resInfo,
      mCoeffs, numMCoeffs, ctx];

    (* If unified solver fails, fall back to Default strategy *)
    If[fParticular === $Failed,
      DiffExp`Utilities`PrintInfo[
        "Unified particular solver failed for integrals ", ctx["Label"],
        ". Falling back to Default strategy."][3];
      Module[{defaultCache, defaultResult, defaultCIndices, defaultFGeneral},
        defaultCache = Lookup[cache, "DefaultFallback", <||>];
        If[!KeyExistsQ[defaultCache, "FMat"],
          defaultResult = SolveDefault[ctx,
            ConstantArray[SeriesData[DiffExp`Symbols`x, 0, {}, 0, maxOrd + 1, 1], systemSize],
            0, defaultCache];
          defaultCache = defaultResult[[3]];
        ];
        defaultResult = SolveDefault[ctx, bVec, epsord, defaultCache];
        cache["DefaultFallback"] = defaultResult[[3]];
        {defaultCIndices, defaultFGeneral} = defaultResult[[{1, 2}]];
        fParticular = defaultFGeneral /. Thread[defaultCIndices -> 0];
        fParticular = DiffExp`SeriesOps`SExpand[fParticular];
      ];
    ];
  ];

  (* General solution: f = particular + sum_i c_i * FMat_column_i *)
  cIndices = Table[Subscript[c, i], {i, systemSize}];
  fGeneral = fParticular + Sum[cIndices[[i]] * FMat[[All, i]], {i, systemSize}];
  fGeneral = DiffExp`SeriesOps`SExpand[fGeneral];

  {cIndices, fGeneral, cache}
];
