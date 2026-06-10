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

ExpandedMatrixMinOrder[AMatExpanded_, systemSize_Integer] := Module[{minOrders},
  minOrders = Flatten[Table[
    If[Head[AMatExpanded[[i, j]]] === SeriesData,
      AMatExpanded[[i, j]][[4]] / AMatExpanded[[i, j]][[6]],
      0
    ],
    {i, systemSize}, {j, systemSize}
  ]];
  Min[minOrders]
];

TrimNonFiniteSeriesTails[expr_, ctx_Association, epsord_, source_String] :=
Module[{trimmed = {}, failed = False, result, trimOne},
  trimOne[s_SeriesData] := Module[
    {coeffs = s[[3]], badFlags, badIndices, firstBad, tailIndices,
     newMax},

    badFlags = DiffExp`Utilities`NonFiniteExpressionQ /@ coeffs;
    badIndices = Flatten@Position[badFlags, True];
    If[badIndices === {}, Return[s]];

    firstBad = First[badIndices];
    tailIndices = Range[firstBad, Length[coeffs]];
    If[badIndices =!= tailIndices,
      failed = True;
      DiffExp`State`LastErrorContext = <|
        "Context" -> ctx,
        "EpsilonOrder" -> epsord,
        "Source" -> source,
        "Series" -> s,
        "BadCoefficientIndices" -> badIndices
      |>;
      Return[s];
    ];

    newMax = s[[4]] + firstBad - 1;
    AppendTo[
      trimmed,
      <|
        "Source" -> source,
        "EpsilonOrder" -> epsord,
        "Label" -> Lookup[ctx, "Label", Missing["Unknown"]],
        "FirstDroppedPower" -> newMax/s[[6]],
        "DroppedCount" -> Length[badIndices]
      |>
    ];
    SeriesData[s[[1]], s[[2]], Take[coeffs, firstBad - 1],
      s[[4]], newMax, s[[6]]]
  ];

  result = expr /. s_SeriesData :> trimOne[s];
  If[failed,
    DiffExp`Utilities`ReportError[
      "General singular recurrence produced non-finite local coefficients ",
      "before the truncation tail for integral(s) ",
      Lookup[ctx, "Label", Missing["Unknown"]],
      " at epsilon order ",
      epsord,
      "."
    ];
  ];
  If[Length[trimmed] > 0,
    DiffExp`Utilities`PrintWarning[
      "General singular recurrence dropped trailing non-finite local-series ",
      "coefficient tail(s): ",
      trimmed
    ];
  ];
  result
];

FuchsianizedSingularRecurrenceApplicableQ[ctx_Association] := Quiet[Check[
  Module[{minOrder},
    If[MissingQ[ctx["AMatFactored"]], Return[False]];
    minOrder = ExpandedMatrixMinOrder[ctx["AMatExpanded"], ctx["SystemSize"]];
    minOrder < -1
  ],
  False
]];

BuildFuchsianizedRecurrenceData[ctx_Association] := Module[
  {x = DiffExp`Symbols`x, maxOrd, systemSize, thetaMat, fb, T, B,
   TInv, AMatFactoredG, AMatExpandedG, ctxG},

  If[MissingQ[ctx["AMatFactored"]],
    Return[$Failed]
  ];

  maxOrd = ctx["ExpansionOrder"];
  systemSize = ctx["SystemSize"];
  thetaMat = Together[x ctx["AMatFactored"]];

  (* The lattice saturation and the trim pass decide pole cancellations with
     exact zero tests (RatCancel pivots, minimum-valuation checks).  On
     numeric coefficients those cancellations are never exact, so the trim
     test keeps "succeeding" on noise and multiplies lattice columns by x
     once per pass, ratcheting T up to x^50-scale entries while the working
     precision collapses.  Rationalize the local matrix (its entries are
     exact data evaluated at the numeric segment center) and run the whole
     reduction in exact arithmetic. *)
  thetaMat = Rationalize[
    thetaMat,
    10^(-Max[20, ctx["WorkingPrecision"] - 50])
  ];

  fb = DiffExp`LocalSeries`FuchsianizeLocal[
    thetaMat,
    x,
    "MaxSteps" -> 200
  ];

  If[MatchQ[fb, _Failure],
    Return[$Failed]
  ];

  {T, B} = fb;
  {T, B} = DiffExp`LocalSeries`TrimFuchsianLattice[thetaMat, T, x];
  TInv = Together[Inverse[T]];
  AMatFactoredG = Together[B/x];
  AMatExpandedG = Table[
    Quiet[Series[AMatFactoredG[[i, j]], {x, 0, maxOrd}]],
    {i, systemSize}, {j, systemSize}
  ];

  If[ExpandedMatrixMinOrder[AMatExpandedG, systemSize] =!= -1,
    Return[$Failed]
  ];

  (* Guard against degenerate lattices: a local balance only needs power
     shifts of the order of the system size plus the pole depth.  Anything
     much larger means the reduction went numerically astray. *)
  Module[{maxDeg},
    maxDeg = Max[Map[
      Function[entry,
        If[TrueQ[PossibleZeroQ[entry]],
          0,
          Module[{tog = Together[entry]},
            Max[
              Exponent[Numerator[tog], x],
              Exponent[Denominator[tog], x]
            ]
          ]
        ]
      ],
      T, {2}
    ]];
    If[!TrueQ[maxDeg <= 4 systemSize + 8],
      DiffExp`Utilities`PrintWarning[
        "BuildFuchsianizedRecurrenceData: local Fuchsianizing lattice has ",
        "entries of x-degree ", maxDeg,
        "; refusing the degenerate balance."
      ];
      Return[$Failed]
    ];
  ];

  ctxG = Join[
    ctx,
    <|
      "AMatFactored" -> AMatFactoredG,
      "AMatExpanded" -> AMatExpandedG,
      "Label" -> ctx["Label"]
    |>
  ];

  <|"T" -> T, "TInv" -> TInv, "Ctx" -> ctxG, "Subcache" -> <||>|>
];

(* Toleranced Jordan data for a numeric residue matrix with a degenerate
   snapped spectrum.  Builds generalized-eigenvector chains per snapped
   eigenvalue from toleranced nested null spaces of (M - lambda)^k, so that
   float-level splitting of degenerate eigenvalues cannot leak ill-scaled
   exact eigenvectors into the fundamental system.  Returns {S, J} in the
   same convention as JordanDecomposition (chains ordered eigenvector
   first; J carries the snapped eigenvalues with 1's on the superdiagonal
   inside each chain), or $Failed. *)
NumericJordanData[residueMat_, snappedSpectrum_, wp_, tolerance_] := Quiet[Check[
  Catch[
    Module[{M, d, clusters, columns = {}, jdiag = {}, jsuper = {},
            projectOut, svdNullBasis, S, J},

      M = N[residueMat, wp];
      d = Length[M];
      clusters = Tally[snappedSpectrum];

      (* Orthonormal toleranced null basis via the right singular vectors;
         NullSpace[..., Tolerance] is unreliable for high-precision input
         (it can underflow to machine arithmetic and return nothing). *)
      svdNullBasis[mat_] := Module[{u, s, v, svals, smax},
        {u, s, v} = SingularValueDecomposition[N[mat, wp]];
        svals = Diagonal[s];
        smax = Max[Append[Abs[svals], 0]];
        If[TrueQ[smax == 0],
          IdentityMatrix[d],
          Select[
            Transpose[v],
            Function[col,
              TrueQ[Norm[mat . col] <= Max[tolerance smax, tolerance]]
            ]
          ]
        ]
      ];

      (* Component of vec orthogonal to span(vs); vs need not be
         orthonormal (orthonormalized on the fly). *)
      projectOut[vec_, vs_] := Module[{work = vec, basis},
        basis = Orthogonalize[vs, Method -> "GramSchmidt"];
        basis = DeleteCases[basis, v_ /; Norm[v] < 1/2];
        Do[
          work = work - (Conjugate[b] . work) b,
          {b, basis}
        ];
        work
      ];

      Do[
        Module[{lambda = cluster[[1]], mult = cluster[[2]], Mshift,
                powers, nullBases, depth, chains = {}, levelMembers,
                tops, cand, k},
          Mshift = M - lambda IdentityMatrix[d];

          (* Nested toleranced null spaces of Mshift^k up to multiplicity. *)
          nullBases = {};
          powers = IdentityMatrix[d];
          depth = 0;
          While[depth < mult + 1,
            depth++;
            powers = powers . Mshift;
            AppendTo[nullBases, svdNullBasis[powers]];
            If[Length[Last[nullBases]] >= mult, Break[]];
          ];
          If[Length[Last[nullBases]] != mult,
            Throw[$Failed, "numericJordan"]
          ];
          depth = Length[nullBases];

          (* levelMembers[k] collects vectors already claimed at level k by
             deeper chains; new level-k tops must be independent of these
             and of Null(Mshift^(k-1)). *)
          levelMembers = ConstantArray[{}, depth];

          Do[
            Do[
              cand = projectOut[
                candidate,
                Join[
                  If[k > 1, nullBases[[k - 1]], {}],
                  levelMembers[[k]]
                ]
              ];
              If[Norm[cand] > tolerance^(1/3),
                Module[{top = cand / Norm[cand], chain},
                  chain = NestList[Mshift . # &, top, k - 1];
                  (* chain = {v_k, v_(k-1), ..., v_1}; the deepest entry is
                     the eigenvector.  Reject defective chains whose lower
                     members degenerate. *)
                  If[Min[Norm /@ chain] < tolerance^(1/3),
                    Throw[$Failed, "numericJordan"]
                  ];
                  Do[
                    levelMembers[[k - j + 1]] =
                      Append[levelMembers[[k - j + 1]], chain[[j]]],
                    {j, k}
                  ];
                  AppendTo[chains, {k, Reverse[chain]}];
                ];
              ],
              {candidate, nullBases[[k]]}
            ],
            {k, depth, 1, -1}
          ];

          If[Total[First /@ chains] != mult,
            Throw[$Failed, "numericJordan"]
          ];

          Do[
            Module[{k = ch[[1]], cols = ch[[2]]},
              columns = Join[columns, cols];
              jdiag = Join[jdiag, ConstantArray[lambda, k]];
              jsuper = Join[jsuper, Append[ConstantArray[1, k - 1], 0]];
            ],
            {ch, chains}
          ];
        ],
        {cluster, clusters}
      ];

      If[Length[columns] != d, Throw[$Failed, "numericJordan"]];

      S = Transpose[columns];
      If[Abs[Det[S]] < tolerance, Throw[$Failed, "numericJordan"]];

      J = DiagonalMatrix[jdiag] + DiagonalMatrix[Most[jsuper], 1];

      {S, J}
    ],
    "numericJordan"
  ],
  $Failed
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
  {eigenvalues, rawEigenvalues, jordanS, jordanJ, jordanSInv, jordanData,
   residueForJordan, snappedSpectrum, degenerateSnapped,
   blockStarts, blockSizes, blockEigenvalues,
   resonanceClasses, visited, class, classIdx,
   solutionExponents, maxLogPowers, initialLogPowers, initialVectors,
   solIdx, blockIdx, posInBlock, currentEV,
   nullVecs, genEigVecs, tolerance},

  tolerance = 10^(-ctx["ChopPrecision"]/2);

  (* The residue of a local system evaluated at a numeric segment center
     carries float-level noise that splits degenerate eigenvalues by roughly
     the square root of the noise.  Exact JordanDecomposition on (a
     rationalization of) such a matrix sees distinct eigenvalues with
     ill-scaled eigenvectors and returns a useless similarity matrix.
     Detect spectral degeneracy from the snapped spectrum and build the
     Jordan data with toleranced generalized-nullspace chains instead. *)
  snappedSpectrum = Quiet[Check[
    Rationalize[
      Eigenvalues[N[residueMat, ctx["WorkingPrecision"]]],
      tolerance
    ],
    $Failed
  ]];
  degenerateSnapped = ListQ[snappedSpectrum] &&
    Max[Last /@ Tally[snappedSpectrum]] > 1;

  jordanData = If[degenerateSnapped,
    NumericJordanData[residueMat, snappedSpectrum,
      ctx["WorkingPrecision"], tolerance],
    $Failed
  ];

  If[!(ListQ[jordanData] && Length[jordanData] === 2 &&
       MatrixQ[jordanData[[1]]] && MatrixQ[jordanData[[2]]]),
    (* Compute Jordan decomposition. Numeric JordanDecomposition can fail for
       exactly rational residues at high precision, so try the rationalized
       residue first and only then fall back to a numerical decomposition. *)
    residueForJordan = Rationalize[residueMat, tolerance];
    jordanData = Quiet[Check[JordanDecomposition[residueForJordan], $Failed]];
  ];
  If[!(ListQ[jordanData] && Length[jordanData] === 2 &&
       MatrixQ[jordanData[[1]]] && MatrixQ[jordanData[[2]]]),
    jordanData = Quiet[Check[
      JordanDecomposition[N[residueMat, ctx["WorkingPrecision"]]],
      $Failed
    ]];
  ];
  If[!(ListQ[jordanData] && Length[jordanData] === 2 &&
       MatrixQ[jordanData[[1]]] && MatrixQ[jordanData[[2]]]),
    DiffExp`State`LastErrorContext = {ctx, residueMat, residueForJordan};
    DiffExp`Utilities`ReportError[
      "General singular recurrence could not compute a Jordan decomposition of the residue matrix for integral(s) ",
      ctx["Label"],
      ". Refusing to fall back to the default Wronskian/Frobenius path."
    ];
  ];
  {jordanS, jordanJ} = jordanData;
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

SafeNumericLinearSolve[mat_, rhs_, wp_, tolerance_] := Module[
  {matN, rhsN, sol, residual, scale},

  matN = N[mat, wp];
  rhsN = N[rhs, wp];

  (* Catastrophic cancellations upstream leave zero-accuracy zeros (0``a)
     in the data, which drag the overall Precision of the input to 0.
     LinearSolve and LeastSquares then FAIL or - worse - return a zero
     "solution" for a manifestly nonzero rhs, silently dropping a source
     term of the recurrence (the banana x = 1/2 eps^-1 impurity).  The
     recurrence data is exact at the working precision by construction,
     so re-fix the precision instead of letting significance arithmetic
     veto the solve. *)
  If[Precision[matN] < wp, matN = SetPrecision[matN, wp]];
  If[Precision[rhsN] < wp, rhsN = SetPrecision[rhsN, wp]];

  sol = Quiet@Check[
    LinearSolve[matN, rhsN],
    $Failed,
    {LinearSolve::nosol, LinearSolve::luc, LinearSolve::sing}
  ];

  If[sol === $Failed || !FreeQ[sol, LinearSolve],
    sol = Quiet@Check[
      LeastSquares[matN, rhsN],
      $Failed
    ];
  ];

  If[sol === $Failed || !FreeQ[sol, LinearSolve | LeastSquares],
    sol = Quiet@Check[
      PseudoInverse[matN, Tolerance -> tolerance] . rhsN,
      $Failed
    ];
  ];

  If[sol === $Failed || !FreeQ[sol, LinearSolve | LeastSquares | PseudoInverse],
    DiffExp`State`LastErrorContext = {mat, rhs, matN, rhsN, sol};
    DiffExp`Utilities`ReportError[
      "General singular recurrence failed to solve a finite-width linear system."
    ];
  ];

  residual = Quiet@Check[matN . sol - rhsN, ConstantArray[0, Dimensions[rhsN]]];
  scale = 1 + Norm[Flatten[N[rhsN, wp]]];
  If[NumericQ[scale] && scale > 0 &&
     NumericQ[Norm[Flatten[N[residual, wp]]]] &&
     Norm[Flatten[N[residual, wp]]]/scale > 100 tolerance,
    DiffExp`Utilities`PrintWarning[
      "SafeNumericLinearSolve: linear solve residual ",
      N[Norm[Flatten[N[residual, wp]]]/scale, 6],
      " exceeds the tolerance scale; the returned solution may be a ",
      "least-squares projection rather than an exact solve."
    ];
  ];

  sol
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
    solution = SafeNumericLinearSolve[Ln, rhs, ctx["WorkingPrecision"], tolerance];
    Return[{solution, {}, False}]
  ];

  (* Singular case *)
  nullVecs = NullSpace[Ln, Tolerance -> tolerance];
  leftNullVecs = NullSpace[Transpose[Ln], Tolerance -> tolerance];

  If[Length[nullVecs] == 0,
    (* Numerically singular but no null space found - treat as non-singular *)
    solution = SafeNumericLinearSolve[Ln, rhs, ctx["WorkingPrecision"], tolerance];
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

    (* For n = 0, iterate k from initK - 1 down to 0.  L_0 is singular by
       construction, so every chain step leaves kernel freedom; the
       minimal-norm pseudo-inverse choice strips kernel components that the
       true generalized-eigenvector chain needs, making the NEXT step's
       right-hand side inconsistent (which SolveRecurrenceStep would
       silently project away).  Resolve each step's freedom from the next
       step's solvability, exactly like the n >= 1 loop below. *)
    prevFreeParams = {};
    Do[
      rhs = -(k + 1) * fCoeffs[[1, k + 2]];
      {fCoeffs[[1, k + 1]], nullVecs, isSingular} =
        SolveRecurrenceStep[Ln, rhs, systemSize, ctx];

      If[isSingular && Length[prevFreeParams] > 0,
        Module[{leftNull, fullRhs, alpha, prevNull, overlapMat,
                overlapMagnitude, rhsVec},
          leftNull = NullSpace[Transpose[Ln], Tolerance -> tolerance];
          prevNull = prevFreeParams;
          If[Length[leftNull] > 0 && Length[prevNull] > 0,
            fullRhs = -(k + 1) * fCoeffs[[1, k + 2]];
            overlapMat = Table[
              leftNull[[a]] . prevNull[[b]],
              {a, Length[leftNull]}, {b, Length[prevNull]}
            ] * (-(k + 1));
            rhsVec = Table[leftNull[[a]] . fullRhs, {a, Length[leftNull]}];
            overlapMagnitude = DiffExp`Utilities`FiniteAbsMax[overlapMat];
            If[overlapMagnitude > tolerance,
              alpha = PseudoInverse[
                N[overlapMat, wp], Tolerance -> tolerance
              ] . (-rhsVec);
              Do[
                fCoeffs[[1, k + 2]] += alpha[[j]] * prevNull[[j]];
                , {j, Length[prevNull]}
              ];
              rhs = -(k + 1) * fCoeffs[[1, k + 2]];
              {fCoeffs[[1, k + 1]], nullVecs, isSingular} =
                SolveRecurrenceStep[Ln, rhs, systemSize, ctx];
            ];
          ];
        ];
      ];

      prevFreeParams = If[isSingular, nullVecs, {}];
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
                  prevNull, overlapMat, overlapMagnitude, rhsVec},
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

              overlapMagnitude = DiffExp`Utilities`FiniteAbsMax[overlapMat];
              If[overlapMagnitude > tolerance,
                (* Solve for alpha *)
                alpha = Quiet[Check[
                  SafeNumericLinearSolve[overlapMat, -rhsVec, wp, tolerance],
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
      If[TrueQ[DiffExp`State`$DebugBlockResidualSeries],
        Print["SINGSTEP sigma=", sigma, " n=", n, " |det|=",
          ScientificForm[N[Abs[det], 3]]];
      ];
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
        fCoeffs[[n + 1, k + 1]] =
          SafeNumericLinearSolve[Ln, rhs, wp, tolerance];
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
  {maxResidual = 0, Ln, rhs, residual, residualMagnitude, n, k},

  (* Check EVERY step, not only the potentially singular ones: a
     non-singular step whose linear solve silently degraded (e.g. on
     precision-collapsed input) is otherwise invisible here, and a
     dropped source term corrupts everything downstream (the banana
     x = 1/2 eps^-1 impurity). *)
  Do[
    Ln = N[(sigma + n) * IdentityMatrix[systemSize] - M0, wp];
    Do[
      rhs = ComputeRecurrenceRHS[fCoeffs, bCoeffs, n, k, kMax,
        bMaxLogK, mCoeffs, numMCoeffs, systemSize, maxOrd, wp];
      residual = Ln . fCoeffs[[n + 1, k + 1]] - rhs;
      (* Analytic-continuation theta symbols may survive in the source;
         DiffExp`Utilities`FiniteAbsMax would silently treat such symbolic residuals as 0,
         blinding this check.  Evaluate both branches instead. *)
      residualMagnitude = If[
        FreeQ[residual, DiffExp`Symbols`\[Theta]p | DiffExp`Symbols`\[Theta]m],
        DiffExp`Utilities`FiniteAbsMax[residual],
        Max[
          DiffExp`Utilities`FiniteAbsMax[residual /. {
            DiffExp`Symbols`\[Theta]p -> 1, DiffExp`Symbols`\[Theta]m -> 0}],
          DiffExp`Utilities`FiniteAbsMax[residual /. {
            DiffExp`Symbols`\[Theta]p -> 0, DiffExp`Symbols`\[Theta]m -> 1}]
        ]
      ];
      maxResidual = DiffExp`Utilities`FiniteAbsMax[{maxResidual, residualMagnitude}];
      , {k, kMax, 0, -1}
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

NormalizeSourcePower[p_, ctx_Association] := Module[
  {tol = 10^(-ctx["ChopPrecision"]/2)},
  Quiet@Rationalize[N[p, ctx["WorkingPrecision"]], tol]
];

SourcePowerClass[p_, ctx_Association] := Mod[NormalizeSourcePower[p, ctx], 1];

MergeSourcePower[p_, groups_Association, ctx_Association] := Module[
  {power = NormalizeSourcePower[p, ctx], key},
  key = SourcePowerClass[power, ctx];
  If[!KeyExistsQ[groups, key] || power < groups[key],
    Join[groups, <|key -> power|>],
    groups
  ]
];

SourcePowersInExpression[expr_, ctx_Association] := Module[
  {x = DiffExp`Symbols`x, e, head, coeffs, ns, factorPowers},

  e = DiffExp`Utilities`PChop[expr];

  Which[
    e === 0,
    {},

    Head[e] === SeriesData,
    coeffs = e[[3]];
    ns = e[[4]] + Range[0, Length[coeffs] - 1];
    NormalizeSourcePower[#, ctx] & /@
      Pick[ns/e[[6]], !TrueQ[DiffExp`Utilities`PChop[#] === 0] & /@ coeffs],

    Head[e] === Plus,
    DeleteDuplicates[Flatten[SourcePowersInExpression[#, ctx] & /@ (List @@ e)]],

    MatchQ[e, Power[x, _]],
    {NormalizeSourcePower[e[[2]], ctx]},

    Head[e] === Times,
    factorPowers = SourcePowersInExpression[#, ctx] & /@ (List @@ e);
    If[MemberQ[factorPowers, {}],
      {},
      DeleteDuplicates[NormalizeSourcePower[Total[#], ctx] & /@ Tuples[factorPowers]]
    ],

    FreeQ[e, x],
    {0},

    True,
    DeleteDuplicates[Flatten[SourcePowersInExpression[#, ctx] & /@ (List @@ Expand[e])]]
  ]
];

SourceBasePowers[bVec_, bMaxLogK_, ctx_Association] := Module[
  {groups = <||>, powers},

  powers = DeleteDuplicates[Flatten[Table[
    SourcePowersInExpression[
      DiffExp`SeriesOps`LogxCoeff[bVec[[comp]], logk],
      ctx
    ],
    {comp, Length[bVec]}, {logk, 0, bMaxLogK}
  ]]];

  Do[
    groups = MergeSourcePower[p, groups, ctx],
    {p, powers}
  ];

  Sort[Values[groups]]
];

PowerCoefficient[expr_, power_, ctx_Association] := Module[
  {x = DiffExp`Symbols`x, p = NormalizeSourcePower[power, ctx], coeff},

  If[DiffExp`Utilities`PChop[expr] === 0,
    Return[0]
  ];

  coeff = Quiet@Check[
    SeriesCoefficient[expr, {x, 0, p}],
    $Failed
  ];

  If[coeff === $Failed || coeff === Indeterminate ||
      !FreeQ[coeff, SeriesCoefficient | Indeterminate | DirectedInfinity | ComplexInfinity],
    0,
    DiffExp`Utilities`PChop[coeff]
  ]
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
   bLeadPow, bMaxLogK, bCoeffs, sourceBases, sigmas,
   sigma, kMax, kMaxInitial, kMaxGrowthAttempts,
   fCoeffs, maxResidual, fParticular, fParticularPieces},

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

  sourceBases = SourceBasePowers[bVec, bMaxLogK, ctx];

  If[sourceBases === {},
    Return[ConstantArray[
      SeriesData[DiffExp`Symbols`x, 0, {}, 0, maxOrd + 1, 1],
      systemSize
    ]]
  ];

  sigmas = NormalizeSourcePower[# + 1, ctx] & /@ sourceBases;

  fParticularPieces = Table[
    sigma = sigmas[[sigmaIndex]];
    bLeadPow = sourceBases[[sigmaIndex]];

    (* Initial K_max: log depth of source + max Jordan block size among resonating eigenvalues.
       MaxLogPowers = blockSize - 1, so we add +1 to get blockSize.
       This ensures the block system at singular n is fully solvable. *)
    kMaxInitial = bMaxLogK + Max[resInfo["MaxLogPowers"]] + 1;
    kMax = kMaxInitial;

    (* Extract beta coefficients for this source sector:
       beta_{n,k,comp} = coeff of x^{sigma-1+n} (Logx)^k in bVec. *)
    bCoeffs = Table[
      PowerCoefficient[
        DiffExp`SeriesOps`LogxCoeff[bVec[[comp]], logk],
        sigma - 1 + n,
        ctx
      ],
      {n, 0, maxOrd}, {logk, 0, bMaxLogK}, {comp, systemSize}
    ];

    If[DiffExp`Utilities`PChop[bCoeffs] === ConstantArray[0, Dimensions[bCoeffs]],
      Nothing,

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

        If[TrueQ[DiffExp`State`$DebugFuchsianizedCheck],
          Print["UNIPART sigma=", sigma, " bMaxLogK=", bMaxLogK,
            " kMax=", kMax, " maxResidual=", ScientificForm[N[maxResidual, 3]],
            " tol=", ScientificForm[N[tolerance, 2]]];
        ];

        (* Accept unless the residual is PROVEN above tolerance: residuals
           of steps fed by accuracy-limited data are significance zeros
           (0``a) whose comparison against the tolerance is undetermined,
           and they must not reject an exact solve. *)
        If[!TrueQ[maxResidual > tolerance],
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

      BuildParticularFromCoeffs[fCoeffs, sigma, kMax, systemSize, maxOrd, wp]
    ]
    ,
    {sigmaIndex, Length[sigmas]}
  ];

  fParticularPieces = DeleteCases[fParticularPieces, Null];

  If[fParticularPieces === {},
    ConstantArray[
      SeriesData[DiffExp`Symbols`x, 0, {}, 0, maxOrd + 1, 1],
      systemSize
    ],
    DiffExp`Utilities`PChop[Total[fParticularPieces]]
  ]
];

(* Extract the M_i coefficients (M(x) = x*A(x) = sum M_i x^i) from the expanded matrix.
   M_0 is the residue, M_1 = A_0, M_2 = A_1, etc.
   Returns {mCoeffs, numMCoeffs} where mCoeffs[[i+1]] = M_i. *)
ExtractMCoefficients[ctx_Association] := Module[
  {maxOrd, mCoeffs},

  maxOrd = ctx["ExpansionOrder"];

  (* M_i = coefficient of x^i in x*A(x) = coefficient of x^{i-1} in A(x)
     So we extract A_{-1}, A_0, A_1, ..., A_{maxOrd-1} *)
  mCoeffs = ExtractAMatCoefficients[ctx, -1, maxOrd - 1];

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

  If[!KeyExistsQ[cache, "GenSingRR"],
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

  If[TrueQ[DiffExp`State`$DebugBlockResidualSeries],
    Print["GENEIG ", ctx["Label"], " eps^", epsord,
      " exponents=", InputForm[N[resInfo["SolutionExponents"], 6]],
      " maxLogPowers=", InputForm[resInfo["MaxLogPowers"]]];
  ];

  (* Compute particular solution using the unified solver *)
  If[DiffExp`Utilities`PChop[bVec] === ConstantArray[0, systemSize],
    fParticular = ConstantArray[
      SeriesData[DiffExp`Symbols`x, 0, {}, 0, maxOrd + 1, 1],
      systemSize
    ];
    ,
    fParticular = Quiet[
      ComputeUnifiedParticular[bVec, resInfo, mCoeffs, numMCoeffs, ctx],
      {Max::nord}
    ];

    (* UseRationalRecurrence is an exclusive recurrence mode. Do not hide
       recurrence failures by using the default Wronskian/Frobenius path. *)
    If[fParticular === $Failed,
      DiffExp`State`LastErrorContext = {ctx, bVec, resInfo, mCoeffs, numMCoeffs, epsord, cache};
      DiffExp`Utilities`ReportError[
        "General singular recurrence failed to construct a particular solution for integral(s) ",
        ctx["Label"],
        " at epsilon order ",
        epsord,
        ". Refusing to fall back to the default Wronskian/Frobenius path."
      ];
    ];
    fParticular = TrimNonFiniteSeriesTails[fParticular, ctx, epsord,
      "particular"];
  ];

  (* General solution: f = particular + sum_i c_i * FMat_column_i *)
  cIndices = Table[Subscript[c, i], {i, systemSize}];
  fGeneral = fParticular + Sum[cIndices[[i]] * FMat[[All, i]], {i, systemSize}];
  fGeneral = DiffExp`SeriesOps`SExpand[fGeneral];
  fGeneral = TrimNonFiniteSeriesTails[fGeneral, ctx, epsord,
    "general"];

  If[TrueQ[DiffExp`State`$DebugFuchsianizedCheck],
    DebugCheckBlockSolution["GENCHECK", ctx, bVec, cIndices, fGeneral, epsord];
  ];

  {cIndices, fGeneral, cache}
];

(* Gated debug self-check shared by the singular solvers: verify each piece
   of a general solution (particular and each fundamental column) against
   the block ODE f' = A f + b around the local origin.  Enable by setting
   DiffExp`State`$DebugFuchsianizedCheck = True. *)
DebugCheckBlockSolution[tag_String, ctx_Association, bVec_, cIndices_,
    fGeneral_, epsord_] := Module[
  {xv, LL, A, toN, ddx, maxAbs, pieces, checkOrd = 8},
  xv = DiffExp`Symbols`x;
  LL = Unique["LL"];
  A = ctx["AMatFactored"];
  If[MissingQ[A], Return[Null, Module]];
  toN[s_SeriesData] := Module[{nmin = s[[4]], den = s[[6]], cf = s[[3]]},
    Sum[cf[[i]] xv^((nmin + i - 1)/den),
      {i, Min[Length[cf], (checkOrd + 3) den - nmin + 1]}] /.
      DiffExp`Symbols`Logx -> LL
  ];
  toN[e_] := e /. DiffExp`Symbols`Logx -> LL;
  ddx[e_] := D[e, xv] + D[e, LL]/xv;
  maxAbs[e_] := Module[{tr, vals},
    tr = Normal[Quiet[Expand[e] + O[xv]^(checkOrd - 2)]];
    vals = Table[
      Max[Flatten[{0, Abs[N[
        (tr /. br /. LL -> Log[xv] /. xv -> SetPrecision[-1/20, 60]),
        20]]}]],
      {br, {
        {DiffExp`Symbols`\[Theta]p -> 1, DiffExp`Symbols`\[Theta]m -> 0},
        {DiffExp`Symbols`\[Theta]p -> 0, DiffExp`Symbols`\[Theta]m -> 1}
      }}
    ];
    Max[vals]
  ];
  pieces = Join[
    {{"particular", fGeneral /. Thread[cIndices -> 0], toN /@ bVec}},
    Table[
      {"hom" <> ToString[i],
       Coefficient[#, cIndices[[i]]] & /@ fGeneral,
       ConstantArray[0, Length[bVec]]},
      {i, Length[cIndices]}
    ]
  ];
  Do[
    Module[{name = piece[[1]], fv = toN /@ piece[[2]], src = piece[[3]], r},
      r = Max[Table[
        maxAbs[ddx[fv[[comp]]] - (A[[comp]] . fv) - src[[comp]]],
        {comp, Length[fv]}
      ]];
      Print[tag, " ", ctx["Label"], " eps^", epsord, " ", name,
        " maxResid=", ScientificForm[N[r, 3]]];
      (* Optional: print the residual series itself to expose WHICH x-powers,
         Logx depths, and theta content escaped the solve.  Enable with
         DiffExp`State`$DebugBlockResidualSeries = True. *)
      If[TrueQ[DiffExp`State`$DebugBlockResidualSeries] && TrueQ[r > 10^-12],
        Do[
          Module[{resid},
            resid = Normal[Quiet[
              Expand[ddx[fv[[comp]]] - (A[[comp]] . fv) - src[[comp]]] +
                O[xv]^(checkOrd - 2)
            ]];
            resid = DiffExp`Utilities`PChop[resid /. LL -> DiffExp`Symbols`Logx];
            If[!TrueQ[resid === 0] && !TrueQ[PossibleZeroQ[resid]],
              Print["RESIDSER ", ctx["Label"], " eps^", epsord, " ", name,
                " comp=", comp, " : ", InputForm[N[Chop[resid, 10^-40], 8]]];
            ];
          ],
          {comp, Length[fv]}
        ];
      ];
    ],
    {piece, pieces}
  ];
];

SolveFuchsianizedSingularRecurrence[ctx_Association, bVec_, epsord_, cacheIn_Association] := Module[
  {cache = cacheIn, data, T, TInv, ctxG, subcache, bVecG,
   cIndices, gGeneral, fGeneral},

  If[!KeyExistsQ[cache, "FuchsianRR"],
    data = BuildFuchsianizedRecurrenceData[ctx];
    If[data === $Failed,
      DiffExp`State`LastErrorContext = {ctx, bVec, epsord};
      DiffExp`Utilities`ReportError[
        "Local Fuchsianization failed for integral(s) ",
        ctx["Label"],
        ". Refusing to fall back to the default Wronskian/Frobenius path."
      ];
    ];

    DiffExp`Utilities`PrintInfo[
      "Using fuchsianized singular recurrence for integrals ",
      ctx["Label"],
      "."
    ][3];

    cache["FuchsianRR"] = data;
  ];

  data = cache["FuchsianRR"];
  T = data["T"];
  TInv = data["TInv"];
  ctxG = data["Ctx"];
  subcache = data["Subcache"];

  bVecG = DiffExp`SeriesOps`SExpand[TInv . bVec];

  {cIndices, gGeneral, subcache} =
    SolveGeneralSingularRecurrence[ctxG, bVecG, epsord, subcache];

  data["Subcache"] = subcache;
  cache["FuchsianRR"] = data;

  fGeneral = DiffExp`SeriesOps`SExpand[T . gGeneral];

  (* Debug self-check: verify each piece of the back-transformed general
     solution against the original-frame block ODE f' = A f + b.  Enable by
     setting DiffExp`State`$DebugFuchsianizedCheck = True. *)
  If[TrueQ[DiffExp`State`$DebugFuchsianizedCheck],
    Print["FUCHSCHECK eps^", epsord,
      " bVec free of thetas: ",
      FreeQ[bVec, DiffExp`Symbols`\[Theta]p | DiffExp`Symbols`\[Theta]m],
      ", free of Logx: ", FreeQ[bVec, DiffExp`Symbols`Logx]];
    Print["FUCHSCHECK eps^", epsord, " bVecG info: ",
      Table[
        Module[{e = bVecG[[i]]},
          Which[
            MatchQ[e, _SeriesData],
            {"ser", e[[4]]/e[[6]], e[[5]]/e[[6]],
              DiffExp`SeriesOps`MaxLogxPower[e]},
            True,
            {Head[e], DiffExp`SeriesOps`MaxLogxPower[e]}
          ]
        ],
        {i, Length[bVecG]}
      ]
    ];
  ];
  If[TrueQ[DiffExp`State`$DebugFuchsianizedCheck],
    DebugCheckBlockSolution["FUCHSCHECK", ctx, bVec, cIndices, fGeneral, epsord];
  ];

  {cIndices, fGeneral, cache}
];
