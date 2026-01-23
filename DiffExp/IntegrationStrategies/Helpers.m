(* IntegrationStrategies/Helpers.m *)
(* Shared helper functions used by multiple integration strategies *)

(* Compute M^(j) matrices via nested differentiation *)
(* Used by VOP and VOPAlt strategies *)
ComputeMatricesMSupj[intind_, line_] := Module[
  {systemSize = Length[intind], MatricesMSupj},

  If[DiffExp`State`FEC["HomogeneousSolve"] === "Expand",
    MatricesMSupj = DiffExp`Utilities`PChop@NestList[
      DiffExp`SeriesOps`SExpand[(DiffExp`SeriesOps`SD[#, DiffExp`Symbols`x] +
        DiffExp`SeriesOps`MatrixMultiplySExpand[#, DiffExp`State`DEqnMatricesExpanded[line][0][[intind, intind]]])] &,
      IdentityMatrix[systemSize], systemSize];
    ,
    MatricesMSupj = NestList[
      Together[(D[#, DiffExp`Symbols`x] + # . DiffExp`State`DEqnMatricesFactored[line][0][[intind, intind]])] &,
      IdentityMatrix[systemSize], systemSize];
  ];

  MatricesMSupj
];

(* Extract nth-order differential equations from MatricesMSupj *)
(* Used by VOP and VOPAlt strategies *)
ComputeNthOrderDiffEqns[MatricesMSupj_, MyN_] := Module[
  {nullSpaceSolutions},

  Table[
    If[DiffExp`State`FEC["HomogeneousSolve"] === "Expand",
      nullSpaceSolutions = MatricesMSupj[[All, ind]] // Transpose //
        NullSpace[#, Method -> "DivisionFreeRowReduction",
          Tolerance -> 10^-DiffExp`State`LinearSolveChopPrecisionVal] &;
      ,
      nullSpaceSolutions = MatricesMSupj[[All, ind]] // Transpose // NullSpace // Together;
    ];
    nullSpaceSolutions = Internal`DeleteTrailingZeros /@ nullSpaceSolutions;
    nullSpaceSolutions = DiffExp`SeriesOps`DiffExpSeries[MinimalBy[nullSpaceSolutions, Length] // First];
    nullSpaceSolutions/(nullSpaceSolutions // Last)
    , {ind, MyN}
  ]
];

(* Compute Wronskian matrix for a set of homogeneous solutions *)
(* Used by VOP and VOPAlt strategies *)
ComputeWronskianMatrix[homogeneousSolutions_] := Module[
  {numSols = Length[homogeneousSolutions]},

  Table[
    DiffExp`SeriesOps`SD[homogeneousSolutions[[mysolind]],
      Sequence @@ ConstantArray[DiffExp`Symbols`x, diffind]]
    , {diffind, 0, numSols - 1}, {mysolind, numSols}
  ]
];

(* Reduce indeterminate constants using vanishing term constraints *)
(* Used by VOP and VOPAlt strategies *)
(* Parameters:
   - fGeneral: current general solution
   - vanishingTerms: terms that must vanish
   - numSolutions: expected number of remaining indeterminates
   Returns: {updatedFGeneral, cIndices} *)
ReduceIndeterminates[fGeneral_, vanishingTerms_, numSolutions_] := Module[
  {currentConstants, IndeterminatesRemaining, seriesOrderOffset, VanishingTermsCurr, OverdeterminedEqns,
   Cmat, Cb, CsPartSol, CsNullVectors, CsGeneralSol, CsReps, c, result},

  IndeterminatesRemaining = currentConstants = DiffExp`Utilities`GetCases[fGeneral, Subscript[c, i_, j_]];
  result = fGeneral;

  If[!DeleteDuplicates[Normal[vanishingTerms]] === {0},
    seriesOrderOffset = -1;

    While[Length[IndeterminatesRemaining] > numSolutions,
      VanishingTermsCurr = Series[vanishingTerms, {DiffExp`Symbols`x, 0, seriesOrderOffset}];

      DiffExp`Utilities`PrintDebug["Reducing number of indeterminate constants, considering terms up to O[x]^", seriesOrderOffset][3];
      OverdeterminedEqns = Flatten[Flatten[DiffExp`SeriesOps`LogxCoeffList /@ VanishingTermsCurr][[All, 3]]];

      If[Length[OverdeterminedEqns] > 0,
        {Cmat, Cb} = {#[[2]], -#[[1]]} &@CoefficientArrays[DeleteCases[OverdeterminedEqns, True], currentConstants];

        CsPartSol = LinearSolve[Cmat, Cb,
          ZeroTest -> (N[DiffExp`Utilities`LSPChop@Expand@Normal[#1], DiffExp`State`LinearSolveChopPrecisionVal] == 0 &)];

        CsNullVectors = Cmat // NullSpace[#, Tolerance -> 10^-DiffExp`State`LinearSolveChopPrecisionVal] &;
        CsGeneralSol = CsPartSol + Sum[CsNullVectors[[ind]] Subscript[c, ind], {ind, Length[CsNullVectors]}] //
          DiffExp`SeriesOps`SExpand;

        CsReps = Thread[currentConstants -> CsGeneralSol] // DiffExp`Utilities`PChop;
        IndeterminatesRemaining = DiffExp`Utilities`GetCases[CsGeneralSol, Subscript[c, __]];
      ];

      seriesOrderOffset += 1;
    ];

    If[Length[IndeterminatesRemaining] < numSolutions,
      DiffExp`State`LastErrorContext = {Cmat, Cb, CsGeneralSol, result};
      DiffExp`Utilities`ReportError["There was an error reducing indeterminates. Try decreasing ChopPrecision or increasing WorkingPrecision."];
    ];

    result = result /. CsReps;
    ,
    (* No constraints needed - just rename constants *)
    result = result /. Table[currentConstants[[ind]] -> Subscript[c, ind], {ind, Length@currentConstants}];
  ];

  result
];

(* Compute general solution matrix GMat from FMat, FMatInv, and inhomogeneous term *)
(* Used by SolveDefault and VOPAlt strategies *)
ComputeGMat[FMat_, FMatInv_, bVec_, intind_] := Module[
  {BMat, cIndices, GMat, c},

  BMat = 1/Length@intind Table[bVec, {iind, intind // Length}] // Transpose;
  cIndices = Table[Subscript[c, i], {i, intind // Length}];

  GMat = DiffExp`SeriesOps`MatrixMultiplySExpand[
    FMat,
    DiffExp`Integration`DiffExpIntegrate[DiffExp`SeriesOps`MatrixMultiplySExpand[FMatInv, BMat]] +
      DiagonalMatrix[cIndices]
  ] // DiffExp`Utilities`PChop;

  {cIndices, GMat}
];
