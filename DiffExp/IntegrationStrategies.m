(* ::Package:: *)

(* DiffExp Integration Strategies Subpackage *)
(* This package provides different integration strategy implementations *)
(* for solving systems of differential equations *)

BeginPackage["DiffExp`IntegrationStrategies`", {
  "DiffExp`Symbols`",
  "DiffExp`State`",
  "DiffExp`Utilities`",
  "DiffExp`SeriesOps`",
  "DiffExp`Integration`",
  "DiffExp`Frobenius`",
  "DiffExp`Wronskian`"
}];

(* Strategy functions *)
SolveSimple::usage = "SolveSimple[intind, bVec, line, epsord] solves a simple integration case where the integral has no homogeneous components.";
SolveDefault::usage = "SolveDefault[intind, bVec, line, epsord, BufferedData] solves using the default Frobenius/Wronskian strategy.";
SolveVOP::usage = "SolveVOP[intind, bVec, line, epsord, BufferedData] solves using the variation of parameters strategy.";
SolveVOPAlt::usage = "SolveVOPAlt[intind, bVec, line, epsord, BufferedData] solves using the alternative variation of parameters strategy.";
DispatchStrategy::usage = "DispatchStrategy[intind, bVec, line, epsord, BufferedData] dispatches to the appropriate integration strategy based on configuration.";

Begin["`Private`"];

(* ============================================================================ *)
(* SHARED HELPER FUNCTIONS                                                      *)
(* ============================================================================ *)

(* Compute M^(j) matrices via nested differentiation *)
(* Used by VOP and VOPAlt strategies *)
ComputeMatricesMSupj[intind_, line_] := Module[
  {MyN = Length[intind], MatricesMSupj},

  If[DiffExp`State`FEC["HomogeneousSolve"] === "Expand",
    MatricesMSupj = DiffExp`Utilities`PChop@NestList[
      DiffExp`SeriesOps`SExpand[(DiffExp`SeriesOps`SD[#, DiffExp`Symbols`x] +
        DiffExp`SeriesOps`MatrixMultiplySExpand[#, DiffExp`State`DEqnMatricesExpanded[line][0][[intind, intind]]])] &,
      IdentityMatrix[MyN], MyN];
    ,
    MatricesMSupj = NestList[
      Together[(D[#, DiffExp`Symbols`x] + # . DiffExp`State`DEqnMatricesFactored[line][0][[intind, intind]])] &,
      IdentityMatrix[MyN], MyN];
  ];

  MatricesMSupj
];

(* Extract nth-order differential equations from MatricesMSupj *)
(* Used by VOP and VOPAlt strategies *)
ComputeNthOrderDiffEqns[MatricesMSupj_, MyN_] := Module[
  {TmpSols},

  Table[
    If[DiffExp`State`FEC["HomogeneousSolve"] === "Expand",
      TmpSols = MatricesMSupj[[All, ind]] // Transpose //
        NullSpace[#, Method -> "DivisionFreeRowReduction",
          Tolerance -> 10^-DiffExp`State`LinearSolveChopPrecisionVal] &;
      ,
      TmpSols = MatricesMSupj[[All, ind]] // Transpose // NullSpace // Together;
    ];
    TmpSols = Internal`DeleteTrailingZeros /@ TmpSols;
    TmpSols = DiffExp`SeriesOps`DiffExpSeries[MinimalBy[TmpSols, Length] // First];
    TmpSols/(TmpSols // Last)
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
  {csCurr, IndeterminatesRemaining, xAdd, VanishingTermsCurr, OverdeterminedEqns,
   Cmat, Cb, CsPartSol, CsNullVectors, CsGeneralSol, CsReps, c, result},

  IndeterminatesRemaining = csCurr = DiffExp`Utilities`GetCases[fGeneral, Subscript[c, i_, j_]];
  result = fGeneral;

  If[!DeleteDuplicates[Normal[vanishingTerms]] === {0},
    xAdd = -1;

    While[Length[IndeterminatesRemaining] > numSolutions,
      VanishingTermsCurr = Series[vanishingTerms, {DiffExp`Symbols`x, 0, xAdd}];

      DiffExp`Utilities`PrintDebug["Reducing number of indeterminate constants, considering terms up to O[x]^", xAdd][3];
      OverdeterminedEqns = Flatten[Flatten[DiffExp`SeriesOps`LogxCoeffList /@ VanishingTermsCurr][[All, 3]]];

      If[Length[OverdeterminedEqns] > 0,
        {Cmat, Cb} = {#[[2]], -#[[1]]} &@CoefficientArrays[DeleteCases[OverdeterminedEqns, True], csCurr];

        CsPartSol = LinearSolve[Cmat, Cb,
          ZeroTest -> (N[DiffExp`Utilities`LSPChop@Expand@Normal[#1], DiffExp`State`LinearSolveChopPrecisionVal] == 0 &)];

        CsNullVectors = Cmat // NullSpace[#, Tolerance -> 10^-DiffExp`State`LinearSolveChopPrecisionVal] &;
        CsGeneralSol = CsPartSol + Sum[CsNullVectors[[ind]] Subscript[c, ind], {ind, Length[CsNullVectors]}] //
          DiffExp`SeriesOps`SExpand;

        CsReps = Thread[csCurr -> CsGeneralSol] // DiffExp`Utilities`PChop;
        IndeterminatesRemaining = DiffExp`Utilities`GetCases[CsGeneralSol, Subscript[c, __]];
      ];

      xAdd += 1;
    ];

    If[Length[IndeterminatesRemaining] < numSolutions,
      DiffExp`State`LastErrorContext = {Cmat, Cb, CsGeneralSol, result};
      DiffExp`Utilities`ReportError["There was an error reducing indeterminates. Try decreasing ChopPrecision or increasing WorkingPrecision."];
    ];

    result = result /. CsReps;
    ,
    (* No constraints needed - just rename constants *)
    result = result /. Table[csCurr[[ind]] -> Subscript[c, ind], {ind, Length@csCurr}];
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

(* ============================================================================ *)
(* INTEGRATION STRATEGIES                                                       *)
(* ============================================================================ *)

(* Simple integration strategy *)
(* Used when there's a single integral without homogeneous components *)
SolveSimple[intind_, bVec_, line_, epsord_] := Module[
  {cIndices, fGeneral, c},

  If[epsord === 0,
    DiffExp`State`BenchmarkData["Segments"][line // N]["HomogeneousSolveAllPreprocessing"]["Integrals"][intind] = 0;
  ];

  cIndices = {Subscript[c, 1]};

  If[DiffExp`Utilities`PChop[bVec] === {0},
    fGeneral = {Subscript[c, 1] + O[DiffExp`Symbols`x]^(DiffExp`State`ExpansionOrderVal + 1)};
    ,
    fGeneral = {DiffExp`Integration`DiffExpIntegrate[bVec[[1]], DiffExp`Symbols`x] + Subscript[c, 1] + O[DiffExp`Symbols`x]^(DiffExp`State`ExpansionOrderVal + 1)};
  ];

  {cIndices, fGeneral}
];

(* Default integration strategy *)
(* Uses Frobenius solutions and Wronskian computation *)
SolveDefault[intind_, bVec_, line_, epsord_, BufferedDataIn_] := Module[
  {ll, HomogeneousEquation, MtildeMat, NMat, Solns, Wronsk, WronskInv, FMat, FMatInv,
   GMat, BMat, cIndices, fGeneral, CrossC, CurrInvWronskSolver,
   HomogeneousEquation2, MtildeMat2, NMat2, Solns2, Wronsk2, WronskInvPrime, WWinvprimeprod, MtildeInv, c,
   BufferedData = BufferedDataIn},

  If[epsord === 0 && !KeyExistsQ[BufferedData, intind],
    DiffExp`State`BenchmarkData["Segments"][line // N]["HomogeneousSolveAllPreprocessing"]["Integrals"][intind] = AbsoluteTime[];

    If[Length[intind] > 1,
      DiffExp`Utilities`PrintInfo["Combining differential equations: ", intind -> intind[[1]], "."][3];
    ];
    ll = intind // Length;
    {HomogeneousEquation, MtildeMat} = DiffExp`Wronskian`CombineDifferentialEquationsHomogeneous[
      If[DiffExp`State`FEC["HomogeneousSolve"] === "Expand",
        DiffExp`State`DEqnMatricesExpanded[line][0][[intind, intind]],
        DiffExp`State`DEqnMatricesFactored[line][0][[intind, intind]]
      ]
      , 1];

    DiffExp`Utilities`PrintDebug["Found homogeneous differential equation: ", HomogeneousEquation + O[DiffExp`Symbols`x]^4 // DiffExp`SeriesOps`SN][3];

    NMat = Table[
      If[ind < ll, UnitVector[ll, ind + 1], -HomogeneousEquation[[Range@ll]]]
      , {ind, ll}];

    DiffExp`Utilities`PrintDebug["Deriving solutions..."][3];
    Solns = DiffExp`Frobenius`FrobeniusSolutions[HomogeneousEquation];
    DiffExp`Utilities`PrintDebug["All homogeneous solutions found..."][3];
    DiffExp`Utilities`PrintDebug["Deriving period matrix..."][3];
    Wronsk = Table[
      DiffExp`SeriesOps`SD[Solns[[iind]], Sequence @@ DiffExp`Utilities`CA[DiffExp`Symbols`x, jind - 1]]
      , {jind, ll}, {iind, ll}
    ] // DiffExp`SeriesOps`SExpand;

    DiffExp`Utilities`PrintDebug["Got Wronskian..."][3];
    DiffExp`Utilities`PrintDebug["Inverting Wronskian..."][3];

    If[(DiffExp`State`FEC["InvWronskSolver"] === "Auto"),
      If[(DiffExp`Utilities`DependsQ[Wronsk, DiffExp`Symbols`Logx]),
        CurrInvWronskSolver = "Frobenius";,
        CurrInvWronskSolver = "Inverse";
      ];
      ,
      CurrInvWronskSolver = DiffExp`State`FEC["InvWronskSolver"];
    ];

    If[CurrInvWronskSolver === "Frobenius",
      DiffExp`Utilities`PrintDebug["Determining inverse Wronskian using the Frobenius method."][1];
      {HomogeneousEquation2, MtildeMat2} = DiffExp`Wronskian`CombineDifferentialEquationsHomogeneous[
        -NMat // Transpose, 1];

      NMat2 = Table[
        If[ind < ll, UnitVector[ll, ind + 1], -HomogeneousEquation2[[Range@ll]]]
        , {ind, ll}];

      Solns2 = DiffExp`Frobenius`FrobeniusSolutions[HomogeneousEquation2];

      DiffExp`Utilities`PrintDebug["Deriving inverse Wronskian."][3];

      Wronsk2 = Table[
        DiffExp`SeriesOps`SD[Solns2[[iind]], Sequence @@ DiffExp`Utilities`CA[DiffExp`Symbols`x, jind - 1]]
        , {jind, ll}, {iind, ll}
      ] // DiffExp`SeriesOps`SExpand;

      (* Cross-checking Wronskians *)
      If[MemberQ[DiffExp`State`CurrCrosscheckFlags, "Wronskians"] === True,
        DiffExp`Utilities`PrintDebug["Cross-checking Wronskians."][1];
        CrossC = (DiffExp`SeriesOps`SD[Wronsk, DiffExp`Symbols`x] - DiffExp`SeriesOps`MatrixMultiplySExpand[NMat, Wronsk] // (DiffExp`Utilities`CPChop@*DiffExp`SeriesOps`SExpand)) + O[DiffExp`Symbols`x]^DiffExp`State`ICrossCheckVerifyResultOrder;
        DiffExp`Utilities`PrintDebug["Found: ", CrossC + O[DiffExp`Symbols`x]^DiffExp`State`ICrossCheckPrintResultOrder // DiffExp`SeriesOps`SN][3];
        If[
          !(SameQ @@ Append[CrossC // Normal // Flatten, 0])
          ,
          DiffExp`Utilities`ReportError["Cross-check failed."];
        ];

        CrossC = (DiffExp`SeriesOps`SD[Wronsk2, DiffExp`Symbols`x] - DiffExp`SeriesOps`MatrixMultiplySExpand[NMat2, Wronsk2] // (DiffExp`Utilities`PChop@*DiffExp`SeriesOps`SExpand)) + O[DiffExp`Symbols`x]^DiffExp`State`ICrossCheckVerifyResultOrder;
        DiffExp`Utilities`PrintDebug["Found: ", CrossC + O[DiffExp`Symbols`x]^DiffExp`State`ICrossCheckPrintResultOrder // DiffExp`SeriesOps`SN][1];
        If[
          !(SameQ @@ Append[CrossC // Normal // Flatten, 0])
          ,
          DiffExp`Utilities`ReportError["Cross-check failed."];
        ]
      ];

      WronskInvPrime = Inverse[MtildeMat2] . Wronsk2 // Transpose // DiffExp`SeriesOps`SExpand;
      DiffExp`Utilities`PrintDebug["Deriving inverse Wronskian.."][3];
      WWinvprimeprod = DiffExp`Utilities`PChop@Normal@DiffExp`SeriesOps`MatrixMultiplySExpand[WronskInvPrime, Wronsk];
      DiffExp`Utilities`PrintDebug["Deriving inverse Wronskian..."][3];
      If[MemberQ[DiffExp`State`CurrCrosscheckFlags, "WronskInv"] === True,
        If[DiffExp`Utilities`DependsQ[(WWinvprimeprod // DiffExp`Utilities`CPChop) + O[DiffExp`Symbols`x]^DiffExp`State`ICrossCheckVerifyResultOrder // Normal, DiffExp`Symbols`x],
          DiffExp`State`LastErrorContext = {WronskInvPrime, Wronsk, WWinvprimeprod};
          DiffExp`Utilities`ReportError["Warning, product of Wronskian inverse times Wronskian does not match identity. Try increasing \"ChopPrecision\" or \"WorkingPrecision\", or try decreasing \"ExpansionOrder\"."][1];
        ];
      ];
      DiffExp`State`LastErrorContext = {WWinvprimeprod};
      WWinvprimeprod = WWinvprimeprod /. DiffExp`Symbols`x^(k_ : 1) :> 0 /. DiffExp`Symbols`Logx -> 0;
      WronskInv = DiffExp`Utilities`PChop@DiffExp`SeriesOps`MatrixMultiplySExpand[(WWinvprimeprod // Inverse), WronskInvPrime];
      DiffExp`Utilities`PrintDebug["Done."][3];
    ];

    If[CurrInvWronskSolver === "Inverse",
      WronskInv = Inverse[Wronsk, Method -> "DivisionFreeRowReduction"];
    ];

    If[CurrInvWronskSolver === "InverseLogx",
      WronskInv = DiffExp`Wronskian`MatrixLogxInverse[Wronsk];
    ];

    DiffExp`Utilities`PrintDebug["Inverting MTilde... "][3];

    Check[
      If[DiffExp`State`FEC["HomogeneousSolve"] === "Expand",
        MtildeInv = Inverse[MtildeMat, Method -> "DivisionFreeRowReduction", ZeroTest -> (Normal[#] == 0 &)];,
        MtildeInv = DiffExp`SeriesOps`DiffExpSeries[Inverse[MtildeMat] // Together];
      ];
      FMat = DiffExp`Utilities`PChop@(DiffExp`SeriesOps`MatrixMultiplySExpand[MtildeInv, Wronsk] // DiffExp`SeriesOps`SExpand);
      ,
      DiffExp`Utilities`PrintWarning["Encountered numerical instability while inverting Mtilde. Turning off DivisionFreeRowReduction and trying again.."];
      FMat = DiffExp`Utilities`PChop@(DiffExp`SeriesOps`MatrixMultiplySExpand[Inverse[MtildeMat], Wronsk] // DiffExp`SeriesOps`SExpand);
    ];

    FMatInv = DiffExp`Utilities`PChop@(DiffExp`SeriesOps`MatrixMultiplySExpand[WronskInv, MtildeMat]);

    (* Cross-checking FMat *)
    If[MemberQ[DiffExp`State`CurrCrosscheckFlags, "PeriodMatrix"] === True,
      DiffExp`Utilities`PrintDebug["Cross-checking period matrix.."][1];
      CrossC = DiffExp`SeriesOps`SD[FMat, DiffExp`Symbols`x] - DiffExp`SeriesOps`MatrixMultiplySExpand[DiffExp`State`DEqnMatricesExpanded[line][0][[intind, intind]], FMat + O[DiffExp`Symbols`x]^DiffExp`State`ICrossCheckVerifyResultOrder] // (DiffExp`Utilities`CPChop@*DiffExp`SeriesOps`SExpand);
      DiffExp`State`LastErrorContext = {FMat, DiffExp`State`DEqnMatricesExpanded[line][0][[intind, intind]], CrossC, MtildeMat, Wronsk};
      DiffExp`Utilities`PrintDebug["Found: ", CrossC + O[DiffExp`Symbols`x]^DiffExp`State`ICrossCheckPrintResultOrder // DiffExp`SeriesOps`SN][3];
      If[
        !(SameQ @@ Append[CrossC // Normal // Flatten, 0])
        ,
        DiffExp`Utilities`ReportError["Cross-check of period matrix failed"];
      ];
    ];

    DiffExp`Utilities`PrintDebug["Period matrix derived."][3];

    BufferedData[intind] = {FMat, FMatInv};
    DiffExp`State`BenchmarkData["Segments"][line // N]["HomogeneousSolveAllPreprocessing"]["Integrals"][intind] = AbsoluteTime[] - DiffExp`State`BenchmarkData["Segments"][line // N]["HomogeneousSolveAllPreprocessing"]["Integrals"][intind];
  ];

  {FMat, FMatInv} = BufferedData[intind];

  (* Use shared helper for GMat computation *)
  DiffExp`Utilities`PrintDebug["Setting up general solution."][3];
  {cIndices, GMat} = ComputeGMat[FMat, FMatInv, bVec, intind];
  BMat = 1/Length@intind Table[bVec, {iind, intind // Length}] // Transpose;

  If[MemberQ[DiffExp`State`CurrCrosscheckFlags, "GeneralSolutionMatrix"] === True,
    DiffExp`Utilities`PrintDebug["Cross-checking GMat with differential equations."][1];
    CrossC = DiffExp`SeriesOps`SD[GMat, DiffExp`Symbols`x] - DiffExp`SeriesOps`MatrixMultiplySExpand[DiffExp`State`DEqnMatricesExpanded[line][0][[intind, intind]], GMat + O[DiffExp`Symbols`x]^DiffExp`State`ICrossCheckVerifyResultOrder] - BMat // (DiffExp`Utilities`CPChop@*DiffExp`SeriesOps`SExpand);
    DiffExp`State`LastErrorContext = {FMat, GMat, FMatInv, BMat, DiffExp`State`DEqnMatricesExpanded[line][0][[intind, intind]], CrossC, MtildeMat, Wronsk};
    DiffExp`Utilities`PrintDebug["Found: ", CrossC + O[DiffExp`Symbols`x]^DiffExp`State`ICrossCheckPrintResultOrder // DiffExp`SeriesOps`SN][3];
    If[
      !(SameQ @@ Append[CrossC + O[DiffExp`Symbols`x]^DiffExp`State`ICrossCheckVerifyResultOrder // Normal // Flatten, 0])
      ,
      DiffExp`Utilities`ReportError["Cross-check of solution matrix failed"];
    ];
  ];

  fGeneral = Total[GMat // Transpose];

  {cIndices, fGeneral, BufferedData}
];

(* VOP integration strategy *)
(* Uses variation of parameters method *)
SolveVOP[intind_, bVec_, line_, epsord_, BufferedDataIn_] := Module[
  {MyN, MatricesMSupj, nthOrderDifferentialEquations, nthOrderDifferentialEquationsSolutions,
   SolveFrom, MtildeMatrix, HighestOrderDifferentialEquationPosition,
   MyHomogeneousSolutions, MyNumberOfSolutions, MyInhomogeneousTerm,
   bSupjVecs, fGeneral, cIndices, VanishingTerms,
   DerivativeVec, BDerVec, MtildeInv, MyWi, MyWronsk, MyWronskDetInv, c, PlaceHolder,
   BufferedData = BufferedDataIn},

  DiffExp`Utilities`PrintDebug["Solving by variation of parameters. ", intind][3];

  MyN = intind // Length;

  (* Initialize MyWronsk and MyWronskDetInv as local associations *)
  MyWronsk = Association[{}];
  MyWronskDetInv = Association[{}];

  If[epsord === 0 && !KeyExistsQ[BufferedData, intind],
    DiffExp`State`BenchmarkData["Segments"][line // N]["HomogeneousSolveAllPreprocessing"]["Integrals"][intind] = AbsoluteTime[];

    (* Use shared helper functions for common computations *)
    MatricesMSupj = ComputeMatricesMSupj[intind, line];
    nthOrderDifferentialEquations = ComputeNthOrderDiffEqns[MatricesMSupj, MyN];

    If[DiffExp`State`FEC["HomogeneousSolve"] === "Expand",
      MatricesMSupj = DiffExp`SeriesOps`DiffExpSeries[MatricesMSupj];
    ];

    HighestOrderDifferentialEquationPosition = FirstPosition[Length /@ nthOrderDifferentialEquations, Max[Length /@ nthOrderDifferentialEquations]][[1]];
    If[Length[nthOrderDifferentialEquations[[HighestOrderDifferentialEquationPosition]]] < MyN + 1,
      SolveFrom = Range@Length@nthOrderDifferentialEquations;
      nthOrderDifferentialEquationsSolutions = DiffExp`Frobenius`FrobeniusSolutions /@ nthOrderDifferentialEquations;
      If[MyN == 1,
        If[DiffExp`State`FEC["HomogeneousSolve"] === "Expand",
          MtildeMatrix = DiffExp`Utilities`PChop@DiffExp`SeriesOps`DiffExpSeries[MatricesMSupj[[All, SolveFrom[[1]]]][[Range[MyN], Range[MyN]]]];
          ,
          MtildeMatrix = MatricesMSupj[[All, SolveFrom[[1]]]][[Range[MyN], Range[MyN]]];
        ];
        ,
        MtildeMatrix = {};
      ];
      ,
      SolveFrom = {HighestOrderDifferentialEquationPosition};

      If[DiffExp`State`FEC["HomogeneousSolve"] === "Expand",
        MtildeMatrix = DiffExp`Utilities`PChop@DiffExp`SeriesOps`DiffExpSeries[MatricesMSupj[[All, SolveFrom[[1]]]][[Range[MyN], Range[MyN]]]];
        ,
        MtildeMatrix = MatricesMSupj[[All, SolveFrom[[1]]]][[Range[MyN], Range[MyN]]];
      ];
      nthOrderDifferentialEquationsSolutions = ConstantArray[Null, Length@nthOrderDifferentialEquations];
      nthOrderDifferentialEquationsSolutions[[HighestOrderDifferentialEquationPosition]] = DiffExp`Frobenius`FrobeniusSolutions@nthOrderDifferentialEquations[[HighestOrderDifferentialEquationPosition]];
    ];

    DiffExp`Utilities`PrintDebug["Determining Wronskian matrices."][3];

    Table[
      MyHomogeneousSolutions = nthOrderDifferentialEquationsSolutions[[myintind]];

      (* Use shared helper for Wronskian matrix computation *)
      MyWronsk[{intind, myintind}] = ComputeWronskianMatrix[MyHomogeneousSolutions];

      (* VOP-specific: compute inverse Wronskian determinant *)
      DiffExp`Utilities`PrintDebug["Determining 1/WronskianDet: ", myintind][3];
      MyWronskDetInv[{intind, myintind}] = DiffExp`SeriesOps`DiffExpSeries[1/(MyWronsk[{intind, myintind}] // Det // DiffExp`SeriesOps`SExpand), DiffExp`State`ExpansionOrderVal];
      DiffExp`Utilities`PrintDebug["Done determining 1/WronskianDet: ", myintind][3];
      , {myintind, SolveFrom}
    ];

    BufferedData[intind] = {MatricesMSupj, nthOrderDifferentialEquations, nthOrderDifferentialEquationsSolutions, SolveFrom, MtildeMatrix, MyWronsk, MyWronskDetInv};

    DiffExp`State`BenchmarkData["Segments"][line // N]["HomogeneousSolveAllPreprocessing"]["Integrals"][intind] = AbsoluteTime[] - DiffExp`State`BenchmarkData["Segments"][line // N]["HomogeneousSolveAllPreprocessing"]["Integrals"][intind];
  ];

  {MatricesMSupj, nthOrderDifferentialEquations, nthOrderDifferentialEquationsSolutions, SolveFrom, MtildeMatrix, MyWronsk, MyWronskDetInv} = BufferedData[intind];

  DiffExp`Utilities`PrintDebug["Determining b's.."][3];
  bSupjVecs = FoldList[DiffExp`SeriesOps`SExpand[#2 . bVec + DiffExp`SeriesOps`SD[#1, DiffExp`Symbols`x]] &, ConstantArray[0, MyN], Delete[MatricesMSupj, -1]];

  DiffExp`Utilities`PrintDebug["Setting up general solutions. ", intind][3];
  fGeneral = Table[
    MyHomogeneousSolutions = nthOrderDifferentialEquationsSolutions[[myintind]];
    MyNumberOfSolutions = MyHomogeneousSolutions // Length;
    MyInhomogeneousTerm = nthOrderDifferentialEquations[[myintind]] . bSupjVecs[[All, myintind]][[Range@(MyNumberOfSolutions + 1)]];

    DiffExp`SeriesOps`SExpand@Sum[
      MyWi = MyWronsk[{intind, myintind}];
      MyWi[[All, ind]] = Append[ConstantArray[0, MyNumberOfSolutions - 1], PlaceHolder];

      (MyHomogeneousSolutions[[ind]] (Subscript[c, myintind, ind] + DiffExp`Integration`DiffExpIntegrate[
        MyWronskDetInv[{intind, myintind}] ((MyWi // Det // DiffExp`SeriesOps`SExpand) /. PlaceHolder -> MyInhomogeneousTerm // DiffExp`SeriesOps`SExpand) // DiffExp`SeriesOps`SExpand
      ]) // DiffExp`SeriesOps`SExpand) // DiffExp`SeriesOps`SExpand

      , {ind, MyHomogeneousSolutions // Length}]

    , {myintind, SolveFrom}
  ];

  DiffExp`Utilities`PrintDebug["Done. ", intind][3];

  If[Length[SolveFrom] > 1,
    VanishingTerms = DiffExp`SeriesOps`SExpand@(DiffExp`SeriesOps`SD[fGeneral, DiffExp`Symbols`x] - DiffExp`State`DEqnMatricesExpanded[line][0][[intind, intind]] . fGeneral - bVec);

    IndeterminatesRemaining = csCurr = DiffExp`Utilities`GetCases[fGeneral, Subscript[c, i_, j_]];

    If[!DeleteDuplicates[Normal[VanishingTerms]] === {0},
      xAdd = -1;

      While[Length[IndeterminatesRemaining] > Length[fGeneral],
        VanishingTermsCurr = Series[VanishingTerms, {DiffExp`Symbols`x, 0, xAdd}];

        DiffExp`Utilities`PrintDebug["Reducing number of indeterminate constants, considering terms up to O[x]^", xAdd][3];
        OverdeterminedEqns = Flatten[Flatten[DiffExp`SeriesOps`LogxCoeffList /@ VanishingTermsCurr][[All, 3]]];

        If[Length[OverdeterminedEqns] > 0,

          {Cmat, Cb} = {#[[2]], -#[[1]]} &@CoefficientArrays[DeleteCases[OverdeterminedEqns, True], csCurr];

          CsPartSol = LinearSolve[Cmat, Cb, ZeroTest -> (N[DiffExp`Utilities`LSPChop@Expand@Normal[#1], DiffExp`State`LinearSolveChopPrecisionVal] == 0 &)];

          CsNullVectors = Cmat // NullSpace[#, Tolerance -> 10^-DiffExp`State`LinearSolveChopPrecisionVal] &;
          CsGeneralSol = CsPartSol + Sum[CsNullVectors[[ind]] Subscript[c, ind], {ind, Length[CsNullVectors]}] // DiffExp`SeriesOps`SExpand;

          CsReps = Thread[csCurr -> CsGeneralSol] // DiffExp`Utilities`PChop;

          IndeterminatesRemaining = DiffExp`Utilities`GetCases[CsGeneralSol, Subscript[c, __]];

        ];

        xAdd += 1;
      ];

      If[Length[IndeterminatesRemaining] < Length[fGeneral],
        DiffExp`State`LastErrorContext = {Cmat, Cb, CsGeneralSol};
        DiffExp`Utilities`ReportError["There was an error in obtaining the solutions for integrals ", intind, " at order ", epsord, " using the \"VOP\" strategy. Try decreasing the value of the option ChopPrecision or increasing the WorkingPrecision."];
      ];

      fGeneral = fGeneral /. CsReps;

      ,

      fGeneral = fGeneral /. Table[csCurr[[ind]] -> Subscript[c, ind], {ind, Length@csCurr}];
    ];

    ,

    (* We can currently put fGeneral[[1]], as SolveFrom is either one element, or all elements. *)
    DerivativeVec = Table[DiffExp`SeriesOps`SD[fGeneral[[1]], Sequence @@ DiffExp`Utilities`CA[DiffExp`Symbols`x, jind - 1]], {jind, intind // Length}] // DiffExp`SeriesOps`SExpand;
    (* But for BDerVec we should choose SolveFrom[[1]] *)
    BDerVec = bSupjVecs[[All, SolveFrom[[1]]]][[Range@MyN]];

    Check[
      If[DiffExp`State`FEC["HomogeneousSolve"] === "Expand",
        MtildeInv = Inverse[MtildeMatrix, Method -> "DivisionFreeRowReduction", ZeroTest -> (Normal[#] == 0 &)];,
        MtildeInv = DiffExp`SeriesOps`DiffExpSeries[Inverse[MtildeMatrix] // Together];
      ];
      fGeneral = DiffExp`Utilities`PChop@(MtildeInv . (DerivativeVec - BDerVec) // DiffExp`SeriesOps`SExpand);
      ,
      DiffExp`Utilities`PrintWarning["Encountered numerical instability while inverting Mtilde. Turning off DivisionFreeRowReduction and trying again.."];
      fGeneral = DiffExp`Utilities`PChop@(Inverse[MtildeMatrix] . (DerivativeVec - BDerVec) // DiffExp`SeriesOps`SExpand);
    ];

    fGeneral = fGeneral /. Thread[DiffExp`Utilities`GetCases[fGeneral, Subscript[c, i_, j_]] -> Table[Subscript[c, ind], {ind, MyN}]];

  ];

  cIndices = DiffExp`Utilities`GetCases[fGeneral, Subscript[c, i_]];

  {cIndices, fGeneral, BufferedData}
];

(* VOPAlt integration strategy *)
(* Alternative variation of parameters method *)
SolveVOPAlt[intind_, bVec_, line_, epsord_, BufferedDataIn_] := Module[
  {MyN, MatricesMSupj, nthOrderDifferentialEquations, nthOrderDifferentialEquationsSolutions,
   SolveFrom, MtildeMatrix, HighestOrderDifferentialEquationPosition,
   MyHomogeneousSolutions, fGeneral, cIndices, VanishingTerms, IndeterminatesRemaining, csCurr,
   xAdd, VanishingTermsCurr, OverdeterminedEqns, Cmat, Cb, CsPartSol, CsNullVectors,
   CsGeneralSol, CsReps, FMat, FMatInv, GMat, BMat, CrossC, MyWronsk, c,
   BufferedData = BufferedDataIn},

  DiffExp`Utilities`PrintDebug["Solving by alternative strategy. ", intind][3];

  MyN = intind // Length;

  If[epsord === 0 && !KeyExistsQ[BufferedData, intind],
    DiffExp`State`BenchmarkData["Segments"][line // N]["HomogeneousSolveAllPreprocessing"]["Integrals"][intind] = AbsoluteTime[];

    (* Use shared helper functions for common computations *)
    MatricesMSupj = ComputeMatricesMSupj[intind, line];
    nthOrderDifferentialEquations = ComputeNthOrderDiffEqns[MatricesMSupj, MyN];

    If[DiffExp`State`FEC["HomogeneousSolve"] === "Expand",
      MatricesMSupj = DiffExp`SeriesOps`DiffExpSeries[MatricesMSupj];
    ];

    HighestOrderDifferentialEquationPosition = FirstPosition[Length /@ nthOrderDifferentialEquations, Max[Length /@ nthOrderDifferentialEquations]][[1]];
    SolveFrom = Range@Length@nthOrderDifferentialEquations;
    nthOrderDifferentialEquationsSolutions = DiffExp`Frobenius`FrobeniusSolutions /@ nthOrderDifferentialEquations;
    If[MyN == 1,
      If[DiffExp`State`FEC["HomogeneousSolve"] === "Expand",
        MtildeMatrix = DiffExp`Utilities`PChop@DiffExp`SeriesOps`DiffExpSeries[MatricesMSupj[[All, SolveFrom[[1]]]][[Range[MyN], Range[MyN]]]];
        ,
        MtildeMatrix = MatricesMSupj[[All, SolveFrom[[1]]]][[Range[MyN], Range[MyN]]];
      ];
      ,
      MtildeMatrix = {};
    ];
    DiffExp`Utilities`PrintDebug["Determining Wronskian matrices."][3];

    (* Initialize MyWronsk as local association *)
    MyWronsk = Association[{}];

    (* Use shared helper for Wronskian matrix computation *)
    Table[
      MyHomogeneousSolutions = nthOrderDifferentialEquationsSolutions[[myintind]];
      MyWronsk[{intind, myintind}] = ComputeWronskianMatrix[MyHomogeneousSolutions];
      , {myintind, SolveFrom}
    ];

    BufferedData[intind] = {MatricesMSupj, nthOrderDifferentialEquations, nthOrderDifferentialEquationsSolutions, SolveFrom, MtildeMatrix};

    DiffExp`State`BenchmarkData["Segments"][line // N]["HomogeneousSolveAllPreprocessing"]["Integrals"][intind] = AbsoluteTime[] - DiffExp`State`BenchmarkData["Segments"][line // N]["HomogeneousSolveAllPreprocessing"]["Integrals"][intind];


    DiffExp`Utilities`PrintDebug["Combining homogeneous solutions. ", intind][3];
    fGeneral = Table[
      MyHomogeneousSolutions = nthOrderDifferentialEquationsSolutions[[myintind]];

      DiffExp`SeriesOps`SExpand@Sum[MyHomogeneousSolutions[[ind]] (Subscript[c, myintind, ind]) // DiffExp`SeriesOps`SExpand, {ind, MyHomogeneousSolutions // Length}]

      , {myintind, SolveFrom}
    ];

    DiffExp`Utilities`PrintDebug["Done. ", intind][3];

    VanishingTerms = DiffExp`SeriesOps`SExpand@(DiffExp`SeriesOps`SD[fGeneral, DiffExp`Symbols`x] - DiffExp`State`DEqnMatricesExpanded[line][0][[intind, intind]] . fGeneral);

    IndeterminatesRemaining = csCurr = DiffExp`Utilities`GetCases[fGeneral, Subscript[c, i_, j_]];

    If[!DeleteDuplicates[Normal[VanishingTerms]] === {0},
      xAdd = -1;

      While[Length[IndeterminatesRemaining] > Length[fGeneral],
        VanishingTermsCurr = Series[VanishingTerms, {DiffExp`Symbols`x, 0, xAdd}];

        DiffExp`Utilities`PrintDebug["Reducing number of indeterminate constants, considering terms up to O[x]^", xAdd][3];
        OverdeterminedEqns = Flatten[Flatten[DiffExp`SeriesOps`LogxCoeffList /@ VanishingTermsCurr][[All, 3]]];

        If[Length[OverdeterminedEqns] > 0,

          {Cmat, Cb} = {#[[2]], -#[[1]]} &@CoefficientArrays[DeleteCases[OverdeterminedEqns, True], csCurr];

          CsPartSol = LinearSolve[Cmat, Cb, ZeroTest -> (N[DiffExp`Utilities`LSPChop@Expand@Normal[#1], DiffExp`State`LinearSolveChopPrecisionVal] == 0 &)];

          CsNullVectors = Cmat // NullSpace[#, Tolerance -> 10^-DiffExp`State`LinearSolveChopPrecisionVal] &;
          CsGeneralSol = CsPartSol + Sum[CsNullVectors[[ind]] Subscript[c, ind], {ind, Length[CsNullVectors]}] // DiffExp`SeriesOps`SExpand;

          CsReps = Thread[csCurr -> CsGeneralSol] // DiffExp`Utilities`PChop;

          IndeterminatesRemaining = DiffExp`Utilities`GetCases[CsGeneralSol, Subscript[c, __]];

        ];

        xAdd += 1;
      ];

      If[Length[IndeterminatesRemaining] < Length[fGeneral],
        DiffExp`State`LastErrorContext = {Cmat, Cb, CsGeneralSol, fGeneral};
        DiffExp`Utilities`ReportError["There was an error in obtaining the solutions for integrals ", intind, " at order ", epsord, " using the \"VOP\" strategy. Try decreasing the value of the option ChopPrecision or increasing the WorkingPrecision."];
      ];

      fGeneral = fGeneral /. CsReps;

      ,

      fGeneral = fGeneral /. Table[csCurr[[ind]] -> Subscript[c, ind], {ind, Length@csCurr}];
    ];


    cIndices = DiffExp`Utilities`GetCases[fGeneral, Subscript[c, i_]];

    FMat = Table[fGeneral /. myc -> 1 /. Subscript[c, _] :> 0, {myc, cIndices}] // Transpose;
    FMatInv = DiffExp`Wronskian`MatrixLogxInverse[FMat];
    BufferedData[intind] = {FMat, FMatInv};


  ];


  {FMat, FMatInv} = BufferedData[intind];

  (* Use shared helper for GMat computation *)
  DiffExp`Utilities`PrintDebug["Setting up general solution."][3];
  {cIndices, GMat} = ComputeGMat[FMat, FMatInv, bVec, intind];
  BMat = 1/Length@intind Table[bVec, {iind, intind // Length}] // Transpose;

  If[MemberQ[DiffExp`State`CurrCrosscheckFlags, "GeneralSolutionMatrix"] === True,
    DiffExp`Utilities`PrintDebug["Cross-checking GMat with differential equations."][1];
    CrossC = DiffExp`SeriesOps`SD[GMat, DiffExp`Symbols`x] - DiffExp`SeriesOps`MatrixMultiplySExpand[DiffExp`State`DEqnMatricesExpanded[line][0][[intind, intind]], GMat + O[DiffExp`Symbols`x]^DiffExp`State`ICrossCheckVerifyResultOrder] - BMat // (DiffExp`Utilities`CPChop@*DiffExp`SeriesOps`SExpand);
    DiffExp`State`LastErrorContext = {FMat, GMat, FMatInv, BMat, DiffExp`State`DEqnMatricesExpanded[line][0][[intind, intind]], CrossC, MtildeMatrix, MyWronsk};
    DiffExp`Utilities`PrintDebug["Found: ", CrossC + O[DiffExp`Symbols`x]^DiffExp`State`ICrossCheckPrintResultOrder // DiffExp`SeriesOps`SN][3];
    If[
      !(SameQ @@ Append[CrossC + O[DiffExp`Symbols`x]^DiffExp`State`ICrossCheckVerifyResultOrder // Normal // Flatten, 0])
      ,
      DiffExp`Utilities`ReportError["Cross-check of solution matrix failed"];
    ];
  ];

  fGeneral = Total[GMat // Transpose];

  {cIndices, fGeneral, BufferedData}
];

(* Dispatch to appropriate strategy based on configuration and problem type *)
DispatchStrategy[intind_, bVec_, line_, epsord_, BufferedData_] := Module[
  {cIndices, fGeneral, result},

  Which[
    (* Simple case: single integral without homogeneous components *)
    Length[intind] === 1 && DiffExp`Utilities`PChop[DiffExp`State`DEqnMatricesExpanded[line][0][[intind, intind]]] == {{0}},
    {cIndices, fGeneral} = SolveSimple[intind, bVec, line, epsord];
    {cIndices, fGeneral, BufferedData}

    (* Default strategy *)
    , DiffExp`State`FEC[IntegrationStrategy] === "Default",
    SolveDefault[intind, bVec, line, epsord, BufferedData]

    (* VOP strategy *)
    , DiffExp`State`FEC[IntegrationStrategy] === "VOP" || DiffExp`State`FEC[IntegrationStrategy] === "VariationOfParameters",
    SolveVOP[intind, bVec, line, epsord, BufferedData]

    (* VOPAlt strategy *)
    , DiffExp`State`FEC[IntegrationStrategy] === "VOPAlt",
    SolveVOPAlt[intind, bVec, line, epsord, BufferedData]
  ]
];

End[];

EndPackage[];
