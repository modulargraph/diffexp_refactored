(* IntegrationStrategies/VOP.m *)
(* SolveVOP and SolveVOPAlt integration strategies *)

(* VOP integration strategy *)
(* Uses variation of parameters method *)
SolveVOP[ctx_Association, bVec_, epsord_, cacheIn_Association] := Module[
  {systemSize, MatricesMSupj, nthOrderDifferentialEquations, nthOrderDifferentialEquationsSolutions,
   SolveFrom, MtildeMatrix, HighestOrderDifferentialEquationPosition,
   homogeneousSolutions, numSolutions, inhomogeneousTerm,
   bSupjVecs, fGeneral, cIndices, VanishingTerms,
   DerivativeVec, boundaryDerivatives, MtildeInv, wronskianColumn, wronskianMatrices, wronskianDetInverse, c, PlaceHolder,
   cache = cacheIn},

  DiffExp`Utilities`PrintDebug["Solving by variation of parameters. ", ctx["Label"]][3];

  systemSize = ctx["SystemSize"];

  (* Initialize wronskianMatrices and wronskianDetInverse as local associations *)
  wronskianMatrices = Association[{}];
  wronskianDetInverse = Association[{}];

  If[epsord === 0 && !KeyExistsQ[cache, "VOP"],
    (* Use shared helper functions for common computations *)
    MatricesMSupj = ComputeMatricesMSupj[ctx];
    nthOrderDifferentialEquations = ComputeNthOrderDiffEqns[MatricesMSupj, systemSize];

    If[ctx["HomogeneousSolve"] === "Expand",
      MatricesMSupj = DiffExp`SeriesOps`DiffExpSeries[MatricesMSupj];
    ];

    HighestOrderDifferentialEquationPosition = FirstPosition[Length /@ nthOrderDifferentialEquations, Max[Length /@ nthOrderDifferentialEquations]][[1]];
    If[Length[nthOrderDifferentialEquations[[HighestOrderDifferentialEquationPosition]]] < systemSize + 1,
      SolveFrom = Range@Length@nthOrderDifferentialEquations;
      nthOrderDifferentialEquationsSolutions = DiffExp`Frobenius`FrobeniusSolutions /@ nthOrderDifferentialEquations;
      If[systemSize == 1,
        If[ctx["HomogeneousSolve"] === "Expand",
          MtildeMatrix = DiffExp`Utilities`PChop@DiffExp`SeriesOps`DiffExpSeries[MatricesMSupj[[All, SolveFrom[[1]]]][[Range[systemSize], Range[systemSize]]]];
          ,
          MtildeMatrix = MatricesMSupj[[All, SolveFrom[[1]]]][[Range[systemSize], Range[systemSize]]];
        ];
        ,
        MtildeMatrix = {};
      ];
      ,
      SolveFrom = {HighestOrderDifferentialEquationPosition};

      If[ctx["HomogeneousSolve"] === "Expand",
        MtildeMatrix = DiffExp`Utilities`PChop@DiffExp`SeriesOps`DiffExpSeries[MatricesMSupj[[All, SolveFrom[[1]]]][[Range[systemSize], Range[systemSize]]]];
        ,
        MtildeMatrix = MatricesMSupj[[All, SolveFrom[[1]]]][[Range[systemSize], Range[systemSize]]];
      ];
      nthOrderDifferentialEquationsSolutions = ConstantArray[Null, Length@nthOrderDifferentialEquations];
      nthOrderDifferentialEquationsSolutions[[HighestOrderDifferentialEquationPosition]] = DiffExp`Frobenius`FrobeniusSolutions@nthOrderDifferentialEquations[[HighestOrderDifferentialEquationPosition]];
    ];

    DiffExp`Utilities`PrintDebug["Determining Wronskian matrices."][3];

    Table[
      homogeneousSolutions = nthOrderDifferentialEquationsSolutions[[myintind]];

      (* Use shared helper for Wronskian matrix computation *)
      wronskianMatrices[myintind] = ComputeWronskianMatrix[homogeneousSolutions];

      (* VOP-specific: compute inverse Wronskian determinant *)
      DiffExp`Utilities`PrintDebug["Determining 1/WronskianDet: ", myintind][3];
      wronskianDetInverse[myintind] = DiffExp`SeriesOps`DiffExpSeries[1/(wronskianMatrices[myintind] // Det // DiffExp`SeriesOps`SExpand), ctx["ExpansionOrder"]];
      DiffExp`Utilities`PrintDebug["Done determining 1/WronskianDet: ", myintind][3];
      , {myintind, SolveFrom}
    ];

    cache["VOP"] = {MatricesMSupj, nthOrderDifferentialEquations, nthOrderDifferentialEquationsSolutions, SolveFrom, MtildeMatrix, wronskianMatrices, wronskianDetInverse};
  ];

  {MatricesMSupj, nthOrderDifferentialEquations, nthOrderDifferentialEquationsSolutions, SolveFrom, MtildeMatrix, wronskianMatrices, wronskianDetInverse} = cache["VOP"];

  DiffExp`Utilities`PrintDebug["Determining b's.."][3];
  bSupjVecs = FoldList[DiffExp`SeriesOps`SExpand[#2 . bVec + DiffExp`SeriesOps`SD[#1, DiffExp`Symbols`x]] &, ConstantArray[0, systemSize], Delete[MatricesMSupj, -1]];

  DiffExp`Utilities`PrintDebug["Setting up general solutions. ", ctx["Label"]][3];
  fGeneral = Table[
    homogeneousSolutions = nthOrderDifferentialEquationsSolutions[[myintind]];
    numSolutions = homogeneousSolutions // Length;
    inhomogeneousTerm = nthOrderDifferentialEquations[[myintind]] . bSupjVecs[[All, myintind]][[Range@(numSolutions + 1)]];

    DiffExp`SeriesOps`SExpand@Sum[
      wronskianColumn = wronskianMatrices[myintind];
      wronskianColumn[[All, ind]] = Append[ConstantArray[0, numSolutions - 1], PlaceHolder];

      (homogeneousSolutions[[ind]] (Subscript[c, myintind, ind] + DiffExp`Integration`DiffExpIntegrate[
        wronskianDetInverse[myintind] ((wronskianColumn // Det // DiffExp`SeriesOps`SExpand) /. PlaceHolder -> inhomogeneousTerm // DiffExp`SeriesOps`SExpand) // DiffExp`SeriesOps`SExpand
      ]) // DiffExp`SeriesOps`SExpand) // DiffExp`SeriesOps`SExpand

      , {ind, homogeneousSolutions // Length}]

    , {myintind, SolveFrom}
  ];

  DiffExp`Utilities`PrintDebug["Done. ", ctx["Label"]][3];

  If[Length[SolveFrom] > 1,
    VanishingTerms = DiffExp`SeriesOps`SExpand@(DiffExp`SeriesOps`SD[fGeneral, DiffExp`Symbols`x] - ctx["AMatExpanded"] . fGeneral - bVec);

    IndeterminatesRemaining = currentConstants = DiffExp`Utilities`GetCases[fGeneral, Subscript[c, i_, j_]];

    If[!DeleteDuplicates[Normal[VanishingTerms]] === {0},
      seriesOrderOffset = -1;

      While[Length[IndeterminatesRemaining] > Length[fGeneral],
        VanishingTermsCurr = Series[VanishingTerms, {DiffExp`Symbols`x, 0, seriesOrderOffset}];

        DiffExp`Utilities`PrintDebug["Reducing number of indeterminate constants, considering terms up to O[x]^", seriesOrderOffset][3];
        OverdeterminedEqns = Flatten[Flatten[DiffExp`SeriesOps`LogxCoeffList /@ VanishingTermsCurr][[All, 3]]];

        If[Length[OverdeterminedEqns] > 0,

          {Cmat, Cb} = {#[[2]], -#[[1]]} &@CoefficientArrays[DeleteCases[OverdeterminedEqns, True], currentConstants];

          CsPartSol = LinearSolve[Cmat, Cb, ZeroTest -> (N[DiffExp`Utilities`LSPChop@Expand@Normal[#1], ctx["LinearSolveChopPrecision"]] == 0 &)];

          CsNullVectors = Cmat // NullSpace[#, Tolerance -> 10^-ctx["LinearSolveChopPrecision"]] &;
          CsGeneralSol = CsPartSol + Sum[CsNullVectors[[ind]] Subscript[c, ind], {ind, Length[CsNullVectors]}] // DiffExp`SeriesOps`SExpand;

          constantReplacements = Thread[currentConstants -> CsGeneralSol] // DiffExp`Utilities`PChop;

          IndeterminatesRemaining = DiffExp`Utilities`GetCases[CsGeneralSol, Subscript[c, __]];

        ];

        seriesOrderOffset += 1;
      ];

      If[Length[IndeterminatesRemaining] < Length[fGeneral],
        DiffExp`State`LastErrorContext = {Cmat, Cb, CsGeneralSol};
        DiffExp`Utilities`ReportError["There was an error in obtaining the solutions for integrals ", ctx["Label"], " at order ", epsord, " using the \"VOP\" strategy. Try decreasing the value of the option ChopPrecision or increasing the WorkingPrecision."];
      ];

      fGeneral = fGeneral /. constantReplacements;

      ,

      fGeneral = fGeneral /. Table[currentConstants[[ind]] -> Subscript[c, ind], {ind, Length@currentConstants}];
    ];

    ,

    (* We can currently put fGeneral[[1]], as SolveFrom is either one element, or all elements. *)
    DerivativeVec = Table[DiffExp`SeriesOps`SD[fGeneral[[1]], Sequence @@ DiffExp`Utilities`CA[DiffExp`Symbols`x, jind - 1]], {jind, systemSize}] // DiffExp`SeriesOps`SExpand;
    (* But for boundaryDerivatives we should choose SolveFrom[[1]] *)
    boundaryDerivatives = bSupjVecs[[All, SolveFrom[[1]]]][[Range@systemSize]];

    Check[
      If[ctx["HomogeneousSolve"] === "Expand",
        MtildeInv = Inverse[MtildeMatrix, Method -> "DivisionFreeRowReduction", ZeroTest -> (Normal[#] == 0 &)];,
        MtildeInv = DiffExp`SeriesOps`DiffExpSeries[Inverse[MtildeMatrix] // Together];
      ];
      fGeneral = DiffExp`Utilities`PChop@(MtildeInv . (DerivativeVec - boundaryDerivatives) // DiffExp`SeriesOps`SExpand);
      ,
      DiffExp`Utilities`PrintWarning["Encountered numerical instability while inverting Mtilde. Turning off DivisionFreeRowReduction and trying again.."];
      fGeneral = DiffExp`Utilities`PChop@(Inverse[MtildeMatrix] . (DerivativeVec - boundaryDerivatives) // DiffExp`SeriesOps`SExpand);
    ];

    fGeneral = fGeneral /. Thread[DiffExp`Utilities`GetCases[fGeneral, Subscript[c, i_, j_]] -> Table[Subscript[c, ind], {ind, systemSize}]];

  ];

  cIndices = DiffExp`Utilities`GetCases[fGeneral, Subscript[c, i_]];

  {cIndices, fGeneral, cache}
];

(* VOPAlt integration strategy *)
(* Alternative variation of parameters method *)
SolveVOPAlt[ctx_Association, bVec_, epsord_, cacheIn_Association] := Module[
  {systemSize, MatricesMSupj, nthOrderDifferentialEquations, nthOrderDifferentialEquationsSolutions,
   SolveFrom, MtildeMatrix, HighestOrderDifferentialEquationPosition,
   homogeneousSolutions, fGeneral, cIndices, VanishingTerms, IndeterminatesRemaining, currentConstants,
   seriesOrderOffset, VanishingTermsCurr, OverdeterminedEqns, Cmat, Cb, CsPartSol, CsNullVectors,
   CsGeneralSol, constantReplacements, FMat, FMatInv, GMat, BMat, crossCheck, wronskianMatrices, c,
   cache = cacheIn},

  DiffExp`Utilities`PrintDebug["Solving by alternative strategy. ", ctx["Label"]][3];

  systemSize = ctx["SystemSize"];

  If[epsord === 0 && !KeyExistsQ[cache, "FMat"],
    (* Use shared helper functions for common computations *)
    MatricesMSupj = ComputeMatricesMSupj[ctx];
    nthOrderDifferentialEquations = ComputeNthOrderDiffEqns[MatricesMSupj, systemSize];

    If[ctx["HomogeneousSolve"] === "Expand",
      MatricesMSupj = DiffExp`SeriesOps`DiffExpSeries[MatricesMSupj];
    ];

    HighestOrderDifferentialEquationPosition = FirstPosition[Length /@ nthOrderDifferentialEquations, Max[Length /@ nthOrderDifferentialEquations]][[1]];
    SolveFrom = Range@Length@nthOrderDifferentialEquations;
    nthOrderDifferentialEquationsSolutions = DiffExp`Frobenius`FrobeniusSolutions /@ nthOrderDifferentialEquations;
    If[systemSize == 1,
      If[ctx["HomogeneousSolve"] === "Expand",
        MtildeMatrix = DiffExp`Utilities`PChop@DiffExp`SeriesOps`DiffExpSeries[MatricesMSupj[[All, SolveFrom[[1]]]][[Range[systemSize], Range[systemSize]]]];
        ,
        MtildeMatrix = MatricesMSupj[[All, SolveFrom[[1]]]][[Range[systemSize], Range[systemSize]]];
      ];
      ,
      MtildeMatrix = {};
    ];
    DiffExp`Utilities`PrintDebug["Determining Wronskian matrices."][3];

    (* Initialize wronskianMatrices as local association *)
    wronskianMatrices = Association[{}];

    (* Use shared helper for Wronskian matrix computation *)
    Table[
      homogeneousSolutions = nthOrderDifferentialEquationsSolutions[[myintind]];
      wronskianMatrices[myintind] = ComputeWronskianMatrix[homogeneousSolutions];
      , {myintind, SolveFrom}
    ];

    DiffExp`Utilities`PrintDebug["Combining homogeneous solutions. ", ctx["Label"]][3];
    fGeneral = Table[
      homogeneousSolutions = nthOrderDifferentialEquationsSolutions[[myintind]];

      DiffExp`SeriesOps`SExpand@Sum[homogeneousSolutions[[ind]] (Subscript[c, myintind, ind]) // DiffExp`SeriesOps`SExpand, {ind, homogeneousSolutions // Length}]

      , {myintind, SolveFrom}
    ];

    DiffExp`Utilities`PrintDebug["Done. ", ctx["Label"]][3];

    VanishingTerms = DiffExp`SeriesOps`SExpand@(DiffExp`SeriesOps`SD[fGeneral, DiffExp`Symbols`x] - ctx["AMatExpanded"] . fGeneral);

    IndeterminatesRemaining = currentConstants = DiffExp`Utilities`GetCases[fGeneral, Subscript[c, i_, j_]];

    If[!DeleteDuplicates[Normal[VanishingTerms]] === {0},
      seriesOrderOffset = -1;

      While[Length[IndeterminatesRemaining] > Length[fGeneral],
        VanishingTermsCurr = Series[VanishingTerms, {DiffExp`Symbols`x, 0, seriesOrderOffset}];

        DiffExp`Utilities`PrintDebug["Reducing number of indeterminate constants, considering terms up to O[x]^", seriesOrderOffset][3];
        OverdeterminedEqns = Flatten[Flatten[DiffExp`SeriesOps`LogxCoeffList /@ VanishingTermsCurr][[All, 3]]];

        If[Length[OverdeterminedEqns] > 0,

          {Cmat, Cb} = {#[[2]], -#[[1]]} &@CoefficientArrays[DeleteCases[OverdeterminedEqns, True], currentConstants];

          CsPartSol = LinearSolve[Cmat, Cb, ZeroTest -> (N[DiffExp`Utilities`LSPChop@Expand@Normal[#1], ctx["LinearSolveChopPrecision"]] == 0 &)];

          CsNullVectors = Cmat // NullSpace[#, Tolerance -> 10^-ctx["LinearSolveChopPrecision"]] &;
          CsGeneralSol = CsPartSol + Sum[CsNullVectors[[ind]] Subscript[c, ind], {ind, Length[CsNullVectors]}] // DiffExp`SeriesOps`SExpand;

          constantReplacements = Thread[currentConstants -> CsGeneralSol] // DiffExp`Utilities`PChop;

          IndeterminatesRemaining = DiffExp`Utilities`GetCases[CsGeneralSol, Subscript[c, __]];

        ];

        seriesOrderOffset += 1;
      ];

      If[Length[IndeterminatesRemaining] < Length[fGeneral],
        DiffExp`State`LastErrorContext = {Cmat, Cb, CsGeneralSol, fGeneral};
        DiffExp`Utilities`ReportError["There was an error in obtaining the solutions for integrals ", ctx["Label"], " at order ", epsord, " using the \"VOP\" strategy. Try decreasing the value of the option ChopPrecision or increasing the WorkingPrecision."];
      ];

      fGeneral = fGeneral /. constantReplacements;

      ,

      fGeneral = fGeneral /. Table[currentConstants[[ind]] -> Subscript[c, ind], {ind, Length@currentConstants}];
    ];


    cIndices = DiffExp`Utilities`GetCases[fGeneral, Subscript[c, i_]];

    FMat = Table[fGeneral /. myc -> 1 /. Subscript[c, _] :> 0, {myc, cIndices}] // Transpose;
    FMatInv = DiffExp`Wronskian`MatrixLogxInverse[FMat];
    cache["FMat"] = FMat;
    cache["FMatInv"] = FMatInv;
  ];


  {FMat, FMatInv} = {cache["FMat"], cache["FMatInv"]};

  (* Use shared helper for GMat computation *)
  DiffExp`Utilities`PrintDebug["Setting up general solution."][3];
  {cIndices, GMat} = ComputeGMat[FMat, FMatInv, bVec];
  BMat = 1/systemSize Table[bVec, {iind, systemSize}] // Transpose;

  If[MemberQ[ctx["CrosscheckFlags"], "GeneralSolutionMatrix"] === True,
    DiffExp`Utilities`PrintDebug["Cross-checking GMat with differential equations."][1];
    crossCheck = DiffExp`SeriesOps`SD[GMat, DiffExp`Symbols`x] - DiffExp`SeriesOps`MatrixMultiplySExpand[ctx["AMatExpanded"], GMat + O[DiffExp`Symbols`x]^ctx["CrossCheckVerifyOrder"]] - BMat // (DiffExp`Utilities`CPChop@*DiffExp`SeriesOps`SExpand);
    DiffExp`State`LastErrorContext = {FMat, GMat, FMatInv, BMat, ctx["AMatExpanded"], crossCheck};
    DiffExp`Utilities`PrintDebug["Found: ", crossCheck + O[DiffExp`Symbols`x]^ctx["CrossCheckPrintOrder"] // DiffExp`SeriesOps`SN][3];
    If[
      !(SameQ @@ Append[crossCheck + O[DiffExp`Symbols`x]^ctx["CrossCheckVerifyOrder"] // Normal // Flatten, 0])
      ,
      DiffExp`Utilities`ReportError["Cross-check of solution matrix failed"];
    ];
  ];

  fGeneral = Total[GMat // Transpose];

  {cIndices, fGeneral, cache}
];
