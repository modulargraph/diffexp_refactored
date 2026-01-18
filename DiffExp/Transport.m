(* ::Package:: *)

(* DiffExp Transport Subpackage *)
(* This package provides the main transport and integration functions *)

BeginPackage["DiffExp`Transport`", {
  "DiffExp`Symbols`",
  "DiffExp`State`",
  "DiffExp`Utilities`",
  "DiffExp`SeriesOps`",
  "DiffExp`Integration`",
  "DiffExp`Pade`",
  "DiffExp`Mobius`",
  "DiffExp`AnalyticContinuation`",
  "DiffExp`LineSegmentation`",
  "DiffExp`Frobenius`",
  "DiffExp`Wronskian`",
  "DiffExp`MatrixLoading`"
}];

(* Transport functions *)
PrepareBoundaryConditions::usage = "PrepareBoundaryConditions[bcs_List, line_List] prepares boundary conditions for use with IntegrateSystem or TransportTo.";
IntegrateSystem::usage = "IntegrateSystem[line_List] obtains general series solutions along a line. IntegrateSystem[bcs_List, line_List] uses boundary conditions.";
TransportTo::usage = "TransportTo[bcs_List, line_List, to_:1, save_:False] transports boundary conditions to arbitrary points.";

Begin["`Private`"];

(* Prepare boundary conditions *)
PrepareBoundaryConditions[bcs_List, line2_Association | line2_List] := Module[
  {line, CoeffList = {}, Coeffs, CoeffSer, bcs1, OneSer = (1 // N[#, DiffExp`State`FEWorkingPrecision] &), Mask, ispoint, LineRat, tmp},

  DiffExp`Utilities`PrintDebug["Preparing boundary conditions."][1];

  If[line2[[0]] === List, line = line2 // Association // KeySort, line = line2 // KeySort];

  If[Accuracy[line] < DiffExp`State`FEWorkingPrecision,
    line = line // Normal // SetPrecision[#, 2 DiffExp`State`FEWorkingPrecision] & // Association;
    DiffExp`Utilities`PrintWarning["Accuracy of the line/point is lower than the working precision. The precision has been artificially increased."];
  ];

  If[Length[bcs] != DiffExp`State`NumIntegrals,
    DiffExp`Utilities`PrintWarning["The number of integrals does not match the size of the system of differential equations."];
  ];

  DiffExp`AnalyticContinuation`PrepareAnalyticContinuation[line];

  ispoint = DiffExp`Utilities`IsPoint[line];

  If[ispoint && Count[DiffExp`Utilities`PChop[(DiffExp`State`MatricesIrreducibleFactors /. line // Together)], 0] > 0,
    DiffExp`Utilities`PrintWarning["The boundary conditions are given as numerical values at a point, but the point lies on a singularity of the differential equations. This only works if the asymptotic limit is finite."];
  ];

  Table[
    Switch[{bcs[[bcind, 0]], bcs[[bcind]]}
      , {List, _},
      DiffExp`Utilities`PrintDebug["Integral ", bcind, ": List. Assuming first entry ~ \!\(\*SuperscriptBox[\(\[Epsilon]\), \(0\)]\)."][1];
      If[Length[bcs[[bcind]]] < DiffExp`State`EpsilonOrderVal + 1,
        DiffExp`Utilities`ReportError["Too few coefficients given for boundary conditions up to order \[Epsilon]^", DiffExp`State`EpsilonOrderVal];
      ];
      Coeffs = bcs[[bcind]][[Range[DiffExp`State`EpsilonOrderVal + 1]]];
      , {String, "?"},
      DiffExp`Utilities`PrintInfo["Integral ", bcind, ": Ignoring boundary conditions."][1];
      Coeffs = DiffExp`Utilities`CA["?", DiffExp`State`EpsilonOrderVal + 1];

      , _,
      DiffExp`Utilities`PrintDebug["Integral ", bcind, ": Assuming closed form expression."][1];
      CoeffSer = DiffExp`SeriesOps`SeriesAlways[bcs[[bcind]], {DiffExp`Symbols`\[Epsilon], 0, DiffExp`State`EpsilonOrderVal}] // DiffExp`SeriesOps`SExpand;

      (* Sanity check *)
      If[(CoeffSer // DiffExp`SeriesOps`SeriesMinPower) < 0, DiffExp`Utilities`ReportError["The boundary conditions should start at finite order."]];

      Coeffs = Table[SeriesCoefficient[CoeffSer, {DiffExp`Symbols`\[Epsilon], 0, ord}], {ord, 0, DiffExp`State`EpsilonOrderVal}];

      (* We put everything below the leading eps coefficient to zero for all orders of x. *)
      If[!ispoint,
        If[!bcs[[bcind]] === 0,
          Do[
            If[!ispoint, DiffExp`Utilities`PrintInfo["Assuming that integral ", bcind, " is exactly zero at epsilon order ", iind - 1, "."][1];];
            Coeffs[[iind]] = SeriesData[DiffExp`Symbols`x, 0, List[], 0, DiffExp`State`ExpansionOrderVal + 1, 1]
            , {iind, CoeffSer[[4]]}];
        ];
      ];
    ];

    AppendTo[CoeffList, Coeffs];
    , {bcind, bcs // Length}
  ];

  bcs1 = If[!ispoint,
    Table[
      If[!(#[[ind]][[0]] === SeriesData),
        (
          DiffExp`SeriesOps`LeadingCoefficientSeries[(#[[ind]] /. line) * OneSer, 2]
        ) /. OneSer "?" + O[DiffExp`Symbols`x]^(1/2) -> "?",
        #[[ind]]
      ]
      , {ind, Length[#]}
    ]
    ,
    If[DiffExp`Utilities`DependsQ[Normal[#], DiffExp`Symbols`x],
      DiffExp`Utilities`ReportError["The boundary terms that are provided depend on the line parameter ", DiffExp`Symbols`x, ", but the line itself does not."];
    ];
    OneSer * # /. line
  ] & /@ CoeffList;

  If[!ispoint,
    Mask = MapAt[Switch[#[[0]],
      SeriesData,
      tmp = #;
      tmp[[3]] = Table["(...)", {ind, #[[3]] // Length}];
      tmp,
      _, #] &, bcs1, {All, All}] // TableForm;

    DiffExp`Utilities`PrintInfo["Prepared boundary conditions in asymptotic limit, of the form:"][1];
    DiffExp`Utilities`PrintInfo[Mask][1];
  ];

  {line, bcs1} /. Log[a_ DiffExp`Symbols`x] /; NumericQ[a] :> Log[a] + DiffExp`Symbols`Logx /. Log[DiffExp`Symbols`x] -> DiffExp`Symbols`Logx
];

(* Main integration function *)
IntegrateSystem[bcs2 : _List : "?", line2_Association | line2_List, opts2_ : {}] := Module[
  {bcs, line, lrln, lrln2, BCSRelevant, relevantinds, IgnorePositions, CurrBlock, HomogeneousEquation, Solns, MtildeMat, GMat, FMat, FMatInv, FMatInvBMat, BMat, CrossC, bVec, IntegrationData, fGeneral, FixAt, BoundaryEqns1, BoundaryEqns2, boundarysols, cIndices, Cmat, Cb, cpivs, csol, NewResults, Wronsk, ll, FMat2, NMat, HomogeneousEquation2, MtildeMat2, NMat2, Solns2, Wronsk2, WronskInvPrime, WWinvprimeprod, WronskInv, opts = opts2, DEqnMatricesExpandedCopy, TurnOffPade, MyN, MatricesMSupj, nthOrderDifferentialEquations, TmpSols, HighestOrderDifferentialEquationPosition, SolveFrom, nthOrderDifferentialEquationsSolutions, MtildeMatrix, bSupjVecs, MyHomogeneousSolutions, MyNumberOfSolutions, MyInhomogeneousTerm, VanishingTerms, csCurr, OverdeterminedEqns, CsPartSol, CsNullVectors, CsGeneralSol, CsReps, DerivativeVec, BDerVec, LineRat, csFreedom, CouldntSolve, xAdd, IndeterminatesRemaining, VanishingTermsCurr, LogsPresent, AlgebraicRootsPresent, TmpSolutionsNormal, MatricesMSupjExp, InhomPlaceHolder, MtildeInv, BufferedData, MyWronsk, MyWronskDetInv, IntegrationDataTab, bVec0, bVec1, bVecRest, jinds, IntegrationDataTab0jind, CurrInvWronskSolver, MyWi, c, PlaceHolder},

  If[line2[[0]] === List, line = line2 // Association // KeySort, line = line2 // KeySort];

  If[!MemberQ[opts, "TransportToCall"],
    DiffExp`State`BenchmarkData = Association[];
    DiffExp`State`BenchmarkData["Segments"] = Association[];
  ];

  If[!KeyExistsQ[DiffExp`State`BenchmarkData["Segments"], line // N],
    DiffExp`State`BenchmarkData["Segments"][line // N] = Association[];
    DiffExp`State`BenchmarkData["Segments"]["ComputationTime"] = AbsoluteTime[];
  ];

  If[Or @@ (!DiffExp`Utilities`DependsQ[Keys@line, #] & /@ DiffExp`State`ExternalScalesVal),
    DiffExp`Utilities`ReportError["The line does not fix all kinematic invariants and masses!"];
  ];

  If[!MemberQ[opts, "TransportToCall"],
    DiffExp`Utilities`PrintInfo["Obtaining series solutions along provided line.."][1];
    If[Accuracy[line] < DiffExp`State`FEWorkingPrecision,
      line = line // Normal // SetPrecision[#, 2 DiffExp`State`FEWorkingPrecision] & // Association;
      DiffExp`Utilities`PrintWarning["Accuracy of the line/point is lower than the working precision. The precision has been artificially increased."];
    ];
  ];

  TurnOffPade[] := If[DiffExp`State`FEC[UsePade] === True && MemberQ[opts, "TransportToCall"],
    DiffExp`Utilities`PrintWarning["Due to the presence of free parameters Pade approximants will be disabled. You will have to manually enable them again in the configuration by setting UsePade -> True!"];
    DiffExp`State`DiffExpConfiguration[UsePade] = False;
  ];

  BufferedData = Association[{}];
  MyWronsk = Association[{}];
  MyWronskDetInv = Association[{}];

  (* To deal with the output of SaveExpansions = True *)
  If[MatchQ[bcs2, {{a_Association, __}, _}],
    bcs = bcs2[[1]];,
    bcs = bcs2;
  ];

  If[!(bcs2 === "?"),
    {bcs, FixAt} = DiffExp`LineSegmentation`CheckBoundaryConditionsAndReparametrize[bcs, line];

    DiffExp`Utilities`PrintDebug["Boundary conditions are given by:"][2];
    DiffExp`Utilities`PrintDebug["Line: ", bcs[[1]] // Normal // N // Association][2];
    DiffExp`Utilities`PrintDebug["Conditions: ", bcs[[2]] // Normal // N][2];
    DiffExp`Utilities`PrintDebug["Fixing boundary conditions at x = ", FixAt][2];

    ,

    TurnOffPade[];
  ];

  If[!Complement[DiffExp`State`ExternalScalesVal, Keys[line]] === {},
    DiffExp`Utilities`ReportError["Line does not fix all kinematic variables!"];
  ];

  DiffExp`AnalyticContinuation`PrepareAnalyticContinuation[line];

  If[!MemberQ[opts, "TransportToCall"],
    (* We do this because the user might have aborted TransportTo while using the option AccuracyGoal. *)
    DiffExp`State`ExpansionOrderVal = DiffExp`State`FEC[ExpansionOrder];
    DiffExp`MatrixLoading`ClearMatrices[];

    DiffExp`State`BenchmarkData["Segments"][line // N]["MatrixExpansion"] = AbsoluteTime[];
  ];

  DiffExp`MatrixLoading`PrepareMatrices[line];

  If[!MemberQ[opts, "TransportToCall"],
    DiffExp`State`BenchmarkData["Segments"][line // N]["MatrixExpansion"] = AbsoluteTime[] - DiffExp`State`BenchmarkData["Segments"][line // N]["MatrixExpansion"];
  ];

  DiffExp`MatrixLoading`InitializeIntegrationSequence[line];

  IntegrationData = Association[{}];

  DEqnMatricesExpandedCopy = DiffExp`State`DEqnMatricesExpanded[line];

  DiffExp`State`BenchmarkData["Segments"][line // N]["ComputationTimes"] = Association[];
  DiffExp`State`BenchmarkData["Segments"][line // N]["ComputationTimes"]["Integrals"] = Association[];

  DiffExp`State`BenchmarkData["Segments"][line // N]["HomogeneousSolveAllPreprocessing"] = Association[];
  DiffExp`State`BenchmarkData["Segments"][line // N]["HomogeneousSolveAllPreprocessing"]["Integrals"] = Association[];

  Do[
    DiffExp`Utilities`PrintInfo["Currently at order: \[Epsilon]^", epsord, "."][2];
    DiffExp`State`BenchmarkData["Segments"][line // N]["ComputationTimes"][epsord] = AbsoluteTime[];

    IntegrationDataTab = Table[IntegrationData[{ind, myeps}], {ind, DiffExp`State`NumIntegrals}, {myeps, 0, epsord - 1}];

    If[epsord > 0,
      bVec1 = DEqnMatricesExpandedCopy[1] . IntegrationDataTab[[All, epsord - 1 + 1]];
    ];
    If[epsord > 1,
      bVecRest = Sum[DEqnMatricesExpandedCopy[lind] . IntegrationDataTab[[All, epsord - lind + 1]] // DiffExp`SeriesOps`SExpand, {lind, 2, epsord}];
    ];

    Do[
      If[!KeyExistsQ[DiffExp`State`BenchmarkData["Segments"][line // N]["ComputationTimes"]["Integrals"], intind],
        DiffExp`State`BenchmarkData["Segments"][line // N]["ComputationTimes"]["Integrals"][intind] = Association[];
      ];
      DiffExp`State`BenchmarkData["Segments"][line // N]["ComputationTimes"]["Integrals"][intind][epsord] = AbsoluteTime[];

      DiffExp`Utilities`PrintInfo["Integrating integral(s) ", intind, " at order \[Epsilon]^", epsord, "."][3];

      DiffExp`Utilities`PrintDebug["Getting inhomogeneous terms."][3];

      jinds = Complement[Range[DiffExp`State`NumIntegrals], intind];
      IntegrationDataTab0jind = Table[IntegrationData[{ind, myeps}], {ind, jinds}, {myeps, 0, epsord}];
      bVec0 = DiffExp`SeriesOps`SExpand[DEqnMatricesExpandedCopy[0][[intind, jinds]] . IntegrationDataTab0jind[[All, epsord + 1]]];
      Which[epsord === 0,
        bVec = bVec0 // DiffExp`SeriesOps`SExpand;,
        epsord === 1,
        bVec = bVec0 + bVec1[[intind]] // DiffExp`SeriesOps`SExpand;,
        epsord > 1,
        bVec = bVec0 + bVec1[[intind]] + bVecRest[[intind]] // DiffExp`SeriesOps`SExpand;
      ];

      DiffExp`Utilities`PrintDebug["Done."][3];

      Which[
        (* Integration of integral without homogeneous components *)
        Length[intind] === 1 && DiffExp`Utilities`PChop[DiffExp`State`DEqnMatricesExpanded[line][0][[intind, intind]]] == {{0}},

        If[epsord === 0,
          DiffExp`State`BenchmarkData["Segments"][line // N]["HomogeneousSolveAllPreprocessing"]["Integrals"][intind] = 0;
        ];

        cIndices = {Subscript[c, 1]};

        If[DiffExp`Utilities`PChop[bVec] === {0},
          fGeneral = {Subscript[c, 1] + O[DiffExp`Symbols`x]^(DiffExp`State`ExpansionOrderVal + 1)};
          ,
          fGeneral = {DiffExp`Integration`DiffExpIntegrate[bVec[[1]], DiffExp`Symbols`x] + Subscript[c, 1] + O[DiffExp`Symbols`x]^(DiffExp`State`ExpansionOrderVal + 1)};
        ];

        , DiffExp`State`FEC[IntegrationStrategy] === "Default",
        (* Default integration strategy *)

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
                Global`DebugData = {WronskInvPrime, Wronsk, WWinvprimeprod};
                DiffExp`Utilities`ReportError["Warning, product of Wronskian inverse times Wronskian does not match identity. Try increasing \"ChopPrecision\" or \"WorkingPrecision\", or try decreasing \"ExpansionOrder\"."][1];
              ];
            ];
            Global`DebugData = {WWinvprimeprod};
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
            Global`DebugData = {FMat, DiffExp`State`DEqnMatricesExpanded[line][0][[intind, intind]], CrossC, MtildeMat, Wronsk};
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

        DiffExp`Utilities`PrintDebug["Setting up general solution."][3];
        BMat = 1/Length@intind Table[bVec, {iind, intind // Length}] // Transpose;
        cIndices = Table[Subscript[c, i], {i, intind // Length}];

        GMat = DiffExp`SeriesOps`MatrixMultiplySExpand[FMat, DiffExp`Integration`DiffExpIntegrate[DiffExp`SeriesOps`MatrixMultiplySExpand[FMatInv, BMat]] + DiagonalMatrix[cIndices]] // DiffExp`Utilities`PChop;

        If[MemberQ[DiffExp`State`CurrCrosscheckFlags, "GeneralSolutionMatrix"] === True,
          DiffExp`Utilities`PrintDebug["Cross-checking GMat with differential equations."][1];
          CrossC = DiffExp`SeriesOps`SD[GMat, DiffExp`Symbols`x] - DiffExp`SeriesOps`MatrixMultiplySExpand[DiffExp`State`DEqnMatricesExpanded[line][0][[intind, intind]], GMat + O[DiffExp`Symbols`x]^DiffExp`State`ICrossCheckVerifyResultOrder] - BMat // (DiffExp`Utilities`CPChop@*DiffExp`SeriesOps`SExpand);
          Global`DebugData = {FMat, GMat, FMatInv, BMat, DiffExp`State`DEqnMatricesExpanded[line][0][[intind, intind]], CrossC, MtildeMat, Wronsk};
          DiffExp`Utilities`PrintDebug["Found: ", CrossC + O[DiffExp`Symbols`x]^DiffExp`State`ICrossCheckPrintResultOrder // DiffExp`SeriesOps`SN][3];
          If[
            !(SameQ @@ Append[CrossC + O[DiffExp`Symbols`x]^DiffExp`State`ICrossCheckVerifyResultOrder // Normal // Flatten, 0])
            ,
            DiffExp`Utilities`ReportError["Cross-check of solution matrix failed"];
          ];
        ];

        fGeneral = Total[GMat // Transpose];

        , DiffExp`State`FEC[IntegrationStrategy] === "VOP" || DiffExp`State`FEC[IntegrationStrategy] === "VariationOfParameters",
        (* VOP integration strategy - variation of parameters *)

        DiffExp`Utilities`PrintDebug["Solving by variation of parameters. ", intind][3];

        MyN = intind // Length;

        If[epsord === 0 && !KeyExistsQ[BufferedData, intind],
          DiffExp`State`BenchmarkData["Segments"][line // N]["HomogeneousSolveAllPreprocessing"]["Integrals"][intind] = AbsoluteTime[];

          If[DiffExp`State`FEC["HomogeneousSolve"] === "Expand",
            MatricesMSupj = DiffExp`Utilities`PChop@NestList[DiffExp`SeriesOps`SExpand[(DiffExp`SeriesOps`SD[#, DiffExp`Symbols`x] + DiffExp`SeriesOps`MatrixMultiplySExpand[#, DiffExp`State`DEqnMatricesExpanded[line][0][[intind, intind]]])] &, IdentityMatrix[MyN], MyN];
            ,
            MatricesMSupj = NestList[Together[(D[#, DiffExp`Symbols`x] + # . DiffExp`State`DEqnMatricesFactored[line][0][[intind, intind]])] &, IdentityMatrix[MyN], MyN];
          ];

          nthOrderDifferentialEquations = Table[
            If[DiffExp`State`FEC["HomogeneousSolve"] === "Expand",
              TmpSols = MatricesMSupj[[All, ind]] // Transpose // NullSpace[#, Method -> "DivisionFreeRowReduction", Tolerance -> 10^-DiffExp`State`LinearSolveChopPrecisionVal] &;
              ,
              TmpSols = MatricesMSupj[[All, ind]] // Transpose // NullSpace // Together;
            ];
            TmpSols = Internal`DeleteTrailingZeros /@ TmpSols;
            TmpSols = DiffExp`SeriesOps`DiffExpSeries[MinimalBy[TmpSols, Length] // First];
            TmpSols/(TmpSols // Last)
            , {ind, MyN}
          ];

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

            MyWronsk[{intind, myintind}] = Table[
              DiffExp`SeriesOps`SD[MyHomogeneousSolutions[[mysolind]], Sequence @@ ConstantArray[DiffExp`Symbols`x, diffind]]
              , {diffind, 0, Length[MyHomogeneousSolutions] - 1}, {mysolind, Length[MyHomogeneousSolutions]}];

            DiffExp`Utilities`PrintDebug["Determining 1/WronskianDet: ", myintind][3];
            MyWronskDetInv[{intind, myintind}] = DiffExp`SeriesOps`DiffExpSeries[1/(MyWronsk[{intind, myintind}] // Det // DiffExp`SeriesOps`SExpand), DiffExp`State`ExpansionOrderVal];
            DiffExp`Utilities`PrintDebug["Done determining 1/WronskianDet: ", myintind][3];
            , {myintind, SolveFrom}
          ];

          BufferedData[intind] = {MatricesMSupj, nthOrderDifferentialEquations, nthOrderDifferentialEquationsSolutions, SolveFrom, MtildeMatrix};

          DiffExp`State`BenchmarkData["Segments"][line // N]["HomogeneousSolveAllPreprocessing"]["Integrals"][intind] = AbsoluteTime[] - DiffExp`State`BenchmarkData["Segments"][line // N]["HomogeneousSolveAllPreprocessing"]["Integrals"][intind];
        ];

        {MatricesMSupj, nthOrderDifferentialEquations, nthOrderDifferentialEquationsSolutions, SolveFrom, MtildeMatrix} = BufferedData[intind];

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
              Global`DebugData = {Cmat, Cb, CsGeneralSol};
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

      ];

      (* Code for fixing boundary conditions. *)
      If[!(bcs2 === "?"),
        DiffExp`Utilities`PrintDebug["General solution found. Fixing the indeterminates from boundary conditions."][2];

        BCSRelevant = bcs[[2]][[intind]][[All, epsord + 1]];
        IgnorePositions = Flatten[Position[BCSRelevant, "?"]];
        relevantinds = Complement[Range@Length[intind], IgnorePositions];
        BCSRelevant = BCSRelevant[[relevantinds]];

        If[FixAt === 0,
          BoundaryEqns1 = (DiffExp`Utilities`PChop@*DiffExp`SeriesOps`SExpand) /@ (fGeneral[[relevantinds]] - BCSRelevant);

          If[
            !SameQ[Append[BoundaryEqns1[[All, 0]], SeriesData]],
            DiffExp`Utilities`ReportError["Internal bug while fixing the solutions from the boundary data."];
          ];

          BoundaryEqns2 = # == 0 & /@ Flatten[{
            #[[3]] /. DiffExp`Symbols`Logx -> 0,
            Table[Coefficient[#[[3]], DiffExp`Symbols`Logx^kind], {kind, 1, DiffExp`State`IMaxLogOrder}]
          } & /@ BoundaryEqns1];

          ,

          BoundaryEqns1 = Table[Sum[
            cIndices[[ind]] DiffExp`Pade`SEval[DiffExp`SeriesOps`SApply[Coefficient[#, cIndices[[ind]]] &, term], FixAt]
            , {ind, cIndices // Length}
          ] + DiffExp`Pade`SEval[term /. (# -> 0 & /@ cIndices), FixAt], {term, fGeneral[[relevantinds]]}] - BCSRelevant;

          BoundaryEqns2 = # == 0 & /@ Flatten[BoundaryEqns1];
        ];

        BoundaryEqns2 = DeleteCases[BoundaryEqns2, True] // Expand;

        If[!BoundaryEqns2 === {},
          DiffExp`Utilities`PrintDebug["Boundary equations are: "][3];
          DiffExp`Utilities`PrintDebug[BoundaryEqns2 // N // ReplaceAll[c -> "c"]][3];

          {Cmat, Cb} = {#[[2]], -#[[1]]} &@CoefficientArrays[BoundaryEqns2, cIndices];

          CouldntSolve = False;
          If[!MemberQ[BoundaryEqns2, False],
            Check[csol = LinearSolve[Cmat, Cb, ZeroTest -> (N[DiffExp`Utilities`LSPChop@Expand@Normal[#1], DiffExp`State`LinearSolveChopPrecisionVal] == 0 &)];
              , CouldntSolve = True;];
            , CouldntSolve = True;];

          If[CouldntSolve,
            Global`DebugData = {BMat, fGeneral[[relevantinds]], BoundaryEqns1, BoundaryEqns2, BCSRelevant, Cmat, Cb, MtildeMat, Wronsk, WronskInv, WronskInvPrime, bVec, FMat, FMatInv, GMat, csol};
            DiffExp`Utilities`ReportError["Boundary conditions cannot be matched to general solution for integral(s): ", intind];
          ];

          DiffExp`Utilities`PrintDebug["Solutions: ", Thread[(ToString /@ cIndices) -> N[csol]]][2];

          csFreedom = Cmat // NullSpace[#, Method -> "DivisionFreeRowReduction", Tolerance -> 10^-DiffExp`State`LinearSolveChopPrecisionVal] &;
          If[Length[csFreedom] > 0,
            DiffExp`Utilities`PrintWarning["Not enough boundary data was provided for integral(s): ", intind, " at epsilon order ", epsord, "."];
            DiffExp`Utilities`PrintWarning["Introducing free parameters: ", Table[Subscript[Global`c, epsord, intind, i], {i, Length@csFreedom}]][1];
            csol = csol + Sum[csFreedom[[nullvecind]] Subscript[Global`c, epsord, intind, nullvecind], {nullvecind, Length@csFreedom}];
            TurnOffPade[];
          ];

          CsReps = Thread[cIndices -> csol] // DiffExp`Utilities`PChop;

          fGeneral = fGeneral /. CsReps // DiffExp`SeriesOps`SExpand;

          ,

          If[Length[cIndices] > 0,
            DiffExp`Utilities`PrintWarning["Not enough boundary data was provided for integral(s): ", intind, " at epsilon order ", epsord, "."];
            DiffExp`Utilities`PrintWarning["Introducing free parameters: ", cIndices /. Subscript[c, i_] :> Subscript[Global`c, epsord, intind, i]][1];

            TurnOffPade[];

            fGeneral = fGeneral /. Subscript[c, i_] :> Subscript[Global`c, epsord, intind, i] // DiffExp`SeriesOps`SExpand;

          ];
        ];

        ,

        (* No boundary terms provided. Relabeling for current epsilon order *)
        fGeneral = fGeneral /. Subscript[c, i_] :> Subscript[Global`c, epsord, intind, i] // DiffExp`SeriesOps`SExpand;
      ];

      If[MemberQ[DiffExp`State`CurrCrosscheckFlags, "GeneralSolution"] === True,
        DiffExp`Utilities`PrintDebug["General solution found. Cross-checking with differential equations."][1];
        CrossC = DiffExp`SeriesOps`SD[fGeneral, DiffExp`Symbols`x] - DiffExp`State`DEqnMatricesExpanded[line][0][[intind, intind]] . fGeneral - bVec + O[DiffExp`Symbols`x]^DiffExp`State`ICrossCheckVerifyResultOrder // (DiffExp`Utilities`CPChop@*DiffExp`SeriesOps`SExpand);
        DiffExp`Utilities`PrintDebug["Found: ", CrossC + O[DiffExp`Symbols`x]^DiffExp`State`ICrossCheckPrintResultOrder // DiffExp`SeriesOps`SN][3];
        If[
          !(SameQ @@ Append[CrossC // Normal // Flatten, 0])
          ,
          DiffExp`Utilities`ReportError["Cross-check failed"];
        ];
      ];

      NewResults = Thread[Table[{iind, epsord}, {iind, intind}] -> (DiffExp`Utilities`PChop@*DiffExp`SeriesOps`SExpand)[fGeneral]];

      AssociateTo[
        IntegrationData,
        NewResults
      ];

      DiffExp`State`BenchmarkData["Segments"][line // N]["ComputationTimes"]["Integrals"][intind][epsord] = AbsoluteTime[] - DiffExp`State`BenchmarkData["Segments"][line // N]["ComputationTimes"]["Integrals"][intind][epsord];

      , {intind, DiffExp`State`IntegrationSequence}
    ];

    DiffExp`State`BenchmarkData["Segments"][line // N]["ComputationTimes"][epsord] = AbsoluteTime[] - DiffExp`State`BenchmarkData["Segments"][line // N]["ComputationTimes"][epsord];

    , {epsord, 0, DiffExp`State`EpsilonOrderVal}];

  GiveMultivaluedError[] := (
    If[(DiffExp`State`FEC["AbortOnAnalyticContinuationFail"] === False),
      DiffExp`Utilities`PrintWarning["After fixing the boundary conditions, the solutions on the current line segment contain multivalued functions. However, the current point is not recognized as a branch point in the configuration, so DiffExp can't perform the analytic continuation! Please add a prescription for the analytic continuation using the option \"DeltaPrescriptions\"."];
      ,
      DiffExp`Utilities`ReportError["After fixing the boundary conditions, the solutions on the current line segment contain multivalued functions. However, the current point is not recognized as a branch point in the configuration, so DiffExp can't perform the analytic continuation! Please add a prescription for the analytic continuation using the option \"DeltaPrescriptions\".", False];

      DiffExp`State`MultivaluedFail = True;
    ];
  );

  If[MemberQ[DiffExp`State`CurrCrosscheckFlags, "SingularityCheck"] === True,
    TmpSolutionsNormal = Values[IntegrationData] // Normal // Chop[#, 10^-DiffExp`State`ICheckMultivaluedChop] &;
    LogsPresent = DiffExp`Utilities`DependsQ[TmpSolutionsNormal, DiffExp`Symbols`Logx | DiffExp`Symbols`Logx^_];
    AlgebraicRootsPresent = DiffExp`Utilities`DependsQ[TmpSolutionsNormal, (DiffExp`Symbols`x^b_ /; Denominator[b] > 1)];

    If[!DiffExp`State`CurrentSingularityHasIDeltaPrescription && (AlgebraicRootsPresent || LogsPresent),
      GiveMultivaluedError[];
    ];

    If[DiffExp`State`CurrentSingularityHasIDeltaPrescription && DiffExp`State`CurrentSingularityWasAddedFromSquareRoot && LogsPresent,
      GiveMultivaluedError[];
    ];
  ];

  If[(DiffExp`State`AnalyticContinuationFailed === True) && !MemberQ[opts, "TransportToCall"],
    DiffExp`Utilities`PrintWarning["Could not transfer i\[Delta]-prescriptions to line parameter on current segment. The results may only be valid in the direction (x > 0 or x < 0) where the boundary conditions were given."];
  ];

  If[MemberQ[opts, "TransportToCall"],
    Table[IntegrationData[{ind, epsord}], {ind, DiffExp`State`NumIntegrals}, {epsord, 0, DiffExp`State`EpsilonOrderVal}]
    ,
    DiffExp`State`BenchmarkData["Segments"]["ComputationTime"] = AbsoluteTime[] - DiffExp`State`BenchmarkData["Segments"]["ComputationTime"];

    DiffExp`MatrixLoading`ClearMatrices[];
    Table[DiffExp`SeriesOps`ApplyAnalyticContinuation[IntegrationData[{ind, epsord}]] // DiffExp`AnalyticContinuation`Project\[Theta]s, {ind, DiffExp`State`NumIntegrals}, {epsord, 0, DiffExp`State`EpsilonOrderVal}]
  ]
];

(* TransportTo function - transports boundary conditions *)
TransportTo[bcs2_List, line2_Association | line2_List, to2 : _?NumericQ : 1, SaveExpansions : _?BooleanQ : False, SampleAtList_List : {}] := Module[
  {line = If[line2[[0]] === List, line2 // Association // KeySort, line2 // KeySort], LineRat, ToRat, to, Tmp, bcs, FixAt, MySingularities, MySingularitiesImaginary, MySingularitiesRelevant, SingularitySegments, PoleIntervals, CurrLine, CurrLineNoMobius = Null, CurrIntegrated, CurrIntervalCurrLine, CurrIntervalLine, Done = False, Tmp2, CurrEvalPoint, CurrEvalPointCurrLine, CurrEval, Currbcs, AllIntegrationData = {}, MyCenter, EvaluateCurrPoint, NextIsPole = False, SegmentCounter = 1, TmpRelateLines, InterSec, CurrIntervalLinePos, CurrIntervalLineNeg, CurrEvalError, CurrIntegratedError, CurrbcsError, PrintError, CurrError, CurrErrorAcc = 0, CurrErrorAccs = ConstantArray[0, {DiffExp`State`NumIntegrals, DiffExp`State`EpsilonOrderVal + 1}], FixWithin, SegmentsToIntegrate, UpdateMatrixExpansionError, TimeStart, TimeStart0, LineReturn, FailedLine, ExpansionsIndeterminates = {}, bcsprev, CurrStatusBackup, BoundaryFixPoint, CurrEvalErrorEx, CurrEvalError1, CurrEvalError2, CurrEvalAtBoundaryFixPoint, CurrEvalEx, AllSegmentsPredivision, TmpFile, Ses, CompressedTermForExport, CompressedTermForExportFN, CurrLineLR, FullLineLR, TmpTmp, PoleIntervals1, RepeatingSegment, LastEvaluation, LastSavedData, LastLine, ExpansionOrders, DigitsNeeded},

  DiffExp`State`BenchmarkData = Association[];
  DiffExp`State`BenchmarkData["TimeStart"] = AbsoluteTime[];

  DiffExp`State`MultivaluedFail = False;

  If[Or @@ (!DiffExp`Utilities`DependsQ[Keys@line, #] & /@ DiffExp`State`ExternalScalesVal),
    DiffExp`Utilities`ReportError["The point/line in the second argument does not fix all kinematic invariants and masses!"];
  ];

  DiffExp`State`ExpansionOrderVal = DiffExp`State`FEC[ExpansionOrder];
  DigitsNeeded = DiffExp`State`FEAccuracyGoal + DiffExp`State`ISafetyDigits;
  TimeStart0 = AbsoluteTime[];

  (* To deal with the output of SaveExpansions = True *)
  If[MatchQ[bcs2, {{a_Association, __}, _}],
    bcs = bcs2[[1]];,
    bcs = bcs2;
  ];

  bcs[[1]] = bcs[[1]] // KeySort;

  (* Increase precision of arguments. *)
  If[Accuracy[line] < DiffExp`State`FEWorkingPrecision,
    line = line // Normal // SetPrecision[#, 2 * DiffExp`State`FEWorkingPrecision] & // Association;
    DiffExp`Utilities`PrintWarning["Accuracy of the line/point is lower than the working precision. The precision has been artificially increased."];
  ];

  If[DiffExp`Utilities`IsPoint[line],
    If[DiffExp`Utilities`IsPoint[bcs[[1]]],
      line = Merge[{bcs[[1]], line}, Expand[#[[1]] (1 - DiffExp`Symbols`x) + DiffExp`Symbols`x #[[2]]] &];
      ,
      Tmp = DiffExp`LineSegmentation`RelateLines[bcs[[1]], line, True];
      If[Tmp === False, DiffExp`Utilities`ReportError["Endpoint does not lie on same line as the boundary conditions."]];
      line = bcs[[1]] /. DiffExp`Symbols`x -> DiffExp`Symbols`x Tmp;
    ];
    to = 1;
    ,
    to = to2;

    If[Accuracy[to] < DiffExp`State`FEWorkingPrecision,
      to = to // SetPrecision[#, 2 * DiffExp`State`FEWorkingPrecision] &;
      DiffExp`Utilities`PrintWarning["Accuracy of the endpoint is lower than the working precision. The precision has been artificially increased."];
    ];
  ];

  ExpansionsIndeterminates = DeleteCases[bcs[[2]] // System`Variables, "?" | DiffExp`Symbols`x | DiffExp`Symbols`Logx];
  If[Length[ExpansionsIndeterminates] > 0 && DiffExp`State`FEC[UsePade] === True,
    DiffExp`Utilities`PrintInfo["The use of Pad\[EAcute] approximants is not possible when transporting boundary conditions with indeterminate coefficients. Pade approximants will be turned off."][1];
    DiffExp`State`DiffExpConfiguration[UsePade] = False;
  ];

  (* Check whether the line lies on a singularity of the differential equations. *)
  If[Length[Tmp = Flatten[Position[Factor[DiffExp`State`MatricesIrreducibleFactors /. line], 0]]] > 0,
    DiffExp`Utilities`ReportError["The line lies on a singularity of the differential equations. The vanishing factors are: ", DiffExp`State`MatricesIrreducibleFactors[[Tmp]], "."];
  ];

  (* Will abort if something is wrong with the boundary conditions. *)
  {bcs, FixAt} = DiffExp`LineSegmentation`CheckBoundaryConditionsAndReparametrize[bcsprev = bcs, line];

  If[Length[bcsprev] > 2,
    CurrErrorAcc = bcsprev[[3]] // Abs // Max;
    CurrErrorAccs = bcsprev[[3]];
    ,
    Quiet[
      CurrErrorAccs = 10^-(MapAt[Accuracy, bcsprev[[2]], {All, All}] // SetPrecision[#, DiffExp`State`FEWorkingPrecision] &);
      CurrErrorAcc = CurrErrorAccs // Max;
    ];
  ];

  If[!((CurrErrorAccs // Dimensions // Last) === DiffExp`State`EpsilonOrderVal + 1),
    CurrErrorAccs = PadRight[#, DiffExp`State`EpsilonOrderVal + 1, 10^-DiffExp`State`FEWorkingPrecision] & /@ CurrErrorAccs;
  ];

  DiffExp`Utilities`PrintInfo["Transporting boundary conditions along ", line // Normal // N // Association, " from x = ", FixAt // N, " to x = ", to // N][1];

  DiffExp`State`BenchmarkData["BasicPreprocessing"] = AbsoluteTime[] - DiffExp`State`BenchmarkData["TimeStart"];
  DiffExp`State`BenchmarkData["MatrixPreprocessing"] = AbsoluteTime[];

  DiffExp`Utilities`PrintInfo["Preparing differential equations along current line.."][1];
  DiffExp`MatrixLoading`PrepareMatricesFactored[line];

  DiffExp`MatrixLoading`InitializeIntegrationSequence[line];

  {MySingularities, MySingularitiesImaginary} = DiffExp`LineSegmentation`FindMatrixSingularities[line, True, {FixAt, to}];

  MySingularitiesRelevant = Select[MySingularities, FixAt <= # <= to || to <= # <= FixAt &];
  If[DiffExp`State`FEC[UseMobius] === True,
    SingularitySegments = {
      #,
      DiffExp`Mobius`GetLineRescaled[line, #, {MySingularities, MySingularitiesImaginary}],
      DiffExp`Mobius`GetLineRescaled[line, #, {MySingularities, MySingularitiesImaginary}, True]
    } & /@ MySingularitiesRelevant;
    ,
    SingularitySegments = {#, DiffExp`Mobius`GetLineRescaled[line, #, {MySingularities, MySingularitiesImaginary}]} & /@ MySingularitiesRelevant;
  ];

  DiffExp`Utilities`PrintInfo["Possible singularities along line at positions ", DeleteCases[SingularitySegments[[All, 1]], \[Infinity] | -\[Infinity]] // N, "."][1];

  DiffExp`State`BenchmarkData["MatrixPreprocessing"] = AbsoluteTime[] - DiffExp`State`BenchmarkData["MatrixPreprocessing"];
  DiffExp`State`BenchmarkData["PoleIntervals"] = AbsoluteTime[];

  If[DiffExp`State`FEC[SegmentationStrategy] === "Dynamic",
    DiffExp`Utilities`PrintInfo["Determining intervals around possible singularities."][1];
  ];

  PoleIntervals = Table[
    If[DiffExp`State`FEC[SegmentationStrategy] === "Dynamic",
      DiffExp`Utilities`PrintInfo["Expanding differential equations around x = ", sline[[1]] // N, "."][1];

      DiffExp`AnalyticContinuation`PrepareAnalyticContinuation[sline[[2]]];
      DiffExp`MatrixLoading`PrepareMatricesFrom[line, sline[[2]]];

      Tmp = {-#, #} &@DiffExp`LineSegmentation`GetMatricesPrecisionDistance[sline[[2]]];
    ];

    If[DiffExp`State`FEC[SegmentationStrategy] === "Predivision",
      Tmp = {-DiffExp`State`RadiusOfConvergenceVal/DiffExp`State`FEC[DivisionOrder], DiffExp`State`RadiusOfConvergenceVal/DiffExp`State`FEC[DivisionOrder]};
    ];

    Tmp = {sline, DiffExp`LineSegmentation`RelateLinesPoint[line, sline[[2]], #] & /@ Tmp};

    DiffExp`Utilities`PrintInfo["Expansion around x = ", sline[[1]] // N, " is valid within region x \[Element] [", Tmp[[2, 1]] // N, ", ", Tmp[[2, 2]] // N, "]."][2];

    Tmp
    , {sline, SingularitySegments}
  ];

  DiffExp`State`BenchmarkData["PoleIntervals"] = AbsoluteTime[] - DiffExp`State`BenchmarkData["PoleIntervals"];

  (* Main loop. *)
  If[DiffExp`Utilities`DependsQ[# === 0 & /@ (SingularitySegments // DiffExp`Utilities`PChop), True],
    DiffExp`Utilities`PrintDebug["First expansion is at singularity."][1];
  ];

  DiffExp`State`BenchmarkData["SegmentCounting"] = AbsoluteTime[];

  (* Duplicate of the code block below, without print statements and integrations. This counts the number of segments. *)
  AllSegmentsPredivision = {};
  If[DiffExp`State`FEC[SegmentationStrategy] === "Predivision",
    DiffExp`Utilities`PrintInfo["Analyzing integration segments."][1];
    PoleIntervals1 = PoleIntervals;
    MyCenter = FixAt;
    CurrLine = DiffExp`Mobius`GetLineRescaled[line, FixAt, {MySingularities, MySingularitiesImaginary}];
    AppendTo[AllSegmentsPredivision, CurrLine];

    Done = False;
    While[!Done,
      PoleIntervals = Select[PoleIntervals, !(#[[1, 1]] === MyCenter) &];

      CurrIntervalCurrLine = {-DiffExp`State`RadiusOfConvergenceVal/DiffExp`State`FEC[DivisionOrder], DiffExp`State`RadiusOfConvergenceVal/DiffExp`State`FEC[DivisionOrder]};
      CurrIntervalLine = DiffExp`LineSegmentation`RelateLinesPoint[line, CurrLine, #] & /@ CurrIntervalCurrLine;

      If[IntervalContainsQ[CurrIntervalLine, to],
        CurrEvalPoint = to;
        CurrEvalPointCurrLine = DiffExp`LineSegmentation`RelateLinesPoint[CurrLine, line, CurrEvalPoint];

        Done = True;

        ,
        SegmentCounter += 1;

        If[to > FixAt,
          Tmp = Select[PoleIntervals, #[[1, 1]] > MyCenter &];
          ,
          Tmp = Select[PoleIntervals, #[[1, 1]] <= MyCenter &];
        ];
        Tmp2 = {#[[1]], IntervalOverlapQ[CurrIntervalLine, #[[2]]]} & /@ Tmp;
        Tmp2 = Flatten[Position[Tmp2[[All, 2]], True]];
        NextIsPole = Length[Tmp2] > 0;

        If[NextIsPole,
          Tmp = Tmp[[Tmp2[[1]]]];

          If[to > FixAt,
            FixWithin = IntervalIntersec[IntervalIntersec[{MyCenter, Tmp[[1, 1]]}, Tmp[[2]]], CurrIntervalLine];
            CurrEvalPoint = (FixWithin[[1]] + FixWithin[[2]])/2;
            ,

            FixWithin = IntervalIntersec[IntervalIntersec[{Tmp[[1, 1]], MyCenter}, Tmp[[2]]], CurrIntervalLine];
            CurrEvalPoint = (FixWithin[[1]] + FixWithin[[2]])/2;
          ];

          CurrEvalPointCurrLine = DiffExp`LineSegmentation`RelateLinesPoint[CurrLine, line, CurrEvalPoint];

          MyCenter = Tmp[[1, 1]];
          CurrLine = Tmp[[1, 2]];
          If[DiffExp`State`FEC[UseMobius] === True, CurrLineNoMobius = Tmp[[1, 3]]];
          ,

          (* Finite point *)
          If[to > FixAt,
            CurrEvalPoint = CurrIntervalLine[[2]],
            CurrEvalPoint = CurrIntervalLine[[1]]
          ];
          CurrEvalPointCurrLine = DiffExp`LineSegmentation`RelateLinesPoint[CurrLine, line, CurrEvalPoint];

          If[DiffExp`State`FEC[SegmentationStrategy] === "Dynamic",
            MyCenter = CurrEvalPoint;
          ];

          If[DiffExp`State`FEC[SegmentationStrategy] === "Predivision",
            If[to > FixAt,
              MyCenter = DiffExp`Mobius`FindNextCenterPointL[CurrEvalPoint, MySingularities],
              MyCenter = DiffExp`Mobius`FindNextCenterPointR[CurrEvalPoint, MySingularities]
            ];
          ];

          CurrLine = DiffExp`Mobius`GetLineRescaled[line, MyCenter, {MySingularities, MySingularitiesImaginary}];
        ];
      ];

      AppendTo[AllSegmentsPredivision, CurrLine];

      PoleIntervals = PoleIntervals1;
    ];

    DiffExp`Utilities`PrintInfo["Segments to integrate: ", SegmentsToIntegrate = SegmentCounter, "."][1];
  ];

  DiffExp`State`BenchmarkData["SegmentCounting"] = AbsoluteTime[] - DiffExp`State`BenchmarkData["SegmentCounting"];
  DiffExp`State`BenchmarkData["AllPreprocessing"] = AbsoluteTime[] - DiffExp`State`BenchmarkData["TimeStart"];

  SegmentCounter = 1;
  ExpansionOrders = {DiffExp`State`ExpansionOrderVal};

  If[DiffExp`State`FEC[SegmentationStrategy] === "Predivision",
    DigitsNeeded = DiffExp`State`FEAccuracyGoal + Ceiling[Log10[SegmentsToIntegrate]] + DiffExp`State`ISafetyDigits;
  ];

  UpdateMatrixExpansionError[] := Block[{CurrErrorTerms, CurrErrorLeft, CurrErrorRight, FixAtLineSegment},
    CurrErrorTerms = DiffExp`LineSegmentation`GetLargestTerm[CurrLine];
    FixAtLineSegment = DiffExp`LineSegmentation`GetMatchingPoint[CurrLine, Currbcs[[1]]];

    If[FixAtLineSegment === 0,
      CurrErrorLeft = 0;,
      CurrErrorLeft = CurrErrorTerms /. (DiffExp`State`AnalyticContinuationReplacements /. If[FixAtLineSegment >= 0, {DiffExp`Symbols`\[Theta]p -> 1, DiffExp`Symbols`\[Theta]m -> 0}, {DiffExp`Symbols`\[Theta]p -> 0, DiffExp`Symbols`\[Theta]m -> 1}]) /. DiffExp`Symbols`Logx -> Log[DiffExp`Symbols`x] /. DiffExp`Symbols`x -> FixAtLineSegment // Abs // N[#, DiffExp`State`FEWorkingPrecision] & // Max;
    ];
    CurrErrorRight = CurrErrorTerms /. (DiffExp`State`AnalyticContinuationReplacements /. If[CurrEvalPointCurrLine >= 0, {DiffExp`Symbols`\[Theta]p -> 1, DiffExp`Symbols`\[Theta]m -> 0}, {DiffExp`Symbols`\[Theta]p -> 0, DiffExp`Symbols`\[Theta]m -> 1}]) /. DiffExp`Symbols`Logx -> Log[DiffExp`Symbols`x] /. DiffExp`Symbols`x -> CurrEvalPointCurrLine // Abs // N[#, DiffExp`State`FEWorkingPrecision] & // Max;
    CurrError = Max[CurrErrorLeft, CurrErrorRight]
  ];

  EvaluateCurrPoint[] := Block[{ErrorSufficient, FixAtLineSegment, NumExpansions, DigitsHave, DeltaDigits, DEqnMatricesExpandedBackup, CurrErrorBackup, MyCounter}, (

    If[DiffExp`State`FEC[SegmentationStrategy] == "Predivision" && KeyExistsQ[DiffExp`State`DiffExpConfiguration, System`AccuracyGoal] && DiffExp`State`DiffExpConfiguration["AccuracyGoalValidate"] === "Before",
      If[NumericQ[DiffExp`State`DiffExpConfiguration[System`AccuracyGoal]],
        ErrorSufficient = False;
        FixAtLineSegment = DiffExp`LineSegmentation`GetMatchingPoint[CurrLine, Currbcs[[1]]];
        DiffExp`State`ExpansionOrderVal = Ceiling[Mean[ExpansionOrders]];

        DiffExp`Utilities`PrintInfo["Determining expansion order for given accuracy goal."][2];
        NumExpansions = 1;
        While[!ErrorSufficient,
          DiffExp`Utilities`PrintInfo["Expanding differential equations at order: ", DiffExp`State`ExpansionOrderVal][2];
          KeyDropFrom[DiffExp`State`DEqnMatricesExpanded, CurrLine];
          DiffExp`MatrixLoading`PrepareMatricesExpanded[CurrLine];

          UpdateMatrixExpansionError[];

          DigitsHave = -Log10[CurrError];
          DeltaDigits = DigitsNeeded - DigitsHave;

          If[DeltaDigits < 0,
            If[NumExpansions == 1 && -DeltaDigits > DiffExp`State`IDigitsSurplusDecreaseExpansionOrder,
              MyCounter = 1;
              While[-DeltaDigits > DiffExp`State`IDigitsSurplusDecreaseExpansionOrder && DiffExp`State`ExpansionOrderVal - DiffExp`State`IExpansionOrderDecrease >= DiffExp`State`IMinExpansionOrder,
                If[MyCounter > 1,
                  DiffExp`Utilities`PrintInfo["Reducing expansion order..: ", DiffExp`State`ExpansionOrderVal][2];
                ];
                DEqnMatricesExpandedBackup = DiffExp`State`DEqnMatricesExpanded[line];
                CurrErrorBackup = CurrError;

                DiffExp`State`ExpansionOrderVal -= DiffExp`State`IExpansionOrderDecrease;
                DiffExp`State`DEqnMatricesExpanded[line] = # + O[DiffExp`Symbols`x]^DiffExp`State`ExpansionOrderVal & /@ DiffExp`State`DEqnMatricesExpanded[line];

                UpdateMatrixExpansionError[];
                DigitsHave = -Log10[CurrError];
                DeltaDigits = DigitsNeeded - DigitsHave;

                MyCounter++;
              ];

              DiffExp`State`ExpansionOrderVal += DiffExp`State`IExpansionOrderDecrease;
              DiffExp`State`DEqnMatricesExpanded[line] = DEqnMatricesExpandedBackup;
              CurrError = CurrErrorBackup;

              DiffExp`Utilities`PrintInfo["Error of matrix expansions: ", CurrError // N][2];
              If[MyCounter == 2,
                DiffExp`Utilities`PrintInfo["Precision reached at order: ", DiffExp`State`ExpansionOrderVal][2];
              ];

              ,

              DiffExp`Utilities`PrintInfo["Error of matrix expansions: ", CurrError // N][2];
              DiffExp`Utilities`PrintInfo["Precision reached at order: ", DiffExp`State`ExpansionOrderVal][2];

            ];
            AppendTo[ExpansionOrders, DiffExp`State`ExpansionOrderVal];
            If[Length[ExpansionOrders] > DiffExp`State`IExpansionOrdersAveraging,
              ExpansionOrders = Delete[ExpansionOrders, 1];
            ];
            ErrorSufficient = True;
            ,
            DiffExp`Utilities`PrintInfo["Error of matrix expansions: ", CurrError // N][2];
            DiffExp`State`ExpansionOrderVal += DiffExp`State`IExpansionOrderIncrease;
            NumExpansions++;
          ];
        ];
      ];
    ];

    CurrIntegrated = IntegrateSystem[Currbcs, CurrLine, {"TransportToCall"}];

    If[SaveExpansions === True,
      TmpRelateLines = DiffExp`Symbols`x -> DiffExp`LineSegmentation`RelateLines[CurrLine, line];
      AppendTo[AllIntegrationData, {
        CurrLine, (* Current line *)
        TmpRelateLines, (* Change of line parameter from CurrLine to line *)
        If[to > FixAt, Identity, Reverse]@{Limit[DiffExp`LineSegmentation`RelateLines[line, Currbcs[[1]]], DiffExp`Symbols`x -> 0], CurrEvalPoint}, (*  Expansions gives results on line between x and y *)
        {Limit[DiffExp`LineSegmentation`RelateLines[CurrLine, Currbcs[[1]]], DiffExp`Symbols`x -> 0], CurrEvalPointCurrLine}, (*  Expansions gives results on CurrLine between x and y *)
        If[DiffExp`State`FEC["SaveExpansionsCompress"] === True,
          CompressedTermForExport = (DiffExp`SeriesOps`ApplyAnalyticContinuation[CurrIntegrated] // DiffExp`AnalyticContinuation`Project\[Theta]s) // If[KeyExistsQ[DiffExp`State`FEC, "SaveExpansionsOrder"], (# + O[DiffExp`Symbols`x]^DiffExp`State`FEC["SaveExpansionsOrder"]) &, Identity] // Compress (* And the expansion data itself *);
          If[!DiffExp`State`FEC["SaveExpansionsCompressDirectory"] === "?",
            CompressedTermForExportFN = FileNameJoin[{DiffExp`State`FEC["SaveExpansionsCompressDirectory"], Hash[CompressedTermForExport, "SHA256", "HexString"] <> ".m"}];
            Export[
              CompressedTermForExportFN,
              CompressedTermForExport
            ];
            CompressedTermForExportFN
            ,
            CompressedTermForExport
          ]
          ,
          (DiffExp`SeriesOps`ApplyAnalyticContinuation[CurrIntegrated] // DiffExp`AnalyticContinuation`Project\[Theta]s) // If[KeyExistsQ[DiffExp`State`FEC, "SaveExpansionsOrder"], (# + O[DiffExp`Symbols`x]^DiffExp`State`FEC["SaveExpansionsOrder"]) &, Identity] (* And the expansion data itself *)
        ]
      }];
    ];

    (* Write out points if provided *)
    If[!SampleAtList === {},
      CurrLineLR = {Limit[DiffExp`LineSegmentation`RelateLines[CurrLine, Currbcs[[1]]], DiffExp`Symbols`x -> 0], CurrEvalPointCurrLine};
      FullLineLR = If[to > FixAt, Identity, Reverse]@{Limit[DiffExp`LineSegmentation`RelateLines[line, Currbcs[[1]]], DiffExp`Symbols`x -> 0], CurrEvalPoint};

      (* SampleAtList has the form {{xval, filename}, ...}*)
      Do[
        DiffExp`Utilities`PrintInfo["Evaluating and exporting at x = ", MyPoint // N, "."][2];

        Export[
          FileNameJoin[{
            SampleAtList[[2]],
            (line /. DiffExp`Symbols`x -> MyPoint // Hash[Normal[KeySort[#]], "SHA256", "HexString"] &) <> ".m"
          }
          ],

          {
            line /. DiffExp`Symbols`x -> MyPoint // Normal // Association,
            DiffExp`Pade`SEval[
              CurrIntegrated,
              DiffExp`LineSegmentation`RelateLinesPoint[CurrLine, line, MyPoint]
            ] // If[Length[SampleAtList] === 3, N[#, SampleAtList[[3]]] &, Identity] // If[DiffExp`State`FEC["SaveExpansionsCompress"] === True, Compress, Identity]
          }
        ];

        , {MyPoint, (* All points lying in the current segment *)Select[SampleAtList[[1]], FullLineLR[[1]] <= # <= FullLineLR[[2]] &]}
      ];
    ];

    (* If MultivaluedFail is True, then multivalued functions were encountered but no singularity was expected. We then abort the computation early. *)
    If[!DiffExp`State`MultivaluedFail,

      DiffExp`Utilities`PrintInfo["Evaluating at x = ", CurrEvalPoint // N][1];
      Global`DebugData = CurrIntegrated;
      (* Either returns the Pade approximant, or the series with Normal applied *)
      CurrEvalEx = DiffExp`Pade`SEval1[CurrIntegrated];
      CurrEval = DiffExp`Pade`SEval2[CurrEvalEx, CurrEvalPointCurrLine];
      If[!DiffExp`State`FEC["EstimateError"] === False,
        Switch[DiffExp`State`FEC["EstimateError"],
          "Fast",
          (* Obtain the point at which the boundary conditions were fixed. *)
          BoundaryFixPoint = DiffExp`LineSegmentation`GetMatchingPoint[CurrLine, Currbcs[[1]]];
          (* We compute the error both at the position where the boundary conditions were fixed, and at the evaluation point of the current line segment.*)
          CurrEvalErrorEx = DiffExp`Pade`SEval1[CurrIntegrated // DiffExp`SeriesOps`DecreaseSeriesOrderBy[#, DiffExp`State`ICurrEvalErrorSeriesDecrease] &];
          CurrEvalError1 = DiffExp`Pade`SEval2[CurrEvalErrorEx, CurrEvalPointCurrLine];
          (* We avoid the point 0, because there might be logarithms, and the error should be manifestly zero. *)
          If[!BoundaryFixPoint === 0,
            CurrEvalError2 = DiffExp`Pade`SEval2[CurrEvalErrorEx, BoundaryFixPoint];
            CurrEvalAtBoundaryFixPoint = DiffExp`Pade`SEval2[CurrEvalEx, BoundaryFixPoint];
            ,
            (* As a small hack, we set the same values as at the endpoint. *)
            CurrEvalError2 = CurrEvalError1;
            CurrEvalAtBoundaryFixPoint = CurrEval;
          ];
        ];
      ];

      DiffExp`Utilities`PrintInfo[
        If[DiffExp`State`FEC[SegmentationStrategy] === "Predivision",
          "Integrated segment " <> ToString[SegmentCounter - 1] <> " out of " <> ToString[SegmentsToIntegrate] <> " in ",
          "Integrated segment in "
        ],
        AbsoluteTime[] - TimeStart // N,
        " seconds."
      ][1];

      ,

      (* Computation should be aborted and last integrated system should be returned *)
      FailedLine = CurrLine;
      CurrEval = CurrIntegrated;
      CurrEvalError = Table[0, {ii, DiffExp`State`NumIntegrals}, {jj, 0, DiffExp`State`EpsilonOrderVal}];

    ];

  )];

  ComputeErrorsPerIndeterminate[aaa_, bbb_, ExpIndets_] := Module[{TmpErrors},
    TmpErrors = Table[DiffExp`SeriesOps`LogxCoeffNS[aaa - bbb, logxord], {logxord, 0, DiffExp`State`IMaxLogOrder}];
    TmpErrors = Table[(

      Flatten[
        Append[
          Table[
            Coefficient[TmpErrors[[All, ii, jj]], var]
            , {var, ExpIndets}]
          ,
          TmpErrors[[All, ii, jj]] /. (# -> 0 & /@ ExpIndets)
        ]
      ]

    ) // Abs // Max, {ii, (TmpErrors // Dimensions)[[2]]}, {jj, (TmpErrors // Dimensions)[[3]]}] // N
  ];

  PrintError[] := Block[{TmpErrors, TmpErrors1, TmpErrors2, ExpIndets},

    If[!DiffExp`State`MultivaluedFail,

      ExpIndets = DeleteCases[Currbcs[[2]] // System`Variables, "?" | DiffExp`Symbols`x | DiffExp`Symbols`Logx];

      TmpErrors1 = ComputeErrorsPerIndeterminate[CurrEvalError1, CurrEval, ExpIndets];
      TmpErrors2 = ComputeErrorsPerIndeterminate[CurrEvalError2, CurrEvalAtBoundaryFixPoint, ExpIndets];
      TmpErrors = Table[Max[{TmpErrors1[[ii, jj]], TmpErrors2[[ii, jj]]}], {ii, First@Dimensions@TmpErrors1}, {jj, Last@Dimensions@TmpErrors1}];

      CurrError = Flatten[TmpErrors] // Abs // Max;

      Switch[DiffExp`State`FEC["EstimateError"],
        "Fast",
        CurrErrorAcc = CurrErrorAcc + CurrError;
        CurrErrorAccs = CurrErrorAccs + TmpErrors;
      ];

      Which[DiffExp`State`FEC["EstimateError"] == "Fast",
        DiffExp`Utilities`PrintInfo["Current segment error estimate: ", CurrError][1];
        DiffExp`Utilities`PrintInfo["Total error estimate: ", CurrErrorAcc][1];
      ];

      If[CurrError > 1,
        DiffExp`Utilities`ReportError["The reported error is very large. This likely indicates a numerical instability."];
      ];

    ];

  ];

  MyCenter = FixAt;
  CurrLine = DiffExp`Mobius`GetLineRescaled[line, FixAt, {MySingularities, MySingularitiesImaginary}];
  If[DiffExp`State`FEC[UseMobius] === True, CurrLineNoMobius = DiffExp`Mobius`GetLineRescaled[line, FixAt, {MySingularities, MySingularitiesImaginary}, True]];
  Currbcs = CurrbcsError = bcs;

  DiffExp`State`BenchmarkData["Segments"] = Association[];

  Done = False;
  While[!Done,
    DiffExp`State`BenchmarkData["Segments"][CurrLine // N] = Association[];
    DiffExp`State`BenchmarkData["Segments"][CurrLine // N]["ComputationTime"] = TimeStart = AbsoluteTime[];

    CurrStatusBackup = {Currbcs, CurrbcsError, bcs, CurrLine, MyCenter, FixAt, CurrErrorAcc, CurrErrorAccs};

    DiffExp`Utilities`PrintInfo["Integrating segment: ", DiffExp`LineSegmentation`PrintMobiusNormalized /@ CurrLine, "."][1];
    If[NextIsPole && MemberQ[MySingularitiesImaginary, MyCenter],
      DiffExp`Utilities`PrintInfo["Current segment is centered at singularity."][1];
    ];

    LastEvaluation = Currbcs;
    LastSavedData = AllIntegrationData;
    LastLine = CurrLine;

    PoleIntervals = Select[PoleIntervals, !(#[[1, 1]] === MyCenter) &];

    DiffExp`State`AnalyticContinuationFailed = False;
    DiffExp`AnalyticContinuation`PrepareAnalyticContinuation[CurrLine];
    If[DiffExp`State`AnalyticContinuationFailed === True && (SegmentCounter > 1 || to < FixAt),
      If[SegmentCounter === 1 && to < FixAt,
        DiffExp`Utilities`ReportError["Can only perform the expansions for x > 0 along the current line segment. Please integrate in the positive line direction."];
        ,
        DiffExp`Utilities`ReportError["Analytic continuation failed. Please separate out the singularities by transporting along a different line."];
      ];
    ];

    DiffExp`State`BenchmarkData["Segments"][LastLine // N]["MatrixExpansion"] = AbsoluteTime[];
    If[DiffExp`State`FEC[SegmentationStrategy] === "Dynamic",
      DiffExp`MatrixLoading`PrepareMatricesFrom[line, CurrLine];
      CurrIntervalCurrLine = {-#, #} &@DiffExp`LineSegmentation`GetMatricesPrecisionDistance[CurrLine];
    ];
    If[DiffExp`State`FEC[SegmentationStrategy] === "Predivision",
      DiffExp`MatrixLoading`PrepareMatricesFrom1[line, CurrLine];
      CurrIntervalCurrLine = {-DiffExp`State`RadiusOfConvergenceVal/DiffExp`State`FEC[DivisionOrder], DiffExp`State`RadiusOfConvergenceVal/DiffExp`State`FEC[DivisionOrder]};
    ];
    DiffExp`State`BenchmarkData["Segments"][LastLine // N]["MatrixExpansion"] = AbsoluteTime[] - DiffExp`State`BenchmarkData["Segments"][LastLine // N]["MatrixExpansion"];

    CurrIntervalLine = DiffExp`LineSegmentation`RelateLinesPoint[line, CurrLine, #] & /@ CurrIntervalCurrLine;
    DiffExp`Utilities`PrintInfo["Current line segment covers x \[Element] [", CurrIntervalLine[[1]] // N, ", ", CurrIntervalLine[[2]] // N, "]."][2];
    SegmentCounter += 1;

    If[IntervalContainsQ[CurrIntervalLine, to],
      CurrEvalPoint = to;
      CurrEvalPointCurrLine = DiffExp`LineSegmentation`RelateLinesPoint[CurrLine, line, CurrEvalPoint];
      EvaluateCurrPoint[];

      Currbcs = {line /. DiffExp`Symbols`x -> CurrEvalPoint, CurrEval};
      If[!DiffExp`State`FEC["EstimateError"] === False,
        CurrbcsError = {line /. DiffExp`Symbols`x -> CurrEvalPoint, CurrEvalError};
        PrintError[];
      ];

      Done = True;

      ,

      If[to > FixAt,
        Tmp = Select[PoleIntervals, #[[1, 1]] > MyCenter &];
        ,
        Tmp = Select[PoleIntervals, #[[1, 1]] <= MyCenter &];
      ];
      Tmp2 = {#[[1]], IntervalOverlapQ[CurrIntervalLine, #[[2]]]} & /@ Tmp;
      Tmp2 = Flatten[Position[Tmp2[[All, 2]], True]];
      NextIsPole = Length[Tmp2] > 0;

      If[NextIsPole,
        Tmp = Tmp[[Tmp2[[1]]]];

        If[to > FixAt,
          (* (Current line center, pole center), (poleinterval) *)
          FixWithin = IntervalIntersec[IntervalIntersec[{MyCenter, Tmp[[1, 1]]}, Tmp[[2]]], CurrIntervalLine];
          CurrEvalPoint = (FixWithin[[1]] + FixWithin[[2]])/2;

          ,

          (* (pole center, current line center), (poleinterval) *)
          FixWithin = IntervalIntersec[IntervalIntersec[{Tmp[[1, 1]], MyCenter}, Tmp[[2]]], CurrIntervalLine];
          CurrEvalPoint = (FixWithin[[1]] + FixWithin[[2]])/2;


        ];


        CurrEvalPointCurrLine = DiffExp`LineSegmentation`RelateLinesPoint[CurrLine, line, CurrEvalPoint];
        EvaluateCurrPoint[];

        MyCenter = Tmp[[1, 1]];
        CurrLine = Tmp[[1, 2]];
        If[DiffExp`State`FEC[UseMobius] === True, CurrLineNoMobius = Tmp[[1, 3]]];
        ,
        (* Finite point *)
        If[to > FixAt,
          CurrEvalPoint = CurrIntervalLine[[2]],
          CurrEvalPoint = CurrIntervalLine[[1]]
        ];
        CurrEvalPointCurrLine = DiffExp`LineSegmentation`RelateLinesPoint[CurrLine, line, CurrEvalPoint];
        EvaluateCurrPoint[];

        If[DiffExp`State`FEC[SegmentationStrategy] === "Dynamic",
          MyCenter = CurrEvalPoint;
        ];
        If[DiffExp`State`FEC[SegmentationStrategy] === "Predivision",
          If[to > FixAt,
            MyCenter = DiffExp`Mobius`FindNextCenterPointL[CurrEvalPoint, MySingularities],
            MyCenter = DiffExp`Mobius`FindNextCenterPointR[CurrEvalPoint, MySingularities]
          ];
        ];

        CurrLine = DiffExp`Mobius`GetLineRescaled[line, MyCenter, {MySingularities, MySingularitiesImaginary}];
        If[DiffExp`State`FEC[UseMobius] === True, CurrLineNoMobius = DiffExp`Mobius`GetLineRescaled[line, MyCenter, {MySingularities, MySingularitiesImaginary}, True]];
      ];

      Currbcs = {
        line /. DiffExp`Symbols`x -> CurrEvalPoint,
        CurrEval
      };

      If[!DiffExp`State`FEC["EstimateError"] === False,
        CurrbcsError = {line /. DiffExp`Symbols`x -> CurrEvalPoint, CurrEvalError};
        PrintError[];
      ];
    ];

    DiffExp`State`BenchmarkData["Segments"][LastLine // N]["ComputationTime"] = AbsoluteTime[] - DiffExp`State`BenchmarkData["Segments"][LastLine // N]["ComputationTime"];

    (* Analytic continuation failed. *)
    If[DiffExp`State`MultivaluedFail,
      DiffExp`Utilities`PrintWarning["The computation will be aborted, and the last line segment will be returned."];
      Done = True;
    ];

    RepeatingSegment = False;
    If[DiffExp`State`FEC[SegmentationStrategy] == "Predivision" && KeyExistsQ[DiffExp`State`DiffExpConfiguration, System`AccuracyGoal] && DiffExp`State`DiffExpConfiguration["AccuracyGoalValidate"] === "After",
      If[CurrError > 10^-DiffExp`State`FEC[System`AccuracyGoal],
        DiffExp`State`ExpansionOrderVal += DiffExp`State`IExpansionOrderIncrease2;
        DiffExp`Utilities`PrintInfo["The estimated error of the results is lower than the requested AccuracyGoal. The expansions will be repeated at the order ", DiffExp`State`ExpansionOrderVal, "."][1];
        (* Reload variables for the computation of the segment *)
        {Currbcs, CurrbcsError, bcs, CurrLine, MyCenter, FixAt, CurrErrorAcc, CurrErrorAccs} = CurrStatusBackup;
        If[SaveExpansions === True,
          AllIntegrationData = Delete[AllIntegrationData, -1];
        ];
        (* Drop the expansion matrices *)
        KeyDropFrom[DiffExp`State`DEqnMatricesExpanded, CurrLine];
        (* Decrease counter *)
        SegmentCounter -= 1;
        (* We are not done yet. *)
        Done = False;
        RepeatingSegment = True;
        ,
        (* If we had reset the value of expansion order, we put it back. *)
        If[DiffExp`State`ExpansionOrderVal != DiffExp`State`FEC[ExpansionOrder],
          DiffExp`State`ExpansionOrderVal = DiffExp`State`FEC[ExpansionOrder];
        ];
      ];
    ];

    If[!RepeatingSegment,
      (* Clear up some memory. *)
      DiffExp`MatrixLoading`ClearMatrices[CurrStatusBackup[[4]]];
    ];

  ];

  DiffExp`Utilities`PrintInfo["Finished integration of ", SegmentCounter - 1, " segments in ", AbsoluteTime[] - TimeStart0 // N, " seconds."][1];
  DiffExp`Utilities`PrintDebug["Performed ", SegmentCounter - SegmentsToIntegrate, " additional expansions due to changing the expansion order."][1];

  If[DiffExp`State`MultivaluedFail,
    LineReturn = FailedLine,
    LineReturn = line /. DiffExp`Symbols`x -> to;
  ];

  DiffExp`State`BenchmarkData["TimeEnd"] = AbsoluteTime[];
  DiffExp`State`BenchmarkData["ComputationTime"] = DiffExp`State`BenchmarkData["TimeEnd"] - DiffExp`State`BenchmarkData["TimeStart"];
  KeyDropFrom[DiffExp`State`BenchmarkData, {"TimeEnd", "TimeStart"}];

  If[SaveExpansions === True,
    If[!DiffExp`State`FEC["EstimateError"] === False,
      {{LineReturn, CurrEval, CurrErrorAccs}, AllIntegrationData},
      {{LineReturn, CurrEval}, AllIntegrationData}
    ]
    ,
    If[!DiffExp`State`FEC["EstimateError"] === False,
      {LineReturn, CurrEval, CurrErrorAccs},
      {LineReturn, CurrEval}
    ]
  ]
];

(* Helper functions *)
IntervalOverlapQ[intv1_, intv2_] := !(IntervalIntersection[Interval[intv1], Interval[intv2]] === Interval[]);
IntervalIntersec[intv1_, intv2_] := IntervalIntersection[Interval[intv1], Interval[intv2]][[1]];
IntervalContainsQ[intv_, point_] := intv[[1]] <= point <= intv[[2]];

End[];

EndPackage[];
