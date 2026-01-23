(* IntegrationStrategies/Default.m *)
(* SolveSimple and SolveDefault integration strategies *)

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
(* Falls back to VOPAlt if no valid pivot is found *)
SolveDefault[intind_, bVec_, line_, epsord_, BufferedDataIn_] := Module[
  {systemSize, HomogeneousEquation, MtildeMat, selectedPivot, pivotResult, NMat, Solns, Wronsk, WronskInv, FMat, FMatInv,
   GMat, BMat, cIndices, fGeneral, crossCheck, CurrInvWronskSolver,
   HomogeneousEquation2, MtildeMat2, NMat2, Solns2, Wronsk2, WronskInvPrime, wronskianProduct, MtildeInv, c,
   BufferedData = BufferedDataIn},

  If[epsord === 0 && !KeyExistsQ[BufferedData, intind],
    DiffExp`State`BenchmarkData["Segments"][line // N]["HomogeneousSolveAllPreprocessing"]["Integrals"][intind] = AbsoluteTime[];

    If[Length[intind] > 1,
      DiffExp`Utilities`PrintInfo["Combining differential equations: ", intind, " with automatic pivot selection."][3];
    ];
    systemSize = intind // Length;

    (* Try to find a valid pivot; fall back to VOPAlt if all pivots fail *)
    pivotResult = DiffExp`Wronskian`CombineDifferentialEquationsWithPivotSelection[
      If[DiffExp`State`FEC["HomogeneousSolve"] === "Expand",
        DiffExp`State`DEqnMatricesExpanded[line][0][[intind, intind]],
        DiffExp`State`DEqnMatricesFactored[line][0][[intind, intind]]
      ]
    ];

    (* Check if fallback to VOPAlt is needed *)
    If[pivotResult === $NeedsFallback,
      DiffExp`State`BenchmarkData["Segments"][line // N]["HomogeneousSolveAllPreprocessing"]["Integrals"][intind] = AbsoluteTime[] - DiffExp`State`BenchmarkData["Segments"][line // N]["HomogeneousSolveAllPreprocessing"]["Integrals"][intind];
      Return[SolveVOPAlt[intind, bVec, line, epsord, BufferedData]]
    ];

    {HomogeneousEquation, MtildeMat, selectedPivot} = pivotResult;
    DiffExp`Utilities`PrintInfo["Using pivot integral ", intind[[selectedPivot]], " (index ", selectedPivot, " of ", systemSize, ")."][3];

    DiffExp`Utilities`PrintDebug["Found homogeneous differential equation: ", HomogeneousEquation + O[DiffExp`Symbols`x]^4 // DiffExp`SeriesOps`SN][3];

    NMat = Table[
      If[ind < systemSize, UnitVector[systemSize, ind + 1], -HomogeneousEquation[[Range@systemSize]]]
      , {ind, systemSize}];

    DiffExp`Utilities`PrintDebug["Deriving solutions..."][3];
    Solns = DiffExp`Frobenius`FrobeniusSolutions[HomogeneousEquation];
    DiffExp`Utilities`PrintDebug["All homogeneous solutions found..."][3];
    DiffExp`Utilities`PrintDebug["Deriving period matrix..."][3];
    Wronsk = Table[
      DiffExp`SeriesOps`SD[Solns[[iind]], Sequence @@ DiffExp`Utilities`CA[DiffExp`Symbols`x, jind - 1]]
      , {jind, systemSize}, {iind, systemSize}
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
        If[ind < systemSize, UnitVector[systemSize, ind + 1], -HomogeneousEquation2[[Range@systemSize]]]
        , {ind, systemSize}];

      Solns2 = DiffExp`Frobenius`FrobeniusSolutions[HomogeneousEquation2];

      DiffExp`Utilities`PrintDebug["Deriving inverse Wronskian."][3];

      Wronsk2 = Table[
        DiffExp`SeriesOps`SD[Solns2[[iind]], Sequence @@ DiffExp`Utilities`CA[DiffExp`Symbols`x, jind - 1]]
        , {jind, systemSize}, {iind, systemSize}
      ] // DiffExp`SeriesOps`SExpand;

      (* Cross-checking Wronskians *)
      If[MemberQ[DiffExp`State`CurrCrosscheckFlags, "Wronskians"] === True,
        DiffExp`Utilities`PrintDebug["Cross-checking Wronskians."][1];
        crossCheck = (DiffExp`SeriesOps`SD[Wronsk, DiffExp`Symbols`x] - DiffExp`SeriesOps`MatrixMultiplySExpand[NMat, Wronsk] // (DiffExp`Utilities`CPChop@*DiffExp`SeriesOps`SExpand)) + O[DiffExp`Symbols`x]^DiffExp`State`ICrossCheckVerifyResultOrder;
        DiffExp`Utilities`PrintDebug["Found: ", crossCheck + O[DiffExp`Symbols`x]^DiffExp`State`ICrossCheckPrintResultOrder // DiffExp`SeriesOps`SN][3];
        If[
          !(SameQ @@ Append[crossCheck // Normal // Flatten, 0])
          ,
          DiffExp`Utilities`ReportError["Cross-check failed."];
        ];

        crossCheck = (DiffExp`SeriesOps`SD[Wronsk2, DiffExp`Symbols`x] - DiffExp`SeriesOps`MatrixMultiplySExpand[NMat2, Wronsk2] // (DiffExp`Utilities`PChop@*DiffExp`SeriesOps`SExpand)) + O[DiffExp`Symbols`x]^DiffExp`State`ICrossCheckVerifyResultOrder;
        DiffExp`Utilities`PrintDebug["Found: ", crossCheck + O[DiffExp`Symbols`x]^DiffExp`State`ICrossCheckPrintResultOrder // DiffExp`SeriesOps`SN][1];
        If[
          !(SameQ @@ Append[crossCheck // Normal // Flatten, 0])
          ,
          DiffExp`Utilities`ReportError["Cross-check failed."];
        ]
      ];

      WronskInvPrime = Inverse[MtildeMat2] . Wronsk2 // Transpose // DiffExp`SeriesOps`SExpand;
      DiffExp`Utilities`PrintDebug["Deriving inverse Wronskian.."][3];
      wronskianProduct = DiffExp`Utilities`PChop@Normal@DiffExp`SeriesOps`MatrixMultiplySExpand[WronskInvPrime, Wronsk];
      DiffExp`Utilities`PrintDebug["Deriving inverse Wronskian..."][3];
      If[MemberQ[DiffExp`State`CurrCrosscheckFlags, "WronskInv"] === True,
        If[DiffExp`Utilities`DependsQ[(wronskianProduct // DiffExp`Utilities`CPChop) + O[DiffExp`Symbols`x]^DiffExp`State`ICrossCheckVerifyResultOrder // Normal, DiffExp`Symbols`x],
          DiffExp`State`LastErrorContext = {WronskInvPrime, Wronsk, wronskianProduct};
          DiffExp`Utilities`ReportError["Warning, product of Wronskian inverse times Wronskian does not match identity. Try increasing \"ChopPrecision\" or \"WorkingPrecision\", or try decreasing \"ExpansionOrder\"."][1];
        ];
      ];
      DiffExp`State`LastErrorContext = {wronskianProduct};
      wronskianProduct = wronskianProduct /. DiffExp`Symbols`x^(k_ : 1) :> 0 /. DiffExp`Symbols`Logx -> 0;
      WronskInv = DiffExp`Utilities`PChop@DiffExp`SeriesOps`MatrixMultiplySExpand[(wronskianProduct // Inverse), WronskInvPrime];
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
      crossCheck = DiffExp`SeriesOps`SD[FMat, DiffExp`Symbols`x] - DiffExp`SeriesOps`MatrixMultiplySExpand[DiffExp`State`DEqnMatricesExpanded[line][0][[intind, intind]], FMat + O[DiffExp`Symbols`x]^DiffExp`State`ICrossCheckVerifyResultOrder] // (DiffExp`Utilities`CPChop@*DiffExp`SeriesOps`SExpand);
      DiffExp`State`LastErrorContext = {FMat, DiffExp`State`DEqnMatricesExpanded[line][0][[intind, intind]], crossCheck, MtildeMat, Wronsk};
      DiffExp`Utilities`PrintDebug["Found: ", crossCheck + O[DiffExp`Symbols`x]^DiffExp`State`ICrossCheckPrintResultOrder // DiffExp`SeriesOps`SN][3];
      If[
        !(SameQ @@ Append[crossCheck // Normal // Flatten, 0])
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
    crossCheck = DiffExp`SeriesOps`SD[GMat, DiffExp`Symbols`x] - DiffExp`SeriesOps`MatrixMultiplySExpand[DiffExp`State`DEqnMatricesExpanded[line][0][[intind, intind]], GMat + O[DiffExp`Symbols`x]^DiffExp`State`ICrossCheckVerifyResultOrder] - BMat // (DiffExp`Utilities`CPChop@*DiffExp`SeriesOps`SExpand);
    DiffExp`State`LastErrorContext = {FMat, GMat, FMatInv, BMat, DiffExp`State`DEqnMatricesExpanded[line][0][[intind, intind]], crossCheck, MtildeMat, Wronsk};
    DiffExp`Utilities`PrintDebug["Found: ", crossCheck + O[DiffExp`Symbols`x]^DiffExp`State`ICrossCheckPrintResultOrder // DiffExp`SeriesOps`SN][3];
    If[
      !(SameQ @@ Append[crossCheck + O[DiffExp`Symbols`x]^DiffExp`State`ICrossCheckVerifyResultOrder // Normal // Flatten, 0])
      ,
      DiffExp`Utilities`ReportError["Cross-check of solution matrix failed"];
    ];
  ];

  fGeneral = Total[GMat // Transpose];

  {cIndices, fGeneral, BufferedData}
];
