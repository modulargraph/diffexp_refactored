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
  "DiffExp`MatrixLoading`",
  "DiffExp`IntegrationStrategies`"
}];

(* Transport functions *)
PrepareBoundaryConditions::usage = "PrepareBoundaryConditions[bcs_List, line_List] prepares boundary conditions for use with IntegrateSystem or TransportTo.";
IntegrateSystem::usage = "IntegrateSystem[line_List] obtains general series solutions along a line. IntegrateSystem[bcs_List, line_List] uses boundary conditions.";
TransportTo::usage = "TransportTo[bcs_List, line_List, to_:1, save_:False] transports boundary conditions to arbitrary points.";
ToPiecewise::usage = "ToPiecewise[segmentdata_List, pade_:False] converts segment data to piecewise functions.";

Begin["`Private`"];

(* Prepare boundary conditions *)
PrepareBoundaryConditions[bcs_List, line2_Association | line2_List] := Module[
  {line, CoeffList = {}, Coeffs, CoeffSer, bcs1, unitSeries = (1 // N[#, DiffExp`State`FEWorkingPrecision] &), Mask, ispoint, LineRat, maskedEntry},

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
          DiffExp`SeriesOps`LeadingCoefficientSeries[(#[[ind]] /. line) * unitSeries, 2]
        ) /. unitSeries "?" + O[DiffExp`Symbols`x]^(1/2) -> "?",
        #[[ind]]
      ]
      , {ind, Length[#]}
    ]
    ,
    If[DiffExp`Utilities`DependsQ[Normal[#], DiffExp`Symbols`x],
      DiffExp`Utilities`ReportError["The boundary terms that are provided depend on the line parameter ", DiffExp`Symbols`x, ", but the line itself does not."];
    ];
    unitSeries * # /. line
  ] & /@ CoeffList;

  If[!ispoint,
    Mask = MapAt[Switch[#[[0]],
      SeriesData,
      maskedEntry = #;
      maskedEntry[[3]] = Table["(...)", {ind, #[[3]] // Length}];
      maskedEntry,
      _, #] &, bcs1, {All, All}] // TableForm;

    DiffExp`Utilities`PrintInfo["Prepared boundary conditions in asymptotic limit, of the form:"][1];
    DiffExp`Utilities`PrintInfo[Mask][1];
  ];

  {line, bcs1} /. Log[a_ DiffExp`Symbols`x] /; NumericQ[a] :> Log[a] + DiffExp`Symbols`Logx /. Log[DiffExp`Symbols`x] -> DiffExp`Symbols`Logx
];

(* Build a SegmentContext association for a given integral block and line *)
BuildSegmentContext[intind_, line_] := <|
  "AMatExpanded" -> DiffExp`State`DEqnMatricesExpanded[line][0][[intind, intind]],
  "AMatFactored" -> If[KeyExistsQ[DiffExp`State`DEqnMatricesFactored, line],
    DiffExp`State`DEqnMatricesFactored[line][0][[intind, intind]],
    Missing["NotAvailable"]
  ],
  "SystemSize" -> Length[intind],
  "ExpansionOrder" -> DiffExp`State`ExpansionOrderVal,
  "WorkingPrecision" -> DiffExp`State`FEWorkingPrecision,
  "ChopPrecision" -> DiffExp`State`ChopPrecisionVal,
  "LinearSolveChopPrecision" -> DiffExp`State`LinearSolveChopPrecisionVal,
  "HomogeneousSolve" -> DiffExp`State`FEC["HomogeneousSolve"],
  "InvWronskSolver" -> DiffExp`State`FEC["InvWronskSolver"],
  "CrosscheckFlags" -> DiffExp`State`CurrCrosscheckFlags,
  "CrossCheckPrintOrder" -> DiffExp`State`ICrossCheckPrintResultOrder,
  "CrossCheckVerifyOrder" -> DiffExp`State`ICrossCheckVerifyResultOrder,
  "IntegrationStrategy" -> DiffExp`State`FEC[IntegrationStrategy],
  "UseRationalRecurrence" -> DiffExp`State`FEC[UseRationalRecurrence],
  "Label" -> intind
|>;

(* Main integration function *)
IntegrateSystem[bcs2 : (_List | _Association) : "?", line2_Association | line2_List, opts2_ : {}] := Module[
  {bcs, line, BCSRelevant, relevantIndices, IgnorePositions, crossCheck, bVec, IntegrationData, fGeneral,
   FixAt, BoundaryEqns1, BoundaryEqns2, cIndices, Cmat, Cb, csol, NewResults, opts = opts2,
   DEqnMatricesExpandedCopy, TurnOffPade, constantReplacements, csFreedom, solveFailed, LogsPresent,
   AlgebraicRootsPresent, normalizedSolutions, segmentCaches, ctx, blockCache, benchStart,
   IntegrationDataTab, bVec0, bVec1, bVecRest, complementIndices, otherIntegralData},

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

  segmentCaches = Association[{}];

  (* Handle Association output from TransportTo/IntegrateSystem *)
  If[AssociationQ[bcs2],
    bcs = {bcs2["KinematicPoint"], bcs2["SeriesValues"]};,
    If[MatchQ[bcs2, {{a_Association, __}, _}],
      bcs = bcs2[[1]];,
      bcs = bcs2;
    ];
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

      complementIndices = Complement[Range[DiffExp`State`NumIntegrals], intind];
      otherIntegralData = Table[IntegrationData[{ind, myeps}], {ind, complementIndices}, {myeps, 0, epsord}];
      bVec0 = DiffExp`SeriesOps`SExpand[DEqnMatricesExpandedCopy[0][[intind, complementIndices]] . otherIntegralData[[All, epsord + 1]]];
      Which[epsord === 0,
        bVec = bVec0 // DiffExp`SeriesOps`SExpand;,
        epsord === 1,
        bVec = bVec0 + bVec1[[intind]] // DiffExp`SeriesOps`SExpand;,
        epsord > 1,
        bVec = bVec0 + bVec1[[intind]] + bVecRest[[intind]] // DiffExp`SeriesOps`SExpand;
      ];

      DiffExp`Utilities`PrintDebug["Done."][3];

      (* Dispatch to appropriate integration strategy *)
      ctx = BuildSegmentContext[intind, line];
      blockCache = Lookup[segmentCaches, Key[intind], <||>];
      If[epsord === 0,
        benchStart = AbsoluteTime[];
      ];
      {cIndices, fGeneral, blockCache} = DiffExp`IntegrationStrategies`DispatchStrategy[ctx, bVec, epsord, blockCache];
      segmentCaches[intind] = blockCache;
      If[epsord === 0,
        DiffExp`State`BenchmarkData["Segments"][line // N]["HomogeneousSolveAllPreprocessing"]["Integrals"][intind] = AbsoluteTime[] - benchStart;
      ];

      (* Code for fixing boundary conditions. *)
      If[!(bcs2 === "?"),
        DiffExp`Utilities`PrintDebug["General solution found. Fixing the indeterminates from boundary conditions."][2];

        BCSRelevant = bcs[[2]][[intind]][[All, epsord + 1]];
        IgnorePositions = Flatten[Position[BCSRelevant, "?"]];
        relevantIndices = Complement[Range@Length[intind], IgnorePositions];
        BCSRelevant = BCSRelevant[[relevantIndices]];

        If[FixAt === 0,
          BoundaryEqns1 = (DiffExp`Utilities`PChop@*DiffExp`SeriesOps`SExpand) /@ (fGeneral[[relevantIndices]] - BCSRelevant);

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
          ] + DiffExp`Pade`SEval[term /. (# -> 0 & /@ cIndices), FixAt], {term, fGeneral[[relevantIndices]]}] - BCSRelevant;

          BoundaryEqns2 = # == 0 & /@ Flatten[BoundaryEqns1];
        ];

        BoundaryEqns2 = DeleteCases[BoundaryEqns2, True] // Expand;

        If[!BoundaryEqns2 === {},
          DiffExp`Utilities`PrintDebug["Boundary equations are: "][3];
          DiffExp`Utilities`PrintDebug[BoundaryEqns2 // N // ReplaceAll[c -> "c"]][3];

          {Cmat, Cb} = {#[[2]], -#[[1]]} &@CoefficientArrays[BoundaryEqns2, cIndices];

          solveFailed = False;
          If[!MemberQ[BoundaryEqns2, False],
            Check[csol = LinearSolve[Cmat, Cb, ZeroTest -> (N[DiffExp`Utilities`LSPChop@Expand@Normal[#1], DiffExp`State`LinearSolveChopPrecisionVal] == 0 &)];
              , solveFailed = True;];
            , solveFailed = True;];

          If[solveFailed,
            DiffExp`State`LastErrorContext = {fGeneral[[relevantIndices]], BoundaryEqns1, BoundaryEqns2, BCSRelevant, Cmat, Cb, bVec, intind, epsord};
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

          constantReplacements = Thread[cIndices -> csol] // DiffExp`Utilities`PChop;

          fGeneral = fGeneral /. constantReplacements // DiffExp`SeriesOps`SExpand;

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
        crossCheck = DiffExp`SeriesOps`SD[fGeneral, DiffExp`Symbols`x] - DiffExp`State`DEqnMatricesExpanded[line][0][[intind, intind]] . fGeneral - bVec + O[DiffExp`Symbols`x]^DiffExp`State`ICrossCheckVerifyResultOrder // (DiffExp`Utilities`CPChop@*DiffExp`SeriesOps`SExpand);
        DiffExp`Utilities`PrintDebug["Found: ", crossCheck + O[DiffExp`Symbols`x]^DiffExp`State`ICrossCheckPrintResultOrder // DiffExp`SeriesOps`SN][3];
        If[
          !(SameQ @@ Append[crossCheck // Normal // Flatten, 0])
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
    normalizedSolutions = Values[IntegrationData] // Normal // Chop[#, 10^-DiffExp`State`ICheckMultivaluedChop] &;
    LogsPresent = DiffExp`Utilities`DependsQ[normalizedSolutions, DiffExp`Symbols`Logx | DiffExp`Symbols`Logx^_];
    AlgebraicRootsPresent = DiffExp`Utilities`DependsQ[normalizedSolutions, (DiffExp`Symbols`x^b_ /; Denominator[b] > 1)];

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
    <|
      "SeriesValues" -> Table[DiffExp`SeriesOps`ApplyAnalyticContinuation[IntegrationData[{ind, epsord}]] // DiffExp`AnalyticContinuation`Project\[Theta]s, {ind, DiffExp`State`NumIntegrals}, {epsord, 0, DiffExp`State`EpsilonOrderVal}],
      "NumIntegrals" -> DiffExp`State`NumIntegrals,
      "EpsilonOrder" -> DiffExp`State`EpsilonOrderVal,
      "ExpansionOrder" -> DiffExp`State`ExpansionOrderVal
    |>
  ]
];

(* TransportTo function - transports boundary conditions *)
TransportTo[bcs2 : (_List | _Association), line2_Association | line2_List, to2 : _?NumericQ : 1, SaveExpansions : _?BooleanQ : False, SampleAtList_List : {}] := Module[
  {line = If[line2[[0]] === List, line2 // Association // KeySort, line2 // KeySort], to, tempResult, bcs, FixAt, singularities, imaginarySingularities, relevantSingularities, SingularitySegments, PoleIntervals, CurrLine, CurrLineNoMobius = Null, CurrIntegrated, CurrIntervalCurrLine, CurrIntervalLine, Done = False, overlapCheck, CurrEvalPoint, CurrEvalPointCurrLine, CurrEval, Currbcs, AllIntegrationData = {}, currentCenter, EvaluateCurrPoint, NextIsPole = False, SegmentCounter = 1, lineRelation, CurrIntervalLinePos, CurrIntervalLineNeg, CurrEvalError, CurrIntegratedError, CurrbcsError, PrintError, CurrError, accumulatedError = 0, accumulatedErrors = ConstantArray[0, {DiffExp`State`NumIntegrals, DiffExp`State`EpsilonOrderVal + 1}], FixWithin, SegmentsToIntegrate, UpdateMatrixExpansionError, TimeStart, TimeStart0, LineReturn, FailedLine, ExpansionsIndeterminates = {}, previousBoundaryConditions, CurrStatusBackup, BoundaryFixPoint, CurrEvalErrorEx, CurrEvalError1, CurrEvalError2, CurrEvalAtBoundaryFixPoint, CurrEvalEx, AllSegmentsPredivision, tempFile, CompressedTermForExport, CompressedTermForExportFN, CurrLineLR, FullLineLR, cachedPoleIntervals, RepeatingSegment, LastEvaluation, LastSavedData, LastLine, ExpansionOrders, DigitsNeeded},

  DiffExp`State`BenchmarkData = Association[];
  DiffExp`State`BenchmarkData["TimeStart"] = AbsoluteTime[];

  DiffExp`State`MultivaluedFail = False;

  If[Or @@ (!DiffExp`Utilities`DependsQ[Keys@line, #] & /@ DiffExp`State`ExternalScalesVal),
    DiffExp`Utilities`ReportError["The point/line in the second argument does not fix all kinematic invariants and masses!"];
  ];

  DiffExp`State`ExpansionOrderVal = DiffExp`State`FEC[ExpansionOrder];
  DigitsNeeded = DiffExp`State`FEAccuracyGoal + DiffExp`State`ISafetyDigits;
  TimeStart0 = AbsoluteTime[];

  (* Handle Association output from TransportTo/IntegrateSystem *)
  If[AssociationQ[bcs2],
    bcs = {bcs2["KinematicPoint"], bcs2["SeriesValues"]};,
    If[MatchQ[bcs2, {{a_Association, __}, _}],
      bcs = bcs2[[1]];,
      bcs = bcs2;
    ];
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
      tempResult = DiffExp`LineSegmentation`RelateLines[bcs[[1]], line, True];
      If[tempResult === False, DiffExp`Utilities`ReportError["Endpoint does not lie on same line as the boundary conditions."]];
      line = bcs[[1]] /. DiffExp`Symbols`x -> DiffExp`Symbols`x tempResult;
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
  If[Length[tempResult = Flatten[Position[Factor[DiffExp`State`MatricesIrreducibleFactors /. line], 0]]] > 0,
    DiffExp`Utilities`ReportError["The line lies on a singularity of the differential equations. The vanishing factors are: ", DiffExp`State`MatricesIrreducibleFactors[[tempResult]], "."];
  ];

  (* Will abort if something is wrong with the boundary conditions. *)
  {bcs, FixAt} = DiffExp`LineSegmentation`CheckBoundaryConditionsAndReparametrize[previousBoundaryConditions = bcs, line];

  If[Length[previousBoundaryConditions] > 2,
    accumulatedError = previousBoundaryConditions[[3]] // Abs // Max;
    accumulatedErrors = previousBoundaryConditions[[3]];
    ,
    Quiet[
      accumulatedErrors = 10^-(MapAt[Accuracy, previousBoundaryConditions[[2]], {All, All}] // SetPrecision[#, DiffExp`State`FEWorkingPrecision] &);
      accumulatedError = accumulatedErrors // Max;
    ];
  ];

  If[!((accumulatedErrors // Dimensions // Last) === DiffExp`State`EpsilonOrderVal + 1),
    accumulatedErrors = PadRight[#, DiffExp`State`EpsilonOrderVal + 1, 10^-DiffExp`State`FEWorkingPrecision] & /@ accumulatedErrors;
  ];

  DiffExp`Utilities`PrintInfo["Transporting boundary conditions along ", line // Normal // N // Association, " from x = ", FixAt // N, " to x = ", to // N][1];

  DiffExp`State`BenchmarkData["BasicPreprocessing"] = AbsoluteTime[] - DiffExp`State`BenchmarkData["TimeStart"];
  DiffExp`State`BenchmarkData["MatrixPreprocessing"] = AbsoluteTime[];

  DiffExp`Utilities`PrintInfo["Preparing differential equations along current line.."][1];
  DiffExp`MatrixLoading`PrepareMatricesFactored[line];

  DiffExp`MatrixLoading`InitializeIntegrationSequence[line];

  {singularities, imaginarySingularities} = DiffExp`LineSegmentation`FindMatrixSingularities[line, True, {FixAt, to}];

  relevantSingularities = Select[singularities, FixAt <= # <= to || to <= # <= FixAt &];
  If[DiffExp`State`FEC[UseMobius] === True,
    SingularitySegments = {
      #,
      DiffExp`Mobius`GetLineRescaled[line, #, {singularities, imaginarySingularities}],
      DiffExp`Mobius`GetLineRescaled[line, #, {singularities, imaginarySingularities}, True]
    } & /@ relevantSingularities;
    ,
    SingularitySegments = {#, DiffExp`Mobius`GetLineRescaled[line, #, {singularities, imaginarySingularities}]} & /@ relevantSingularities;
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

      tempResult = {-#, #} &@DiffExp`LineSegmentation`GetMatricesPrecisionDistance[sline[[2]]];
    ];

    If[DiffExp`State`FEC[SegmentationStrategy] === "Predivision",
      tempResult = {-DiffExp`State`RadiusOfConvergenceVal/DiffExp`State`FEC[DivisionOrder], DiffExp`State`RadiusOfConvergenceVal/DiffExp`State`FEC[DivisionOrder]};
    ];

    tempResult = {sline, DiffExp`LineSegmentation`RelateLinesPoint[line, sline[[2]], #] & /@ tempResult};

    DiffExp`Utilities`PrintInfo["Expansion around x = ", sline[[1]] // N, " is valid within region x \[Element] [", tempResult[[2, 1]] // N, ", ", tempResult[[2, 2]] // N, "]."][2];

    tempResult
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
    cachedPoleIntervals = PoleIntervals;
    currentCenter = FixAt;
    CurrLine = DiffExp`Mobius`GetLineRescaled[line, FixAt, {singularities, imaginarySingularities}];
    AppendTo[AllSegmentsPredivision, CurrLine];

    Done = False;
    While[!Done,
      PoleIntervals = Select[PoleIntervals, !(#[[1, 1]] === currentCenter) &];

      CurrIntervalCurrLine = {-DiffExp`State`RadiusOfConvergenceVal/DiffExp`State`FEC[DivisionOrder], DiffExp`State`RadiusOfConvergenceVal/DiffExp`State`FEC[DivisionOrder]};
      CurrIntervalLine = DiffExp`LineSegmentation`RelateLinesPoint[line, CurrLine, #] & /@ CurrIntervalCurrLine;

      If[IntervalContainsQ[CurrIntervalLine, to],
        CurrEvalPoint = to;
        CurrEvalPointCurrLine = DiffExp`LineSegmentation`RelateLinesPoint[CurrLine, line, CurrEvalPoint];

        Done = True;

        ,
        SegmentCounter += 1;

        If[to > FixAt,
          tempResult = Select[PoleIntervals, #[[1, 1]] > currentCenter &];
          ,
          tempResult = Select[PoleIntervals, #[[1, 1]] <= currentCenter &];
        ];
        overlapCheck = {#[[1]], IntervalOverlapQ[CurrIntervalLine, #[[2]]]} & /@ tempResult;
        overlapCheck = Flatten[Position[overlapCheck[[All, 2]], True]];
        NextIsPole = Length[overlapCheck] > 0;

        If[NextIsPole,
          tempResult = tempResult[[overlapCheck[[1]]]];

          If[to > FixAt,
            FixWithin = IntervalIntersec[IntervalIntersec[{currentCenter, tempResult[[1, 1]]}, tempResult[[2]]], CurrIntervalLine];
            CurrEvalPoint = (FixWithin[[1]] + FixWithin[[2]])/2;
            ,

            FixWithin = IntervalIntersec[IntervalIntersec[{tempResult[[1, 1]], currentCenter}, tempResult[[2]]], CurrIntervalLine];
            CurrEvalPoint = (FixWithin[[1]] + FixWithin[[2]])/2;
          ];

          CurrEvalPointCurrLine = DiffExp`LineSegmentation`RelateLinesPoint[CurrLine, line, CurrEvalPoint];

          currentCenter = tempResult[[1, 1]];
          CurrLine = tempResult[[1, 2]];
          If[DiffExp`State`FEC[UseMobius] === True, CurrLineNoMobius = tempResult[[1, 3]]];
          ,

          (* Finite point *)
          If[to > FixAt,
            CurrEvalPoint = CurrIntervalLine[[2]],
            CurrEvalPoint = CurrIntervalLine[[1]]
          ];
          CurrEvalPointCurrLine = DiffExp`LineSegmentation`RelateLinesPoint[CurrLine, line, CurrEvalPoint];

          If[DiffExp`State`FEC[SegmentationStrategy] === "Dynamic",
            currentCenter = CurrEvalPoint;
          ];

          If[DiffExp`State`FEC[SegmentationStrategy] === "Predivision",
            If[to > FixAt,
              currentCenter = DiffExp`Mobius`FindNextCenterPointL[CurrEvalPoint, singularities],
              currentCenter = DiffExp`Mobius`FindNextCenterPointR[CurrEvalPoint, singularities]
            ];
          ];

          CurrLine = DiffExp`Mobius`GetLineRescaled[line, currentCenter, {singularities, imaginarySingularities}];
        ];
      ];

      AppendTo[AllSegmentsPredivision, CurrLine];

      PoleIntervals = cachedPoleIntervals;
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
      lineRelation = DiffExp`Symbols`x -> DiffExp`LineSegmentation`RelateLines[CurrLine, line];
      AppendTo[AllIntegrationData, {
        CurrLine, (* Current line *)
        lineRelation, (* Change of line parameter from CurrLine to line *)
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
      DiffExp`State`LastErrorContext = CurrIntegrated;
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

  ComputeErrorsPerIndeterminate[aaa_, bbb_, indeterminates_] := Module[{errorsByOrder},
    errorsByOrder = Table[DiffExp`SeriesOps`LogxCoeffNS[aaa - bbb, logxord], {logxord, 0, DiffExp`State`IMaxLogOrder}];
    errorsByOrder = Table[(

      Flatten[
        Append[
          Table[
            Coefficient[errorsByOrder[[All, ii, jj]], var]
            , {var, indeterminates}]
          ,
          errorsByOrder[[All, ii, jj]] /. (# -> 0 & /@ indeterminates)
        ]
      ]

    ) // Abs // Max, {ii, (errorsByOrder // Dimensions)[[2]]}, {jj, (errorsByOrder // Dimensions)[[3]]}] // N
  ];

  PrintError[] := Block[{segmentErrors, endpointErrors, boundaryErrors, indeterminates},

    If[!DiffExp`State`MultivaluedFail,

      indeterminates = DeleteCases[Currbcs[[2]] // System`Variables, "?" | DiffExp`Symbols`x | DiffExp`Symbols`Logx];

      endpointErrors = ComputeErrorsPerIndeterminate[CurrEvalError1, CurrEval, indeterminates];
      boundaryErrors = ComputeErrorsPerIndeterminate[CurrEvalError2, CurrEvalAtBoundaryFixPoint, indeterminates];
      segmentErrors = Table[Max[{endpointErrors[[ii, jj]], boundaryErrors[[ii, jj]]}], {ii, First@Dimensions@endpointErrors}, {jj, Last@Dimensions@endpointErrors}];

      CurrError = Flatten[segmentErrors] // Abs // Max;

      Switch[DiffExp`State`FEC["EstimateError"],
        "Fast",
        accumulatedError = accumulatedError + CurrError;
        accumulatedErrors = accumulatedErrors + segmentErrors;
      ];

      Which[DiffExp`State`FEC["EstimateError"] == "Fast",
        DiffExp`Utilities`PrintInfo["Current segment error estimate: ", CurrError][1];
        DiffExp`Utilities`PrintInfo["Total error estimate: ", accumulatedError][1];
      ];

      If[CurrError > 1,
        DiffExp`Utilities`ReportError["The reported error is very large. This likely indicates a numerical instability."];
      ];

    ];

  ];

  currentCenter = FixAt;
  CurrLine = DiffExp`Mobius`GetLineRescaled[line, FixAt, {singularities, imaginarySingularities}];
  If[DiffExp`State`FEC[UseMobius] === True, CurrLineNoMobius = DiffExp`Mobius`GetLineRescaled[line, FixAt, {singularities, imaginarySingularities}, True]];
  Currbcs = CurrbcsError = bcs;

  DiffExp`State`BenchmarkData["Segments"] = Association[];

  Done = False;
  While[!Done,
    DiffExp`State`BenchmarkData["Segments"][CurrLine // N] = Association[];
    DiffExp`State`BenchmarkData["Segments"][CurrLine // N]["ComputationTime"] = TimeStart = AbsoluteTime[];

    CurrStatusBackup = {Currbcs, CurrbcsError, bcs, CurrLine, currentCenter, FixAt, accumulatedError, accumulatedErrors};

    DiffExp`Utilities`PrintInfo["Integrating segment: ", DiffExp`Utilities`PrintMobiusNormalized /@ CurrLine, "."][1];
    If[NextIsPole && MemberQ[imaginarySingularities, currentCenter],
      DiffExp`Utilities`PrintInfo["Current segment is centered at singularity."][1];
    ];

    LastEvaluation = Currbcs;
    LastSavedData = AllIntegrationData;
    LastLine = CurrLine;

    PoleIntervals = Select[PoleIntervals, !(#[[1, 1]] === currentCenter) &];

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
        tempResult = Select[PoleIntervals, #[[1, 1]] > currentCenter &];
        ,
        tempResult = Select[PoleIntervals, #[[1, 1]] <= currentCenter &];
      ];
      overlapCheck = {#[[1]], IntervalOverlapQ[CurrIntervalLine, #[[2]]]} & /@ tempResult;
      overlapCheck = Flatten[Position[overlapCheck[[All, 2]], True]];
      NextIsPole = Length[overlapCheck] > 0;

      If[NextIsPole,
        tempResult = tempResult[[overlapCheck[[1]]]];

        If[to > FixAt,
          (* (Current line center, pole center), (poleinterval) *)
          FixWithin = IntervalIntersec[IntervalIntersec[{currentCenter, tempResult[[1, 1]]}, tempResult[[2]]], CurrIntervalLine];
          CurrEvalPoint = (FixWithin[[1]] + FixWithin[[2]])/2;

          ,

          (* (pole center, current line center), (poleinterval) *)
          FixWithin = IntervalIntersec[IntervalIntersec[{tempResult[[1, 1]], currentCenter}, tempResult[[2]]], CurrIntervalLine];
          CurrEvalPoint = (FixWithin[[1]] + FixWithin[[2]])/2;


        ];


        CurrEvalPointCurrLine = DiffExp`LineSegmentation`RelateLinesPoint[CurrLine, line, CurrEvalPoint];
        EvaluateCurrPoint[];

        currentCenter = tempResult[[1, 1]];
        CurrLine = tempResult[[1, 2]];
        If[DiffExp`State`FEC[UseMobius] === True, CurrLineNoMobius = tempResult[[1, 3]]];
        ,
        (* Finite point *)
        If[to > FixAt,
          CurrEvalPoint = CurrIntervalLine[[2]],
          CurrEvalPoint = CurrIntervalLine[[1]]
        ];
        CurrEvalPointCurrLine = DiffExp`LineSegmentation`RelateLinesPoint[CurrLine, line, CurrEvalPoint];
        EvaluateCurrPoint[];

        If[DiffExp`State`FEC[SegmentationStrategy] === "Dynamic",
          currentCenter = CurrEvalPoint;
        ];
        If[DiffExp`State`FEC[SegmentationStrategy] === "Predivision",
          If[to > FixAt,
            currentCenter = DiffExp`Mobius`FindNextCenterPointL[CurrEvalPoint, singularities],
            currentCenter = DiffExp`Mobius`FindNextCenterPointR[CurrEvalPoint, singularities]
          ];
        ];

        CurrLine = DiffExp`Mobius`GetLineRescaled[line, currentCenter, {singularities, imaginarySingularities}];
        If[DiffExp`State`FEC[UseMobius] === True, CurrLineNoMobius = DiffExp`Mobius`GetLineRescaled[line, currentCenter, {singularities, imaginarySingularities}, True]];
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
        {Currbcs, CurrbcsError, bcs, CurrLine, currentCenter, FixAt, accumulatedError, accumulatedErrors} = CurrStatusBackup;
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

  <|
    "KinematicPoint" -> LineReturn,
    "SeriesValues" -> CurrEval,
    "ErrorEstimates" -> If[!DiffExp`State`FEC["EstimateError"] === False,
      accumulatedErrors, Missing["NotComputed"]],
    "SegmentData" -> If[SaveExpansions === True,
      AllIntegrationData, Missing["NotRequested"]],
    "ComputationTime" -> DiffExp`State`BenchmarkData["ComputationTime"],
    "NumIntegrals" -> DiffExp`State`NumIntegrals,
    "EpsilonOrder" -> DiffExp`State`EpsilonOrderVal,
    "ExpansionOrder" -> DiffExp`State`FEC[ExpansionOrder]
  |>
];

(* Helper functions *)
IntervalOverlapQ[intv1_, intv2_] := !(IntervalIntersection[Interval[intv1], Interval[intv2]] === Interval[]);
IntervalIntersec[intv1_, intv2_] := IntervalIntersection[Interval[intv1], Interval[intv2]][[1]];
IntervalContainsQ[intv_, point_] := intv[[1]] <= point <= intv[[2]];

(* ToPiecewise - convert saved segment data to piecewise functions *)
ToPiecewise[SavedData2_, Pade : _?BooleanQ : False, Ord_Integer : Null] := Module[
  {SavedData, piecewiseResult, Uncompressed, Counter},

  If[AssociationQ[SavedData2],
    If[MissingQ[SavedData2["SegmentData"]],
      DiffExp`Utilities`ReportError["No segment data. TransportTo was not called with SaveExpansions -> True."];
    ];
    SavedData = SavedData2["SegmentData"],
    If[MatchQ[SavedData2, {{a_Association, _}, {__}}] || MatchQ[SavedData2, {{a_Association, _, _}, {__}}],
      SavedData = SavedData2[[2]],
      SavedData = SavedData2
    ]
  ];

  If[!(MatchQ[SavedData[[0]] === List] && Quiet[Dimensions[SavedData][[2]] === 5]),
    DiffExp`Utilities`ReportError["Could not interpret the argument. Maybe TransportTo[...] was not called with the option save_ set to True?"];
  ];

  Counter = 1;
  If[!$FrontEnd === Null,
    PrintTemporary["Processing ", Dynamic[Counter]];
  ];

  If[DiffExp`State`FEC["SaveExpansionsCompress"] === True,
    If[!DiffExp`State`FEC["SaveExpansionsCompressDirectory"] === "?",
      If[StringJoin[SavedData[[1, 5]] // StringPart[#, -2 ;; -1] &] === ".m",
        Uncompressed[ind_] := Uncompressed[ind] = Uncompress[Import[SavedData[[ind, 5]]]];,
        Uncompressed[ind_] := Uncompressed[ind] = Uncompress[SavedData[[ind, 5]]];
      ];,
      Uncompressed[ind_] := Uncompressed[ind] = Uncompress[SavedData[[ind, 5]]];
    ],
    Uncompressed[ind_] := SavedData[[ind, 5]];
  ];

  Table[
    piecewiseResult = Piecewise@Table[
      Counter = {ind, intind, epsord};
      {
        (If[Pade === True,
          (DiffExp`AnalyticContinuation`Project\[Theta]s[#, DiffExp`Pade`GetPade] &@#) /.
            DiffExp`Symbols`Logx -> Log[DiffExp`Symbols`x] /.
            DiffExp`Symbols`\[Theta]p -> HeavisideTheta[DiffExp`Symbols`x] /.
            DiffExp`Symbols`\[Theta]m -> HeavisideTheta[-DiffExp`Symbols`x] /.
            (SavedData[[ind, 2]]),
          (Normal@#) /.
            DiffExp`Symbols`Logx -> Log[DiffExp`Symbols`x] /.
            DiffExp`Symbols`\[Theta]p -> HeavisideTheta[DiffExp`Symbols`x] /.
            DiffExp`Symbols`\[Theta]m -> HeavisideTheta[-DiffExp`Symbols`x] /.
            (SavedData[[ind, 2]])
        ] &@(Uncompressed[ind][[intind, epsord]] + If[Ord === Null, 0, O[DiffExp`Symbols`x]^Ord])),
        DiffExp`Symbols`x >= SavedData[[ind, 3, 1]] && DiffExp`Symbols`x <= SavedData[[ind, 3, 2]]
      },
      {ind, Length@SavedData}
    ];
    Evaluate[piecewiseResult /. DiffExp`Symbols`x -> #] &,
    {intind, Uncompressed[1] // Dimensions // First},
    {epsord, Uncompressed[1] // Dimensions // Last}
  ]
];

End[];

EndPackage[];
