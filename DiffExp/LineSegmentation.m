(* ::Package:: *)

(* DiffExp LineSegmentation Subpackage *)
(* This package provides line segmentation and interval functions *)

BeginPackage["DiffExp`LineSegmentation`", {"DiffExp`Symbols`", "DiffExp`State`", "DiffExp`Utilities`", "DiffExp`SeriesOps`"}];

(* Line segmentation functions *)
RelateLines::usage = "RelateLines[a_,b_,noerror_] relates line parameters of different segments.";
RelateLinesPoint::usage = "RelateLinesPoint[a_,b_,pointb_] relates lines at a specific point.";
FindMatrixSingularities::usage = "FindMatrixSingularities[line_,getcomplex_,{fixat_,to_}] finds singularities.";
GetLargestTerm::usage = "GetLargestTerm[line_] gives a large term for error estimation.";
GetMatricesPrecisionDistance::usage = "GetMatricesPrecisionDistance[line_] gets intervals in dynamic segmentation.";
CheckBoundaryConditionsAndReparametrize::usage = "CheckBoundaryConditionsAndReparametrize[bcs_,line_] validates and reparametrizes boundary conditions.";
GetMatchingPoint::usage = "GetMatchingPoint[line_,bcsline_] gets the matching point on a line.";

Begin["`Private`"];

(* Relates line parameters of different segments *)
RelateLines[a2_Association, b2_Association, noerror_: False] := Module[
  {lineA = KeySort[a2], lineB = KeySort[b2], lineSolution, allSolutions, CanSolveFrom, solvableKey, suggestedSolution},

  (* Find one component of a that contains x, and use that to solve the equation *)
  CanSolveFrom = Flatten[Position[DiffExp`Utilities`DependsQ[DiffExp`Utilities`PChop@Expand@#, DiffExp`Symbols`x] & /@ lineA, True]];
  If[Length[CanSolveFrom] === 0,
    DiffExp`State`LastErrorContext = {lineA, lineB};
    DiffExp`Utilities`ReportError["Argument does not depend on the line parameter!"];
  ];

  solvableKey = CanSolveFrom[[1, 1]];

  lineSolution = Quiet[Solve[(lineA[solvableKey] /. DiffExp`Symbols`x -> y) - lineB[solvableKey] == 0, y]];

  If[Length[lineSolution] > 1,
    DiffExp`Utilities`ReportError["Multiple solutions encountered while relating lines: ",
      lineA // Normal // N // Association, " and ", lineB // Normal // N // Association,
      ". This is likely because you chose a non-linear line segment."];
  ];

  If[!(Length[lineSolution] === 1),
    DiffExp`State`LastErrorContext = {lineA, lineB};
    If[!noerror, DiffExp`Utilities`ReportError["Could not relate lines."];];
    False,
    (* Crosscheck *)
    suggestedSolution = lineSolution[[1, 1, 2]] // Together // DiffExp`Utilities`PChop;

    Quiet[
      If[!(SameQ @@ Append[Values[(DiffExp`Utilities`PChop@Together@N[#, DiffExp`State`FEC[WorkingPrecision]]) & /@
          ((lineA /. DiffExp`Symbols`x -> suggestedSolution) - lineB)], 0]),
        DiffExp`Utilities`PrintDebug["Could not relate lines:"][1];
        DiffExp`Utilities`PrintDebug[Values[(lineA /. DiffExp`Symbols`x -> suggestedSolution) - lineB] // Factor // N][1];
        DiffExp`State`LastErrorContext = {lineA, lineB};
        If[!noerror, DiffExp`Utilities`ReportError["Could not relate lines."];];
        False,
        suggestedSolution
      ]
    ]
  ]
];

RelateLinesPoint[a_Association, b_Association, pointb_] :=
  RelateLines[a, b] /. DiffExp`Symbols`x -> pointb;

(* Derives singularity positions *)
FindMatrixSingularities[line_, getcomplex_: False, {fixat_, to_}] := Module[
  {AllSingularTerms, RatPoles, Sqrts, singularityPositions, singularityPositionsNumeric, projectedSingularities, SPosProjections, currentSingularity, currentSingularityNumeric, NumPoles, Tmp, Tmp2, Tmp3, Tmp4, leftInterval, rightInterval},

  DiffExp`Utilities`PrintInfo["Determining positions of singularities and branch-cuts."][1];

  AllSingularTerms = Flatten[FactorList /@ Factor[Union[
    DiffExp`State`MatricesIrreducibleFactors /. line,
    Denominator@Together@Values@line
  ]]] // DeleteDuplicates // DeleteCases[#, a_ /; !DiffExp`Utilities`DependsQ[a, DiffExp`Symbols`x]] & //
    DeleteDuplicates[#, Expand[#1] === Expand[#2] || Expand[#1] === Expand[-#2] &] &;

  singularityPositions = Quiet[Solve[# == 0, DiffExp`Symbols`x] & /@ AllSingularTerms];
  singularityPositions = Join[(singularityPositions // Flatten)[[All, 2]] // DeleteDuplicatesBy[#, N[#, DiffExp`State`ChopPrecisionVal] &] &,
    {-\[Infinity], \[Infinity]}] // Sort[#, Re[N[#1, DiffExp`State`ChopPrecisionVal]] < Re[N[#2, DiffExp`State`ChopPrecisionVal]] &] &;
  singularityPositionsNumeric = N[singularityPositions, DiffExp`State`FEWorkingPrecision];

  (* Project complex singularities onto the real line *)
  SPosProjections = (
    currentSingularity = #[[1]]; currentSingularityNumeric = #[[2]];

    leftInterval = Select[singularityPositionsNumeric, Re[currentSingularityNumeric] - Im[currentSingularityNumeric] < Re[#] < Re[currentSingularityNumeric] &];
    rightInterval = Select[singularityPositionsNumeric, Re[currentSingularityNumeric] < Re[#] < Re[currentSingularityNumeric] + Im[currentSingularityNumeric] &];

    {
      If[leftInterval === {}, Re[currentSingularity] - Im[currentSingularity], Null],
      Re[currentSingularity],
      If[rightInterval === {}, Re[currentSingularity] + Im[currentSingularity], Null]
    }
  ) & /@ Transpose[{singularityPositions, singularityPositionsNumeric}];

  Quiet[
    projectedSingularities = Sort[
      DeleteDuplicatesBy[DeleteCases[Flatten[SPosProjections], Null], N[#, DiffExp`State`ChopPrecisionVal] &],
      N[#1, DiffExp`State`FEWorkingPrecision] < N[#2, DiffExp`State`FEWorkingPrecision] &
    ];
  ];

  If[getcomplex === True,
    {projectedSingularities, singularityPositions},
    projectedSingularities
  ]
];

(* Gives a large term which can be used to estimate the error *)
GetLargestTerm[line_] := Block[{},
  (Flatten[Coefficient[DiffExp`State`DEqnMatricesExpanded[line] // Values, DiffExp`Symbols`x,
    DiffExp`State`ExpansionOrderVal - DiffExp`State`ISafetyExpansionSubtract - (DiffExp`State`MaxCouplingOrder - 1)]] // Abs // Max) *
    DiffExp`Symbols`x^(DiffExp`State`ExpansionOrderVal - DiffExp`State`ISafetyExpansionSubtract - (DiffExp`State`MaxCouplingOrder - 1))
];

(* Get intervals in the dynamic segmentation strategy *)
GetMatricesPrecisionDistance[line_Association] := Module[{DiffCoeffs, MaxVariation, SmallestCanGo, MaxDisc, DigitsNeeded},
  If[!NumericQ[DiffExp`State`FEC[System`AccuracyGoal]],
    DiffExp`Utilities`ReportError["Accuracy goal is not given as a number."];
  ];

  DigitsNeeded = DiffExp`State`FEAccuracyGoal + DiffExp`State`ISafetyDigits;

  DiffCoeffs = {GetLargestTerm[line]} /. (a_: 1) DiffExp`Symbols`x^(b_: 1) :> {Abs@a, b};

  If[DiffCoeffs === {-\[Infinity]},
    DiffExp`State`LastErrorContext = {line, DiffExp`State`DEqnMatricesExpanded[line]};
    DiffExp`Utilities`ReportError["Could not determine variation in the expanded matrices."];
  ];

  MaxVariation = First[DiffCoeffs];

  SmallestCanGo = Block[{eps = 10^-DigitsNeeded, Const = MaxVariation[[1]], Ord = MaxVariation[[2]]},
    If[OddQ[DiffExp`State`ExpansionOrderVal],
      {
        Root[eps (-1) + Const #1^MaxVariation[[2]] &, 1],
        Root[eps (-1) + Const #1^MaxVariation[[2]] &, 2]
      },
      {Root[eps (-1) + Const #1^MaxVariation[[2]] &, 1]}
    ]
  ];

  SmallestCanGo = Min[Abs /@ SmallestCanGo];

  SmallestCanGo
];

(* Provides consistency checks and reparametrizes boundary conditions *)
CheckBoundaryConditionsAndReparametrize[bcs3_, line_Association] := Module[
  {bcs2, bcs = {Null, Null}, lineRelation, inverseLineRelation, FixAt, zeroLimit, lineRelationSeries, logTerms},

  bcs2 = bcs3;

  If[DiffExp`Utilities`DependsQ[line, DiffExp`Symbols`x^(k_: 1) /; (k != 1 & k != -1)],
    DiffExp`Utilities`ReportError["Non-linear line segments are currently not supported!"];
  ];
  If[Or @@ (!DiffExp`Utilities`DependsQ[Keys@bcs2[[1]], #] & /@ DiffExp`State`ExternalScalesVal),
    DiffExp`State`LastErrorContext = {Keys@bcs2[[1]], DiffExp`State`ExternalScalesVal};
    DiffExp`Utilities`ReportError["The point/line where the boundary conditions are defined does not fix all kinematic invariants and masses!"];
  ];
  If[(DiffExp`Utilities`IsPoint[bcs2[[1]]]) && DiffExp`Utilities`DependsQ[bcs2[[2]] // Normal, DiffExp`Symbols`x],
    DiffExp`Utilities`ReportError["The boundary terms that are provided depend on the line parameter ",
      DiffExp`Symbols`x, ", but the second argument does not."];
  ];
  If[DiffExp`Utilities`IsLine[bcs2[[1]]],
    If[!(DeleteDuplicates@Flatten@MapAt[Head, bcs2[[2]] /. "?" -> O[DiffExp`Symbols`x], {All, All}]) === {SeriesData},
      DiffExp`Utilities`ReportError["Line depends on x but boundary data is not given as a series."];
    ]
  ];
  If[!(Dimensions[bcs2[[2]]][[1]] === DiffExp`State`NumIntegrals),
    DiffExp`Utilities`ReportError["The number of entries in the boundary conditions does not match the size of the system."];
  ];
  If[!(Dimensions[bcs2[[2]]][[2]] >= DiffExp`State`EpsilonOrderVal + 1),
    DiffExp`Utilities`ReportError["Not enough orders in \[Epsilon] are given in the boundary conditions."];
  ];
  bcs2[[2]] = bcs2[[2]][[Range@DiffExp`State`NumIntegrals, Range@(DiffExp`State`EpsilonOrderVal + 1)]];

  lineRelation = RelateLines[line, bcs2[[1]], True];

  If[lineRelation === False,
    DiffExp`State`LastErrorContext = {bcs2, line};
    DiffExp`Utilities`ReportError["Chosen boundary point does not lie on line."];
  ];

  bcs2[[2]] = MapAt[If[!Accuracy[#] === \[Infinity], SetPrecision[#, DiffExp`State`FEWorkingPrecision], #] &, bcs2[[2]], {All, All}];

  zeroLimit = Limit[lineRelation, DiffExp`Symbols`x -> 0] // DiffExp`Utilities`PChop;
  If[zeroLimit === 0,
    FixAt = 0;

    If[!DiffExp`Utilities`DependsQ[bcs2[[1]], DiffExp`Symbols`x],
      (* Point that happens to lie at origin *)
      bcs[[1]] = line /. DiffExp`Symbols`x -> FixAt;
      bcs[[2]] = MapAt[DiffExp`SeriesOps`SeriesAlways[#, {DiffExp`Symbols`x, 0, 0}, 2] &, bcs2[[2]], {All, All}];
      bcs[[2]] = bcs[[2]] /. (a_: 1) "?" + O[DiffExp`Symbols`x]^(1/2) -> "?";,

      (* Line segment around origin *)
      inverseLineRelation = RelateLines[bcs2[[1]], line];

      If[(Sign@Limit[D[inverseLineRelation, DiffExp`Symbols`x], DiffExp`Symbols`x -> 0]) === -1,
        DiffExp`Utilities`ReportError["Asymptotic boundary conditions should be oriented in the same direction as the integration line."];
      ];

      bcs[[1]] = bcs2[[1]] /. DiffExp`Symbols`x -> inverseLineRelation;
      lineRelationSeries = DiffExp`Symbols`x -> (Series[inverseLineRelation, {DiffExp`Symbols`x, 0, DiffExp`State`ISeriesChangeCoefficient DiffExp`State`ExpansionOrderVal},
        Assumptions -> DiffExp`Symbols`x > 0]);
      bcs[[2]] = Assuming[DiffExp`Symbols`x > 0,
        bcs2[[2]] /. lineRelationSeries /. DiffExp`Symbols`Logx -> Log[(lineRelationSeries[[2]] // Normal) /. DiffExp`Symbols`x -> yy]];

      logTerms = bcs[[2]] // DiffExp`Utilities`GetCases[#, Log[_]] &;
      logTerms = # -> Normal[Series[#, {yy, 0, DiffExp`State`ExpansionOrderVal}, Assumptions -> yy > 0]] & /@ logTerms;
      logTerms[[All, 2]] = (logTerms[[All, 2]] /. {
        Log[a_ yy] /; a > 0 :> Log[a] + Log[yy],
        Log[a_ yy] /; a < 0 :> Log[-a] + Log[yy]
      } /. Log[yy] -> DiffExp`Symbols`Logx /. yy -> DiffExp`Symbols`x);

      bcs[[2]] = DiffExp`SeriesOps`SafeReplaceSeries11[bcs[[2]], logTerms] //
        Quiet[N[#, DiffExp`State`FEC[WorkingPrecision]]] &;

      If[DiffExp`Utilities`DependsQ[bcs[[2]], Log[_]],
        DiffExp`State`LastErrorContext = DiffExp`SeriesOps`SExpand@bcs[[2]];
        DiffExp`Utilities`ReportError["Error rescaling boundary conditions."];
      ];
    ];,

    If[DiffExp`Utilities`DependsQ[Normal[lineRelation], DiffExp`Symbols`x],
      DiffExp`Utilities`ReportError["Boundary conditions given on asymptotic limit that is not centered at current line."];
    ];

    (* From hereon we can assume we are fixing at a point *)
    bcs[[1]] = line /. DiffExp`Symbols`x -> lineRelation;
    bcs[[2]] = bcs2[[2]];
    FixAt = lineRelation;
  ];

  {bcs, FixAt}
];

(* Gets matching point on a line *)
GetMatchingPoint[line_Association, bcsline_] := Module[{lineRelation, zeroLimit, FixAt},
  lineRelation = RelateLines[line, bcsline, True];

  If[lineRelation === False,
    DiffExp`State`LastErrorContext = {bcsline, line};
    DiffExp`Utilities`ReportError["GetMatchingPoint: Internal error"];
  ];

  If[DiffExp`Utilities`DependsQ[lineRelation, DiffExp`Symbols`x],
    zeroLimit = Limit[lineRelation, DiffExp`Symbols`x -> 0] // DiffExp`Utilities`PChop;

    If[zeroLimit === 0,
      FixAt = 0;,
      DiffExp`Utilities`ReportError["GetMatchingPoint: Internal error."];
    ];,
    FixAt = lineRelation;
  ];

  FixAt
];

End[];

EndPackage[];
