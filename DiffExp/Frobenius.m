(* ::Package:: *)

(* DiffExp Frobenius Subpackage *)
(* This package provides Frobenius solution methods *)

BeginPackage["DiffExp`Frobenius`", {"DiffExp`Symbols`", "DiffExp`State`", "DiffExp`Utilities`", "DiffExp`SeriesOps`", "DiffExp`Integration`"}];

(* Frobenius functions *)
Frobenius1::usage = "Frobenius1[DEqn_] finds one Frobenius solution.";
FrobeniusSolutions::usage = "FrobeniusSolutions[DEqn_] finds all Frobenius solutions.";

Begin["`Private`"];

(* Solution from the Frobenius ansatz *)
Frobenius1[DEqn_] := Module[
  {DEqnSer, equationLength = Length[DEqn], IndicialEquation, rMax, Ansatz, c, Eqns, cUnknowns, cSols, r, validatedSolution, crossCheckResult, maxCoefficientIndex, LeadingPowers},

  LeadingPowers = (#[[4]]/#[[6]] &@DiffExp`SeriesOps`LeadingCoefficientSeries[#, 2]) & /@ DEqn;
  LeadingPowers = LeadingPowers - Last[LeadingPowers];

  Do[
    (* Check that equation is Fuchsian *)
    If[LeadingPowers[[ind]] < ind - Length[LeadingPowers],
      DiffExp`Utilities`ReportError["Obtained a higher order differential equation that is not in Fuchsian form."];
    ];,
    {ind, Length@LeadingPowers}
  ];

  DiffExp`Utilities`PrintDebug["Frobenius method: Series expanding.."][3];
  DEqnSer = DiffExp`SeriesOps`DiffExpSeries[DEqn];

  IndicialEquation = Sum[
    DiffExp`SeriesOps`LeadingCoefficientSeries[DEqnSer[[ind]]] *
      (D[DiffExp`Symbols`x^r, Sequence @@ ConstantArray[DiffExp`Symbols`x, ind - 1]] DiffExp`Symbols`x^-r),
    {ind, equationLength}
  ] // Normal // Rationalize[#, 10^-DiffExp`State`ChopPrecisionVal] &;

  DiffExp`Utilities`PrintDebug["Indicial equation: ", IndicialEquation /. r -> "r"][3];
  rMax = Max[IndicialEquation == 0 // Solve[#, r] & // Part[#, All, 1, 2] &];
  DiffExp`Utilities`PrintDebug["Considering largest root: r = ", rMax][3];

  If[Denominator[rMax] > 2 && !(DiffExp`State`FEC["IgnoreIndicialCheck"] === True),
    DiffExp`Utilities`PrintWarning["The root of the indicial equation is of degree greater than two: ", rMax,
      ". This case is not thoroughly tested."];
  ];

  cUnknowns = Table[Subscript[c, iind], {iind, Ceiling[DiffExp`State`ExpansionOrderVal - rMax]}];
  Ansatz = Expand[DiffExp`Symbols`x^rMax Sum[
    Subscript[c, iind] DiffExp`Symbols`x^iind /. Subscript[c, 0] -> 1,
    {iind, 0, Ceiling[DiffExp`State`ExpansionOrderVal - rMax]}
  ]] + O[DiffExp`Symbols`x]^Ceiling[DiffExp`State`ExpansionOrderVal];

  Eqns = DiffExp`Utilities`PChop[# == 0 & /@ Sum[
    DEqnSer[[ind]] (D[Ansatz, Sequence @@ ConstantArray[DiffExp`Symbols`x, ind - 1]]),
    {ind, equationLength}
  ][[3]]];

  Check[
    cSols = N[Thread[cUnknowns -> LinearSolve[
      Sequence @@ {#[[2]], -#[[1]]} &@CoefficientArrays[DeleteCases[Eqns, True], cUnknowns],
      ZeroTest -> (N[DiffExp`Utilities`LSPChop@Expand@Normal[#1], DiffExp`State`LinearSolveChopPrecisionVal] == 0 &)
    ]], DiffExp`State`FEWorkingPrecision];,

    DiffExp`State`LastErrorContext = {DEqn, DEqnSer, IndicialEquation, rMax, Ansatz, cUnknowns, Eqns, validatedSolution, cSols};

    If[DiffExp`Utilities`DependsQ[cSols, LinearSolve],
      DiffExp`Utilities`ReportError["Something went wrong while applying the Frobenius method."];,
      DiffExp`Utilities`PrintWarning["Encountered possible instability during evaluation of the Frobenius method. Cross-checking the result."];
      crossCheckResult = Ansatz /. (DiffExp`Utilities`PChop[cSols]);
      crossCheckResult = Normal[DiffExp`SeriesOps`SExpand@Sum[
        DEqnSer[[ord]] D[crossCheckResult, Sequence @@ ConstantArray[DiffExp`Symbols`x, ord - 1]],
        {ord, Length[DEqnSer]}
      ]];
      If[crossCheckResult === 0,
        DiffExp`Utilities`PrintInfo["Solution is valid. Continuing."][1];,
        DiffExp`Utilities`ReportError["Result is incorrect. Aborting."];
      ];
    ];
  ];

  DiffExp`Utilities`PrintDebug["Frobenius solution found."][3];

  maxCoefficientIndex = Max[Cases[Eqns, Subscript[c, i_] :> i, Infinity]];
  Series[Ansatz, {DiffExp`Symbols`x, 0, Floor[maxCoefficientIndex + rMax]}] /. (DiffExp`Utilities`PChop[cSols])
];

(* All Frobenius solutions *)
FrobeniusSolutions[DEqn_] := Module[
  {DEqnOrder = Length[DEqn] - 1, Solns, DEqnReduced, DEqnSer, DEqnZeros, Checks},

  DEqnSer = DiffExp`SeriesOps`DiffExpSeries[DEqn, DiffExp`State`ExpansionOrderVal];

  DEqnZeros = If[Normal[#] === 0, 0, 1] & /@ DEqnSer;

  Solns = {Frobenius1[DEqnSer]};

  If[DEqnOrder > 1,
    DEqnReduced = Table[
      DiffExp`SeriesOps`SExpand[Sum[
        Binomial[iind - 1, kind - 1] DiffExp`SeriesOps`SMultiply[
          DEqnSer[[iind]],
          DiffExp`SeriesOps`SD[Solns[[1]], Sequence @@ DiffExp`Utilities`CA[DiffExp`Symbols`x, iind - kind]]
        ],
        {iind, kind, DEqnOrder + 1}
      ]],
      {kind, DEqnOrder + 1}
    ];

    (* Cross-check *)
    If[MemberQ[DiffExp`State`CurrCrosscheckFlags, "FrobeniusSolutions"] === True,
      If[!Normal[(DiffExp`Utilities`CPChop@DEqnReduced[[1]])] === 0,
        DiffExp`Utilities`PrintDebug[DEqnReduced[[1]] + O[DiffExp`Symbols`x]^DiffExp`State`ExpansionOrderVal // DiffExp`SeriesOps`SN][3];
        DiffExp`Utilities`ReportError["Something went wrong while deriving a lower order differential equation from a multiplicative ansatz."];
      ];
    ];

    DEqnReduced = Delete[DEqnReduced, 1];
    Solns = Join[
      Solns,
      DiffExp`SeriesOps`SMultiply[Solns[[1]], #] & /@ DiffExp`Integration`DiffExpIntegrate[FrobeniusSolutions[DEqnReduced]]
    ];
  ];

  If[MemberQ[DiffExp`State`CurrCrosscheckFlags, "FrobeniusSolutions"] === True,
    DiffExp`Utilities`PrintDebug["Running cross-checks on Frobenius solutions.."][1];
    Checks = (DiffExp`Utilities`CPChop@*DiffExp`SeriesOps`SExpand)[Sum[
      DiffExp`SeriesOps`SMultiply[DEqnSer[[iind]], DiffExp`SeriesOps`SD[#, Sequence @@ DiffExp`Utilities`CA[DiffExp`Symbols`x, iind - 1]]],
      {iind, DEqnOrder + 1}
    ]] & /@ Solns;
    DiffExp`Utilities`PrintDebug["Cross-checks: ", Checks // DiffExp`SeriesOps`SN][3];
    If[!SameQ[Append[Normal[Checks], 0]],
      DiffExp`Utilities`ReportError["Error while validating the Frobenius solutions. Aborting.."];
    ];
  ];

  Solns
];

End[];

EndPackage[];
