(* ::Package:: *)

(* DiffExp Pade Subpackage *)
(* This package provides Pade approximant functions *)

BeginPackage["DiffExp`Pade`", {"DiffExp`Symbols`", "DiffExp`State`", "DiffExp`Utilities`", "DiffExp`SeriesOps`"}];

(* Pade functions *)
GetPade::usage = "GetPade[a_] returns the Pade approximant of a series.";
SEval1::usage = "SEval1[a_] evaluates series, optionally using Pade approximant.";
SEval2::usage = "SEval2[a_,at_] evaluates expression at a point with analytic continuation.";
SEval::usage = "SEval[a_,at_] evaluates series at a point.";
ToPiecewise::usage = "ToPiecewise[segmentdata_List, pade_:False] converts segment data to piecewise functions.";

Begin["`Private`"];

(* GetPade for various inputs *)
GetPade[0] = 0;
GetPade[a_?NumericQ] := a;
GetPade[a_SeriesData] := Block[
  {
    $MinPrecision = DiffExp`State`FEWorkingPrecision,
    maxPadeOrder = {Floor[((a[[5]] - a[[4]])/a[[6]] + 1)/2], Floor[((a[[5]] - a[[4]])/a[[6]] + 1)/2]},
    MaxLogOrder = DiffExp`SeriesOps`MaxLogxPower[a],
    tempPadeApprox
  },

  Sum[
    DiffExp`Symbols`Logx^ind (
      (Quiet[
        tempPadeApprox = PadeApproximant[Chop[#, 10^-DiffExp`State`ChopPrecisionVal], {DiffExp`Symbols`x, 0, maxPadeOrder}];
        If[Length[DiffExp`Utilities`GetCases[tempPadeApprox, PadeApproximant[__]]] > 0,
          DiffExp`Utilities`PrintWarning["Error in Pade approximant of ", # // N, ". Evaluating normally."][1];
          DiffExp`Utilities`PChop@#,
          tempPadeApprox
        ]
      ]) &@(Normal@DiffExp`SeriesOps`LogxCoeff[a, ind])
    ),
    {ind, 0, MaxLogOrder}
  ]
];

(* SEval1 - get evaluable form of series *)
SetAttributes[SEval1, Listable];
SEval1[a_SeriesData] := If[DiffExp`State`FEC[UsePade] === True, GetPade[a], Normal[a]];
SEval1[0] := 0;
SEval1[a_ /; NumericQ[a]] := a;

(* SEval2 - evaluate at a point with analytic continuation *)
SetAttributes[SEval2, Listable];
SEval2[a_, at_] := Block[{$MinPrecision = DiffExp`State`FEWorkingPrecision},
  a /. (DiffExp`State`AnalyticContinuationReplacements /.
    If[at >= 0,
      {DiffExp`Symbols`\[Theta]p -> 1, DiffExp`Symbols`\[Theta]m -> 0},
      {DiffExp`Symbols`\[Theta]p -> 0, DiffExp`Symbols`\[Theta]m -> 1}
    ]) /. DiffExp`Symbols`Logx -> Log[DiffExp`Symbols`x] /. DiffExp`Symbols`x -> at //
    SetPrecision[#, DiffExp`State`FEWorkingPrecision] & // Expand
];

(* SEval - evaluate series at a point *)
SEval[a_SeriesData, at_] := Block[{$MinPrecision = DiffExp`State`FEWorkingPrecision},
  If[DiffExp`State`FEC[UsePade] === True, GetPade[a], Normal[a]] /.
    (DiffExp`State`AnalyticContinuationReplacements /.
      If[at >= 0,
        {DiffExp`Symbols`\[Theta]p -> 1, DiffExp`Symbols`\[Theta]m -> 0},
        {DiffExp`Symbols`\[Theta]p -> 0, DiffExp`Symbols`\[Theta]m -> 1}
      ]) /. DiffExp`Symbols`Logx -> Log[DiffExp`Symbols`x] /. DiffExp`Symbols`x -> at //
    SetPrecision[#, DiffExp`State`FEWorkingPrecision] & // Expand
];
SEval[0, at_] := 0;
SEval[a_ /; NumericQ[a], at_] := a;

(* ToPiecewise - convert saved segment data to piecewise functions *)
ToPiecewise[SavedData2_, Pade : _?BooleanQ : False, Ord_Integer : Null] := Module[
  {SavedData, piecewiseResult, Uncompressed, Counter},

  If[MatchQ[SavedData2, {{a_Association, _}, {__}}] || MatchQ[SavedData2, {{a_Association, _, _}, {__}}],
    SavedData = SavedData2[[2]],
    SavedData = SavedData2
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
          (DiffExp`AnalyticContinuation`Project\[Theta]s[#, GetPade] &@#) /.
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
