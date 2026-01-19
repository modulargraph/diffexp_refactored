(* ::Package:: *)

(* DiffExp AnalyticContinuation Subpackage *)
(* This package handles analytic continuation logic *)

BeginPackage["DiffExp`AnalyticContinuation`", {"DiffExp`Symbols`", "DiffExp`State`", "DiffExp`Utilities`", "DiffExp`SeriesOps`"}];

(* Analytic continuation functions *)
PrepareAnalyticContinuation::usage = "PrepareAnalyticContinuation[Line_] prepares replacement rules for analytic continuation.";
Project\[Theta]s::usage = "Project\[Theta]s[Expr_,f_] removes multiple occurrences of theta functions.";

Begin["`Private`"];

(* Set attributes as in original *)
SetAttributes[Project\[Theta]s, Listable];

(* Prepares replacement rules for analytic continuation along the given line segment *)
PrepareAnalyticContinuation[Line_] := Module[
  {leadingCoeffs, indexedTerms, SignsNeeded, uniqueSigns, CurrSignChoices, FindVanishingFactors, AllPrescriptions, VanishingFactors},

  FindVanishingFactors[Factors_] := (
    leadingCoeffs = {Normal[DiffExp`SeriesOps`LeadingCoefficientSeries[#[[1]]]], #[[2]]} & /@
      (Factors /. Line // Expand // DiffExp`Utilities`PChop);
    indexedTerms = Table[{iind, leadingCoeffs[[iind]]}, {iind, Length[leadingCoeffs]}];
    indexedTerms = Select[indexedTerms, DiffExp`Utilities`DependsQ[#[[2, 1]], DiffExp`Symbols`x] &];
    VanishingFactors = indexedTerms[[All, 2]];
  );

  (* Check whether the singularity coincides with an automatically added +I\[Delta] coming from a root in the differential equations *)
  DiffExp`State`CurrentSingularityWasAddedFromSquareRoot = False;
  FindVanishingFactors[DiffExp`State`SquareRootPrescriptionsAdded[]];
  If[Length[VanishingFactors] > 0,
    FindVanishingFactors[Complement[DiffExp`State`DeltaPrescriptionsVal, DiffExp`State`SquareRootPrescriptionsAdded[]]];
    If[Length[VanishingFactors] === 0,
      DiffExp`State`CurrentSingularityWasAddedFromSquareRoot = True;
    ];
  ];

  FindVanishingFactors[DiffExp`State`DeltaPrescriptionsVal];
  If[Length[VanishingFactors] > 0,
    DiffExp`State`CurrentSingularityHasIDeltaPrescription = True;,
    DiffExp`State`CurrentSingularityHasIDeltaPrescription = False
  ];

  SignsNeeded = (
    {Exponent[#[[1]], DiffExp`Symbols`x], #[[2]]/(Sign[#[[1]]] /. Sign[DiffExp`Symbols`x] -> 1)} /. {
      {k_, 1} :> (k /. {1 -> 1, a_ -> "?"}),
      {k_, -1} :> (k /. {1 -> -1, a_ -> "?"})
    }
  ) & /@ VanishingFactors;
  uniqueSigns = DeleteDuplicates@SignsNeeded;

  CurrSignChoices = Thread[List[DiffExp`State`DeltaPrescriptionsVal[[indexedTerms[[All, 1]]]][[All, 1]], SignsNeeded]];

  If[(DiffExp`Utilities`DependsQ[uniqueSigns, "?"]) || Length[uniqueSigns] > 1,
    If[!KeyExistsQ[DiffExp`State`AnalyticContinuationReplacementsAssociation, Line],
      DiffExp`Utilities`PrintInfo["Singularity => Sign[Im[x]]:"][2];
      DiffExp`Utilities`PrintInfo[
        Rule @@ # & /@ (Thread[List[
          DiffExp`State`DeltaPrescriptionsVal[[indexedTerms[[All, 1]]]][[All, 1]] +
            I "\[Delta]" DiffExp`State`DeltaPrescriptionsVal[[indexedTerms[[All, 1]]]][[All, 2]],
          SignsNeeded
        ]])
      ][2];
    ];
    DiffExp`State`AnalyticContinuationFailed = True;
    uniqueSigns = {1};
  ];

  DiffExp`State`AnalyticContinuationReplacements = {};
  If[uniqueSigns === {-1},
    DiffExp`State`AnalyticContinuationReplacements = {
      DiffExp`Symbols`Logx -> (DiffExp`Symbols`\[Theta]p + DiffExp`Symbols`\[Theta]m) DiffExp`Symbols`Logx - 2 \[Pi] I DiffExp`Symbols`\[Theta]m,
      DiffExp`Symbols`x^b_ /; Denominator[b] == 2 :>
        DiffExp`Symbols`x^(b - 1/2) (DiffExp`Symbols`\[Theta]p - DiffExp`Symbols`\[Theta]m) Sqrt[DiffExp`Symbols`x],
      DiffExp`Symbols`x^b_ /; Denominator[b] > 2 :>
        (DiffExp`Symbols`\[Theta]p + Exp[-2 \[Pi] I b] DiffExp`Symbols`\[Theta]m) DiffExp`Symbols`x^b
    };
  ];

  If[!KeyExistsQ[DiffExp`State`AnalyticContinuationReplacementsAssociation, Line],
    If[uniqueSigns === {},
      DiffExp`Utilities`PrintDebug["Line not centered at a singularity."][1];,
      DiffExp`Utilities`PrintDebug["Analytic continuation: x carries ", uniqueSigns[[1]], "*i\[Delta]"][1];
      DiffExp`Utilities`PrintDebug["Using replacement rules: ", DiffExp`State`AnalyticContinuationReplacements][1];
    ];
  ];

  DiffExp`State`AnalyticContinuationReplacementsAssociation[Line] = DiffExp`State`AnalyticContinuationReplacements;
];

(* Removes multiple occurrences of theta functions *)
Project\[Theta]s[Expr_, f_: Expand] := Module[
  {\[Theta]mPart, \[Theta]pPart},

  If[Length[DiffExp`Utilities`GetCases[Expr, DiffExp`Symbols`\[Theta]p | DiffExp`Symbols`\[Theta]m]] === 0,
    f[Expr] // DiffExp`SeriesOps`SExpand,
    \[Theta]mPart = (Expr /. DiffExp`Symbols`\[Theta]p -> 0 /. DiffExp`Symbols`\[Theta]m -> 1);
    \[Theta]pPart = (Expr /. DiffExp`Symbols`\[Theta]m -> 0 /. DiffExp`Symbols`\[Theta]p -> 1);
    (f[\[Theta]pPart + \[Theta]mPart] // (DiffExp`Utilities`PChop@*DiffExp`SeriesOps`SExpand)) (1/2) +
      (f[\[Theta]pPart - \[Theta]mPart] // (DiffExp`Utilities`PChop@*DiffExp`SeriesOps`SExpand)) (1/2 DiffExp`Symbols`\[Theta]p - 1/2 DiffExp`Symbols`\[Theta]m)
  ]
];

End[];

EndPackage[];
