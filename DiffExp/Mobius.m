(* ::Package:: *)

(* DiffExp Mobius Subpackage *)
(* This package provides Mobius transformation functions *)

BeginPackage["DiffExp`Mobius`", {"DiffExp`Symbols`", "DiffExp`State`", "DiffExp`Utilities`"}];

(* Mobius transformation functions *)
GetMobius::usage = "GetMobius[{zmin_,zmid_,zmax_}] returns Mobius transformation.";
GetLineRescaled::usage = "GetLineRescaled[line_,at_,{signsproj_,signsim_},nomobius_] rescales line around a point.";
GetMobiusCPL::usage = "GetMobiusCPL[...] gets Mobius center point left.";
GetMobiusCPR::usage = "GetMobiusCPR[...] gets Mobius center point right.";
GetCPLRep::usage = "GetCPLRep[MyEq_] gets center point replacement.";
GetCPL::usage = "GetCPL[...] gets center point left.";
GetCPR::usage = "GetCPR[...] gets center point right.";
FindNextCenterPointL::usage = "FindNextCenterPointL[xbc_,singsproj_] finds next center point to the left.";
FindNextCenterPointR::usage = "FindNextCenterPointR[xbc_,singsproj_] finds next center point to the right.";

Begin["`Private`"];

(* GetMobius for different singularity configurations *)
GetMobius[{-\[Infinity], zmid_, \[Infinity]}] := DiffExp`Symbols`x - zmid;

GetMobius[{zmin_, zmid_, \[Infinity]}] /; zmin != -\[Infinity] :=
  (zmid + DiffExp`Symbols`x zmid - 2 DiffExp`Symbols`x zmin)/(1 - DiffExp`Symbols`x) // Simplify;

GetMobius[{-\[Infinity], zmid_, zmax_}] /; zmax != \[Infinity] :=
  (2 DiffExp`Symbols`x zmax + zmid - DiffExp`Symbols`x zmid)/(1 + DiffExp`Symbols`x) // Simplify;

GetMobius[{zmin_, zmid_, zmax_}] /; zmax != \[Infinity] && zmin != -\[Infinity] :=
  -((-1 + DiffExp`Symbols`x) zmid zmin + zmax (zmid + DiffExp`Symbols`x zmid - 2 DiffExp`Symbols`x zmin))/
   ((-1 + DiffExp`Symbols`x) zmax + zmin + DiffExp`Symbols`x (-2 zmid + zmin)) // Simplify;

GetMobius[{zmin_, \[Infinity], zmax_}] /; zmax != \[Infinity] && zmin != -\[Infinity] :=
  (zmin (-1 + DiffExp`Symbols`x) + zmax (1 + DiffExp`Symbols`x))/(2 DiffExp`Symbols`x);

(* GetLineRescaled - get a line segment centered at at_ *)
GetLineRescaled[line_Association, at_, {signsproj_, signsim_}, nomobius_: False] :=
  Block[{$MaxExtraPrecision = 1000},
    If[DiffExp`State`UseMobiusVal === True && nomobius === False,
      Module[{leftBound, rightBound},
        leftBound = Select[signsproj, # < at &] // Last;
        rightBound = Select[signsproj, # > at &] // First;

        (Collect[Numerator[#], DiffExp`Symbols`x]/Collect[Denominator[#], DiffExp`Symbols`x]) & /@
          ((line /. DiffExp`Symbols`x -> GetMobius[{leftBound, at, rightBound}])) /.
          DiffExp`Symbols`x -> DiffExp`Symbols`x/DiffExp`State`FEC[RadiusOfConvergence] //
          SetPrecision[#, 2 DiffExp`State`FEWorkingPrecision] &
      ],
      Module[{minDistance, leftBound, rightBound},
        leftBound = Select[signsproj, # < at &] // Last;
        rightBound = Select[signsproj, # > at &] // First;
        minDistance = Min[at - leftBound, rightBound - at];
        (* Handle infinite bounds: mirror the finite side, or use unit distance *)
        If[!TrueQ[minDistance < Infinity],
          minDistance = Which[
            leftBound =!= -Infinity, at - leftBound,
            rightBound =!= Infinity, rightBound - at,
            True, 1
          ];
        ];

        (Collect[Numerator[#], DiffExp`Symbols`x]/Collect[Denominator[#], DiffExp`Symbols`x]) & /@
          (((Normal[line] /. DiffExp`Symbols`x -> at + DiffExp`Symbols`x # /.
            DiffExp`Symbols`x -> DiffExp`Symbols`x/DiffExp`State`FEC[RadiusOfConvergence] //
            SetPrecision[#, 2 DiffExp`State`FEWorkingPrecision] & // Expand // Association) &@minDistance))
      ]
    ]
  ];

(* Functions for deriving center points in the predivision segmentation strategy *)
GetMobiusCPL[{-\[Infinity], zbound_, zmax_}] :=
  (-zbound + k zbound + 2 zmax)/(1 + k) /. k -> DiffExp`State`FEC[DivisionOrder];

GetMobiusCPL[{zmin_, zbound_, \[Infinity]}] :=
  (zbound + k zbound - 2 zmin)/(-1 + k) /. k -> DiffExp`State`FEC[DivisionOrder];

GetMobiusCPL[{zmin_, zbound_, zmax_}] :=
  (zbound zmax + k zbound zmax + zbound zmin - k zbound zmin - 2 zmax zmin)/
   (2 zbound - zmax + k zmax - zmin - k zmin) /. k -> DiffExp`State`FEC[DivisionOrder];

GetMobiusCPR[{-\[Infinity], zbound_, zmax_}] :=
  (zbound + k zbound - 2 zmax)/(-1 + k) /. k -> DiffExp`State`FEC[DivisionOrder];

GetMobiusCPR[{zmin_, zbound_, \[Infinity]}] :=
  (-zbound + k zbound + 2 zmin)/(1 + k) /. k -> DiffExp`State`FEC[DivisionOrder];

GetMobiusCPR[{zmin_, zbound_, zmax_}] :=
  (zbound zmax - k zbound zmax + zbound zmin + k zbound zmin - 2 zmax zmin)/
   (2 zbound - zmax - k zmax - zmin + k zmin) /. k -> DiffExp`State`FEC[DivisionOrder];

(* Helper function for center point calculation *)
GetCPLRep[MyEq_] := Block[{ruleSet},
  ruleSet = Cases[MyEq, xnew == a_ | s == a_] /. Equal -> Rule;
  xnew //. ruleSet
];

GetCPL[{-\[Infinity], zbound_, zmax_}, k2_: Null] := Block[
  {MyEq, k = If[k2 === Null, DiffExp`State`DivisionOrderVal, k2]},
  MyEq = zmax \[Element] Reals && zbound < zmax && k >= 1 && s == (-k zbound + k zmax)/(1 + k) && xnew == -s + zmax;
  GetCPLRep[MyEq]
];

GetCPL[{zmin_, zbound_, \[Infinity]}, k2_: Null] := Block[
  {MyEq, k = If[k2 === Null, DiffExp`State`DivisionOrderVal, k2]},
  MyEq = zmin \[Element] Reals && zbound > zmin && k > 1 && s == (k zbound - k zmin)/(-1 + k) && xnew == s + zmin;
  GetCPLRep[MyEq]
];

GetCPL[{zmin_, zbound_, zmax_}, k2_: Null] := Block[
  {MyEq, k = If[k2 === Null, DiffExp`State`DivisionOrderVal, k2]},
  MyEq = zmin \[Element] Reals && zmax > zmin &&
    ((zmin < zbound < (zmax + zmin)/2 &&
      ((1 <= k <= (zmax - zmin)/(-2 zbound + zmax + zmin) && s == (-k zbound + k zmax)/(1 + k)) ||
       (k > (zmax - zmin)/(-2 zbound + zmax + zmin) && s == (k zbound - k zmin)/(-1 + k)))) ||
     ((zmax + zmin)/2 <= zbound < zmax && k >= 1 && s == (-k zbound + k zmax)/(1 + k))) &&
    xnew == (s + k zbound)/k;
  GetCPLRep[MyEq]
];

GetCPR[{-\[Infinity], zbound_, zmax_}, k2_: Null] := Block[
  {MyEq, k = If[k2 === Null, DiffExp`State`DivisionOrderVal, k2]},
  MyEq = zmax \[Element] Reals && zbound < zmax && k > 1 && s == (-k zbound + k zmax)/(-1 + k) && xnew == -s + zmax;
  GetCPLRep[MyEq]
];

GetCPR[{zmin_, zbound_, \[Infinity]}, k2_: Null] := Block[
  {MyEq, k = If[k2 === Null, DiffExp`State`DivisionOrderVal, k2]},
  MyEq = zmin \[Element] Reals && zbound > zmin && k >= 1 && s == (k zbound - k zmin)/(1 + k) && xnew == s + zmin;
  GetCPLRep[MyEq]
];

GetCPR[{zmin_, zbound_, zmax_}, k2_: Null] := Block[
  {MyEq, k = If[k2 === Null, DiffExp`State`DivisionOrderVal, k2]},
  MyEq = zmin \[Element] Reals && zmax > zmin &&
    ((zmin < zbound <= (zmax + zmin)/2 && k >= 1 && s == (k zbound - k zmin)/(1 + k)) ||
     ((zmax + zmin)/2 < zbound < zmax &&
      ((1 <= k <= (-zmax + zmin)/(-2 zbound + zmax + zmin) && s == (k zbound - k zmin)/(1 + k)) ||
       (k > (-zmax + zmin)/(-2 zbound + zmax + zmin) && s == (-k zbound + k zmax)/(-1 + k))))) &&
    xnew == -((s - k zbound)/k);
  GetCPLRep[MyEq]
];

(* Find next center points *)
FindNextCenterPointL[xbc_, singsproj_] := Module[{leftBound, rightBound},
  leftBound = Select[singsproj, # < xbc &] // Last;
  rightBound = Select[singsproj, # > xbc &] // First;
  If[DiffExp`State`FEC[UseMobius] === True, GetMobiusCPL[{leftBound, xbc, rightBound}], GetCPL[{leftBound, xbc, rightBound}]]
];

FindNextCenterPointR[xbc_, singsproj_] := Module[{leftBound, rightBound},
  leftBound = Select[singsproj, # < xbc &] // Last;
  rightBound = Select[singsproj, # > xbc &] // First;
  If[DiffExp`State`FEC[UseMobius] === True, GetMobiusCPR[{leftBound, xbc, rightBound}], GetCPR[{leftBound, xbc, rightBound}]]
];

End[];

EndPackage[];
