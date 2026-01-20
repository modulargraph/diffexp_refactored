(* ::Package:: *)

(* DiffExp - Main Package *)
(* This is the refactored version that loads modular subpackages *)

(* Copyright (C) 2024 Martijn Hidding *)

(* DiffExp is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>. *)

(* This software is maintained on https://gitlab.com/hiddingm/diffexp *)

(* Load all subpackages explicitly before BeginPackage *)
Get["Symbols.m"];
Get["State.m"];
Get["Utilities.m"];
Get["SeriesOps.m"];
Get["Integration.m"];
Get["Pade.m"];
Get["Mobius.m"];
Get["AnalyticContinuation.m"];
Get["LineSegmentation.m"];
Get["Frobenius.m"];
Get["Wronskian.m"];
Get["MatrixLoading.m"];
Get["IntegrationStrategies.m"];
Get["Transport.m"];
Get["SingularityDecomposition.m"];
Get["RegularizedIntegration.m"];

(* Begin package with imported subpackage contexts *)
BeginPackage["DiffExp`", {
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
  "DiffExp`IntegrationStrategies`",
  "DiffExp`Transport`",
  "DiffExp`SingularityDecomposition`",
  "DiffExp`RegularizedIntegration`"
}];

Print["Loading DiffExp version 1.1 (refactored)"];
Print["For questions, email: martijnhidding@outlook.com"];
Print["For the latest version, see: https://gitlab.com/hiddingm/diffexp"];

(* Note: Symbols like \[Epsilon], eps, Logx, \[Theta]p, \[Theta]m are available
   from DiffExp`Symbols` context which is imported above.
   Functions like TransportTo, IntegrateSystem, PrepareBoundaryConditions are
   available from DiffExp`Transport` context.
   Configuration functions CurrentConfiguration, LoadConfiguration, UpdateConfiguration
   are defined below but use DiffExp`State` for storage.
   ToPiecewise is available from DiffExp`Pade` context.
*)

(* Configuration functions - these are defined here in DiffExp` context *)
CurrentConfiguration::usage = "CurrentConfiguration[] returns the current configuration values.";
LoadConfiguration::usage = "LoadConfiguration[l_List] can be used to parse configuration options to DiffExp.";
UpdateConfiguration::usage = "UpdateConfiguration[l_List] or UpdateConfiguration[l__Rule] can be used to update configuration options.";

Begin["`Private`"];

(* CurrentConfiguration *)
CurrentConfiguration[] := KeySelect[DiffExp`State`DiffExpConfiguration, MemberQ[{
  System`AccuracyGoal, ChopPrecision, DeltaPrescriptions,
  DivisionOrder, EpsilonOrder, ExpansionOrder,
  LineParameter, MatrixDirectory, RadiusOfConvergence,
  SegmentationStrategy, IntegrationStrategy, UseMobius,
  UsePade, System`Variables, Verbosity, System`WorkingPrecision
}, #] &];

(* LoadConfiguration *)
LoadConfiguration[a__] := (
  DiffExp`State`DiffExpConfiguration = DiffExp`State`DefaultConfiguration;
  #[[1]][] & /@ DiffExp`State`DiffExpExtensions;
  UpdateConfiguration[a]
);

(* UpdateConfiguration *)
UpdateConfiguration[a__Rule] := UpdateConfiguration[{a}];
UpdateConfiguration[l_List] := UpdateConfiguration[l // Association];
UpdateConfiguration[assoc_Association] := Module[{PSFL, DeltaPrescriptionsAdjusted},
  DiffExp`State`DiffExpConfiguration = Merge[{DiffExp`State`DiffExpConfiguration, assoc}, Last];

  If[KeyExistsQ[assoc, LogFile],
    If[ValueQ[DiffExp`State`LogStream],
      Close[DiffExp`State`LogStream];
    ];
    DiffExp`State`LogStream = OpenAppend[assoc[LogFile]];
    AppendTo[$Output, DiffExp`State`LogStream];
    AppendTo[$Messages, DiffExp`State`LogStream];
  ];
  If[KeyExistsQ[DiffExp`State`FEC, "CrosscheckLevel"],
    DiffExp`State`CurrCrosscheckFlags = Select[DiffExp`State`CrosscheckFlags // Normal, #[[2]] <= DiffExp`State`FEC["CrosscheckLevel"] &][[All, 1]];
  ];
  If[KeyExistsQ[DiffExp`State`FEC, "CrosscheckFlags"],
    DiffExp`State`CurrCrosscheckFlags = Union[DiffExp`State`CurrCrosscheckFlags, DiffExp`State`FEC["CrosscheckFlags"]];
  ];
  If[KeyExistsQ[assoc, System`WorkingPrecision] || KeyExistsQ[assoc, ChopPrecision],
    If[DiffExp`State`FEC[ChopPrecision] >= DiffExp`State`FEC[WorkingPrecision],
      DiffExp`Utilities`ReportError["The value of ChopPrecision should be smaller than the value of WorkingPrecision."];
    ];
  ];
  DiffExp`State`DiffExpConfiguration["LinearSolveChopPrecision"] = DiffExp`State`FEC[ChopPrecision];

  If[KeyExistsQ[assoc, "LinearSolveChopPrecision"],
    DiffExp`State`DiffExpConfiguration["LinearSolveChopPrecision"] = assoc["LinearSolveChopPrecision"];
  ];
  If[KeyExistsQ[assoc, LineParameter] && KeyExistsQ[DiffExp`State`FEC, System`Variables],
    If[MemberQ[DiffExp`State`FEC[System`Variables], assoc[LineParameter]],
      DiffExp`Utilities`ReportError["The symbol for the line parameter can't be equal to one of the kinematic variables or masses."];
    ];
  ];
  If[KeyExistsQ[assoc, LineParameter],
    DiffExp`State`LineParameterVal = assoc[LineParameter];
    DiffExp`Symbols`x = assoc[LineParameter];
  ];
  If[KeyExistsQ[assoc, DeltaPrescriptions],
    DeltaPrescriptionsAdjusted = If[#[[0]] === List, #,
      {# /. Global`\[Delta] -> 0, Coefficient[#, Global`\[Delta]]/I}
    ] & /@ assoc[DeltaPrescriptions];
    DiffExp`State`DiffExpConfiguration[DeltaPrescriptions] = DeltaPrescriptionsAdjusted;

    (* We separately keep track of the delta prescriptions added by the user *)
    DiffExp`State`UserDeltaPrescriptions = DeltaPrescriptionsAdjusted;
  ];

  If[KeyExistsQ[assoc, "EstimateError"],
    If[assoc["EstimateError"] === "False",
      DiffExp`State`DiffExpConfiguration["EstimateError"] = False;
      ,
      If[!MemberQ[{True, "Fast"}, assoc["EstimateError"]],
        DiffExp`State`DiffExpConfiguration["EstimateError"] = False;
      ];
    ];
  ];

  If[KeyExistsQ[assoc, ExpansionOrder],
    DiffExp`State`ExpansionOrderVal = DiffExp`State`DiffExpConfiguration[ExpansionOrder];
  ];

  If[KeyExistsQ[assoc, DeltaPrescriptions],
    PSFL = FactorList /@ DiffExp`State`FEC[DeltaPrescriptions][[All, 1]];
    If[DiffExp`Utilities`DependsQ[Length[#] > 2 & /@ PSFL, True],
      DiffExp`Utilities`ReportError["Physical singularities should be irreducible polynomials!"];
    ];
  ];

  If[KeyExistsQ[assoc, MatrixDirectory] || (KeyExistsQ[DiffExp`State`FEC, MatrixDirectory] && KeyExistsQ[assoc, EpsilonOrder]),
    If[(!KeyExistsQ[assoc, System`Variables]) && KeyExistsQ[assoc, MatrixDirectory], DiffExp`State`DiffExpConfiguration[System`Variables] = {};];

    DiffExp`MatrixLoading`LoadMatrices[DiffExp`State`MatrixDirectoryVal];
    ,
    If[!KeyExistsQ[DiffExp`State`FEC, MatrixDirectory],
      DiffExp`Utilities`PrintWarning["No differential equations are loaded!"];
    ];
  ];

  DiffExp`State`IntegrationSequence = {};
  DiffExp`State`DEqnMatricesFactored = Association[{}];
  DiffExp`State`DEqnMatricesFactoredClosedForm = Association[{}];
  DiffExp`State`DEqnMatricesExpanded = Association[{}];
  DiffExp`State`AlphabetLogRulesFactored = Association[{}];
  DiffExp`State`AlphabetLogRulesExpanded = Association[{}];

  DiffExp`State`AnalyticContinuationReplacementsAssociation = Association[{}];

  DiffExp`Integration`UpdateIntReps[DiffExp`State`IMaxLogOrderDefault];

  #[[2]][assoc] & /@ DiffExp`State`DiffExpExtensions;

  DiffExp`Utilities`PrintInfo["Configuration updated."][1];

  CurrentConfiguration[]
];

End[];

EndPackage[];
