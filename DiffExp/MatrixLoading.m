(* ::Package:: *)

(* DiffExp MatrixLoading Subpackage *)
(* This package handles matrix loading and preparation *)

BeginPackage["DiffExp`MatrixLoading`", {"DiffExp`Symbols`", "DiffExp`State`", "DiffExp`Utilities`", "DiffExp`SeriesOps`"}];

(* Matrix loading functions *)
LoadMatrices::usage = "LoadMatrices[Folder_] loads matrices from the specified folder.";
PrepareMatrices::usage = "PrepareMatrices[line_] prepares matrices for a line.";
PrepareMatricesFrom1::usage = "PrepareMatricesFrom1[lineorig_,linenew_] prepares matrices reusing previous results.";
PrepareMatricesFrom::usage = "PrepareMatricesFrom[lineorig_,linenew_] prepares and expands matrices.";
PrepareMatricesFactored::usage = "PrepareMatricesFactored[line_] prepares factored matrices.";
PrepareMatricesExpanded::usage = "PrepareMatricesExpanded[line_] prepares expanded matrices.";
ClearMatrices::usage = "ClearMatrices[line_] clears matrix caches.";
InitializeIntegrationSequence::usage = "InitializeIntegrationSequence[line_] detects integration sequence.";

Begin["`Private`"];

(* Loads files in the location of MatrixDirectory *)
LoadMatrices[Folder_] := Module[
  {matrixFiles, CurrFileName, CurrMatrix, dimensionCheck, MatrixHeads, unsupportedHeads, primeFactorList, NewSquareRootSingularities,
   parsedVariables, SqrtFlips, VarsPartialDerivatives, ExtraVars, AlphabetLogs, FilePattern, ExtraVarsEncountered},

  DiffExp`Utilities`PrintInfo["Loading matrices."][1];

  VarsPartialDerivatives = DiffExp`State`FEC[System`Variables];
  If[DiffExp`State`FEC[System`Variables] === {},
    matrixFiles = FileNames[FileNameJoin[{DiffExp`State`FEC[MatrixDirectory], "d*_*.m"}]];

    If[Length[matrixFiles] === 0,
      DiffExp`Utilities`ReportError["No partial derivative matrices found in the given directory."];
    ];

    FilePattern = Longest[pre__] ~~ $PathnameSeparator ~~ "d" ~~ Shortest[a___] ~~ "_" ~~ Shortest[b__] ~~ ".m" /;
      NumericQ[ToExpression[b]] || b === "d";
    DiffExp`Utilities`PrintInfo["Found files: ",
      Flatten[StringCases[#, FilePattern :> "d" <> a <> "_" <> b <> ".m"] & @matrixFiles]][1];
    parsedVariables = Flatten[StringCases[#, FilePattern :> {a, b}], 1] & /@ matrixFiles;
    parsedVariables = DeleteCases[parsedVariables, {}];

    If[DeleteDuplicates[parsedVariables[[All, 2]]] === {"d"},
      DiffExp`Utilities`PrintInfo["Found differential equations in closed form."][1];
      DiffExp`State`UsingClosedFormMatrix = True;,
      DiffExp`State`UsingClosedFormMatrix = False;
    ];

    parsedVariables = ToExpression /@ (parsedVariables[[All, 1]] // DeleteDuplicates);
    parsedVariables = DeleteCases[parsedVariables, Null];
    VarsPartialDerivatives = parsedVariables;

    If[MemberQ[parsedVariables, DiffExp`State`LineParameterVal],
      DiffExp`Utilities`ReportError["Some of the kinematic invariants or masses are named the same as the line parameter."];
    ];
  ];

  DiffExp`State`ExpansionMatricesClosedForm = Association[];

  (* Handling of special canonical matrix file *)
  If[(FileExistsQ[FileNameJoin[{DiffExp`State`FEC[MatrixDirectory], "d_1.m"}]]) &&
      !DiffExp`State`UsingClosedFormMatrix,
    DiffExp`State`ExpansionMatricesCanonical1 = Import[FileNameJoin[{DiffExp`State`FEC[MatrixDirectory], "d_1.m"}]] // Expand;

    (* Checking consistency of d_1.m *)
    If[!(And @@ Flatten[MapAt[
        (And @@ (MatchQ[#, (a_: 1) Log[b_] /; NumericQ[a]] & /@ DeleteCases[DiffExp`Utilities`SplitSum[#], 0])) &,
        DiffExp`State`ExpansionMatricesCanonical1, {All, All}]]) === True,
      DiffExp`Utilities`ReportError["The matrix d_1.m is not of the right form."];
    ];

    AlphabetLogs = DiffExp`Utilities`GetCases[DiffExp`State`ExpansionMatricesCanonical1, Log[a_]];

    ExtraVars = AlphabetLogs /. Log[a_] :> a // System`Variables;
    VarsPartialDerivatives = Union[VarsPartialDerivatives, ExtraVars];

    DiffExp`State`AlphabetLogRules = Table[Log[Subscript[l, ind]] -> AlphabetLogs[[ind]], {ind, Length[AlphabetLogs]}];

    DiffExp`State`ExpansionMatricesCanonical1 = DiffExp`State`ExpansionMatricesCanonical1 /. (Reverse /@ DiffExp`State`AlphabetLogRules);

    DiffExp`Utilities`PrintInfo["Loaded canonical matrix."][1];,

    DiffExp`State`ExpansionMatricesCanonical1 = "ZeroM";
    AlphabetLogs = {};
    DiffExp`State`AlphabetLogRules = {};
  ];

  DiffExp`Utilities`PrintInfo["Kinematic invariants and masses: ", VarsPartialDerivatives][1];
  DiffExp`State`DiffExpConfiguration[System`Variables] = VarsPartialDerivatives;

  If[!DiffExp`State`UsingClosedFormMatrix,
    DiffExp`State`ExpansionMatrices = Table[
      CurrFileName = FileNameJoin[{DiffExp`State`MatrixDirectoryVal, "d" <> ToString[var] <> "_" <> ToString[ord] <> ".m"}];

      If[FileExistsQ[CurrFileName],
        CurrMatrix = Import[CurrFileName];

        If[!(ExtraVarsEncountered = Complement[System`Variables[CurrMatrix], DiffExp`State`ExternalScalesVal]) === {},
          DiffExp`Utilities`ReportError["The file ", "d" <> ToString[var] <> "_" <> ToString[ord] <> ".m",
            " contains the variables: ", ExtraVarsEncountered, ", but no partial derivatives are given for these."];
        ];

        {var, ord} -> CurrMatrix,

        If[MemberQ[VarsPartialDerivatives, var] && AlphabetLogs === {},
          DiffExp`Utilities`PrintInfo["Assuming M[", var, "][", ord, "] is zero."][1];
        ];

        {var, ord} -> "ZeroM"
      ],
      {var, DiffExp`State`ExternalScalesVal}, {ord, Range[0, DiffExp`State`EpsilonOrderVal]}
    ] // Flatten[#, 1] & // Association;,

    DiffExp`State`ExpansionMatrices = Association[];
    DiffExp`State`ExpansionMatricesClosedForm = Table[
      CurrFileName = FileNameJoin[{DiffExp`State`MatrixDirectoryVal, "d" <> ToString[var] <> "_d.m"}];

      CurrMatrix = Import[CurrFileName];
      If[!(ExtraVarsEncountered = Complement[System`Variables[CurrMatrix], Join[DiffExp`State`ExternalScalesVal, {DiffExp`Symbols`eps, DiffExp`Symbols`\[Epsilon]}]]) === {},
        DiffExp`Utilities`ReportError["The file contains variables without partial derivatives."];
      ];

      var -> CurrMatrix,
      {var, DiffExp`State`ExternalScalesVal}
    ];
  ];

  (* Sanity check on matrix dimensions *)
  dimensionCheck = (If[!# === "ZeroM", # // Dimensions, Null] & /@
    Join[Values[DiffExp`State`ExpansionMatrices], {DiffExp`State`ExpansionMatricesCanonical1},
      Values[DiffExp`State`ExpansionMatricesClosedForm]]) // DeleteCases[#, Null] &;

  If[!(dimensionCheck // SameQ),
    DiffExp`Utilities`ReportError["Loaded matrices are of different dimensions: ", dimensionCheck];
  ];

  DiffExp`Utilities`PrintInfo["Loaded system of size ",
    ToString[#[[1]]] <> " x " <> ToString[#[[2]]] &@First@
      DeleteCases[Dimensions /@ Join[Values[DiffExp`State`ExpansionMatrices],
        {DiffExp`State`ExpansionMatricesCanonical1}, Values[DiffExp`State`ExpansionMatricesClosedForm]], {}]][1];

  DiffExp`State`ExpansionMatrices = DiffExp`State`ExpansionMatrices /.
    "ZeroM" -> DiffExp`Utilities`CA[0, {dimensionCheck[[1, 1]], dimensionCheck[[1, 1]]}];
  DiffExp`State`ExpansionMatricesCanonical1 = DiffExp`State`ExpansionMatricesCanonical1 /.
    "ZeroM" -> DiffExp`Utilities`CA[0, {dimensionCheck[[1, 1]], dimensionCheck[[1, 1]]}];

  DiffExp`State`NumIntegrals = dimensionCheck[[1, 1]];

  MatrixHeads = DiffExp`Utilities`GetCases[
    Join[DiffExp`State`ExpansionMatrices // Values, DiffExp`State`ExpansionMatricesClosedForm // Values,
      AlphabetLogs /. Log[a_] :> a] // Flatten,
    a_ :> a[[0]]
  ];
  unsupportedHeads = Complement[MatrixHeads, {Association, List, Complex, Integer, Plus, Power, Rational, Symbol, Times}];

  If[Length[unsupportedHeads] > 0,
    DiffExp`Utilities`ReportError["The differential equation matrices contain unsupported functions: ", unsupportedHeads];
  ];

  If[Length[DiffExp`Utilities`GetCases[
      Join[DiffExp`State`ExpansionMatrices // Values, DiffExp`State`ExpansionMatricesClosedForm // Values,
        AlphabetLogs /. Log[a_] :> a],
      Power[a_, b_] /; Denominator[b] > 2]] > 0,
    DiffExp`Utilities`ReportError["Differential equations contain higher order roots."];
  ];

  (* Detecting square roots in the differential equations *)
  DiffExp`State`DEqnSquareRoots = DeleteCases[
    DiffExp`Utilities`GetCases[{DiffExp`State`ExpansionMatrices // Values, DiffExp`State`AlphabetLogRules,
      DiffExp`State`ExpansionMatricesClosedForm // Values},
      Power[a_, b_] /; Denominator[b] === 2 :> Expand[a]],
    !DiffExp`Utilities`DependsQ[#, Alternatives @@ DiffExp`State`ExternalScalesVal] &
  ];
  DiffExp`State`DEqnSquareRoots = DeleteDuplicates[DiffExp`State`DEqnSquareRoots,
    Expand[#1] === Expand[#2] || Expand[#1] === Expand[-#2] &];

  primeFactorList = FactorList /@ DiffExp`State`DEqnSquareRoots;
  If[DiffExp`Utilities`DependsQ[Length[#] > 2 & /@ primeFactorList, True],
    DiffExp`Utilities`ReportError["Matrices contain square roots which are not irreducible!"];
  ];

  NewSquareRootSingularities = DiffExp`State`SquareRootPrescriptionsAdded[];

  If[Length[NewSquareRootSingularities] > 0,
    DiffExp`Utilities`PrintInfo["Additional square roots encountered in the partial derivative matrices:"][1];
    DiffExp`Utilities`PrintInfo[NewSquareRootSingularities[[All, 1]] // TableForm][1];
    DiffExp`Utilities`PrintInfo["Assigning these roots +i\[Delta]."][1];
  ];

  DiffExp`State`DiffExpConfiguration[DeltaPrescriptions] =
    Union[DiffExp`State`FEC[DeltaPrescriptions], DiffExp`State`SquareRootPrescriptionsAdded[]];

  (* Flipping square roots which have -I\[Delta] prescription *)
  SqrtFlips = Table[
    If[\[Sigma][[2]] === 1,
      Sqrt[-\[Sigma][[1]] // Expand] -> -I Sqrt[\[Sigma][[1]]],
      Sqrt[\[Sigma][[1]]] -> -I Sqrt[-\[Sigma][[1]] // Expand]
    ],
    {\[Sigma], DiffExp`State`DiffExpConfiguration[DeltaPrescriptions]}
  ];

  If[Length[DiffExp`State`DEqnSquareRoots] > 0,
    DiffExp`Utilities`PrintDebug["Flipping roots in matrices according to i\[Delta]-prescriptions given: ",
      Select[{Sqrt[#] -> (Sqrt[#] /. SqrtFlips), Sqrt[#] - (Sqrt[#] /. SqrtFlips) === 0} & /@
        DiffExp`State`DEqnSquareRoots, #[[2]] === False &][[All, 1]]][1];
  ];

  DiffExp`State`ExpansionMatrices = Association[Normal[DiffExp`State`ExpansionMatrices] /.
    a_^b_ /; Denominator[b] == 2 :> a^(b - 1/2) (Sqrt[a // Expand] /. SqrtFlips)];
  DiffExp`State`ExpansionMatricesClosedForm = Association[Normal[DiffExp`State`ExpansionMatricesClosedForm] /.
    a_^b_ /; Denominator[b] == 2 :> a^(b - 1/2) (Sqrt[a // Expand] /. SqrtFlips)];
  DiffExp`State`AlphabetLogRules = DiffExp`State`AlphabetLogRules /.
    a_^b_ /; Denominator[b] == 2 :> a^(b - 1/2) (Sqrt[a // Expand] /. SqrtFlips);

  (* Loading all the factors which may yield singularities *)
  DiffExp`Utilities`PrintInfo["Getting irreducible factors.."][1];

  DiffExp`State`MatricesIrreducibleFactors = Replace[
    Join[
      (FactorList /@ ({
        DiffExp`State`ExpansionMatrices // Values // Flatten //
          DiffExp`Utilities`GetCases[#, a_^b_ /; IntegerQ[b] && b < 0 :> a] &,
        DiffExp`State`ExpansionMatricesClosedForm // Values // Flatten //
          DiffExp`Utilities`GetCases[#, a_^b_ /; IntegerQ[b] && b < 0 :> a] &,
        AlphabetLogs /. Log[a_] :> Factor[a]
      } // Flatten)) // Flatten,
      Cases[{DiffExp`State`ExpansionMatrices, DiffExp`State`ExpansionMatricesClosedForm, AlphabetLogs},
        a_^b_ /; Denominator[b] == 2 :> a, Infinity]
    ] /. DiffExp`Symbols`\[Epsilon] -> 0 // DeleteDuplicates //
    DeleteCases[#, a_ /; !DiffExp`Utilities`DependsQ[a, Alternatives @@ DiffExp`State`FEC[System`Variables]]] &,
    a_^b_ /; Denominator[b] == 2 :> a, 1
  ] // DeleteDuplicates[#, Expand[#1] === Expand[#2] || Expand[#1] === Expand[-#2] &] &;
];

(* Clear matrix caches *)
ClearMatrices[line_] := Module[{},
  If[!DiffExp`State`FEC["KeepMatrixExpansions"] === True,
    KeyDropFrom[DiffExp`State`DEqnMatricesFactored, line];
    KeyDropFrom[DiffExp`State`DEqnMatricesFactoredClosedForm, line];
    KeyDropFrom[DiffExp`State`DEqnMatricesExpanded, line];
  ];
];

ClearMatrices[] := Module[{},
  If[!DiffExp`State`FEC["KeepMatrixExpansions"] === True,
    DiffExp`State`DEqnMatricesFactored = Association[];
    DiffExp`State`DEqnMatricesFactoredClosedForm = Association[];
    DiffExp`State`DEqnMatricesExpanded = Association[];
  ];
];

(* Prepare matrices along a line *)
PrepareMatrices[line_Association] := (
  If[!KeyExistsQ[DiffExp`State`DEqnMatricesExpanded, line],
    PrepareMatricesFactored[line];
    PrepareMatricesExpanded[line];
  ];
);

(* Prepare matrices by reusing previously factored matrices *)
PrepareMatricesFrom1[lineorig_Association, linenew_Association] := Module[{ParRelns},
  If[(!KeyExistsQ[DiffExp`State`DEqnMatricesFactored, linenew]) &&
     (!KeyExistsQ[DiffExp`State`DEqnMatricesExpanded, linenew]),
    DiffExp`Utilities`PrintInfo["Preparing differential equations along current line segment."][2];
  ];

  ParRelns = DiffExp`LineSegmentation`RelateLines[lineorig, linenew];

  If[!KeyExistsQ[DiffExp`State`DEqnMatricesFactored, lineorig],
    DiffExp`Utilities`PrintDebug["Warning!: PrepareMatricesFrom is asked to re-use results that don't exist."][1];
    PrepareMatricesFactored[linenew];,

    If[!KeyExistsQ[DiffExp`State`DEqnMatricesFactored, linenew],
      DiffExp`State`DEqnMatricesFactored[linenew] = Association[
        Table[
          epsord -> (DiffExp`Utilities`PChop[
            If[DiffExp`State`FEC[UseMobius] === True, Identity, DiffExp`Utilities`FactorOrTogether[lineorig, linenew]][
              (DiffExp`State`DEqnMatricesFactored[lineorig][epsord] /. DiffExp`Symbols`x -> ParRelns) D[ParRelns, DiffExp`Symbols`x]]]),
          {epsord, 0, DiffExp`State`EpsilonOrderVal}
        ]
      ];

      If[!DiffExp`State`AlphabetLogRules === {},
        DiffExp`State`AlphabetLogRulesFactored[linenew] = MapAt[
          (# /. Log[a_] :> Log[a /. linenew // DiffExp`Utilities`FactorOrTogether[linenew] // DiffExp`Utilities`PChop]) &,
          DiffExp`State`AlphabetLogRules, {All, 2}
        ];
      ];
    ];
  ];
];

PrepareMatricesFrom[lineorig_Association, linenew_Association] := (
  PrepareMatricesFrom1[lineorig, linenew];
  DiffExp`Utilities`PrintInfo["Expanding partial derivative matrices."][2];
  PrepareMatricesExpanded[linenew];
);

(* Prepare factored matrices *)
PrepareMatricesFactored[line_Association] := (
  If[!KeyExistsQ[DiffExp`State`DEqnMatricesFactored, line],

    If[!DiffExp`State`UsingClosedFormMatrix,

      DiffExp`State`DEqnMatricesFactored[line] = Association[
        Table[
          epsord -> (DiffExp`Utilities`PChop[DiffExp`Utilities`FactorOrTogether[line][
            Sum[(DiffExp`State`ExpansionMatrices[{v, epsord}] /. line) D[line[v], DiffExp`Symbols`x],
              {v, DiffExp`State`ExternalScalesVal}]]]),
          {epsord, 0, DiffExp`State`EpsilonOrderVal}
        ]
      ];

      If[!DiffExp`State`AlphabetLogRules === {},
        DiffExp`State`AlphabetLogRulesFactored[line] = MapAt[
          (# /. Log[a_] :> Log[a /. line // DiffExp`Utilities`FactorOrTogether[line] // DiffExp`Utilities`PChop]) &,
          DiffExp`State`AlphabetLogRules, {All, 2}
        ];
      ];,

      DiffExp`State`DEqnMatricesFactoredClosedForm[line] = (DiffExp`Utilities`PChop[DiffExp`Utilities`FactorOrTogether[line][
        Sum[(DiffExp`State`ExpansionMatricesClosedForm[v] /. line) D[line[v], DiffExp`Symbols`x],
          {v, DiffExp`State`ExternalScalesVal}]]]);
      DiffExp`State`DEqnMatricesFactored[line] = Series[DiffExp`State`DEqnMatricesFactoredClosedForm[line],
        {DiffExp`Symbols`\[Epsilon], 0, DiffExp`State`EpsilonOrderVal}];
      DiffExp`State`DEqnMatricesFactored[line] = Association[
        Table[
          epsord -> DiffExp`Utilities`PChop[DiffExp`Utilities`FactorOrTogether[line][
            SeriesCoefficient[DiffExp`State`DEqnMatricesFactored[line], {DiffExp`Symbols`\[Epsilon], 0, epsord}]]],
          {epsord, 0, DiffExp`State`EpsilonOrderVal}
        ]
      ];
    ];
  ];
);

(* Prepare expanded matrices *)
PrepareMatricesExpanded[line_Association] := (
  If[!KeyExistsQ[DiffExp`State`DEqnMatricesExpanded, line],
    DiffExp`State`DEqnMatricesExpanded[line] = Series[# // N[#, DiffExp`State`FEWorkingPrecision] &,
      {DiffExp`Symbols`x, 0, DiffExp`State`ExpansionOrderVal}, Assumptions -> DiffExp`Symbols`x > 0] & /@
      DiffExp`State`DEqnMatricesFactored[line];

    If[!DiffExp`State`AlphabetLogRules === {},
      DiffExp`State`AlphabetLogRulesExpanded[line] = MapAt[
        D[Series[# // N[#, DiffExp`State`FEWorkingPrecision] &,
          {DiffExp`Symbols`x, 0, DiffExp`State`ExpansionOrderVal}, Assumptions -> DiffExp`Symbols`x > 0], DiffExp`Symbols`x] &,
        DiffExp`State`AlphabetLogRulesFactored[line], {All, 2}
      ];

      DiffExp`State`DEqnMatricesExpanded[line][1] = DiffExp`State`DEqnMatricesExpanded[line][1] +
        (DiffExp`State`ExpansionMatricesCanonical1 /. DiffExp`State`AlphabetLogRulesExpanded[line]) //
        DiffExp`SeriesOps`SExpand;
    ];
  ];
);

(* Detects the integration sequence *)
InitializeIntegrationSequence[line_] := Module[
  {HomogeneousMask, dependencyEdges, CurrIndex, IntegrationDependencies, redundantKeys},

  DiffExp`Utilities`PrintDebug["Analyzing integration sequence on current line."][1];
  HomogeneousMask = MapAt[DiffExp`Utilities`ZeroQ, DiffExp`State`DEqnMatricesFactored[line][0], {All, All}];

  (* A graph describing which integrals couple together in the differential equations *)
  dependencyEdges = Flatten[MapIndexed[
    (CurrIndex = #2[[1]]; CurrIndex -> # &) /@ #1 &,
    (Flatten[Position[#, False]]) & /@ HomogeneousMask
  ]];
  dependencyEdges = Join[dependencyEdges, Table[DirectedEdge[i, i], {i, DiffExp`State`NumIntegrals}]] // DeleteDuplicates;

  IntegrationDependencies = (# -> VertexOutComponent[Graph[dependencyEdges], {#}]) & /@ Range[DiffExp`State`NumIntegrals];

  (* Integration sequence and coupled integrals *)
  DiffExp`State`IntegrationSequence = SortBy[IntegrationDependencies, Length[#[[2]]] &][[All, 1]];
  DiffExp`State`IntegrationSequence = Association[
    # -> Union[DiffExp`Utilities`GetCases[ConnectedComponents[dependencyEdges // Graph, #], _Integer], {#}] & /@
      DiffExp`State`IntegrationSequence
  ];
  redundantKeys = (DeleteDuplicates@*Flatten)[Delete[#, -1] & /@ (DiffExp`State`IntegrationSequence // Values)];
  DiffExp`State`IntegrationSequence = KeyDrop[DiffExp`State`IntegrationSequence, redundantKeys] // Values;

  DiffExp`Utilities`PrintDebug["Integration sequence is ", DiffExp`State`IntegrationSequence][1];
  DiffExp`State`MaxCouplingOrder = (Length /@ DiffExp`State`IntegrationSequence) // Max;
  DiffExp`Utilities`PrintDebug["Maximum coupling order is ", DiffExp`State`MaxCouplingOrder][1];
];

End[];

EndPackage[];
