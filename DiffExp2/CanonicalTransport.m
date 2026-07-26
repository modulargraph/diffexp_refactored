(* DiffExp2/CanonicalTransport.m

   Transport for canonical dlog systems

       d f = eps Sum_i C_i dlog(W_i) f

   along a straight line in kinematic space.  This is deliberately separate
   from the exact rational-matrix solver: algebraic letters may contain square
   roots even though LoadSystem cannot accept an irrational line matrix.
   Epsilon triangularity is retained explicitly, including at a
   regular-singular starting point. *)

BeginPackage[
  "DiffExp2`CanonicalTransport`",
  {"DiffExp2`Config`", "DiffExp2`Tolerances`"}
];

LoadCanonicalSystem::usage =
  "LoadCanonicalSystem[spec] validates constant matrices, alphabet letters, and kinematic variables for a canonical dlog system.";
TransportCanonicalLine::usage =
  "TransportCanonicalLine[sys,boundary,from,to,opts] transports a finite epsilon-coefficient boundary through explicitly supplied canonical charts.";

Begin["`Private`"];

err[id_String, detail_, payload_:<||>] :=
  DiffExp2`Tolerances`DE2Error[id, Join[
    <|"Module" -> "CanonicalTransport", "Message" -> detail|>, payload]];

squareMatrixDimension[matrix_] := Module[{dimensions = Dimensions[matrix]},
  If[Length[dimensions] === 2 && dimensions[[1]] === dimensions[[2]],
    dimensions[[1]], Missing["NotSquare"]]];

LoadCanonicalSystem[spec_Association] := Module[
  {matrices, letters, variables, dimensions, letterSymbols,
    unknownSymbols},
  matrices = Lookup[spec, "ConstantMatrices", Missing["NotAvailable"]];
  letters = Lookup[spec, "Letters", Missing["NotAvailable"]];
  variables = Lookup[spec, "Variables", Missing["NotAvailable"]];
  If[!ListQ[matrices] || matrices === {},
    err["E6", "\"ConstantMatrices\" must be a nonempty list"]];
  If[!ListQ[letters] || Length[letters] =!= Length[matrices],
    err["E6",
      "\"Letters\" must contain one expression per constant matrix",
      <|"MatrixCount" -> Length[matrices],
        "LetterCount" -> If[ListQ[letters], Length[letters], Missing[]]|>]];
  If[!ListQ[variables] || variables === {} ||
      !AllTrue[variables, MatchQ[#, _Symbol] &] ||
      DuplicateFreeQ[variables] =!= True,
    err["E6",
      "\"Variables\" must be a nonempty duplicate-free list of symbols",
      <|"Variables" -> variables|>]];
  dimensions = squareMatrixDimension /@ matrices;
  If[MemberQ[dimensions, _Missing] ||
      Length[DeleteDuplicates[dimensions]] =!= 1,
    err["E6",
      "all canonical matrices must be square and have one common dimension",
      <|"Dimensions" -> (Dimensions /@ matrices)|>]];
  If[!AllTrue[matrices,
      ArrayQ[Normal[#], 2, NumericQ] &],
    err["E6",
      "canonical matrices must contain only constant numeric entries"]];
  If[!FreeQ[matrices, Alternatives @@ variables],
    err["E6",
      "canonical matrices must not depend on kinematic variables"]];
  If[AnyTrue[letters, TrueQ[PossibleZeroQ[#]] &],
    err["E6", "alphabet letters must not be identically zero"]];
  letterSymbols = DeleteDuplicates@Cases[
    letters,
    symbol_Symbol /; Context[symbol] =!= "System`",
    Infinity
  ];
  unknownSymbols = Select[
    letterSymbols, !MemberQ[variables, #] &];
  If[unknownSymbols =!= {},
    err["E6",
      "every non-System symbol in the alphabet must be declared as a variable",
      <|"Symbols" -> unknownSymbols|>]];
  <|
    "Schema" -> "DiffExp2.CanonicalSystem/v1",
    "ConstantMatrices" -> (SparseArray /@ matrices),
    "Letters" -> letters,
    "Variables" -> variables,
    "Dimension" -> First[dimensions],
    "LetterCount" -> Length[letters]
  |>
];

LoadCanonicalSystem[x___] :=
  err["E6", "LoadCanonicalSystem expects one Association",
    <|"Arguments" -> {x}|>];

Options[TransportCanonicalLine] = {
  "ExpansionOrder" -> Automatic,
  "EpsilonOrder" -> Automatic,
  "WorkingPrecision" -> Automatic,
  "UsePade" -> Automatic,
  "ChartCenters" -> Automatic,
  "ChartBoundaries" -> Automatic
};

pointAssociation[point_Association] := point;
pointAssociation[point : {___Rule}] := Association[point];
pointAssociation[point_] :=
  err["E6", "line endpoints must be Associations or lists of rules",
    <|"Point" -> point|>];

pointValues[variables_List, point_] := Module[{association, values},
  association = pointAssociation[point];
  values = Lookup[association, variables, Missing["NotAvailable"]];
  If[MemberQ[values, _Missing],
    err["E6", "line endpoint does not assign every canonical variable",
      <|"Variables" -> variables, "Point" -> point|>]];
  If[!AllTrue[values, realNumericQ],
    err["E6", "canonical line endpoint values must be real numeric values",
      <|"Values" -> values|>]];
  values
];

realNumericQ[value_] :=
  NumericQ[value] && TrueQ[PossibleZeroQ[Im[N[value, 50]]]];

resolvedInteger[value_, key_String, minimum_Integer] := Module[{answer},
  answer = Replace[value, Automatic :> DiffExp2`Config`CFG[key]];
  If[!IntegerQ[answer] || answer < minimum,
    err["E6", key <> " must be an integer at least " <> ToString[minimum],
      <|"Value" -> answer|>]];
  answer
];

resolvedBoolean[value_, key_String] := Module[{answer},
  answer = Replace[value, Automatic :> DiffExp2`Config`CFG[key]];
  If[!BooleanQ[answer],
    err["E6", key <> " must be True or False", <|"Value" -> answer|>]];
  answer
];

chartGeometry[centersOption_, boundariesOption_] := Module[
  {centers, boundaries},
  Which[
    centersOption === Automatic && boundariesOption === Automatic,
      err["E6",
        "canonical transport requires explicit chart centers and boundaries"],
    centersOption === Automatic && ListQ[boundariesOption],
      err["E6",
        "\"ChartCenters\" is required when \"ChartBoundaries\" is explicit"],
    ListQ[centersOption] && boundariesOption === Automatic,
      err["E6",
        "\"ChartBoundaries\" is required when \"ChartCenters\" is explicit"],
    ListQ[centersOption] && ListQ[boundariesOption],
      centers = centersOption; boundaries = boundariesOption,
    True,
      err["E6", "invalid canonical chart geometry"]
  ];
  If[Length[boundaries] =!= Length[centers] + 1 ||
      Length[centers] === 0,
    err["E6",
      "canonical charts require one more boundary than center",
      <|"Centers" -> centers, "Boundaries" -> boundaries|>]];
  If[!AllTrue[Join[centers, boundaries], realNumericQ],
    err["E6",
      "canonical chart centers and boundaries must be real numeric values"]];
  If[!TrueQ[PossibleZeroQ[boundaries[[1]]]] ||
      !TrueQ[PossibleZeroQ[boundaries[[-1]] - 1]] ||
      !And @@ Thread[Differences[N[boundaries, 50]] > 0],
    err["E6",
      "canonical chart boundaries must increase strictly from 0 to 1",
      <|"Boundaries" -> boundaries|>]];
  If[!And @@ MapThread[
      TrueQ[N[#2 - #1, 50] <= 0 && N[#1 - #3, 50] <= 0] &,
      {centers, Most[boundaries], Rest[boundaries]}],
    err["E6", "each canonical chart center must lie in its boundary interval",
      <|"Centers" -> centers, "Boundaries" -> boundaries|>]];
  {centers, boundaries}
];

finiteNumericQ[value_] :=
  NumericQ[value] &&
  FreeQ[value, Indeterminate | ComplexInfinity | DirectedInfinity[_]];

materialize[values_, workingPrecision_Integer] :=
  N[values, workingPrecision];

matrixCombination[values_List, matrices_List, dimension_Integer,
    workingPrecision_Integer] := Module[{combined},
  combined = Total[MapThread[#1 #2 &, {values, matrices}]];
  If[combined === 0,
    combined = ConstantArray[0, {dimension, dimension}],
    combined = Normal[combined]];
  materialize[combined, workingPrecision]
];

letterCoefficientTable[dlogs_List, center_, order_Integer,
    workingPrecision_Integer, lineParameter_Symbol] := Module[
  {z = Unique["canonicalLocal$"], series, coefficients},
  series = Map[
    Refine[
      Normal@Series[
        N[# /. lineParameter -> center + z, workingPrecision],
        {z, 0, order - 1}
      ],
      Element[z, Reals]
    ] &,
    dlogs
  ];
  coefficients = Table[
    Map[Coefficient[#, z, power] &, series],
    {power, 0, order - 1}
  ];
  coefficients = materialize[coefficients, workingPrecision];
  If[!ArrayQ[coefficients, 2, finiteNumericQ] ||
      !FreeQ[coefficients, Piecewise | ConditionalExpression | z],
    err["E7",
      "an algebraic letter branch remained unresolved in a regular chart",
      <|"Center" -> center|>]];
  coefficients
];

polynomialEvaluate[table_List, step_, order_Integer,
    workingPrecision_Integer] :=
  materialize[
    Total[
      MapIndexed[#1 step^(First[#2] - 1) &, table]
    ],
    workingPrecision
  ];

padeEvaluate[table_List, step_, order_Integer,
    workingPrecision_Integer, usePade_] := Module[
  {z = Unique["canonicalPade$"], numeratorDegree, denominatorDegree,
    fallback = 0, values},
  If[!TrueQ[usePade],
    Return[<|"Value" ->
      polynomialEvaluate[table, step, order, workingPrecision],
      "Fallbacks" -> 0|>]];
  numeratorDegree = Floor[order/2];
  denominatorDegree = Ceiling[order/2];
  values = Map[
    Function[coefficients,
      Module[{polynomial, approximant, value, fallbackValue,
        nonzeroPositions, effectiveDegree},
        fallbackValue = Sum[
          coefficients[[power + 1]] step^power,
          {power, 0, order}
        ];
        polynomial = Sum[
          coefficients[[power + 1]] z^power,
          {power, 0, order}
        ];
        nonzeroPositions = Flatten@Position[
          coefficients, coefficient_ /;
            !TrueQ[PossibleZeroQ[coefficient]]];
        effectiveDegree = If[
          nonzeroPositions === {}, -Infinity,
          Max[nonzeroPositions] - 1
        ];
        If[effectiveDegree <= numeratorDegree,
          Return[materialize[fallbackValue, workingPrecision], Module]
        ];
        approximant = Quiet[
          Check[
            PadeApproximant[
              polynomial,
              {z, 0, {numeratorDegree, denominatorDegree}}
            ],
            $Failed
          ]
        ];
        value = If[approximant === $Failed ||
            Head[Unevaluated[approximant]] === PadeApproximant,
          $Failed,
          Quiet[Check[approximant /. z -> step, $Failed]]
        ];
        If[value === $Failed || !finiteNumericQ[value],
          fallback++;
          value = fallbackValue
        ];
        materialize[value, workingPrecision]
      ]
    ],
    Transpose[table]
  ];
  <|"Value" -> values, "Fallbacks" -> fallback|>
];

regularChart[inputBlocks_List, left_, center_, right_, dlogs_List,
    matrices_List, dimension_Integer, epsilonOrder_Integer,
    expansionOrder_Integer, workingPrecision_Integer, usePade_,
    lineParameter_Symbol] := Module[
  {letterCoefficients, coefficientMatrices, coefficients, forcing,
    leftStep, rightStep, evaluated, leftTail, output, fallbackCount = 0,
    timing},
  timing = AbsoluteTiming[
    letterCoefficients = letterCoefficientTable[
      dlogs, center, expansionOrder, workingPrecision, lineParameter];
    coefficientMatrices = matrixCombination[
        #, matrices, dimension, workingPrecision] & /@
      letterCoefficients;
    coefficients = Table[
      ConstantArray[0, {expansionOrder + 1, dimension}],
      {epsilonIndex, 0, epsilonOrder}
    ];
    coefficients[[1, 1]] = inputBlocks[[1]];
    leftStep = left - center;
    rightStep = right - center;
    Do[
      Do[
        forcing = Total[
          Table[
            coefficientMatrices[[matrixPower + 1]] .
              coefficients[[
                epsilonIndex, power - matrixPower + 1]],
            {matrixPower, 0, power}
          ]
        ];
        coefficients[[epsilonIndex + 1, power + 2]] =
          materialize[forcing/(power + 1), workingPrecision],
        {power, 0, expansionOrder - 1}
      ];
      evaluated = padeEvaluate[
        coefficients[[epsilonIndex + 1]], leftStep,
        expansionOrder, workingPrecision, usePade];
      fallbackCount += evaluated["Fallbacks"];
      leftTail = evaluated["Value"];
      coefficients[[epsilonIndex + 1, 1]] =
        materialize[
          inputBlocks[[epsilonIndex + 1]] - leftTail,
          workingPrecision
        ],
      {epsilonIndex, 1, epsilonOrder}
    ];
    output = Table[
      evaluated = padeEvaluate[
        coefficients[[epsilonIndex + 1]], rightStep,
        expansionOrder, workingPrecision, usePade];
      fallbackCount += evaluated["Fallbacks"];
      evaluated["Value"],
      {epsilonIndex, 0, epsilonOrder}
    ];
  ][[1]];
  {
    output,
    <|"Center" -> center, "Domain" -> {left, right},
      "SingularStart" -> False, "TimingSeconds" -> timing,
      "PadeFallbacks" -> fallbackCount|>
  }
];

singularStartChart[inputBlocks_List, center_, right_, dlogs_List,
    matrices_List, dimension_Integer, epsilonOrder_Integer,
    expansionOrder_Integer, workingPrecision_Integer, usePade_,
    lineParameter_Symbol] := Module[
  {z = Unique["canonicalSingular$"], residues, residueMatrix,
    compatibility, regularExpressions, series, regularValues,
    coefficientMatrices, coefficients, forcing, output, evaluated,
    fallbackCount = 0, timing},
  timing = AbsoluteTiming[
    residues = Table[
      Quiet[
        Check[
          Limit[
            (lineParameter - center) dlogs[[index]],
            lineParameter -> center,
            Direction -> "FromAbove"
          ],
          $Failed
        ]
      ],
      {index, Length[dlogs]}
    ];
    If[MemberQ[residues, $Failed] ||
        !AllTrue[residues, NumericQ],
      err["E7", "could not determine the starting dlog residue",
        <|"Center" -> center|>]];
    residueMatrix = matrixCombination[
      residues, matrices, dimension, workingPrecision];
    compatibility = Table[
      residueMatrix . inputBlocks[[epsilonIndex]],
      {epsilonIndex, 1, epsilonOrder}
    ];
    If[AnyTrue[Flatten[compatibility], !TrueQ[PossibleZeroQ[#]] &],
      err["E7",
        "the supplied boundary is not a finite regular solution at the canonical singular point",
        <|"Center" -> center,
          "ResidueActions" ->
            (Max[Abs[#]] & /@ compatibility)|>]];
    regularExpressions = MapThread[
      Together[#1 - #2/(lineParameter - center)] &,
      {dlogs, residues}
    ];
    series = Map[
      Refine[
        Normal@Series[
          N[# /. lineParameter -> center + z, workingPrecision],
          {z, 0, expansionOrder - 1}
        ],
        z > 0
      ] &,
      regularExpressions
    ];
    regularValues = Table[
      Map[Coefficient[#, z, power] &, series],
      {power, 0, expansionOrder - 1}
    ];
    regularValues = materialize[regularValues, workingPrecision];
    If[!ArrayQ[regularValues, 2, finiteNumericQ] ||
        !FreeQ[regularValues, Piecewise | ConditionalExpression | z],
      err["E7",
        "an algebraic letter branch remained unresolved at the singular start",
        <|"Center" -> center|>]];
    coefficientMatrices = matrixCombination[
        #, matrices, dimension, workingPrecision] & /@
      regularValues;
    coefficients = Table[
      ConstantArray[0, {expansionOrder + 1, dimension}],
      {epsilonIndex, 0, epsilonOrder}
    ];
    Do[
      coefficients[[epsilonIndex + 1, 1]] =
        inputBlocks[[epsilonIndex + 1]],
      {epsilonIndex, 0, epsilonOrder}
    ];
    Do[
      Do[
        forcing =
          residueMatrix .
            coefficients[[epsilonIndex, power + 2]] +
          Total[
            Table[
              coefficientMatrices[[matrixPower + 1]] .
                coefficients[[
                  epsilonIndex, power - matrixPower + 1]],
              {matrixPower, 0, power}
            ]
          ];
        coefficients[[epsilonIndex + 1, power + 2]] =
          materialize[forcing/(power + 1), workingPrecision],
        {epsilonIndex, 1, epsilonOrder}
      ],
      {power, 0, expansionOrder - 1}
    ];
    output = Table[
      evaluated = padeEvaluate[
        coefficients[[epsilonIndex + 1]], right - center,
        expansionOrder, workingPrecision, usePade];
      fallbackCount += evaluated["Fallbacks"];
      evaluated["Value"],
      {epsilonIndex, 0, epsilonOrder}
    ];
  ][[1]];
  {
    output,
    <|"Center" -> center, "Domain" -> {center, right},
      "SingularStart" -> True, "Residues" -> residues,
      "ResidueRank" -> MatrixRank[residueMatrix],
      "TimingSeconds" -> timing, "PadeFallbacks" -> fallbackCount|>
  }
];

TransportCanonicalLine[system_Association, boundary_List, from_, to_,
    OptionsPattern[]] := Module[
  {dimension, matrices, letters, variables, epsilonOrder, expansionOrder,
    workingPrecision, usePade, centers, boundaries, fromValues, toValues,
    lineParameter = Unique["canonicalLine$"], lineRules, lineLetters, dlogs,
    blocks, chartResult, chartRecords = {}, totalTiming},
  If[Lookup[system, "Schema", None] =!= "DiffExp2.CanonicalSystem/v1",
    err["E6", "first argument is not a canonical-system record"]];
  dimension = system["Dimension"];
  matrices = system["ConstantMatrices"];
  letters = system["Letters"];
  variables = system["Variables"];
  If[Dimensions[boundary][[1]] =!= dimension ||
      Length[Dimensions[boundary]] =!= 2,
    err["E6",
      "boundary must be a rectangular dimension-by-epsilon array",
      <|"ExpectedDimension" -> dimension,
        "BoundaryDimensions" -> Dimensions[boundary]|>]];
  epsilonOrder = resolvedInteger[
    Replace[OptionValue["EpsilonOrder"],
      Automatic :> Dimensions[boundary][[2]] - 1],
    "EpsilonOrder", 0
  ];
  If[Dimensions[boundary][[2]] < epsilonOrder + 1,
    err["E6", "boundary does not contain the requested epsilon order",
      <|"BoundaryColumns" -> Dimensions[boundary][[2]],
        "EpsilonOrder" -> epsilonOrder|>]];
  expansionOrder = resolvedInteger[
    OptionValue["ExpansionOrder"], "ExpansionOrder",
    DiffExp2`Tolerances`$MinExpansionOrder];
  workingPrecision = resolvedInteger[
    OptionValue["WorkingPrecision"], "WorkingPrecision", 20];
  usePade = resolvedBoolean[OptionValue["UsePade"], "UsePade"];
  {centers, boundaries} = chartGeometry[
    OptionValue["ChartCenters"], OptionValue["ChartBoundaries"]];
  fromValues = pointValues[variables, from];
  toValues = pointValues[variables, to];
  lineRules = Thread[
    variables ->
      (fromValues + lineParameter (toValues - fromValues))
  ];
  lineLetters = Together[letters /. lineRules];
  dlogs = Map[Together[D[Log[#], lineParameter]] &, lineLetters];
  blocks = materialize[
    Transpose[boundary[[All, 1 ;; epsilonOrder + 1]]],
    workingPrecision
  ];
  totalTiming = AbsoluteTiming[
    Do[
      chartResult = If[chartIndex === 1 &&
          TrueQ[PossibleZeroQ[centers[[1]] - boundaries[[1]]]],
        singularStartChart[
          blocks, centers[[1]], boundaries[[2]], dlogs, matrices,
          dimension, epsilonOrder, expansionOrder, workingPrecision,
          usePade, lineParameter
        ],
        regularChart[
          blocks, boundaries[[chartIndex]], centers[[chartIndex]],
          boundaries[[chartIndex + 1]], dlogs, matrices, dimension,
          epsilonOrder, expansionOrder, workingPrecision, usePade,
          lineParameter
        ]
      ];
      blocks = chartResult[[1]];
      AppendTo[chartRecords,
        Join[<|"Index" -> chartIndex|>, chartResult[[2]]]],
      {chartIndex, Length[centers]}
    ];
  ][[1]];
  <|
    "Schema" -> "DiffExp2.CanonicalTransportResult/v1",
    "Value" -> Transpose[blocks],
    "From" -> from,
    "To" -> to,
    "Dimension" -> dimension,
    "EpsilonOrder" -> epsilonOrder,
    "ExpansionOrder" -> expansionOrder,
    "WorkingPrecision" -> workingPrecision,
    "UsePade" -> usePade,
    "ChartCenters" -> centers,
    "ChartBoundaries" -> boundaries,
    "Charts" -> chartRecords,
    "TimingSeconds" -> totalTiming,
    "PadeFallbacks" -> Total[Lookup[chartRecords, "PadeFallbacks", 0]]
  |>
];

TransportCanonicalLine[x___] :=
  err["E6",
    "TransportCanonicalLine expects a canonical system, boundary, and two endpoints",
    <|"Arguments" -> {x}|>];

End[];
EndPackage[];
