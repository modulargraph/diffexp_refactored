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
CanonicalLineChartGeometry::usage =
  "CanonicalLineChartGeometry[sys,from,to,opts] constructs clearance-certified charts from the finite algebraic singularities of the active canonical alphabet.";
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
      ArrayQ[#, 2, NumericQ] &],
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
  "ChartBoundaries" -> Automatic,
  "ImaginaryDetour" -> None,
  "BranchTrackingSubsteps" -> 32
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

resolvedDetour[None | Automatic, variableCount_Integer] :=
  <|"Amplitude" -> 0, "Directions" -> ConstantArray[0, variableCount]|>;
resolvedDetour[amplitude_?realNumericQ, variableCount_Integer] /;
    TrueQ[amplitude > 0] :=
  <|"Amplitude" -> amplitude,
    "Directions" -> ConstantArray[1, variableCount]|>;
resolvedDetour[spec_Association, variableCount_Integer] := Module[
  {amplitude, directions},
  amplitude = Lookup[spec, "Amplitude", Missing["NotAvailable"]];
  directions = Lookup[
    spec, "Directions", ConstantArray[1, variableCount]];
  If[!realNumericQ[amplitude] || !TrueQ[amplitude > 0],
    err["E6",
      "imaginary-detour amplitude must be a positive real number",
      <|"Amplitude" -> amplitude|>]];
  If[!ListQ[directions] ||
      Length[directions] =!= variableCount ||
      !AllTrue[directions, realNumericQ],
    err["E6",
      "imaginary-detour directions must contain one real number per variable",
      <|"Directions" -> directions,
        "VariableCount" -> variableCount|>]];
  <|"Amplitude" -> amplitude, "Directions" -> directions|>
];
resolvedDetour[value_, variableCount_Integer] :=
  err["E6",
    "\"ImaginaryDetour\" must be None, a positive amplitude, or an Association",
    <|"Value" -> value, "VariableCount" -> variableCount|>];

zeroCanonicalMatrixQ[matrix_] := If[
  Head[matrix] === SparseArray,
  Length[ArrayRules[matrix]] === 1,
  AllTrue[Flatten[matrix], TrueQ[PossibleZeroQ[#]] &]
];

algebraicNorm[expression_] := Fold[
  Function[{current, root},
    Together[current (current /. root -> -root)]
  ],
  expression,
  DeleteDuplicates@Cases[
    expression, Power[_, Rational[1, 2]], Infinity]
];

canonicalSingularityPolynomials[lineLetters_List,
    lineParameter_Symbol] := Module[
  {radicands, norms, expressions, polynomials},
  radicands = DeleteDuplicates@Cases[
    lineLetters,
    Power[base_, Rational[1, 2]] :> base,
    Infinity
  ];
  norms = algebraicNorm /@ lineLetters;
  If[!FreeQ[norms, Power[_, Rational[1, 2]]],
    err["E7",
      "could not eliminate all square roots while deriving canonical chart clearance"]];
  expressions = Join[
    radicands,
    Flatten[
      {Numerator[#], Denominator[#]} & /@ (Together /@ norms)
    ]
  ];
  polynomials = DeleteCases[
    DeleteDuplicates[expressions],
    value_ /; FreeQ[value, lineParameter]
  ];
  If[!AllTrue[polynomials, PolynomialQ[#, lineParameter] &],
    err["E7",
      "canonical chart clearance currently requires algebraic letters with polynomial line singularities"]];
  polynomials
];

canonicalSingularities[polynomials_List, lineParameter_Symbol,
    workingPrecision_Integer] := Module[
  {rootLists, roots, mergeTolerance},
  rootLists = Table[
    Quiet[
      Check[
        lineParameter /. NSolve[
          polynomials[[index]] == 0,
          lineParameter,
          WorkingPrecision -> workingPrecision
        ],
        $Failed
      ]
    ],
    {index, Length[polynomials]}
  ];
  If[MemberQ[rootLists, $Failed] ||
      !AllTrue[Flatten[rootLists], finiteNumericQ],
    err["E7",
      "could not determine every finite canonical line singularity"]];
  roots = SortBy[Flatten[rootLists], {Re, Im}];
  mergeTolerance = 10^-Floor[workingPrecision/2];
  DeleteDuplicates[
    roots,
    Abs[#1 - #2] < mergeTolerance &
  ]
];

clearanceCharts[singularities_List, safety_, workingPrecision_Integer,
    searchSubdivisions_Integer] := Module[
  {clearance, nextChart, centers = {}, boundaries = {0},
    left, chart, ratios, realTolerance},
  If[singularities === {},
    Return[<|
      "Centers" -> {1/2},
      "Boundaries" -> {0, 1},
      "MaximumClearanceRatio" -> 0
    |>]
  ];
  realTolerance = 10^-Floor[workingPrecision/3];
  If[AnyTrue[
      singularities,
      -realTolerance <= Re[#] <= 1 + realTolerance &&
        Abs[Im[#]] <= realTolerance &
    ],
    err["E7",
      "the proposed canonical contour contains a real singularity; choose an imaginary detour whose directions move every active letter singularity off the path",
      <|"Singularities" -> Select[
        singularities,
        -realTolerance <= Re[#] <= 1 + realTolerance &&
          Abs[Im[#]] <= realTolerance &
      ]|>]];
  clearance[point_?NumericQ] :=
    Min[Abs[N[point, workingPrecision] - singularities]];
  nextChart[leftPoint_?NumericQ] := Module[
    {midpoint, function, samples, bracket, lo, hi, mid,
      center, right},
    midpoint = (leftPoint + 1)/2;
    If[midpoint - leftPoint <= safety clearance[midpoint],
      Return[{midpoint, 1}]
    ];
    function[point_?NumericQ] :=
      point - leftPoint - safety clearance[point];
    samples = N[
      Subdivide[leftPoint, 1, searchSubdivisions],
      workingPrecision
    ];
    bracket = SelectFirst[
      Partition[samples, 2, 1],
      function[#[[1]]] <= 0 && function[#[[2]]] >= 0 &,
      Missing["NotFound"]
    ];
    If[MissingQ[bracket],
      err["E7",
        "could not bracket the next clearance-certified canonical chart",
        <|"LeftBoundary" -> leftPoint|>]];
    {lo, hi} = bracket;
    Do[
      mid = (lo + hi)/2;
      If[function[mid] >= 0, hi = mid, lo = mid],
      {Ceiling[Log2[10] workingPrecision]}
    ];
    center = (lo + hi)/2;
    right = Min[1, center + safety clearance[center]];
    {center, right}
  ];
  While[Last[boundaries] < 1,
    left = Last[boundaries];
    chart = nextChart[left];
    If[chart[[2]] <= left,
      err["E7",
        "clearance-certified canonical chart construction stopped making progress",
        <|"LeftBoundary" -> left, "Chart" -> chart|>]];
    AppendTo[centers, chart[[1]]];
    AppendTo[boundaries, chart[[2]]];
  ];
  ratios = MapThread[
    Max[#1 - #2[[1]], #2[[2]] - #1]/clearance[#1] &,
    {centers, Partition[boundaries, 2, 1]}
  ];
  <|
    "Centers" -> centers,
    "Boundaries" -> boundaries,
    "MaximumClearanceRatio" -> Max[ratios]
  |>
];

Options[CanonicalLineChartGeometry] = {
  "WorkingPrecision" -> Automatic,
  "ImaginaryDetour" -> None,
  "ClearanceFactor" -> 1/3,
  "SearchSubdivisions" -> 2000
};

CanonicalLineChartGeometry[system_Association, from_, to_,
    OptionsPattern[]] := Module[
  {workingPrecision, safety, searchSubdivisions, variables,
    fromValues, toValues, detour,
    lineParameter = Unique["canonicalGeometry$"], lineRules,
    activeIndices, activeLetters, lineLetters, polynomials,
    singularities, charts},
  If[Lookup[system, "Schema", None] =!= "DiffExp2.CanonicalSystem/v1",
    err["E6", "first argument is not a canonical-system record"]];
  workingPrecision = resolvedInteger[
    OptionValue["WorkingPrecision"], "WorkingPrecision", 20];
  safety = OptionValue["ClearanceFactor"];
  If[!realNumericQ[safety] || !TrueQ[0 < safety < 1],
    err["E6", "\"ClearanceFactor\" must be a real number between 0 and 1",
      <|"Value" -> safety|>]];
  searchSubdivisions = OptionValue["SearchSubdivisions"];
  If[!IntegerQ[searchSubdivisions] || searchSubdivisions < 10,
    err["E6", "\"SearchSubdivisions\" must be an integer at least 10",
      <|"Value" -> searchSubdivisions|>]];
  variables = system["Variables"];
  fromValues = pointValues[variables, from];
  toValues = pointValues[variables, to];
  detour = resolvedDetour[
    OptionValue["ImaginaryDetour"], Length[variables]];
  lineRules = Thread[
    variables ->
      (fromValues + lineParameter (toValues - fromValues) +
        I detour["Amplitude"] 4 lineParameter (1 - lineParameter) *
          detour["Directions"])
  ];
  activeIndices = Select[
    Range[system["LetterCount"]],
    !zeroCanonicalMatrixQ[system["ConstantMatrices"][[#]]] &
  ];
  activeLetters = system["Letters"][[activeIndices]];
  lineLetters = Together[activeLetters /. lineRules];
  polynomials = canonicalSingularityPolynomials[
    lineLetters, lineParameter];
  singularities = canonicalSingularities[
    polynomials, lineParameter, workingPrecision];
  charts = clearanceCharts[
    singularities, safety, workingPrecision, searchSubdivisions];
  Join[
    <|
      "Schema" -> "DiffExp2.CanonicalChartGeometry/v1",
      "ImaginaryDetour" -> detour,
      "ClearanceFactor" -> safety,
      "ActiveLetterIndices" -> activeIndices,
      "SingularityPolynomials" -> polynomials,
      "Singularities" -> singularities
    |>,
    charts
  ]
];

CanonicalLineChartGeometry[x___] :=
  err["E6",
    "CanonicalLineChartGeometry expects a canonical system and two endpoints",
    <|"Arguments" -> {x}|>];

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

squareRootBase[Power[base_, Rational[1, 2]]] := base;

trackedSquareRootSigns[roots_List, centers_List,
    lineParameter_Symbol, workingPrecision_Integer,
    substeps_Integer] := Module[
  {bases, previousPoint = 0, previousBases, continuousRoots,
    signs = {}, target, grid, nextBases, principalRoots, chartSigns},
  If[roots === {},
    Return[ConstantArray[{}, Length[centers]]]];
  bases = squareRootBase /@ roots;
  previousBases = materialize[
    bases /. lineParameter -> previousPoint, workingPrecision];
  continuousRoots = materialize[
    roots /. lineParameter -> previousPoint, workingPrecision];
  If[!AllTrue[Join[previousBases, continuousRoots], finiteNumericQ] ||
      AnyTrue[previousBases, TrueQ[PossibleZeroQ[#]] &],
    err["E7",
      "a square-root branch cannot be initialized at the path start",
      <|"RootValues" -> continuousRoots,
        "RadicandValues" -> previousBases|>]];
  Do[
    target = centers[[chartIndex]];
    grid = Table[
      previousPoint +
        (target - previousPoint) stepIndex/substeps,
      {stepIndex, 1, substeps}
    ];
    Do[
      nextBases = materialize[
        bases /. lineParameter -> point, workingPrecision];
      If[!AllTrue[nextBases, finiteNumericQ] ||
          AnyTrue[nextBases, TrueQ[PossibleZeroQ[#]] &],
        err["E7",
          "a square-root radicand vanished during branch tracking",
          <|"Point" -> point, "RadicandValues" -> nextBases|>]];
      continuousRoots = materialize[
        continuousRoots *
          Exp[Log[nextBases/previousBases]/2],
        workingPrecision
      ];
      previousBases = nextBases,
      {point, grid}
    ];
    principalRoots = materialize[
      roots /. lineParameter -> target, workingPrecision];
    chartSigns = MapThread[
      If[Abs[#1 - #2] <= Abs[#1 + #2], 1, -1] &,
      {continuousRoots, principalRoots}
    ];
    AppendTo[signs, chartSigns];
    continuousRoots = materialize[
      chartSigns principalRoots, workingPrecision];
    previousPoint = target,
    {chartIndex, Length[centers]}
  ];
  signs
];

chartDlogs[lineLetters_List, centers_List, lineParameter_Symbol,
    workingPrecision_Integer, substeps_Integer] := Module[
  {roots, signs, replacements, signedLetters},
  roots = DeleteDuplicates@Cases[
    lineLetters, Power[_, Rational[1, 2]], Infinity];
  signs = trackedSquareRootSigns[
    roots, centers, lineParameter, workingPrecision, substeps];
  <|
    "Roots" -> roots,
    "Signs" -> signs,
    "Dlogs" -> Table[
      replacements = MapThread[
        Rule, {roots, signs[[chartIndex]] roots}];
      signedLetters = lineLetters /. replacements;
      Map[Together[D[Log[#], lineParameter]] &, signedLetters],
      {chartIndex, Length[centers]}
    ]
  |>
];

matrixCombination[values_List, matrices_List, dimension_Integer,
    workingPrecision_Integer] := Module[{combined},
  combined = Total[MapThread[#1 #2 &, {values, matrices}]];
  If[combined === 0,
    combined = SparseArray[{}, {dimension, dimension}]];
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

seriesForcing[coefficientMatrices_List, previousCoefficients_List,
    order_Integer, dimension_Integer,
    workingPrecision_Integer] := Module[
  {forcing = ConstantArray[0, {order, dimension}], contribution},
  Do[
    contribution = Transpose[
      coefficientMatrices[[matrixPower + 1]] .
        Transpose[
          previousCoefficients[[1 ;; order - matrixPower]]
        ]
    ];
    forcing[[matrixPower + 1 ;; order]] += contribution,
    {matrixPower, 0, order - 1}
  ];
  materialize[forcing, workingPrecision]
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
      forcing = seriesForcing[
        coefficientMatrices, coefficients[[epsilonIndex]],
        expansionOrder, dimension, workingPrecision
      ];
      coefficients[[epsilonIndex + 1, 2 ;; expansionOrder + 1]] =
        materialize[
          forcing/Range[expansionOrder],
          workingPrecision
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
    coefficientMatrices, coefficients, forcing, residueForcing,
    output, evaluated,
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
      forcing = seriesForcing[
        coefficientMatrices, coefficients[[epsilonIndex]],
        expansionOrder, dimension, workingPrecision
      ];
      residueForcing = Transpose[
        residueMatrix .
          Transpose[
            coefficients[[
              epsilonIndex, 2 ;; expansionOrder + 1]]
          ]
      ];
      coefficients[[epsilonIndex + 1, 2 ;; expansionOrder + 1]] =
        materialize[
          (forcing + residueForcing)/Range[expansionOrder],
          workingPrecision
        ],
      {epsilonIndex, 1, epsilonOrder}
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
    detour, branchTrackingSubsteps,
    lineParameter = Unique["canonicalLine$"], lineRules,
    lineLetters, branchData, blocks, chartResult, chartRecords = {},
    totalTiming},
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
  branchTrackingSubsteps = resolvedInteger[
    OptionValue["BranchTrackingSubsteps"],
    "BranchTrackingSubsteps", 2];
  {centers, boundaries} = chartGeometry[
    OptionValue["ChartCenters"], OptionValue["ChartBoundaries"]];
  fromValues = pointValues[variables, from];
  toValues = pointValues[variables, to];
  detour = resolvedDetour[
    OptionValue["ImaginaryDetour"], Length[variables]];
  lineRules = Thread[
    variables ->
      (fromValues + lineParameter (toValues - fromValues) +
        I detour["Amplitude"] 4 lineParameter (1 - lineParameter) *
          detour["Directions"])
  ];
  lineLetters = Together[letters /. lineRules];
  branchData = chartDlogs[
    lineLetters, centers, lineParameter, workingPrecision,
    branchTrackingSubsteps];
  blocks = materialize[
    Transpose[boundary[[All, 1 ;; epsilonOrder + 1]]],
    workingPrecision
  ];
  totalTiming = AbsoluteTiming[
    Do[
      DiffExp2`Config`PrintInfo[
        2, "DiffExp2 canonical chart ", chartIndex, "/",
        Length[centers], " centered at ", centers[[chartIndex]]];
      chartResult = If[chartIndex === 1 &&
          TrueQ[PossibleZeroQ[centers[[1]] - boundaries[[1]]]],
        singularStartChart[
          blocks, centers[[1]], boundaries[[2]],
          branchData["Dlogs"][[chartIndex]], matrices,
          dimension, epsilonOrder, expansionOrder, workingPrecision,
          usePade, lineParameter
        ],
        regularChart[
          blocks, boundaries[[chartIndex]], centers[[chartIndex]],
          boundaries[[chartIndex + 1]],
          branchData["Dlogs"][[chartIndex]], matrices, dimension,
          epsilonOrder, expansionOrder, workingPrecision, usePade,
          lineParameter
        ]
      ];
      blocks = chartResult[[1]];
      AppendTo[chartRecords,
        Join[<|"Index" -> chartIndex|>, chartResult[[2]]]];
      DiffExp2`Config`PrintInfo[
        2, "DiffExp2 canonical chart ", chartIndex,
        " completed in ", chartResult[[2, "TimingSeconds"]],
        " seconds"],
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
    "ImaginaryDetour" -> detour,
    "SquareRootCount" -> Length[branchData["Roots"]],
    "SquareRootSigns" -> branchData["Signs"],
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
