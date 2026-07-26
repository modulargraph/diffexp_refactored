(* DiffExp 2 reproduction of Examples/MultiplePolylogarithms.nb.

   A weight-w Goncharov polylogarithm is embedded as the eps^w coefficient
   of an eps-canonical triangular dlog system. This makes ordinary GPL log
   chains exactly match DiffExp 2's (eps Log[t])^p sector normalization. *)

repoRoot = ParentDirectory[
  ParentDirectory[DirectoryName[$InputFileName]]];
Get[FileNameJoin[{repoRoot, "DiffExp2.m"}]];

environmentValue[name_String] := Quiet[
  Check[Environment[name], $Failed]];
booleanEnvironment[name_String, default_] := Replace[
  environmentValue[name],
  {
    "0" | "false" | "False" -> False,
    "1" | "true" | "True" -> True,
    _ -> default
  }
];
includeWeight20 = booleanEnvironment[
  "MPL_INCLUDE_WEIGHT20", True];

t = Global`t;
ClearAll[G, GFB, ListShuffleProduct, ShuffleG];

(* These shuffle-regularization rules are the rules used by the original
   notebook to remove trailing-zero base-point logarithms. *)
ListShuffleProduct[
    u : {a_, x___}, v : {b_, y___}, c___] :=
  Join[
    ListShuffleProduct[{x}, v, c, a],
    ListShuffleProduct[u, {y}, c, b]
  ];
ListShuffleProduct[{x___}, {y___}, c___] := {{c, x, y}};
ShuffleG[G[a__, x_], G[b__, x_]] :=
  Total[
    G[Sequence @@ #, x] & /@
      ListShuffleProduct[{a}, {b}]
  ];
GFB[a___, a1_, b___, 0, x_] /;
    !a1 === 0 && b === 0 :=
  G[a, a1, b, 0, x] +
    G[a, a1, x] G[b, 0, x] -
    ShuffleG[G[a, a1, x], G[b, 0, x]] /. G -> GFB;
GFB[a___, a1_, x_] /; !a1 === 0 := G[a, a1, x];
GFB[a___, x_] /; a === 0 := Module[
  {length = Length[{a}]},
  Log[x]^length/length!
];

$MPLExpansionOrder = 50;
$MPLWorkingPrecision = 100;
$MPLClearanceFactor = 1/3;
$MPLGeometrySeconds = 0;
$MPLTransportSeconds = 0;
$MPLCharts = 0;

ClearAll[mplFiniteBasepoint];
mplFiniteBasepoint[arguments__] := Module[
  {all = {arguments}, indices, endpoint, weight, dimension,
    constantMatrices, system, boundary, detour, geometry,
    geometrySeconds, result, transportSeconds},
  endpoint = Last[all];
  indices = Most[all];
  weight = Length[indices];
  dimension = weight + 1;
  constantMatrices = Table[
    SparseArray[
      {{index, index + 1} -> 1},
      {dimension, dimension}
    ],
    {index, weight}
  ];
  system = DiffExp2`LoadCanonicalSystem[<|
    "ConstantMatrices" -> constantMatrices,
    "Letters" -> (t - # & /@ indices),
    "Variables" -> {t}
  |>];
  boundary = ConstantArray[0, {dimension, weight + 1}];
  boundary[[-1, 1]] = 1;

  (* A lower-half-plane path gives t-a-I0 for real positive letters,
     hence Mathematica's principal +I Pi value after crossing. *)
  detour = <|"Amplitude" -> 1/10, "Directions" -> {-1}|>;
  geometrySeconds = AbsoluteTiming[
    geometry = DiffExp2`CanonicalLineChartGeometry[
      system,
      {t -> 0},
      {t -> endpoint},
      "WorkingPrecision" -> $MPLWorkingPrecision,
      "ImaginaryDetour" -> detour,
      "ClearanceFactor" -> $MPLClearanceFactor
    ];
  ][[1]];
  transportSeconds = AbsoluteTiming[
    result = DiffExp2`TransportCanonicalLine[
      system,
      boundary,
      {t -> 0},
      {t -> endpoint},
      "ExpansionOrder" -> $MPLExpansionOrder,
      "EpsilonOrder" -> weight,
      "WorkingPrecision" -> $MPLWorkingPrecision,
      "UsePade" -> False,
      "ImaginaryDetour" -> detour,
      "ChartCenters" -> geometry["Centers"],
      "ChartBoundaries" -> geometry["Boundaries"]
    ];
  ][[1]];
  $MPLGeometrySeconds += geometrySeconds;
  $MPLTransportSeconds += transportSeconds;
  $MPLCharts += Length[geometry["Centers"]];
  result["Value"][[1, weight + 1]]
];

ClearAll[mplEvaluate];
mplEvaluate[a__, endpoint_] /; Im[endpoint] === 0 :=
  Expand[
    GFB[a, endpoint] /.
      G[a1__, b_] /; SameQ[a1] :>
        Log[1 - b/{a1}[[1]]]^Length[{a1}]/
          Length[{a1}]! /.
      G[args__] :> mplFiniteBasepoint[args]
  ];

shortCases = {
  <|
    "Name" -> "G[1,0,1;4]",
    "Expression" -> G[1, 0, 1, 4],
    "Reference" ->
      -6.7782180257804207212554826775005988168291802221955692129682 +
      0.9250147943833369547396749852220309435917997631163983727603 I
  |>,
  <|
    "Name" -> "G[1,-10,0;4]",
    "Expression" -> G[1, -10, 0, 4],
    "Reference" ->
      -0.0191508840720296721365611597236750922866172732200324064383 -
      0.3066358899483403657463434439014286049874538907865239005438 I
  |>,
  <|
    "Name" -> "G[10,-10+I,-1/2,-50;1]",
    "Expression" -> G[10, -10 + I, -2^(-1), -50, 1],
    "Reference" ->
      -9.8802442781507281548895360764863423574760704710022738 10^-6 -
      9.352314628872620198852585457725779647560856726857223 10^-7 I
  |>
};
historicalSeconds = <|
  50 -> {0.169930, 0.504031, 0.141181},
  75 -> {0.306044, 0.657636, 0.138740},
  100 -> {38.701777}
|>;

allPassed = True;
Do[
  $MPLExpansionOrder = order;
  $MPLClearanceFactor = 1/3;
  Do[
    $MPLGeometrySeconds = 0;
    $MPLTransportSeconds = 0;
    $MPLCharts = 0;
    wallSeconds = AbsoluteTiming[
      value = shortCases[[index, "Expression"]] /.
        G -> mplEvaluate;
    ][[1]];
    error = Abs[value - shortCases[[index, "Reference"]]];
    tolerance = If[order === 50, 10^-23, 10^-34];
    Print[
      "MPL case=", shortCases[[index, "Name"]],
      " expansionOrder=", order,
      " wallSeconds=", N[wallSeconds, 8],
      " geometrySeconds=", N[$MPLGeometrySeconds, 8],
      " transportSeconds=", N[$MPLTransportSeconds, 8],
      " charts=", $MPLCharts,
      " originalDiffExp1Seconds=",
        historicalSeconds[order][[index]],
      " speedup=",
        N[historicalSeconds[order][[index]]/wallSeconds, 6],
      " error=", InputForm[N[error, 8]],
      " value=", InputForm[value]
    ];
    If[!TrueQ[error < tolerance],
      allPassed = False;
      Print[
        "MPL FAIL case=", shortCases[[index, "Name"]],
        " expansionOrder=", order,
        " tolerance=", tolerance
      ],
      Print[
        "MPL PASS case=", shortCases[[index, "Name"]],
        " expansionOrder=", order
      ]
    ],
    {index, Length[shortCases]}
  ],
  {order, {50, 75}}
];

If[includeWeight20,
  $MPLExpansionOrder = 100;
  $MPLClearanceFactor = 1/2;
  $MPLGeometrySeconds = 0;
  $MPLTransportSeconds = 0;
  $MPLCharts = 0;
  weight20Reference =
    5.13066731719907533179589918813462949766948546641803107466616580097 10^-12;
  wallSeconds = AbsoluteTiming[
    value = G[Sequence @@ Range[20], 21] /. G -> mplEvaluate;
  ][[1]];
  error = Abs[value - weight20Reference];
  Print[
    "MPL case=weight20",
    " expansionOrder=100",
    " wallSeconds=", N[wallSeconds, 8],
    " geometrySeconds=", N[$MPLGeometrySeconds, 8],
    " transportSeconds=", N[$MPLTransportSeconds, 8],
    " charts=", $MPLCharts,
    " originalDiffExp1Seconds=", historicalSeconds[100][[1]],
    " speedup=",
      N[historicalSeconds[100][[1]]/wallSeconds, 6],
    " error=", InputForm[N[error, 8]]
  ];
  If[!TrueQ[error < 10^-40],
    allPassed = False;
    Print["MPL FAIL case=weight20 tolerance=1e-40"],
    Print["MPL PASS case=weight20"]
  ]
];

If[allPassed,
  Print["MPL PASS"];
  Exit[0],
  Exit[1]
];
