(* DiffExp 2 reproduction of Examples/5pNonPlanar.nb from DiffExp 1.

   Physics source:
     D. Chicherina, T. Gehrmann, J. M. Henn, P. Wasser, Y. Zhang,
     S. Zoia, "All Master Integrals for Three-Jet Production at NNLO",
     arXiv:1812.11160.

   The ancillary matrix and boundary files are fetched separately by
     Scripts/fetch_henn_nonplanar_data.sh
   so their upstream provenance and hashes remain explicit. *)

repoRoot = ParentDirectory[
  ParentDirectory[DirectoryName[$InputFileName]]];
Get[FileNameJoin[{repoRoot, "DiffExp2.m"}]];

environmentValue[name_String] := Quiet[
  Check[Environment[name], $Failed]];
integerEnvironment[name_String, default_Integer] := Module[{value},
  value = environmentValue[name];
  If[StringQ[value] && StringMatchQ[value, DigitCharacter ..],
    FromDigits[value], default]
];

dataDirectory = Replace[
  environmentValue["HENN_NONPLANAR_DATA"],
  Except[_String] :> FileNameJoin[{
    DirectoryName[$InputFileName], "Data", "HennNonplanar"}]
];
requiredFiles = {
  "XB_Atilde.txt",
  "XB_Boundary_values_X0.txt",
  "XB_Boundary_values_X1.txt"
};
missingFiles = Select[
  FileNameJoin[{dataDirectory, #}] & /@ requiredFiles,
  !FileExistsQ[#] &
];
If[missingFiles =!= {},
  Print["Missing Henn nonplanar ancillary data: ", missingFiles];
  Print["Run Scripts/fetch_henn_nonplanar_data.sh or set ",
    "HENN_NONPLANAR_DATA."];
  Exit[2]
];

workingPrecision = integerEnvironment[
  "HENN_NONPLANAR_WORKING_PRECISION", 50];
expansionOrder = integerEnvironment[
  "HENN_NONPLANAR_EXPANSION_ORDER", 25];
epsilonOrder = 4;

{v1, v2, v3, v4, v5} =
  {Global`v1, Global`v2, Global`v3, Global`v4, Global`v5};
gramDelta =
  v1^2 (v2 - v5)^2 +
  (v2 v3 + v4 (-v3 + v5))^2 +
  2 v1 (
    -v2^2 v3 +
    v4 (v3 - v5) v5 +
    v2 (v3 v4 + (v3 + v4) v5)
  );
gramRoot = Sqrt[gramDelta];
oddLetterBases = {
  v1 v2 - v2 v3 + v3 v4 - v1 v5 - v4 v5,
  -v1 v2 + v2 v3 - v3 v4 - v1 v5 + v4 v5,
  -v1 v2 - v2 v3 + v3 v4 + v1 v5 - v4 v5,
  v1 v2 - v2 v3 - v3 v4 - v1 v5 + v4 v5,
  -v1 v2 + v2 v3 - v3 v4 + v1 v5 - v4 v5
};
pentagonLetters = Join[
  {
    v1, v2, v3, v4, v5,
    v3 + v4, v4 + v5, v1 + v5, v1 + v2, v2 + v3,
    v1 - v4, v2 - v5, -v1 + v3, -v2 + v4, -v3 + v5,
    v1 + v2 - v4,
    v2 + v3 - v5,
    -v1 + v3 + v4,
    -v2 + v4 + v5,
    v1 - v3 + v5,
    -v1 - v2 + v3 + v4,
    -v2 - v3 + v4 + v5,
    v1 - v3 - v4 + v5,
    v1 + v2 - v4 - v5,
    -v1 + v2 + v3 - v5
  },
  (# - gramRoot)/(# + gramRoot) & /@ oddLetterBases,
  {gramRoot}
];

preparationTiming = AbsoluteTiming[
  potential = ToExpression[
    Import[FileNameJoin[{dataDirectory, requiredFiles[[1]]}]]];
  constantMatrices = Table[
    SparseArray[
      Map[
        Coefficient[#, Log[Global`W[index]]] &,
        potential,
        {2}
      ]
    ],
    {index, Length[pentagonLetters]}
  ];
  canonicalSystem = DiffExp2`LoadCanonicalSystem[<|
    "ConstantMatrices" -> constantMatrices,
    "Letters" -> pentagonLetters,
    "Variables" -> {v1, v2, v3, v4, v5}
  |>];
  boundary = ToExpression[
    Import[FileNameJoin[{dataDirectory, requiredFiles[[2]]}]]];
  reference = ToExpression[
    Import[FileNameJoin[{dataDirectory, requiredFiles[[3]]}]]];
][[1]];

x0 = {
  v1 -> 3,
  v2 -> -1,
  v3 -> 1,
  v4 -> 1,
  v5 -> -1
};
x1 = {
  v1 -> 4,
  v2 -> -(113/47),
  v3 -> 281/149,
  v4 -> 349/257,
  v5 -> -(863/541)
};

result = DiffExp2`TransportCanonicalLine[
  canonicalSystem,
  boundary,
  x0,
  x1,
  "ExpansionOrder" -> expansionOrder,
  "EpsilonOrder" -> epsilonOrder,
  "WorkingPrecision" -> workingPrecision,
  "UsePade" -> True,
  "ChartCenters" -> {0, 1/3, 7/10, 19/20},
  "ChartBoundaries" -> {
    0,
    0.19698903994874387,
    0.5085056419319689,
    0.8409884239904127,
    1
  }
];

errorsByOrder = Table[
  Max[Abs[
    result["Value"][[All, order + 1]] -
    reference[[All, order + 1]]
  ]],
  {order, 0, epsilonOrder}
];
maxError = Max[errorsByOrder];
historicalDiffExp1Seconds = 49.765221;
speedup = historicalDiffExp1Seconds/result["TimingSeconds"];

Print["HENN_NONPLANAR dimension=", result["Dimension"],
  " letters=", canonicalSystem["LetterCount"],
  " epsilonOrder=", epsilonOrder,
  " expansionOrder=", expansionOrder];
Print["HENN_NONPLANAR preparationSeconds=",
  N[preparationTiming, 8],
  " transportSeconds=", N[result["TimingSeconds"], 8],
  " totalSeconds=", N[preparationTiming + result["TimingSeconds"], 8]];
Print["HENN_NONPLANAR originalDiffExp1TransportSeconds=",
  historicalDiffExp1Seconds,
  " transportSpeedup=", N[speedup, 6]];
Print["HENN_NONPLANAR errorsByOrder=",
  InputForm[N[errorsByOrder, 12]],
  " maxError=", InputForm[N[maxError, 12]]];
Print["HENN_NONPLANAR padeFallbacks=", result["PadeFallbacks"]];

If[!TrueQ[maxError < 10^-9],
  Print["HENN_NONPLANAR FAIL: reference error exceeds 1e-9"];
  Exit[1]
];

Print["HENN_NONPLANAR PASS"];
Exit[0];
