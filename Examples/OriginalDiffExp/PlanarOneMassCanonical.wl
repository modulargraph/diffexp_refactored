(* DiffExp 2 reproduction of the four canonical systems in
   Examples/5pPlanar1Mass.nb from DiffExp 1.

   Physics source:
     D. Chicherina, V. Sotnikov, S. Zoia,
     "Pentagon Functions for One-Mass Planar Scattering Amplitudes",
     arXiv:2005.04195.

   Fetch the hash-pinned ancillary files with
     Scripts/fetch_planar_one_mass_data.sh

   PLANAR_ONE_MASS_FAMILY may be 1loop, zmz, mzz, zzz, or all. *)

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
  environmentValue["PLANAR_ONE_MASS_DATA"],
  Except[_String] :> FileNameJoin[{
    DirectoryName[$InputFileName], "Data", "PlanarOneMass"}]
];
familySelection = Replace[
  environmentValue["PLANAR_ONE_MASS_FAMILY"],
  Except[_String] -> "all"
];
families = Switch[
  familySelection,
  "all", {"1loop", "zmz", "mzz", "zzz"},
  "1loop" | "zmz" | "mzz" | "zzz", {familySelection},
  _,
    Print["PLANAR_ONE_MASS invalid family: ", familySelection];
    Exit[2]
];
workingPrecision = integerEnvironment[
  "PLANAR_ONE_MASS_WORKING_PRECISION", 50];
expansionOrder = integerEnvironment[
  "PLANAR_ONE_MASS_EXPANSION_ORDER", 25];
epsilonOrder = 4;

requiredFiles = Join[
  {"alphabet.m"},
  Flatten[
    {
      FileNameJoin[{#, "diffEq-" <> # <> ".m"}],
      FileNameJoin[{#, "numIntegrals-" <> # <> ".m"}]
    } & /@ families
  ]
];
missingFiles = Select[
  FileNameJoin[{dataDirectory, #}] & /@ requiredFiles,
  !FileExistsQ[#] &
];
If[missingFiles =!= {},
  Print["Missing planar one-mass ancillary data: ", missingFiles];
  Print["Run Scripts/fetch_planar_one_mass_data.sh or set ",
    "PLANAR_ONE_MASS_DATA."];
  Exit[2]
];

{s12, s15, s23, s34, s45, p1s} =
  {Global`s12, Global`s15, Global`s23, Global`s34, Global`s45,
    Global`p1s};
variables = {s12, s15, s23, s34, s45, p1s};
tr5Expression = Sqrt[
  (-s12 s15 + s12 s23 + p1s s34 + s15 s45 -
      s34 s45 - s23 s34)^2 -
  4 s23 s34 s45 (p1s - s12 - s15 + s34)];
sqrtG3Expression = Sqrt[
  p1s^2 + (s23 - s45)^2 - 2 p1s (s23 + s45)];
sqrtG3ncExpression = Sqrt[
  (s12 + s15)^2 - 4 p1s s34];
letters = Import[
  FileNameJoin[{dataDirectory, "alphabet.m"}]] /.
  {
    Global`tr5 -> tr5Expression,
    Global`sqrtG3 -> sqrtG3Expression,
    Global`sqrtG3nc -> sqrtG3ncExpression
  };

(* All directions are positive, so this is in the published +i0 homotopy.
   The small generic separation prevents cancellations from leaving an
   active singularity exactly on the real line parameter. *)
detour = <|
  "Amplitude" -> 1/10,
  "Directions" -> {1, 11/10, 6/5, 13/10, 7/5, 3/2}
|>;
historicalDiffExp1Seconds = <|
  "1loop" -> 85.556919,
  "zmz" -> 1391.141062,
  "mzz" -> 992.459149,
  "zzz" -> 1670.664508
|>;

allPassed = True;
Do[
  preparationSeconds = AbsoluteTiming[
    matrixTensor = Import[
      FileNameJoin[{
        dataDirectory, family, "diffEq-" <> family <> ".m"}]];
    samples = Import[
      FileNameJoin[{
        dataDirectory, family, "numIntegrals-" <> family <> ".m"}]];
    canonicalSystem = DiffExp2`LoadCanonicalSystem[<|
      "ConstantMatrices" -> Table[
        matrixTensor[[index]],
        {index, First[Dimensions[matrixTensor]]}
      ],
      "Letters" -> letters,
      "Variables" -> variables
    |>];
    boundary = Transpose[samples[[1, 2]]];
    reference = Transpose[samples[[6, 2]]];
  ][[1]];

  geometrySeconds = AbsoluteTiming[
    geometry = DiffExp2`CanonicalLineChartGeometry[
      canonicalSystem,
      samples[[1, 1]],
      samples[[6, 1]],
      "WorkingPrecision" -> workingPrecision,
      "ImaginaryDetour" -> detour,
      "ClearanceFactor" -> 1/3
    ];
  ][[1]];

  transportWallSeconds = AbsoluteTiming[
    result = DiffExp2`TransportCanonicalLine[
      canonicalSystem,
      boundary,
      samples[[1, 1]],
      samples[[6, 1]],
      "ExpansionOrder" -> expansionOrder,
      "EpsilonOrder" -> epsilonOrder,
      "WorkingPrecision" -> workingPrecision,
      "UsePade" -> False,
      "ImaginaryDetour" -> detour,
      "ChartCenters" -> geometry["Centers"],
      "ChartBoundaries" -> geometry["Boundaries"]
    ];
  ][[1]];

  errorsByOrder = Table[
    Max[Abs[
      result["Value"][[All, order + 1]] -
      reference[[All, order + 1]]
    ]],
    {order, 0, epsilonOrder}
  ];
  maxError = Max[errorsByOrder];
  comparableSeconds = geometrySeconds + transportWallSeconds;
  speedup =
    historicalDiffExp1Seconds[family]/comparableSeconds;

  Print[
    "PLANAR_ONE_MASS family=", family,
    " dimension=", result["Dimension"],
    " letters=", canonicalSystem["LetterCount"],
    " activeLetters=", Length[geometry["ActiveLetterIndices"]],
    " charts=", Length[geometry["Centers"]],
    " preparationSeconds=", N[preparationSeconds, 8],
    " geometrySeconds=", N[geometrySeconds, 8],
    " transportSeconds=", N[transportWallSeconds, 8],
    " recurrenceSeconds=", N[result["TimingSeconds"], 8],
    " comparableSeconds=", N[comparableSeconds, 8],
    " originalDiffExp1Seconds=",
      historicalDiffExp1Seconds[family],
    " speedup=", N[speedup, 6]
  ];
  Print[
    "PLANAR_ONE_MASS family=", family,
    " errorsByOrder=", InputForm[N[errorsByOrder, 12]],
    " maxError=", InputForm[N[maxError, 12]],
    " padeFallbacks=", result["PadeFallbacks"]
  ];
  If[!TrueQ[maxError < 10^-8],
    allPassed = False;
    Print[
      "PLANAR_ONE_MASS FAIL family=", family,
      ": reference error exceeds 1e-8"
    ],
    Print["PLANAR_ONE_MASS PASS family=", family]
  ],
  {family, families}
];

If[allPassed,
  Print["PLANAR_ONE_MASS PASS"];
  Exit[0],
  Exit[1]
];
