(* DiffExp 2 reproduction of the unequal-mass part of DiffExp 1's
   Examples/Banana.nb.

   The public notebook computes

     (mm1,mm2,mm3,mm4) = (2,3/2,4/3,1),  psq = 50

   by two independent contours.  This runnable example uses its faster
   crosscheck contour: transport the equal-mass solution to psq=50 first,
   lift the four equal-mass masters into the 15-master unequal family, then
   deform the masses at fixed momentum.  The separately qualified main
   contour (mass deformation at psq=1/2 followed by momentum transport)
   agrees with this endpoint through eps^7.

   Exact matrices are fetched by Scripts/fetch_original_banana_data.sh from
   DiffExp revision 784c8229bf92369a03f011a48e161522c8c54bbd. *)

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
SetAttributes[catchDE2, HoldFirst];
catchDE2[expression_] := Quiet[
  Catch[expression, "DiffExp2Error"]];

dataDirectory = Replace[
  environmentValue["ORIGINAL_BANANA_DATA"],
  Except[_String] :> FileNameJoin[{
    DirectoryName[$InputFileName], "Data", "Banana"}]
];
equalMassDirectory = FileNameJoin[{dataDirectory, "EqualMass"}];
unequalMassDirectory = FileNameJoin[{dataDirectory, "UnequalMass"}];
referenceDirectory = FileNameJoin[{
  DirectoryName[$InputFileName], "Reference"}];

equalFiles = {"dt_0.m", "dt_1.m"};
unequalFiles = {
  "dmm1_0.m", "dmm1_1.m",
  "dmm2_0.m", "dmm2_1.m",
  "dmm3_0.m", "dmm3_1.m"
};
requiredFiles = Join[
  FileNameJoin[{equalMassDirectory, #}] & /@ equalFiles,
  FileNameJoin[{unequalMassDirectory, #}] & /@ unequalFiles,
  {
    FileNameJoin[{
      referenceDirectory, "BananaBoundaryAtMinusOneEps7.m"}],
    FileNameJoin[{
      referenceDirectory, "BananaUnequalMassAt50.m"}]
  }
];
missingFiles = Select[requiredFiles, !FileExistsQ[#] &];
If[missingFiles =!= {},
  Print["Missing original unequal-mass Banana data: ", missingFiles];
  Print["Run Scripts/fetch_original_banana_data.sh or set ",
    "ORIGINAL_BANANA_DATA."];
  Exit[2]
];

workingPrecision = integerEnvironment[
  "ORIGINAL_BANANA_UNEQUAL_WORKING_PRECISION", 100];
expansionOrder = integerEnvironment[
  "ORIGINAL_BANANA_UNEQUAL_EXPANSION_ORDER", 50];
epsilonOrder = 7;
referenceEpsilonOrder = 4;

DiffExp2`LoadConfiguration[{
  "RecurrenceBackend" -> "Cpp",
  "WorkingPrecision" -> workingPrecision,
  "LinearSolveChopPrecision" -> 10,
  "ExpansionOrder" -> expansionOrder,
  "EpsilonOrder" -> epsilonOrder,
  "DivisionOrder" -> 2,
  "UsePade" -> True,
  "Verbosity" -> 0
}];

{t, x, eps, psq, mm1, mm2, mm3, mm4} = {
  Global`t, Global`x, Global`eps, Global`psq,
  Global`mm1, Global`mm2, Global`mm3, Global`mm4
};

boundaryAtMinusOne = Get[FileNameJoin[{
  referenceDirectory, "BananaBoundaryAtMinusOneEps7.m"}]];
referenceAt50 = Get[FileNameJoin[{
  referenceDirectory, "BananaUnequalMassAt50.m"}]];

equalMatrix =
  Get[FileNameJoin[{equalMassDirectory, equalFiles[[1]]}]] +
  eps Get[FileNameJoin[{equalMassDirectory, equalFiles[[2]]}]];

dmm1 =
  Get[FileNameJoin[{unequalMassDirectory, unequalFiles[[1]]}]] +
  eps Get[FileNameJoin[{unequalMassDirectory, unequalFiles[[2]]}]];
dmm2 =
  Get[FileNameJoin[{unequalMassDirectory, unequalFiles[[3]]}]] +
  eps Get[FileNameJoin[{unequalMassDirectory, unequalFiles[[4]]}]];
dmm3 =
  Get[FileNameJoin[{unequalMassDirectory, unequalFiles[[5]]}]] +
  eps Get[FileNameJoin[{unequalMassDirectory, unequalFiles[[6]]}]];

plainValue[value_] := Transpose[
  Table[
    DiffExp2`EpsilonCoefficient[value, order],
    {order, 0, epsilonOrder}
  ]
];

transportLeg[boundary_, matrix_, variable_Symbol, from_, to_,
    name_String] := Module[
  {line, pulled, system, result, wallSeconds},
  line = from + (to - from) x;
  pulled = Map[
    Together,
    (matrix /. variable -> line) (to - from),
    {2}
  ];
  system = DiffExp2`LoadSystem[<|
    "Matrix" -> pulled,
    "Variable" -> x
  |>];
  wallSeconds = AbsoluteTiming[
    result = catchDE2[
      DiffExp2`TransportEndpoint[system, boundary, 0, 1]];
  ][[1]];
  If[FailureQ[result],
    Print["ORIGINAL_BANANA_UNEQUAL FAIL leg=", name, " ", result];
    Exit[1]
  ];
  Print[
    "ORIGINAL_BANANA_UNEQUAL leg=", name,
    " seconds=", N[wallSeconds, 8],
    " segments=", result["SegmentCount"]
  ];
  {plainValue[result["Value"]], wallSeconds}
];

(* The upper-half-plane contour avoids the equal-mass thresholds at
   psq=0,4,16 and has the same physical continuation as the notebook's
   Mobius-segmented real path. *)
equalUp = transportLeg[
  boundaryAtMinusOne, equalMatrix, t,
  -1, -1 + 10 I, "equal-up"];
equalAcross = transportLeg[
  equalUp[[1]], equalMatrix, t,
  -1 + 10 I, 50 + 10 I, "equal-across"];
equalDown = transportLeg[
  equalAcross[[1]], equalMatrix, t,
  50 + 10 I, 50, "equal-down"];

(* At equal masses the 15 unequal-family masters collapse onto these four
   equal-family masters in the exact order used by Banana.nb. *)
equalMassLift = {
  1, 1, 1, 1, 1, 1,
  2, 2, 2, 2,
  3,
  4, 4, 4, 4
};
unequalBoundaryAt50 = equalDown[[1]][[equalMassLift]];

(* mm1=1+x, mm2=1+x/2, mm3=1+x/3, mm4=1 at fixed psq=50. *)
massLineMatrix = Map[
  Together,
  (dmm1 + dmm2/2 + dmm3/3) /. {
    psq -> 50,
    mm1 -> 1 + x,
    mm2 -> 1 + x/2,
    mm3 -> 1 + x/3,
    mm4 -> 1
  },
  {2}
];
massSystem = DiffExp2`LoadSystem[<|
  "Matrix" -> massLineMatrix,
  "Variable" -> x
|>];
massSeconds = AbsoluteTiming[
  massResult = catchDE2[
    DiffExp2`TransportEndpoint[
      massSystem, unequalBoundaryAt50, 0, 1]];
][[1]];
If[FailureQ[massResult],
  Print["ORIGINAL_BANANA_UNEQUAL FAIL leg=mass-deformation ",
    massResult];
  Exit[1]
];
Print[
  "ORIGINAL_BANANA_UNEQUAL leg=mass-deformation",
  " seconds=", N[massSeconds, 8],
  " segments=", massResult["SegmentCount"]
];

valueAt50 = plainValue[massResult["Value"]];
errorsByOrder = Table[
  Max[Abs[
    valueAt50[[All, order + 1]] -
    referenceAt50[[All, order + 1]]
  ]],
  {order, 0, referenceEpsilonOrder}
];
maxError = Max[errorsByOrder];

equalTransportSeconds = Total[{
  equalUp[[2]], equalAcross[[2]], equalDown[[2]]
}];
totalSeconds = equalTransportSeconds + massSeconds;

(* Like-for-like saved DiffExp 1 crosscheck timings.  Its separate
   asymptotic-boundary-to-t=-1 call (33.272064 s) is excluded because this
   example, like BananaEqualMass.wl, starts from the numerical t=-1 seed. *)
historicalDiffExp1EqualMinusOneTo50Seconds = 76.018792;
historicalDiffExp1MassDeformationSeconds = 244.617557;
historicalDiffExp1ComparableSeconds =
  historicalDiffExp1EqualMinusOneTo50Seconds +
  historicalDiffExp1MassDeformationSeconds;
speedup = historicalDiffExp1ComparableSeconds/totalSeconds;

Print[
  "ORIGINAL_BANANA_UNEQUAL dimension=15",
  " epsilonOrder=", epsilonOrder,
  " checkedEpsilonOrder=", referenceEpsilonOrder,
  " expansionOrder=", expansionOrder,
  " equalTransportSeconds=", N[equalTransportSeconds, 8],
  " massDeformationSeconds=", N[massSeconds, 8],
  " totalSeconds=", N[totalSeconds, 8]
];
Print[
  "ORIGINAL_BANANA_UNEQUAL originalDiffExp1ComparableSeconds=",
  historicalDiffExp1ComparableSeconds,
  " speedup=", N[speedup, 6]
];
Print[
  "ORIGINAL_BANANA_UNEQUAL errorsByOrder=",
  InputForm[N[errorsByOrder, 12]],
  " maxError=", InputForm[N[maxError, 12]]
];

If[!TrueQ[maxError < 10^-8],
  Print[
    "ORIGINAL_BANANA_UNEQUAL FAIL: psq=50 reference error exceeds 1e-8"];
  Exit[1]
];

Print["ORIGINAL_BANANA_UNEQUAL PASS"];
Exit[0];
