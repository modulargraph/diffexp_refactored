(* Reconstruct the X0 boundary of the 108-integral canonical double-pentagon
   basis from arXiv:1812.11160 with the public Feynman-trick pipeline.

   The paper defines

     XB[a1,...,a11] =
       Integral d^D k1 d^D k2/(i Pi)^D Product_i D_i^(-ai)

   with D_i = q_i^2. DiffExp 2 uses D_i' = -q_i^2, so each raw result is
   multiplied by (-1)^Total[{a1,...,a11}] before the published canonical
   linear combinations are formed. The paper's final normalization is
   eps^4 Exp[2 eps EulerGamma].

   The maintained default computes canonical component 1, which exercises a
   genuine four-level two-loop Feynman-trick ladder and is compared through
   epsilon^4.  --all requests the complete 108-component reconstruction.

   Usage:

     wolframscript -script Examples/FeynmanTrick/HennDoublePentagonBoundary.wl
     wolframscript -script Examples/FeynmanTrick/HennDoublePentagonBoundary.wl --plan
     wolframscript -script Examples/FeynmanTrick/HennDoublePentagonBoundary.wl --component 1
     wolframscript -script Examples/FeynmanTrick/HennDoublePentagonBoundary.wl --all *)

repoRoot = ParentDirectory[
  ParentDirectory[DirectoryName[ExpandFileName[$InputFileName]]]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "FeynmanTrick.m"}]];

environmentOrDefault[name_String, default_String] := Module[{value},
  value = Environment[name];
  If[StringQ[value] && StringLength[StringTrim[value]] > 0,
    StringTrim[value], default]
];

positiveIntegerEnvironment[name_String, default_Integer] := Module[
  {text, value},
  text = environmentOrDefault[name, ToString[default, InputForm]];
  If[!StringMatchQ[text, DigitCharacter ..],
    Print[name, " must be a positive integer; received ", InputForm[text]];
    Exit[2]
  ];
  value = FromDigits[text];
  If[value < 1,
    Print[name, " must be a positive integer; received ", value];
    Exit[2]
  ];
  value
];

arguments = Rest[$ScriptCommandLine];
planOnly = MemberQ[arguments, "--plan"];
arguments = DeleteCases[arguments, "--plan"];
writeManifestPosition = FirstPosition[arguments, "--write-manifest"];
compareManifestPosition = FirstPosition[arguments, "--compare-manifest"];
If[writeManifestPosition =!= Missing["NotFound"] &&
    compareManifestPosition =!= Missing["NotFound"],
  Print["--write-manifest and --compare-manifest are mutually exclusive"];
  Exit[2]
];
executionMode = Which[
  writeManifestPosition =!= Missing["NotFound"], "WriteManifest",
  compareManifestPosition =!= Missing["NotFound"], "CompareManifest",
  True, "Direct"
];
manifestPath = Switch[executionMode,
  "WriteManifest",
    If[writeManifestPosition[[1]] >= Length[arguments],
      Print["--write-manifest requires a destination path"];
      Exit[2]
    ];
    arguments[[writeManifestPosition[[1]] + 1]],
  "CompareManifest",
    If[compareManifestPosition[[1]] >= Length[arguments],
      Print["--compare-manifest requires a manifest path"];
      Exit[2]
    ];
    arguments[[compareManifestPosition[[1]] + 1]],
  _, None
];
If[executionMode =!= "Direct",
  modePosition = If[executionMode === "WriteManifest",
    writeManifestPosition[[1]], compareManifestPosition[[1]]];
  arguments = Delete[arguments, {
    {modePosition}, {modePosition + 1}}]
];
manifestInput = If[executionMode === "CompareManifest",
  Quiet@Check[Import[ExpandFileName[manifestPath], "RawJSON"], $Failed],
  None
];
If[executionMode === "CompareManifest" &&
    (!AssociationQ[manifestInput] ||
      Lookup[manifestInput, "Schema", None] =!=
        "DiffExp2.HennDoublePentagonBoundaryRun/v1"),
  Print["Invalid Henn FT run manifest: ", InputForm[manifestPath]];
  Exit[2]
];
componentSelection = Which[
  executionMode === "CompareManifest",
    Lookup[manifestInput, "Components", $Failed],
  arguments === {}, {1},
  arguments === {"--all"}, All,
  MatchQ[arguments, {"--component", _String}] &&
      StringMatchQ[arguments[[2]], DigitCharacter ..],
    {FromDigits[arguments[[2]]]},
  True,
    Print["Usage: wolframscript -script ", $InputFileName,
      " [--plan] [--all | --component 1..108]"];
    Print["       wolframscript -script ", $InputFileName,
      " --write-manifest path [--all | --component 1..108]"];
    Print["       wolframscript -script ", $InputFileName,
      " --compare-manifest path"];
    Exit[2]
];

dataDirectory = Replace[
  Environment["HENN_NONPLANAR_DATA"],
  Except[_String] :> FileNameJoin[{
    repoRoot, "Examples", "OriginalDiffExp", "Data", "HennNonplanar"}]
];
requiredFiles = {
  "dlogBasisXB.txt",
  "XB_Boundary_values_X0.txt"
};
missingFiles = Select[
  FileNameJoin[{dataDirectory, #}] & /@ requiredFiles,
  !FileExistsQ[#] &
];
If[missingFiles =!= {},
  Print["Missing Henn double-pentagon ancillary data: ", missingFiles];
  Print["Run Scripts/fetch_henn_nonplanar_data.sh or set ",
    "HENN_NONPLANAR_DATA."];
  Exit[2]
];

canonicalBasis = ToExpression[
  Import[FileNameJoin[{dataDirectory, requiredFiles[[1]]}]]];
paperBoundary = ToExpression[
  Import[FileNameJoin[{dataDirectory, requiredFiles[[2]]}]]];
If[Dimensions[paperBoundary] =!= {108, 5} ||
    Length[canonicalBasis] =!= 108,
  Print["Unexpected Henn ancillary shape: basis=", Length[canonicalBasis],
    " boundary=", Dimensions[paperBoundary]];
  Exit[2]
];
If[componentSelection === All,
  componentSelection = Range[Length[canonicalBasis]]];
If[!AllTrue[componentSelection, 1 <= # <= Length[canonicalBasis] &],
  Print["Henn component must lie between 1 and 108; received ",
    componentSelection];
  Exit[2]
];

selectedBasis = canonicalBasis[[componentSelection]];
selectedReference = paperBoundary[[componentSelection]];
selectionSlug = If[
  componentSelection === Range[Length[canonicalBasis]],
  "all",
  "c" <> StringRiffle[ToString /@ componentSelection, "_"]];
scalarTargets = DeleteDuplicates[
  Cases[selectedBasis, Global`XB[a__] :> {a}, Infinity]];

{k1, k2, p1, p2, p3, p4} =
  {Global`k1, Global`k2, Global`p1, Global`p2, Global`p3, Global`p4};
p5 = -p1 - p2 - p3 - p4;

(* X0 = {s12,s23,s34,s45,s15} = {3,-1,1,1,-1}.
   The remaining nonadjacent invariants are all -1. *)
x0MomentumRules = {
  p1^2 -> 0, p2^2 -> 0, p3^2 -> 0, p4^2 -> 0,
  p1 p2 -> 3/2,
  p2 p3 -> -1/2,
  p3 p4 -> 1/2,
  p1 p3 -> -1/2,
  p2 p4 -> -1/2,
  p1 p4 -> -1/2
};
x0InvariantRules = {
  Global`s12 -> 3,
  Global`s23 -> -1,
  Global`s34 -> 1,
  Global`s45 -> 1,
  Global`s15 -> -1,
  Global`eps5 -> I Sqrt[3]
};

familyDefinition = <|
  "LoopMomenta" -> {k1, k2},
  "ExternalMomenta" -> {p1, p2, p3, p4},
  "Propagators" -> -{
    k1^2,
    (-p1 + k1)^2,
    (-p1 - p2 + k1)^2,
    k2^2,
    (p4 + p5 + k2)^2,
    (p5 + k2)^2,
    (k1 - k2)^2,
    (p3 + k1 - k2)^2,
    (p5 + k1)^2,
    (-p1 + k2)^2,
    (-p1 - p2 + k2)^2
  },
  "Replacements" -> x0MomentumRules,
  "Dimension" -> 4 - 2 FeynmanTrick`FTeps
|>;

(* D1,...,D8 are physical propagators in the top sector, but a negative
   power of an absent line acts as an irreducible numerator in a subsector.
   A zero-power line must not be introduced into the Feynman-trick ladder and
   later removed through a singular endpoint limit: its beta exponent is
   exactly zero, so it is absent from the scalar identity.  Group the scalar
   terms by both numerator pattern and positive denominator support.  Each
   exact family request then merges precisely the lines occurring in its
   targets. *)
negativeDenominatorPositions[index_List] :=
  Flatten@Position[Take[index, 8], _?(# < 0 &)];
activeDenominatorPositions[index_List] :=
  Flatten@Position[Take[index, 8], _?(# > 0 &)];
groupPatterns = DeleteDuplicates[
  (<|
      "NumeratorPattern" -> negativeDenominatorPositions[#],
      "ActiveDenominators" -> activeDenominatorPositions[#]
    |> &) /@ scalarTargets
];
targetGroups = Map[
  Function[pattern, Join[pattern, <|
    "Targets" -> Select[
      scalarTargets,
      negativeDenominatorPositions[#] === pattern["NumeratorPattern"] &&
        activeDenominatorPositions[#] ===
          pattern["ActiveDenominators"] &]
  |>]],
  groupPatterns
];

familyForGroup[group_Association] := Module[
  {pattern, activeDenominators, anchor, suffix},
  pattern = group["NumeratorPattern"];
  activeDenominators = group["ActiveDenominators"];
  anchor = First[activeDenominators];
  suffix = If[pattern === {}, "none",
    StringRiffle[ToString /@ pattern, "_"]];
  Join[familyDefinition, <|
    "Name" -> "henn_double_pentagon_x0_n" <> suffix <> "_s" <>
      StringRiffle[ToString /@ activeDenominators, "_"],
    "EliminatedPositions" ->
      Complement[Range[8], Join[activeDenominators, pattern]],
    "NumeratorPositions" -> Join[pattern, {9, 10, 11}],
    "CombinationSequence" ->
      ({anchor, #} & /@ Rest[activeDenominators])
  |>]
];

workingPrecision = positiveIntegerEnvironment[
  "HENN_FT_WORKING_PRECISION", 80];
(* The maintained paper comparison is an eight-digit gate.  The runner
   independently propagates guarded producer accuracy through each
   support-specific ladder and caches the exact per-level profile learned for
   that request. *)
matchingDigits = positiveIntegerEnvironment["HENN_FT_MATCH_DIGITS", 8];
matchingCertificationDigits = positiveIntegerEnvironment[
  "HENN_FT_MATCHING_CERTIFICATION_DIGITS", 6];
expansionOrder = positiveIntegerEnvironment[
  "HENN_FT_EXPANSION_ORDER", 25];
boundaryExtraOrder = positiveIntegerEnvironment[
  "HENN_FT_BOUNDARY_EXTRA_ORDER", 4];
divisionOrder = positiveIntegerEnvironment[
  "HENN_FT_DIVISION_ORDER", 5];
fixedParameterValues = {
  1/5, 3/10, 2/5, 1/2, 3/5, 7/10, 4/5};
cppThreads = positiveIntegerEnvironment[
  "DE2_CPP_THREADS", Min[8, Max[1, $ProcessorCount]]];
(* The X0 point is reached on the lower Feynman rim.  The final x1
   deformation crosses the time-like s12 channel in the opposite projected
   orientation; the remaining parameter contours stay on the lower rim.
   These are physical sheet data, not numerical tuning parameters. *)
deltaPrescriptionSign = -1;
cacheRoot = ExpandFileName[environmentOrDefault[
  "DIFFEXP2_CACHE_DIR",
  FileNameJoin[{$HomeDirectory, ".cache", "diffexp2"}]]];
checkpointRoot = ExpandFileName[environmentOrDefault[
  "HENN_FT_CHECKPOINT_DIR",
  FileNameJoin[{cacheRoot, "ladder"}]]];
firePath = environmentOrDefault["FT_FIRE_PATH", ""];
wolframScript = environmentOrDefault["HENN_FT_WOLFRAMSCRIPT", ""];

pipelineOptionsForGroup[group_Association] := Module[
  {levelCount, levelDeltaPrescriptionSigns},
  levelCount = Length[group["ActiveDenominators"]] - 1;
  levelDeltaPrescriptionSigns = Join[
    If[levelCount >= 1, {1}, {}],
    ConstantArray[-1, Max[0, levelCount - 1]]];
  {
  "WorkingPrecision" -> workingPrecision,
  "MatchingDigits" -> matchingDigits,
  "ExpansionOrder" -> expansionOrder,
  "EpsilonOrder" -> 0,
  "BoundaryExtraOrder" -> boundaryExtraOrder,
  (* Public halos request additional published coefficients.  This boundary
     needs only epsilon^0 after its known two-loop pole normalization; native
     intrinsic-loss planning and private reservoir retries supply internal
     rows without turning them into public matching obligations. *)
  "LevelEpsilonHalos" -> ConstantArray[0, levelCount],
  (* A single repeated anchor would put level 2 on the diagonal
     epsilon-zero singularity x2=x1.  Distinct exact interior anchors keep
     every ladder start regular while preserving real 0..1 contours. *)
  "FixedParameterValues" -> Take[fixedParameterValues, levelCount],
  (* This family crosses several closely spaced real and projected algebraic
     singularities.  Fifth-radius overlaps keep the post-singular regular
     basis well conditioned without asking the Taylor retry loop to repair a
     geometrically marginal one-third-radius handoff. *)
  "DivisionOrder" -> divisionOrder,
  "RadiusOfConvergence" -> 1,
  "RecurrenceBackend" -> "Cpp",
  "CppThreads" -> cppThreads,
  "ValueTransport" -> True,
  (* The direct-value shortcut must first certify a physical ODE disk at the
     long level-3 singular handoff.  On this family that optional proof is
     slower than the already-required basis match (more than 15 minutes
     versus about five) and does not tighten the terminal line enclosure. *)
  "NativeValueHopExecution" -> False,
  "BatchEndpointArms" -> True,
  "DeltaPrescriptionSign" -> deltaPrescriptionSign,
  "LevelDeltaPrescriptionSigns" -> levelDeltaPrescriptionSigns,
  "PreparedCacheDirectory" -> FileNameJoin[{cacheRoot, "fire"}],
  "CheckpointDirectory" ->
    FileNameJoin[{checkpointRoot,
      "henn-double-pentagon-x0-" <> selectionSlug <> "-n" <>
        If[group["NumeratorPattern"] === {}, "none",
          StringRiffle[
            ToString /@ group["NumeratorPattern"], "_"]] <> "-s" <>
        StringRiffle[
          ToString /@ group["ActiveDenominators"], "_"]}],
  (* Full Acb state propagation and endpoint publication remain at eight
     digits.  The measured level-3 chart-to-chart consistency diagnostic has
     a stable 1.7*10^-7 midpoint floor, so only that internal diagnostic is
     capped at six digits.  Upstream producer levels retain their guarded
     certification targets, and the final paper oracle remains the independent
     10^-8 numerical comparison. *)
  "ExtraEnvironment" -> <|
    "FT_MATCHING_CERTIFICATION_DIGITS_BY_LEVEL" ->
      "3:" <> ToString[matchingCertificationDigits, InputForm]|>,
  "WorkingDirectory" -> repoRoot,
  "EchoOutput" -> True
  }
];

plans = Map[
  Function[group,
    Module[{options = pipelineOptionsForGroup[group]},
      If[firePath =!= "",
        options = Append[options, "FIREPath" -> firePath]];
      If[wolframScript =!= "",
        options = Append[options, "WolframScript" -> wolframScript]];
      FeynmanTrick`PipelinePlan[
        familyForGroup[group], group["Targets"],
        Sequence @@ options]
    ]
  ],
  targetGroups
];
If[AnyTrue[plans, FailureQ],
  Print["Could not construct Henn double-pentagon FT plan: ",
    InputForm[Select[plans, FailureQ]]];
  Exit[2]
];

Print["HENN_FT_BOUNDARY components=", Length[componentSelection],
  " scalarTargets=", Length[scalarTargets],
  " numeratorTargets=", Count[scalarTargets, _?(Min[#] < 0 &)],
  " familyGroups=", Length[targetGroups],
  " prescriptionSign=", deltaPrescriptionSign];

If[planOnly,
  MapThread[
    Function[{group, groupPlan},
      Print["HENN_FT_PLAN numeratorPattern=",
        group["NumeratorPattern"],
        " activeDenominators=", group["ActiveDenominators"],
        " targets=", Length[group["Targets"]]];
      Print[InputForm[KeyTake[groupPlan, {
        "FamilyName", "FamilyID", "RequestID", "Settings",
        "PreparedCacheDirectory", "CheckpointDirectory"}]]]
    ],
    {targetGroups, plans}
  ];
  Exit[0]
];

If[executionMode === "WriteManifest",
  manifestDirectory = DirectoryName[ExpandFileName[manifestPath]];
  If[!DirectoryQ[manifestDirectory],
    CreateDirectory[manifestDirectory,
      CreateIntermediateDirectories -> True]];
  requestWrites = Map[
    FeynmanTrick`PipelineRequest`WritePipelineRequest[
      #["Request"], #["RequestFile"]] &,
    plans
  ];
  If[AnyTrue[requestWrites, FailureQ],
    Print["Could not publish Henn FT pipeline request(s): ",
      InputForm[Select[requestWrites, FailureQ]]];
    Exit[2]
  ];
  manifestRuns = MapIndexed[
    Function[{groupPlan, ordinal},
      <|
        "Ordinal" -> First[ordinal],
        "Command" -> groupPlan["Command"],
        "WorkingDirectory" -> groupPlan["WorkingDirectory"],
        "Environment" -> groupPlan["Environment"],
        "LogFile" -> FileNameJoin[{manifestDirectory,
          "henn-ft-run-" <> ToString[First[ordinal]] <> ".log"}],
        "ExpectedMasters" ->
          Lookup[groupPlan["OutputRequests"], "IndexVector"]
      |>
    ],
    plans
  ];
  manifestOutput = <|
    "Schema" -> "DiffExp2.HennDoublePentagonBoundaryRun/v1",
    "Components" -> componentSelection,
    "ScalarTargets" -> scalarTargets,
    "PrescriptionSign" -> deltaPrescriptionSign,
    "Runs" -> manifestRuns
  |>;
  If[Quiet@Check[
      Export[ExpandFileName[manifestPath], manifestOutput, "RawJSON"];
      True, False] =!= True,
    Print["Could not write Henn FT run manifest: ", manifestPath];
    Exit[2]
  ];
  Print["HENN_FT_MANIFEST ", ExpandFileName[manifestPath]];
  Exit[0]
];

If[executionMode === "CompareManifest",
  manifestRuns = Lookup[manifestInput, "Runs", {}];
  parsedOutputs = Map[
    Function[run,
      Module[{logFile, lines, stepwiseLines, finalLines, stepwiseRows,
          finalRows, rows, expected},
        logFile = Lookup[run, "LogFile", ""];
        If[!FileExistsQ[logFile],
          Print["Missing Henn FT runner log: ", logFile];
          Exit[1]
        ];
        lines = Import[logFile, "Lines"];
        stepwiseLines = Select[
          lines, StringStartsQ[#, "STEPWISE "] &];
        finalLines = Select[lines, StringStartsQ[#, "FINAL "] &];
        stepwiseRows = Quiet@Check[
          ImportString[StringDrop[#, StringLength["STEPWISE "]],
            "RawJSON"] & /@ stepwiseLines,
          $Failed
        ];
        finalRows = Quiet@Check[
          ImportString[StringDrop[#, StringLength["FINAL "]],
            "RawJSON"] & /@ finalLines,
          $Failed
        ];
        rows = If[ListQ[stepwiseRows],
          Select[stepwiseRows, Lookup[#, "Level", None] === 0 &],
          $Failed
        ];
        expected = Lookup[run, "ExpectedMasters", {}];
        If[rows === $Failed || finalRows === $Failed ||
            Length[rows] =!= Length[expected] ||
            Length[finalRows] =!= Length[expected] ||
            Lookup[rows, "Master"] =!= expected ||
            Lookup[finalRows, "Master"] =!= expected,
          Print["Malformed or incomplete Henn FT runner log: ", logFile,
            " expected=", Length[expected],
            " actual=", If[ListQ[rows], Length[rows], $Failed]];
          Exit[1]
        ];
        rows
      ]
    ],
    manifestRuns
  ];
  allOutputs = Join @@ parsedOutputs,
  results = FeynmanTrick`RunIntegrationPipeline /@ plans;
  If[AnyTrue[results,
      !AssociationQ[#] || Lookup[#, "Status", None] =!= "Succeeded" &],
    Print["HENN_FT_BOUNDARY pipeline failed: ",
      InputForm[Select[results,
        !AssociationQ[#] || Lookup[#, "Status", None] =!= "Succeeded" &]]];
    Exit[1]
  ];
  allOutputs = Join @@ Lookup[results, "Outputs"]
];

jsonNumber[value_?NumberQ] := value;
jsonNumber[value_Association] /;
    Sort[Keys[value]] === Sort[{"Re", "Im"}] :=
  value["Re"] + I value["Im"];

outputSeries[row_Association, variable_] := Total[
  (#[[2]] // jsonNumber) variable^#[[1]] & /@ row["Coefficients"]];
targetKey[index_List] := ToString[index, InputForm];

rawSeries = AssociationThread[
  targetKey /@ Lookup[allOutputs, "Master"],
  outputSeries[#, Global`hennEps] & /@ allOutputs
];

missingOutputs = Select[
  scalarTargets, !KeyExistsQ[rawSeries, targetKey[#]] &];
If[missingOutputs =!= {},
  Print["FT result omitted requested Henn scalar integrals: ",
    missingOutputs];
  Exit[1]
];

rawRule = Global`XB[a__] :> Module[{index = {a}},
  (-1)^Total[index] rawSeries[targetKey[index]]
];
canonicalRaw = selectedBasis /. x0InvariantRules /. rawRule;
computedSeries = Map[
  Series[
    Global`hennEps^4 Exp[2 Global`hennEps EulerGamma] #,
    {Global`hennEps, 0, 4}] &,
  canonicalRaw
];
computedBoundary = Table[
  SeriesCoefficient[series, order],
  {series, computedSeries}, {order, 0, 4}
];
seriesMinimumPower[series_SeriesData] := series[[4]]/series[[6]];
seriesMinimumPower[_] := 0;
unexpectedPoleCoefficients = Map[
  Function[series,
    With[{minimum = seriesMinimumPower[series]},
      If[IntegerQ[minimum] && minimum < 0,
        Table[SeriesCoefficient[series, order],
          {order, minimum, -1}],
        {}
      ]
    ]
  ],
  computedSeries
];
unexpectedPoleErrors = Map[
  If[# === {}, 0, Max[Abs[N[#, 30]]]] &,
  unexpectedPoleCoefficients
];
maxUnexpectedPoleError = Max[unexpectedPoleErrors];

errors = Abs[N[computedBoundary - selectedReference, 30]];
errorsByOrder = Max /@ Transpose[errors];
errorsByComponent = Max /@ errors;
maxPaperError = Max[errorsByOrder];
maxError = Max[maxPaperError, maxUnexpectedPoleError];
tolerance = 10^-8;
worstOrdinal = First@First@Position[
  errorsByComponent, Max[errorsByComponent]];

Print["HENN_FT_BOUNDARY errorsByOrder=",
  InputForm[N[errorsByOrder, 12]],
  " maxUnexpectedPoleError=",
  InputForm[N[maxUnexpectedPoleError, 12]],
  " maxError=", InputForm[N[maxError, 12]]];
Print["HENN_FT_BOUNDARY worstComponent=",
  componentSelection[[worstOrdinal]]];

If[!TrueQ[maxError < tolerance],
  Print["HENN_FT_BOUNDARY computedWorst=",
    InputForm[N[computedBoundary[[worstOrdinal]], 16]]];
  Print["HENN_FT_BOUNDARY referenceWorst=",
    InputForm[N[selectedReference[[worstOrdinal]], 16]]];
  Print["HENN_FT_BOUNDARY FAIL: maximum paper-boundary error exceeds ",
    tolerance];
  Exit[1]
];

Print["HENN_FT_BOUNDARY PASS"];
Exit[0];
