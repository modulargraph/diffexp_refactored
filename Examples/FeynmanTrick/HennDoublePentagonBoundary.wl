(* Reconstruct the X0 boundary of the 108-integral canonical double-pentagon
   basis from arXiv:1812.11160 with the public Feynman-trick pipeline.

   The paper defines

     XB[a1,...,a11] =
       Integral d^D k1 d^D k2/(i Pi)^D Product_i D_i^(-ai)

   with D_i = q_i^2. DiffExp 2 uses D_i' = -q_i^2, so each raw result is
   multiplied by (-1)^Total[{a1,...,a11}] before the published canonical
   linear combinations are formed. The paper's final normalization is
   eps^4 Exp[2 eps EulerGamma].

   One seven-level two-loop Feynman-trick ladder transports the complete
   108-master L0 basis.  The selected canonical component or all 108
   components are then reconstructed by one exact L0 FIRE reduction and
   compared through epsilon^4.

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
(* Basis discovery must span every observable the universal 108-component
   family promises, independent of a particular --component selection.
   FIRE's generic corner/first-shell scan misses one rank-two dotted master
   needed by this span. *)
basisProbeTargets = DeleteDuplicates[
  Cases[canonicalBasis, Global`XB[a__] :> {a}, Infinity]];

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

(* One complete denominator family owns every subsector.  The Feynman-trick
   recursion transports the complete FIRE master vector at each level; it
   must not be repeated for every scalar term in the canonical definitions.
   D9,D10,D11 are the genuine irreducible-numerator slots.  Negative powers
   of D1,...,D8 are ordinary subsector numerators and are reduced onto the L0
   master basis after that basis has been discovered. *)
universalFamily = Join[familyDefinition, <|
  "Name" -> "henn_double_pentagon_x0",
  "EliminatedPositions" -> {},
  "NumeratorPositions" -> {9, 10, 11},
  "BasisProbeIntegrals" -> basisProbeTargets,
  "CombinationSequence" -> ({1, #} & /@ Range[2, 8])
|>];

workingPrecision = positiveIntegerEnvironment[
  "HENN_FT_WORKING_PRECISION", 80];
(* The maintained paper comparison is an eight-digit gate.  The runner
   independently propagates guarded producer accuracy through the shared
   master-basis ladder and caches the exact per-level profile learned for
   that request. *)
matchingDigits = positiveIntegerEnvironment["HENN_FT_MATCH_DIGITS", 8];
matchingCertificationDigits = positiveIntegerEnvironment[
  "HENN_FT_MATCHING_CERTIFICATION_DIGITS", 6];
expansionOrder = positiveIntegerEnvironment[
  "HENN_FT_EXPANSION_ORDER", 25];
masterEpsilonOrder = positiveIntegerEnvironment[
  "HENN_FT_MASTER_EPSILON_ORDER", 8];
boundaryExtraOrder = positiveIntegerEnvironment[
  "HENN_FT_BOUNDARY_EXTRA_ORDER", 4];
divisionOrder = positiveIntegerEnvironment[
  "HENN_FT_DIVISION_ORDER", 5];
fireTimeoutSeconds = positiveIntegerEnvironment[
  "HENN_FT_FIRE_TIMEOUT_SECONDS", 1800];
l0ReductionBatchSize = positiveIntegerEnvironment[
  "HENN_FT_L0_REDUCTION_BATCH_SIZE", 4];
l0FireThreads = positiveIntegerEnvironment[
  "HENN_FT_L0_FIRE_THREADS", Min[8, Max[1, $ProcessorCount]]];
l0FireBackend = environmentOrDefault[
  "HENN_FT_L0_FIRE_BACKEND", "Classical"];
If[!MemberQ[{"Classical", "Modular"}, l0FireBackend],
  Print["HENN_FT_L0_FIRE_BACKEND must be Classical or Modular; received ",
    InputForm[l0FireBackend]];
  Exit[2]
];
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

pipelineOptions[] := Module[
  {levelCount, levelDeltaPrescriptionSigns},
  levelCount = Length[universalFamily["CombinationSequence"]];
  levelDeltaPrescriptionSigns = Join[
    If[levelCount >= 1, {1}, {}],
    ConstantArray[-1, Max[0, levelCount - 1]]];
  {
  "WorkingPrecision" -> workingPrecision,
  "MatchingDigits" -> matchingDigits,
  "ExpansionOrder" -> expansionOrder,
  (* L0 scalar observables are reconstructed from exact FIRE reductions after
     the single master-basis ladder.  Retain positive master orders so poles
     in those reduction coefficients cannot consume the requested raw
     epsilon^0 target coefficient. *)
  "EpsilonOrder" -> masterEpsilonOrder,
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
  "FIRETimeoutSeconds" -> fireTimeoutSeconds,
  "PreparedCacheDirectory" -> FileNameJoin[{cacheRoot, "fire"}],
  "CheckpointDirectory" ->
    FileNameJoin[{checkpointRoot,
      "henn-double-pentagon-x0-master-basis"}],
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

pipelineOptionList = pipelineOptions[];
If[firePath =!= "",
  pipelineOptionList = Append[pipelineOptionList, "FIREPath" -> firePath]];
If[wolframScript =!= "",
  pipelineOptionList =
    Append[pipelineOptionList, "WolframScript" -> wolframScript]];
plan = FeynmanTrick`PipelinePlan[
  universalFamily, All, Sequence @@ pipelineOptionList];
If[FailureQ[plan],
  Print["Could not construct Henn double-pentagon FT plan: ",
    InputForm[plan]];
  Exit[2]
];
plans = {plan};

Print["HENN_FT_BOUNDARY components=", Length[componentSelection],
  " scalarTargets=", Length[scalarTargets],
  " basisProbes=", Length[basisProbeTargets],
  " numeratorTargets=", Count[scalarTargets, _?(Min[#] < 0 &)],
  " familyGroups=1",
  " masterSelection=All",
  " prescriptionSign=", deltaPrescriptionSign];

If[planOnly,
  Print["HENN_FT_PLAN universalDenominators=", Range[8],
    " scalarObservables=", Length[scalarTargets],
    " recursionLevels=", Length[universalFamily["CombinationSequence"]]];
  Print[InputForm[KeyTake[plan, {
    "FamilyName", "FamilyID", "RequestID", "OutputRequests", "Settings",
    "PreparedCacheDirectory", "CheckpointDirectory"}]]];
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
        "OutputMode" -> groupPlan["Request", "OutputMode"]
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
          finalRows, resolutionLines, resolutionRows, resolution,
          rows, expected},
        logFile = Lookup[run, "LogFile", ""];
        If[!FileExistsQ[logFile],
          Print["Missing Henn FT runner log: ", logFile];
          Exit[1]
        ];
        lines = Import[logFile, "Lines"];
        stepwiseLines = Select[
          lines, StringStartsQ[#, "STEPWISE "] &];
        finalLines = Select[lines, StringStartsQ[#, "FINAL "] &];
        resolutionLines = Select[
          lines, StringStartsQ[#, "OUTPUT_RESOLUTION "] &];
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
        resolutionRows = Quiet@Check[
          ImportString[
            StringDrop[#, StringLength["OUTPUT_RESOLUTION "]],
            "RawJSON"] & /@ resolutionLines,
          $Failed
        ];
        rows = If[ListQ[stepwiseRows],
          Select[stepwiseRows, Lookup[#, "Level", None] === 0 &],
          $Failed
        ];
        If[!ListQ[resolutionRows] || Length[resolutionRows] =!= 1,
          Print["Missing or ambiguous Henn FT All-master resolution: ",
            logFile];
          Exit[1]
        ];
        resolution = First[resolutionRows];
        expected = Lookup[resolution, "Masters", $Failed];
        If[!ListQ[expected] || expected === {} ||
            rows === $Failed || finalRows === $Failed ||
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

masterIntegrals = Lookup[allOutputs, "Master", $Failed];
If[!ListQ[masterIntegrals] || masterIntegrals === {} ||
    !DuplicateFreeQ[masterIntegrals],
  Print["FT result did not publish a unique nonempty L0 master basis"];
  Exit[1]
];
masterSeries = AssociationThread[
  targetKey /@ Lookup[allOutputs, "Master"],
  outputSeries[#, Global`hennEps] & /@ allOutputs
];

(* Reduce all canonical scalar terms onto the one transported L0 master
   basis.  This is one sector-scheduled exact FIRE reduction, not another
   Feynman-trick ladder. *)
prepareScalarReductions[masters_List] := Module[
  {family, topology, setupTopology, detailed, reductions, reportedMasters,
   coefficientRows, residuals, reductionDirectory, cacheDirectory,
   targetBatches, support, complexity, batchContract, batchKey, batchFile,
   batchPayload, batchResults, loadBatch, saveBatch, validBatchQ},
  family = FeynmanTrick`CreateFamily[universalFamily, All];
  If[!AssociationQ[family],
    Return[Failure["HennL0Reduction", <|
      "Detail" -> "could not canonicalize the universal Henn family"|>],
      Module]
  ];
  topology = family["Topology"];
  topology["Name"] = "henn_double_pentagon_x0_l0_observables_" <>
    ToLowerCase[l0FireBackend];
  reductionDirectory = FileNameJoin[{
    cacheRoot, "fire", "henn-double-pentagon-x0-l0-observables-" <>
      ToLowerCase[l0FireBackend]}];
  cacheDirectory = FileNameJoin[{
    cacheRoot, "fire", "henn-double-pentagon-x0-l0-reductions"}];
  FeynmanTrick`SetFTOption[
    "DimensionExpression", universalFamily["Dimension"]];
  FeynmanTrick`SetFTOption["WorkDirectory", reductionDirectory];
  FeynmanTrick`SetFTOption["FIRETimeoutSeconds", fireTimeoutSeconds];
  FeynmanTrick`SetFTOption["FIREBackend", l0FireBackend];
  FeynmanTrick`SetFTOption["Threads", l0FireThreads];
  FeynmanTrick`SetFTOption["FThreads", l0FireThreads];
  If[firePath =!= "",
    FeynmanTrick`SetFTOption["FIREPath", firePath]];
  setupTopology =
    FeynmanTrick`FIREInterface`SetupFIRE[topology, reductionDirectory];
  If[!AssociationQ[setupTopology],
    Return[Failure["HennL0Reduction", <|
      "Detail" -> "FIRE setup failed for the universal L0 family"|>],
      Module]
  ];
  If[setupTopology["NumPropagators"] =!= Length[First[scalarTargets]] ||
      setupTopology["NumeratorPositions"] =!= {9, 10, 11},
    Return[Failure["HennL0Reduction", <|
      "Detail" ->
        "FIRE changed the Henn L0 arity or irreducible-numerator contract",
      "NumPropagators" -> setupTopology["NumPropagators"],
      "NumeratorPositions" -> setupTopology["NumeratorPositions"]|>],
      Module]
  ];
  (* Pin the exact basis already discovered and numerically transported by
     the single ladder.  Otherwise FIRE starts a fresh master-selection
     problem for each observable batch. *)
  setupTopology["Masters"] = masters;
  If[!DirectoryQ[cacheDirectory],
    CreateDirectory[cacheDirectory, CreateIntermediateDirectories -> True]];
  support[index_List] :=
    Flatten@Position[Take[index, 8], _?(# > 0 &)];
  complexity[index_List] := Total[Abs[Min[index, 0]]];
  (* One monolithic 257-target modular reconstruction exceeded both the disk
     and time guards.  Keep the exact FIRE scheduling resumable by L0 sector.
     Lower sectors use small batches, but all top-sector numerators must share
     one invocation: splitting them would repeat the same expensive descent
     through every lower sector for each observable. *)
  targetBatches = Flatten[
    (Partition[
        SortBy[#, complexity],
        UpTo[If[Length[support[First[#]]] === 8,
          Length[#], l0ReductionBatchSize]]] &) /@
      SortBy[GatherBy[scalarTargets, support],
        {Max[complexity /@ #] &, Length}],
    1
  ];
  validBatchQ[payload_, contract_] :=
    AssociationQ[payload] &&
      Lookup[payload, "Schema", None] ===
        "DiffExp2.HennL0ReductionBatch/v1" &&
      Lookup[payload, "Contract", None] === contract &&
      AssociationQ[Lookup[payload, "Reductions", None]] &&
      Sort[Keys[payload["Reductions"]]] ===
        Sort[contract["Targets"]] &&
      ListQ[Lookup[payload, "Masters", None]] &&
      Complement[payload["Masters"], masters] === {};
  loadBatch[file_, contract_] := Module[{payload},
    If[!FileExistsQ[file], Return[$Failed, Module]];
    payload = Quiet@Check[Import[file, "WXF"], $Failed];
    If[validBatchQ[payload, contract], payload, $Failed]
  ];
  saveBatch[file_, payload_, contract_] := Module[
    {temporary, written, reloaded},
    If[!validBatchQ[payload, contract], Return[$Failed, Module]];
    temporary = file <> ".tmp-" <> ToString[$ProcessID];
    If[FileExistsQ[temporary], Quiet[DeleteFile[temporary]]];
    written = Quiet@Check[Export[temporary, payload, "WXF"], $Failed];
    If[written === $Failed || !FileExistsQ[temporary],
      Return[$Failed, Module]];
    reloaded = Quiet@Check[Import[temporary, "WXF"], $Failed];
    If[!validBatchQ[reloaded, contract],
      Quiet[DeleteFile[temporary]];
      Return[$Failed, Module]
    ];
    If[!Quiet@Check[
        RenameFile[temporary, file, OverwriteTarget -> True]; True, False],
      If[FileExistsQ[temporary], Quiet[DeleteFile[temporary]]];
      Return[$Failed, Module]
    ];
    file
  ];
  batchResults = Catch[
    MapIndexed[
      Function[{targets, ordinal},
      batchContract = <|
        "Schema" -> "DiffExp2.HennL0ReductionBatchContract/v1",
        "FamilyID" -> family["FamilyID"],
        "SetupFingerprintRecord" -> setupTopology["SetupFingerprintRecord"],
        "Backend" -> l0FireBackend,
        "TransportedMasters" -> masters,
        "Targets" -> targets|>;
      batchKey = IntegerString[Hash[batchContract, "SHA256"], 16, 64];
      batchFile = FileNameJoin[{cacheDirectory,
        "batch-" <> batchKey <> ".wxf"}];
      batchPayload = loadBatch[batchFile, batchContract];
      If[AssociationQ[batchPayload],
        Print["HENN_FT_L0_REDUCTION CACHE HIT batch=", First[ordinal],
          "/", Length[targetBatches], " targets=", Length[targets]],
        Print["HENN_FT_L0_REDUCTION CACHE MISS batch=", First[ordinal],
          "/", Length[targetBatches], " targets=", Length[targets],
          " support=", support[First[targets]]];
        detailed = FeynmanTrick`FIREInterface`ReduceIntegralsDetailed[
          setupTopology, targets];
        If[!AssociationQ[detailed],
          Throw[Failure["HennL0Reduction", <|
            "Detail" -> "FIRE failed to reduce an L0 observable batch",
            "Batch" -> First[ordinal], "Targets" -> targets|>],
            "HennL0ReductionAbort"]
        ];
        batchPayload = <|
          "Schema" -> "DiffExp2.HennL0ReductionBatch/v1",
          "Contract" -> batchContract,
          "Reductions" -> detailed["Reductions"],
          "Masters" -> DeleteDuplicates[detailed["Masters"]]|>;
        If[saveBatch[batchFile, batchPayload, batchContract] === $Failed,
          Throw[Failure["HennL0Reduction", <|
            "Detail" -> "could not publish an exact L0 reduction batch",
            "Batch" -> First[ordinal], "File" -> batchFile|>],
            "HennL0ReductionAbort"]
        ]
      ];
      batchPayload
      ],
      targetBatches
    ],
    "HennL0ReductionAbort"
  ];
  If[FailureQ[batchResults], Return[batchResults, Module]];
  reductions = Join @@ Lookup[batchResults, "Reductions"];
  reportedMasters =
    DeleteDuplicates[Join @@ Lookup[batchResults, "Masters"]];
  If[Sort[Keys[reductions]] =!= Sort[scalarTargets],
    Return[Failure["HennL0Reduction", <|
      "Detail" -> "batched L0 reductions did not cover every scalar target",
      "MissingTargets" -> Complement[scalarTargets, Keys[reductions]]|>],
      Module]
  ];
  If[Complement[reportedMasters, masters] =!= {},
    Return[Failure["HennL0Reduction", <|
      "Detail" ->
        "observable reduction reported a master outside the transported basis",
      "UnexpectedMasters" -> Complement[reportedMasters, masters]|>],
      Module]
  ];
  coefficientRows = AssociationMap[
    Function[target,
      Table[
        Together[Coefficient[
          reductions[target], Global`G[1, masters[[j]]]]],
        {j, Length[masters]}]
    ],
    scalarTargets
  ];
  residuals = AssociationMap[
    Function[target, Together[
      reductions[target] -
        Sum[
          coefficientRows[target][[j]]*
            Global`G[1, masters[[j]]],
          {j, Length[masters]}]
    ]],
    scalarTargets
  ];
  If[AnyTrue[Values[residuals], !TrueQ[PossibleZeroQ[#]] &],
    Return[Failure["HennL0Reduction", <|
      "Detail" ->
        "observable reduction was not linear in the transported master basis",
      "Residuals" ->
        Select[residuals, !TrueQ[PossibleZeroQ[#]] &]|>],
      Module]
  ];
  coefficientRows
];

reductionCoefficients = prepareScalarReductions[masterIntegrals];
If[FailureQ[reductionCoefficients],
  Print["HENN_FT_L0_REDUCTION FAIL ", InputForm[reductionCoefficients]];
  Exit[1]
];

rationalMinimumPower[coefficient_] := Module[{expanded, minimum},
  If[TrueQ[PossibleZeroQ[coefficient]], Return[Infinity, Module]];
  expanded = Quiet@Check[
    Series[
      coefficient /.
        Global`d -> 4 - 2 Global`hennEps /.
        FeynmanTrick`FTeps -> Global`hennEps,
      {Global`hennEps, 0, 0}],
    $Failed
  ];
  If[expanded === $Failed, Return[$Failed, Module]];
  minimum = If[Head[expanded] === SeriesData,
    expanded[[4]]/expanded[[6]], 0];
  If[IntegerQ[minimum], minimum, $Failed]
];
reductionMinimumPowers =
  rationalMinimumPower /@ Flatten[Values[reductionCoefficients]];
If[MemberQ[reductionMinimumPowers, $Failed],
  Print["Could not determine the epsilon pole depth of the L0 reductions"];
  Exit[1]
];
finiteReductionMinimumPowers =
  DeleteCases[reductionMinimumPowers, Infinity];
requiredMasterEpsilonOrder = If[finiteReductionMinimumPowers === {},
  0, Max[0, -Min[finiteReductionMinimumPowers]]];
If[requiredMasterEpsilonOrder > masterEpsilonOrder,
  Print["HENN_FT_MASTER_EPSILON_ORDER=", masterEpsilonOrder,
    " is insufficient for the exact L0 reduction; rerun with at least ",
    requiredMasterEpsilonOrder];
  Exit[1]
];
Print["HENN_FT_L0_REDUCTION masters=", Length[masterIntegrals],
  " observables=", Length[scalarTargets],
  " requiredMasterEpsilonOrder=", requiredMasterEpsilonOrder];

epsilonReductionCoefficients = AssociationMap[
  Function[row,
    Together[# /. Global`d -> 4 - 2 Global`hennEps /.
        FeynmanTrick`FTeps -> Global`hennEps] & /@ row],
  reductionCoefficients
];
rawSeries = AssociationThread[
  targetKey /@ scalarTargets,
  Map[
    Function[target,
      Series[
        Sum[
          epsilonReductionCoefficients[target][[j]]*
            masterSeries[targetKey[masterIntegrals[[j]]]],
          {j, Length[masterIntegrals]}],
        {Global`hennEps, 0, 0}]
    ],
    scalarTargets
  ]
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
  Coefficient[Normal[series], Global`hennEps, order],
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
