(* ::Package:: *)

(* ==========================================================================
   Scripts/dump_transport_checkpoints.m  —  M0 task (15), RewritePlan.md sec 6.

   THE MISSING M4 PARITY HARNESS: runs the OLD (frozen-oracle) DiffExp on a
   configured classic transport line and emits machine-readable values of the
   full solution vector at EVERY segment boundary (matching point) and at the
   final point.  DiffExp2's Transport.m (M4 gate, RewritePlan sec 6 "classic
   parity on bubble/sunrise/2F1/banana lines vs M0 checkpoint dumps, 1e-25")
   is compared against this output with ZERO old-code kernel time.

   Modeled on Scripts/run_ft_stepwise.m (env-var configuration, Print of
   compact RawJSON rows, Exit codes).

   MECHANISM
   - Classic presets call DiffExp`Transport`TransportTo[bcs, line, to, True]
     with SaveExpansions = True (4th positional arg, DiffExp/Transport.m:514).
     The result association carries "SegmentData" (DiffExp/Transport.m:1240),
     a list of 5-element records appended per segment at
     DiffExp/Transport.m:847-867 (and :1063-1080 for the singular-endpoint
     final segment):
       [[1]] CurrLine        - the segment chart: assoc var -> expr(x_local)
       [[2]] lineRelation    - x -> RelateLines[CurrLine, line]
       [[3]] global-parameter interval covered (ascending)
       [[4]] {matchLocal, evalLocal} - the boundary-fix (matching) point and
             the evaluation (hand-off) point IN THE LOCAL CHART COORDINATE
       [[5]] the expansion matrix [integral, epsSlot] of x_local-series,
             already passed through ApplyAnalyticContinuation // Project0s
             at save time (DiffExp/Transport.m:853/865)
   - bubble / sunrise are the FT LEVEL LINES (the M4 "parity examples run
     WITHOUT Mobius where Integrate is involved" case): we reproduce
     Scripts/run_ft_stepwise.m up to each level's TransportLevel call
     (FeynmanTrick/DiffExpIntegration.m:267-525), which itself calls
     DiffExp`Transport`TransportTo[..., True] twice (lines 440-446, 474-480;
     UseMobius -> False, UsePade -> False per :352-353) and returns the two
     full TransportTo results under "LowerResult"/"UpperResult" (:509-510).
   - Per-segment evaluation copies the canonical replay path of
     ToPiecewise (DiffExp/Transport.m:1289-1312): the saved series is
     evaluated by
        (Normal@series)                              [UsePade == False]
        (Project0s[series, GetPade])                 [UsePade == True ]
        /. Logx -> Log[x] /. 0p -> HeavisideTheta[x]
        /. 0m -> HeavisideTheta[-x] /. x -> localPoint
     inside Block[{$MinPrecision = FEWorkingPrecision}] (the SEval2
     convention, DiffExp/Pade.m:70-77).  Substituting the LOCAL point
     directly is equivalent to ToPiecewise's substitute-lineRelation-then-
     evaluate-globally route, because lineRelation(global) == local by
     construction (DiffExp/Transport.m:846).
   - The kinematic Point of a checkpoint is computed EXACTLY from the chart:
     CurrLine /. x -> localPoint.

   OUTPUT (one Print per row, prefix + compact RawJSON):
     HARNESS    {run metadata}
     CONFIG     {Example, Line, pinned old-config snapshot}
     CHECKPOINT {"Example", "Line", "SegmentIndex", "Position":"Start"|"End",
                 "Point":  exact value of the line's variable as a string
                           (single-variable lines) or the full exact
                           kinematic rule list as a string (multi-variable),
                 "KinematicPoint":  {varName -> exact-string},
                 "KinematicPointN": {varName -> 40-digit float},
                 "LocalPoint": exact local chart coordinate as string,
                 "Values": [[ [re,im], ... per eps slot ] ... per integral],
                 "EpsOrders": [0, ..., EpsilonOrder]  (transport eps SLOTS),
                 "EpsPrefactors": per-master FT prefactors (FT presets only:
                           true eps power of slot k of master i is
                           k - EpsPrefactors[i], run_ft_stepwise convention),
                 "NumIntegrals", "NumSegments", "Note" (only when values
                 could not be evaluated, e.g. singular chart center)}
     ERRORS     {Example, Line, ErrorEstimates [integral][epsSlot]} — the old
                library's accumulated two-point error probe
                (DiffExp/Transport.m:1238-1239), when EstimateError != False.
     TRANSPORT  {Example, Line, NumSegments, EndpointIsSingularity,
                 ComputationTime, ...} — per-line summary.
     STEPWISE   (FT presets only) run_ft_stepwise-compatible boundary rows.
     FATAL      loud failure description; the script then Exit[1]s.

   Per checkpoint, "Start" is the segment's boundary-matching point (where
   its boundary data was fixed = previous segment's hand-off) and "End" is
   its evaluation point (the value handed to the next segment; the last
   segment's End is the final point).  Start(i+1) == End(i) numerically — a
   free cross-segment consistency probe for the parity gate.

   NUMERIC HYGIENE: every emitted value is N[..., 40] split into [re, im]
   float pairs.  KNOWN QUIRK: Mathematica RawJSON prints subnormal zeros as
   "0.e-63"-style tokens; downstream parsers must regex-fix
   ([0-9])\.e -> \1.0e before json.loads (same fix compare_stepwise_log.py
   applies to STEPWISE rows).

   CONFIGURATION (environment variables)
     DTC_EXAMPLE          comma list of presets, default "2f1".  Presets:
        2f1               resonant 2F1 (a=1/4,b=1/3,c=2), z: 0 -> 1/2,
                          config EXACTLY Tests/test_resonant_2f1.m:53-92,
                          matrices Tests/Hypergeometric2F1_Resonant_Matrices
        2f1_regular       2F1 (a=1/4,b=1/3,c=3/2), z: 1/100 -> 1/2,
                          config EXACTLY Tests/test_hypergeometric2f1.m:48-95,
                          matrices Tests/Hypergeometric2F1_Matrices
        banana            unequal-mass banana: equal-mass BC chain
                          (Tests/test_unequal_mass_full.m:24-59) then the
                          UnequalMassConfiguration mass line
                          psq=1/2, mm_i = {1+x, 1+x/2, 1+x/3, 1}, x: 0 -> 1
                          (ChopPrecision 500, DivisionOrder 4,
                          ExpansionOrder 70, RadiusOfConvergence 10,
                          UseMobius True, UsePade True, WP 1000 — the M4
                          pinned config), matrices Tests/Banana_Matrices;
                          then (DTC_BANANA_PSQ_LINE=1, default) the psq line
                          1/2 -> 37/10 (canonical Banana_example.m:140-146
                          uses machine 3.7; we pin the EXACT rational 37/10
                          so the oracle endpoint is reproducible).
        banana_equalmass  equal-mass banana t-line -1 -> 10 with
                          DeltaPrescriptions {t-16+I delta}
                          (Reference/Examples/Banana_example.m:15-53),
                          matrices Tests/Banana_EqualMass_Matrices
        bubble | sunrise  FT level lines: full FT iteration as in
                          Scripts/run_ft_stepwise.m, then per level dump the
                          lower (11/23 -> 0) and upper (11/23 -> 1)
                          TransportTo results (both endpoints singular; the
                          singular End row is emitted with Values:null).
     DTC_OUTPUT_DIGITS    output precision for values (default 40)
     DTC_VERBOSITY        DiffExp Verbosity for classic presets (default 1;
                          print-only, does not affect values)
     DTC_FORCE_NORMAL_EVAL  "1": evaluate checkpoints with Normal[] even if
                          the config used Pade (default "0": follow the
                          loaded config, i.e. reproduce the hand-off values
                          the old library actually propagated)
     DTC_EMIT_START_ROWS  "1" (default): emit Start AND End rows per segment;
                          "0": End rows only (plus segment 1 Start)
     DTC_USE_RATIONAL_RECURRENCE  "1" (default, campaign standard) for the
                          banana presets (test_unequal_mass_full.m:121)
     DTC_BANANA_PSQ_LINE  "1" (default): also dump the banana psq line
     DTC_FT_WORKING_PRECISION (500), DTC_FT_EPS_ORDER (0),
     DTC_FT_EXPANSION_ORDER (50), DTC_FT_DIVISION_ORDER (4),
     DTC_FT_BOUNDARY_EXTRA_ORDER (4), DTC_FT_POLE_ALLOWANCE (4)
                          FT-preset knobs, identical to run_ft_stepwise.m

   HARD RULES HONORED: never run this yourself if you are an implementation
   agent — the orchestrator owns the single kernel.  This script never
   deletes anything and never debugs the old core (RewritePlan R8); any
   anomaly is a loud FATAL + Exit[1].
   ========================================================================== *)

(* ---------- bootstrap (parse order matters: packages must be loaded ------
   ---------- BEFORE any code referencing their symbols is parsed) -------- *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
$Path = DeleteDuplicates[Prepend[$Path, repoRoot]];

dtcEnv[name_, default_] := Module[{value = Environment[name]},
  If[StringQ[value] && StringLength[StringTrim[value]] > 0, value, default]
];

dtcRequested = StringTrim /@ StringSplit[dtcEnv["DTC_EXAMPLE", "2f1"], ","];
dtcKnownPresets = {"2f1", "2f1_regular", "banana", "banana_equalmass",
  "bubble", "sunrise"};
If[!And @@ (MemberQ[dtcKnownPresets, #] & /@ dtcRequested),
  Print["FATAL unknown preset(s) in DTC_EXAMPLE=", dtcRequested,
    "; known: ", dtcKnownPresets];
  Exit[1];
];
dtcNeedsFT = Or @@ (MemberQ[{"bubble", "sunrise"}, #] & /@ dtcRequested);

(* Load order copied from Scripts/run_ft_stepwise.m:7-9.  FeynmanTrick is
   loaded ONLY for the FT presets so the classic presets see exactly the
   canonical-test package environment.  All cross-package references below
   are FULLY QUALIFIED (memory: unexported symbols silently no-op
   cross-package calls), so parse-time resolution is deterministic in both
   load modes. *)
If[dtcNeedsFT,
  Quiet[Get["FeynmanTrick/FeynmanTrick.m"], {General::shdw, Symbol::shdw}];
];
Quiet[Get["DiffExp.m"], {General::shdw, Symbol::shdw}];
If[dtcNeedsFT,
  Get[FileNameJoin[{repoRoot, "Scripts", "FTExamples.m"}]];
];

(* ---------- env configuration ---------- *)

dtcOutputDigits = ToExpression[dtcEnv["DTC_OUTPUT_DIGITS", "40"]];
dtcVerbosity = ToExpression[dtcEnv["DTC_VERBOSITY", "1"]];
dtcForceNormalEval = dtcEnv["DTC_FORCE_NORMAL_EVAL", "0"] === "1";
dtcEmitStartRows = dtcEnv["DTC_EMIT_START_ROWS", "1"] === "1";
dtcUseRationalRecurrence =
  dtcEnv["DTC_USE_RATIONAL_RECURRENCE", "1"] === "1";
dtcBananaPsqLine = dtcEnv["DTC_BANANA_PSQ_LINE", "1"] === "1";

dtcFTWorkingPrecision =
  ToExpression[dtcEnv["DTC_FT_WORKING_PRECISION", "500"]];
dtcFTEpsOrder = ToExpression[dtcEnv["DTC_FT_EPS_ORDER", "0"]];
dtcFTExpansionOrder = ToExpression[dtcEnv["DTC_FT_EXPANSION_ORDER", "50"]];
dtcFTDivisionOrder = ToExpression[dtcEnv["DTC_FT_DIVISION_ORDER", "4"]];
dtcFTBoundaryExtraOrder =
  ToExpression[dtcEnv["DTC_FT_BOUNDARY_EXTRA_ORDER", "4"]];

(* FT options exactly as Scripts/run_ft_stepwise.m:30-37 (only evaluated for
   FT presets; the symbols parse harmlessly when FT is not loaded). *)
If[dtcNeedsFT,
  FeynmanTrick`SetFTOption["Threads", 1];
  FeynmanTrick`SetFTOption["FThreads", 1];
  FeynmanTrick`SetFTOption["Verbosity", 0];
  FeynmanTrick`SetFTOption["WorkingPrecision", dtcFTWorkingPrecision];
  FeynmanTrick`SetFTOption["ReductionCache", False];
  FeynmanTrick`SetFTOption["FIRETimeoutSeconds", 1800];
  FeynmanTrick`SetFTOption["IntegrationPoleAllowance",
    ToExpression[dtcEnv["DTC_FT_POLE_ALLOWANCE", "4"]]];
];

(* ---------- generic helpers ---------- *)

ClearAll[dtcFail, dtcJSON, dtcCleanPair, dtcCleanReal, dtcExactString,
  dtcAssertMatrixDir, dtcFindVar, dtcUncompressExpansions, dtcEvalExpansions,
  dtcConfigSnapshot, dtcDumpTransport, dtcCheckTransportResult];

(* Loud failure: print FATAL row, throw out of the current preset. *)
dtcFail[msg__] := (
  Print["FATAL ", msg];
  Throw[$Failed, "DTCFatal"]
);

dtcJSON[assoc_] := ExportString[assoc, "RawJSON", "Compact" -> True];

(* N[..., 40] (DTC_OUTPUT_DIGITS) and split into a simple [re, im] float
   pair.  Keep outputs SIMPLE FLOATS: downstream parsers regex-fix the
   RawJSON "(d).e" -> "(d).0e" quirk (see header).  Non-numeric leaves are
   stringified so the row stays valid JSON (and flagged via "Note"). *)
dtcCleanPair[v_] := Module[{n = N[v, dtcOutputDigits]},
  If[NumericQ[n], {Re[n], Im[n]}, ToString[v // InputForm]]
];
dtcCleanReal[v_] := Module[{n = N[v, dtcOutputDigits]},
  If[NumericQ[n], Re[n], ToString[v // InputForm]]
];

dtcExactString[v_] := ToString[InputForm[v]];

(* Matrix-directory guard: directory must exist and contain the per-variable
   derivative slices d<var>_<k>.m the old loader consumes
   (DiffExp/MatrixLoading.m:36-50 parses variables from these names). *)
dtcAssertMatrixDir[dir_String] := Module[{files},
  If[!DirectoryQ[dir],
    dtcFail["matrix directory missing: ", dir];
  ];
  files = FileNames["d*_*.m", dir];
  If[Length[files] === 0,
    dtcFail["matrix directory has no d<var>_<k>.m slices: ", dir];
  ];
  dir
];

(* Resolve a kinematic variable BY NAME from the loaded system
   (DiffExp`State`FEC[System`Variables], populated by LoadConfiguration from
   the matrix file names).  This sidesteps every symbol-context trap: the
   harness never has to guess which context the parsed variable landed in. *)
dtcFindVar[name_String] := Module[
  {vars = DiffExp`State`FEC[System`Variables], hits},
  If[!ListQ[vars] || Length[vars] === 0,
    dtcFail["no variables loaded (LoadConfiguration not run / failed?)"];
  ];
  hits = Select[vars,
    (Head[#] === Symbol && SymbolName[#] === name) || ToString[#] === name &];
  If[Length[hits] =!= 1,
    dtcFail["variable '", name, "' not uniquely found among ", vars];
  ];
  First[hits]
];

(* SaveExpansionsCompress replay, copied from ToPiecewise
   (DiffExp/Transport.m:1278-1287).  Default config has compression off
   (DiffExp/State.m:124); handle the compressed forms defensively anyway. *)
dtcUncompressExpansions[data_] := Which[
  StringQ[data] && FileExistsQ[data], Uncompress[Import[data]],
  StringQ[data], Uncompress[data],
  True, data
];

(* Evaluate the saved per-segment expansion matrix at a local chart point.
   EXACT copy of the ToPiecewise evaluation chain
   (DiffExp/Transport.m:1293-1304), with the local point substituted
   directly (equivalent, see header), wrapped in the SEval2 $MinPrecision
   convention (DiffExp/Pade.m:70).  usePade -> Project0s[#, GetPade]
   (DiffExp/AnalyticContinuation.m:93, DiffExp/Pade.m:32). *)
dtcEvalExpansions[expansions_List, localPoint_, usePade_] := Module[{evalOne},
  evalOne[s_] := Module[{body},
    body = If[usePade === True,
      DiffExp`AnalyticContinuation`Project\[Theta]s[s, DiffExp`Pade`GetPade],
      Normal[s]
    ];
    body /. DiffExp`Symbols`Logx -> Log[DiffExp`Symbols`x] /.
      DiffExp`Symbols`\[Theta]p -> HeavisideTheta[DiffExp`Symbols`x] /.
      DiffExp`Symbols`\[Theta]m -> HeavisideTheta[-DiffExp`Symbols`x] /.
      DiffExp`Symbols`x -> localPoint
  ];
  Quiet[
    Block[{$MinPrecision = DiffExp`State`FEWorkingPrecision},
      Map[evalOne, expansions, {2}]
    ]
  ]
];

(* Pin the old config inside the dump itself (M4: "old config pinned per
   example").  Accessors per DiffExp/State.m:147-162. *)
dtcConfigSnapshot[] := <|
  "WorkingPrecision" -> DiffExp`State`FEWorkingPrecision,
  "ChopPrecision" -> DiffExp`State`ChopPrecisionVal,
  "ExpansionOrder" ->
    dtcExactString[DiffExp`State`FEC[DiffExp`State`ExpansionOrder]],
  "DivisionOrder" ->
    dtcExactString[DiffExp`State`FEC[DiffExp`State`DivisionOrder]],
  "EpsilonOrder" -> DiffExp`State`EpsilonOrderVal,
  "RadiusOfConvergence" -> dtcExactString[DiffExp`State`RadiusOfConvergenceVal],
  "UseMobius" -> TrueQ[DiffExp`State`FEC[DiffExp`State`UseMobius]],
  "UsePade" -> TrueQ[DiffExp`State`FEC[DiffExp`State`UsePade]],
  "UseRationalRecurrence" ->
    TrueQ[DiffExp`State`FEC[DiffExp`State`UseRationalRecurrence]],
  "SegmentationStrategy" ->
    ToString[DiffExp`State`FEC[DiffExp`State`SegmentationStrategy]],
  "EstimateError" -> ToString[DiffExp`State`FEC["EstimateError"]],
  "DeltaPrescriptions" ->
    dtcExactString[DiffExp`State`FEC[DiffExp`State`DeltaPrescriptions]],
  "Variables" ->
    (ToString /@ Flatten[{DiffExp`State`FEC[System`Variables]}])
|>;

dtcCheckTransportResult[example_, label_, result_] := (
  If[!AssociationQ[result],
    dtcFail[example, "/", label, ": TransportTo did not return an ",
      "association (got head ", Head[result], ")"];
  ];
  If[TrueQ[DiffExp`State`MultivaluedFail],
    dtcFail[example, "/", label, ": old transport aborted with ",
      "MultivaluedFail (analytic continuation failure mid-line); the dump ",
      "would be incomplete"];
  ];
  If[MissingQ[result["SegmentData"]] || !ListQ[result["SegmentData"]],
    dtcFail[example, "/", label, ": result has no SegmentData; ",
      "TransportTo must be called with SaveExpansions == True ",
      "(DiffExp/Transport.m:1240)"];
  ];
  If[Length[result["SegmentData"]] === 0,
    dtcFail[example, "/", label, ": SegmentData is empty"];
  ];
);

(* The central dump: one CHECKPOINT row per segment boundary, ERRORS and
   TRANSPORT summary rows per line.  epsPrefactors is FT-layer metadata
   (Null for classic presets). *)
dtcDumpTransport[example_String, label_String, result_, epsPrefactors_] :=
Module[
  {segs, numInts, epsSlots, usePade, nSegs, numericRowCount = 0,
   emitRow, seg, chartRules, expansions, matchLocal, evalLocal, dims},

  dtcCheckTransportResult[example, label, result];
  segs = result["SegmentData"];
  nSegs = Length[segs];
  numInts = result["NumIntegrals"];
  epsSlots = result["EpsilonOrder"] + 1;
  usePade = !dtcForceNormalEval &&
    TrueQ[DiffExp`State`FEC[DiffExp`State`UsePade]];

  Print["CONFIG ", dtcJSON[Join[
    <|"Example" -> example, "Line" -> label,
      "NumIntegrals" -> numInts, "NumSegments" -> nSegs,
      "PadeEvaluation" -> usePade|>,
    dtcConfigSnapshot[]
  ]]];

  emitRow[segIndex_, position_, localPoint_, chartRules2_, expansions2_] :=
  Module[{vals, cleanVals, kinRules, note = Missing[], row, allNumeric},
    vals = dtcEvalExpansions[expansions2, localPoint, usePade];
    cleanVals = Map[dtcCleanPair, vals, {2}];
    allNumeric = And @@ Flatten[Map[ListQ, cleanVals, {2}]];
    If[!allNumeric,
      note = "values not numeric at this point (singular chart " <>
        "center / divergent endpoint); Values set to null";
      cleanVals = Null;
      ,
      numericRowCount += 1;
    ];
    kinRules = chartRules2 /. DiffExp`Symbols`x -> localPoint;
    kinRules = MapAt[Together, kinRules, {All, 2}];
    row = <|
      "Example" -> example,
      "Line" -> label,
      "SegmentIndex" -> segIndex,
      "Position" -> position,
      "Point" -> If[Length[kinRules] === 1,
        dtcExactString[kinRules[[1, 2]]],
        dtcExactString[kinRules]
      ],
      "KinematicPoint" -> Association[
        (ToString[#[[1]]] -> dtcExactString[#[[2]]]) & /@ kinRules],
      "KinematicPointN" -> Association[
        (ToString[#[[1]]] -> dtcCleanReal[#[[2]]]) & /@ kinRules],
      "LocalPoint" -> dtcExactString[localPoint],
      "Values" -> cleanVals,
      "EpsOrders" -> Range[0, epsSlots - 1],
      "NumIntegrals" -> numInts,
      "NumSegments" -> nSegs
    |>;
    If[ListQ[epsPrefactors],
      row["EpsPrefactors"] = epsPrefactors;
    ];
    If[!MissingQ[note], row["Note"] = note];
    Print["CHECKPOINT ", dtcJSON[row]];
  ];

  Do[
    seg = segs[[ind]];
    If[Length[seg] =!= 5,
      dtcFail[example, "/", label, ": segment record ", ind,
        " does not have 5 elements (DiffExp/Transport.m:847-867 layout)"];
    ];
    chartRules = Normal[seg[[1]]];      (* var -> expr(x_local) *)
    matchLocal = seg[[4, 1]];
    evalLocal = seg[[4, 2]];
    expansions = dtcUncompressExpansions[seg[[5]]];
    If[!ListQ[expansions] || Length[expansions] =!= numInts,
      dtcFail[example, "/", label, ": segment ", ind,
        " expansion matrix has ", Length[expansions],
        " rows, expected NumIntegrals = ", numInts];
    ];
    dims = {Length[expansions], Length[expansions[[1]]]};
    If[dims[[2]] =!= epsSlots,
      Print["WARNING ", example, "/", label, " segment ", ind,
        ": eps slots in data = ", dims[[2]],
        " differ from EpsilonOrder+1 = ", epsSlots, "; using data dims"];
      epsSlots = dims[[2]];
    ];

    (* Start row: the segment's boundary-matching point (segment 1's Start
       is the initial point of the line).  End row: the hand-off point;
       the last segment's End is the final point (or the singular endpoint,
       emitted with Values:null). *)
    If[dtcEmitStartRows || ind === 1,
      emitRow[ind, "Start", matchLocal, chartRules, expansions];
    ];
    emitRow[ind, "End", evalLocal, chartRules, expansions];
    ,
    {ind, nSegs}
  ];

  If[numericRowCount === 0,
    dtcFail[example, "/", label,
      ": no checkpoint produced numeric values — dump is useless"];
  ];

  (* Old library's accumulated error estimates, per integral per eps slot
     (DiffExp/Transport.m:1238-1239). *)
  If[!MissingQ[result["ErrorEstimates"]] && ListQ[result["ErrorEstimates"]],
    Print["ERRORS ", dtcJSON[<|
      "Example" -> example,
      "Line" -> label,
      "ErrorEstimates" -> Map[dtcCleanReal, result["ErrorEstimates"], {2}]
    |>]];
    ,
    Print["WARNING ", example, "/", label,
      ": ErrorEstimates not available (EstimateError == False?)"];
  ];

  Print["TRANSPORT ", dtcJSON[<|
    "Example" -> example,
    "Line" -> label,
    "NumSegments" -> nSegs,
    "EndpointIsSingularity" -> TrueQ[result["EndpointIsSingularity"]],
    "ComputationTime" -> dtcCleanReal[result["ComputationTime"]],
    "EpsilonOrder" -> result["EpsilonOrder"],
    "ExpansionOrder" -> dtcExactString[result["ExpansionOrder"]],
    "FinalKinematicPoint" -> dtcExactString[Normal[result["KinematicPoint"]]]
  |>]];
];

(* ==========================================================================
   CLASSIC PRESETS
   ========================================================================== *)

(* ---- 2f1 (resonant): config copied from Tests/test_resonant_2f1.m:53-92.
   a=1/4, b=1/3, c=2; residue eigenvalues {0,-1} at z=0 (resonant); the
   regular solution is transported from the singular point z=0 to z=1/2. *)
dtcRun2F1[variant_String] := Module[
  {matrixDir, config, zv, bcExpressions, prepared, result, dz,
   hyA, hyB, hyC, label},

  If[variant === "resonant",
    matrixDir = dtcAssertMatrixDir[
      FileNameJoin[{repoRoot, "Tests", "Hypergeometric2F1_Resonant_Matrices"}]
    ] <> "/";
    {hyA, hyB, hyC} = {1/4, 1/3, 2};
    label = "z0_to_half";
    (* Tests/test_resonant_2f1.m:53-64 *)
    config = {
      DiffExp`State`MatrixDirectory -> matrixDir,
      DiffExp`State`Verbosity -> dtcVerbosity,
      (* DEC-18: oracle re-baselined without Mobius; RoC rescaling kept *)
      DiffExp`State`UseMobius -> False,
      DiffExp`State`UsePade -> True,
      DiffExp`State`UseRationalRecurrence -> True,
      System`WorkingPrecision -> 200,
      DiffExp`State`ExpansionOrder -> 60,
      DiffExp`State`DivisionOrder -> 4,
      DiffExp`State`ChopPrecision -> 150,
      DiffExp`State`EpsilonOrder -> 0
    };
    ,
    matrixDir = dtcAssertMatrixDir[
      FileNameJoin[{repoRoot, "Tests", "Hypergeometric2F1_Matrices"}]
    ] <> "/";
    {hyA, hyB, hyC} = {1/4, 1/3, 3/2};
    label = "z001_to_half";
    (* Tests/test_hypergeometric2f1.m:48-58 *)
    config = {
      DiffExp`State`MatrixDirectory -> matrixDir,
      DiffExp`State`Verbosity -> dtcVerbosity,
      (* DEC-18: oracle re-baselined without Mobius; RoC rescaling kept *)
      DiffExp`State`UseMobius -> False,
      DiffExp`State`UsePade -> True,
      System`WorkingPrecision -> 200,
      DiffExp`State`ExpansionOrder -> 60,
      DiffExp`State`DivisionOrder -> 4,
      DiffExp`State`ChopPrecision -> 150,
      DiffExp`State`EpsilonOrder -> 0
    };
  ];

  DiffExp`LoadConfiguration[config];
  If[DiffExp`State`NumIntegrals =!= 2,
    dtcFail["2f1: expected 2 integrals, got ", DiffExp`State`NumIntegrals];
  ];
  zv = dtcFindVar["z"];

  If[variant === "resonant",
    (* BCs select the regular solution at the singular point z=0
       (Tests/test_resonant_2f1.m:78-85): series expansion along z = x. *)
    bcExpressions = {
      Hypergeometric2F1[hyA, hyB, hyC, zv],
      D[Hypergeometric2F1[hyA, hyB, hyC, dz], dz] /. dz -> zv
    };
    prepared = DiffExp`Transport`PrepareBoundaryConditions[
      bcExpressions, Association[zv -> DiffExp`Symbols`x]];
    (* Tests/test_resonant_2f1.m:92 + SaveExpansions for the dump *)
    result = DiffExp`Transport`TransportTo[
      prepared, Association[zv -> 1/2], 1, True];
    ,
    (* numeric BCs near z=0 (Tests/test_hypergeometric2f1.m:67-88) *)
    prepared = {
      Association[zv -> 1/100],
      {
        {N[Hypergeometric2F1[hyA, hyB, hyC, 1/100], 200]},
        {N[D[Hypergeometric2F1[hyA, hyB, hyC, dz], dz] /. dz -> 1/100, 200]}
      }
    };
    (* Tests/test_hypergeometric2f1.m:93 + SaveExpansions *)
    result = DiffExp`Transport`TransportTo[
      prepared, Association[zv -> 1/2], 1, True];
  ];

  dtcDumpTransport[
    If[variant === "resonant", "2f1", "2f1_regular"], label, result, Null];
];

(* ---- equal-mass banana boundary-condition chain, EXACTLY
   Tests/test_unequal_mass_full.m:24-59 (= Reference/Examples/
   Banana_example.m:15-91).  Returns the TransportTo result at t = 1/2. *)
dtcEqualMassChain[] := Module[
  {matrixDir, config, tv, ep, bcs, prepared, tmpAsymptotic, atHalf},

  matrixDir = dtcAssertMatrixDir[
    FileNameJoin[{repoRoot, "Tests", "Banana_EqualMass_Matrices"}]] <> "/";

  (* Tests/test_unequal_mass_full.m:24-30.  NOTE Global`t / Global`\[Delta]:
     the matrix variable is parsed at runtime into Global` (nothing exports
     a symbol named t), matching the canonical top-level test scripts.  The
     config must be loaded BEFORE we can resolve t, so DeltaPrescriptions is
     installed with the variable looked up afterwards via
     UpdateConfiguration (values identical; LoadConfiguration parses the
     variable from the d t_<k>.m file names first). *)
  config = {
    DiffExp`State`MatrixDirectory -> matrixDir,
    DiffExp`State`Verbosity -> dtcVerbosity,
    (* DEC-18: oracle re-baselined without Mobius; RoC rescaling kept *)
      DiffExp`State`UseMobius -> False,
    DiffExp`State`UsePade -> True
  };
  DiffExp`LoadConfiguration[config];
  If[DiffExp`State`NumIntegrals =!= 4,
    dtcFail["banana equal-mass: expected 4 integrals, got ",
      DiffExp`State`NumIntegrals];
  ];
  tv = dtcFindVar["t"];
  DiffExp`UpdateConfiguration[{
    DiffExp`State`DeltaPrescriptions -> {tv - 16 + I*Global`\[Delta]}
  }];

  ep = DiffExp`Symbols`\[Epsilon];
  (* Tests/test_unequal_mass_full.m:36-46 *)
  bcs = {
    "?",
    "?",
    ep (1 + 3 ep) (1 + 4 ep) * (
      -4 E^(3 EulerGamma ep) Gamma[ep]^3/tv +
      6 E^(3 EulerGamma ep) (-1/tv)^(1 + ep) ep Gamma[-ep]^2 Gamma[ep]^3 /
        Gamma[-2 ep] +
      8 E^(3 EulerGamma ep) (-1/tv)^(1 + 2 ep) ep Gamma[-ep]^3 Gamma[ep] *
        Gamma[2 ep]/Gamma[-3 ep] +
      3 E^(3 EulerGamma ep) (-1/tv)^(1 + 3 ep) ep Gamma[-ep]^4 *
        Gamma[3 ep]/Gamma[-4 ep]
    ),
    E^(3 EulerGamma ep) ep^3 Gamma[ep]^3
  };

  prepared = DiffExp`Transport`PrepareBoundaryConditions[
    bcs, Association[tv -> -1/DiffExp`Symbols`x]];

  (* Tests/test_unequal_mass_full.m:52 *)
  DiffExp`UpdateConfiguration[{
    DiffExp`State`DivisionOrder -> 4,
    DiffExp`State`ExpansionOrder -> 70
  }];

  (* Tests/test_unequal_mass_full.m:55-58 *)
  Print["INFO banana: transporting equal-mass BCs to the asymptotic limit"];
  tmpAsymptotic = DiffExp`Transport`TransportTo[
    prepared, Association[tv -> -1/DiffExp`Symbols`x]];
  If[!AssociationQ[tmpAsymptotic],
    dtcFail["banana: asymptotic equal-mass transport failed"]];
  Print["INFO banana: transporting equal-mass BCs to t = 1/2"];
  atHalf = DiffExp`Transport`TransportTo[
    tmpAsymptotic, Association[tv -> 1/2]];
  If[!AssociationQ[atHalf] || !ListQ[atHalf["SeriesValues"]] ||
      Length[atHalf["SeriesValues"]] =!= 4,
    dtcFail["banana: equal-mass boundary values at t = 1/2 malformed"]];
  atHalf
];

(* ---- banana: the M4 pinned parity line.  UnequalMassConfiguration with
   UsePade -> True, UseMobius -> True, RadiusOfConvergence -> 10
   (Tests/test_unequal_mass_full.m:63-73, Reference/Examples/
   Banana_example.m:70-80; RewritePlan M4).  Mass line x: 0 -> 1, then
   (optional) psq line 1/2 -> 37/10. *)
dtcRunBanana[] := Module[
  {atHalf, matrixDir, config, psqv, mm1v, mm2v, mm3v, mm4v, sv, ubcs,
   lineAssoc, massResult, psqLine, psqResult, xs = DiffExp`Symbols`x},

  atHalf = dtcEqualMassChain[];

  matrixDir = dtcAssertMatrixDir[
    FileNameJoin[{repoRoot, "Tests", "Banana_Matrices"}]] <> "/";

  (* Tests/test_unequal_mass_full.m:63-73 *)
  config = {
    DiffExp`State`ChopPrecision -> 500,
    DiffExp`State`DivisionOrder -> 4,
    DiffExp`State`ExpansionOrder -> 70,
    DiffExp`State`MatrixDirectory -> matrixDir,
    DiffExp`State`RadiusOfConvergence -> 10,
    (* DEC-18: oracle re-baselined without Mobius; RoC rescaling kept *)
      DiffExp`State`UseMobius -> False,
    DiffExp`State`UsePade -> True,
    System`WorkingPrecision -> 1000,
    DiffExp`State`Verbosity -> dtcVerbosity
  };
  DiffExp`LoadConfiguration[config];
  (* campaign standard (Tests/test_unequal_mass_full.m:121) *)
  DiffExp`UpdateConfiguration[{
    DiffExp`State`UseRationalRecurrence -> dtcUseRationalRecurrence
  }];
  If[DiffExp`State`NumIntegrals =!= 15,
    dtcFail["banana: expected 15 integrals, got ",
      DiffExp`State`NumIntegrals];
  ];
  {psqv, mm1v, mm2v, mm3v, mm4v} =
    dtcFindVar /@ {"psq", "mm1", "mm2", "mm3", "mm4"};

  sv = atHalf["SeriesValues"];
  (* Master mapping, Tests/test_unequal_mass_full.m:78-97: 6 copies of
     integral 1, 4 of integral 2, 1 of integral 3, 4 of integral 4. *)
  ubcs = {
    Association[psqv -> 1/2, mm1v -> 1, mm2v -> 1, mm3v -> 1, mm4v -> 1],
    sv[[#]] & /@ {1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 4, 4, 4, 4}
  };

  (* The mass line, Tests/test_unequal_mass_full.m:124-128 /
     Banana_example.m:128-132, plus SaveExpansions for the dump. *)
  lineAssoc = Association[
    psqv -> 1/2, mm1v -> 1 + xs, mm2v -> 1 + xs/2,
    mm3v -> 1 + xs/3, mm4v -> 1];
  Print["INFO banana: transporting along the mass line x: 0 -> 1"];
  massResult = DiffExp`Transport`TransportTo[ubcs, lineAssoc, 1, True];
  dtcDumpTransport["banana", "mass_line", massResult, Null];

  If[dtcBananaPsqLine,
    (* Banana_example.m:140-146 with the endpoint pinned to the EXACT
       rational 37/10 (canonical uses machine 3.7; see header). *)
    psqLine = Association[
      psqv -> xs, mm1v -> 2, mm2v -> 3/2, mm3v -> 4/3, mm4v -> 1];
    Print["INFO banana: transporting along the psq line 1/2 -> 37/10"];
    psqResult = DiffExp`Transport`TransportTo[
      massResult, psqLine, 37/10, True];
    dtcDumpTransport["banana", "psq_line", psqResult, Null];
  ];
];

(* ---- banana_equalmass: the canonical equal-mass t-line -1 -> 10 with
   SaveExpansions, exactly Reference/Examples/Banana_example.m:15-53
   (Results2 is the example's own SaveExpansions == True call). *)
dtcRunBananaEqualMass[] := Module[
  {matrixDir, config, tv, ep, bcs, prepared, results1, results2},

  matrixDir = dtcAssertMatrixDir[
    FileNameJoin[{repoRoot, "Tests", "Banana_EqualMass_Matrices"}]] <> "/";

  config = {
    DiffExp`State`MatrixDirectory -> matrixDir,
    DiffExp`State`Verbosity -> dtcVerbosity,
    (* DEC-18: oracle re-baselined without Mobius; RoC rescaling kept *)
      DiffExp`State`UseMobius -> False,
    DiffExp`State`UsePade -> True
  };
  DiffExp`LoadConfiguration[config];
  If[DiffExp`State`NumIntegrals =!= 4,
    dtcFail["banana_equalmass: expected 4 integrals, got ",
      DiffExp`State`NumIntegrals];
  ];
  tv = dtcFindVar["t"];
  DiffExp`UpdateConfiguration[{
    DiffExp`State`DeltaPrescriptions -> {tv - 16 + I*Global`\[Delta]}
  }];

  ep = DiffExp`Symbols`\[Epsilon];
  bcs = {
    "?",
    "?",
    ep (1 + 3 ep) (1 + 4 ep) * (
      -4 E^(3 EulerGamma ep) Gamma[ep]^3/tv +
      6 E^(3 EulerGamma ep) (-1/tv)^(1 + ep) ep Gamma[-ep]^2 Gamma[ep]^3 /
        Gamma[-2 ep] +
      8 E^(3 EulerGamma ep) (-1/tv)^(1 + 2 ep) ep Gamma[-ep]^3 Gamma[ep] *
        Gamma[2 ep]/Gamma[-3 ep] +
      3 E^(3 EulerGamma ep) (-1/tv)^(1 + 3 ep) ep Gamma[-ep]^4 *
        Gamma[3 ep]/Gamma[-4 ep]
    ),
    E^(3 EulerGamma ep) ep^3 Gamma[ep]^3
  };
  prepared = DiffExp`Transport`PrepareBoundaryConditions[
    bcs, Association[tv -> -1/DiffExp`Symbols`x]];

  (* Banana_example.m:47 *)
  Print["INFO banana_equalmass: transporting to t = -1"];
  results1 = DiffExp`Transport`TransportTo[prepared, Association[tv -> -1]];
  If[!AssociationQ[results1],
    dtcFail["banana_equalmass: transport to t = -1 failed"]];

  (* Banana_example.m:52 — the canonical SaveExpansions call *)
  Print["INFO banana_equalmass: transporting along t: -1 -> 10"];
  results2 = DiffExp`Transport`TransportTo[
    results1, Association[tv -> DiffExp`Symbols`x], 10, True];
  dtcDumpTransport["banana_equalmass", "t_line_m1_to_10", results2, Null];
];

(* ==========================================================================
   FT PRESETS (bubble | sunrise): level-line transports, structure copied
   from Scripts/run_ft_stepwise.m:110-245 with checkpoint dumps inserted
   after each TransportLevel call.  TransportLevel pins the old config for
   these lines itself: UseMobius -> False, UsePade -> False, Predivision,
   ChopPrecision = WP - 50 (FeynmanTrick/DiffExpIntegration.m:346-361).
   ========================================================================== *)

(* run_ft_stepwise.m:47-108 helpers, kept verbatim-compatible so the
   STEPWISE rows remain parseable by Scripts/compare_stepwise_log.py. *)
dtcLaurentCoefficient[laur_Association, power_Integer] := Module[
  {idx = power - laur["MinPower"] + 1},
  If[idx >= 1 && idx <= Length[laur["Coefficients"]],
    laur["Coefficients"][[idx]],
    0
  ]
];

dtcCleanNumberStepwise[value_] := Module[{n = N[value, 50], re, im},
  If[!NumericQ[n], Return[value]];
  re = Re[n];
  im = Im[n];
  If[Abs[im] < 10^-80,
    re,
    <|"Re" -> re, "Im" -> im|>
  ]
];
dtcCleanExprStepwise[expr_Association] :=
  AssociationThread[Keys[expr], dtcCleanExprStepwise /@ Values[expr]];
dtcCleanExprStepwise[expr_List] := dtcCleanExprStepwise /@ expr;
dtcCleanExprStepwise[expr_?NumericQ] := dtcCleanNumberStepwise[expr];
dtcCleanExprStepwise[expr_] := expr;

dtcCoeffTable[raw_Association, minPower_Integer, maxPower_Integer] :=
  Table[{p, dtcCleanNumberStepwise[dtcLaurentCoefficient[raw, p]]},
    {p, minPower, maxPower}];

dtcPrintBoundaryRows[example_String, level_Integer, masters_List,
    boundary_Association] := Module[
  {rawValues, prefactors, rawMin, rows, rowMin},
  rawValues = boundary["RawBoundaryValues"];
  prefactors = boundary["EpsPrefactors"];
  rawMin = boundary["RawMinPower"];
  rows = Table[
    rowMin = Lookup[rawValues[[i]], "MinPower", rawMin];
    <|
      "Example" -> example,
      "Level" -> level,
      "Master" -> masters[[i]],
      "EpsPrefactor" -> prefactors[[i]],
      "RawMinPower" -> rowMin,
      "Coefficients" -> dtcCoeffTable[rawValues[[i]], rowMin, 0]
    |>,
    {i, Length[masters]}
  ];
  Do[
    Print["STEPWISE ",
      ExportString[dtcCleanExprStepwise[row], "RawJSON", "Compact" -> True]],
    {row, rows}
  ];
];

dtcRunFTExample[name_String] := Module[
  {
    topology, ftData, outputDir, nLevels, boundaryOrder, deepBoundary,
    currentBCs, currentPrefactors, transportOrder, matrixDir,
    extraSingularFactors, transportResult, levelBoundary, levelIBPBatch,
    requiredOrder
  },
  Print["EXAMPLE ", name];
  FeynmanTrick`SetFTOption["DimensionExpression",
    Global`FTExampleDimension[name]];
  topology = Global`FTExampleTopology[name, "step"];
  If[topology === $Failed,
    dtcFail[name, ": FTExampleTopology failed"]];
  ftData = FeynmanTrick`FeynmanTrickIteration`DefineFTIteration[
    topology, Global`FTExampleSequence[name], {}
  ];
  outputDir = FileNameJoin[{$TemporaryDirectory,
    "DTC_" <> name <> "_" <> ToString[$ProcessID]}];
  If[DirectoryQ[outputDir],
    DeleteDirectory[outputDir, DeleteContents -> True]];
  CreateDirectory[outputDir, CreateIntermediateDirectories -> True];
  ftData = FeynmanTrick`FeynmanTrickIteration`RunFullIteration[
    ftData, outputDir];
  If[ftData === $Failed, dtcFail[name, ": RunFullIteration failed"]];

  nLevels = ftData["NumLevels"];
  boundaryOrder = dtcFTEpsOrder + nLevels + dtcFTBoundaryExtraOrder;
  deepBoundary = FeynmanTrick`BoundaryConditions`DeepestLevelBoundary[
    ftData, boundaryOrder
  ];
  If[!AssociationQ[deepBoundary],
    dtcFail[name, ": DeepestLevelBoundary failed"]];

  dtcPrintBoundaryRows[
    name, nLevels, ftData["Levels"][nLevels]["Masters"],
    <|
      "RawBoundaryValues" -> Map[
        <|
          "MinPower" -> -deepBoundary["EpsPrefactors"][[#]],
          "Coefficients" -> deepBoundary["BoundaryValues"][[#]]
        |> &,
        Range[Length[deepBoundary["BoundaryValues"]]]
      ],
      "EpsPrefactors" -> deepBoundary["EpsPrefactors"],
      "RawMinPower" -> Min[-deepBoundary["EpsPrefactors"]]
    |>
  ];

  currentBCs = deepBoundary["BoundaryValues"];
  currentPrefactors = deepBoundary["EpsPrefactors"];
  Do[
    levelIBPBatch = FeynmanTrick`LevelReduction`PrepareLevelIBPBatch[
      ftData, level];
    If[levelIBPBatch === $Failed,
      dtcFail[name, ": FIRE boundary reduction batch failed at level ", level]];
    (* Carry the full boundary depth, run_ft_stepwise.m:154-172 (the
       c0b24f3 budget lesson: single-level RequiredTransportEpsilonOrder
       underestimates deep towers). *)
    requiredOrder = FeynmanTrick`LevelReduction`RequiredTransportEpsilonOrder[
      ftData, level, dtcFTEpsOrder, currentPrefactors, levelIBPBatch];
    If[requiredOrder === $Failed,
      dtcFail[name, ": stale FIRE boundary reduction batch at level ", level]];
    transportOrder = Length[First[currentBCs]] - 1;
    If[transportOrder < requiredOrder,
      Print["WARNING level ", level, " transport depth ", transportOrder,
        " is below the single-level requirement ", requiredOrder,
        "; increase DTC_FT_BOUNDARY_EXTRA_ORDER."];
    ];
    FeynmanTrick`FeynmanTrickIteration`ExportLevel[
      ftData, level, outputDir, "diffexp", transportOrder
    ];
    matrixDir = FileNameJoin[{outputDir,
      "Level_" <> ToString[level] <> "_Matrices"}];
    dtcAssertMatrixDir[matrixDir];
    extraSingularFactors =
      FeynmanTrick`LevelReduction`CollectLevelIBPSingularFactors[
        ftData, level, levelIBPBatch
      ];
    If[extraSingularFactors === $Failed,
      dtcFail[name, ": stale FIRE boundary reduction batch at level ", level]];
    Print["EXTRA_SINGULAR_FACTORS level ", level, ": ",
      InputForm[extraSingularFactors]];
    transportResult = FeynmanTrick`DiffExpIntegration`TransportLevel[
      matrixDir, currentBCs, transportOrder,
      "WorkingPrecision" -> dtcFTWorkingPrecision,
      "ExpansionOrder" -> dtcFTExpansionOrder,
      "DivisionOrder" -> dtcFTDivisionOrder,
      "Verbosity" -> dtcVerbosity,
      "EpsPrefactors" -> currentPrefactors,
      "ExtraSingularFactors" -> extraSingularFactors,
      "UseRationalRecurrence" -> True
    ];
    If[transportResult === $Failed,
      dtcFail[name, ": TransportLevel failed at level ", level]];

    (* THE CHECKPOINT DUMPS: TransportLevel returns the two raw TransportTo
       results (SaveExpansions == True hardwired,
       FeynmanTrick/DiffExpIntegration.m:440-446 and 474-480). *)
    If[AssociationQ[transportResult["LowerResult"]],
      dtcDumpTransport[name, "level" <> ToString[level] <> "_lower",
        transportResult["LowerResult"], currentPrefactors];
    ];
    If[AssociationQ[transportResult["UpperResult"]],
      dtcDumpTransport[name, "level" <> ToString[level] <> "_upper",
        transportResult["UpperResult"], currentPrefactors];
    ];

    transportResult["BoundaryValuesAbove"] = currentBCs;
    transportResult["EpsPrefactorsAbove"] = currentPrefactors;
    levelBoundary = FeynmanTrick`DiffExpIntegration`ComputeLevelBoundary[
      ftData, level - 1, transportResult, dtcFTEpsOrder, levelIBPBatch
    ];
    If[!AssociationQ[levelBoundary],
      dtcFail[name, ": ComputeLevelBoundary failed at level ", level - 1]];
    dtcPrintBoundaryRows[
      name, level - 1, ftData["Levels"][level - 1]["Masters"],
      levelBoundary
    ];
    currentBCs = levelBoundary["BoundaryValues"];
    currentPrefactors = levelBoundary["EpsPrefactors"];
    ,
    {level, nLevels, 1, -1}
  ];
];

(* ==========================================================================
   DISPATCH
   ========================================================================== *)

Print["HARNESS ", dtcJSON[<|
  "Script" -> "dump_transport_checkpoints",
  "Examples" -> dtcRequested,
  "OutputDigits" -> dtcOutputDigits,
  "EmitStartRows" -> dtcEmitStartRows,
  "ForceNormalEval" -> dtcForceNormalEval,
  "Date" -> DateString["ISODateTime"]
|>]];

dtcRunPreset[name_String] := Switch[name,
  "2f1", dtcRun2F1["resonant"],
  "2f1_regular", dtcRun2F1["regular"],
  "banana", dtcRunBanana[],
  "banana_equalmass", dtcRunBananaEqualMass[],
  "bubble" | "sunrise", dtcRunFTExample[name],
  _, dtcFail["unknown preset ", name]
];

Do[
  Module[{outcome},
    Print["EXAMPLE_BEGIN ", name];
    outcome = CheckAbort[
      Catch[dtcRunPreset[name]; "ok", "DTCFatal"],
      $Failed
    ];
    If[outcome =!= "ok",
      Print["FAILED ", name];
      Exit[1];
    ];
    Print["EXAMPLE_DONE ", name];
  ],
  {name, dtcRequested}
];
Exit[0];
