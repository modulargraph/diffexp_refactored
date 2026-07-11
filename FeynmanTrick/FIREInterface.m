(* ::Package:: *)
(* FIREInterface - Wrapper for FIRE6 IBP reduction *)
(* Handles topology definition, FIRE setup, basis finding, *)
(* integral reduction, and differential matrix computation *)
(*
   ARCHITECTURE: Multi-topology batch processing

   Instead of clearing and restarting FIRE for each topology, we:
   1. Generate all .start files first (SetupFIRE or SetupFIREBatch)
   2. Load all topologies with unique problem numbers
   3. Run FIRE6 once with all problems in a single config
   4. Parse combined results

   This avoids fermat subprocess race conditions and is how FIRE6 is
   designed to work with multiple related topologies.
*)

BeginPackage["FeynmanTrick`FIREInterface`", {"FeynmanTrick`", "FeynmanTrick`PropagatorAlgebra`"}];

DefineTopology::usage =
  "DefineTopology[name, loopMomenta, externalMomenta, propagators, replacements] \
defines a topology for FIRE reduction. Propagator convention: D_j = -q_j^2 + m_j^2. \
Returns a Topology association.";

SetupFIRE::usage =
  "SetupFIRE[topology, workDir] generates .start file for FIRE6. \
Does NOT run FIRE - call FindBasis or ReduceIntegrals after setup.";

SetupFIREBatch::usage =
  "SetupFIREBatch[{topo1, topo2, ...}, workDir] sets up multiple topologies \
at once, assigning each a unique problem number. Returns list of updated topologies.";

FindBasis::usage =
  "FindBasis[topology] runs FIRE6 to determine master integrals. \
Returns the topology with Masters field populated.";

FindBasisBatch::usage =
  "FindBasisBatch[{topo1, topo2, ...}] finds master integrals for all topologies \
in a single FIRE6 run. Returns list of updated topologies.";

ReduceIntegrals::usage =
  "ReduceIntegrals[topology, integrals] reduces a list of integrals \
(given as index vectors) to the master basis using FIRE6. \
Returns an association: integral -> linear combination of masters.";

ReduceIntegralsDetailed::usage =
  "ReduceIntegralsDetailed[topology, integrals] reduces integrals with FIRE6 and \
returns an association containing reductions and the masters reported by FIRE.";

ReduceIntegralsBatch::usage =
  "ReduceIntegralsBatch[{{topo1, ints1}, {topo2, ints2}, ...}] reduces \
integrals for multiple topologies in a single FIRE6 run.";

ComputeDiffMatrix::usage =
  "ComputeDiffMatrix[topology, variable] computes the differential equation \
matrix A such that d/dx Masters = A . Masters. \
The variable is typically the Feynman parameter xx.";

GetFIREResult::usage =
  "GetFIREResult[topology, integral] returns the reduction of a single integral.";

TopologyQ::usage = "TopologyQ[expr] returns True if expr is a valid topology association.";

ClearFIREState::usage = "ClearFIREState[] clears all FIRE internal state. Call before setting up new topologies.";

Begin["`Private`"];

(* Global state: track which topologies have been set up *)
$SetupTopologies = <||>;  (* name -> topology association *)
$NextProblemNumber = 1;
$ReductionCache = <||>;

(* ============================================================ *)
(* DefineTopology                                                *)
(* ============================================================ *)

DefineTopology[name_String, loopMomenta_List, externalMomenta_List,
               propagators_List, replacements_List:{}] :=
  <|
    "Name" -> name,
    "LoopMomenta" -> loopMomenta,
    "ExternalMomenta" -> externalMomenta,
    "Propagators" -> propagators,
    "Replacements" -> replacements,
    "NumPropagators" -> Length[propagators],
    "EliminatedPositions" -> {},
    "WorkDirectory" -> "",
    "ProblemNumber" -> 0,  (* Assigned during SetupFIRE *)
    "Masters" -> {},
    "MasterRules" -> {},
    "StartFileReady" -> False,
    "Variables" -> extractVariables[propagators, replacements, loopMomenta, externalMomenta]
  |>;

TopologyQ[t_Association] := KeyExistsQ[t, "Propagators"] && KeyExistsQ[t, "LoopMomenta"];
TopologyQ[_] := False;

topologyEliminatedPositions[topology_Association] :=
Module[{nP, eliminated},
  nP = topology["NumPropagators"];
  eliminated = Lookup[topology, "EliminatedPositions", {}];
  Sort[DeleteDuplicates[
    Select[eliminated, IntegerQ[#] && 1 <= # <= nP &]
  ]]
];

feynmanTrickRestrictions[topology_Association] :=
Module[{nP, eliminated},
  nP = topology["NumPropagators"];
  eliminated = topologyEliminatedPositions[topology];
  Table[
    ReplacePart[ConstantArray[0, nP], pos -> 1],
    {pos, eliminated}
  ]
];

allowedBasisSectorQ[sector_List, eliminated_List] :=
  AllTrue[eliminated, sector[[#]] === 0 &];

basisSectors[topology_Association] :=
Module[{nP, eliminated},
  nP = topology["NumPropagators"];
  eliminated = topologyEliminatedPositions[topology];
  Select[Rest[Tuples[{0, 1}, nP]], allowedBasisSectorQ[#, eliminated] &]
];

coefficientMatrixForPropagators[props_List, squares_List] :=
  Table[
    Coefficient[props[[i]], squares[[j]]],
    {i, Length[props]}, {j, Length[squares]}
  ];

completePropagators[props_List, squares_List] := Module[
  {completed = props, currentRank, targetRank, candidates, improved},
  targetRank = Length[squares];
  currentRank = MatrixRank[coefficientMatrixForPropagators[completed, squares]];
  candidates = DeleteCases[squares, _?(MemberQ[completed, #] &)];

  While[currentRank < targetRank,
    improved = SelectFirst[
      candidates,
      MatrixRank[coefficientMatrixForPropagators[Append[completed, #], squares]] >
        currentRank &,
      Missing["NotFound"]
    ];
    If[MissingQ[improved], Break[]];
    AppendTo[completed, improved];
    candidates = DeleteCases[candidates, improved];
    currentRank = MatrixRank[coefficientMatrixForPropagators[completed, squares]];
  ];

  completed
];

normalizeIntegralIndex[topology_Association, integral_List] := Module[
  {nP = topology["NumPropagators"]},
  Which[
    Length[integral] == nP,
      integral,
    Length[integral] < nP,
      PadRight[integral, nP, 0],
    True,
      Print["Error: integral index has length ", Length[integral],
        " but topology ", topology["Name"], " has ", nP, " propagators."];
      $Failed
  ]
];

normalizeIntegralIndices[topology_Association, integrals_List] :=
  normalizeIntegralIndex[topology, #] & /@ integrals;

$reductionCacheSchema = "FeynmanTrick.ReductionCache/v2";
$fireSetupFingerprintSchema = "FeynmanTrick.FIRESetup/v1";

fileSHA256[path_String] := If[FileExistsQ[path],
  IntegerString[FileHash[path, "SHA256"], 16, 64], Missing["NotFound"]];

currentFIRERuntimeFingerprintRecord[] := Module[{firePath},
  firePath = ExpandFileName[FeynmanTrick`Private`$FTConfig["FIREPath"]];
  <|
    "FIREPath" -> firePath,
    "FIREBinarySHA256" -> fileSHA256[
      FileNameJoin[{firePath, "bin", "FIRE6"}]],
    "FIRESourceSHA256" -> fileSHA256[
      FileNameJoin[{firePath, "FIRE6.m"}]],
    "SystemID" -> $SystemID,
    "DimensionVariable" ->
      FeynmanTrick`Private`$FTConfig["DimensionVariable"]
  |>
];

fallbackSetupFingerprintRecord[topology_Association] := <|
  "Schema" -> $fireSetupFingerprintSchema,
  "VerifiedStartFile" -> False,
  (* Unverified/test topology objects stay scoped to their exact setup
     instance as well as their mathematical content. *)
  "LegacyScope" -> Lookup[topology, {
    "WorkDirectory", "Name", "ProblemNumber"}, Missing["Absent"]],
  "Topology" -> Lookup[topology, {
    "LoopMomenta", "ExternalMomenta", "Propagators", "Replacements",
    "NumPropagators", "OriginalNumPropagators", "EliminatedPositions",
    "NumeratorPositions"}, Missing["Absent"]],
  "AutoDetectRestrictions" -> Lookup[
    FeynmanTrick`Private`$FTConfig, "AutoDetectRestrictions",
    Missing["Unset"]],
  "DimensionVariable" -> Lookup[
    FeynmanTrick`Private`$FTConfig, "DimensionVariable", Missing["Unset"]]
|>;

reductionSetupFingerprintRecord[topology_Association] :=
  Lookup[topology, "SetupFingerprintRecord",
    fallbackSetupFingerprintRecord[topology]];

cachedTopologySemanticCompatibleQ[topology_Association] := Module[
  {stored, storedDimension, currentDimension},
  stored = Lookup[topology, "SetupFingerprintRecord", Missing["Absent"]];
  (* Unverified topology keys already embed the current dimension through the
     fallback setup record, so a changed dimension cannot produce a hit. *)
  If[!AssociationQ[stored], Return[True, Module]];
  storedDimension = Lookup[stored, "DimensionVariable", Missing["Absent"]];
  currentDimension =
    FeynmanTrick`Private`$FTConfig["DimensionVariable"];
  If[storedDimension =!= currentDimension,
    Print["Error: cached FIRE reduction dimension variable ",
      storedDimension, " does not match current ", currentDimension,
      ". Restore the setup configuration or rerun SetupFIRE."];
    Return[False, Module]];
  True
];

preparedTopologyCompatibleQ[topology_Association] := Module[
  {stored, current, startFile, keys, mismatches},
  stored = Lookup[topology, "SetupFingerprintRecord", Missing["Absent"]];
  If[!AssociationQ[stored],
    Print["Error: prepared topology has no setup fingerprint; rerun SetupFIRE."];
    Return[False, Module]];
  startFile = FileNameJoin[{topology["WorkDirectory"],
    topology["Name"] <> ".start"}];
  current = Join[currentFIRERuntimeFingerprintRecord[], <|
    "StartFileSHA256" -> fileSHA256[startFile]|>];
  keys = {"FIREPath", "FIREBinarySHA256", "FIRESourceSHA256", "SystemID",
    "DimensionVariable", "StartFileSHA256"};
  mismatches = Select[keys,
    Lookup[stored, #, Missing["StoredAbsent"]] =!=
      Lookup[current, #, Missing["CurrentAbsent"]] &];
  If[!TrueQ[Lookup[stored, "VerifiedStartFile", False]] ||
      !FileExistsQ[startFile],
    mismatches = DeleteDuplicates[Append[mismatches, "StartFile"]]];
  If[mismatches =!= {},
    Print["Error: prepared FIRE setup is incompatible with current runtime: ",
      mismatches, ". Rerun SetupFIRE."];
    Return[False, Module]];
  True
];

reductionCacheKey[topology_Association, fireIntegral_List] := {
  $reductionCacheSchema,
  reductionSetupFingerprintRecord[topology],
  fireIntegral
};

(* The exact setup record is part of the key, not merely its digest.  This
   keeps cache hits deterministic and collision-safe while allowing identical
   prepared content to reuse exact reductions even if its transient FIRE
   problem number differs. *)

cacheReduction[topology_Association, fireIntegral_List, reduction_, masters_List] :=
Module[{key},
  If[!TrueQ[Lookup[FeynmanTrick`Private`$FTConfig, "ReductionCache", True]],
    Return[Null, Module]
  ];
  key = reductionCacheKey[topology, fireIntegral];
  $ReductionCache[key] = <|
    "Reduction" -> reduction,
    "Masters" -> masters
  |>;
];

cachedReduction[topology_Association, fireIntegral_List] :=
Module[{key},
  If[!TrueQ[Lookup[FeynmanTrick`Private`$FTConfig, "ReductionCache", True]],
    Return[Missing["Disabled"], Module]
  ];
  key = reductionCacheKey[topology, fireIntegral];
  Lookup[$ReductionCache, Key[key], Missing["NotCached"]]
];


(* Extract symbolic variables from propagators and replacements *)
(* Returns Global` symbols suitable for FIRE/Fermat *)
extractVariables[props_List, repls_List, loops_List, exts_List] :=
Module[{allSyms, momenta, vars},
  allSyms = Cases[{props, repls}, _Symbol, Infinity] // DeleteDuplicates;
  momenta = Join[loops, exts];
  vars = Complement[allSyms, momenta];
  (* Return as Global symbols for FIRE compatibility *)
  Symbol[SymbolName[#]] & /@ vars
];

(* Build substitution rules mapping contexted symbols to Global equivalents *)
(* This ensures Fermat (which FIRE uses internally) gets simple variable names *)
buildFIRESubstitution[topology_Association] :=
Module[{allSyms, momenta, contextedSyms},
  allSyms = Cases[topology["Propagators"], _Symbol, Infinity] // DeleteDuplicates;
  momenta = Join[topology["LoopMomenta"], topology["ExternalMomenta"]];
  (* Find symbols not in Global` context that aren't momenta *)
  contextedSyms = Select[
    Complement[allSyms, momenta],
    (Context[#] =!= "Global`") &
  ];
  (* Map each to its Global equivalent *)
  Rule[#, Symbol[SymbolName[#]]] & /@ contextedSyms
];


(* ============================================================ *)
(* ClearFIREState - reset all FIRE state                        *)
(* ============================================================ *)

ClearFIREState[] := Module[{},
  $SetupTopologies = <||>;
  $NextProblemNumber = 1;
  $ReductionCache = <||>;

  (* Load FIRE6.m if needed *)
  ensureFIRELoaded[];

  (* Clear FIRE internal state *)
  Quiet[
    Unprotect[FIRE`Internal, FIRE`External, FIRE`Propagators, FIRE`Replacements,
              FIRE`RESTRICTIONS];
    Clear[FIRE`Internal, FIRE`External, FIRE`Propagators, FIRE`Replacements,
          FIRE`RESTRICTIONS];
    (* Clear all problem-indexed data *)
    Clear[FIRE`ExampleDimension, FIRE`SBasis0L, FIRE`SBasis0D, FIRE`SBasis0C,
          FIRE`SBasisL, FIRE`SBasisS, FIRE`SBasisR, FIRE`SBasisRL, FIRE`SBasisM, FIRE`HPI];
    FIRE`Burning = False;
  , {Unprotect::ssym, Clear::ssym}];
];


(* ============================================================ *)
(* ensureFIRELoaded - load FIRE6.m if not already               *)
(* ============================================================ *)

ensureFIRELoaded[] := Module[{firePath},
  (* Check if FIRE is actually loaded by testing for a key function, not just $Packages *)
  (* $Packages can be misleading if the context was created without loading the package *)
  If[!ValueQ[FIRE`SaveStart] || Head[FIRE`SaveStart] === Symbol,
    firePath = FeynmanTrick`Private`$FTConfig["FIREPath"];
    If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 2,
      Print["Loading FIRE6.m from: ", firePath];
    ];
    Block[{$Path = Prepend[$Path, firePath]},
      Get["FIRE6.m"];
    ];
    (* Verify FIRE was actually loaded *)
    If[!ValueQ[FIRE`SaveStart] || Head[FIRE`SaveStart] === Symbol,
      Print["Error: Failed to load FIRE6.m. Check FIREPath: ", firePath];
      Return[$Failed];
    ];
  ];
];


(* ============================================================ *)
(* FIRE runner helpers                                           *)
(* ============================================================ *)

shellQuote[s_String] := "\"" <> StringReplace[s, {
    "\\" -> "\\\\",
    "\"" -> "\\\""
  }] <> "\"";

safeReadString[path_String] := Module[{txt},
  If[!FileExistsQ[path], Return["", Module]];
  txt = Quiet[ReadString[path], {ReadString::read, General::stop}];
  If[StringQ[txt], txt, ""]
];

cleanupFIREProcesses[fireBin_String, configName_String] :=
Module[{fireDir, flameBin, fermatBin},
  fireDir = DirectoryName[fireBin];
  flameBin = FileNameJoin[{fireDir, "FLAME6"}];
  fermatBin = FileNameJoin[{fireDir, "..", "extra", "fuel", "extra", "ferm64", "fer64"}];

  (* Kill workers tied to this config first, then any orphan Fermat from this FIRE install. *)
  Quiet[Run["pkill -9 -f " <> shellQuote[fireBin <> " -c " <> configName] <> " 2>/dev/null || true"]];
  Quiet[Run["pkill -9 -f " <> shellQuote[flameBin <> " -c " <> configName] <> " 2>/dev/null || true"]];
  Quiet[Run["pkill -9 -f " <> shellQuote[fermatBin] <> " 2>/dev/null || true"]];
];

cleanupFIREWorkDir[dir_String] := Module[{tempDir},
  tempDir = FileNameJoin[{dir, "temp"}];
  If[DirectoryQ[tempDir],
    Quiet[DeleteDirectory[tempDir, DeleteContents -> True]];
  ];
  (* Fallback for stubborn temp dirs left by aborted workers. *)
  If[DirectoryQ[tempDir],
    Quiet[Run["rm -rf " <> shellQuote[tempDir] <> " 2>/dev/null || true"]];
  ];

  (* Remove stale per-attempt run directories. *)
  Scan[
    Quiet[DeleteDirectory[#, DeleteContents -> True]] &,
    Select[FileNames["run_*", dir], DirectoryQ]
  ];
];


(* ============================================================ *)
(* runFIRE6 - helper to run FIRE6 binary robustly               *)
(* ============================================================ *)

runFIRE6[fireBin_String, dir_String, configName_String] :=
Module[{exitCode, maxRetries = 3, attempt},
  Do[
    cleanupFIREProcesses[fireBin, configName];
    cleanupFIREWorkDir[dir];
    exitCode = runFIRE6Once[fireBin, dir, configName];
    If[exitCode === 0, Return[0, Module]];
    If[attempt < maxRetries,
      If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
        Print["FIRE6 failed (exit ", exitCode, "), retrying (attempt ", attempt+1, "/", maxRetries, ")..."];
      ];
      (* Additional cleanup before retry *)
      cleanupFIREProcesses[fireBin, configName];
      cleanupFIREWorkDir[dir];
      Pause[2];
    ];
  , {attempt, maxRetries}];
  If[exitCode =!= 0,
    Print["FIRE6 failed after ", maxRetries, " attempts (last exit code: ", exitCode, ")."];
  ];
  exitCode
];

runFIRE6Once[fireBin_String, dir_String, configName_String] :=
Module[{exitCode, result, logFile, cmd, stdoutFile, stderrFile,
        timeoutSeconds, fireBinQ, configQ, stdoutQ, stderrQ,
        runDir, inputFiles, tableFiles},
  logFile = FileNameJoin[{dir, "fire_stdout.log"}];
  stdoutFile = FileNameJoin[{dir, "fire_stdout.tmp"}];
  stderrFile = FileNameJoin[{dir, "fire_stderr.tmp"}];

  timeoutSeconds = Lookup[FeynmanTrick`Private`$FTConfig, "FIRETimeoutSeconds", 600];
  If[!IntegerQ[timeoutSeconds] || timeoutSeconds <= 0,
    timeoutSeconds = 600;
  ];

  fireBinQ = shellQuote[fireBin];
  configQ = shellQuote[configName];
  stdoutQ = shellQuote[stdoutFile];
  stderrQ = shellQuote[stderrFile];

  (* Isolate each run in its own working directory to avoid stale temp/db locks
     from interrupted FIRE workers. *)
  runDir = FileNameJoin[{dir, "run_" <> StringReplace[CreateUUID[], "-" -> ""]}];
  CreateDirectory[runDir, CreateIntermediateDirectories -> True];
  inputFiles = Join[
    FileNames["*.start", dir],
    FileNames["*.config", dir],
    FileNames["*.m", dir]
  ] // DeleteDuplicates;
  Scan[
    Quiet[CopyFile[#, FileNameJoin[{runDir, FileNameTake[#]}], OverwriteTarget -> True]] &,
    inputFiles
  ];

  cmd = StringJoin[
    (* Keep the FIRE binary itself as $!, rather than backgrounding the whole
       `cd && FIRE` compound command and monitoring an extra subshell. *)
    "cd ", shellQuote[runDir], " || exit 125; ",
    fireBinQ, " -c ", configQ, " > ", stdoutQ, " 2> ", stderrQ, " & ",
    "fire_pid=$!; ",
    "start=$(date +%s); timeout=", ToString[timeoutSeconds], "; ",
    (* A two-second poll made every successful short reduction take at least
       two seconds.  Poll quickly while startup-sized FIRE jobs normally
       finish, then back off so long reductions do not wake the shell often. *)
    "poll_interval=0.05; poll_count=0; ",
    "while kill -0 $fire_pid 2>/dev/null; do ",
      "now=$(date +%s); ",
      "if [ $((now-start)) -ge $timeout ]; then ",
        (* Capture the complete worker tree before its parent is terminated
           and descendants are reparented. Descendants are signalled first,
           giving FIRE a short chance to reap them before its own TERM. *)
        "collect_descendants() ( ",
          "for descendant in $(pgrep -P \"$1\" 2>/dev/null); do ",
            "collect_descendants \"$descendant\"; echo \"$descendant\"; ",
          "done; ",
        "); ",
        "child_pids=$(collect_descendants \"$fire_pid\"); ",
        "if [ -n \"$child_pids\" ]; then kill -TERM $child_pids 2>/dev/null || true; fi; ",
        "sleep 0.1; ",
        "kill -TERM $fire_pid 2>/dev/null || true; ",
        "sleep 2; ",
        "if [ -n \"$child_pids\" ]; then kill -KILL $child_pids 2>/dev/null || true; fi; ",
        "kill -KILL $fire_pid 2>/dev/null || true; ",
        (* Reap the monitored parent before reporting the timeout. *)
        "wait $fire_pid 2>/dev/null || true; ",
        "echo ", shellQuote["FIRE6 timeout after " <> ToString[timeoutSeconds] <> "s"], " >> ", stderrQ, "; ",
        "exit 124; ",
      "fi; ",
      "sleep $poll_interval; ",
      "poll_count=$((poll_count+1)); ",
      "if [ $poll_count -eq 40 ]; then poll_interval=0.25; fi; ",
    "done; ",
    "wait $fire_pid; ",
    "exit $?"
  ];

  exitCode = Run[cmd];
  (* Run returns exit code * 256 on Unix *)
  exitCode = If[IntegerQ[exitCode], BitShiftRight[exitCode, 8], exitCode];

  (* Copy generated table outputs back to the canonical topology directory. *)
  tableFiles = FileNames["*.tables", runDir];
  Scan[
    Quiet[CopyFile[#, FileNameJoin[{dir, FileNameTake[#]}], OverwriteTarget -> True]] &,
    tableFiles
  ];

  (* Read output *)
  result = <|
    "ExitCode" -> exitCode,
    "StandardOutput" -> safeReadString[stdoutFile],
    "StandardError" -> safeReadString[stderrFile]
  |>;

  (* Clean up temp files *)
  Quiet[DeleteFile[stdoutFile]];
  Quiet[DeleteFile[stderrFile]];

  (* Save output to log file *)
  Export[logFile, result["StandardOutput"] <> result["StandardError"], "Text"];

  exitCode = result["ExitCode"];

  (* Ensure no orphan workers leak across retries/runs. *)
  If[exitCode =!= 0,
    cleanupFIREProcesses[fireBin, configName];
  ];

  (* Best-effort cleanup of this attempt's isolated run directory. *)
  Quiet[DeleteDirectory[runDir, DeleteContents -> True]];

  (* Print output on verbose or failure *)
  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 2 || exitCode =!= 0,
    If[StringLength[result["StandardOutput"]] > 0,
      Print["FIRE6 output (tail): "];
      Print[StringTake[result["StandardOutput"], -Min[500, StringLength[result["StandardOutput"]]]]];
    ];
    If[exitCode =!= 0,
      Print["FIRE6 exit code: ", exitCode];
      If[StringLength[result["StandardError"]] > 0,
        Print["FIRE6 stderr: ", StringTake[result["StandardError"], -Min[300, StringLength[result["StandardError"]]]]];
      ];
    ];
  ];

  exitCode
];


(* ============================================================ *)
(* SetupFIRE - Generate .start file for a single topology       *)
(* ============================================================ *)

SetupFIRE[topology_Association, workDir_String:""] :=
Module[{dir, name, result, fireSubst, fireProps, fireRepls, pn, restrictions,
        squares, completedFireProps, originalFirePropCount, numeratorPositions,
        autoDetectRestrictions, startFile, fireVariables,
        runtimeFingerprintRecord, setupFingerprintRecord},
  name = topology["Name"];

  (* Determine working directory *)
  dir = If[workDir === "",
    FileNameJoin[{FeynmanTrick`Private`$FTConfig["WorkDirectory"], name}],
    workDir
  ];
  If[!DirectoryQ[dir], CreateDirectory[dir, CreateIntermediateDirectories -> True]];

  (* Load FIRE6.m if needed *)
  ensureFIRELoaded[];

  (* Assign problem number *)
  pn = $NextProblemNumber++;

  (* Substitute contexted symbols with Global equivalents for FIRE/Fermat *)
  fireSubst = buildFIRESubstitution[topology];
  fireProps = topology["Propagators"] /. fireSubst;
  fireRepls = topology["Replacements"] /. fireSubst;
  originalFirePropCount = Length[fireProps];

  (* CRITICAL: Completely clear all FIRE state before setting up new topology *)
  (* Without this, FIRE reuses stale IBP relations and hangs *)
  Quiet[
    (* Unprotect everything that FIRE protects *)
    Unprotect[FIRE`Internal, FIRE`External, FIRE`Propagators, FIRE`Replacements,
              FIRE`RESTRICTIONS,
              FIRE`PrepareIBPd, FIRE`BackMatrix, FIRE`Squares, FIRE`startinglist];
    (* Clear the "already prepared" flag and cached IBP data - this is critical! *)
    Clear[FIRE`PrepareIBPd, FIRE`BackMatrix, FIRE`Squares, FIRE`startinglist];
    Clear[FIRE`Internal, FIRE`External, FIRE`Propagators, FIRE`Replacements,
          FIRE`RESTRICTIONS];
    (* Clear IBP and sector data *)
    Clear[FIRE`ExampleDimension, FIRE`SBasis0L, FIRE`SBasis0D, FIRE`SBasis0C,
          FIRE`SBasisL, FIRE`SBasisD, FIRE`SBasisA, FIRE`SBasisH, FIRE`SBasisO,
          FIRE`SBasisC, FIRE`SBasisS, FIRE`SBasisR, FIRE`SBasisRL, FIRE`SBasisM,
          FIRE`SBasisN, FIRE`HPI, FIRE`LRules];
    FIRE`Burning = False;
    FIRE`ProblemNumber = 0;
  , {Unprotect::ssym, Clear::ssym, Set::write}];

  (* Set FIRE variables *)
  FIRE`Internal = topology["LoopMomenta"];
  FIRE`External = topology["ExternalMomenta"];
  If[fireRepls =!= {}, FIRE`Replacements = fireRepls];
  FIRE`Propagators = fireProps;

  (* FIRE requires a complete basis of loop-momentum scalar products. Add
     irreducible numerator slots explicitly instead of relying on Prepare[] to
     limp on after PrepareIBP[] reports "Not enough propagators". *)
  squares = FIRE`SquaresEv[];
  completedFireProps = completePropagators[fireProps, squares];
  If[Length[completedFireProps] > Length[fireProps] &&
     FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["  Added ", Length[completedFireProps] - Length[fireProps],
      " irreducible numerator(s) for FIRE."]
  ];
  fireProps = completedFireProps;
  FIRE`Propagators = fireProps;

  restrictions = feynmanTrickRestrictions[topology];
  restrictions = PadRight[#, Length[fireProps], 0] & /@ restrictions;
  numeratorPositions = Range[originalFirePropCount + 1, Length[fireProps]];
  restrictions = Join[
    restrictions,
    {ConstantArray[-1, Length[fireProps]]},
    Table[
      ReplacePart[ConstantArray[0, Length[fireProps]], pos -> 1],
      {pos, numeratorPositions}
    ]
  ];
  FIRE`RESTRICTIONS = restrictions;

  (* Capture the option actually used to build this .start.  Reading the
     mutable global later would make cache identity depend on current state
     rather than the prepared FIRE problem. *)
  autoDetectRestrictions =
    FeynmanTrick`Private`$FTConfig["AutoDetectRestrictions"];

  (* Prepare IBP relations *)
  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["Setting up FIRE for topology: ", name, " (problem ", pn, ")"];
    If[Length[restrictions] > 0,
      Print["  Sector restrictions: ", restrictions];
    ];
  ];
  FIRE`PrepareIBP[];

  (* Prepare sector basis - this should now say "Prepared" not "Already prepared" *)
  FIRE`Prepare[
    AutoDetectRestrictions -> autoDetectRestrictions
  ];

  (* Save start file *)
  Module[{oldDir = Directory[]},
    SetDirectory[dir];
    FIRE`SaveStart[name];
    SetDirectory[oldDir];
  ];

  (* Return updated topology *)
  result = topology;
  result["WorkDirectory"] = dir;
  result["ProblemNumber"] = pn;
  result["StartFileReady"] = True;
  result["OriginalNumPropagators"] = Lookup[
    topology, "OriginalNumPropagators", topology["NumPropagators"]];
  result["OriginalPropagators"] = Lookup[
    topology, "OriginalPropagators", topology["Propagators"]];
  result["Propagators"] = fireProps;
  result["NumPropagators"] = Length[fireProps];
  result["NumeratorPositions"] = numeratorPositions;
  fireVariables = extractVariables[
    fireProps, fireRepls, topology["LoopMomenta"], topology["ExternalMomenta"]];
  result["Variables"] = fireVariables;

  (* Content-address the exact prepared problem.  The full record is retained
     and used in cache keys; the digest is only a compact diagnostic. *)
  startFile = FileNameJoin[{dir, name <> ".start"}];
  runtimeFingerprintRecord = currentFIRERuntimeFingerprintRecord[];
  setupFingerprintRecord = Join[<|
    "Schema" -> $fireSetupFingerprintSchema,
    "VerifiedStartFile" -> FileExistsQ[startFile],
    "StartFileSHA256" -> fileSHA256[startFile],
    "AutoDetectRestrictions" -> autoDetectRestrictions
  |>, runtimeFingerprintRecord, <|
    "LoopMomenta" -> topology["LoopMomenta"],
    "ExternalMomenta" -> topology["ExternalMomenta"],
    "Propagators" -> fireProps,
    "Replacements" -> fireRepls,
    "OriginalNumPropagators" -> result["OriginalNumPropagators"],
    "NumPropagators" -> result["NumPropagators"],
    "EliminatedPositions" -> topologyEliminatedPositions[topology],
    "NumeratorPositions" -> numeratorPositions,
    "Restrictions" -> restrictions,
    "Variables" -> fireVariables
  |>];
  result["SetupFingerprintRecord"] = setupFingerprintRecord;
  result["SetupFingerprint"] = IntegerString[
    Hash[setupFingerprintRecord, "SHA256"], 16, 64];

  (* Store in global registry *)
  $SetupTopologies[name] = result;

  result
];


(* ============================================================ *)
(* SetupFIREBatch - Set up multiple topologies at once          *)
(* ============================================================ *)

SetupFIREBatch[topologies_List, workDir_String:""] :=
Module[{results},
  (* Clear state before batch setup *)
  ClearFIREState[];

  (* Set up each topology *)
  results = SetupFIRE[#, workDir] & /@ topologies;

  results
];


(* ============================================================ *)
(* writeMultiProblemConfig - config file for multiple problems  *)
(* ============================================================ *)

writeMultiProblemConfig[topologies_List, dir_String, configName_String,
                        integralsFile_String, outputFile_String] :=
Module[{vars, content, threads, fthreads, problemLines, allVars},
  threads = FeynmanTrick`Private`$FTConfig["Threads"];
  fthreads = FeynmanTrick`Private`$FTConfig["FThreads"];

  (* Collect all variables from all topologies *)
  allVars = Join[
    {FeynmanTrick`Private`$FTConfig["DimensionVariable"]},
    Flatten[#["Variables"] & /@ topologies]
  ] // DeleteDuplicates;

  (* Build #problem lines *)
  problemLines = StringJoin[
    "#problem           ", ToString[#["ProblemNumber"]], " ",
    #["Name"], ".start\n"
  ] & /@ topologies;

  content = StringJoin[
    "#threads           ", ToString[threads], "\n",
    "#fthreads          ", ToString[fthreads], "\n",
    "#variables         ", StringRiffle[SymbolName /@ allVars, ","], "\n",
    "#start\n",
    Sequence @@ problemLines,
    "#integrals         ", integralsFile, "\n",
    "#output            ", outputFile, "\n"
  ];

  Export[FileNameJoin[{dir, configName}], content, "Text"];
];


(* Write single-problem config (for backwards compatibility) *)
writeSingleProblemConfig[topology_Association, dir_String, configName_String,
                         integralsFile_String:"", outputFile_String:""] :=
Module[{vars, content, threads, fthreads, intFile, outFile, name, pn},
  name = topology["Name"];
  pn = topology["ProblemNumber"];
  threads = FeynmanTrick`Private`$FTConfig["Threads"];
  fthreads = FeynmanTrick`Private`$FTConfig["FThreads"];

  vars = Join[
    {FeynmanTrick`Private`$FTConfig["DimensionVariable"]},
    topology["Variables"]
  ] // DeleteDuplicates;

  intFile = If[integralsFile === "", name <> ".m", integralsFile];
  outFile = If[outputFile === "", name <> ".tables", outputFile];

  content = StringJoin[
    "#threads           ", ToString[threads], "\n",
    "#fthreads          ", ToString[fthreads], "\n",
    "#variables         ", StringRiffle[SymbolName /@ vars, ","], "\n",
    "#start\n",
    "#problem           ", ToString[pn], " ", name, ".start\n",
    "#integrals         ", intFile, "\n",
    "#output            ", outFile, "\n"
  ];

  Export[FileNameJoin[{dir, configName}], content, "Text"];
];


(* ============================================================ *)
(* FindBasis - find master integrals for a single topology      *)
(* ============================================================ *)

FindBasis[topology_Association] :=
Module[{dir, name, fireBin, intFile, result, masters, tablesFile, pn},

  If[!topology["StartFileReady"],
    Print["Error: Must call SetupFIRE before FindBasis."];
    Return[$Failed];
  ];

  dir = topology["WorkDirectory"];
  name = topology["Name"];
  pn = topology["ProblemNumber"];
  fireBin = FileNameJoin[{FeynmanTrick`Private`$FTConfig["FIREPath"], "bin", "FIRE6"}];

  (* Write integrals file with all sector corners *)
  Module[{allSectors, intContent},
    allSectors = basisSectors[topology];
    intContent = StringJoin["{",
      StringRiffle[("{" <> ToString[pn] <> "," <> ToString[#, InputForm] <> "}") & /@ allSectors, ",\n"],
      "}\n"];
    intFile = FileNameJoin[{dir, name <> ".m"}];
    Export[intFile, intContent, "Text"];
  ];

  (* Write config *)
  writeSingleProblemConfig[topology, dir, name <> ".config"];

  (* Run FIRE6 *)
  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["Running FIRE6 to find basis for ", name, "..."];
  ];
  If[runFIRE6[fireBin, dir, name] =!= 0,
    Print["Error: FIRE6 basis discovery failed."];
    Return[$Failed];
  ];

  (* Load tables and extract masters *)
  tablesFile = FileNameJoin[{dir, name <> ".tables"}];
  If[!FileExistsQ[tablesFile],
    Print["Error: FIRE6 did not produce output tables."];
    Return[$Failed];
  ];

  masters = FIRE`Tables2Masters[tablesFile];
  (* Filter for this problem number and extract indices *)
  masters = Cases[masters, {pn, indices_List} :> indices];

  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["Found ", Length[masters], " master integrals."];
  ];

  result = topology;
  result["Masters"] = masters;
  result
];


(* ============================================================ *)
(* FindBasisBatch - find masters for multiple topologies        *)
(* ============================================================ *)

FindBasisBatch[topologies_List] :=
Module[{dir, fireBin, intContent, allSectors, tablesFile, masters, results},

  If[Length[topologies] == 0, Return[{}]];

  (* Use first topology's work directory *)
  dir = topologies[[1]]["WorkDirectory"];
  fireBin = FileNameJoin[{FeynmanTrick`Private`$FTConfig["FIREPath"], "bin", "FIRE6"}];

  (* Build combined integrals file *)
  intContent = StringJoin["{",
    StringRiffle[
      Flatten[
        Table[
          With[{topo = topologies[[i]], pn = topologies[[i]]["ProblemNumber"]},
            allSectors = basisSectors[topo];
            ("{" <> ToString[pn] <> "," <> ToString[#, InputForm] <> "}") & /@ allSectors
          ],
          {i, Length[topologies]}
        ],
        1
      ],
      ",\n"
    ],
    "}\n"
  ];
  Export[FileNameJoin[{dir, "batch_basis.m"}], intContent, "Text"];

  (* Write multi-problem config *)
  writeMultiProblemConfig[topologies, dir, "batch_basis.config",
                          "batch_basis.m", "batch_basis.tables"];

  (* Run FIRE6 once *)
  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["Running FIRE6 to find basis for ", Length[topologies], " topologies..."];
  ];
  If[runFIRE6[fireBin, dir, "batch_basis"] =!= 0,
    Print["Error: FIRE6 batch basis discovery failed."];
    Return[$Failed];
  ];

  (* Load results *)
  tablesFile = FileNameJoin[{dir, "batch_basis.tables"}];
  If[!FileExistsQ[tablesFile],
    Print["Error: FIRE6 did not produce output tables."];
    Return[$Failed];
  ];

  masters = FIRE`Tables2Masters[tablesFile];

  (* Update each topology with its masters *)
  results = Table[
    With[{topo = topologies[[i]], pn = topologies[[i]]["ProblemNumber"]},
      Module[{topoMasters, result},
        topoMasters = Cases[masters, {pn, indices_List} :> indices];
        If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
          Print["  ", topo["Name"], ": ", Length[topoMasters], " masters"];
        ];
        result = topo;
        result["Masters"] = topoMasters;
        result
      ]
    ],
    {i, Length[topologies]}
  ];

  results
];


(* ============================================================ *)
(* ReduceIntegrals - reduce integrals for a single topology     *)
(* ============================================================ *)

ReduceIntegralsDetailed[topology_Association, integrals_List] :=
Module[{dir, name, fireBin, intFile, tablesFile, rules, intContent, result,
        pn, masters = {}, fireIntegrals, cacheEntries, missingPositions,
        missingFireIntegrals, newReductions = {},
        computedReductions, allEntries, allMasters, cacheEnabled},

  If[!topology["StartFileReady"],
    Print["Error: Must call SetupFIRE before ReduceIntegrals."];
    Return[$Failed];
  ];

  If[integrals === {},
    Return[<|"Reductions" -> <||>, "Masters" -> {}, "Rules" -> {},
      "TablesFile" -> Missing["NoIntegrals"]|>]
  ];

  fireIntegrals = normalizeIntegralIndices[topology, integrals];
  If[MemberQ[fireIntegrals, $Failed], Return[$Failed]];

  cacheEnabled = TrueQ[
    Lookup[FeynmanTrick`Private`$FTConfig, "ReductionCache", True]
  ];
  cacheEntries = (cachedReduction[topology, #] &) /@ fireIntegrals;
  missingPositions = If[cacheEnabled,
    Flatten[Position[cacheEntries, _Missing, {1}]],
    Range[Length[fireIntegrals]]
  ];

  If[missingPositions === {},
    (* Cached rational functions still carry the setup-time dimension symbol.
       Reject a semantic mismatch before a caller can expand, for example,
       1/(d-4) as though the active regulator lived in a different symbol. *)
    If[!TrueQ[cachedTopologySemanticCompatibleQ[topology]],
      Return[$Failed, Module]];
    result = Association[
      MapThread[Rule, {integrals, ((#["Reduction"] &) /@ cacheEntries)}]
    ];
    allMasters = DeleteDuplicates[
      Flatten[Cases[cacheEntries, assoc_Association :> assoc["Masters"]], 1]
    ];
    Return[
      <|
        "Reductions" -> result,
        "Masters" -> allMasters,
        "Rules" -> {},
        "TablesFile" -> Missing["Cached"]
      |>,
      Module
    ];
  ];

  (* A cache miss launches the current FIRE binary against the prepared
     .start.  Refuse that launch if the start file, FIRE runtime, or dimension
     variable no longer matches the setup-time fingerprint carried by the
     topology.  Exact cache hits need no FIRE runtime and remain usable. *)
  If[!TrueQ[preparedTopologyCompatibleQ[topology]], Return[$Failed, Module]];

  ensureFIRELoaded[];

  dir = topology["WorkDirectory"];
  name = topology["Name"];
  pn = topology["ProblemNumber"];
  fireBin = FileNameJoin[{FeynmanTrick`Private`$FTConfig["FIREPath"], "bin", "FIRE6"}];

  (* Different public indices can become identical after zero-padding, and a
     caller may also repeat an index verbatim.  FIRE accepts both but repeats
     the same reduction work.  Stable de-duplication preserves the first-seen
     request order and the result is still expanded back to every original
     key below. *)
  missingFireIntegrals = DeleteDuplicates[
    fireIntegrals[[missingPositions]]
  ];

  (* Write integrals file with problem number *)
  intContent = StringJoin[
    "{",
    StringRiffle[
      ("{" <> ToString[pn] <> "," <> ToString[#, InputForm] <> "}") & /@ missingFireIntegrals,
      ",\n"
    ],
    "}\n"
  ];
  intFile = name <> "_reduce.m";
  Export[FileNameJoin[{dir, intFile}], intContent, "Text"];

  (* Write config *)
  writeSingleProblemConfig[topology, dir, name <> "_reduce.config",
                           intFile, name <> "_reduce.tables"];

  (* Run FIRE6 *)
  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 2,
    Print["Reducing ", Length[missingFireIntegrals], " integrals..."];
  ];
  If[runFIRE6[fireBin, dir, name <> "_reduce"] =!= 0,
    Print["Error: FIRE6 integral reduction failed."];
    Return[$Failed];
  ];

  (* Load results *)
  tablesFile = FileNameJoin[{dir, name <> "_reduce.tables"}];
  If[!FileExistsQ[tablesFile],
    Print["Error: FIRE6 reduction did not produce output."];
    Return[$Failed];
  ];

  rules = FIRE`Tables2Rules[tablesFile];
  masters = Cases[FIRE`Tables2Masters[tablesFile],
    {pn, indices_List} :> indices
  ];

  (* Convert to our format *)
  (* FIRE uses G[problemNumber, indices] *)
  newReductions = ((Global`G[pn, #] /. rules) &) /@ missingFireIntegrals;
  (* Replace G[pn, ...] with G[1, ...] for consistency *)
  newReductions = newReductions /. Global`G[pn, idx_] :> Global`G[1, idx];

  (* FIRE has already proved that each returned master reduces to itself.
     Cache those identities immediately: the adaptive derivative-closure loop
     often asks for one of them in its next wave, and launching FIRE again for
     that request is pure overhead.  This uses only the existing in-kernel,
     exact-key cache; no FIRE database is retained or shared across runs. *)
  Scan[
    cacheReduction[topology, #, Global`G[1, #], masters] &,
    DeleteDuplicates[masters]
  ];

  MapThread[
    cacheReduction[topology, #1, #2, masters] &,
    {missingFireIntegrals, newReductions}
  ];

  computedReductions = Association[
    MapThread[Rule, {missingFireIntegrals, newReductions}]
  ];
  allEntries = (cachedReduction[topology, #] &) /@ fireIntegrals;

  result = Association[
    MapThread[
      Function[{orig, fire, entry},
        orig -> If[AssociationQ[entry],
          entry["Reduction"],
          computedReductions[fire]
        ]
      ],
      {integrals, fireIntegrals, allEntries}
    ]
  ];

  allMasters = DeleteDuplicates[
    Join[
      masters,
      Flatten[Cases[allEntries, assoc_Association :> assoc["Masters"]], 1]
    ]
  ];

  <|
    "Reductions" -> result,
    "Masters" -> allMasters,
    "Rules" -> rules,
    "TablesFile" -> tablesFile
  |>
];

ReduceIntegrals[topology_Association, integrals_List] :=
Module[{detailed},
  detailed = ReduceIntegralsDetailed[topology, integrals];
  If[detailed === $Failed, Return[$Failed]];
  detailed["Reductions"]
];

(* ============================================================ *)
(* ReduceIntegralsBatch - reduce integrals for multiple topos   *)
(* ============================================================ *)

ReduceIntegralsBatch[topoIntegralPairs_List] :=
Module[{dir, fireBin, intContent, tablesFile, rules, results, allIntegrals},

  If[Length[topoIntegralPairs] == 0, Return[{}]];

  (* Use first topology's work directory *)
  dir = topoIntegralPairs[[1, 1]]["WorkDirectory"];
  fireBin = FileNameJoin[{FeynmanTrick`Private`$FTConfig["FIREPath"], "bin", "FIRE6"}];

  (* Build combined integrals file *)
  allIntegrals = Flatten[
    Table[
      With[{topo = topoIntegralPairs[[i, 1]],
            ints = topoIntegralPairs[[i, 2]],
            pn = topoIntegralPairs[[i, 1]]["ProblemNumber"]},
        Module[{fireInts = normalizeIntegralIndices[topo, ints]},
          If[MemberQ[fireInts, $Failed], Return[$Failed, Module]];
          MapThread[{pn, #1, #2} &, {ints, fireInts}]
        ]
      ],
      {i, Length[topoIntegralPairs]}
    ],
    1
  ];

  intContent = StringJoin[
    "{",
    StringRiffle[
      ("{" <> ToString[#[[1]]] <> "," <> ToString[#[[3]], InputForm] <> "}") & /@ allIntegrals,
      ",\n"
    ],
    "}\n"
  ];
  Export[FileNameJoin[{dir, "batch_reduce.m"}], intContent, "Text"];

  (* Write multi-problem config *)
  writeMultiProblemConfig[topoIntegralPairs[[All, 1]], dir, "batch_reduce.config",
                          "batch_reduce.m", "batch_reduce.tables"];

  (* Run FIRE6 *)
  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["Reducing integrals for ", Length[topoIntegralPairs], " topologies..."];
  ];
  If[runFIRE6[fireBin, dir, "batch_reduce"] =!= 0,
    Print["Error: FIRE6 batch integral reduction failed."];
    Return[$Failed];
  ];

  (* Load results *)
  tablesFile = FileNameJoin[{dir, "batch_reduce.tables"}];
  If[!FileExistsQ[tablesFile],
    Print["Error: FIRE6 reduction did not produce output."];
    Return[$Failed];
  ];

  rules = FIRE`Tables2Rules[tablesFile];

  (* Parse results for each topology *)
  results = Table[
    With[{topo = topoIntegralPairs[[i, 1]],
          ints = topoIntegralPairs[[i, 2]],
          pn = topoIntegralPairs[[i, 1]]["ProblemNumber"]},
      Module[{fireInts = normalizeIntegralIndices[topo, ints]},
        If[MemberQ[fireInts, $Failed], Return[$Failed, Module]];
        Association[
          (Global`G[pn, #] /. rules) & /@ fireInts //
          (* Normalize G[pn, ...] -> G[1, ...] *)
          (# /. Global`G[pn, idx_] :> Global`G[1, idx]) & //
          MapThread[Rule, {ints, #}] &
        ]
      ]
    ],
    {i, Length[topoIntegralPairs]}
  ];

  results
];


(* ============================================================ *)
(* ComputeDiffMatrix                                             *)
(* ============================================================ *)

ComputeDiffMatrix[topology_Association, variable_, precomputedDecomp_:{Automatic}] :=
Module[{masters, nMasters, coeffMat, constVec, shiftedPerMaster,
        uniqueIntegrals, reductions, diffMatrix, i, j, expr},

  masters = topology["Masters"];
  If[masters === {},
    Print["Error: No masters found. Call FindBasis first."];
    Return[$Failed];
  ];
  nMasters = Length[masters];

  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["Computing differential matrix in ", variable,
          " for ", nMasters, " masters..."];
  ];

  (* Step 1: Decompose propagator derivatives *)
  If[precomputedDecomp =!= {Automatic} && Length[precomputedDecomp] == 2,
    {coeffMat, constVec} = precomputedDecomp;
  ,
    {coeffMat, constVec} = DecomposePropagatorDerivative[
      topology["Propagators"],
      topology["LoopMomenta"],
      variable,
      topology["Replacements"]
    ];
  ];

  (* Step 2: For each master, find shifted integrals in d/dx M_i *)
  shiftedPerMaster = Table[
    DifferentiatedIntegrals[masters[[i]], coeffMat, constVec],
    {i, nMasters}
  ];

  (* Step 3: Collect all unique integrals needing reduction *)
  uniqueIntegrals = DeleteDuplicates[
    Flatten[shiftedPerMaster[[All, All, 1]], 1]
  ];

  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["Need to reduce ", Length[uniqueIntegrals], " integrals."];
  ];

  (* Step 4: Reduce all integrals *)
  reductions = ReduceIntegrals[topology, uniqueIntegrals];

  If[reductions === $Failed,
    Print["Error: FIRE reduction failed."];
    Return[$Failed];
  ];

  (* Map FIRE's Global variables back to the requested variable symbol *)
  If[Context[variable] =!= "Global`",
    Module[{globalVar = Symbol[SymbolName[variable]]},
      reductions = Map[(# /. globalVar -> variable) &, reductions];
    ];
  ];

  (* Step 5: Assemble the matrix *)
  diffMatrix = Table[0, {nMasters}, {nMasters}];

  Do[
    expr = Total[
      (#[[2]] * (reductions[#[[1]]])) & /@ shiftedPerMaster[[i]]
    ];
    Do[
      diffMatrix[[i, j]] = Coefficient[expr, Global`G[1, masters[[j]]]];
    , {j, nMasters}];
  , {i, nMasters}];

  diffMatrix = Map[Together, diffMatrix, {2}];

  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["Differential matrix computed (", nMasters, "x", nMasters, ")."];
  ];

  diffMatrix
];


(* ============================================================ *)
(* GetFIREResult                                                 *)
(* ============================================================ *)

GetFIREResult[topology_Association, integral_List] :=
Module[{result},
  result = ReduceIntegrals[topology, {integral}];
  If[result === $Failed, Return[$Failed]];
  result[integral]
];


End[];
EndPackage[];
