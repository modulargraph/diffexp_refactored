(* ::Package:: *)
(* FIREInterface - Wrapper for FIRE 7 IBP reduction *)
(* Handles topology definition, FIRE setup, basis finding, *)
(* integral reduction, and differential matrix computation *)
(*
   ARCHITECTURE: Multi-topology batch processing

   Instead of clearing and restarting FIRE for each topology, we:
   1. Generate all .start files first (SetupFIRE or SetupFIREBatch)
   2. Load all topologies with unique problem numbers
   3. Run FIRE 7 once with all problems in a single config
   4. Parse combined results

   This avoids reducer subprocess races and is how FIRE 7 is
   designed to work with multiple related topologies.
*)

Get[FileNameJoin[{DirectoryName[$InputFileName], "FIRE7Runner.m"}]];

BeginPackage["FeynmanTrick`FIREInterface`", {"FeynmanTrick`", "FeynmanTrick`PropagatorAlgebra`"}];

DefineTopology::usage =
  "DefineTopology[name, loopMomenta, externalMomenta, propagators, replacements] \
defines a topology for FIRE reduction. Propagator convention: D_j = -q_j^2 + m_j^2. \
Returns a Topology association.";

SetupFIRE::usage =
  "SetupFIRE[topology, workDir] generates a .start file for FIRE 7. \
Does NOT run FIRE - call FindBasis or ReduceIntegrals after setup.";

SetupFIREBatch::usage =
  "SetupFIREBatch[{topo1, topo2, ...}, workDir] sets up multiple topologies \
at once, assigning each a unique problem number. Returns list of updated topologies.";

FindBasis::usage =
  "FindBasis[topology] runs FIRE 7 to determine master integrals. \
Returns the topology with Masters field populated.";

FindBasisBatch::usage =
  "FindBasisBatch[{topo1, topo2, ...}] finds master integrals for all topologies \
in a single FIRE 7 run. Returns list of updated topologies.";

ReduceIntegrals::usage =
  "ReduceIntegrals[topology, integrals] reduces a list of integrals \
(given as index vectors) to the master basis using FIRE 7. \
Returns an association: integral -> linear combination of masters.";

ReduceIntegralsDetailed::usage =
  "ReduceIntegralsDetailed[topology, integrals] reduces integrals with FIRE 7 and \
returns an association containing reductions and the masters reported by FIRE.";

ReduceIntegralsBatch::usage =
  "ReduceIntegralsBatch[{{topo1, ints1}, {topo2, ints2}, ...}] reduces \
integrals for multiple topologies in a single FIRE 7 run.";

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
If[!ValueQ[$LoadedFIRE7Source], $LoadedFIRE7Source = None];
If[!ValueQ[$LoadedFIRE7SourceSHA256], $LoadedFIRE7SourceSHA256 = None];

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
Module[{nP, eliminated, numerators, inactive, active},
  nP = topology["NumPropagators"];
  eliminated = topologyEliminatedPositions[topology];
  numerators = Select[Lookup[topology, "NumeratorPositions", {}],
    IntegerQ[#] && 1 <= # <= nP &];
  inactive = DeleteDuplicates[Join[eliminated, numerators]];
  active = Complement[Range[nP], inactive];
  If[active === {}, Return[{}, Module]];
  (* Added irreducible numerators are never positive denominator corners.
     Enumerating them used to turn a five-denominator four-loop family into
     as many as 2^14-1 basis requests. *)
  ReplacePart[ConstantArray[0, nP], Thread[active -> #]] & /@
    Rest[Tuples[{0, 1}, Length[active]]]
];

(* FIRE's reported masters are request-dependent: reducing only a sector
   corner need not expose additional masters represented by a dot or an
   irreducible numerator in that same sector.  Probe the complete first
   complexity shell while keeping numerator growth linear (rather than
   enumerating 2^n spurious positive-numerator sectors). *)
basisSeedIntegrals[topology_Association] := Module[
  {corners, nP, eliminated, numerators, denominatorPositions},
  corners = basisSectors[topology];
  If[corners === {}, Return[{}, Module]];
  nP = topology["NumPropagators"];
  eliminated = topologyEliminatedPositions[topology];
  numerators = Complement[
    Select[Lookup[topology, "NumeratorPositions", {}],
      IntegerQ[#] && 1 <= # <= nP &], eliminated];
  denominatorPositions = Complement[Range[nP], Join[eliminated, numerators]];
  DeleteDuplicates[Flatten[Map[Function[corner,
    Join[
      {corner},
      (ReplacePart[corner, # -> 2] & /@
        Select[denominatorPositions, corner[[#]] === 1 &]),
      (ReplacePart[corner, # -> -1] & /@ numerators)
    ]], corners], 1]]
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

$reductionCacheSchema = "FeynmanTrick.ReductionCache/v4";
$fireSetupFingerprintSchema = "FeynmanTrick.FIRESetup/v2";

fileSHA256[path_String] := If[FileExistsQ[path],
  IntegerString[FileHash[path, "SHA256"], 16, 64], Missing["NotFound"]];

fireConfigValue[key_String, default_] := Lookup[
  FeynmanTrick`Private`$FTConfig, key, default];

resolvedFIREMPIExecutable[] := Module[{configured},
  configured = fireConfigValue["FIREMPIExecutable", Automatic];
  If[configured === Automatic, "mpirun", configured]
];

resolvedExecutablePath[command_String] := Module[
  {explicit, candidates, found},
  explicit = !MemberQ[{"", "."}, DirectoryName[command]];
  If[explicit,
    Return[If[FileExistsQ[ExpandFileName[command]],
      ExpandFileName[command], Missing["NotFound"]], Module]];
  candidates = FileNameJoin[{#, command}] & /@
    StringSplit[Environment["PATH"],
      If[$OperatingSystem === "Windows", ";", ":"]];
  found = SelectFirst[candidates, FileExistsQ, Missing["NotFound"]];
  If[StringQ[found], ExpandFileName[found], found]
];

resolvedFIREMPIExecutablePath[] :=
  resolvedExecutablePath[resolvedFIREMPIExecutable[]];

currentFIRESourceSHA256[firePath_String] := Module[
  {source, diskHash},
  source = FileNameJoin[{firePath, "FIRE7.m"}];
  diskHash = fileSHA256[source];
  If[$LoadedFIRE7Source === source &&
      StringQ[$LoadedFIRE7SourceSHA256],
    If[diskHash === $LoadedFIRE7SourceSHA256,
      $LoadedFIRE7SourceSHA256, Missing["LoadedSourceChanged"]],
    diskHash]
];

commandVersion[command_String] := Module[{result, text},
  result = Quiet[Check[RunProcess[{command, "--version"}], $Failed],
    {RunProcess::pnfd, RunProcess::posix}];
  If[result === $Failed || !AssociationQ[result] ||
      Lookup[result, "ExitCode", 1] =!= 0,
    Return[Missing["Unavailable"], Module]];
  text = StringTrim[Lookup[result, "StandardOutput", ""] <> "\n" <>
    Lookup[result, "StandardError", ""]];
  If[text === "", Missing["Unavailable"], First[StringSplit[text, "\n"]]]
];

currentFIRERuntimeFingerprintRecord[] := Module[
  {firePath, backend, common, mpiPath},
  firePath = ExpandFileName[FeynmanTrick`Private`$FTConfig["FIREPath"]];
  backend = fireConfigValue["FIREBackend", "Modular"];
  mpiPath = If[backend === "Modular",
    resolvedFIREMPIExecutablePath[], Missing["Inactive"]];
  common = <|
    "FIREPath" -> firePath,
    "FIRESourceSHA256" -> currentFIRESourceSHA256[firePath],
    "FIREBackend" -> backend,
    "FIREVariableOrderPolicy" ->
      "kinematics-first-dimension-last/v1",
    "SystemID" -> $SystemID,
    "DimensionVariable" ->
      FeynmanTrick`Private`$FTConfig["DimensionVariable"]
  |>;
  Join[common, Switch[backend,
    "Classical", <|
      "FIREBinarySHA256" -> fileSHA256[
        FileNameJoin[{firePath, "bin", "FIRE7"}]]|>,
    "Modular", <|
      "FIREPrimeBinarySHA256" -> fileSHA256[
        FileNameJoin[{firePath, "bin", "FIRE7p"}]],
      "FIREMultiprimeBinarySHA256" -> fileSHA256[
        FileNameJoin[{firePath, "bin", "FIRE7mp"}]],
      "FIREMPIBinarySHA256" -> fileSHA256[
        FileNameJoin[{firePath, "bin", "FIRE7_MPI"}]],
      "FIREReconstructBinarySHA256" -> fileSHA256[
        FileNameJoin[{firePath, "bin", "reconstruct"}]],
      "FIREMPILauncherPath" -> mpiPath,
      "FIREMPILauncherSHA256" -> If[StringQ[mpiPath],
        fileSHA256[mpiPath], Missing["NotFound"]],
      "FIREPathsIncSHA256" -> fileSHA256[
        FileNameJoin[{firePath, "paths.inc"}]],
      "FIRECalc" -> fireConfigValue["FIRECalc", "flint"],
      "FIRESamplingPolicy" -> "geometric-zippel/v1",
      "FIREVariableBasePolicy" -> "big-primes/v1",
      "FIREUseMultiprime" -> fireConfigValue["FIREUseMultiprime", True],
      "FIREDimensionSeparated" ->
        fireConfigValue["FIREDimensionSeparated", False],
      "FIREMultiprimeWidth" ->
        fireConfigValue["FIREMultiprimeWidth", 16]|>,
    _, <|"InvalidFIREBackend" -> backend|>
  ]]
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
  keys = Keys[current];
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

preparedStartHashes[topologies_List] := Module[{pairs, grouped},
  pairs = Map[Function[topology,
    With[{record = Lookup[topology, "SetupFingerprintRecord", <||>]},
      topology["Name"] <> ".start" ->
        Lookup[record, "StartFileSHA256", Missing["Absent"]]]], topologies];
  If[!AllTrue[Last /@ pairs, StringQ[#] && StringLength[#] === 64 &],
    Return[$Failed, Module]];
  grouped = GatherBy[pairs, First];
  If[AnyTrue[grouped, Length[DeleteDuplicates[Last /@ #]] =!= 1 &],
    Return[$Failed, Module]];
  Association[pairs]
];

preparedRunnerRuntimeHashes[topologies_List, operation_String] := Module[
  {records, first, backend, useMultiprime, keyMap, expected},
  records = Lookup[topologies, "SetupFingerprintRecord", Missing["Absent"]];
  If[records === {} || !AllTrue[records, AssociationQ],
    Return[$Failed, Module]];
  first = First[records];
  backend = Lookup[first, "FIREBackend", Missing["Absent"]];
  useMultiprime = TrueQ[Lookup[first, "FIREUseMultiprime", False]];
  keyMap = Switch[{backend, operation},
    {"Classical", _}, <|
      "FIREBinarySHA256" -> "FIREBinarySHA256"|>,
    {"Modular", "Basis"}, <|
      "FIRE7SourceSHA256" -> "FIRESourceSHA256",
      "FIRE7PrimeSHA256" -> "FIREPrimeBinarySHA256"|>,
    {"Modular", "Reduction"}, Join[<|
      "FIRE7SourceSHA256" -> "FIRESourceSHA256",
      "FIRE7PrimeSHA256" -> "FIREPrimeBinarySHA256",
      "FIRE7MPISHA256" -> "FIREMPIBinarySHA256",
      "ReconstructSHA256" -> "FIREReconstructBinarySHA256"|>,
      <|"MPILauncherSHA256" -> "FIREMPILauncherSHA256"|>,
      If[useMultiprime, <|
        "FIRE7MultiprimeSHA256" -> "FIREMultiprimeBinarySHA256",
        "PathsIncSHA256" -> "FIREPathsIncSHA256"|>, <||>]],
    _, Return[$Failed, Module]];
  expected = Association@KeyValueMap[
    #1 -> Lookup[first, #2, Missing["Absent"]] &, keyMap];
  If[!AllTrue[Values[expected], StringQ[#] && StringLength[#] === 64 &] ||
      !AllTrue[Rest[records], Function[record,
        And @@ KeyValueMap[
          Lookup[record, keyMap[#1], Missing["Absent"]] === #2 &,
          expected]]],
    Return[$Failed, Module]];
  expected
];

reductionCacheKey[topology_Association, fireIntegral_List] := {
  $reductionCacheSchema,
  reductionSetupFingerprintRecord[topology],
  Lookup[topology, "Masters", {}],
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


globalFIRESymbol[s_Symbol] := Symbol["Global`" <> SymbolName[s]];

symbolNameCollisions[symbols_List] := Select[
  GatherBy[DeleteDuplicates[symbols], SymbolName], Length[#] > 1 &];

(* Extract the coefficients that remain after the replacement closure.
   Preserve first occurrence order: it becomes the modular sampling order. *)
extractVariables[props_List, repls_List, loops_List, exts_List] :=
Module[{allInputSymbols, collisions, subst, normalizedProps,
        normalizedRepls, rhs, reduced, allSyms, momenta},
  allInputSymbols = DeleteDuplicates[Cases[
    {props, repls, loops, exts}, _Symbol, Infinity]];
  collisions = symbolNameCollisions[allInputSymbols];
  If[collisions =!= {},
    Print["Error: FIRE symbol basenames collide across contexts: ",
      SymbolName /@ # & /@ collisions];
    Return[$Failed, Module]];
  subst = Rule[#, globalFIRESymbol[#]] & /@ Select[allInputSymbols,
    Context[#] =!= "Global`" &];
  normalizedProps = props /. subst;
  normalizedRepls = repls /. subst;
  rhs = Cases[normalizedRepls, (Rule | RuleDelayed)[_, value_] :> value];
  reduced = Quiet[Check[
    Join[normalizedProps, rhs] //. normalizedRepls,
    Join[normalizedProps, rhs]]];
  allSyms = DeleteDuplicates[Cases[reduced, _Symbol, Infinity]];
  momenta = globalFIRESymbol /@ Join[loops, exts];
  Select[allSyms, !MemberQ[momenta, #] &]
];

(* Map contexted coefficients from both propagators and replacement rules to
   Global`, which is the namespace serialized by FIRE. *)
buildFIRESubstitution[topology_Association] :=
Module[{allSyms, collisions, dimension},
  allSyms = DeleteDuplicates[Cases[{
    topology["Propagators"], topology["Replacements"],
    topology["LoopMomenta"], topology["ExternalMomenta"]},
    _Symbol, Infinity]];
  dimension = FeynmanTrick`Private`$FTConfig["DimensionVariable"];
  collisions = symbolNameCollisions[Append[allSyms, dimension]];
  If[MemberQ[allSyms, dimension] || collisions =!= {},
    Print["Error: FIRE symbol basenames collide across contexts: ",
      InputForm[If[collisions === {}, {{dimension}}, collisions]]];
    Return[$Failed, Module]];
  Rule[#, globalFIRESymbol[#]] & /@
    Select[allSyms, Context[#] =!= "Global`" &]
];

exactFIREInputQ[expression_] := FreeQ[expression,
  _Real | Indeterminate | ComplexInfinity | _DirectedInfinity];

validFIRENumericalPointQ[rules_] :=
  ListQ[rules] && AllTrue[rules, MatchQ[#, _Rule] &] &&
  DuplicateFreeQ[First /@ rules] &&
  AllTrue[First /@ rules, Head[#] === Symbol &] &&
  exactFIREInputQ[rules];

(* Apply a numerical point to a fixed point so chained exact assignments such
   as {s -> q, q -> -1} are resolved deterministically.  Refuse cyclic or
   self-expanding assignments instead of relying on the kernel's iteration
   limit. *)
applyFIRENumericalPoint[expression_, rules_List] := Module[
  {current = expression, next, seen = {}, steps = 0,
   maxSteps = Length[rules] + 1},
  While[steps <= maxSteps,
    next = current /. rules;
    If[SameQ[next, current], Return[current, Module]];
    If[AnyTrue[seen, SameQ[#, next] &], Return[$Failed, Module]];
    AppendTo[seen, current];
    current = next;
    steps++;
  ];
  $Failed
];

(* NumericalPoint fixes coefficients, not the scalar-product patterns on the
   left of kinematic replacement rules. *)
applyFIRENumericalPointToReplacements[replacements_List, rules_List] :=
Module[{rhs},
  rhs = applyFIRENumericalPoint[Last /@ replacements, rules];
  If[rhs === $Failed, $Failed,
    MapThread[Rule, {First /@ replacements, rhs}]]
];


(* ============================================================ *)
(* ClearFIREState - reset all FIRE state                        *)
(* ============================================================ *)

ClearFIREState[] := Module[{},
  $SetupTopologies = <||>;
  $NextProblemNumber = 1;
  $ReductionCache = <||>;

  If[ensureFIRELoaded[] === $Failed, Return[$Failed, Module]];

  If[Length[DownValues[FIRE`ClearIBP]] > 0,
    Quiet[Check[FIRE`ClearIBP[], Null]]];

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
(* ensureFIRELoaded - load the configured FIRE7.m exactly once *)
(* ============================================================ *)

fireSaveStartDefinedQ[] := Total[Length /@ {
  OwnValues[FIRE`SaveStart], DownValues[FIRE`SaveStart],
  SubValues[FIRE`SaveStart], UpValues[FIRE`SaveStart]}] > 0;

ensureFIRELoaded[] := Module[
  {firePath, source, sourceHash, loadedQ, result, afterHash},
  firePath = ExpandFileName[FeynmanTrick`Private`$FTConfig["FIREPath"]];
  source = FileNameJoin[{firePath, "FIRE7.m"}];
  If[!FileExistsQ[source],
    Print["Error: FIRE7.m does not exist at ", source];
    Return[$Failed, Module]];
  sourceHash = fileSHA256[source];
  If[!StringQ[sourceHash] || StringLength[sourceHash] =!= 64,
    Print["Error: could not fingerprint FIRE7.m at ", source];
    Return[$Failed, Module]];
  loadedQ = fireSaveStartDefinedQ[];
  If[loadedQ,
    If[$LoadedFIRE7Source === source &&
        $LoadedFIRE7SourceSHA256 === sourceHash,
      Return[True, Module]];
    Print["Error: a FIRE package is already loaded from an unverified source. Start a fresh Wolfram kernel before changing FIREPath."];
    Return[$Failed, Module]];
  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 2,
    Print["Loading FIRE7.m from: ", firePath]];
  result = Quiet[Check[Get[source], $Failed]];
  If[result === $Failed || !fireSaveStartDefinedQ[],
    Print["Error: failed to load FIRE7.m from ", source];
    Return[$Failed, Module]];
  afterHash = fileSHA256[source];
  If[afterHash =!= sourceHash,
    (* FIRE definitions may now reflect mixed file contents.  Poison this
       kernel's provenance state permanently; restoring the old file must not
       make those definitions trusted. *)
    $LoadedFIRE7Source = source;
    $LoadedFIRE7SourceSHA256 = Missing["ChangedDuringLoad"];
    Print["Error: FIRE7.m changed while it was being loaded. Start a fresh Wolfram kernel."];
    Return[$Failed, Module]];
  $LoadedFIRE7Source = source;
  $LoadedFIRE7SourceSHA256 = sourceHash;
  True
];


(* ============================================================ *)
(* FIRE runner helpers                                           *)
(* ============================================================ *)

safeReadString[path_String] := Module[{txt},
  If[!FileExistsQ[path], Return["", Module]];
  txt = Quiet[ReadString[path], {ReadString::read, General::stop}];
  If[StringQ[txt], txt, ""]
];

fire7RunnerSettings[] := Module[{workers, mpi},
  workers = fireConfigValue["FIREModularWorkers", Automatic];
  workers = Replace[workers, Automatic :> Max[1, Min[10, $ProcessorCount]]];
  mpi = resolvedFIREMPIExecutable[];
  <|
    "Calc" -> fireConfigValue["FIRECalc", "flint"],
    "Workers" -> workers,
    "UseMultiprime" -> fireConfigValue["FIREUseMultiprime", True],
    "PrimeLimit" -> fireConfigValue["FIREPrimeLimit", 127],
    "DimensionSeparated" ->
      fireConfigValue["FIREDimensionSeparated", False],
    "MultiprimeWidth" -> fireConfigValue["FIREMultiprimeWidth", 16],
    "MPIExecutable" -> mpi,
    "BasisProbeCount" -> fireConfigValue["FIREBasisProbeCount", 2],
    "KeepTables" -> fireConfigValue["FIREKeepModularTables", True],
    "CacheDirectory" -> fireConfigValue["FIREModularCacheDirectory",
      FileNameJoin[{$TemporaryDirectory, "FeynmanTrick_FIRE7_Modular"}]],
    "TimeoutSeconds" -> fireConfigValue["FIRETimeoutSeconds", 1800]
  |>
];

publishFIRETable[source_String, destination_String] := Module[{temp, result},
  temp = destination <> ".tmp-" <> StringReplace[CreateUUID[], "-" -> ""];
  result = Quiet[Check[
    CopyFile[source, temp, OverwriteTarget -> True];
    RenameFile[temp, destination, OverwriteTarget -> True];
    True, False]];
  If[FileExistsQ[temp], Quiet[DeleteFile[temp]]];
  result
];

canonicalMasterList[table_String] := Module[{masters},
  masters = Quiet[Check[FIRE`Tables2Masters[table], $Failed]];
  If[!ListQ[masters] || !AllTrue[masters,
      MatchQ[#, {_Integer, {_Integer ...}}] &], Return[$Failed, Module]];
  SortBy[DeleteDuplicates[masters], ToString[#, InputForm] &]
];

validFIREPairsForRequestsQ[pairs_List, requested_List] := Module[
  {arities},
  If[!AllTrue[Join[pairs, requested],
      MatchQ[#, {_Integer, {_Integer ...}}] &], Return[False, Module]];
  arities = GroupBy[DeleteDuplicates[requested],
    First -> (Length[Last[#]] &)];
  If[AnyTrue[Values[arities], Length[DeleteDuplicates[#]] =!= 1 &],
    Return[False, Module]];
  AllTrue[pairs, Function[pair,
    KeyExistsQ[arities, First[pair]] &&
      MemberQ[arities[First[pair]], Length[Last[pair]]]]]
];

(* Inspect only FIRE's relation-ID map, never modular coefficients.  This
   proves that every requested basis seed is structurally present in a prime
   table without feeding finite-field residues to Tables2Rules. *)
rawFIRETableCoverage[table_String] := Module[
  {raw, relations, idMap, relationIDs, mapping},
  raw = Quiet[Check[Get[table], $Failed]];
  If[!MatchQ[raw, {_List, _List}], Return[$Failed, Module]];
  {relations, idMap} = raw;
  If[!AllTrue[relations, MatchQ[#, {_Integer, _List}] &] ||
      !AllTrue[idMap,
        MatchQ[#, {_Integer, {_Integer, {_Integer ...}}}] &],
    Return[$Failed, Module]];
  relationIDs = First /@ relations;
  If[!DuplicateFreeQ[First /@ idMap], Return[$Failed, Module]];
  mapping = Association[Rule @@@ idMap];
  If[!AllTrue[relationIDs, KeyExistsQ[mapping, #] &],
    Return[$Failed, Module]];
  DeleteDuplicates[Lookup[mapping, relationIDs]]
];

parseAndValidateFIRETable[table_String, requested_List,
    allowedProblems_List] := Module[
  {rules, masters, rulePairs, rhsIntegrals, rhsPairs, covered, malformedRules,
   requestedArities, validPairQ},
  If[!FileExistsQ[table] || Quiet[Check[FileByteCount[table], 0]] <= 0,
    Return[$Failed, Module]];
  rules = CheckAbort[Quiet[Check[
    FIRE`Tables2Rules[table, Identity, True,
      "ResultVar" -> Global`G], $Failed]], $Failed];
  masters = CheckAbort[Quiet[Check[FIRE`Tables2Masters[table], $Failed]],
    $Failed];
  If[!ListQ[rules] || !ListQ[masters] ||
      !AllTrue[masters, MatchQ[#, {_Integer, {_Integer ...}}] &],
    Return[$Failed, Module]];
  If[!FreeQ[rules, _Real | Indeterminate | ComplexInfinity |
      DirectedInfinity[_]], Return[$Failed, Module]];
  malformedRules = Select[rules,
    !MatchQ[#, HoldPattern[Rule[Global`G[_Integer, {_Integer ...}], _]]] &];
  If[malformedRules =!= {}, Return[$Failed, Module]];
  (* FIRE7 reconstruction can retain distinct internal relation IDs for the
     same requested master.  Tables2Rules then emits the same exact identity
     rule more than once.  Canonicalize those harmless duplicates before the
     functional-LHS check below; distinct right-hand sides for one integral
     remain a malformed table and are still rejected. *)
  rules = DeleteDuplicates[rules];
  rulePairs = Cases[rules,
    HoldPattern[Rule[Global`G[p_Integer, idx_List], _]] :> {p, idx}];
  If[!DuplicateFreeQ[rulePairs], Return[$Failed, Module]];
  rhsIntegrals = Cases[Last /@ rules,
    HoldPattern[Global`G[___]], Infinity];
  If[!AllTrue[rhsIntegrals,
      MatchQ[#, HoldPattern[Global`G[_Integer, {_Integer ...}]]] &],
    Return[$Failed, Module]];
  rhsPairs = DeleteDuplicates[Cases[Last /@ rules,
    Global`G[p_Integer, idx_List] :> {p, idx}, Infinity]];
  requestedArities = GroupBy[DeleteDuplicates[requested], First ->
    (Length[Last[#]] &)];
  If[AnyTrue[Values[requestedArities],
      Length[DeleteDuplicates[#]] =!= 1 &],
    Return[$Failed, Module]];
  validPairQ[pair_] := MemberQ[allowedProblems, First[pair]] &&
    KeyExistsQ[requestedArities, First[pair]] &&
    MemberQ[requestedArities[First[pair]], Length[Last[pair]]];
  If[!AllTrue[Join[rulePairs, masters, rhsPairs], validPairQ],
    Return[$Failed, Module]];
  If[!SubsetQ[masters, rhsPairs], Return[$Failed, Module]];
  covered = DeleteDuplicates[Join[rulePairs, masters]];
  If[!SubsetQ[covered, DeleteDuplicates[requested]],
    Return[$Failed, Module]];
  <|"Rules" -> rules, "Masters" -> DeleteDuplicates[masters],
    "RuleIntegrals" -> rulePairs|>
];

runFIRE6BasisProbes[fireRoot_String, dir_String, configName_String,
    expectedInputHashes_:<||>, expectedRuntimeHashes_:<||>] :=
Module[{probeResult, tables, masterSets, manifest, canonical, allowedProblems,
        settings, requested},
  manifest = FeynmanTrick`FIRE7Runner`ParseConfig[dir, configName];
  If[FailureQ[manifest],
    Print["FIRE7 basis config failed validation: ", manifest];
    Return[1, Module]];
  allowedProblems = Lookup[manifest["Problems"], "ProblemNumber"];
  requested = Quiet[Check[
    Get[FileNameJoin[{dir, manifest["IntegralsFile"]}]], $Failed]];
  If[!ListQ[requested] ||
      !AllTrue[requested, MatchQ[#, {_Integer, {_Integer ...}}] &],
    Print["FIRE7 basis integral list is malformed."];
    Return[1, Module]];
  settings = Join[fire7RunnerSettings[], <|
    "ExpectedInputHashes" -> expectedInputHashes,
    "ExpectedRuntimeHashes" -> expectedRuntimeHashes,
    "BasisSetValidator" -> With[
      {allowed = allowedProblems, wanted = DeleteDuplicates[requested]},
      Function[probeTables, Module[{sets, coverage},
        sets = canonicalMasterList /@ probeTables;
        coverage = rawFIRETableCoverage /@ probeTables;
        !MemberQ[sets, $Failed] && SameQ @@ sets &&
          AllTrue[sets, validFIREPairsForRequestsQ[#, wanted] &] &&
          !MemberQ[coverage, $Failed] &&
          AllTrue[coverage, SubsetQ[#, wanted] &] &&
          AllTrue[Flatten[coverage, 1], Function[pair,
            MemberQ[allowed, First[pair]] &&
              AnyTrue[wanted, First[#] === First[pair] &&
                Length[Last[#]] === Length[Last[pair]] &]]] &&
          AllTrue[Flatten[sets, 1],
            MatchQ[#, {p_Integer, {_Integer ...}} /; MemberQ[allowed, p]] &]
      ]]]|>];
  probeResult = FeynmanTrick`FIRE7Runner`RunBasisProbes[
    fireRoot, dir, configName, settings];
  If[FailureQ[probeResult],
    Print["FIRE7 basis probes failed: ", probeResult];
    Return[1, Module]];
  tables = probeResult["Tables"];
  masterSets = canonicalMasterList /@ tables;
  If[MemberQ[masterSets, $Failed] || !SameQ @@ masterSets,
    Print["Error: independent FIRE7 finite-field basis probes disagree or are malformed."];
    Return[2, Module]];
  manifest = probeResult["Manifest"];
  canonical = FileNameJoin[{dir, manifest["OutputFile"]}];
  If[!publishFIRETable[First[tables], canonical],
    Print["Error: could not publish FIRE7 basis table."];
    Return[3, Module]];
  0
];

(* Keep the historical private name as a test/extension seam.  It now
   dispatches explicitly to FIRE 7 classical or modular execution and never
   falls back from one backend to the other. *)
runFIRE6[fireBin_String, dir_String, configName_String,
    expectedInputHashes_:<||>, operation_:Automatic,
    expectedRuntimeHashes_:<||>] :=
Module[{backend, fireRoot, outcome, manifest, requested, allowed, settings,
        resolvedOperation},
  backend = fireConfigValue["FIREBackend", "Modular"];
  fireRoot = ExpandFileName[FeynmanTrick`Private`$FTConfig["FIREPath"]];
  resolvedOperation = Replace[operation, Automatic :>
    If[StringEndsQ[configName, "_reduce"], "Reduction", "Basis"]];
  If[!MemberQ[{"Basis", "Reduction"}, resolvedOperation], Return[2, Module]];
  Switch[backend,
    "Classical", runFIRE6Once[fireBin, dir, configName,
      expectedInputHashes, expectedRuntimeHashes],
    "Modular",
      If[resolvedOperation === "Basis",
        runFIRE6BasisProbes[fireRoot, dir, configName,
          expectedInputHashes, expectedRuntimeHashes],
        manifest = FeynmanTrick`FIRE7Runner`ParseConfig[dir, configName];
        If[FailureQ[manifest],
          Print["FIRE7 reduction config failed validation: ", manifest];
          Return[1, Module]];
        requested = Quiet[Check[Get[FileNameJoin[
          {dir, manifest["IntegralsFile"]}]], $Failed]];
        allowed = Lookup[manifest["Problems"], "ProblemNumber"];
        If[!ListQ[requested] ||
            !AllTrue[requested, MatchQ[#, {_Integer, {_Integer ...}}] &],
          Print["FIRE7 reduction integral list is malformed."];
          Return[1, Module]];
        settings = Join[fire7RunnerSettings[], <|
          "ExpectedInputHashes" -> expectedInputHashes,
          "ExpectedRuntimeHashes" -> expectedRuntimeHashes,
          "FinalValidator" -> With[
            {wanted = requested, problems = allowed},
            Function[table, AssociationQ[
              parseAndValidateFIRETable[table, wanted, problems]]]]|>];
        outcome = FeynmanTrick`FIRE7Runner`RunModular[
          fireRoot, dir, configName, settings];
        If[FailureQ[outcome],
          Print["FIRE7 modular reduction failed: ", outcome]; 1,
          If[TrueQ[Lookup[outcome, "Success", False]], 0, 1]]],
    _, Print["Error: FIREBackend must be \"Modular\" or \"Classical\"."]; 2
  ]
];

runFIRE6Once[fireBin_String, dir_String, configName_String,
    expectedInputHashes_:<||>, expectedRuntimeHashes_:<||>] :=
Module[{exitCode, result, logFile, timeoutSeconds, runDir, inputFiles,
        stdout, stderr, command, manifest, generatedTable, canonicalTable,
        configTarget, configText, configLines, artifactsDir, sourceHashes,
        copiedHashes},
  logFile = FileNameJoin[{dir, "fire_stdout.log"}];

  manifest = FeynmanTrick`FIRE7Runner`ParseConfig[dir, configName];
  If[FailureQ[manifest],
    Print["FIRE7 Classical config failed validation: ", manifest];
    Return[126, Module]];
  If[!AssociationQ[expectedInputHashes], Return[126, Module]];
  If[!AssociationQ[expectedRuntimeHashes] ||
      Lookup[expectedRuntimeHashes, "FIREBinarySHA256", fileSHA256[fireBin]] =!=
        fileSHA256[fireBin],
    Return[126, Module]];

  timeoutSeconds = Lookup[FeynmanTrick`Private`$FTConfig, "FIRETimeoutSeconds", 600];
  If[!IntegerQ[timeoutSeconds] || timeoutSeconds <= 0,
    timeoutSeconds = 600;
  ];

  (* Isolate each run in its own working directory to avoid stale temp/db locks
     from interrupted FIRE workers. *)
  runDir = FileNameJoin[{dir, "run_" <> StringReplace[CreateUUID[], "-" -> ""]}];
  CreateDirectory[runDir, CreateIntermediateDirectories -> True];
  inputFiles = DeleteDuplicates[Join[
    {manifest["ConfigPath"],
      FileNameJoin[{dir, manifest["IntegralsFile"]}]},
    FileNameJoin[{dir, #}] & /@ Lookup[manifest, "MasterFiles", {}],
    FileNameJoin[{dir, #}] & /@ Lookup[manifest["Problems"], "StartFile"]
  ]];
  sourceHashes = AssociationThread[FileNameTake /@ inputFiles,
    fileSHA256 /@ inputFiles];
  If[!And @@ KeyValueMap[
      Lookup[sourceHashes, #1, Missing["Absent"]] === #2 &,
      expectedInputHashes],
    Quiet[DeleteDirectory[runDir, DeleteContents -> True]];
    Return[126, Module]];
  Scan[
    Quiet[CopyFile[#, FileNameJoin[{runDir, FileNameTake[#]}], OverwriteTarget -> True]] &,
    inputFiles
  ];
  copiedHashes = AssociationThread[FileNameTake /@ inputFiles,
    fileSHA256[FileNameJoin[{runDir, FileNameTake[#]}]] & /@ inputFiles];
  If[sourceHashes =!= AssociationThread[FileNameTake /@ inputFiles,
        fileSHA256 /@ inputFiles] || copiedHashes =!= sourceHashes,
    Quiet[DeleteDirectory[runDir, DeleteContents -> True]];
    Return[126, Module]];

  (* FIRE 7.1 calls create_directories on the output parent even when the
     configured path is a bare filename.  Give the isolated Classical run a
     real parent exactly as the modular runner does, then publish only that
     rewritten config's declared output. *)
  artifactsDir = FileNameJoin[{runDir, "artifacts"}];
  CreateDirectory[artifactsDir, CreateIntermediateDirectories -> True];
  configTarget = FileNameJoin[{runDir, manifest["ConfigFile"]}];
  configText = safeReadString[configTarget];
  configLines = StringSplit[configText, {"\r\n", "\n", "\r"}];
  configLines = Replace[configLines, line_ /;
      StringStartsQ[StringTrim[line], "#output"] :>
        "#output            artifacts/" <> manifest["OutputFile"], {1}];
  If[Export[configTarget, StringRiffle[configLines, "\n"] <> "\n",
      "Text"] === $Failed,
    Quiet[DeleteDirectory[runDir, DeleteContents -> True]];
    Return[126, Module]];

  command = {ExpandFileName[fireBin], "-c", manifest["ConfigStem"]};
  result = FeynmanTrick`FIRE7Runner`RunSupervised[command, runDir,
    timeoutSeconds, FileNameJoin[{runDir, "fire7_classical"}]];
  exitCode = If[FailureQ[result], 125, Lookup[result, "ExitCode", 125]];

  generatedTable = FileNameJoin[{artifactsDir, manifest["OutputFile"]}];
  canonicalTable = FileNameJoin[{dir, manifest["OutputFile"]}];
  If[KeyExistsQ[expectedRuntimeHashes, "FIREBinarySHA256"] &&
      expectedRuntimeHashes["FIREBinarySHA256"] =!= fileSHA256[fireBin],
    exitCode = 126];
  (* FIRE7 process exit alone is not success.  Only this isolated attempt's
     exact config-declared, nonempty output may be published; unrelated tables
     and stale caller-side canonicals are never evidence of completion. *)
  If[exitCode === 0 &&
      (!FileExistsQ[generatedTable] ||
       Quiet[Check[FileByteCount[generatedTable], 0]] <= 0),
    exitCode = 126];
  If[exitCode === 0 &&
      !publishFIRETable[generatedTable, canonicalTable],
    exitCode = 127];

  stdout = If[AssociationQ[result],
    safeReadString[result["StandardOutputFile"]], ""];
  stderr = If[AssociationQ[result],
    safeReadString[result["StandardErrorFile"]], ToString[result, InputForm]];

  (* Save output to log file *)
  Export[logFile, stdout <> stderr, "Text"];

  (* Best-effort cleanup of this attempt's isolated run directory. *)
  Quiet[DeleteDirectory[runDir, DeleteContents -> True]];

  (* Print output on verbose or failure *)
  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 2 || exitCode =!= 0,
    If[StringLength[stdout] > 0,
      Print["FIRE7 output (tail): "];
      Print[StringTake[stdout, -Min[500, StringLength[stdout]]]];
    ];
    If[exitCode =!= 0,
      Print["FIRE7 exit code: ", exitCode];
      If[StringLength[stderr] > 0,
        Print["FIRE7 stderr: ", StringTake[stderr, -Min[300, StringLength[stderr]]]];
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
        declaredNumeratorPositions, autoNumeratorPositions,
        autoDetectRestrictions, startFile, fireVariables,
        runtimeFingerprintRecord, setupFingerprintRecord, setupTopology,
        numericalPoint, fixedProps, fixedRepls,
        prepareIBPResult, prepareResult, saveResult, oldDir,
        fireInternal, fireExternal, dimensionVariable},
  name = topology["Name"];

  (* Determine working directory *)
  dir = If[workDir === "",
    FileNameJoin[{FeynmanTrick`Private`$FTConfig["WorkDirectory"], name}],
    workDir
  ];
  If[!DirectoryQ[dir], CreateDirectory[dir, CreateIntermediateDirectories -> True]];

  If[ensureFIRELoaded[] === $Failed, Return[$Failed, Module]];

  If[!MemberQ[{"Modular", "Classical"},
      fireConfigValue["FIREBackend", "Modular"]],
    Print["Error: FIREBackend must be \"Modular\" or \"Classical\"."];
    Return[$Failed, Module]];

  (* Assign problem number *)
  pn = $NextProblemNumber++;

  numericalPoint = Lookup[topology, "NumericalPoint", {}];
  If[!validFIRENumericalPointQ[numericalPoint],
    Print["Error: NumericalPoint must be finite exact assignments from unique symbols."];
    Return[$Failed, Module]];
  If[Intersection[First /@ numericalPoint,
      Join[topology["LoopMomenta"], topology["ExternalMomenta"]]] =!= {},
    Print["Error: NumericalPoint cannot assign loop or external momentum symbols."];
    Return[$Failed, Module]];
  fixedProps = applyFIRENumericalPoint[
    topology["Propagators"], numericalPoint];
  fixedRepls = applyFIRENumericalPointToReplacements[
    topology["Replacements"], numericalPoint];
  If[fixedProps === $Failed || fixedRepls === $Failed,
    Print["Error: NumericalPoint substitutions are cyclic or do not reach a fixed point."];
    Return[$Failed, Module]];
  setupTopology = topology;
  setupTopology["Propagators"] = fixedProps;
  setupTopology["Replacements"] = fixedRepls;

  (* Substitute contexted symbols with Global equivalents for FIRE/Fermat. *)
  fireSubst = buildFIRESubstitution[setupTopology];
  If[fireSubst === $Failed, Return[$Failed, Module]];
  fireProps = setupTopology["Propagators"] /. fireSubst;
  fireRepls = setupTopology["Replacements"] /. fireSubst;
  fireInternal = topology["LoopMomenta"] /. fireSubst;
  fireExternal = topology["ExternalMomenta"] /. fireSubst;
  If[!exactFIREInputQ[{fireProps, fireRepls}],
    Print["Error: FIRE setup requires exact propagators and kinematics; rationalize inexact numerical input explicitly."];
    Return[$Failed, Module]];
  originalFirePropCount = Length[fireProps];
  declaredNumeratorPositions = Sort[DeleteDuplicates[Select[
    Lookup[topology, "NumeratorPositions", {}],
    IntegerQ[#] && 1 <= # <= originalFirePropCount &]]];
  If[Intersection[declaredNumeratorPositions,
      topologyEliminatedPositions[topology]] =!= {},
    Print["Error: a FIRE slot cannot be both eliminated and a declared numerator."];
    Return[$Failed, Module]];

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
  FIRE`Internal = fireInternal;
  FIRE`External = fireExternal;
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
  autoNumeratorPositions =
    Range[originalFirePropCount + 1, Length[fireProps]];
  numeratorPositions = Sort[DeleteDuplicates[
    Join[declaredNumeratorPositions, autoNumeratorPositions]]];
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
  prepareIBPResult = Quiet[Check[FIRE`PrepareIBP[], $Failed]];
  If[prepareIBPResult =!= True,
    Print["Error: FIRE 7 PrepareIBP failed for topology ", name, "."];
    Return[$Failed, Module]];

  dimensionVariable = FeynmanTrick`Private`$FTConfig["DimensionVariable"];
  prepareResult = Quiet[Check[
    FIRE`Prepare[
      "AutoDetectRestrictions" -> autoDetectRestrictions,
      "DVar" -> dimensionVariable
    ]; True,
    $Failed]];
  If[prepareResult =!= True ||
      !IntegerQ[Quiet[Check[FIRE`ExampleDimension[0], $Failed]]],
    Print["Error: FIRE 7 Prepare failed for topology ", name, "."];
    Return[$Failed, Module]];

  (* Save start file *)
  oldDir = Directory[];
  saveResult = Quiet[Check[
    SetDirectory[dir];
    FIRE`SaveStart[name];
    True,
    $Failed]];
  SetDirectory[oldDir];
  startFile = FileNameJoin[{dir, name <> ".start"}];
  If[saveResult =!= True || !FileExistsQ[startFile] ||
      Quiet[Check[FileByteCount[startFile], 0]] <= 0,
    Print["Error: FIRE 7 did not create a nonempty start file for ", name, "."];
    Return[$Failed, Module]];

  (* Return updated topology *)
  result = topology;
  result["WorkDirectory"] = dir;
  result["ProblemNumber"] = pn;
  result["StartFileReady"] = True;
  result["OriginalNumPropagators"] = Lookup[
    topology, "OriginalNumPropagators", topology["NumPropagators"]];
  result["OriginalPropagators"] = Lookup[
    topology, "OriginalPropagators", topology["Propagators"]];
  result["OriginalLoopMomenta"] = Lookup[
    topology, "OriginalLoopMomenta", topology["LoopMomenta"]];
  result["OriginalExternalMomenta"] = Lookup[
    topology, "OriginalExternalMomenta", topology["ExternalMomenta"]];
  result["OriginalReplacements"] = Lookup[
    topology, "OriginalReplacements", topology["Replacements"]];
  result["Propagators"] = fireProps;
  result["LoopMomenta"] = fireInternal;
  result["ExternalMomenta"] = fireExternal;
  result["Replacements"] = fireRepls;
  result["NumPropagators"] = Length[fireProps];
  result["NumeratorPositions"] = numeratorPositions;
  fireVariables = extractVariables[
    fireProps, fireRepls, fireInternal, fireExternal];
  If[fireVariables === $Failed ||
      !AllTrue[SymbolName /@ Append[fireVariables, dimensionVariable],
        safeFIREVariableNameQ],
    Print["Error: topology contains a variable name unsupported by FIRE 7."];
    Return[$Failed, Module]];
  result["Variables"] = fireVariables;

  (* Content-address the exact prepared problem.  The full record is retained
     and used in cache keys; the digest is only a compact diagnostic. *)
  If[ensureFIRELoaded[] === $Failed,
    Print["Error: FIRE7.m changed during topology preparation."];
    Return[$Failed, Module]];
  runtimeFingerprintRecord = currentFIRERuntimeFingerprintRecord[];
  setupFingerprintRecord = Join[<|
    "Schema" -> $fireSetupFingerprintSchema,
    "VerifiedStartFile" -> True,
    "StartFileSHA256" -> fileSHA256[startFile],
    "AutoDetectRestrictions" -> autoDetectRestrictions
  |>, runtimeFingerprintRecord, <|
    "LoopMomenta" -> fireInternal,
    "ExternalMomenta" -> fireExternal,
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
Module[{results, cleared},
  (* Clear state before batch setup *)
  cleared = ClearFIREState[];
  If[cleared === $Failed, Return[$Failed, Module]];

  (* Set up each topology *)
  results = SetupFIRE[#, workDir] & /@ topologies;
  If[MemberQ[results, $Failed], $Failed, results]
];


(* ============================================================ *)
(* writeMultiProblemConfig - config file for multiple problems  *)
(* ============================================================ *)

safeFIREVariableNameQ[name_String] :=
  StringMatchQ[name, RegularExpression["[a-z][a-z0-9]*"]];

orderedFIREVariables[topologies_List] := Module[
  {dimension, variables, ordered, names, collisions},
  dimension = FeynmanTrick`Private`$FTConfig["DimensionVariable"];
  variables = DeleteDuplicates[
    Flatten[Lookup[topologies, "Variables", {}]]];
  If[!AllTrue[Append[variables, dimension], Head[#] === Symbol &] ||
      dimension =!= globalFIRESymbol[dimension],
    Print["Error: FIRE variables must be symbols and the dimension variable must be in Global`."];
    Return[$Failed, Module]];
  ordered = DeleteDuplicates[
    Append[DeleteCases[variables, dimension], dimension]];
  names = SymbolName /@ ordered;
  collisions = Select[GatherBy[ordered, SymbolName], Length[#] > 1 &];
  If[collisions =!= {},
    Print["Error: FIRE variables have colliding basenames: ",
      InputForm[collisions]];
    Return[$Failed, Module]];
  If[!AllTrue[names, safeFIREVariableNameQ],
    Print["Error: FIRE 7 variables require lowercase alphanumeric names: ",
      Select[names, !safeFIREVariableNameQ[#] &]];
    Return[$Failed, Module]];
  (* FIRE 7.1's modular reconstruction stores monomial exponents in fixed
     exp[16] arrays.  Classical execution does not use that code path. *)
  If[fireConfigValue["FIREBackend", "Modular"] === "Modular" &&
      Length[ordered] > 16,
    Print["Error: FIRE 7 modular reconstruction supports at most 16 variables including the dimension; received ",
      Length[ordered], "."];
    Return[$Failed, Module]];
  ordered
];

writeMultiProblemConfig[topologies_List, dir_String, configName_String,
                        integralsFile_String, outputFile_String] :=
Module[{content, threads, fthreads, problemLines, allVars, path},
  threads = FeynmanTrick`Private`$FTConfig["Threads"];
  fthreads = FeynmanTrick`Private`$FTConfig["FThreads"];

  (* FIRE7 reconstruction requires kinematics first and the dimension last. *)
  allVars = orderedFIREVariables[topologies];
  If[allVars === $Failed, Return[$Failed, Module]];

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

  path = FileNameJoin[{dir, configName}];
  If[Export[path, content, "Text"] === $Failed,
    $Failed,
    <|"Path" -> path, "OrderedVariables" -> allVars,
      "OutputFile" -> outputFile, "IntegralsFile" -> integralsFile|>]
];


(* Write single-problem config (for backwards compatibility) *)
writeSingleProblemConfig[topology_Association, dir_String, configName_String,
                         integralsFile_String:"", outputFile_String:""] :=
Module[{vars, content, threads, fthreads, intFile, outFile, name, pn, path,
        masters, mastersFile, mastersLine, mastersContent},
  name = topology["Name"];
  pn = topology["ProblemNumber"];
  threads = FeynmanTrick`Private`$FTConfig["Threads"];
  fthreads = FeynmanTrick`Private`$FTConfig["FThreads"];

  vars = orderedFIREVariables[{topology}];
  If[vars === $Failed, Return[$Failed, Module]];

  intFile = If[integralsFile === "", name <> ".m", integralsFile];
  outFile = If[outputFile === "", name <> ".tables", outputFile];
  masters = Lookup[topology, "Masters", {}];
  If[!ListQ[masters] ||
      !AllTrue[masters,
        ListQ[#] && Length[#] === topology["NumPropagators"] &&
          AllTrue[#, IntegerQ] &],
    Return[$Failed, Module]];
  mastersFile = If[masters === {}, "",
    FileBaseName[configName] <> ".preferred"];
  mastersLine = If[masters === {}, "",
    "#preferred         " <> mastersFile <> "\n"];
  If[masters =!= {},
    mastersContent = ToString[
      ({pn, #} & /@ masters), InputForm] <> "\n";
    If[Export[FileNameJoin[{dir, mastersFile}], mastersContent, "Text"] ===
        $Failed,
      Return[$Failed, Module]]
  ];

  content = StringJoin[
    "#threads           ", ToString[threads], "\n",
    "#fthreads          ", ToString[fthreads], "\n",
    "#variables         ", StringRiffle[SymbolName /@ vars, ","], "\n",
    "#start\n",
    "#problem           ", ToString[pn], " ", name, ".start\n",
    mastersLine,
    "#integrals         ", intFile, "\n",
    "#output            ", outFile, "\n"
  ];

  path = FileNameJoin[{dir, configName}];
  If[Export[path, content, "Text"] === $Failed,
    $Failed,
    <|"Path" -> path, "OrderedVariables" -> vars,
      "OutputFile" -> outFile, "IntegralsFile" -> intFile,
      "MasterFile" -> If[mastersFile === "", Missing["NotUsed"],
        mastersFile]|>]
];

fireInvocationStem[prefix_String, kind_String, identity_] :=
  prefix <> "_" <> StringTake[
    IntegerString[Hash[{"FeynmanTrick.FIREInvocation/v1", identity},
      "SHA256"], 16, 64], 16] <> "_" <> kind;


(* ============================================================ *)
(* FindBasis - find master integrals for a single topology      *)
(* ============================================================ *)

FindBasis[topology_Association] :=
Module[{dir, name, fireBin, intFile, result, masters, tablesFile, pn,
        allSectors, intContent, configManifest, parsedMasters, parsedCoverage,
        requestedPairs, runStem, expectedInputs, expectedRuntime},

  If[!TrueQ[Lookup[topology, "StartFileReady", False]],
    Print["Error: Must call SetupFIRE before FindBasis."];
    Return[$Failed, Module];
  ];
  If[!TrueQ[preparedTopologyCompatibleQ[topology]], Return[$Failed, Module]];
  If[ensureFIRELoaded[] === $Failed, Return[$Failed, Module]];

  dir = topology["WorkDirectory"];
  name = topology["Name"];
  pn = topology["ProblemNumber"];
  fireBin = FileNameJoin[{FeynmanTrick`Private`$FTConfig["FIREPath"], "bin", "FIRE7"}];

  (* Reduce the corner plus first-complexity shell of every allowed sector. *)
  allSectors = basisSeedIntegrals[topology];
  If[allSectors === {},
    Print["Error: topology has no active denominator sectors."];
    Return[$Failed, Module]];
  requestedPairs = {pn, #} & /@ allSectors;
  intContent = StringJoin["{",
    StringRiffle[("{" <> ToString[pn] <> "," <>
      ToString[#, InputForm] <> "}") & /@ allSectors, ",\n"], "}\n"];
  runStem = fireInvocationStem[name, "basis", {
    Lookup[topology, "SetupFingerprint", Missing["Absent"]], allSectors}];
  intFile = FileNameJoin[{dir, runStem <> ".m"}];
  If[Export[intFile, intContent, "Text"] === $Failed,
    Return[$Failed, Module]];

  configManifest = writeSingleProblemConfig[
    topology, dir, runStem <> ".config", runStem <> ".m",
    runStem <> ".tables"];
  If[configManifest === $Failed, Return[$Failed, Module]];
  expectedInputs = preparedStartHashes[{topology}];
  expectedRuntime = preparedRunnerRuntimeHashes[{topology}, "Basis"];
  If[expectedInputs === $Failed || expectedRuntime === $Failed,
    Return[$Failed, Module]];

  (* Run FIRE7 *)
  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["Running FIRE7 to find basis for ", name, "..."];
  ];
  If[runFIRE6[fireBin, dir, runStem, expectedInputs, "Basis",
      expectedRuntime] =!= 0,
    Print["Error: FIRE7 basis discovery failed."];
    Return[$Failed, Module];
  ];
  If[!TrueQ[preparedTopologyCompatibleQ[topology]],
    Return[$Failed, Module]];

  (* Load tables and extract masters *)
  tablesFile = FileNameJoin[{dir, runStem <> ".tables"}];
  If[!FileExistsQ[tablesFile],
    Print["Error: FIRE7 did not produce output tables."];
    Return[$Failed, Module];
  ];

  parsedMasters = canonicalMasterList[tablesFile];
  parsedCoverage = rawFIRETableCoverage[tablesFile];
  If[parsedMasters === $Failed ||
      parsedCoverage === $Failed ||
      !SubsetQ[parsedCoverage, requestedPairs] ||
      !validFIREPairsForRequestsQ[parsedMasters, requestedPairs] ||
      !AllTrue[Join[parsedMasters, parsedCoverage], First[#] === pn &],
    Print["Error: FIRE 7 basis table is malformed, incomplete, or contains a foreign problem."];
    Return[$Failed, Module]];
  masters = Cases[parsedMasters, {pn, indices_List} :> indices];

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
Module[{dir, fireBin, intContent, allSectors, tablesFile, masters, results,
        configManifest, problemNumbers, runStem, parsedCoverage,
        requestedPairs, expectedInputs, expectedRuntime},

  If[Length[topologies] == 0, Return[{}, Module]];
  If[!AllTrue[topologies,
      TrueQ[Lookup[#, "StartFileReady", False]] &],
    Print["Error: Must call SetupFIREBatch before FindBasisBatch."];
    Return[$Failed, Module]];
  If[!AllTrue[topologies, TrueQ[preparedTopologyCompatibleQ[#]] &],
    Return[$Failed, Module]];
  If[ensureFIRELoaded[] === $Failed, Return[$Failed, Module]];

  (* Use first topology's work directory *)
  dir = topologies[[1]]["WorkDirectory"];
  If[!AllTrue[topologies, #["WorkDirectory"] === dir &],
    Print["Error: batch FIRE topologies must share one work directory."];
    Return[$Failed, Module]];
  fireBin = FileNameJoin[{FeynmanTrick`Private`$FTConfig["FIREPath"], "bin", "FIRE7"}];
  problemNumbers = Lookup[topologies, "ProblemNumber"];

  (* Build combined integrals file *)
  intContent = StringJoin["{",
    StringRiffle[
      Flatten[
        Table[
          With[{topo = topologies[[i]], pn = topologies[[i]]["ProblemNumber"]},
            allSectors = basisSeedIntegrals[topo];
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
  runStem = fireInvocationStem["batch", "basis", {
    Lookup[topologies, "SetupFingerprint"], intContent}];
  If[Export[FileNameJoin[{dir, runStem <> ".m"}], intContent, "Text"] === $Failed,
    Return[$Failed, Module]];

  (* Write multi-problem config *)
  configManifest = writeMultiProblemConfig[topologies, dir,
    runStem <> ".config", runStem <> ".m", runStem <> ".tables"];
  If[configManifest === $Failed, Return[$Failed, Module]];
  expectedInputs = preparedStartHashes[topologies];
  expectedRuntime = preparedRunnerRuntimeHashes[topologies, "Basis"];
  If[expectedInputs === $Failed || expectedRuntime === $Failed,
    Return[$Failed, Module]];

  (* Run FIRE7 once *)
  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["Running FIRE7 to find basis for ", Length[topologies], " topologies..."];
  ];
  If[runFIRE6[fireBin, dir, runStem, expectedInputs, "Basis",
      expectedRuntime] =!= 0,
    Print["Error: FIRE7 batch basis discovery failed."];
    Return[$Failed, Module];
  ];
  If[!AllTrue[topologies, TrueQ[preparedTopologyCompatibleQ[#]] &],
    Return[$Failed, Module]];

  (* Load results *)
  tablesFile = FileNameJoin[{dir, runStem <> ".tables"}];
  If[!FileExistsQ[tablesFile],
    Print["Error: FIRE7 did not produce output tables."];
    Return[$Failed, Module];
  ];

  requestedPairs = Quiet[Check[Get[FileNameJoin[{dir, runStem <> ".m"}]],
    $Failed]];
  masters = canonicalMasterList[tablesFile];
  parsedCoverage = rawFIRETableCoverage[tablesFile];
  If[!ListQ[requestedPairs] || masters === $Failed ||
      parsedCoverage === $Failed ||
      !SubsetQ[parsedCoverage, requestedPairs] ||
      !validFIREPairsForRequestsQ[masters, requestedPairs] ||
      !AllTrue[Join[masters, parsedCoverage],
        MemberQ[problemNumbers, First[#]] &],
    Print["Error: FIRE 7 batch basis table is malformed, incomplete, or contains a foreign problem."];
    Return[$Failed, Module]];

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
        computedReductions, allEntries, allMasters, cacheEnabled,
        configManifest, parsedTable, knownMasters, runStem, expectedInputs,
        expectedRuntime},

  If[!TrueQ[Lookup[topology, "StartFileReady", False]],
    Print["Error: Must call SetupFIRE before ReduceIntegrals."];
    Return[$Failed, Module];
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
  knownMasters = Lookup[topology, "Masters", {}];
  (* A prepared master basis is already an exact set of identity reductions.
     Seed those identities before looking for misses so a boundary-observable
     batch never asks FIRE to rediscover its own transported masters. *)
  If[cacheEnabled && knownMasters =!= {},
    Scan[
      cacheReduction[
        topology, #, Global`G[1, #], knownMasters] &,
      DeleteDuplicates[knownMasters]]
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

  If[ensureFIRELoaded[] === $Failed, Return[$Failed, Module]];

  dir = topology["WorkDirectory"];
  name = topology["Name"];
  pn = topology["ProblemNumber"];
  fireBin = FileNameJoin[{FeynmanTrick`Private`$FTConfig["FIREPath"], "bin", "FIRE7"}];

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
  runStem = fireInvocationStem[name, "reduce", {
    Lookup[topology, "SetupFingerprint", Missing["Absent"]],
    Lookup[topology, "Masters", {}],
    missingFireIntegrals}];
  intFile = runStem <> ".m";
  If[Export[FileNameJoin[{dir, intFile}], intContent, "Text"] === $Failed,
    Return[$Failed, Module]];

  (* Write config *)
  configManifest = writeSingleProblemConfig[topology, dir,
    runStem <> ".config", intFile, runStem <> ".tables"];
  If[configManifest === $Failed, Return[$Failed, Module]];
  expectedInputs = preparedStartHashes[{topology}];
  expectedRuntime = preparedRunnerRuntimeHashes[{topology}, "Reduction"];
  If[expectedInputs === $Failed || expectedRuntime === $Failed,
    Return[$Failed, Module]];

  (* Run FIRE7 *)
  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 2,
    Print["Reducing ", Length[missingFireIntegrals], " integrals..."];
  ];
  If[runFIRE6[fireBin, dir, runStem, expectedInputs, "Reduction",
      expectedRuntime] =!= 0,
    Print["Error: FIRE7 integral reduction failed."];
    Return[$Failed, Module];
  ];
  If[!TrueQ[preparedTopologyCompatibleQ[topology]],
    Return[$Failed, Module]];

  (* Load results *)
  tablesFile = FileNameJoin[{dir, runStem <> ".tables"}];
  If[!FileExistsQ[tablesFile],
    Print["Error: FIRE7 reduction did not produce output."];
    Return[$Failed, Module];
  ];

  parsedTable = parseAndValidateFIRETable[tablesFile,
    ({pn, #} & /@ missingFireIntegrals), {pn}];
  If[parsedTable === $Failed,
    Print["Error: FIRE 7 reduction table is malformed or incomplete."];
    Return[$Failed, Module]];
  rules = parsedTable["Rules"];
  masters = Cases[parsedTable["Masters"],
    {pn, indices_List} :> indices
  ];
  If[knownMasters =!= {} && !SubsetQ[knownMasters, masters],
    Print["Error: FIRE 7 reduction reported a master outside the prepared basis."];
    Return[$Failed, Module]];

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
Module[{dir, fireBin, intContent, tablesFile, rules, results, allIntegrals,
        topologies, configManifest, parsedTable, problemNumbers,
        reportedMasters, knownMasters, runStem, expectedInputs,
        expectedRuntime},

  If[Length[topoIntegralPairs] == 0, Return[{}, Module]];
  If[!AllTrue[topoIntegralPairs,
      MatchQ[#, {_Association, _List}] &], Return[$Failed, Module]];
  topologies = topoIntegralPairs[[All, 1]];
  If[!AllTrue[topologies,
      TrueQ[Lookup[#, "StartFileReady", False]] &],
    Print["Error: Must prepare every topology before ReduceIntegralsBatch."];
    Return[$Failed, Module]];
  If[!AllTrue[topologies, TrueQ[preparedTopologyCompatibleQ[#]] &],
    Return[$Failed, Module]];
  If[ensureFIRELoaded[] === $Failed, Return[$Failed, Module]];

  (* Use first topology's work directory *)
  dir = topoIntegralPairs[[1, 1]]["WorkDirectory"];
  If[!AllTrue[topologies, #["WorkDirectory"] === dir &],
    Print["Error: batch FIRE topologies must share one work directory."];
    Return[$Failed, Module]];
  fireBin = FileNameJoin[{FeynmanTrick`Private`$FTConfig["FIREPath"], "bin", "FIRE7"}];
  problemNumbers = Lookup[topologies, "ProblemNumber"];

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
  If[allIntegrals === {},
    Return[ConstantArray[<||>, Length[topologies]], Module]];

  intContent = StringJoin[
    "{",
    StringRiffle[
      ("{" <> ToString[#[[1]]] <> "," <> ToString[#[[3]], InputForm] <> "}") & /@ allIntegrals,
      ",\n"
    ],
    "}\n"
  ];
  runStem = fireInvocationStem["batch", "reduce", {
    Lookup[topologies, "SetupFingerprint"], allIntegrals}];
  If[Export[FileNameJoin[{dir, runStem <> ".m"}], intContent, "Text"] === $Failed,
    Return[$Failed, Module]];

  (* Write multi-problem config *)
  configManifest = writeMultiProblemConfig[topologies, dir,
    runStem <> ".config", runStem <> ".m", runStem <> ".tables"];
  If[configManifest === $Failed, Return[$Failed, Module]];
  expectedInputs = preparedStartHashes[topologies];
  expectedRuntime = preparedRunnerRuntimeHashes[topologies, "Reduction"];
  If[expectedInputs === $Failed || expectedRuntime === $Failed,
    Return[$Failed, Module]];

  (* Run FIRE7 *)
  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["Reducing integrals for ", Length[topoIntegralPairs], " topologies..."];
  ];
  If[runFIRE6[fireBin, dir, runStem, expectedInputs, "Reduction",
      expectedRuntime] =!= 0,
    Print["Error: FIRE7 batch integral reduction failed."];
    Return[$Failed, Module];
  ];
  If[!AllTrue[topologies, TrueQ[preparedTopologyCompatibleQ[#]] &],
    Return[$Failed, Module]];

  (* Load results *)
  tablesFile = FileNameJoin[{dir, runStem <> ".tables"}];
  If[!FileExistsQ[tablesFile],
    Print["Error: FIRE7 reduction did not produce output."];
    Return[$Failed, Module];
  ];

  parsedTable = parseAndValidateFIRETable[tablesFile,
    ({#[[1]], #[[3]]} & /@ allIntegrals), problemNumbers];
  If[parsedTable === $Failed,
    Print["Error: FIRE 7 batch reduction table is malformed or incomplete."];
    Return[$Failed, Module]];
  rules = parsedTable["Rules"];
  reportedMasters = parsedTable["Masters"];
  Do[
    knownMasters = Lookup[topologies[[i]], "Masters", {}];
    If[knownMasters =!= {} && !SubsetQ[knownMasters,
        Cases[reportedMasters,
          {p_, idx_List} /; p === problemNumbers[[i]] :> idx]],
      Print["Error: FIRE 7 batch reduction reported an unexpected master."];
      Return[$Failed, Module]],
    {i, Length[topologies]}];

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
