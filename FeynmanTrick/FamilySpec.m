(* ::Package:: *)
(* Exact, process-free input normalization for user-defined FT families. *)

BeginPackage["FeynmanTrick`FamilySpec`", {
  "FeynmanTrick`",
  "FeynmanTrick`FIREInterface`"
}];

CreateFamily::usage =
  "CreateFamily[family, opts] validates an exact raw family or an unprepared " <>
  "FIRE topology and returns a canonical FeynmanTrick family specification. " <>
  "CreateFamily[family, integrals, opts] sets the default output integrals; " <>
  "integrals may be one integer index vector, an ordered list of vectors, " <>
  "or All. All is recorded for later master discovery and is not discovered " <>
  "by CreateFamily.";

NormalizeOutputIntegrals::usage =
  "NormalizeOutputIntegrals[spec, n, eliminated] converts Automatic, one index vector, " <>
  "an ordered list of vectors, or All into the canonical output-target form " <>
  "for an n-propagator family and enforces zero indices at eliminated positions.";

ValidateOutputIntegralsForSequence::usage =
  "ValidateOutputIntegralsForSequence[integrals, sequence, eliminated] verifies that " <>
  "every explicit target is supported by each ordered Feynman-trick merge. " <>
  "Negative indices are allowed only at positions that never participate " <>
  "in the sequence, and eliminated positions cannot be merged. The third " <>
  "argument defaults to an empty list.";

CreateFamily::input =
  "The first argument must be a raw family association or an existing topology association.";
CreateFamily::missing =
  "The family association is missing required key(s) `1`.";
CreateFamily::name =
  "Family name `1` must be a FIRE-safe identifier beginning with a letter and containing only ASCII letters, digits, and underscores.";
CreateFamily::field =
  "Family field `1` has an unsupported shape: `2`.";
CreateFamily::momenta =
  "LoopMomenta and ExternalMomenta must be duplicate-free, disjoint lists of symbols; received `1` and `2`.";
CreateFamily::rules =
  "Family field `1` must be a list of exact immediate Rules with unique left-hand sides.";
CreateFamily::inexact =
  "Family field `1` contains an inexact number. Family definitions must be exact.";
CreateFamily::point =
  "NumericalPoint must contain finite exact assignments from non-momentum symbols; received `1`.";
CreateFamily::nprop =
  "NumPropagators `1` does not agree with the propagator-list length `2`.";
CreateFamily::dimension =
  "Dimension `1` must resolve to an exact expression. Supply an exact Dimension field or configure an exact DimensionExpression.";
CreateFamily::eliminated =
  "EliminatedPositions `1` must be a duplicate-free list of integer positions between 1 and `2`.";
CreateFamily::sequence =
  "CombinationSequence `1` must be Automatic or a list of distinct active-position pairs.";
CreateFamily::pair =
  "Combination pair `1` at step `2` is invalid for the currently active positions `3`.";
CreateFamily::prepared =
  "Prepared FIRE topology state cannot be ingested as a family definition. Pass the original unprepared DefineTopology association instead.";
CreateFamily::canonical =
  "Canonical FamilySpec input contains inconsistent duplicated mathematical data in Definition, Topology, or TopTopology.";

NormalizeOutputIntegrals::shape =
  "OutputIntegrals `1` must be Automatic, All, one integer index vector, or a nonempty list of integer index vectors.";
NormalizeOutputIntegrals::length =
  "Output integral `1` has length `2`; this family requires exactly `3` indices.";
NormalizeOutputIntegrals::integer =
  "Output integral `1` is unsupported: every index must be an exact integer.";
NormalizeOutputIntegrals::nprop =
  "The propagator count `1` must be a positive integer.";
NormalizeOutputIntegrals::eliminated =
  "EliminatedPositions `1` must be a duplicate-free list of integer positions between 1 and `2`.";
NormalizeOutputIntegrals::sector =
  "Output integral `1` is incompatible with EliminatedPositions `2`; every eliminated index must be zero.";
ValidateOutputIntegralsForSequence::merge =
  "Output integral `1` is unsupported at merge step `2` (`3`): merged powers `4` must both be nonnegative integers. Negative indices are supported only at positions that never participate in a merge.";
ValidateOutputIntegralsForSequence::sequence =
  "CombinationSequence `1` is not a list of integer position pairs valid for output vectors of length `2`.";

Begin["`Private`"];

Options[CreateFamily] = {
  "Name" -> Automatic,
  "CombinationSequence" -> Automatic,
  "OutputIntegrals" -> Automatic
};

requiredFamilyKeys = {
  "LoopMomenta", "ExternalMomenta", "Propagators", "Replacements"
};

mathematicalOptionalKeys = {
  "AnalyticPrescription", "Prescriptions", "KinematicAssumptions"
};

inexactNumberPresentQ[expr_] := !FreeQ[Unevaluated[expr], _Real];

validRuleListQ[rules_] :=
  ListQ[rules] && AllTrue[rules, MatchQ[#, _Rule] &] &&
  DuplicateFreeQ[First /@ rules];

validNumericalPointRuleListQ[rules_] :=
  validRuleListQ[rules] && AllTrue[First /@ rules, Head[#] === Symbol &];

validDimensionQ[dimension_] :=
  dimension =!= Automatic &&
  !inexactNumberPresentQ[dimension] &&
  FreeQ[dimension,
    _Missing | Indeterminate | _DirectedInfinity | ComplexInfinity] &&
  !MatchQ[dimension,
    _List | _Association | _Rule | _RuleDelayed | _String |
      True | False | Null];

validFamilyNameQ[name_] :=
  StringQ[name] &&
  StringMatchQ[name, RegularExpression["[A-Za-z][A-Za-z0-9_]*"]];

stableDigest[expr_] :=
  IntegerString[Hash[expr, "SHA256"], 16, 64];

automaticLeftToRightSequence[n_Integer, eliminated_List] := Module[
  {active = Complement[Range[n], eliminated]},
  If[Length[active] <= 1, {},
    Thread[{ConstantArray[First[active], Length[active] - 1], Rest[active]}]
  ]
];

normalizeCombinationSequence[Automatic, n_Integer, eliminated_List] :=
  automaticLeftToRightSequence[n, eliminated];

normalizeCombinationSequence[sequence_, n_Integer,
    eliminated_List] := Module[
  {active = Complement[Range[n], eliminated], pair, step},
  If[!ListQ[sequence] || Length[sequence] > Max[0, Length[active] - 1] ||
      !AllTrue[sequence, MatchQ[#, {_Integer, _Integer}] &],
    Message[CreateFamily::sequence, sequence];
    Return[$Failed, Module]
  ];
  Do[
    pair = sequence[[step]];
    If[pair[[1]] === pair[[2]] ||
        !MemberQ[active, pair[[1]]] || !MemberQ[active, pair[[2]]],
      Message[CreateFamily::pair, pair, step, active];
      Return[$Failed, Module]
    ];
    active = DeleteCases[active, pair[[2]]],
    {step, Length[sequence]}
  ];
  sequence
];

validateTargetVector[vector_List, n_Integer] := Module[{},
  If[Length[vector] =!= n,
    Message[NormalizeOutputIntegrals::length, vector, Length[vector], n];
    Return[$Failed, Module]
  ];
  If[!AllTrue[vector, IntegerQ],
    Message[NormalizeOutputIntegrals::integer, vector];
    Return[$Failed, Module]
  ];
  vector
];

validEliminatedPositionsQ[positions_, n_Integer] :=
  ListQ[positions] && DuplicateFreeQ[positions] &&
  AllTrue[positions, IntegerQ[#] && 1 <= # <= n &];

canonicalDuplicatedMathematicsQ[input_Association] := Module[
  {definition, topology, topTopology, keys},
  If[Lookup[input, "Schema", None] =!= "FeynmanTrick.FamilySpec/v1",
    Return[True, Module]];
  definition = Lookup[input, "Definition", None];
  topology = Lookup[input, "Topology", None];
  topTopology = Lookup[input, "TopTopology", None];
  keys = Join[
    requiredFamilyKeys,
    {"NumericalPoint", "Dimension", "EliminatedPositions"},
    mathematicalOptionalKeys
  ];
  AssociationQ[definition] && AssociationQ[topology] &&
    AssociationQ[topTopology] && topTopology === topology &&
    AllTrue[keys,
      Lookup[topology, #, Missing["Absent"]] ===
          Lookup[definition, #, Missing["Absent"]] &]
];

NormalizeOutputIntegrals[spec_, n_, eliminatedPositions_:{}] := Module[
  {vectors, checked, result, incompatible},
  If[!IntegerQ[n] || n <= 0,
    Message[NormalizeOutputIntegrals::nprop, n];
    Return[$Failed, Module]
  ];
  If[!validEliminatedPositionsQ[eliminatedPositions, n],
    Message[NormalizeOutputIntegrals::eliminated, eliminatedPositions, n];
    Return[$Failed, Module]
  ];
  result = Which[
    spec === Automatic,
      {ReplacePart[ConstantArray[1, n],
        Thread[eliminatedPositions -> 0]]},
    spec === All,
      All,
    ListQ[spec] && Length[spec] === n && !AnyTrue[spec, ListQ],
      checked = validateTargetVector[spec, n];
      If[checked === $Failed, $Failed, {checked}],
    ListQ[spec] && spec =!= {} && AllTrue[spec, ListQ],
      vectors = validateTargetVector[#, n] & /@ spec;
      If[MemberQ[vectors, $Failed], $Failed, vectors],
    ListQ[spec] && !AnyTrue[spec, ListQ],
      (* A flat list is an attempted single vector even when its length is
         wrong.  Give the useful arity diagnostic instead of a shape error. *)
      checked = validateTargetVector[spec, n];
      If[checked === $Failed, $Failed, {checked}],
    True,
      Message[NormalizeOutputIntegrals::shape, spec];
      $Failed
  ];
  If[result === $Failed || result === All, Return[result, Module]];
  incompatible = Select[
    result,
    Function[vector,
      AnyTrue[eliminatedPositions,
        Function[position, vector[[position]] =!= 0]]]
  ];
  If[incompatible =!= {},
    Message[NormalizeOutputIntegrals::sector,
      First[incompatible], eliminatedPositions];
    Return[$Failed, Module]
  ];
  result
];

ValidateOutputIntegralsForSequence[All, _List, _List:{}] := All;

ValidateOutputIntegralsForSequence[integrals_List,
    sequence_List, eliminatedPositions_List:{}] := Module[
  {n, result = integrals, vector, pair, step, vi, vj, target, active},
  If[integrals === {}, Return[integrals, Module]];
  If[!AllTrue[integrals, ListQ] ||
      !SameQ @@ (Length /@ integrals),
    Message[ValidateOutputIntegralsForSequence::sequence,
      sequence, Missing["InconsistentTargetLengths"]];
    Return[$Failed, Module]
  ];
  n = Length[First[integrals]];
  If[!validEliminatedPositionsQ[eliminatedPositions, n],
    Message[NormalizeOutputIntegrals::eliminated, eliminatedPositions, n];
    Return[$Failed, Module]
  ];
  If[!AllTrue[sequence,
      MatchQ[#, {_Integer, _Integer}] &&
      #[[1]] =!= #[[2]] &&
      1 <= #[[1]] <= n && 1 <= #[[2]] <= n &],
    Message[ValidateOutputIntegralsForSequence::sequence, sequence, n];
    Return[$Failed, Module]
  ];
  (* A syntactically valid pair can still be impossible after an earlier
     merge removed its second position.  Validate the ordered active-position
     evolution here as well as in CreateFamily so direct DefineFTIteration
     callers cannot bypass the canonical sequence contract. *)
  active = Complement[Range[n], eliminatedPositions];
  Do[
    pair = sequence[[step]];
    If[!MemberQ[active, pair[[1]]] || !MemberQ[active, pair[[2]]],
      Message[ValidateOutputIntegralsForSequence::sequence, sequence, n];
      Return[$Failed, Module]
    ];
    active = DeleteCases[active, pair[[2]]],
    {step, Length[sequence]}
  ];
  Do[
    vector = result[[target]];
    Do[
      pair = sequence[[step]];
      {vi, vj} = vector[[pair]];
      If[!IntegerQ[vi] || !IntegerQ[vj] || vi < 0 || vj < 0,
        Message[ValidateOutputIntegralsForSequence::merge,
          integrals[[target]], step, pair, {vi, vj}];
        Return[$Failed, Module]
      ];
      (* This is the exact index map used by FeynmanTrickNeededIntegral.
         With nonnegative powers its four branches are exhaustive; the
         direct branch is reached only for the genuine {0,0} case. *)
      Which[
        vi > 0 && vj > 0,
          vector[[pair[[1]]]] = vi + vj,
        vi > 0 && vj === 0,
          vector[[pair[[1]]]] = vi,
        vi === 0 && vj > 0,
          vector[[pair[[1]]]] = vj,
        vi === 0 && vj === 0,
          vector[[pair[[1]]]] = 0
      ];
      vector[[pair[[2]]]] = 0,
      {step, Length[sequence]}
    ];
    result[[target]] = vector,
    {target, Length[result]}
  ];
  integrals
];

ValidateOutputIntegralsForSequence[targets_, sequence_, eliminated_:{}] := (
  Message[ValidateOutputIntegralsForSequence::sequence,
    sequence, Missing["UnknownTargetLength"]];
  $Failed
);

familyDefinition[source_Association, dimension_, numericalPoint_List,
    eliminatedPositions_List] := Association @ Join[
  {
    "LoopMomenta" -> source["LoopMomenta"],
    "ExternalMomenta" -> source["ExternalMomenta"],
    "Propagators" -> source["Propagators"],
    "Replacements" -> source["Replacements"],
    "NumericalPoint" -> numericalPoint,
    "Dimension" -> dimension,
    "EliminatedPositions" -> eliminatedPositions
  },
  Cases[
    mathematicalOptionalKeys,
    key_ /; KeyExistsQ[source, key] :> key -> source[key]
  ]
];

validateFamilySource[source_Association] := Module[
  {missing, loops, externals, propagators, replacements, n, key,
   eliminatedPositions},
  missing = Select[requiredFamilyKeys, !KeyExistsQ[source, #] &];
  If[missing =!= {},
    Message[CreateFamily::missing, missing];
    Return[$Failed, Module]
  ];

  If[Lookup[source, "StartFileReady", False] =!= False ||
      Lookup[source, "ProblemNumber", 0] =!= 0 ||
      Lookup[source, "WorkDirectory", ""] =!= "" ||
      Lookup[source, "Masters", {}] =!= {} ||
      Lookup[source, "MasterRules", {}] =!= {} ||
      AnyTrue[{
        "SetupFingerprint", "SetupFingerprintRecord",
        "NumeratorPositions", "OriginalPropagators",
        "OriginalNumPropagators"
      }, KeyExistsQ[source, #] &],
    Message[CreateFamily::prepared];
    Return[$Failed, Module]
  ];

  loops = source["LoopMomenta"];
  externals = source["ExternalMomenta"];
  propagators = source["Propagators"];
  replacements = source["Replacements"];

  If[!ListQ[loops] || loops === {},
    Message[CreateFamily::field, "LoopMomenta", loops];
    Return[$Failed, Module]
  ];
  If[!ListQ[externals],
    Message[CreateFamily::field, "ExternalMomenta", externals];
    Return[$Failed, Module]
  ];
  If[!ListQ[propagators] || propagators === {},
    Message[CreateFamily::field, "Propagators", propagators];
    Return[$Failed, Module]
  ];
  If[!VectorQ[loops, MatchQ[#, _Symbol] &] ||
      !VectorQ[externals, MatchQ[#, _Symbol] &] ||
      !DuplicateFreeQ[loops] || !DuplicateFreeQ[externals] ||
      Intersection[loops, externals] =!= {},
    Message[CreateFamily::momenta, loops, externals];
    Return[$Failed, Module]
  ];
  If[!validRuleListQ[replacements],
    Message[CreateFamily::rules, "Replacements"];
    Return[$Failed, Module]
  ];

  Do[
    If[inexactNumberPresentQ[source[key]],
      Message[CreateFamily::inexact, key];
      Return[$Failed, Module]
    ],
    {key, Join[requiredFamilyKeys,
      Select[mathematicalOptionalKeys, KeyExistsQ[source, #] &]]}
  ];

  n = Length[propagators];
  If[KeyExistsQ[source, "NumPropagators"] &&
      source["NumPropagators"] =!= n,
    Message[CreateFamily::nprop, source["NumPropagators"], n];
    Return[$Failed, Module]
  ];
  eliminatedPositions = Lookup[source, "EliminatedPositions", {}];
  If[!validEliminatedPositionsQ[eliminatedPositions, n],
    Message[CreateFamily::eliminated, eliminatedPositions, n];
    Return[$Failed, Module]
  ];
  <|
    "NumPropagators" -> n,
    "EliminatedPositions" -> Sort[eliminatedPositions]
  |>
];

resolveSource[input_Association] :=
  If[Lookup[input, "Schema", None] === "FeynmanTrick.FamilySpec/v1" &&
      AssociationQ[Lookup[input, "Topology", None]],
    input["Topology"],
    input
  ];

resolveEmbedded[input_Association, source_Association, key_String,
    optionValue_] := Which[
  optionValue =!= Automatic, optionValue,
  KeyExistsQ[input, key], input[key],
  KeyExistsQ[source, key], source[key],
  True, Automatic
];

requestRecord[familyID_String, vector_List, ordinal_Integer] := Module[
  {identity},
  identity = "ft-l0-" <>
    stableDigest[{"FeynmanTrick.L0OutputRequest/v1", familyID, vector}];
  <|
    "Schema" -> "FeynmanTrick.L0OutputRequest/v1",
    "RequestID" -> identity,
    "RequestIdentity" -> identity,
    "PhysicalIntegralID" -> identity,
    "RequestOrdinal" -> ordinal,
    "IndexVector" -> vector,
    "Resolved" -> True
  |>
];

allRequestRecord[familyID_String] := Module[{identity},
  identity = "ft-l0-all-" <>
    stableDigest[{"FeynmanTrick.L0OutputRequest/v1", familyID, All}];
  <|
    "Schema" -> "FeynmanTrick.L0OutputRequest/v1",
    "RequestID" -> identity,
    "RequestIdentity" -> identity,
    "PhysicalIntegralID" -> identity,
    "RequestOrdinal" -> 1,
    "IndexVector" -> All,
    "Selection" -> All,
    "Resolved" -> False
  |>
];

createFamily[input_Association, positionalTargets_, nameOption_,
    sequenceOption_, outputOption_] := Module[
  {source, inputKind, validated, n, eliminatedPositions,
   numericalPoint, dimensionSpec, dimension, definition, familyID, name,
   sequenceSpec,
   sequence, outputSpec, normalizedTargets, topology, requests,
   outputSource},

  If[!canonicalDuplicatedMathematicsQ[input],
    Message[CreateFamily::canonical];
    Return[$Failed, Module]
  ];
  source = resolveSource[input];
  inputKind = If[KeyExistsQ[source, "NumPropagators"],
    "Topology", "RawFamily"];
  validated = validateFamilySource[source];
  If[validated === $Failed, Return[$Failed, Module]];
  n = validated["NumPropagators"];
  eliminatedPositions = validated["EliminatedPositions"];

  numericalPoint = Lookup[source, "NumericalPoint", {}];
  If[!validNumericalPointRuleListQ[numericalPoint],
    Message[CreateFamily::rules, "NumericalPoint"];
    Return[$Failed, Module]
  ];
  If[inexactNumberPresentQ[numericalPoint],
    Message[CreateFamily::inexact, "NumericalPoint"];
    Return[$Failed, Module]
  ];
  If[!FreeQ[numericalPoint,
      Indeterminate | ComplexInfinity | _DirectedInfinity] ||
      Intersection[First /@ numericalPoint,
        Join[source["LoopMomenta"], source["ExternalMomenta"]]] =!= {},
    Message[CreateFamily::point, numericalPoint];
    Return[$Failed, Module]
  ];
  If[FeynmanTrick`FIREInterface`Private`applyFIRENumericalPoint[
        First /@ numericalPoint, numericalPoint] === $Failed,
    Message[CreateFamily::point, numericalPoint];
    Return[$Failed, Module]
  ];

  dimensionSpec = Lookup[source, "Dimension", Automatic];
  dimension = If[dimensionSpec === Automatic,
    FeynmanTrick`Private`DimensionExpression[], dimensionSpec];
  If[!validDimensionQ[dimension],
    Message[CreateFamily::dimension, dimension];
    Return[$Failed, Module]
  ];

  definition = familyDefinition[
    source, dimension, numericalPoint, eliminatedPositions];
  familyID = "ft-family-" <>
    stableDigest[{"FeynmanTrick.FamilyDefinition/v1", definition}];

  name = Which[
    nameOption =!= Automatic, nameOption,
    KeyExistsQ[source, "Name"] && source["Name"] =!= Automatic,
      source["Name"],
    True, "family_" <> StringTake[familyID, -12]
  ];
  If[!validFamilyNameQ[name],
    Message[CreateFamily::name, name];
    Return[$Failed, Module]
  ];

  sequenceSpec = resolveEmbedded[
    input, source, "CombinationSequence", sequenceOption];
  sequence = normalizeCombinationSequence[
    sequenceSpec, n, eliminatedPositions];
  If[sequence === $Failed, Return[$Failed, Module]];

  outputSpec = If[positionalTargets =!= Automatic,
    positionalTargets,
    resolveEmbedded[input, source, "OutputIntegrals", outputOption]
  ];
  outputSource = Which[
    positionalTargets =!= Automatic, "Positional",
    outputOption =!= Automatic, "Option",
    KeyExistsQ[input, "OutputIntegrals"] ||
      KeyExistsQ[source, "OutputIntegrals"], "Embedded",
    True, "Automatic"
  ];
  normalizedTargets = NormalizeOutputIntegrals[
    outputSpec, n, eliminatedPositions];
  If[normalizedTargets === $Failed, Return[$Failed, Module]];
  normalizedTargets = ValidateOutputIntegralsForSequence[
    normalizedTargets, sequence, eliminatedPositions];
  If[normalizedTargets === $Failed, Return[$Failed, Module]];

  topology = FeynmanTrick`FIREInterface`DefineTopology[
    name,
    source["LoopMomenta"],
    source["ExternalMomenta"],
    source["Propagators"],
    source["Replacements"]
  ];
  (* Existing topology metadata and exact family metadata survive ingestion;
     canonical structural fields and the validated arity remain authoritative. *)
  topology = Join[
    topology,
    KeyDrop[source, {"CombinationSequence", "OutputIntegrals"}],
    <|
      "Name" -> name,
      "NumPropagators" -> n,
      "EliminatedPositions" -> eliminatedPositions,
      "NumericalPoint" -> numericalPoint,
      "Dimension" -> dimension
    |>
  ];

  requests = If[normalizedTargets === All,
    {allRequestRecord[familyID]},
    MapIndexed[requestRecord[familyID, #1, First[#2]] &,
      normalizedTargets]
  ];

  <|
    "Schema" -> "FeynmanTrick.FamilySpec/v1",
    "FamilyID" -> familyID,
    "FamilyIdentity" -> familyID,
    "InputKind" -> inputKind,
    "Name" -> name,
    "Definition" -> definition,
    "Topology" -> topology,
    "TopTopology" -> topology,
    "NumPropagators" -> n,
    "EliminatedPositions" -> eliminatedPositions,
    "NumericalPoint" -> numericalPoint,
    "Dimension" -> dimension,
    "CombinationSequence" -> sequence,
    "OutputIntegralSource" -> outputSource,
    "OutputIntegralMode" ->
      If[normalizedTargets === All, "AllPendingDiscovery", "Explicit"],
    "OutputIntegrals" -> normalizedTargets,
    "L0OutputRequests" -> requests
  |>
];

CreateFamily[input_Association, opts : OptionsPattern[]] :=
  createFamily[
    input,
    Automatic,
    OptionValue["Name"],
    OptionValue["CombinationSequence"],
    OptionValue["OutputIntegrals"]
  ];

CreateFamily[input_Association, targets : (All | _List),
    opts : OptionsPattern[]] :=
  createFamily[
    input,
    targets,
    OptionValue["Name"],
    OptionValue["CombinationSequence"],
    OptionValue["OutputIntegrals"]
  ];

CreateFamily[___] := (Message[CreateFamily::input]; $Failed);

End[];
EndPackage[];
