(* ::Package:: *)
(* Exact per-level FIRE reductions shared by the DiffExp2 ladder consumers. *)

BeginPackage["FeynmanTrick`LevelReduction`", {"FeynmanTrick`"}];

BoundaryRequestRecords::usage =
  "BoundaryRequestRecords[masters, combinedPositions] classifies each lower-level master and returns the exact upper-level integral it needs.";
LevelIBPBatchSpec::usage =
  "LevelIBPBatchSpec[ftData, upperLevel] returns the content-addressed specification for one level's complete boundary-reduction batch.";
PrepareLevelIBPBatch::usage =
  "PrepareLevelIBPBatch[ftData, upperLevel] performs the level's exact FIRE reductions once and returns a validated reusable bundle.";
CollectLevelIBPSingularFactors::usage =
  "CollectLevelIBPSingularFactors[ftData, level, batch] extracts Feynman-parameter denominator factors from a prepared level bundle.";
RequiredTransportEpsilonOrder::usage =
  "RequiredTransportEpsilonOrder[ftData, level, epsilonOrder, prefactors, batch] computes the transport epsilon depth required by exact IBP coefficients.";

Begin["`Private`"];

$levelIBPBatchSchema = "FeynmanTrick.LevelIBPBatch/v1";

configuration[] := FeynmanTrick`FTConfiguration[];

dimensionExpression[] := Module[{cfg = configuration[], expr, eps},
  eps = Lookup[cfg, "EpsilonSymbol", FeynmanTrick`FTeps];
  expr = Lookup[cfg, "DimensionExpression", Automatic];
  If[expr === Automatic, 4 - 2 eps, expr]
];

verbosity[] := Module[{value = Lookup[configuration[], "Verbosity", 1]},
  If[IntegerQ[value], value, 1]
];

(* FIRE reduction coefficients are exact.  In particular, a tiny exact
   Laurent coefficient still fixes the epsilon pole depth and must never be
   chopped merely because of its magnitude. *)
zeroCoefficientQ[coefficient_] := TrueQ[PossibleZeroQ[coefficient]];

laurentZero[minPower_Integer, maxPower_Integer] := <|
  "MinPower" -> minPower,
  "Coefficients" -> Table[0, {Max[0, maxPower - minPower + 1]}]
|>;

laurentTrim[laurent_Association] := Module[
  {minPower = laurent["MinPower"], coefficients = laurent["Coefficients"]},
  While[Length[coefficients] > 0 && zeroCoefficientQ[First[coefficients]],
    coefficients = Rest[coefficients];
    minPower++
  ];
  If[coefficients === {},
    <|"MinPower" -> 0, "Coefficients" -> {0}|>,
    <|"MinPower" -> minPower, "Coefficients" -> coefficients|>
  ]
];

expandIBPCoefficientLaurent[coefficient_, maxPower_Integer] := Module[
  {cfg, dimensionVariable, epsilonSymbol, expanded, minPower, coefficients},
  cfg = configuration[];
  dimensionVariable = Lookup[cfg, "DimensionVariable", Global`d];
  epsilonSymbol = Lookup[cfg, "EpsilonSymbol", FeynmanTrick`FTeps];

  If[coefficient === 0, Return[laurentZero[0, maxPower], Module]];

  expanded = Quiet[Check[
    Series[coefficient /. dimensionVariable -> dimensionExpression[],
      {epsilonSymbol, 0, maxPower}],
    $Failed
  ]];
  If[expanded === $Failed, Return[$Failed, Module]];

  minPower = If[Head[expanded] === SeriesData,
    expanded[[4]]/expanded[[6]],
    0
  ];
  (* Rational functions of the dimension have integer Laurent powers.  A
     fractional/unknown epsilon frame is unsupported here; flooring it would
     silently under-budget downstream transport. *)
  If[!IntegerQ[minPower], Return[$Failed, Module]];

  coefficients = Quiet[Check[
    Table[
      Together[Normal[If[Head[expanded] === SeriesData,
        SeriesCoefficient[expanded, power],
        SeriesCoefficient[
          expanded + O[epsilonSymbol]^(maxPower + 1),
          {epsilonSymbol, 0, power}]
      ]]],
      {power, minPower, maxPower}
    ],
    $Failed
  ]];
  If[coefficients === $Failed, Return[$Failed, Module]];
  laurentTrim[<|
    "MinPower" -> minPower,
    "Coefficients" -> coefficients
  |>]
];

BoundaryRequestRecords[mastersAtLevel_List, combinedPositions_List] := Module[
  {positionI, positionJ},
  If[Length[combinedPositions] =!= 2, Return[$Failed, Module]];
  {positionI, positionJ} = combinedPositions;
  Table[
    Module[{masterVector, vi, vj, neededVector, case},
      masterVector = mastersAtLevel[[masterIndex]];
      If[!And @@ (1 <= # <= Length[masterVector] & /@
          {positionI, positionJ}), Return[$Failed, Module]];
      vi = masterVector[[positionI]];
      vj = masterVector[[positionJ]];
      case = Which[
        vi > 0 && vj > 0, "integrate",
        vi > 0 && vj == 0, "limitUpper",
        vi == 0 && vj > 0, "limitLower",
        True, "direct"
      ];
      neededVector = masterVector;
      Switch[case,
        "integrate",
          neededVector[[positionI]] = vi + vj;
          neededVector[[positionJ]] = 0,
        "limitUpper",
          neededVector[[positionI]] = vi;
          neededVector[[positionJ]] = 0,
        "limitLower",
          neededVector[[positionI]] = vj;
          neededVector[[positionJ]] = 0,
        "direct",
          neededVector[[positionJ]] = 0
      ];
      <|
        "MasterIndex" -> masterIndex,
        "MasterVec" -> masterVector,
        "Vi" -> vi,
        "Vj" -> vj,
        "Case" -> case,
        "NeededVec" -> neededVector
      |>
    ],
    {masterIndex, Length[mastersAtLevel]}
  ]
];

levelIBPBatchPayloadKey[reductions_Association,
    coefficientVectors_Association] := IntegerString[
  Hash[{reductions, coefficientVectors}, "SHA256"], 16, 64];

LevelIBPBatchSpec[ftData_Association, upperLevel_Integer] := Module[
  {levels, levelAbove, levelBelow, topology, mastersAbove, mastersBelow,
   combinedPositions, requests, neededIntegrals, normalizedNeededIntegrals,
   setupFingerprintRecord, keyRecord, key},
  levels = Lookup[ftData, "Levels", Missing["NotAvailable"]];
  If[!AssociationQ[levels] || upperLevel <= 0 ||
      !KeyExistsQ[levels, upperLevel] ||
      !KeyExistsQ[levels, upperLevel - 1],
    Return[$Failed, Module]];
  levelAbove = levels[upperLevel];
  levelBelow = levels[upperLevel - 1];
  topology = Lookup[levelAbove, "Topology", Missing["NotAvailable"]];
  mastersAbove = Lookup[levelAbove, "Masters", Missing["NotAvailable"]];
  mastersBelow = Lookup[levelBelow, "Masters", Missing["NotAvailable"]];
  combinedPositions = Lookup[
    levelAbove, "CombinedPositions", Missing["NotAvailable"]];
  If[!AssociationQ[topology] || !ListQ[mastersAbove] ||
      !ListQ[mastersBelow] || !MatchQ[combinedPositions, {_Integer, _Integer}],
    Return[$Failed, Module]];
  requests = BoundaryRequestRecords[mastersBelow, combinedPositions];
  If[requests === $Failed || MemberQ[requests, $Failed],
    Return[$Failed, Module]];
  neededIntegrals = DeleteDuplicates[#["NeededVec"] & /@ requests];
  normalizedNeededIntegrals =
    FeynmanTrick`FIREInterface`Private`normalizeIntegralIndices[
      topology, neededIntegrals];
  If[normalizedNeededIntegrals === $Failed ||
      MemberQ[normalizedNeededIntegrals, $Failed], Return[$Failed, Module]];
  normalizedNeededIntegrals = Sort[DeleteDuplicates[normalizedNeededIntegrals]];
  setupFingerprintRecord =
    FeynmanTrick`FIREInterface`Private`reductionSetupFingerprintRecord[
      topology];
  If[setupFingerprintRecord === $Failed, Return[$Failed, Module]];
  keyRecord = {
    $levelIBPBatchSchema, upperLevel,
    setupFingerprintRecord,
    Lookup[levelAbove, "FeynmanParameter", Missing["NotAvailable"]],
    combinedPositions, mastersBelow, mastersAbove, normalizedNeededIntegrals
  };
  key = IntegerString[Hash[keyRecord, "SHA256"], 16, 64];
  <|
    "Schema" -> $levelIBPBatchSchema,
    "Key" -> key,
    "KeyRecord" -> keyRecord,
    "UpperLevel" -> upperLevel,
    "Topology" -> topology,
    "MastersAbove" -> mastersAbove,
    "BoundaryRequests" -> requests,
    "NeededIntegrals" -> neededIntegrals
  |>
];

PrepareLevelIBPBatch[ftData_Association, upperLevel_Integer] := Module[
  {spec, reductions, coefficientVectors, mastersAbove, payloadKey},
  spec = LevelIBPBatchSpec[ftData, upperLevel];
  If[spec === $Failed, Return[$Failed, Module]];
  reductions = If[spec["NeededIntegrals"] === {}, <||>,
    FeynmanTrick`FIREInterface`ReduceIntegrals[
      spec["Topology"], spec["NeededIntegrals"]]];
  If[reductions === $Failed || !AssociationQ[reductions] ||
      !AllTrue[spec["NeededIntegrals"], KeyExistsQ[reductions, #] &],
    Return[$Failed, Module]];
  mastersAbove = spec["MastersAbove"];
  coefficientVectors = AssociationMap[
    Function[integral, Table[
      Coefficient[reductions[integral], Global`G[1, mastersAbove[[j]]]],
      {j, Length[mastersAbove]}]],
    spec["NeededIntegrals"]];
  payloadKey = levelIBPBatchPayloadKey[reductions, coefficientVectors];
  Join[spec, <|
    "Reductions" -> reductions,
    "CoefficientVectors" -> coefficientVectors,
    "PayloadKey" -> payloadKey
  |>]
];

validLevelIBPBatchQ[batch_, ftData_Association, upperLevel_Integer] := Module[
  {spec = LevelIBPBatchSpec[ftData, upperLevel], needed},
  If[spec === $Failed || !AssociationQ[batch], Return[False, Module]];
  needed = spec["NeededIntegrals"];
  TrueQ[
    Lookup[batch, "Schema", None] === $levelIBPBatchSchema &&
    Lookup[batch, "Key", None] === spec["Key"] &&
    Lookup[batch, "KeyRecord", None] === spec["KeyRecord"] &&
    Lookup[batch, "UpperLevel", None] === upperLevel &&
    Lookup[batch, "BoundaryRequests", None] === spec["BoundaryRequests"] &&
    Lookup[batch, "NeededIntegrals", None] === needed &&
    AssociationQ[Lookup[batch, "Reductions", None]] &&
    AssociationQ[Lookup[batch, "CoefficientVectors", None]] &&
    Lookup[batch, "PayloadKey", None] === levelIBPBatchPayloadKey[
      batch["Reductions"], batch["CoefficientVectors"]] &&
    AllTrue[needed,
      KeyExistsQ[batch["Reductions"], #] &&
      KeyExistsQ[batch["CoefficientVectors"], #] &&
      ListQ[batch["CoefficientVectors"][#]] &&
      Length[batch["CoefficientVectors"][#]] ===
        Length[spec["MastersAbove"]] &]
  ]
];

resolveLevelIBPBatch[ftData_Association, upperLevel_Integer, batch_] := Which[
  batch === Automatic,
    PrepareLevelIBPBatch[ftData, upperLevel],
  validLevelIBPBatchQ[batch, ftData, upperLevel],
    batch,
  True,
    If[verbosity[] >= 1,
      Print["Error: stale or mismatched level IBP batch for level ",
        upperLevel]];
    $Failed
];

CollectLevelIBPSingularFactors[ftData_Association, level_Integer,
    suppliedBatch_:Automatic] := Module[
  {levels, levelBelow, levelAbove, mastersBelow, feynmanParameter,
   variableName, factors = {}, batch, boundaryRequests, reductions},
  levels = Lookup[ftData, "Levels", Missing["NotAvailable"]];
  If[!AssociationQ[levels] || level <= 0 ||
      !KeyExistsQ[levels, level] || !KeyExistsQ[levels, level - 1],
    Return[{}, Module]];
  levelBelow = levels[level - 1];
  levelAbove = levels[level];
  mastersBelow = Lookup[levelBelow, "Masters", {}];
  feynmanParameter = Lookup[
    levelAbove, "FeynmanParameter", Missing["NotAvailable"]];
  If[Head[feynmanParameter] =!= Symbol, Return[{}, Module]];
  variableName = SymbolName[feynmanParameter];
  batch = resolveLevelIBPBatch[ftData, level, suppliedBatch];
  If[batch === $Failed,
    Return[If[suppliedBatch === Automatic, {}, $Failed], Module]];
  boundaryRequests = batch["BoundaryRequests"];
  reductions = batch["Reductions"];

  Do[
    Module[{request, neededVector, ibpCoefficients},
      request = boundaryRequests[[masterIndex]];
      neededVector = request["NeededVec"];
      If[!KeyExistsQ[reductions, neededVector], Continue[]];
      ibpCoefficients = batch["CoefficientVectors"][neededVector];
      factors = Join[factors, Flatten[Table[
        Module[{denominator, factorList},
          denominator = Denominator[Together[coefficient]];
          factorList = If[denominator === 1, {},
            First /@ FactorList[Factor[denominator]]];
          Select[factorList,
            !FreeQ[#, symbol_Symbol /;
              SymbolName[symbol] === variableName] &]
        ],
        {coefficient, ibpCoefficients}
      ]]]
    ],
    {masterIndex, Length[mastersBelow]}
  ];
  DeleteDuplicates[
    DeleteCases[Factor /@ factors, 0 | 1 | -1],
    TrueQ[PossibleZeroQ[Expand[#1 - #2]]] ||
      TrueQ[PossibleZeroQ[Expand[#1 + #2]]] &
  ]
];

RequiredTransportEpsilonOrder[ftData_Association, level_Integer,
    epsilonOrder_Integer, epsilonPrefactors_List:{},
    suppliedBatch_:Automatic] := Module[
  {levels, levelBelow, levelAbove, mastersBelow, mastersAbove,
   combinedPositions, prefactors, required = epsilonOrder, maxProbe,
   boundaryRequests, integrationPoleAllowance, batch, reductions,
   budgetStatus},
  levels = Lookup[ftData, "Levels", Missing["NotAvailable"]];
  If[!AssociationQ[levels] || level <= 0 ||
      !KeyExistsQ[levels, level] || !KeyExistsQ[levels, level - 1],
    Return[epsilonOrder + If[epsilonPrefactors === {}, 0,
      Max[epsilonPrefactors]], Module]];
  levelBelow = levels[level - 1];
  levelAbove = levels[level];
  mastersBelow = Lookup[levelBelow, "Masters", {}];
  mastersAbove = Lookup[levelAbove, "Masters", {}];
  combinedPositions = Lookup[levelAbove, "CombinedPositions", {}];
  prefactors = If[Length[epsilonPrefactors] === Length[mastersAbove],
    epsilonPrefactors,
    ConstantArray[0, Length[mastersAbove]]
  ];
  maxProbe = epsilonOrder + If[prefactors === {}, 0, Max[prefactors]] + 20;
  boundaryRequests = BoundaryRequestRecords[mastersBelow, combinedPositions];
  If[boundaryRequests === $Failed, Return[$Failed, Module]];
  integrationPoleAllowance = If[
    AnyTrue[boundaryRequests, #["Case"] === "integrate" &],
    Module[{raw, parsed, configured},
      raw = Environment["FT_INTEGRATION_POLE_ALLOWANCE"];
      parsed = If[StringQ[raw], Quiet[Check[ToExpression[raw], None]], None];
      configured = Lookup[configuration[], "IntegrationPoleAllowance", 4];
      Which[
        IntegerQ[parsed] && parsed >= 0, parsed,
        IntegerQ[configured] && configured >= 0, configured,
        True, 4
      ]
    ],
    0
  ];
  required = Max[required, epsilonOrder + integrationPoleAllowance];
  batch = resolveLevelIBPBatch[ftData, level, suppliedBatch];
  If[batch === $Failed,
    Return[If[suppliedBatch === Automatic,
      Max[0, Ceiling[required]], $Failed], Module]];
  boundaryRequests = batch["BoundaryRequests"];
  reductions = batch["Reductions"];

  budgetStatus = Catch[Do[
    Module[{request, neededVector, ibpCoefficients, coefficientLaurent},
      request = boundaryRequests[[masterIndex]];
      neededVector = request["NeededVec"];
      If[!KeyExistsQ[reductions, neededVector], Continue[]];
      ibpCoefficients = batch["CoefficientVectors"][neededVector];
      Do[
        If[ibpCoefficients[[j]] =!= 0,
          coefficientLaurent = expandIBPCoefficientLaurent[
            ibpCoefficients[[j]], maxProbe];
          If[coefficientLaurent === $Failed,
            Throw[$Failed, epsilonBudgetFailure]];
          required = Max[required,
            epsilonOrder + prefactors[[j]] -
              coefficientLaurent["MinPower"] + integrationPoleAllowance]
        ],
        {j, Length[mastersAbove]}
      ]
    ],
    {masterIndex, Length[mastersBelow]}
  ], epsilonBudgetFailure];
  If[budgetStatus === $Failed, Return[$Failed, Module]];
  Max[0, Ceiling[required]]
];

End[];
EndPackage[];
