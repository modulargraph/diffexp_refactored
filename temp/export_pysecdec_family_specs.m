(* Export exact FeynmanTrick family/master FP specs for pySecDec. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
$Path = DeleteDuplicates[Prepend[$Path, repoRoot]];

Quiet[Get["FeynmanTrick/FeynmanTrick.m"], {General::shdw, Symbol::shdw}];

FeynmanTrick`SetFTOption["Threads", 1];
FeynmanTrick`SetFTOption["FThreads", 1];
FeynmanTrick`SetFTOption["Verbosity", 0];
FeynmanTrick`SetFTOption["DimensionExpression", 2 - 2*FeynmanTrick`FTeps];
FeynmanTrick`SetFTOption["ReductionCache", False];
FeynmanTrick`SetFTOption["FIRETimeoutSeconds", 1800];

envOrDefault[name_, default_] := Module[{value = Environment[name]},
  If[StringQ[value] && StringLength[StringTrim[value]] > 0, value, default]
];

makeTopology[name_] := Switch[name,
  "bubble",
    FeynmanTrick`FIREInterface`DefineTopology[
      "pysecdec_bubble",
      {Global`l1},
      {Global`p},
      {1 - Global`l1^2, 1 - (Global`l1 - Global`p)^2},
      {Global`p^2 -> -1}
    ],
  "sunrise",
    FeynmanTrick`FIREInterface`DefineTopology[
      "pysecdec_sunrise",
      {Global`l1, Global`l2},
      {Global`p},
      {
        1 - Global`l1^2,
        1 - Global`l2^2,
        1 - (-Global`l1 - Global`l2 + Global`p)^2
      },
      {Global`p^2 -> -1}
    ],
  "banana",
    FeynmanTrick`FIREInterface`DefineTopology[
      "pysecdec_banana",
      {Global`l1, Global`l2, Global`l3},
      {Global`p},
      {
        1 - Global`l1^2,
        1 - Global`l2^2,
        1 - Global`l3^2,
        1 - (-Global`l1 - Global`l2 - Global`l3 + Global`p)^2
      },
      {Global`p^2 -> -1}
    ]
];

sequence[name_] := Switch[name,
  "bubble", {{1, 2}},
  "sunrise", {{1, 2}, {1, 3}},
  "banana", {{1, 2}, {1, 3}, {1, 4}}
];

pythonString[expr_] := Module[{s},
  s = ToString[InputForm[expr]];
  StringReplace[s, {
    "^" -> "**",
    " " -> "",
    "FTeps" -> "eps",
    "Eps" -> "eps",
    "FeynmanTrick`FTeps" -> "eps",
    "Global`" -> "",
    "Gamma[" -> "gamma(",
    "[" -> "(",
    "]" -> ")"
  }]
];

unitSimplexMap[n_Integer] := Module[{params, rescaled, jac},
  If[n <= 1,
    Return[{{1}, {}, 1}, Module]
  ];
  params = Table[ToExpression["Global`z" <> ToString[i - 1]], {i, n - 1}];
  rescaled = Table[
    Which[
      j == 1,
        Product[params[[i]], {i, 1, n - 1}],
      j <= n - 1,
        (1 - params[[j - 1]]) * Product[params[[i]], {i, j, n - 1}],
      j == n,
        1 - params[[n - 1]]
    ],
    {j, n}
  ];
  jac = FullSimplify[
    Abs[Det[D[rescaled[[1 ;; n - 1]], {params}]]],
    And @@ Thread[0 < params < 1]
  ];
  {rescaled, params, jac}
];

specForMaster[example_String, level_Integer, master_List, ftData_Association,
    fixedValue_] := Module[
  {levelData, topo, levelParam, fixedRules, active, powers, activeProps,
   uf, U, F, fvars, n, L, v, simplex, rescaled, params, jac, subRules,
   Usub, Fsub, remainder, gammaArg, uPow, fPow, prefactor, name},

  levelData = ftData["Levels"][level];
  topo = levelData["Topology"];
  levelParam = Lookup[levelData, "FeynmanParameter", None];
  fixedRules = If[levelParam === None || MissingQ[levelParam],
    {},
    {levelParam -> fixedValue}
  ];

  active = Flatten@Position[master, _?(# > 0 &), {1}, Heads -> False];
  powers = master[[active]];
  activeProps = topo["Propagators"][[active]] /. fixedRules;
  n = Length[activeProps];
  L = Length[topo["LoopMomenta"]];
  v = Total[powers];

  uf = FeynmanTrick`BoundaryConditions`ComputeSymanzikPolynomials[
    activeProps, topo["LoopMomenta"], topo["Replacements"]
  ];
  If[uf === $Failed, Return[$Failed, Module]];
  {U, F, fvars} = uf;

  simplex = unitSimplexMap[n];
  {rescaled, params, jac} = simplex;
  subRules = Thread[fvars -> rescaled];
  Usub = FullSimplify[Together[U /. subRules], And @@ Thread[0 < params < 1]];
  Fsub = FullSimplify[Together[F /. subRules], And @@ Thread[0 < params < 1]];
  remainder = FullSimplify[
    jac * Product[rescaled[[i]]^(powers[[i]] - 1), {i, n}],
    And @@ Thread[0 < params < 1]
  ];

  gammaArg = v - L + L*FeynmanTrick`FTeps;
  uPow = v - (L + 1) + (L + 1)*FeynmanTrick`FTeps;
  fPow = -v + L - L*FeynmanTrick`FTeps;
  prefactor = Gamma[gammaArg] / Times @@ (Gamma /@ powers);

  name = StringJoin[
    example, "_L", ToString[level], "_",
    StringRiffle[ToString /@ master, "x"]
  ];

  <|
    "Name" -> name,
    "Example" -> example,
    "Level" -> level,
    "Master" -> master,
    "ActivePositions" -> active,
    "ActivePowers" -> powers,
    "LoopCount" -> L,
    "Variables" -> (SymbolName /@ params),
    "U" -> pythonString[Usub],
    "F" -> pythonString[Fsub],
    "UPower" -> pythonString[uPow],
    "FPower" -> pythonString[fPow],
    "Remainder" -> pythonString[remainder],
    "Prefactor" -> pythonString[prefactor],
    "RequestedOrders" -> {0},
    "FixedValue" -> pythonString[fixedValue]
  |>
];

example = envOrDefault["FT_EXAMPLE", "bubble"];
outFile = envOrDefault[
  "PYSECDEC_SPEC_FILE",
  FileNameJoin[{$TemporaryDirectory, "pysecdec_" <> example <> "_specs.json"}]
];
includeNeeded = envOrDefault["INCLUDE_NEEDED", "1"] =!= "0";
fixedValue = Rationalize[ToExpression[envOrDefault["FT_FIXED_VALUE", "11/23"]], 0];

topology = makeTopology[example];
ftData = FeynmanTrick`FeynmanTrickIteration`DefineFTIteration[
  topology, sequence[example], {}
];
outDir = FileNameJoin[{$TemporaryDirectory,
  "pysecdec_specs_" <> example <> "_" <> ToString[$ProcessID]}];
If[DirectoryQ[outDir], DeleteDirectory[outDir, DeleteContents -> True]];
CreateDirectory[outDir, CreateIntermediateDirectories -> True];
ftData = FeynmanTrick`FeynmanTrickIteration`RunFullIteration[ftData, outDir];
If[ftData === $Failed, Exit[10]];

specs = {};
Do[
  masters = ftData["Levels"][level]["Masters"];
  Do[
    spec = specForMaster[example, level, master, ftData, fixedValue];
    If[AssociationQ[spec], AppendTo[specs, spec]],
    {master, masters}
  ];
  ,
  {level, 0, ftData["NumLevels"]}
];

If[includeNeeded,
  Do[
    If[level < ftData["NumLevels"],
      requests = FeynmanTrick`DiffExpIntegration`Private`BoundaryRequestRecords[
        ftData["Levels"][level]["Masters"],
        ftData["Levels"][level + 1]["CombinedPositions"]
      ];
      needed = DeleteDuplicates[#["NeededVec"] & /@ requests];
      Do[
        spec = specForMaster[example, level + 1, master, ftData, fixedValue];
        If[AssociationQ[spec],
          spec = Association[spec, "Name" -> spec["Name"] <> "_needed"];
          AppendTo[specs, spec]
        ],
        {master, needed}
      ];
    ],
    {level, 0, ftData["NumLevels"] - 1}
  ];
];

Export[outFile, specs, "RawJSON"];
Print["SpecFile=", outFile];
Print["SpecCount=", Length[specs]];
Do[
  Print[
    spec["Name"], " level=", spec["Level"],
    " master=", InputForm[spec["Master"]],
    " vars=", spec["Variables"],
    " pref=", spec["Prefactor"],
    " U=", spec["U"],
    " F=", spec["F"]
  ],
  {spec, specs}
];
