(* IntegrationStrategies/Dispatch.m *)
(* Strategy dispatch logic *)

(* Dispatch to appropriate strategy based on configuration and problem type *)
DispatchStrategy[ctx_Association, bVec_, epsord_, cacheIn_Association] := Module[
  {cIndices, fGeneral, result, cache = cacheIn,
   useRationalRecurrence, useSingularRecurrence, useGeneralSingularRecurrence,
   useFuchsianizedSingularRecurrence, singularEigenData},

  If[StringQ[Environment["DEBUG_DUMP_DISPATCH_DIR"]],
    Module[{ddir = Environment["DEBUG_DUMP_DISPATCH_DIR"], nfile},
      nfile = Length[FileNames["dispatch_*.m", ddir]] + 1;
      Put[<|"ctx" -> ctx, "bVec" -> bVec, "epsord" -> epsord|>,
        FileNameJoin[{ddir, "dispatch_" <> IntegerString[nfile, 10, 4] <> ".m"}]];
    ];
  ];

    (* Check if rational recurrence should be used (non-singular points) *)
    useRationalRecurrence = (
      ctx["UseRationalRecurrence"] === True &&
      RationalRecurrenceApplicableQ[ctx]
    );

  (* Check if singular recurrence should be used (simple pole, non-resonant) *)
  (* Uses PrepareSingularRecurrence which combines the check and eigenvalue computation *)
  useSingularRecurrence = False;
  If[!useRationalRecurrence && ctx["UseRationalRecurrence"] === True &&
     !DiffExp`Utilities`DependsQ[bVec, DiffExp`Symbols`Logx],
    (* If the solver cache already exists, we know this strategy was previously selected *)
    If[KeyExistsQ[cache, "SingRR"],
      useSingularRecurrence = True;
      ,
      singularEigenData = PrepareSingularRecurrence[ctx];
      If[singularEigenData =!= $Failed,
        useSingularRecurrence = True;
        (* Store eigenvalue data in cache so the solver doesn't recompute it *)
        If[!KeyExistsQ[cache, "SingularEigenData"],
          cache["SingularEigenData"] = singularEigenData;
        ];
      ];
    ];
  ];

  (* Check if general singular recurrence should be used (simple pole, any eigenvalues) *)
  (* This handles resonant eigenvalues, Jordan blocks, and Logx in bVec *)
  useGeneralSingularRecurrence = (
    !useRationalRecurrence &&
    !useSingularRecurrence &&
    ctx["UseRationalRecurrence"] === True &&
    GeneralSingularRecurrenceApplicableQ[ctx]
  );

  useFuchsianizedSingularRecurrence = (
    !useRationalRecurrence &&
    !useSingularRecurrence &&
    !useGeneralSingularRecurrence &&
    ctx["UseRationalRecurrence"] === True &&
    FuchsianizedSingularRecurrenceApplicableQ[ctx]
  );

  (* Optional dispatch trace for debugging local solves: enable by setting
     DiffExp`State`$LogStrategyDispatch = True (survives LoadConfiguration). *)
  If[TrueQ[DiffExp`State`$LogStrategyDispatch],
    If[!ListQ[DiffExp`State`StrategyDispatchLog],
      DiffExp`State`StrategyDispatchLog = {}
    ];
    AppendTo[
      DiffExp`State`StrategyDispatchLog,
      <|
        "Label" -> ctx["Label"],
        "EpsilonOrder" -> epsord,
        "SystemSize" -> ctx["SystemSize"],
        "Strategy" -> Which[
          ctx["SystemSize"] === 1 &&
            Normal[DiffExp`Utilities`PChop[ctx["AMatExpanded"]]] === {{0}},
          "Simple",
          useRationalRecurrence, "RationalRecurrence",
          useSingularRecurrence, "SingularRecurrence",
          useGeneralSingularRecurrence, "GeneralSingularRecurrence",
          useFuchsianizedSingularRecurrence, "FuchsianizedSingularRecurrence",
          True, "Other"
        ]
      |>
    ];
  ];

  Which[
    (* Simple case: single integral without homogeneous components *)
    ctx["SystemSize"] === 1 && Normal[DiffExp`Utilities`PChop[ctx["AMatExpanded"]]] === {{0}},
    {cIndices, fGeneral} = SolveSimple[ctx, bVec, epsord];
    {cIndices, fGeneral, cache}

    (* Rational recurrence method for non-singular points with rational matrices *)
    , useRationalRecurrence,
    SolveRationalRecurrence[ctx, bVec, epsord, cache]

    (* Singular recurrence method for non-resonant regular singular points *)
    , useSingularRecurrence,
    SolveSingularRecurrence[ctx, bVec, epsord, cache]

    (* General singular recurrence for resonant/non-diagonalizable regular singular points *)
    , useGeneralSingularRecurrence,
    SolveGeneralSingularRecurrence[ctx, bVec, epsord, cache]

    (* Higher-pole singular systems that can be locally fuchsianized *)
    , useFuchsianizedSingularRecurrence,
    SolveFuchsianizedSingularRecurrence[ctx, bVec, epsord, cache]

    (* Recurrence was requested, but no recursive strategy accepted this block. *)
    , ctx["UseRationalRecurrence"] === True,
    DiffExp`State`LastErrorContext = {ctx, bVec, epsord, cache};
    DiffExp`Utilities`ReportError[
      "UseRationalRecurrence -> True was requested, but no recursive finite-width strategy accepted integral(s) ",
      ctx["Label"],
      " at epsilon order ",
      epsord,
      ". Refusing to fall back to the default Wronskian/Frobenius path."
    ]

    (* Default strategy *)
    , ctx["IntegrationStrategy"] === "Default",
    SolveDefault[ctx, bVec, epsord, cache]

    (* VOP strategy *)
    , ctx["IntegrationStrategy"] === "VOP" || ctx["IntegrationStrategy"] === "VariationOfParameters",
    SolveVOP[ctx, bVec, epsord, cache]

    (* VOPAlt strategy *)
    , ctx["IntegrationStrategy"] === "VOPAlt",
    SolveVOPAlt[ctx, bVec, epsord, cache]
  ]
];
