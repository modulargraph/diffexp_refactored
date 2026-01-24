(* IntegrationStrategies/Dispatch.m *)
(* Strategy dispatch logic *)

(* Dispatch to appropriate strategy based on configuration and problem type *)
DispatchStrategy[ctx_Association, bVec_, epsord_, cache_Association] := Module[
  {cIndices, fGeneral, result,
   useRationalRecurrence, useSingularRecurrence, useGeneralSingularRecurrence},

  (* Check if rational recurrence should be used (non-singular points) *)
  useRationalRecurrence = (
    ctx["UseRationalRecurrence"] === True &&
    RationalRecurrenceApplicableQ[ctx] &&
    (* bVec must not contain Logx (non-singular point solutions are Logx-free) *)
    !DiffExp`Utilities`DependsQ[bVec, DiffExp`Symbols`Logx]
  );

  (* Check if singular recurrence should be used (simple pole, non-resonant) *)
  useSingularRecurrence = (
    !useRationalRecurrence &&
    ctx["UseRationalRecurrence"] === True &&
    SingularRecurrenceApplicableQ[ctx] &&
    (* Solutions at non-resonant singular points are Logx-free *)
    !DiffExp`Utilities`DependsQ[bVec, DiffExp`Symbols`Logx]
  );

  (* Check if general singular recurrence should be used (simple pole, any eigenvalues) *)
  (* This handles resonant eigenvalues, Jordan blocks, and Logx in bVec *)
  useGeneralSingularRecurrence = (
    !useRationalRecurrence &&
    !useSingularRecurrence &&
    ctx["UseRationalRecurrence"] === True &&
    GeneralSingularRecurrenceApplicableQ[ctx]
  );

  Which[
    (* Simple case: single integral without homogeneous components *)
    ctx["SystemSize"] === 1 && DiffExp`Utilities`PChop[ctx["AMatExpanded"]] == {{0}},
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
