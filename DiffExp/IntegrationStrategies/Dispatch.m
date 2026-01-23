(* IntegrationStrategies/Dispatch.m *)
(* Strategy dispatch logic *)

(* Dispatch to appropriate strategy based on configuration and problem type *)
DispatchStrategy[intind_, bVec_, line_, epsord_, BufferedData_] := Module[
  {cIndices, fGeneral, result, useRationalRecurrence, useSingularRecurrence},

  (* Check if rational recurrence should be used (non-singular points) *)
  useRationalRecurrence = (
    DiffExp`State`FEC[UseRationalRecurrence] === True &&
    RationalRecurrenceApplicableQ[intind, line] &&
    (* bVec must not contain Logx (non-singular point solutions are Logx-free) *)
    !DiffExp`Utilities`DependsQ[bVec, DiffExp`Symbols`Logx]
  );

  (* Check if singular recurrence should be used (simple pole, non-resonant) *)
  useSingularRecurrence = (
    !useRationalRecurrence &&
    DiffExp`State`FEC[UseRationalRecurrence] === True &&
    SingularRecurrenceApplicableQ[intind, line] &&
    (* Solutions at non-resonant singular points are Logx-free *)
    !DiffExp`Utilities`DependsQ[bVec, DiffExp`Symbols`Logx]
  );

  Which[
    (* Simple case: single integral without homogeneous components *)
    Length[intind] === 1 && DiffExp`Utilities`PChop[DiffExp`State`DEqnMatricesExpanded[line][0][[intind, intind]]] == {{0}},
    {cIndices, fGeneral} = SolveSimple[intind, bVec, line, epsord];
    {cIndices, fGeneral, BufferedData}

    (* Rational recurrence method for non-singular points with rational matrices *)
    , useRationalRecurrence,
    SolveRationalRecurrence[intind, bVec, line, epsord, BufferedData]

    (* Singular recurrence method for regular singular points *)
    , useSingularRecurrence,
    SolveSingularRecurrence[intind, bVec, line, epsord, BufferedData]

    (* Default strategy *)
    , DiffExp`State`FEC[IntegrationStrategy] === "Default",
    SolveDefault[intind, bVec, line, epsord, BufferedData]

    (* VOP strategy *)
    , DiffExp`State`FEC[IntegrationStrategy] === "VOP" || DiffExp`State`FEC[IntegrationStrategy] === "VariationOfParameters",
    SolveVOP[intind, bVec, line, epsord, BufferedData]

    (* VOPAlt strategy *)
    , DiffExp`State`FEC[IntegrationStrategy] === "VOPAlt",
    SolveVOPAlt[intind, bVec, line, epsord, BufferedData]
  ]
];
