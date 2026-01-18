(* Multiple Polylogarithms evaluation using DiffExp *)
(* Demonstrates evaluation of Goncharov polylogarithms G[a1,...,an,x] *)

(* Set the directory to the script location *)
SetDirectory[DirectoryName[$InputFileName]];
scriptDir = Directory[];
Print["Working directory: ", scriptDir];

(* Load the DiffExp package *)
Print["Loading DiffExp package..."];
Get[FileNameJoin[{ParentDirectory[scriptDir], "DiffExp.m"}]];
Print["Package loaded successfully!"];

(* ========================================================================= *)
(* Functions for shuffling out basepoint divergences                        *)
(* ========================================================================= *)

ListShuffleProduct[u:{a_, x___}, v:{b_, y___}, c___] := Join[
  ListShuffleProduct[{x}, v, c, a],
  ListShuffleProduct[u, {y}, c, b]
];
ListShuffleProduct[{x___}, {y___}, c___] := {{c, x, y}};

ShuffleG[G[a__, x_], G[b__, x_]] := Total[
  G[Sequence @@ #, x] & /@ ListShuffleProduct[{a}, {b}]
];

ClearAll[GFB];
GFB[a___, a1_, b___, 0, x_] /; (!SameQ[a1, 0] && SameQ[b, 0]) := (
  G[a, a1, b, 0, x] + G[a, a1, x] G[b, 0, x] -
  ShuffleG[G[a, a1, x], G[b, 0, x]]
) /. G -> GFB;

GFB[a___, a1_, x_] /; !SameQ[a1, 0] := G[a, a1, x];

GFB[a___, x_] /; SameQ[a, 0] := Module[{l = Length@{a}},
  1/l! Log[x]^l
];

ExtractLogarithms[ex_] := ex /. G -> GFB;

(* ========================================================================= *)
(* DiffExp configuration                                                    *)
(* ========================================================================= *)

(* Create a temporary directory for the matrices *)
If[!ValueQ[TmpDirectory], TmpDirectory = CreateDirectory[]];
Print["Temporary directory: ", TmpDirectory];

(* Configuration for MPL evaluation - using LOW expansion order for quick tests *)
DiffExpConfiguration = Association[{
  EpsilonOrder -> 0,
  WorkingPrecision -> 250,
  ChopPrecision -> 225,
  ExpansionOrder -> 20,  (* LOW for quick testing - increase for more precision *)
  UseMobius -> True,
  UsePade -> False,
  DivisionOrder -> 3,
  Verbosity -> 1
}];

Print["Configuration: ", DiffExpConfiguration];

(* ========================================================================= *)
(* Evaluation functions                                                     *)
(* ========================================================================= *)

GEvaluateFiniteBasepoint[a__, endpoint_] /; Im[endpoint] === 0 := Module[
  {MPLMatrix, l = Length[{a}], MyConfiguration, BoundaryData, Res},

  (* Write out differential matrix *)
  MPLMatrix = Append[
    PadLeft[#, l + 1, 0] & /@ DiagonalMatrix[
      Table[1/(t - ind), {ind, {a}}]
    ],
    ConstantArray[0, l + 1]
  ];
  Export[TmpDirectory <> $PathnameSeparator <> "dt_0.m", MPLMatrix];

  (* Update the configuration *)
  MyConfiguration = Join[
    DiffExpConfiguration // Normal,
    {
      DeltaPrescriptions -> Table[
        t - ind + If[ind > 0, -I \[Delta], +I \[Delta]],
        {ind, Select[{a}, Im[#] === 0 &]}
      ],
      MatrixDirectory -> TmpDirectory
    }
  ];
  LoadConfiguration[MyConfiguration];

  (* Boundary data at t = 0 *)
  BoundaryData = Append[
    ConstantArray[0 + O[x]^(1/2), l],
    1 + O[x]^(1/2)
  ] // PrepareBoundaryConditions[#, {t -> endpoint x}] &;

  (* Obtain the numerical results *)
  Res = TransportTo[BoundaryData, {t -> endpoint}];
  (Res[[2, 1, 1]] + pm Res[[3, 1, 1]]) // N[#,
    Round[-Log10[Res[[3, 1, 1]]] + Log10[Abs@Res[[2, 1, 1]]]]
  ] &
];

GEvaluate[a__, endpoint_] /; Im[endpoint] === 0 :=
  Expand[
    (GFB[a, endpoint] /.
      G[a1__, b_] /; SameQ[a1] :>
        1/Length[{a1}]! Log[1 - b/{a1}[[1]]]^Length[{a1}]
    ) /. G[args__] :> GEvaluateFiniteBasepoint[args]
  ];

(* ========================================================================= *)
(* Examples                                                                 *)
(* ========================================================================= *)

Print["\n========================================================================"];
Print["=== MULTIPLE POLYLOGARITHM EVALUATION EXAMPLES ==="];
Print["========================================================================\n"];

Print["Using low expansion order (20) for quick testing."];
Print["Increase ExpansionOrder for higher precision results.\n"];

(* Example 1: Simple G function *)
Print["=== Example 1: G[1,0,1,4] ==="];
result1 = G[1, 0, 1, 4] /. G -> GEvaluate // AbsoluteTiming;
Print["Time: ", result1[[1]], " seconds"];
Print["Result: ", result1[[2]]];
Print["(The 'pm' indicates the error estimate)\n"];

(* Example 2: G function with negative argument *)
Print["=== Example 2: G[1,-10,0,4] ==="];
result2 = G[1, -10, 0, 4] /. G -> GEvaluate // AbsoluteTiming;
Print["Time: ", result2[[1]], " seconds"];
Print["Result: ", result2[[2]]];
Print[];

(* Example 3: G function with complex argument *)
Print["=== Example 3: G[10,-10+I,-1/2,-50,1] ==="];
result3 = G[10, -10 + I, -1/2, -50, 1] /. G -> GEvaluate // AbsoluteTiming;
Print["Time: ", result3[[1]], " seconds"];
Print["Result: ", result3[[2]]];
Print[];

Print["========================================================================"];
Print["=== ALL EXAMPLES COMPLETED SUCCESSFULLY! ==="];
Print["========================================================================"];
