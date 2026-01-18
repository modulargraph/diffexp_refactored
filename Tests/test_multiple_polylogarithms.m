(* Test: Multiple Polylogarithms evaluation using DiffExp *)
(* Tests evaluation of Goncharov polylogarithms G[a1,...,an,x] *)

(* Set directory and paths *)
SetDirectory[DirectoryName[$InputFileName]];
scriptDir = Directory[];
parentDir = ParentDirectory[scriptDir];

(* Add DiffExp to path *)
$Path = Prepend[$Path, FileNameJoin[{parentDir, "DiffExp"}]];

Print["==========================================="];
Print["Test: Multiple Polylogarithms Evaluation"];
Print["===========================================\n"];

(* Load the refactored DiffExp package *)
Print["Loading DiffExp package..."];
Get["DiffExp.m"];
Print["Package loaded!\n"];

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
  ExpansionOrder -> 20,  (* LOW for quick testing *)
  UseMobius -> True,
  UsePade -> False,
  DivisionOrder -> 3,
  Verbosity -> 1
}];

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
(* Tests                                                                    *)
(* ========================================================================= *)

testsPassed = 0;
testsTotal = 0;

(* Test 1: Simple G function - G[1,0,1,4] *)
Print["=== Test 1: G[1,0,1,4] ==="];
testsTotal++;
result1 = Quiet[G[1, 0, 1, 4] /. G -> GEvaluate];
(* Check that result is numerical and has expected form *)
If[NumericQ[result1 /. pm -> 0],
  Print["  [PASS] G[1,0,1,4] evaluated to numerical result"];
  Print["  Result: ", result1];
  testsPassed++;
  ,
  Print["  [FAIL] G[1,0,1,4] did not evaluate to numerical result"];
];

(* Test 2: G function with negative argument - G[1,-10,0,4] *)
Print["\n=== Test 2: G[1,-10,0,4] ==="];
testsTotal++;
result2 = Quiet[G[1, -10, 0, 4] /. G -> GEvaluate];
If[NumericQ[result2 /. pm -> 0],
  Print["  [PASS] G[1,-10,0,4] evaluated to numerical result"];
  Print["  Result: ", result2];
  testsPassed++;
  ,
  Print["  [FAIL] G[1,-10,0,4] did not evaluate to numerical result"];
];

(* Cleanup temporary directory *)
Print["\nCleaning up temporary directory..."];
DeleteDirectory[TmpDirectory, DeleteContents -> True];

(* Summary *)
Print["\n==========================================="];
Print["Results: ", testsPassed, " / ", testsTotal, " tests passed"];
If[testsPassed === testsTotal,
  Print["All tests PASSED!"];
  ,
  Print["Some tests FAILED!"];
];
Print["==========================================="];
