(* Test singularity decomposition for equal mass banana *)
(* See SINGULARITY_DECOMPOSITION_PROBLEM.md for problem description *)

SetDirectory[DirectoryName[$InputFileName]];
scriptDir = Directory[];
parentDir = ParentDirectory[scriptDir];
$Path = Prepend[$Path, FileNameJoin[{parentDir, "DiffExp"}]];

Print["=== Singularity Decomposition Test ===\n"];
Print["Loading DiffExp..."];
Get["DiffExp.m"];

(* Configuration: singularity at t = 16 *)
LoadConfiguration[{
  DeltaPrescriptions -> {t - 16 + I * \[Delta]},
  MatrixDirectory -> FileNameJoin[{scriptDir, "Banana_EqualMass_Matrices"}] <> "/",
  Verbosity -> 1,
  UseMobius -> True,
  UsePade -> False,
  ExpansionOrder -> 25,
  EpsilonOrder -> 4
}];

(* Boundary conditions for equal mass banana *)
EqualMassBoundaryConditions = {
  "?",
  "?",
  \[Epsilon] (1 + 3 \[Epsilon]) (1 + 4 \[Epsilon]) * (
    -4 E^(3 EulerGamma \[Epsilon]) Gamma[\[Epsilon]]^3/t +
    6 E^(3 EulerGamma \[Epsilon]) (-1/t)^(1 + \[Epsilon]) \[Epsilon] Gamma[-\[Epsilon]]^2 Gamma[\[Epsilon]]^3/Gamma[-2 \[Epsilon]] +
    8 E^(3 EulerGamma \[Epsilon]) (-1/t)^(1 + 2 \[Epsilon]) \[Epsilon] Gamma[-\[Epsilon]]^3 Gamma[\[Epsilon]] Gamma[2 \[Epsilon]]/Gamma[-3 \[Epsilon]] +
    3 E^(3 EulerGamma \[Epsilon]) (-1/t)^(1 + 3 \[Epsilon]) \[Epsilon] Gamma[-\[Epsilon]]^4 Gamma[3 \[Epsilon]]/Gamma[-4 \[Epsilon]]
  ),
  E^(3 EulerGamma \[Epsilon]) \[Epsilon]^3 Gamma[\[Epsilon]]^3
};

(* ============================================================ *)
(* Step 1: Transport to t = 15 (near singularity at t = 16) *)
(* ============================================================ *)
Print["\n=== Step 1: Transport to t = 15 ===\n"];

PreparedBCs = PrepareBoundaryConditions[EqualMassBoundaryConditions, <|t -> -1/x|>];
Results1 = TransportTo[PreparedBCs, <|t -> -1|>];
Results15 = TransportTo[Results1, <|t -> x|>, 15];

Print["Transported to t = 15.\n"];

(* ============================================================ *)
(* Step 2: IntegrateSystem centered at singularity t = 16 *)
(* Line: t = 16 + x, so x = 0 is the singularity *)
(* ============================================================ *)
Print["=== Step 2: IntegrateSystem around singularity t = 16 ===\n"];
Print["Line: t = 16 + x"];
Print["  x = 0  -> t = 16 (singularity)"];
Print["  x = -1 -> t = 15 (where we have BCs)\n"];

singLine = <|t -> 16 + x|>;
intResult = IntegrateSystem[Results15, singLine];

Print["\nIntegrateSystem complete.\n"];

(* ============================================================ *)
(* Step 3: Examine the structure of the result *)
(* ============================================================ *)
Print["=== Step 3: Structure of IntegrateSystem Output ===\n"];

Do[
  Print["Integral ", i, ":"];
  Do[
    ser = intResult["SeriesValues"][[i, j]];
    If[MatchQ[ser, _SeriesData],
      minPow = ser[[4]] / ser[[6]];
      maxPow = (ser[[5]] - 1) / ser[[6]];
      step = 1 / ser[[6]];
      maxLogPow = DiffExp`SeriesOps`MaxLogxPower[ser];
      Print["  eps^", j-1, ": x^(", minPow, ") to x^(", maxPow, "), step=", step, ", MaxLogx=", maxLogPow];
    ,
      Print["  eps^", j-1, ": ", Head[ser], " = ", Short[ser, 2]];
    ];
  , {j, intResult["EpsilonOrder"] + 1}];
  Print[""];
, {i, intResult["NumIntegrals"]}];

(* ============================================================ *)
(* Step 4: Apply DecomposeSingularity *)
(* ============================================================ *)
Print["=== Step 4: Decomposition ===\n"];

(* Test on each integral *)
decompositions = DecomposeSingularityAll[intResult];

Do[
  Print["--- Integral ", i, " ---"];
  PrintDecomposition[decompositions[[i]]];
  Print[""];
, {i, Length[decompositions]}];

(* Detailed output for integral 1 *)
Print["\n=== Detailed Analysis: Integral 1 ===\n"];
decomp1 = decompositions[[1]];
If[Length[decomp1] > 0,
  Print["a = ", decomp1[[1]]["a"]];
  Print["b = ", decomp1[[1]]["b"]];
  Print["\nThis means the singular behavior is: x^(", decomp1[[1]]["a"],
        If[decomp1[[1]]["b"] != 0, " + " <> ToString[decomp1[[1]]["b"]] <> "*eps", ""], ")"];
  Print["\ng(x, eps) is finite (starts at x^0):"];
  Do[
    Print["  g at eps^", j-1, ": ", Short[decomp1[[1]]["g"][[j]], 4]];
  , {j, Min[3, Length[decomp1[[1]]["g"]]]}];
];

(* Verification: check that g starts at x^0 or higher *)
Print["\n=== Verification ===\n"];
Do[
  decomp = decompositions[[i]];
  If[Length[decomp] > 0,
    gLeadingPowers = Table[
      If[MatchQ[decomp[[1]]["g"][[j]], _SeriesData],
        decomp[[1]]["g"][[j]][[4]] / decomp[[1]]["g"][[j]][[6]],
        0  (* constant = x^0 *)
      ],
      {j, Length[decomp[[1]]["g"]]}
    ];
    minGPower = Min[gLeadingPowers];
    If[minGPower >= 0,
      Print["Integral ", i, ": PASS - g starts at x^", minGPower, " (non-negative)"];
    ,
      Print["Integral ", i, ": FAIL - g starts at x^", minGPower, " (negative!)"];
    ];
  ,
    Print["Integral ", i, ": Empty decomposition"];
  ];
, {i, Length[decompositions]}];

Print["\n=== Done ==="];
