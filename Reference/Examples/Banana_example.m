(* Simple test script for DiffExp package - Banana example *)

(* Set the directory to the script location *)
SetDirectory[DirectoryName[$InputFileName]];
scriptDir = Directory[];
Print["Working directory: ", scriptDir];

(* Load the DiffExp package *)
Print["Loading DiffExp package..."];
Get[FileNameJoin[{ParentDirectory[scriptDir], "DiffExp.m"}]];
Print["Package loaded successfully!"];

(* Simple configuration for equal mass banana graph *)
Print["\n=== Equal Mass Banana Configuration ==="];
EqualMassConfiguration = {
  DeltaPrescriptions -> {t - 16 + I * \[Delta]},
  MatrixDirectory -> scriptDir <> "/Banana_EqualMass_Matrices/",
  Verbosity -> 2,
  UseMobius -> True,
  UsePade -> True
};

Print["Loading configuration..."];
LoadConfiguration[EqualMassConfiguration];

Print["\nCurrent configuration:"];
Print[CurrentConfiguration[]];

(* Boundary conditions *)
Print["\n=== Setting up boundary conditions ==="];
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
} // PrepareBoundaryConditions[#, <|t -> -1/x|>] &;

Print["Boundary conditions prepared"];

(* Transport to t = -1 *)
Print["\n=== Transporting to t = -1 ==="];
Results1 = TransportTo[EqualMassBoundaryConditions, <|t -> -1|>];
Print["Transport to t = -1 complete"];

(* Transport from t = -1 to t = 10 *)
Print["\n=== Transporting from t = -1 to t = 10 ==="];
Results2 = TransportTo[Results1, <|t -> x|>, 10, True];
Print["Transport complete!"];

(* Display results at t = 5 *)
Print["\n=== Results at t = 5 ==="];
ResultsForEval = ToPiecewise[Results2];
Print["Master integral values at t = 5:"];
Table[ResultsForEval[[i, j]][5] // N[#, 10] &, {i, 1, 4}, {j, 1, 4}] // TableForm // Print;

(* ========================================================================= *)
(* Unequal Mass (3-mass) Banana Graph Example                              *)
(* ========================================================================= *)

Print["\n\n========================================================================"];
Print["=== UNEQUAL MASS BANANA CONFIGURATION ==="];
Print["========================================================================\n"];

(* Configuration for unequal mass banana graph *)
UnequalMassConfiguration = {
  ChopPrecision -> 500,
  DivisionOrder -> 4,
  ExpansionOrder -> 70,
  MatrixDirectory -> scriptDir <> "/Banana_Matrices/",
  RadiusOfConvergence -> 10, (* Without setting this option, the intermediate expansions blow up. *)
  UseMobius -> True,
  UsePade -> True,
  WorkingPrecision -> 1000,
  Verbosity -> 2
};

(* First we generate a set of boundary conditions using the equal mass differential equations *)
Print["=== Generating boundary conditions at t = 1/2 using equal mass system ==="];
LoadConfiguration[EqualMassConfiguration];
(* Generate results at high precision *)
UpdateConfiguration[{DivisionOrder -> 4, ExpansionOrder -> 70}];
Print["Transporting to asymptotic limit..."];
TmpAsymptotic = TransportTo[EqualMassBoundaryConditions, <|t -> -1/x|>];
Print["Transporting to t = 1/2..."];
BoundaryConditionsAtHalf = TransportTo[TmpAsymptotic, <|t -> 1/2|>];
Print["Boundary conditions at t = 1/2 generated successfully!\n"];

(* The boundary conditions for the unequal mass family *)
Print["=== Setting up unequal mass boundary conditions ==="];
UnequalMassBoundaryConditions = {
  <|psq -> 1/2, mm1 -> 1, mm2 -> 1, mm3 -> 1, mm4 -> 1|>,
  {
    BoundaryConditionsAtHalf["SeriesValues"][[1]],
    BoundaryConditionsAtHalf["SeriesValues"][[1]],
    BoundaryConditionsAtHalf["SeriesValues"][[1]],
    BoundaryConditionsAtHalf["SeriesValues"][[1]],
    BoundaryConditionsAtHalf["SeriesValues"][[1]],
    BoundaryConditionsAtHalf["SeriesValues"][[1]],
    BoundaryConditionsAtHalf["SeriesValues"][[2]],
    BoundaryConditionsAtHalf["SeriesValues"][[2]],
    BoundaryConditionsAtHalf["SeriesValues"][[2]],
    BoundaryConditionsAtHalf["SeriesValues"][[2]],
    BoundaryConditionsAtHalf["SeriesValues"][[3]],
    BoundaryConditionsAtHalf["SeriesValues"][[4]],
    BoundaryConditionsAtHalf["SeriesValues"][[4]],
    BoundaryConditionsAtHalf["SeriesValues"][[4]],
    BoundaryConditionsAtHalf["SeriesValues"][[4]]
  }
};
Print["Unequal mass boundary conditions set up successfully!\n"];

(* Load the unequal mass configuration *)
Print["=== Loading unequal mass configuration ==="];
LoadConfiguration[UnequalMassConfiguration];
UpdateConfiguration[Verbosity -> 1];
Print["Configuration loaded\n"];

(* Example: obtain results in the mass configuration mm1 = 2, mm2 = 3/2, mm3 = 4/3, mm4 = 1 *)
Print["=== Transporting to mass configuration: mm1 = 2, mm2 = 3/2, mm3 = 4/3, mm4 = 1 ==="];
Print["Starting from: psq = 1/2, all masses = 1"];
Print["Target: psq = 1/2, mm1 = 1 + x, mm2 = 1 + x/2, mm3 = 1 + x/3, mm4 = 1, with x = 1\n"];

UnequalMassResults1 = TransportTo[
  UnequalMassBoundaryConditions,
  <|psq -> 1/2, mm1 -> 1 + x, mm2 -> 1 + x/2, mm3 -> 1 + x/3, mm4 -> 1|>,
  1
];
Print["Transport to intermediate point complete!\n"];

(* Now transport to a different momentum point *)
Print["=== Transporting to psq = 3.7 (random point) ==="];
Print["From: psq = 1/2, mm1 = 2, mm2 = 3/2, mm3 = 4/3, mm4 = 1"];
Print["To: psq = 3.7, mm1 = 2, mm2 = 3/2, mm3 = 4/3, mm4 = 1\n"];

UnequalMassResults2 = TransportTo[
  UnequalMassResults1,
  <|psq -> x, mm1 -> 2, mm2 -> 3/2, mm3 -> 4/3, mm4 -> 1|>,
  3.7,
  True
];
Print["Transport to psq = 3.7 complete!\n"];

(* Display results *)
Print["=== Unequal mass results at psq = 3.7 ==="];
Print["Mass configuration: mm1 = 2, mm2 = 3/2, mm3 = 4/3, mm4 = 1"];
Print["Evaluating first 5 master integrals at epsilon orders 0 through 4:\n"];
UnequalMassForEval = ToPiecewise[UnequalMassResults2];
Table[
  UnequalMassForEval[[i, j]][3.7] // N[#, 15] &,
  {i, 1, 5}, {j, 1, 5}
] // TableForm // Print;

Print["\n========================================================================"];
Print["=== ALL TESTS COMPLETED SUCCESSFULLY! ==="];
Print["========================================================================"];
