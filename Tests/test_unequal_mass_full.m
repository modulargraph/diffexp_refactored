(* Full unequal mass banana test - tests actual transport computation *)
(* This mirrors the unequal mass section of Examples/Banana_example.m *)

(* Set directory and paths *)
SetDirectory[DirectoryName[$InputFileName]];
scriptDir = Directory[];
parentDir = ParentDirectory[scriptDir];

(* Add subpackages to path *)
$Path = Prepend[$Path, FileNameJoin[{parentDir, "DiffExp"}]];

Print["==========================================="];
Print["Full Unequal Mass Banana Test"];
Print["===========================================\n"];

(* Load the refactored DiffExp package *)
Print["Loading refactored DiffExp package..."];
Get["DiffExp.m"];
Print["Package loaded!\n"];

(* ====== Step 1: Generate boundary conditions from equal mass system ====== *)
Print["=== Step 1: Generate boundary conditions at t = 1/2 using equal mass ==="];
EqualMassConfiguration = {
  DeltaPrescriptions -> {t - 16 + I * \[Delta]},
  MatrixDirectory -> FileNameJoin[{scriptDir, "Banana_EqualMass_Matrices"}] <> "/",
  Verbosity -> 1,
  UseMobius -> True,
  UsePade -> True
};

Print["Loading equal mass configuration..."];
LoadConfiguration[EqualMassConfiguration];

(* Equal mass boundary conditions *)
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

Print["Preparing boundary conditions..."];
PreparedBCs = PrepareBoundaryConditions[EqualMassBoundaryConditions, <|t -> -1/x|>];

(* Generate results at high precision for unequal mass *)
Print["Updating configuration for high precision..."];
UpdateConfiguration[{DivisionOrder -> 4, ExpansionOrder -> 70}];

Print["Transporting to asymptotic limit..."];
TmpAsymptotic = TransportTo[PreparedBCs, <|t -> -1/x|>];

Print["Transporting to t = 1/2..."];
BoundaryConditionsAtHalf = TransportTo[TmpAsymptotic, <|t -> 1/2|>];
Print["Equal mass boundary conditions at t = 1/2 generated.\n"];

(* ====== Step 2: Load unequal mass configuration ====== *)
Print["=== Step 2: Load unequal mass configuration ==="];
UnequalMassConfiguration = {
  ChopPrecision -> 500,
  DivisionOrder -> 4,
  ExpansionOrder -> 70,
  MatrixDirectory -> FileNameJoin[{scriptDir, "Banana_Matrices"}] <> "/",
  RadiusOfConvergence -> 10,
  UseMobius -> True,
  UsePade -> True,
  WorkingPrecision -> 1000,
  Verbosity -> 1
};

Print["Loading unequal mass configuration..."];
LoadConfiguration[UnequalMassConfiguration];
Print["Unequal mass matrices loaded: ", DiffExp`State`NumIntegrals, " integrals\n"];

(* ====== Step 3: Prepare unequal mass boundary conditions ====== *)
Print["=== Step 3: Prepare unequal mass boundary conditions ==="];
(* The mapping for the unequal mass family:
   6 copies of integral 1, 4 copies of integral 2, 1 copy of integral 3, 4 copies of integral 4 *)
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
Print["Unequal mass boundary conditions prepared.\n"];

(* ====== Step 4: Transport in unequal mass system ====== *)
Print["=== Step 4: Transport in unequal mass system ==="];
Print["Transporting to mass configuration: mm1 = 2, mm2 = 3/2, mm3 = 4/3, mm4 = 1"];
Print["Starting from: psq = 1/2, all masses = 1"];
Print["Target: psq = 1/2, mm1 = 1 + x, mm2 = 1 + x/2, mm3 = 1 + x/3, mm4 = 1, with x = 1\n"];

UnequalMassResults1 = TransportTo[
  UnequalMassBoundaryConditions,
  <|psq -> 1/2, mm1 -> 1 + x, mm2 -> 1 + x/2, mm3 -> 1 + x/3, mm4 -> 1|>,
  1
];

Print["Transport complete!"];
Print["Result point: ", UnequalMassResults1["KinematicPoint"]];

Print["\n==========================================="];
Print["Test completed successfully!"];
Print["==========================================="];
