(* Minimal test for TransportLevel with pre-computed Feynman trick matrices *)
(* Tests the DiffExp configuration and transport from the FeynmanTrick package *)

Print["=== TransportLevel Diagnostic Test ===\n"];

(* Setup *)
SetDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$Path = Prepend[$Path, Directory[]];

(* Load packages *)
Print["Loading FeynmanTrick..."];
Quiet[Get["FeynmanTrick/FeynmanTrick.m"], {General::shdw, Symbol::shdw}];
Print["Loading DiffExp..."];
Quiet[Get["DiffExp.m"], {General::shdw, Symbol::shdw}];
Print["Packages loaded.\n"];

(* Use pre-computed matrices from the pipeline test *)
matrixDir = "/private/var/folders/r9/t81xsghn0r198pktdqz_0v700000gn/T/FT_pipeline_test_32517/Level_3_Matrices";

If[!DirectoryQ[matrixDir],
  Print["Error: Matrix directory not found: ", matrixDir];
  Print["Run test_feynmantrick_pipeline.m first to generate matrices."];
  Exit[1];
];

Print["Matrix directory: ", matrixDir];
Print["Files: ", FileNames["*.m", matrixDir]];

(* Read the first matrix to check dimensions *)
dxx0 = Get[FileNameJoin[{matrixDir, "dxx_0.m"}]];
Print["Matrix size: ", Dimensions[dxx0]];
Print["Matrix variables: ", Variables[dxx0]];
Print[""];

(* Create dummy boundary conditions (7 masters, 3 eps orders) *)
numMasters = 7;
epsOrder = 2;
precision = 200;

(* Simple boundary: identity-like values at each eps order *)
boundaryValues = Table[
  SetPrecision[If[i == j, 1, 0], precision],
  {i, numMasters}, {j, epsOrder + 1}
];
Print["Boundary dimensions: ", Dimensions[boundaryValues]];
Print[""];

(* --- Diagnostic: Check what TransportLevel does step by step --- *)
Print["=== Step-by-step TransportLevel diagnosis ===\n"];

fixedVal = 11/23;

Print["Step 1: Load DiffExp configuration"];
diffExpConfig = {
  DiffExp`State`MatrixDirectory -> matrixDir,
  System`WorkingPrecision -> precision,
  DiffExp`State`ChopPrecision -> precision - 50,
  DiffExp`State`ExpansionOrder -> 30,
  DiffExp`State`EpsilonOrder -> epsOrder,
  DiffExp`State`UseMobius -> False,
  DiffExp`State`UsePade -> False,
  DiffExp`State`Verbosity -> 1,
  DiffExp`State`SegmentationStrategy -> "Predivision"
};

DiffExp`LoadConfiguration[diffExpConfig];

Print["\nStep 2: Check DiffExp state"];
Print["  FEC[Variables] = ", DiffExp`State`FEC[System`Variables]];
Print["  ExternalScalesVal = ", DiffExp`State`ExternalScalesVal];
Print["  NumIntegrals = ", DiffExp`State`NumIntegrals];

detectedVar = First[DiffExp`State`FEC[System`Variables]];
Print["  detectedVar = ", detectedVar];
Print["  Context[detectedVar] = ", Context[detectedVar]];
Print["  SymbolName[detectedVar] = ", SymbolName[detectedVar]];
Print["  detectedVar === xx: ", detectedVar === xx];
Print["  detectedVar === Global`xx: ", detectedVar === Global`xx];
Print[""];

Print["Step 3: Create boundary conditions"];
startPoint = Association[detectedVar -> SetPrecision[fixedVal, precision]];
Print["  startPoint = ", startPoint];
Print["  Keys[startPoint] = ", Keys[startPoint]];
Print["  Keys match ExternalScalesVal: ",
  DiffExp`Utilities`DependsQ[Keys[startPoint], detectedVar]];
Print[""];

Print["Step 4: Create endpoint"];
endpoint = Association[detectedVar -> 0];
Print["  endpoint = ", endpoint];
Print["  Keys[endpoint] = ", Keys[endpoint]];
Print["  DependsQ check: ",
  !Or @@ (!DiffExp`Utilities`DependsQ[Keys[endpoint], #] & /@ DiffExp`State`ExternalScalesVal)];
Print[""];

Print["Step 5: Call TransportTo"];
Print["  bcs format: {startPoint, boundaryValues}"];
Print["  startPoint: ", startPoint];
Print["  Dimensions[boundaryValues]: ", Dimensions[boundaryValues]];

(* Call TransportTo directly *)
result = DiffExp`Transport`TransportTo[
  {startPoint, boundaryValues},
  endpoint,
  1,
  True
];

Print["\nStep 6: Check result"];
If[AssociationQ[result],
  Print["  PASS: TransportTo returned Association"];
  Print["  Keys: ", Keys[result]];
  Print["  EndpointIsSingularity: ", result["EndpointIsSingularity"]];
  If[KeyExistsQ[result, "SegmentData"],
    Print["  SegmentData length: ", Length[result["SegmentData"]]];
  ];
,
  Print["  FAIL: TransportTo returned: ", Head[result]];
  If[Head[result] === Symbol, Print["  Value: ", result]];
];

Print["\n=== Test Complete ==="];
