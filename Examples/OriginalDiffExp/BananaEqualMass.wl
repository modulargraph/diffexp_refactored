(* DiffExp 2 reproduction of the equal-mass part of DiffExp 1's
   Examples/Banana.nb.

   The exact matrices are fetched separately with
     Scripts/fetch_original_banana_data.sh
   from DiffExp revision 784c8229bf92369a03f011a48e161522c8c54bbd.

   DiffExp 2 does not yet expose the old notebook's generic asymptotic
   boundary constructor.  The t=-1 seed below was therefore evaluated once
   with that pinned DiffExp 1 notebook at working precision 1000.  Its
   displayed 80 digits are more than the 25-digit saved endpoint oracle can
   test. *)

repoRoot = ParentDirectory[
  ParentDirectory[DirectoryName[$InputFileName]]];
Get[FileNameJoin[{repoRoot, "DiffExp2.m"}]];

environmentValue[name_String] := Quiet[
  Check[Environment[name], $Failed]];
integerEnvironment[name_String, default_Integer] := Module[{value},
  value = environmentValue[name];
  If[StringQ[value] && StringMatchQ[value, DigitCharacter ..],
    FromDigits[value], default]
];

dataDirectory = Replace[
  environmentValue["ORIGINAL_BANANA_DATA"],
  Except[_String] :> FileNameJoin[{
    DirectoryName[$InputFileName], "Data", "Banana"}]
];
equalMassDirectory = FileNameJoin[{dataDirectory, "EqualMass"}];
requiredFiles = {"dt_0.m", "dt_1.m"};
missingFiles = Select[
  FileNameJoin[{equalMassDirectory, #}] & /@ requiredFiles,
  !FileExistsQ[#] &
];
If[missingFiles =!= {},
  Print["Missing original equal-mass Banana data: ", missingFiles];
  Print["Run Scripts/fetch_original_banana_data.sh or set ",
    "ORIGINAL_BANANA_DATA."];
  Exit[2]
];

workingPrecision = integerEnvironment[
  "ORIGINAL_BANANA_WORKING_PRECISION", 100];
expansionOrder = integerEnvironment[
  "ORIGINAL_BANANA_EXPANSION_ORDER", 50];
epsilonOrder = 4;

DiffExp2`LoadConfiguration[{
  "RecurrenceBackend" -> "Cpp",
  "WorkingPrecision" -> workingPrecision,
  "ExpansionOrder" -> expansionOrder,
  "EpsilonOrder" -> epsilonOrder,
  "DivisionOrder" -> 2,
  "UsePade" -> True,
  "Verbosity" -> 0
}];

(* Rows are the four equal-mass masters and columns are eps^0,...,eps^4.
   This is an external numerical boundary, not a rigorous ball oracle. *)
boundaryAtMinusOne = {
  {
    0,
    0.6098142615973314758171843477675682786476851268581215492083784369709609591742826667825`80,
    0.2246831620574879673425088336760305200821496864053600751025356545995343050926440178729`80,
    1.4892677286998821044617316488935897021650734458050394436656984112396821052863861395151`80,
    -0.1565473321967297789519358641742159477818933654598919829576794790581059025101518181084`80
  },
  {
    0,
    2.0315611364240859983912125025668199346884866545759549624453331661622642306612357858196`80,
    2.5015987338179334195798206613354159669923686891906532609670009656329301206317588584070`80,
    1.6155825238660982227573351635736612126312522743982385587051728903653984165141949092109`80,
    10.8388645335762554167769955752532358451281534617665196169009649283486984509159155851137`80
  },
  {
    0,
    8.2681045358689687315430153454799888687286184838449789674473628501495693113214418501594`80,
    18.5951889391375558214192352312249597137513127021191037704195486584712737119302679940955`80,
    -7.8357482657303532985002807519183681338088655455091335817317062707456695632681624875081`80,
    94.2784399247556444094551889948288307388305879159143090744908559987787836077066963154757`80
  },
  {
    1,
    0,
    2.4674011002723396547086227499690377838284248518101976566033373440550112056048013107504`80,
    -1.2020569031595942853997381615114499907649862923404988817922715553418382057863130901865`80,
    3.8557765200959298072757631689279106536350502662104646086206027675196796951660765158590`80
  }
};

(* Saved numerical table from the original notebook at t=20. *)
referenceAt20 = {
  {
    0,
    -0.43111218711011440753058887570625234775`25 -
      1.46531981358168928811420549968675918393`25 I,
    15.89653222034845789132057998270246251081`25 -
      7.72801500720397251672883169498739962242`25 I,
    59.84058102236571360802300785386791558357`25 +
      59.33639678015538383601290922158884583655`25 I,
    -78.14760301949654720648161753194522302586`25 +
      241.6134185032753295924503092430610577872`25 I
  },
  {
    0,
    -0.42124100906276135590317328422352941217`25 +
      4.92898116976737860164588734632119973156`25 I,
    -28.97194748338290349682545273302885619591`25 +
      7.21693641303080076639080813479168920472`25 I,
    -72.65714043196471770166032229848029993149`25 -
      93.52413764607496274354097722014159017438`25 I,
    146.42486282105256234814942282521463486864`25 -
      293.54086971255504083340992568984902777595`25 I
  },
  {
    0,
    14.24555780317066670763531703106371433047`25 +
      7.89418052687550709137844634758758100669`25 I,
    41.01245222029950705727032503507939152056`25 +
      80.81033195741338610644944738680159019006`25 I,
    -130.47105599852056541944638101010685847346`25 +
      281.8170328545998371239286773067131426978`25 I,
    -652.20254238193932474842015963911016301467`25 +
      272.28910592681830856443755861144977786314`25 I
  },
  {
    1,
    0,
    2.46740110027233965470862274996903778383`25,
    -1.20205690315959428539973816151144999076`25,
    3.85577652009592980727576316892791065364`25
  }
};

{t, x, eps} = {Global`t, Global`x, Global`eps};
matrix0 = Get[FileNameJoin[{equalMassDirectory, requiredFiles[[1]]}]];
matrix1 = Get[FileNameJoin[{equalMassDirectory, requiredFiles[[2]]}]];
matrix = matrix0 + eps matrix1;

plainValue[value_] := Transpose[
  Table[
    DiffExp2`EpsilonCoefficient[value, order],
    {order, 0, epsilonOrder}
  ]
];

transportLeg[boundary_, from_, to_, name_String] := Module[
  {line, pulled, system, result, wallSeconds},
  line = from + (to - from) x;
  pulled = Map[
    Together,
    (matrix /. t -> line) (to - from),
    {2}
  ];
  system = DiffExp2`LoadSystem[<|
    "Matrix" -> pulled,
    "Variable" -> x
  |>];
  wallSeconds = AbsoluteTiming[
    result = Catch[
      DiffExp2`TransportEndpoint[system, boundary, 0, 1],
      "DiffExp2Error"
    ];
  ][[1]];
  If[FailureQ[result],
    Print["ORIGINAL_BANANA FAIL leg=", name, " ", result];
    Exit[1]
  ];
  Print[
    "ORIGINAL_BANANA leg=", name,
    " seconds=", N[wallSeconds, 8],
    " segments=", result["SegmentCount"]
  ];
  {plainValue[result["Value"]], wallSeconds}
];

(* The upper-half-plane route avoids the real regular-singular points
   t=0,4,16 while retaining the notebook's physical continuation. *)
leg1 = transportLeg[
  boundaryAtMinusOne, -1, -1 + 5 I, "vertical-up"];
leg2 = transportLeg[
  leg1[[1]], -1 + 5 I, 20 + 5 I, "horizontal"];
leg3 = transportLeg[
  leg2[[1]], 20 + 5 I, 20, "vertical-down"];

valueAt20 = leg3[[1]];
transportSeconds = Total[{leg1[[2]], leg2[[2]], leg3[[2]]}];
errorsByOrder = Table[
  Max[Abs[
    valueAt20[[All, order + 1]] -
    referenceAt20[[All, order + 1]]
  ]],
  {order, 0, epsilonOrder}
];
maxError = Max[errorsByOrder];

(* The original 41.237399-second notebook call marched from t=-1 to t=32
   and produced the saved t=20 table on that path.  It is useful historical
   context, but not a strict endpoint-for-endpoint timing comparison. *)
historicalDiffExp1PathSeconds = 41.237399;

Print[
  "ORIGINAL_BANANA equalMass dimension=4 epsilonOrder=", epsilonOrder,
  " expansionOrder=", expansionOrder,
  " transportSeconds=", N[transportSeconds, 8],
  " originalDiffExp1PathTo32Seconds=", historicalDiffExp1PathSeconds
];
Print[
  "ORIGINAL_BANANA errorsByOrder=",
  InputForm[N[errorsByOrder, 12]],
  " maxError=", InputForm[N[maxError, 12]]
];

If[!TrueQ[maxError < 10^-10],
  Print["ORIGINAL_BANANA FAIL: t=20 reference error exceeds 1e-10"];
  Exit[1]
];

Print["ORIGINAL_BANANA PASS"];
Exit[0];
