(* Inspect the singularity decomposition and residual towers of one segment
   of a Laurent dump: terms (a, b), subtracted local powers, and the
   epsilon/Logx tower at each subtracted power (branch-resolved). *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
$Path = DeleteDuplicates[Prepend[$Path, repoRoot]];

Quiet[Get["DiffExp.m"], {General::shdw, Symbol::shdw}];

envOrDefault[name_, default_] := Module[{value = Environment[name]},
  If[StringQ[value] && StringLength[StringTrim[value]] > 0, value, default]
];

Get[envOrDefault["DUMP_FILE", "/tmp/diffexp_banana_l0_dumps/laurent_integral_0009.m"]];
dump = DiffExp`RegularizedIntegration`Private`laurentIntegralDump;
segIndex = ToExpression[envOrDefault["SEG_INDEX", "1"]];
seg = dump["SavedData"]["SegmentData"][[segIndex]];
epsMin = dump["EpsMinPower"];
tol = DiffExp`State`FEC[DiffExp`State`RationalizationTolerance];
xLocal = DiffExp`Symbols`x;

bounds = DiffExp`RegularizedIntegration`Private`segmentActualBounds[seg];
overlap = {Max[Min[bounds], dump["Bounds"][[1]]], Min[Max[bounds], dump["Bounds"][[2]]]};
localA = DiffExp`RegularizedIntegration`Private`segmentLocalCoordinateForValue[seg, overlap[[1]]];
localB = DiffExp`RegularizedIntegration`Private`segmentLocalCoordinateForValue[seg, overlap[[2]]];

Print["segment=", segIndex, " mainBounds=", N[bounds, 10],
  " localA=", N[localA, 10], " localB=", N[localB, 10]];

ser = DiffExp`RegularizedIntegration`Private`uncompressSeriesData[seg[[5]]][[1]];
decomp = DiffExp`SingularityDecomposition`DecomposeSingularity[ser];
Print["decomposition terms: ", Length[decomp]];
Do[
  Module[{t = decomp[[ti]], a, b, g, branchDir, branchPoint, branchRules,
          branchB, subtracted},
    a = t["a"]; b = t["b"]; g = t["g"];
    (* replicate the integration-time branch resolution: the singular point
       local 0 lies at one end of [localA, localB] after the flip logic *)
    branchPoint = If[Abs[N[localA]] < 10^-12, localB, localA];
    branchDir = Sign[Re[N[branchPoint]]];
    branchRules = DiffExp`RegularizedIntegration`Private`thetaRulesAtPoint[
      branchPoint, branchDir
    ];
    branchB = DiffExp`Utilities`PChop[Expand[b /. branchRules]];
    Print["TERM ", ti, ": a=", a, " rawB=", InputForm[N[b, 20] /. branchRules],
      " branchB=", InputForm[N[branchB, 25]]];
    subtracted = Select[
      DeleteDuplicates[Flatten[Table[
        If[MatchQ[g[[idx]], _SeriesData],
          Table[(g[[idx, 4]] + k - 1)/g[[idx, 6]], {k, Length[g[[idx, 3]]]}],
          {0}
        ],
        {idx, Length[g]}
      ]]],
      TrueQ[N[a + # + 1] <= 10^-12] &
    ];
    Print["  subtracted local powers: ", subtracted];
    Do[
      Module[{tower},
        tower = Table[
          Module[{c = If[MatchQ[g[[idx]], _SeriesData],
              Module[{i2 = lp*g[[idx, 6]] - g[[idx, 4]] + 1},
                If[IntegerQ[i2] && 1 <= i2 <= Length[g[[idx, 3]]], g[[idx, 3, i2]], 0]
              ],
              If[lp === 0, g[[idx]], 0]
            ]},
            DiffExp`Utilities`PChop[Expand[c /. branchRules]]
          ],
          {idx, Length[g]}
        ];
        Print["  lp=", lp, " tower logs per offset:"];
        Do[
          Module[{c = tower[[idx]], lps},
            lps = Select[DiffExp`SeriesOps`LogxPowerRange[c],
              !DiffExp`RegularizedIntegration`Private`EffectiveZeroExprQ[
                DiffExp`SeriesOps`LogxCoeffNS[c, #], tol] &];
            Print["    off=", idx - 1, " eps^", epsMin + idx - 1, " logs=",
              Table[{k, InputForm[N[DiffExp`SeriesOps`LogxCoeffNS[c, k], 15]]}, {k, lps}]
            ];
          ],
          {idx, Min[Length[tower], 7]}
        ];
      ],
      {lp, subtracted}
    ];
  ],
  {ti, Length[decomp]}
];

Quit[0];
