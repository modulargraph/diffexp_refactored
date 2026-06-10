(* Print the raw combined series of one dump segment at selected epsilon
   orders: leading x-coefficients with their Logx content. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
$Path = DeleteDuplicates[Prepend[$Path, repoRoot]];

Quiet[Get["DiffExp.m"], {General::shdw, Symbol::shdw}];

envOrDefault[name_, default_] := Module[{value = Environment[name]},
  If[StringQ[value] && StringLength[StringTrim[value]] > 0, value, default]
];

Get[envOrDefault["DUMP_FILE", "/tmp/diffexp_banana_l0_dumps/laurent_integral_0009.m"]];
dump = DiffExp`RegularizedIntegration`Private`laurentIntegralDump;
segIndex = ToExpression[envOrDefault["SEG_INDEX", "33"]];
maxOrders = ToExpression[envOrDefault["MAX_ORDERS", "4"]];
maxPowers = ToExpression[envOrDefault["MAX_POWERS", "5"]];
seg = dump["SavedData"]["SegmentData"][[segIndex]];
epsMin = dump["EpsMinPower"];

ser = DiffExp`RegularizedIntegration`Private`uncompressSeriesData[seg[[5]]][[1]];
Print["segment=", segIndex, " epsMin=", epsMin, " orders=", Length[ser]];

Do[
  Module[{s = ser[[k]]},
    Print["ORDER eps^", epsMin + k - 1, " head=", Head[s]];
    If[MatchQ[s, _SeriesData],
      Module[{nmin = s[[4]], den = s[[6]], coeffs = s[[3]]},
        Do[
          Module[{c = coeffs[[idx]], lps},
            If[!TrueQ[PossibleZeroQ[c]],
              lps = DiffExp`SeriesOps`LogxPowerRange[c];
              Print["  x^", (nmin + idx - 1)/den, ": ",
                Table[{lp, InputForm[Chop[N[
                  DiffExp`SeriesOps`LogxCoeffNS[c, lp] /. {
                    DiffExp`Symbols`\[Theta]p -> "TP",
                    DiffExp`Symbols`\[Theta]m -> "TM"
                  } /. {"TP" -> 1, "TM" -> 0},
                  12], 10^-25]]}, {lp, lps}]
              ];
            ];
          ],
          {idx, Min[maxPowers, Length[coeffs]]}
        ];
      ],
      Print["  value=", InputForm[N[s, 12]]];
    ];
  ],
  {k, Min[maxOrders, Length[ser]]}
];

Quit[0];
