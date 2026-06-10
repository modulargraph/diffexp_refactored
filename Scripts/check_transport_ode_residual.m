(* Check, per master and epsilon order, whether the saved transported local
   series of a chosen segment satisfy the level ODE
     d/dx J_n = sum_k m'(x) A_k(main(x)) J_{n-k}
   as truncated symbolic expansions around the segment center (Logx handled
   as Log[x], theta symbols checked on both branches).  A failing
   master/order localizes a wrong local solve. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
$Path = DeleteDuplicates[Prepend[$Path, repoRoot]];

Quiet[Get["FeynmanTrick/FeynmanTrick.m"], {General::shdw, Symbol::shdw}];
Quiet[Get["DiffExp.m"], {General::shdw, Symbol::shdw}];

envOrDefault[name_, default_] := Module[{value = Environment[name]},
  If[StringQ[value] && StringLength[StringTrim[value]] > 0, value, default]
];

saved = Get[envOrDefault["TRANSPORT_FILE", "/tmp/ft_transport_save/transport_level_1.m"]];
segIndex = ToExpression[envOrDefault["SEG_INDEX", "33"]];
checkXOrder = ToExpression[envOrDefault["CHECK_X_ORDER", "10"]];
maxEpsCheck = ToExpression[envOrDefault["MAX_EPS_CHECK", "6"]];

transport = saved["TransportResult"];
matrixDir = saved["MatrixDir"];
segs = transport["SegmentData"];
seg = segs[[segIndex]];
numMasters = transport["NumIntegrals"];
xLocal = DiffExp`Symbols`x;
LL = Unique["LogSym"];

Print["level=", saved["Level"], " segments=", Length[segs],
  " masters=", numMasters, " segIndex=", segIndex];

xMainExpr = DiffExp`Utilities`PChop[Expand[
  DiffExp`RegularizedIntegration`Private`segmentMainExpression[seg]]];
xMainExpr = FeynmanTrick`DiffExpIntegration`Private`snapMainExpression[
  xMainExpr, seg[[4]], {0, 1, 1/2}
];

matFull = Get[First[FileNames["*_full.m", matrixDir]]];
epsSym = FeynmanTrick`FTeps;
mainVar = SelectFirst[Variables[Level[matFull, {-1}]],
  StringMatchQ[SymbolName[#], "xx" ~~ DigitCharacter ...] &];

series = DiffExp`RegularizedIntegration`Private`uncompressSeriesData[seg[[5]]];
numEps = Min[Length[series[[1]]], maxEpsCheck + 1];

(* Truncated normal form of a saved series: polynomial in xLocal (possibly
   negative powers), LL (= Log[x]), theta symbols. *)
toNormal[ser_SeriesData] := Module[{nmin = ser[[4]], den = ser[[6]], coeffs = ser[[3]], top},
  top = Min[Length[coeffs], (checkXOrder + 2) * den - nmin + 1];
  Sum[
    coeffs[[i]] * xLocal^((nmin + i - 1)/den),
    {i, top}
  ] /. DiffExp`Symbols`Logx -> LL
];
toNormal[expr_] := expr /. DiffExp`Symbols`Logx -> LL;

(* Matrix in local coordinates as truncated Laurent polynomials in xLocal. *)
jacExpr = D[xMainExpr, xLocal];
epsLocal = Unique["ee"];
matLocal = jacExpr * (matFull /. mainVar -> xMainExpr) /. epsSym -> epsLocal;
matEpsOrders = Table[
  Map[
    Function[entry,
      Module[{c = SeriesCoefficient[entry + O[epsLocal]^(numEps + 1), {epsLocal, 0, k}]},
        If[PossibleZeroQ[c],
          0,
          Normal[Quiet[Series[c, {xLocal, 0, checkXOrder + 2}]]]
        ]
      ]
    ],
    matLocal, {2}
  ],
  {k, 0, numEps - 1}
];

normals = Table[toNormal[series[[m, n]]], {m, numMasters}, {n, numEps}];

(* d/dx with LL = Log[x]: D[f] + (coefficientwise dLL/dx = 1/x) *)
ddx[expr_] := D[expr, xLocal] + D[expr, LL]/xLocal;

thetaBranches = {
  {DiffExp`Symbols`\[Theta]p -> 1, DiffExp`Symbols`\[Theta]m -> 0},
  {DiffExp`Symbols`\[Theta]p -> 0, DiffExp`Symbols`\[Theta]m -> 1}
};

Do[
  Module[{maxRel = 0, badOrders = {}},
    Do[
      Module[{lhs, rhs, resid, scale, residTrunc, vals},
        lhs = ddx[normals[[m, n + 1]]];
        rhs = Sum[
          Sum[
            If[matEpsOrders[[k + 1, m, j]] === 0 || normals[[j, n - k + 1]] === 0,
              0,
              matEpsOrders[[k + 1, m, j]] * normals[[j, n - k + 1]]
            ],
            {j, numMasters}
          ],
          {k, 0, n}
        ];
        resid = Expand[lhs - rhs];
        (* keep x powers <= checkXOrder - 2 (both lhs/rhs truncations safe) *)
        residTrunc = Normal[Quiet[resid + O[xLocal]^(checkXOrder - 2)]];
        scale = Max[Flatten[{1, Abs[N[
          (normals[[m, n + 1]] /. LL -> 1 /. Thread[thetaBranches[[1]]]) /.
            xLocal -> 1/10, 30]]}]];
        vals = Table[
          Module[{r = residTrunc /. branch /. LL -> Log[xLocal]},
            Max[Flatten[{0, Abs[N[
              Table[r /. xLocal -> SetPrecision[pt, 80], {pt, {-1/10, -1/20}}],
            30]]}]]
          ],
          {branch, thetaBranches}
        ];
        Module[{rr = Max[vals]/scale},
          If[rr > maxRel, maxRel = rr];
          If[rr > 10^-12, AppendTo[badOrders, {n, ScientificForm[N[rr, 2]]}]];
        ];
      ],
      {n, 0, numEps - 1}
    ];
    Print["MASTER ", m, " maxRelResidual=", ScientificForm[N[maxRel, 3]],
      If[badOrders =!= {}, "  bad: " <> ToString[badOrders], ""]];
  ],
  {m, numMasters}
];

Quit[0];
