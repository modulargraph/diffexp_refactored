(* ::Package:: *)
(* EpsPrefactors - Find and apply epsilon prefactors to remove poles *)
(* Strategy: find D = diag(eps^k1,...,eps^kn) such that D*A*D^{-1} *)
(* has no eps poles. Extensible to handle more complex pole types. *)

BeginPackage["FeynmanTrick`EpsPrefactors`", {"FeynmanTrick`"}];

FindEpsPrefactors::usage =
  "FindEpsPrefactors[matrix, epsSymbol] finds integers {k1,...,kn} such that \
the transformed matrix D*A*D^{-1} (D = Diag[eps^k1,...,eps^kn]) \
has no eps poles. Returns {k1,...,kn} with k1=0 normalization. \
epsSymbol defaults to FTConfiguration[\"EpsilonSymbol\"].";

ApplyEpsPrefactors::usage =
  "ApplyEpsPrefactors[matrix, prefactors, epsSymbol] applies the transformation \
D*A*D^{-1} where D = Diag[eps^k1,...,eps^kn]. Returns the transformed matrix.";

MatrixPoleOrders::usage =
  "MatrixPoleOrders[matrix, epsSymbol] returns a matrix of pole orders in eps \
for each entry. Positive values indicate poles.";

CheckEpsPoles::usage =
  "CheckEpsPoles[matrix, epsSymbol] returns True if the matrix has any eps poles.";

Begin["`Private`"];

(* ============================================================ *)
(* MatrixPoleOrders                                              *)
(* For each entry A_ij, find the negative of the minimum power   *)
(* of eps. A positive value means there is a pole.               *)
(* ============================================================ *)

MatrixPoleOrders[matrix_, epsSymbol_:Automatic] :=
Module[{eps, n},
  eps = If[epsSymbol === Automatic,
    FeynmanTrick`Private`$FTConfig["EpsilonSymbol"],
    epsSymbol
  ];
  n = Length[matrix];

  Table[
    If[matrix[[i, j]] === 0,
      0,
      -poleOrder[matrix[[i, j]], eps]
    ],
    {i, n}, {j, n}
  ]
];

(* Find the minimum power of eps in an expression *)
(* Returns negative number if there are poles *)
poleOrder[expr_, eps_] :=
Module[{expanded, terms, powers, minPow},
  If[expr === 0, Return[0]];

  (* Try to find the leading power via Series *)
  (* This handles rational functions in eps *)
  Module[{ser},
    ser = Quiet[Series[expr, {eps, 0, 0}]];
    If[Head[ser] === SeriesData,
      Return[ser[[4]]]; (* The minimum exponent *)
    ];
  ];

  (* Fallback: try Exponent with Min *)
  minPow = Quiet[Exponent[expr, eps, Min]];
  If[NumericQ[minPow], Return[minPow]];

  (* If all else fails, assume no pole *)
  0
];


(* ============================================================ *)
(* CheckEpsPoles                                                 *)
(* Returns True if any entry has a pole in eps                   *)
(* ============================================================ *)

CheckEpsPoles[matrix_, epsSymbol_:Automatic] :=
Module[{orders},
  orders = MatrixPoleOrders[matrix, epsSymbol];
  Max[Flatten[orders]] > 0
];


(* ============================================================ *)
(* FindEpsPrefactors                                             *)
(* Finds integer shifts k_i such that eps^{k_i-k_j} * A_ij      *)
(* has no poles for all (i,j).                                   *)
(*                                                               *)
(* Constraint: k_i - k_j >= p_ij for all nonzero A_ij            *)
(* where p_ij is the pole order of A_ij.                         *)
(* Normalization: k_1 = 0.                                       *)
(* Minimize: sum of |k_i| (to keep prefactors small).            *)
(* ============================================================ *)

FindEpsPrefactors[matrix_, epsSymbol_:Automatic] :=
Module[{n, poleOrders, constraints, vars, k, solution,
        objectiveVars, allConstraints, result},
  n = Length[matrix];
  poleOrders = MatrixPoleOrders[matrix, epsSymbol];

  If[Max[Flatten[poleOrders]] <= 0,
    (* No poles - return zeros *)
    Return[Table[0, {n}]];
  ];

  (* Set up integer linear program *)
  vars = Table[k[i], {i, n}];

  (* Constraints: k_i - k_j >= p_ij for all (i,j) with A_ij != 0 and p_ij > 0 *)
  constraints = Flatten[Table[
    If[matrix[[i, j]] =!= 0 && poleOrders[[i, j]] > 0,
      k[i] - k[j] >= poleOrders[[i, j]],
      Nothing
    ],
    {i, n}, {j, n}
  ]];

  (* Normalization: k_1 = 0 *)
  AppendTo[constraints, k[1] == 0];

  (* Integer constraints *)
  Do[AppendTo[constraints, k[i] \[Element] Integers], {i, n}];

  (* Objective: minimize sum of |k_i| *)
  (* Use auxiliary variables for absolute values *)
  Module[{absVars, absConstraints, objective},
    absVars = Table[Global`abs[i], {i, n}];
    absConstraints = Flatten[Table[
      {Global`abs[i] >= k[i], Global`abs[i] >= -k[i],
       Global`abs[i] \[Element] Integers},
      {i, n}
    ]];

    objective = Total[absVars];
    allConstraints = Join[constraints, absConstraints];

    solution = Quiet[
      Minimize[{objective, allConstraints}, Join[vars, absVars], Integers]
    ];

    If[Head[solution] === Minimize || solution[[1]] === Infinity,
      (* Fallback: try a simpler greedy approach *)
      result = greedySolve[poleOrders, matrix, n];
    ,
      result = vars /. solution[[2]];
    ];
  ];

  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1 && Max[Abs[result]] > 0,
    Print["Eps prefactors: ", result];
  ];

  result
];


(* Greedy fallback solver for the constraint system *)
(* Uses a shortest-path-like approach *)
greedySolve[poleOrders_, matrix_, n_] :=
Module[{k, changed, iterations},
  k = Table[0, {n}];

  (* Iteratively satisfy constraints *)
  changed = True;
  iterations = 0;
  While[changed && iterations < 100,
    changed = False;
    iterations++;
    Do[
      Do[
        If[matrix[[i, j]] =!= 0 && poleOrders[[i, j]] > 0,
          If[k[[i]] - k[[j]] < poleOrders[[i, j]],
            k[[i]] = k[[j]] + poleOrders[[i, j]];
            changed = True;
          ];
        ];
      , {j, n}];
    , {i, n}];
  ];

  (* Normalize: shift so k[[1]] = 0 *)
  k - k[[1]]
];


(* ============================================================ *)
(* ApplyEpsPrefactors                                            *)
(* Applies D * A * D^{-1} where D = diag(eps^k_i)               *)
(* Result: A_new_ij = eps^{k_i - k_j} * A_old_ij                *)
(* ============================================================ *)

ApplyEpsPrefactors[matrix_, prefactors_List, epsSymbol_:Automatic] :=
Module[{eps, n},
  eps = If[epsSymbol === Automatic,
    FeynmanTrick`Private`$FTConfig["EpsilonSymbol"],
    epsSymbol
  ];
  n = Length[matrix];

  Table[
    If[prefactors[[i]] - prefactors[[j]] == 0,
      matrix[[i, j]],
      eps^(prefactors[[i]] - prefactors[[j]]) * matrix[[i, j]] // Together
    ],
    {i, n}, {j, n}
  ]
];


End[];
EndPackage[];
