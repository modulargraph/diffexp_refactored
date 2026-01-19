(* ::Package:: *)

(* DiffExp Singularity Decomposition Subpackage *)
(*
   Decomposes IntegrateSystem output near singularities into canonical form:

   f(x, eps) = sum_i x^(a_i + b_i*eps) * g_i(x, eps)

   where:
   - a_i in Q : singular power (can be negative)
   - b_i in Q : epsilon-dependent exponent (determined from Logx structure)
   - g_i(x, eps) : finite series in x, starting at x^0 or higher (no negative powers)

   See Tests/SINGULARITY_DECOMPOSITION_PROBLEM.md for full specification.
*)

BeginPackage["DiffExp`SingularityDecomposition`", {
  "DiffExp`Symbols`",
  "DiffExp`State`",
  "DiffExp`Utilities`",
  "DiffExp`SeriesOps`"
}];

(* Main function *)
DecomposeSingularity::usage = "DecomposeSingularity[seriesList] decomposes a list of series (one per eps order) from IntegrateSystem into canonical singularity form. Returns a list of terms {<|\"a\"->..., \"b\"->..., \"g\"->...|>, ...}. Uses RationalizationTolerance from configuration for determining a and b.";

DecomposeSingularityAll::usage = "DecomposeSingularityAll[integrateSystemOutput] applies DecomposeSingularity to all integrals.";

PrintDecomposition::usage = "PrintDecomposition[decomp] prints the decomposition in readable form.";

Begin["`Private`"];

(* ============================================================ *)
(* Helper functions *)
(* ============================================================ *)

(* Get the leading (minimum) power of x in a SeriesData *)
GetLeadingPower[ser_SeriesData] := ser[[4]] / ser[[6]];
GetLeadingPower[0] := Infinity;
GetLeadingPower[n_?NumericQ] := 0;

(* Get the coefficient at a specific power of x in a SeriesData *)
(* Power k means x^k; returns the full coefficient (may contain Logx) *)
GetCoefficientAtPower[ser_SeriesData, k_] := Module[
  {nmin, den, coeffs, idx},
  nmin = ser[[4]];
  den = ser[[6]];
  coeffs = ser[[3]];
  (* Index for power k: k = (nmin + idx - 1) / den, so idx = k*den - nmin + 1 *)
  idx = k * den - nmin + 1;
  If[idx < 1 || idx > Length[coeffs], 0, coeffs[[idx]]]
];
GetCoefficientAtPower[0, k_] := 0;
GetCoefficientAtPower[n_?NumericQ, k_] := If[k == 0, n, 0];

(* Multiply a series by x^(-a) to shift powers *)
(* After this, leading power becomes 0 *)
ShiftSeriesByPower[ser_SeriesData, a_] := Module[
  {nmin, nmax, den, newNmin, newNmax},
  nmin = ser[[4]];
  nmax = ser[[5]];
  den = ser[[6]];
  newNmin = nmin - a * den;
  newNmax = nmax - a * den;
  ReplacePart[ser, {4 -> newNmin, 5 -> newNmax}]
];
ShiftSeriesByPower[0, a_] := 0;
ShiftSeriesByPower[n_?NumericQ, a_] := n * DiffExp`Symbols`x^(-a);

(* Multiply series by x^(-b*eps) expanded to given epsilon order *)
(* x^(-b*eps) = exp(-b*eps*Logx) = Sum[(-b*Logx)^n/n! * eps^n, {n, 0, Infinity}] *)
(* This modifies the epsilon expansion coefficients *)
MultiplyByXMinusBEps[seriesList_List, b_] := Module[
  {epsOrder, result, n, k, factor},
  epsOrder = Length[seriesList] - 1;

  (* For each eps order n, we need to account for contributions from
     lower orders k via (-b*Logx)^(n-k)/(n-k)! *)
  result = Table[
    Sum[
      (* When n = k, the factor is 1 (avoid 0^0 issue) *)
      factor = If[n == k, 1, (-b * DiffExp`Symbols`Logx)^(n - k) / (n - k)!];
      factor * seriesList[[k + 1]],
      {k, 0, n}
    ] // DiffExp`SeriesOps`SExpand,
    {n, 0, epsOrder}
  ];
  result
];

(* Check if all series have non-negative leading powers *)
AllNonNegativePowers[seriesList_List] := Module[{minPower},
  minPower = Min[GetLeadingPower /@ seriesList];
  minPower >= 0
];

(* Check if all series are effectively zero (within numerical precision) *)
(* Uses RationalizationTolerance since this is for determining decomposition structure *)
AllEffectivelyZero[seriesList_List] := Module[{},
  And @@ (EffectivelyZero /@ seriesList)
];

EffectivelyZero[0] := True;
EffectivelyZero[n_?NumericQ] := Abs[n] < DiffExp`State`FEC[RationalizationTolerance];
EffectivelyZero[ser_SeriesData] := Module[{coeffs, maxAbs},
  coeffs = Flatten[{ser[[3]]}];
  If[Length[coeffs] == 0, Return[True]];
  (* Check if all coefficients (including Logx parts) are small *)
  maxAbs = Max[Abs[coeffs /. DiffExp`Symbols`Logx -> 0]];
  maxAbs < DiffExp`State`FEC[RationalizationTolerance]
];

(* ============================================================ *)
(* Main decomposition algorithm *)
(* ============================================================ *)

DecomposeSingularity[seriesList_List] := Module[
  {terms, current, a, b, g, c0, c1Logx, maxIter, iter, ratTol},

  (* Initialize *)
  terms = {};
  current = seriesList;
  maxIter = 20; (* Safety limit *)
  iter = 0;

  While[!AllEffectivelyZero[current] && iter < maxIter,
    iter++;

    (* Step 1: Find leading power a (most negative across all eps orders) *)
    a = Min[GetLeadingPower /@ current];

    (* If already non-negative, we're extracting a finite term *)
    If[a >= 0, a = 0];

    (* Handle case where current is effectively constant *)
    If[a === Infinity,
      Break[];
    ];

    (* Step 2: Determine b from Logx structure *)
    (* We need to find the first non-zero epsilon order to get a reference coefficient,
       then look at the next order to determine b from the Logx structure.
       x^(b*eps) = 1 + b*eps*Logx + ... so:
       - If eps^n has coeff c_n at x^a (without Logx)
       - Then eps^(n+1) has Logx coeff = b * c_n at x^a
       - So b = (Logx coeff at eps^(n+1)) / (non-Logx coeff at eps^n) *)
    b = 0;
    ratTol = DiffExp`State`FEC[RationalizationTolerance];
    Do[
      c0 = GetCoefficientAtPower[current[[refOrder]], a] /. DiffExp`Symbols`Logx -> 0;
      (* Use rationalization tolerance to check if c0 is effectively non-zero *)
      If[Abs[c0] > ratTol && Length[current] > refOrder,
        c1Logx = Coefficient[GetCoefficientAtPower[current[[refOrder + 1]], a], DiffExp`Symbols`Logx];
        b = c1Logx / c0;
        (* Clean up numerical noise: take real part if imaginary part is tiny *)
        If[NumericQ[b] && Abs[Im[b]] < ratTol,
          b = Re[b];
        ];
        (* Rationalize b - it should be an integer or simple fraction *)
        If[NumericQ[b],
          b = Rationalize[b, ratTol];
        ];
        Break[];
      ];
    , {refOrder, 1, Length[current] - 1}];

    (* Step 3: Extract g(x, eps) by factoring out x^(a + b*eps) *)
    (* g = x^(-a) * x^(-b*eps) * f *)
    (* Note: g coefficients keep their full precision from the original series *)
    g = ShiftSeriesByPower[#, a] & /@ current;
    g = MultiplyByXMinusBEps[g, b];

    (* Step 4: Compute residual = current - x^(a + b*eps) * g *)
    (* But we want g to be finite, so we need to truncate it to start at x^0 *)
    (* The finite part of g is the part with non-negative powers *)
    (* The residual is what's left: current - x^(a + b*eps) * g_finite *)

    (* For now, we assume that extracting one term per iteration works.
       The g we computed should be the full contribution from this singularity. *)

    (* Store the term *)
    AppendTo[terms, <|"a" -> a, "b" -> b, "g" -> g|>];

    (* If g is finite (starts at x^0), we're done with this singular structure *)
    (* Check if there's remaining singular structure by looking at g's leading power *)
    If[AllNonNegativePowers[g],
      (* g is fully finite; compute residual to check for more singular terms *)
      (* For a clean decomposition, we extract the leading coefficient and subtract *)
      (* For now, break after first term to avoid infinite loops *)
      Break[];
    ];

    (* If g still has negative powers, continue extracting *)
    current = g;
  ];

  If[iter >= maxIter,
    Print["Warning: DecomposeSingularity reached maximum iterations"];
  ];

  terms
];

DecomposeSingularityAll[integrateSystemOutput_List] :=
  Table[DecomposeSingularity[integrateSystemOutput[[i]]], {i, Length[integrateSystemOutput]}];

(* ============================================================ *)
(* Output formatting *)
(* ============================================================ *)

PrintDecomposition[terms_List] := Module[{},
  If[Length[terms] == 0,
    Print["Empty decomposition"];
    Return[];
  ];
  Print["Decomposition has ", Length[terms], " term(s):"];
  Do[
    Print["  Term ", i, ": x^(", terms[[i]]["a"],
          If[terms[[i]]["b"] != 0,
             " + " <> ToString[terms[[i]]["b"]] <> "*eps",
             ""],
          ") * g_", i, "(x, eps)"];
    (* Print first few coefficients of g at eps^0 *)
    If[MatchQ[terms[[i]]["g"][[1]], _SeriesData],
      Print["    g_", i, " at eps^0: ", Short[terms[[i]]["g"][[1]], 3]];
    ];
  , {i, Length[terms]}];
];

End[];

EndPackage[];
