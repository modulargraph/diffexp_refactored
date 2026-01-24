(* Test script for FeynmanTrick PropagatorAlgebra *)
(* Tests the core algebraic functions without requiring FIRE *)

SetDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$Path = Prepend[$Path, Directory[]];

Print["=== Testing FeynmanTrick PropagatorAlgebra ===\n"];

(* Load the package *)
Print["Loading FeynmanTrick package..."];
Get["FeynmanTrick.m"];
Print["Package loaded.\n"];

passed = 0;
failed = 0;

test[name_, expr_, expected_] := Module[{result},
  result = (expr === expected);
  If[result,
    Print["PASS: ", name];
    passed++;
  ,
    Print["FAIL: ", name];
    Print["  Expected: ", expected];
    Print["  Got:      ", expr];
    failed++;
  ];
];

testApprox[name_, expr_, expected_] := Module[{result},
  result = (Simplify[expr - expected] === 0);
  If[result,
    Print["PASS: ", name];
    passed++;
  ,
    Print["FAIL: ", name];
    Print["  Expected: ", expected];
    Print["  Got:      ", expr];
    Print["  Diff:     ", Simplify[expr - expected]];
    failed++;
  ];
];

(* ============================================================ *)
(* Test 1: CombinePropagatorSymbolic *)
(* ============================================================ *)

Print["\n--- Test: CombinePropagatorSymbolic ---"];

Module[{Di, Dj, result},
  Di = -k^2 + m1^2;
  Dj = -k^2 + 2*k*p + m2^2;
  result = FeynmanTrick`PropagatorAlgebra`CombinePropagatorSymbolic[Di, Dj, xx];
  testApprox["Combine two propagators",
    result,
    -k^2 + 2*(1-xx)*k*p + xx*m1^2 + (1-xx)*m2^2
  ];
];

(* ============================================================ *)
(* Test 2: FeynmanTrickDecomposition *)
(* ============================================================ *)

Print["\n--- Test: FeynmanTrickDecomposition ---"];

Module[{coeffMat, constVec, n = 5, cp = 2, op = 3},
  {coeffMat, constVec} = FeynmanTrick`PropagatorAlgebra`FeynmanTrickDecomposition[n, cp, op, xx];

  test["CoeffMatrix shape", Dimensions[coeffMat], {5, 5}];
  test["ConstVector shape", Length[constVec], 5];
  test["CoeffMatrix[cp,cp] = 1/xx", coeffMat[[cp, cp]], 1/xx];
  test["CoeffMatrix[cp,op] = -1/xx", coeffMat[[cp, op]], -1/xx];
  test["CoeffMatrix[1,1] = 0", coeffMat[[1, 1]], 0];
  test["ConstVector all zero", constVec, {0, 0, 0, 0, 0}];
];

(* ============================================================ *)
(* Test 3: DifferentiatedIntegrals with FeynmanTrick decomp *)
(* ============================================================ *)

Print["\n--- Test: DifferentiatedIntegrals (Feynman trick) ---"];

Module[{coeffMat, constVec, result, masterIdx},
  (* 4 propagators, combined at position 1, other at position 2 *)
  {coeffMat, constVec} = FeynmanTrick`PropagatorAlgebra`FeynmanTrickDecomposition[4, 1, 2, xx];

  (* Master with index {1,1,1,0} *)
  masterIdx = {1, 1, 1, 0};
  result = FeynmanTrick`PropagatorAlgebra`DifferentiatedIntegrals[masterIdx, coeffMat, constVec];

  Print["  Master: ", masterIdx];
  Print["  Shifted integrals: ", result];

  (* Expected: from k=1 (v_1=1):
     j=1: coeff = -1*(1/xx), shift k1+1, j1-1 → {1,1,1,0} (no change!) with coeff -1/xx
     j=2: coeff = -1*(-1/xx) = 1/xx, shift k1+1, j2-1 → {2,0,1,0} with coeff 1/xx
     from k=2 (v_2=1):
     coeffMatrix row 2 is all zeros, so no contribution
     from k=3 (v_3=1):
     coeffMatrix row 3 is all zeros, so no contribution *)
  test["Number of shifted integrals", Length[result], 2];

  (* Check the {1,1,1,0} term (diagonal, coeff -1/xx) *)
  Module[{diagTerm, offDiagTerm},
    diagTerm = Select[result, #[[1]] === {1, 1, 1, 0} &];
    offDiagTerm = Select[result, #[[1]] === {2, 0, 1, 0} &];

    If[Length[diagTerm] > 0,
      testApprox["Diagonal term coefficient", diagTerm[[1, 2]], -1/xx];
    ,
      Print["FAIL: Diagonal term not found"];
      failed++;
    ];

    If[Length[offDiagTerm] > 0,
      testApprox["Off-diagonal term coefficient", offDiagTerm[[1, 2]], 1/xx];
    ,
      Print["FAIL: Off-diagonal term not found"];
      failed++;
    ];
  ];
];

(* ============================================================ *)
(* Test 4: DecomposePropagatorDerivative - simple 1-loop case *)
(* ============================================================ *)

Print["\n--- Test: DecomposePropagatorDerivative (1-loop) ---"];

Module[{props, loops, result, coeffMat, constVec},
  (* Simple 1-loop with 2 propagators:
     D1 = -k^2 (massless tadpole)
     D2 = -(k-p)^2 = -k^2 + 2*k*p - p^2
     Combined: D_comb = xx*D1 + (1-xx)*D2 = -k^2 + 2*(1-xx)*k*p - (1-xx)*p^2
  *)
  props = {
    -k^2 + 2*(1-xx)*k*p - (1-xx)*p^2,  (* combined propagator *)
    -k^2 + 2*k*p - p^2                  (* D2 stays *)
  };
  loops = {k};

  {coeffMat, constVec} = FeynmanTrick`PropagatorAlgebra`DecomposePropagatorDerivative[
    props, loops, xx
  ];

  Print["  CoeffMatrix: ", coeffMat];
  Print["  ConstVector: ", constVec];

  (* Expected: d(D_comb)/dxx = -2*k*p + p^2 = D1_orig - D2
     D1_orig = -k^2, D2 = -k^2 + 2*k*p - p^2
     D1_orig - D2 = -2*k*p + p^2

     We need to express this in terms of current props:
     D_comb/xx - D2/xx: works if D_comb = xx*D1 + (1-xx)*D2
       D_comb/xx - D2/xx = D1 + (1-xx)/xx*D2 - D2/xx = D1 - D2/xx + (1-xx)/xx*D2
       = D1 + D2*((1-xx)-1)/xx = D1 - D2/xx ... hmm

     Actually: derivative of prop 1 w.r.t. xx:
     d/dxx (-k^2 + 2*(1-xx)*k*p - (1-xx)*p^2) = -2*k*p + p^2

     This should decompose as c1*prop1 + c2*prop2 + const where:
     loop-dependent part of derivative: -2*k*p (coefficient of k*p is -2, coeff of k^2 is 0)
     loop-dependent part of prop1: k^2 coeff = -1, k*p coeff = 2*(1-xx)
     loop-dependent part of prop2: k^2 coeff = -1, k*p coeff = 2

     System: c1*(-1) + c2*(-1) = 0  (from k^2)
             c1*2*(1-xx) + c2*2 = -2  (from k*p)
     From first eq: c1 = -c2
     Substituting: -c2*2*(1-xx) + c2*2 = -2 → c2*(2-2+2*xx) = -2 → c2*2*xx = -2 → c2 = -1/xx
     So c1 = 1/xx

     This matches FeynmanTrickDecomposition!
  *)

  testApprox["CoeffMat[1,1] = 1/xx", coeffMat[[1, 1]], 1/xx];
  testApprox["CoeffMat[1,2] = -1/xx", coeffMat[[1, 2]], -1/xx];
  test["CoeffMat[2,*] all zero", coeffMat[[2]], {0, 0}];

  (* Verify: the constant part should also match *)
  (* constVec[1] = derivConst - c1*propConst1 - c2*propConst2
     derivConst = p^2
     propConst1 = -(1-xx)*p^2
     propConst2 = -p^2
     constVec[1] = p^2 - (1/xx)*(-(1-xx)*p^2) - (-1/xx)*(-p^2)
                 = p^2 + (1-xx)/xx*p^2 - p^2/xx
                 = p^2 * (1 + (1-xx)/xx - 1/xx)
                 = p^2 * (1 + (1-xx-1)/xx)
                 = p^2 * (1 - 1/xx + 1/xx - 1/xx)
     Hmm let me redo: = p^2 + (1-xx)*p^2/xx - p^2/xx = p^2(1 + (1-xx-1)/xx) = p^2*(1-1/xx+1/xx-1/xx)
     Actually: = p^2 + p^2*(1-xx)/xx - p^2/xx = p^2*(1 + (1-xx-1)/xx) = p^2*(1 + (-xx)/xx) = p^2*(1-1) = 0
     Wait: (1-xx)/xx - 1/xx = (1-xx-1)/xx = -xx/xx = -1
     So: p^2 + p^2*(-1) = 0. Great! constVec[1] = 0.
  *)
  testApprox["ConstVec[1] = 0", Simplify[constVec[[1]]], 0];
  test["ConstVec[2] = 0", constVec[[2]], 0];
];

(* ============================================================ *)
(* Test 5: Full round-trip verification *)
(* ============================================================ *)

Print["\n--- Test: Full round-trip (derivative = Sum c_j D_j + const) ---"];

Module[{props, loops, coeffMat, constVec, reconstructed, deriv},
  props = {
    -k^2 + 2*(1-xx)*k*p - (1-xx)*p^2,
    -k^2 + 2*k*p - p^2
  };
  loops = {k};

  {coeffMat, constVec} = FeynmanTrick`PropagatorAlgebra`DecomposePropagatorDerivative[
    props, loops, xx
  ];

  (* Verify: Sum_j coeffMat[[1,j]]*props[[j]] + constVec[[1]] == D[props[[1]], xx] *)
  reconstructed = Expand[Sum[coeffMat[[1, j]] * props[[j]], {j, 2}] + constVec[[1]]];
  deriv = Expand[D[props[[1]], xx]];

  testApprox["Reconstruction matches derivative",
    Simplify[reconstructed - deriv], 0
  ];
];

(* ============================================================ *)
(* Test 6: 2-loop example *)
(* ============================================================ *)

Print["\n--- Test: 2-loop propagator decomposition ---"];

Module[{props, loops, coeffMat, constVec, deriv1, recon1},
  (* Simplified 2-loop sunrise-like with combined propagator *)
  props = {
    -k1^2 + 2*xx*k1*p - xx^2*p^2 + xx*(1-xx)*p^2,   (* combined: xx*(-k1^2) + (1-xx)*(-(k1-p)^2) *)
    -(k1-p)^2,    (* = -k1^2 + 2*k1*p - p^2 *)
    -k2^2,
    -(k1-k2)^2    (* = -k1^2 + 2*k1*k2 - k2^2 *)
  };
  loops = {k1, k2};

  (* Simplify props[1] *)
  props[[1]] = Expand[xx*(-k1^2) + (1-xx)*(-(k1-p)^2)];
  props[[2]] = Expand[-(k1-p)^2];
  props[[4]] = Expand[-(k1-k2)^2];

  Print["  props[1] = ", props[[1]]];
  Print["  props[2] = ", props[[2]]];

  {coeffMat, constVec} = FeynmanTrick`PropagatorAlgebra`DecomposePropagatorDerivative[
    props, loops, xx
  ];

  Print["  CoeffMatrix row 1: ", coeffMat[[1]]];
  Print["  ConstVector[1]: ", constVec[[1]]];

  (* Verify reconstruction *)
  deriv1 = Expand[D[props[[1]], xx]];
  recon1 = Expand[Sum[coeffMat[[1, j]] * props[[j]], {j, 4}] + constVec[[1]]];

  testApprox["2-loop reconstruction matches derivative",
    Simplify[recon1 - deriv1], 0
  ];

  (* Expected: coeffMat[1,1] = 1/xx, coeffMat[1,2] = -1/xx *)
  testApprox["2-loop CoeffMat[1,1] = 1/xx", coeffMat[[1, 1]], 1/xx];
  testApprox["2-loop CoeffMat[1,2] = -1/xx", coeffMat[[1, 2]], -1/xx];
];

(* ============================================================ *)
(* Test 7: EpsPrefactors *)
(* ============================================================ *)

Print["\n--- Test: EpsPrefactors ---"];

Module[{testMatrix, prefactors, transformed, eps = FTeps},
  (* Matrix with an eps pole in entry [1,2] *)
  testMatrix = {
    {1/xx, 1/(eps*xx)},
    {0, -1/xx}
  };

  Print["  eps symbol: ", eps, " (context: ", Context[eps], ")"];
  Print["  Has pole: ", FeynmanTrick`EpsPrefactors`CheckEpsPoles[testMatrix, eps]];

  prefactors = FeynmanTrick`EpsPrefactors`FindEpsPrefactors[testMatrix, eps];
  Print["  Prefactors: ", prefactors];

  (* Expected: k1-k2 >= 1, k1=0 → k2 <= -1, so k2=-1 works *)
  test["Prefactor removes pole", prefactors[[1]] - prefactors[[2]] >= 1, True];

  transformed = FeynmanTrick`EpsPrefactors`ApplyEpsPrefactors[testMatrix, prefactors, eps];
  Print["  Transformed matrix: ", transformed];

  (* Check no poles remain *)
  test["No poles after transform",
    FeynmanTrick`EpsPrefactors`CheckEpsPoles[transformed, eps], False
  ];
];

(* ============================================================ *)
(* Summary *)
(* ============================================================ *)

Print["\n============================"];
Print["Results: ", passed, " passed, ", failed, " failed."];
Print["============================"];

If[failed > 0, Exit[1]];
