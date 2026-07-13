(* Shared FeynmanTrick example definitions for the Scripts/ tools.
   Load FeynmanTrick/FeynmanTrick.m before Get-ing this file.

   All kinematics are baked in numerically (Euclidean region):
     - bubble/sunrise/banana/banana4/kite: massive lines (m=1),
       p^2 = -1, d = 2-2eps
     - banana_unequal: squared masses {2,3/2,4/3,1},
       p^2 = -1, d = 2-2eps
     - banana4_unequal: squared masses {2,3/2,4/3,5/4,1},
       p^2 = -1, d = 2-2eps
     - box family + pentagon: massless on-shell legs, d = 4-2eps
       box/box_bubble/box_triangle/double_box_planar: s = -1, t = -1/3
       pentagon: {s12,s23,s34,s45,s15} = {-1,-2,-3,-5,-7}
     - pentagon_massive: squared internal masses {1,3/2,4/3,5/4,6/5},
       p_i^2 = -1 and p_i.p_j = 1/4 (i != j), d = 4-2eps

   box_bubble drops double-box propagators 1 and 2; p1 and p2 then enter at
   the same vertex, so it is kinematically a 2-loop THREE-point function
   (depends on s only).  box_triangle drops double-box propagator 1 and
   remains a genuine 4-point function of s and t. *)

FTExampleNames[] := FeynmanTrick`SupportedExamples[];

FTExampleDoubleBoxPropagators[] := {
  -Global`l1^2,
  -(Global`l1 + Global`p1)^2,
  -(Global`l1 + Global`p1 + Global`p2)^2,
  -(Global`l1 - Global`l2)^2,
  -Global`l2^2,
  -(Global`l2 + Global`p1 + Global`p2)^2,
  -(Global`l2 + Global`p1 + Global`p2 + Global`p3)^2
};

(* s = -1, t = -1/3:  p1.p2 = s/2, p2.p3 = t/2, p1.p3 = -(s+t)/2 *)
FTExampleBoxReplacements[] := {
  Global`p1^2 -> 0, Global`p2^2 -> 0, Global`p3^2 -> 0,
  Global`p1 Global`p2 -> -1/2,
  Global`p2 Global`p3 -> -1/6,
  Global`p1 Global`p3 -> 2/3
};

(* {s12,s23,s34,s45,s15} = {-1,-2,-3,-5,-7}; p5 = -(p1+p2+p3+p4) is
   on-shell by p1.p4 = 5 (sum of all off-diagonal dot products vanishes). *)
FTExamplePentagonReplacements[] := {
  Global`p1^2 -> 0, Global`p2^2 -> 0, Global`p3^2 -> 0, Global`p4^2 -> 0,
  Global`p1 Global`p2 -> -1/2,
  Global`p2 Global`p3 -> -1,
  Global`p3 Global`p4 -> -3/2,
  Global`p1 Global`p3 -> -1,
  Global`p2 Global`p4 -> -1,
  Global`p1 Global`p4 -> 5
};

(* Symmetric, strictly Euclidean four-vector point.  With
   p5=-(p1+p2+p3+p4), p5^2=-1 as well.  The Euclidean Gram matrix -p_i.p_j
   has eigenvalues {1/4,5/4,5/4,5/4}. *)
FTExampleMassivePentagonReplacements[] := {
  Global`p1^2 -> -1, Global`p2^2 -> -1,
  Global`p3^2 -> -1, Global`p4^2 -> -1,
  Global`p1 Global`p2 -> 1/4,
  Global`p1 Global`p3 -> 1/4,
  Global`p1 Global`p4 -> 1/4,
  Global`p2 Global`p3 -> 1/4,
  Global`p2 Global`p4 -> 1/4,
  Global`p3 Global`p4 -> 1/4
};

FTExampleSpec[name_String] := Switch[name,
  "bubble",
    <|
      "LoopMomenta" -> {Global`l1},
      "ExternalMomenta" -> {Global`p},
      "Propagators" -> {1 - Global`l1^2, 1 - (Global`l1 - Global`p)^2},
      "Replacements" -> {Global`p^2 -> -1},
      "Dimension" -> 2 - 2*FeynmanTrick`FTeps
    |>,
  "sunrise",
    <|
      "LoopMomenta" -> {Global`l1, Global`l2},
      "ExternalMomenta" -> {Global`p},
      "Propagators" -> {
        1 - Global`l1^2,
        1 - Global`l2^2,
        1 - (-Global`l1 - Global`l2 + Global`p)^2
      },
      "Replacements" -> {Global`p^2 -> -1},
      "Dimension" -> 2 - 2*FeynmanTrick`FTeps
    |>,
  "banana",
    <|
      "LoopMomenta" -> {Global`l1, Global`l2, Global`l3},
      "ExternalMomenta" -> {Global`p},
      "Propagators" -> {
        1 - Global`l1^2,
        1 - Global`l2^2,
        1 - Global`l3^2,
        1 - (-Global`l1 - Global`l2 - Global`l3 + Global`p)^2
      },
      "Replacements" -> {Global`p^2 -> -1},
      "Dimension" -> 2 - 2*FeynmanTrick`FTeps
    |>,
  (* The constants are squared masses: FIRE and the FT algebra use the
     propagator convention D_j = m_j^2 - q_j^2.  This is the Euclidean
     point mm={2,3/2,4/3,1}, psq=-1 of the legacy unequal-mass family. *)
  "banana_unequal",
    <|
      "LoopMomenta" -> {Global`l1, Global`l2, Global`l3},
      "ExternalMomenta" -> {Global`p},
      "Propagators" -> {
        2 - Global`l1^2,
        3/2 - Global`l2^2,
        4/3 - Global`l3^2,
        1 - (-Global`l1 - Global`l2 - Global`l3 + Global`p)^2
      },
      "Replacements" -> {Global`p^2 -> -1},
      "Dimension" -> 2 - 2*FeynmanTrick`FTeps
    |>,
  "banana4",
    <|
      "LoopMomenta" -> {Global`l1, Global`l2, Global`l3, Global`l4},
      "ExternalMomenta" -> {Global`p},
      "Propagators" -> {
        1 - Global`l1^2,
        1 - Global`l2^2,
        1 - Global`l3^2,
        1 - Global`l4^2,
        1 - (-Global`l1 - Global`l2 - Global`l3 - Global`l4 + Global`p)^2
      },
      "Replacements" -> {Global`p^2 -> -1},
      "Dimension" -> 2 - 2*FeynmanTrick`FTeps
    |>,
  (* Four-loop unequal-mass Euclidean showcase.  As above, these constants
     are squared masses in D_j = m_j^2-q_j^2. *)
  "banana4_unequal",
    <|
      "LoopMomenta" -> {Global`l1, Global`l2, Global`l3, Global`l4},
      "ExternalMomenta" -> {Global`p},
      "Propagators" -> {
        2 - Global`l1^2,
        3/2 - Global`l2^2,
        4/3 - Global`l3^2,
        5/4 - Global`l4^2,
        1 - (-Global`l1 - Global`l2 - Global`l3 - Global`l4 + Global`p)^2
      },
      "Replacements" -> {Global`p^2 -> -1},
      "Dimension" -> 2 - 2*FeynmanTrick`FTeps
    |>,
  (* Fully massive equal-mass two-loop kite.  This is the symmetric
     Euclidean benchmark, rather than the separate three-mass/two-massless
     four-dimensional kite convention. *)
  "kite",
    <|
      "LoopMomenta" -> {Global`l1, Global`l2},
      "ExternalMomenta" -> {Global`p},
      "Propagators" -> {
        1 - Global`l1^2,
        1 - (Global`l1 - Global`p)^2,
        1 - Global`l2^2,
        1 - (Global`l2 - Global`p)^2,
        1 - (Global`l1 - Global`l2)^2
      },
      "Replacements" -> {Global`p^2 -> -1},
      "Dimension" -> 2 - 2*FeynmanTrick`FTeps
    |>,
  "box",
    <|
      "LoopMomenta" -> {Global`l1},
      "ExternalMomenta" -> {Global`p1, Global`p2, Global`p3},
      "Propagators" -> {
        -Global`l1^2,
        -(Global`l1 + Global`p1)^2,
        -(Global`l1 + Global`p1 + Global`p2)^2,
        -(Global`l1 + Global`p1 + Global`p2 + Global`p3)^2
      },
      "Replacements" -> FTExampleBoxReplacements[],
      "Dimension" -> 4 - 2*FeynmanTrick`FTeps
    |>,
  "pentagon",
    <|
      "LoopMomenta" -> {Global`l1},
      "ExternalMomenta" -> {Global`p1, Global`p2, Global`p3, Global`p4},
      "Propagators" -> {
        -Global`l1^2,
        -(Global`l1 + Global`p1)^2,
        -(Global`l1 + Global`p1 + Global`p2)^2,
        -(Global`l1 + Global`p1 + Global`p2 + Global`p3)^2,
        -(Global`l1 + Global`p1 + Global`p2 + Global`p3 + Global`p4)^2
      },
      "Replacements" -> FTExamplePentagonReplacements[],
      "Dimension" -> 4 - 2*FeynmanTrick`FTeps
    |>,
  "pentagon_massive",
    <|
      "LoopMomenta" -> {Global`l1},
      "ExternalMomenta" -> {
        Global`p1, Global`p2, Global`p3, Global`p4},
      "Propagators" -> {
        1 - Global`l1^2,
        3/2 - (Global`l1 + Global`p1)^2,
        4/3 - (Global`l1 + Global`p1 + Global`p2)^2,
        5/4 - (Global`l1 + Global`p1 + Global`p2 + Global`p3)^2,
        6/5 - (Global`l1 + Global`p1 + Global`p2 +
          Global`p3 + Global`p4)^2
      },
      "Replacements" -> FTExampleMassivePentagonReplacements[],
      "Dimension" -> 4 - 2*FeynmanTrick`FTeps
    |>,
  "box_bubble",
    <|
      "LoopMomenta" -> {Global`l1, Global`l2},
      "ExternalMomenta" -> {Global`p1, Global`p2, Global`p3},
      "Propagators" -> FTExampleDoubleBoxPropagators[][[{3, 4, 5, 6, 7}]],
      "Replacements" -> FTExampleBoxReplacements[],
      "Dimension" -> 4 - 2*FeynmanTrick`FTeps
    |>,
  "box_triangle",
    <|
      "LoopMomenta" -> {Global`l1, Global`l2},
      "ExternalMomenta" -> {Global`p1, Global`p2, Global`p3},
      "Propagators" -> FTExampleDoubleBoxPropagators[][[{2, 3, 4, 5, 6, 7}]],
      "Replacements" -> FTExampleBoxReplacements[],
      "Dimension" -> 4 - 2*FeynmanTrick`FTeps
    |>,
  "double_box_planar",
    <|
      "LoopMomenta" -> {Global`l1, Global`l2},
      "ExternalMomenta" -> {Global`p1, Global`p2, Global`p3},
      "Propagators" -> FTExampleDoubleBoxPropagators[],
      "Replacements" -> FTExampleBoxReplacements[],
      "Dimension" -> 4 - 2*FeynmanTrick`FTeps
    |>,
  _,
    Print["FTExampleSpec: unknown example ", name]; $Failed
];

FTExampleTopology[name_String, prefix_String:"ft"] := Module[
  {spec = FTExampleSpec[name]},
  If[spec === $Failed, Return[$Failed, Module]];
  FeynmanTrick`FIREInterface`DefineTopology[
    prefix <> "_" <> name,
    spec["LoopMomenta"],
    spec["ExternalMomenta"],
    spec["Propagators"],
    spec["Replacements"]
  ]
];

(* Combine everything into position 1, left to right - the convention used
   by every verified example so far. *)
FTExampleSequence[name_String] := Module[{spec = FTExampleSpec[name]},
  If[spec === $Failed, Return[$Failed, Module]];
  Table[{1, k}, {k, 2, Length[spec["Propagators"]]}]
];

FTExampleDimension[name_String] := Module[{spec = FTExampleSpec[name]},
  If[spec === $Failed, Return[$Failed, Module]];
  spec["Dimension"]
];
