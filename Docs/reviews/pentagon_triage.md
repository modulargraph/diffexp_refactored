# Pentagon triage (RewritePlan M0 task 16; risks R8/R9)

Status: COMPLETE, 2026-06-11.  Inputs: `/tmp/f_pentagon.log` (completed
pentagon run, FAILS pin at all orders; copy frozen at
`Tests/refs/oracle_logs/f_pentagon.log`), `Scripts/FTExamples.m`,
`Scripts/run_ft_stepwise.m`, `FeynmanTrick/DiffExpIntegration.m`,
`DiffExp/Transport.m`, `DiffExp/AnalyticContinuation.m`,
`DiffExp/RegularizedIntegration.m`, `Tests/PINS.md`,
`Docs/reviews/rewrite_plan_review_3lens.json` (key `result.execution`,
finding 1), and the surviving level-1 integral replay dump
`/tmp/diffexp_pentagon_dumps/laurent_integral_0006.m` (smoking-gun
evidence; vendor it — see section 6).  No Wolfram kernel was run for this
triage; every numerical claim below is verified to >= 150 digits against
the dump files with exact arithmetic.

VERDICT IN ONE PARAGRAPH: the "current point is not recognized as a branch
point" warning (fired twice, log lines 72-73) is caused by ONE missing
prescription factor at pentagon level 1: the three-mass-triangle Kallen
factor `16 xx1^2 - 19 xx1 + 4`, whose roots `(19 +- Sqrt[105])/32`
(~0.273533 and ~0.913967) are sqrt-type branch points INTERIOR to the
integration interval [0,1] — one on each transport leg, hence exactly two
warnings.  The factor lives only in the level-1 DE-matrix denominators;
the FT pipeline's prescription auto-generation scans only IBP-coefficient
denominators and so never sees it.  Without it DiffExp crosses two genuine
branch points with NO analytic-continuation rule, corrupting the level-1
transport at ALL eps orders (the observed spurious `Im(eps^0) = -0.319` in
a real Euclidean quantity is its fingerprint).  This is the
prescription-config failure class: it survives the DiffExp2 rewrite
unchanged unless the example data is patched (section 2).  Independently,
the run contains 20 silent zero-regulator drops at level 2 and 2x19
omitted endpoint coefficients at level 1 (D2 class — the rewrite fixes
these by construction).  The R9 sign protocol confirms the recorded pin
conventions are internally consistent; the failure is NOT a pin-sign
artifact (section 3).

---

## 1. Root cause of the unrecognized-branch-point warning

### 1.1 The warning's exact source and trigger chain

- Warning text emitted at `DiffExp/Transport.m:472` — the `PrintWarning`
  branch of `GiveMultivaluedError[]` (defined `DiffExp/Transport.m:470-478`).
  The warning (not abort) branch is taken because the FT layer sets
  `"AbortOnAnalyticContinuationFail" -> False`
  (`FeynmanTrick/DiffExpIntegration.m:411-414`, re-applied at `:466-469`).
- `GiveMultivaluedError[]` is invoked from the `SingularityCheck` block at
  `DiffExp/Transport.m:480-492`, specifically the first condition at
  `:485-487`: the boundary-fixed solutions on the current segment contain
  `Logx` (`:482`) or algebraic roots `x^b` with `Denominator[b] > 1`
  (`:483`), AND `CurrentSingularityHasIDeltaPrescription === False`.
- That flag is computed in
  `DiffExp/AnalyticContinuation.m:39-43` (`PrepareAnalyticContinuation`,
  called per segment from `DiffExp/Transport.m:211`; also `:47`, `:639`,
  `:1025`): a prescription "recognizes" the current chart iff one of the
  `DeltaPrescriptions` polynomials vanishes at the segment center
  (`FindVanishingFactors`, `AnalyticContinuation.m:21-27`).  If none
  vanishes, `VanishingFactors === {}`, the flag is `False`, AND
  `AnalyticContinuationReplacements = {}` (`AnalyticContinuation.m:70-79`
  never populates it) — i.e. the segment is crossed with NO `Logx`/
  `x^(p/q)` continuation rules at all.  The crossing is not merely
  warned about; it is performed wrongly (principal-branch values on the
  far side, no theta-function bookkeeping).
- The check runs inside `IntegrateSystem` (`DiffExp/Transport.m:146`),
  which `TransportTo` calls once per segment
  (`DiffExp/Transport.m:843`, `:1058`).

### 1.2 The transport line and the two singular points

The warnings fire during the LEVEL-1 stage of the run (log line 71 prints
`EXTRA_SINGULAR_FACTORS level 1` immediately before them).
`Scripts/run_ft_stepwise.m:154-231` iterates levels 4 -> 1; for each level
`TransportLevel` (`FeynmanTrick/DiffExpIntegration.m:267-`) performs TWO
`TransportTo` calls in the detected variable `xx1`:

- leg A: `<|xx1 -> 11/23|>` -> `<|xx1 -> 0|>`
  (`FeynmanTrick/DiffExpIntegration.m:438-446`),
- leg B: `<|xx1 -> 11/23|>` -> `<|xx1 -> 1|>`
  (`FeynmanTrick/DiffExpIntegration.m:472-480`).

The level-1 integral replay dump
`/tmp/diffexp_pentagon_dumps/laurent_integral_0006.m` (`"Variable" -> xx1`)
records every segment chart of that line.  Among its chart centers are two
quadratic irrationals, matching to 159 digits (verified with 165-digit
decimal arithmetic):

- `xx1 = (19 - Sqrt[105])/32 = 0.27353278856376255052...`  (crossed by leg A)
- `xx1 = (19 + Sqrt[105])/32 = 0.91396721143623744947...`  (crossed by leg B)

These are the roots of `16 xx1^2 - 19 xx1 + 4`.  Each leg crosses exactly
one of them => exactly two warnings (log lines 72-73).  All other interior
chart centers in the dump are rational predivision points; `xx1 = 1/2` also
appears as a chart center (see 1.4) but produced no warning.

### 1.3 Why that factor: full reconstruction of the level systems

Pentagon (`Scripts/FTExamples.m:95-108`): 1-loop, 5 massless on-shell legs,
propagators `-l1^2, -(l1+p1)^2, ..., -(l1+p1+p2+p3+p4)^2`
(`FTExamples.m:99-105`), kinematics `{s12,s23,s34,s45,s15} =
{-1,-2,-3,-5,-7}` via the dot products at `FTExamples.m:38-46`.
Combination sequence (`FTExamples.m:151-154`): `{{1,2},{1,3},{1,4},{1,5}}`,
i.e. level k combines the running position-1 propagator with propagator
k+1 using Feynman parameter `xx_k`; deeper-but-unintegrated parameters are
frozen at 11/23 in each level's matrix export.  The empirical weight
convention (pinned below by the IBP factors) is
`P_combined = xx_k * P_1 + (1 - xx_k) * P_(k+1)`.

Level-by-level singular factors (EXTRA_SINGULAR_FACTORS = output of
`CollectLevelIBPSingularFactors`, which scans ONLY the denominators of the
boundary-request IBP reduction coefficients,
`FeynmanTrick/DiffExpIntegration.m:1326-1384`, line `:1369` —
never the exported DE matrices):

| level | var | IBP factors (log line) | roots | interior (0,1)? | warnings at this level |
|---|---|---|---|---|---|
| 4 | xx4 | `{xx4, -38617 + 27958 xx4}` (line 25) | 0, 38617/27958 ~ 1.3813 | none | none |
| 3 | xx3 | `{xx3, -1817 + 1685 xx3}` (line 28) | 0, 1817/1685 ~ 1.0783 | none | 3x `General::munfl` (cosmetic underflow, lines 31-39) |
| 2 | xx2 | `{xx2, -69 + 64 xx2, 23 + 5 xx2, -1 + xx2}` (line 43) | 0, 69/64, -23/5, 1 | none | 20x zero-regulator drop + 1x omitted (>= eps^17) (lines 44-64) |
| 1 | xx1 | `{2 + 3 xx1, -1 + xx1, xx1}` (line 71) | -2/3, 1, 0 | none from IBP | 2x branch-point config (lines 72-73), 2x omitted-19 (>= eps^0) (lines 74-75) |

(The 23-laden integers at levels 2-4 are the frozen `xx_j = 11/23`
residues; the outside-[0,1] roots 38617/27958 and 1817/1685 are visible as
chart-map tokens in dumps 0001/0002 — consistency check passed.  Chart
centers at levels 2-4, extracted from dumps 0001-0005, contain NO
irrational interior points, consistent with zero branch-point warnings
there.)

Level-1 system, derived from `FTExamples.m` kinematics.  With
`P12 = xx1 * (-l1^2) + (1-xx1) * (-(l1+p1)^2) = -(l1 + (1-xx1) p1)^2`
(exact since `p1^2 = 0`), shift `l1' = l1 + (1-xx1) p1`: the level-1 family
is a one-loop TWO-MASS-HARD BOX with external legs

```
q1 = xx1 p1 + p2          m1^2 = q1^2 = xx1 s12      = -xx1
q2 = p3                    massless
q3 = p4                    massless
q4 = (1-xx1) p1 + p5      m4^2 = q4^2 = (1-xx1) s15  = 7 xx1 - 7
S  = (q1+q2)^2 = -3 xx1 - 2          T = (q2+q3)^2 = s34 = -3
```

This assignment is PINNED by the IBP factor list itself: the level-1
masters (log lines 65-70) are `[1,0,1,1,1]` (the 2mh box), `[1,0,1,0,1]`
(three-mass triangle on legs `{m1^2, T, m4^2}`), and bubbles
`[1,0,1,0,0]` (scale `m1^2 = -xx1` -> factor `xx1`), `[1,0,0,0,1]`
(scale `m4^2 = 7(xx1-1)` -> factor `-1+xx1`), `[1,0,0,1,0]` (scale
`S = -(2+3 xx1)` -> factor `2+3 xx1`), `[0,0,1,0,1]` (scale `T = -3`,
constant).  Every IBP factor is a bubble scale; the alternative weight
convention (`(1-xx1)` on prop 1) is excluded because it predicts a factor
`5 - 3 xx1` instead of `2 + 3 xx1`.

The complete Landau alphabet of this family adds the DE-matrix-only
factors:

- THE MISSING FACTOR — three-mass-triangle Kallen determinant
  (`[1,0,1,0,1]` subsector):
  `lambda(m1^2, T, m4^2) = lambda(-xx1, -3, 7 xx1 - 7)
   = 64 xx1^2 - 76 xx1 + 16 = 4 (16 xx1^2 - 19 xx1 + 4)`,
  roots `(19 +- Sqrt[105])/32 ~= 0.273533, 0.913967` — BOTH interior to
  [0,1], sqrt-type branch points (the local solution space at a Kallen
  zero carries half-integer exponents -> `AlgebraicRootsPresent` at
  `Transport.m:483`).  Identity `lambda = 4(16x^2-19x+4)` and the root
  values verified exactly.
- `S - m4^2 = 5 - 10 xx1`, root `xx1 = 1/2` (degenerate two-mass-triangle
  Landau point, squared Kallen `(S-m4^2)^2`): present as a chart center
  in dump 0006 but the fixed solutions there are single-valued (the
  two-mass triangle is analytic at `S = m4^2`) -> no warning, no
  prescription needed.
- `S - m1^2 = -2 xx1 - 2` (root -1) and box leading singularity
  `S T - m1^2 m4^2 = 7 xx1^2 + 2 xx1 + 6` (complex roots — the
  complex-singularity class RewritePlan section 3.1 flags for pentagon
  lines): never on [0,1].

### 1.4 Why the prescription list misses it

The effective `DeltaPrescriptions` for the level-1 transports are
auto-generated by `deltaPrescriptionsForFactors`
(`FeynmanTrick/DiffExpIntegration.m:224-238`, applied at `:398-408`):
endpoints `{detectedVar, 1 - detectedVar}` plus the IBP factors from
`CollectLevelIBPSingularFactors`, all with sign +1.  For pentagon level 1
that yields exactly

```
{{xx1, 1}, {1 - xx1, 1}, {2 + 3 xx1, 1}}
```

(`-1 + xx1` is folded into `1 - xx1` by the sign-blind dedup at
`DiffExpIntegration.m:234-237`; see section 6c).  None of these vanishes at
`(19 +- Sqrt[105])/32` -> `CurrentSingularityHasIDeltaPrescription = False`
-> warning, and the crossing happens with empty replacement rules.

Two structural reasons the factor cannot currently get in:

1. `CollectLevelIBPSingularFactors` scans IBP boundary-reduction
   coefficient denominators only (`DiffExpIntegration.m:1359-1380`); the
   Kallen factor appears only in the d/dxx1 DE-matrix denominators.
   DiffExp itself DOES know the factor — `MatricesIrreducibleFactors`
   is built from the loaded matrix denominators
   (`DiffExp/MatrixLoading.m:217-231`) and segmentation correctly centers
   charts on its roots (`DiffExp/LineSegmentation.m:65-71`, hence the
   irrational chart centers in the dump) — but matrix factors are NOT
   propagated to prescriptions unless `"PrescribeMatrixFactors" -> True`
   (default `False`, `DiffExpIntegration.m:281`, used at `:400-406`),
   which `Scripts/run_ft_stepwise.m:183-192` does not pass.
2. The old auto-prescription for sqrt singularities
   (`SquareRootPrescriptionsAdded`, `DiffExp/State.m:228-230`, applied at
   `MatrixLoading.m:189-190`) keys on literal `Sqrt[...]` in matrix
   ENTRIES (`DEqnSquareRoots`).  The FIRE-exported FT matrices are
   rational in `xx1`; the sqrt lives only in the SOLUTIONS, so this
   mechanism is empty here.

NOTE / premise correction: `Scripts/FTExamples.m` contains NO
`DeltaPrescriptions` data today (grep confirms; the only DeltaPrescriptions
in the FT path are the auto-generated ones above).  The task brief's
"pentagon's DeltaPrescriptions in FTExamples.m" therefore maps to "the
absence of any per-example prescription surface in FTExamples.m" — the
patch below creates it.

### 1.5 Consequence chain (why this corrupts ALL orders)

Crossing a genuine sqrt branch point with
`AnalyticContinuationReplacements = {}` is an O(1) RELATIVE error on every
sector transported through the chart — branch errors are not
eps-suppressed.  Both legs of the level-1 line cross one such point, so
all level-1 segment data outboard of `0.2735` / `0.9140` (in particular
the data feeding the `[0,1]` integration and both endpoint limits) is on
an uncontrolled branch at every eps order, including the ones that build
the final `eps^-2`.  Fingerprint: the final result (log line 76) has
`Im(eps^0) = -0.3193` although the Euclidean-region pentagon (all
`s_ij < 0`) is real; the level-1 BCs as produced by level 2 (log line 65)
still had `Im ~ 1.5e-19`.  The imaginary contamination enters exactly at
the level-1 stage, where the two warnings fire.  A pure-D2 explanation
does not produce this: the drop rules discard real Laurent content but do
not inject O(1) imaginary parts, and the honest omission warnings (lines
74-75) only disclaim orders >= eps^0 while the run is wrong already at
eps^-2 (0.25214 vs pin-implied +0.30952 = +13/42).

---

## 2. Proposed patch (NOT applied — for the M5 pentagon retest)

Two hunks: the example-data surface in `Scripts/FTExamples.m` (primary,
as tasked) and the 3-line consumption hook in `Scripts/run_ft_stepwise.m`
(without which the new function is dead code).  Routing the factor through
`"ExtraSingularFactors"` is deliberate: `TransportLevel` then both adds it
to the segmentation alphabet (`appendMatrixFactors`,
`FeynmanTrick/DiffExpIntegration.m:204-222` — harmless, the factor is
already in the matrix alphabet) and emits a `{factor, +1}` prescription
(`deltaPrescriptionsForFactors`, `:398-408`), with the name-based variable
remap at `:374-378` making the `Global`xx1` spelling context-safe (the
level-1 detected variable is literally named `xx1`, per the dump).

### 2.1 `Scripts/FTExamples.m` (append after `FTExampleDimension`, line 159)

```diff
--- a/Scripts/FTExamples.m
+++ b/Scripts/FTExamples.m
@@ -156,3 +156,32 @@
 FTExampleDimension[name_String] := Module[{spec = FTExampleSpec[name]},
   If[spec === $Failed, Return[$Failed, Module]];
   spec["Dimension"]
 ];
+
+(* Per-level DE-matrix singular factors that need explicit i*delta
+   prescriptions but are invisible to CollectLevelIBPSingularFactors
+   (FeynmanTrick/DiffExpIntegration.m:1326-1384 scans only IBP
+   boundary-reduction coefficient denominators, never the exported
+   d<var> matrices).  Join these into the "ExtraSingularFactors" passed
+   to TransportLevel: they then enter both the segmentation alphabet
+   (appendMatrixFactors) and the auto-generated DeltaPrescriptions
+   (deltaPrescriptionsForFactors, DiffExpIntegration.m:398-408, sign +1).
+   Factors must be written in the level's Feynman parameter symbol NAME
+   (xx<level>); TransportLevel remaps by SymbolName to the detected
+   matrix variable (DiffExpIntegration.m:374-378).
+
+   pentagon, level 1: the level-1 system is the one-loop two-mass-hard
+   box with legs q1 = xx1 p1 + p2, q2 = p3, q3 = p4, q4 = (1-xx1) p1 + p5
+   (props 1,2 combined with weight xx1 on prop 1).  Its three-mass
+   triangle subsector [1,0,1,0,1] carries the Kallen factor
+   lambda(q1^2, s34, q4^2) = lambda(-xx1, -3, 7 xx1 - 7)
+   = 4 (16 xx1^2 - 19 xx1 + 4), roots (19 +- Sqrt[105])/32
+   ~= 0.273533 and 0.913967 — sqrt-type branch points INTERIOR to the
+   integration interval [0,1].  Without a prescription DiffExp cannot
+   continue across them ("current point is not recognized as a branch
+   point", DiffExp/Transport.m:472) and the level-1 transport is wrong
+   at ALL eps orders.  Full derivation: Docs/reviews/pentagon_triage.md. *)
+FTExampleLevelDeltaFactors[name_String, level_Integer] := Switch[
+  {name, level},
+  {"pentagon", 1}, {16*Global`xx1^2 - 19*Global`xx1 + 4},
+  _, {}
+];
```

### 2.2 `Scripts/run_ft_stepwise.m` (companion hook, lines 177-180)

```diff
--- a/Scripts/run_ft_stepwise.m
+++ b/Scripts/run_ft_stepwise.m
@@ -174,10 +174,13 @@
     ];
     matrixDir = FileNameJoin[{outputDir, "Level_" <> ToString[level] <> "_Matrices"}];
-    extraSingularFactors =
-      FeynmanTrick`DiffExpIntegration`CollectLevelIBPSingularFactors[
-        ftData, level
-      ];
+    extraSingularFactors = Join[
+      FeynmanTrick`DiffExpIntegration`CollectLevelIBPSingularFactors[
+        ftData, level
+      ],
+      FTExampleLevelDeltaFactors[name, level]
+    ];
     Print["EXTRA_SINGULAR_FACTORS level ", level, ": ",
       InputForm[extraSingularFactors]];
```

### 2.3 Patch caveats and acceptance criteria

- SIGN: `deltaPrescriptionsForFactors` assigns the pipeline-wide default
  `+1` (`DiffExpIntegration.m:280`).  For an interior Euclidean
  pseudo-threshold the physical answer is real, which gives a sharp
  acceptance test rather than an a-priori sign argument.  Accept the
  patch iff, on rerun: (a) the two `Transport.m:472` warnings are gone;
  (b) `|Im(eps^0)|` of FINAL collapses from 0.319 to ~0; (c) `eps^-2`
  moves toward `+13/42 = +0.30952380952...` (the exact closed-form
  anchor, `Tests/PINS.md:117-118`).  If a substantial imaginary part
  SURVIVES, the prescription sign is wrong for this factor: flip it
  (requires a small extension, since `deltaPrescriptionsForFactors`
  takes one sign for all factors — add a per-factor `{poly, sign}`
  passthrough, or set `"DeltaPrescriptionSign" -> -1` for a diagnostic
  run, noting that flips ALL prescriptions).
- A single prescription on the irreducible quadratic covers both roots
  with mutually consistent chart Im-signs (DiffExp derives the per-chart
  sign from the leading coefficient, `AnalyticContinuation.m:45-50`);
  do NOT split it into two linear sqrt factors.
- The blanket alternative `"PrescribeMatrixFactors" -> True` (one
  option, no FTExamples change) would also prescribe `5 - 10 xx1` (an
  analytic point — harmless in principle but a behavior change) and
  every other matrix factor at every level FOR ALL EXAMPLES if set
  globally.  Rejected: the targeted data patch is auditable per example
  and cannot disturb the five validated ladder entries.
- DiffExp2 forward-compatibility: the same factor must reach DiffExp2's
  `LocalSolution["Prescriptions"]` for the pentagon M5 run.  The
  FTExamples surface is engine-agnostic data, which is the point of
  putting it there.

---

## 3. R9 sign-protocol result

Rule under test (`Tests/PINS.md:42`):
`FT = (-1)^(sum of propagator powers) x loop_package`
(= `(-1)^nprops` here; all ladder examples use corner integrals, powers 1).
Origin of the factor: FT propagators are Euclidean `-(q^2)`
(`Scripts/FTExamples.m:18-26,99-105`) while the pySecDec loop_package
route uses `+(q^2)` (e.g.
`Tests/refs/generators/pysecdec_pentagon.py:20-23`); each propagator flip
contributes one sign.  Cross-check: pySecDec's printed pentagon prefactor
`-gamma(eps+3)` = `(-1)^5 Gamma[3+eps]`
(`Tests/refs/pysecdec/pentagon_pin_result.txt:1`) carries the same
counting.

| example | props (source) | predicted sign | recorded status | consistent? |
|---|---|---|---|---|
| box | 4 (`FTExamples.m:86-91`) | +1 | VALIDATED end-to-end vs `box1l` pin (`Tests/PINS.md:46,64-77`; anchor protocol `Tests/refs/generators/pysecdec_box1l_pin.py`) | YES |
| box_bubble | 5 (`FTExamples.m:113`, dbox props `[[{3,4,5,6,7}]]`) | -1 | VALIDATED at all 4 orders as `FT = -1 x loop_package` (`Tests/PINS.md:47,79-93`; memory brief) | YES |
| pentagon | 5 (`FTExamples.m:99-105`) | **-1** | recorded as `-1`, PROVISIONAL (`Tests/PINS.md:48,105-128`; memory brief stores loop_package values `{-0.309523809524, +0.0164223473, +0.0254117770, -0.2797210931}` to be flipped) | YES |
| box_triangle | 6 (`FTExamples.m:121`, dbox props `[[{2,3,4,5,6,7}]]`) | **+1** | recorded as `+1` (`Tests/PINS.md:49,130`; memory brief "box_triangle loop_package (sign +1)") | YES |
| double_box_planar | 7 (`FTExamples.m:129`) | -1 | recorded as `-1`, never integrated (`Tests/PINS.md:50`) | YES |

FINDINGS:

1. NO INCONSISTENCY between the rule, the propagator counts in
   `FTExamples.m`, the memory brief, and `Tests/PINS.md`.  Derived
   pentagon sign: `(-1)^5 = -1`; FT pin = `{+0.309523809524,
   -0.0164223473, -0.0254117770, +0.2797210931}` with the exact
   leading-pole anchor `FT eps^-2 = +13/42`.  Derived box_triangle sign:
   `(-1)^6 = +1`.
2. The pentagon failure is NOT a sign/normalization artifact.  Two
   independent arguments: (a) the run's `eps^-2 = +0.25214` agrees in
   SIGN with the `-1`-rule expectation `+0.30952` and disagrees with the
   raw loop_package value `-0.30952` — corroborating the `-1`; (b) a
   wrong global sign or normalization scales all orders by one constant,
   but the observed ratios are order-dependent (`eps^-2`: 0.815;
   `eps^-1`: `-0.17471 / -0.01642` = 10.6; `eps^0` complex) — genuine
   order-dependent content corruption, exactly what sections 1 and 4
   describe.
3. WHAT REMAINS OPEN (why the pin stays PROVISIONAL): the `-1` for
   pentagon is predicted, not yet empirically pinned, because the FT
   side has never produced a correct value to take a ratio with.  The
   closing step is the PINS.md section-1 protocol run
   (`Tests/PINS.md:52-58`): compute the same integral through
   `Scripts/export_pysecdec_family_specs.m` +
   `Scripts/pysecdec_family_driver.py` (exact FT normalization, no sign
   ambiguity) and confirm the ratio — kernel + pySecDec work for the
   orchestrator, ~minutes replay with the cached pentagon package
   (`Tests/refs/generators/README.md`).  Schedule it BEFORE any M5
   pentagon judgment (it is independent of, and parallel to, the
   section-2 patch).  Same protocol run for box_triangle (also
   PREDICTED ONLY).

---

## 4. M5 triage-taxonomy classification of the pentagon failure

Per the M5 gate taxonomy (RewritePlan section 6: D2-class /
prescription-config-class / pin-class):

### 4.1 Warning inventory of `/tmp/f_pentagon.log`

| class | count | log lines | source | content impact |
|---|---|---|---|---|
| D2: zero-regulator drop ("resonant endpoint coefficient with zero epsilon regulator; dropping contribution") | 20 | 44-63 (level-2 boundary integration) | `DiffExp/RegularizedIntegration.m:1467-1477` (drop at `:1476` when `b ~= 0` AND `a+n+1 ~= 0`) | SILENT loss with NO order bookkeeping — can corrupt any eps order of the level-1 BCs, poles included |
| D2: omitted endpoint coefficients, in-window ("omitted 19 ... starting at epsilon order 0 ... NOT trustworthy") | 2 (19 coeffs each) | 74-75 (level-1 -> final integration) | `DiffExp/RegularizedIntegration.m:1668-1679` | honestly disclaimed: final orders >= eps^0 untrusted; does NOT cover the poles |
| D2: omitted endpoint coefficient, near-window | 1 (1 coeff, >= eps^17) | 64 (level 2) | same | negligible at the needed orders |
| prescription-config: unrecognized branch point | 2 | 72-73 (level-1 transport, one per leg) | `DiffExp/Transport.m:472` via `:485-487` | O(1) relative branch error at ALL eps orders of all level-1 data crossing `(19 +- Sqrt[105])/32`; injects the spurious `Im(eps^0) = -0.319` |
| cosmetic: `General::munfl` underflow | 3 + suppression | 31-39 (level-3 transport) | WL numerics (~1e-308 magnitudes) | none |

Calibration: the five VALIDATED ladder runs show ZERO occurrences of
either D2 warning (checked: all box_bubble logs incl.
`Tests/refs/oracle_dumps/m0_oracle_boxbubble.log`).  These warnings are
not benign background; pentagon (and box_triangle) are the only emitters.

### 4.2 Attribution

- The leading-pole error (`eps^-2`: 0.25214 vs +13/42, ~18.5% low) is
  attributable to BOTH live classes: the 20 level-2 silent drops (order
  bookkeeping absent by construction) and the level-1 branch
  misconfiguration (order-blind by nature).  The omission warnings
  (lines 74-75) explicitly do NOT cover the poles.
- The complex finite part is attributable to the prescription-config
  class specifically (section 1.5).
- Pin-class is EXCLUDED as a cause (section 3, finding 2).

### 4.3 What the rewrite fixes vs what persists

- D2 share — FIXED BY CONSTRUCTION in DiffExp2: exact
  `x^(a + b eps) (eps Logx)^p` tags mean there are no collapsed towers to
  resum and no zero-regulator denominators to drop
  (`RegularizedIntegration.m`'s drop/salvage rules have no DiffExp2
  counterpart; RewritePlan sections 2/3.3, Integrate.m case table).  If a
  DiffExp2 pentagon run still fails at the leading pole AFTER the
  section-2 patch, that is an M5 BLOCKER-class finding.
- Prescription-config share — PERSISTS in DiffExp2 if `FTExamples.m` is
  not patched: prescriptions are irreducibly per-example input data.
  Behavior difference: DiffExp2's contract escalates the silent-warn to
  a LOUD ERROR (RewritePlan section 3.1 "Prescriptions"), so an
  unpatched pentagon will refuse to produce a number rather than produce
  a wrong one.  SPEC GAP flagged while verifying this: the plan's
  loud-error trigger is worded "missing prescription at a chart with
  b!=0/p>0 sectors", but the pentagon charts at the Kallen roots have
  HALF-INTEGER `a` with `b = 0, p = 0` — multivalued via fractional `a`.
  The Prescriptions error contract must also trigger on
  `Denominator[a] > 1` (the exact analogue of old `Transport.m:483`
  `AlgebraicRootsPresent`).  Spec agents for SectorSeries/Transport
  (M0 tasks 4-11) must fold this in, else DiffExp2 reproduces this very
  silent failure.

Bottom-line split: of the THREE observable defects (wrong poles, untrusted
finite window, spurious imaginary part) — the untrusted-window defect is
pure D2; the imaginary part is pure prescription-config; the wrong poles
are jointly caused, with the exact split only measurable by the
section-5 step-3 diagnostic rerun.

---

## 5. Recommended M5 retest order (pentagon)

1. ZERO-KERNEL (now): land the section-2 patch on a branch (it is inert
   for all currently-validated examples: the new function returns `{}`
   for every `{example, level}` except `{"pentagon", 1}`).
2. IDLE-KERNEL, INDEPENDENT (before M5 judgment; can run during M0
   oracle generation): the R9 protocol run for pentagon AND box_triangle
   (section 3, finding 3).  Promote both pins from PROVISIONAL/PREDICTED
   in `Tests/PINS.md`.  Cheap: pySecDec packages are cached, replays are
   minutes.
3. OLD-CORE DIAGNOSTIC RERUN (recommended, ~1 pentagon pass): rerun
   `run_ft_stepwise.m` pentagon WITH the patch on the frozen old core.
   This does not violate R8 (no old-core code is debugged; only
   example/config data changes in Scripts/).  Expected reading:
   warnings 72-73 gone; `Im(eps^0) -> ~0`; `eps^-2` moves toward +13/42
   but plausibly still off while the 20 level-2 drops remain.  The
   residual error IS the measured D2 share — the calibration number the
   M5 gate needs to judge the DiffExp2 pentagon run, and the empirical
   confirmation (or refutation) of the `+1` prescription sign per
   section 2.3.  Budget note: also raise the epsilon lookahead
   (`FT_POLE_ALLOWANCE`, default 4, `run_ft_stepwise.m:36-37`) if the
   omitted-19 warnings persist — they said "rerun with more epsilon
   lookahead".
4. M5 LADDER POSITION (unchanged from the plan: ...box_bubble,
   box_triangle, pentagon...): run DiffExp2 pentagon only after (a) the
   budget-formula unit reproduces box_bubble's 9/11 (RewritePlan
   section 3.4 validation), since pentagon's level-1 `EpsPrefactor = 2`
   /`RawMinPower = -2` tower plus the level-2 needs is exactly the
   cumulative-budget shape that starved box_bubble pre-c0b24f3, and
   (b) steps 1-2 are done.
5. M5 VERDICT TABLE for pentagon (decision rule at the gate):
   - passes all orders -> close R8 item; promote pin to VALIDATED.
   - leading pole still wrong, warnings clean -> D2-residual/BLOCKER
     (rewrite was supposed to fix it; do not reclassify).
   - only `Im(eps^0)` nonzero -> prescription SIGN: flip per section 2.3.
   - all orders off by one constant factor -> pin-class; re-run R9
     protocol before touching any code.
   - DiffExp2 loud-errors "missing prescription" naming a chart NOT at
     the Kallen roots -> a further pentagon-level factor was missed at
     levels 2-4 (none is predicted by section 1.3's reconstruction, but
     levels 2-4 alphabets were only checked against IBP factors + chart
     centers, not exhaustively): extend
     `FTExampleLevelDeltaFactors` with the named factor, rerun.

---

## 6. Collateral findings (for other M0 agents; none blocks this triage)

(a) VENDOR THE DUMPS: `/tmp/diffexp_pentagon_dumps/laurent_integral_0006.m`
    is the only artifact proving the chart centers (sections 1.2) and
    /tmp is not storage (power-loss lesson).  Copy all six dump files to
    `Tests/refs/oracle_dumps/pentagon_laurent/` alongside the already
    vendored `f_pentagon.log` (pins agent, M0 task 14 scope).
(b) SPEC GAP (math/minimalism spec reviewers, M0 tasks 12-13): the
    DiffExp2 Prescriptions loud-error must trigger on fractional-`a`
    charts too — see section 4.3.
(c) HYGIENE (FT hardening batch, A1/M6): the dedup in
    `deltaPrescriptionsForFactors` (`DiffExpIntegration.m:234-237`)
    identifies factors up to OVERALL SIGN, so `{-1+xx1, +1}` silently
    becomes `{1-xx1, +1}` — which denotes the OPPOSITE iδ side for the
    same zero.  Benign for the validated examples (endpoint factors) but
    it is exactly the class of silent convention flip A1 exists to kill;
    make the dedup sign-aware or normalize factors to a canonical sign
    BEFORE attaching the prescription sign.
(d) PREMISE CORRECTION for RewritePlan section 2/6 wording: the pentagon
    DeltaPrescriptions gap is in the AUTO-GENERATION path
    (`FeynmanTrick/DiffExpIntegration.m`) + missing example data, not in
    an existing `FTExamples.m` prescription list (none exists).  The M5
    gate text "prescription/config-class (FTExamples fix)" remains
    accurate once the section-2 patch creates that surface.
