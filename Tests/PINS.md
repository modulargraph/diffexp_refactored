# Validation pin suite (M0 deliverable — execution contract)

Status: consolidated 2026-06-11 (M0 task 14 of Docs/RewritePlan.md).  Every
pin the DiffExp2 rewrite is judged against, with value, provenance, and
precision.  All reference data is vendored IN-REPO under `Tests/refs/`
(power-loss lesson: /tmp is not storage).  See `Tests/refs/MANIFEST.md` for
the file-by-file inventory and `Tests/refs/generators/README.md` for how to
regenerate any pin from scratch.

Conventions used throughout:

- Kinematics: box family at s = -1, t = -1/3; pentagon at
  {s12,s23,s34,s45,s15} = {-1,-2,-3,-5,-7}; d = 4-2eps; Feynman-trick
  parameter fixed at 11/23 for boundary-level (per-master) references.
  Examples defined in `Scripts/FTExamples.m`.
- "FT value" = the FT-family normalization produced by the pipeline
  (`Scripts/run_ft_stepwise.m`):
  prefactor = gamma(v - L + L*eps)/prod_i gamma(v_i),
  U power = v-(L+1)+(L+1)eps, F power = -v+L-L*eps
  (Docs/FeynmanTrickBananaStatus.md:27-33).  STEPWISE log rows print Laurent
  coefficient pairs [k, c_k] of the full FT-normalized master (prefactor
  included); they are DIRECTLY comparable to family-route pySecDec values
  and to the FT-side of the sign rule below.
- "loop_package value" = pySecDec `loop_package` with Minkowski propagators
  (l+...)^2.  The FT family uses Euclidean propagators -(...)^2, so:

      FT = (-1)^(total propagator power) x loop_package
         = (-1)^nprops for corner integrals (all powers 1).

- Comparator: `Scripts/compare_stepwise_log.py <log> [level] [--refs f.json]`,
  relative tolerance TOL = 5e-6 (compare_stepwise_log.py:31).  Banana
  references are hardcoded there (lines 18-29); driver `results.json` files
  load via `--refs` (caveat: entries whose source_name ends `_needed` are
  skipped at compare_stepwise_log.py:65-66 — see Tests/refs/MANIFEST.md).
- Quoted "+/-" errors are pySecDec Cuhre error estimates; at low maxeval they
  can be optimistic by up to ~10x (measured on box_bubble, section 2.2).

---

## 1. The sign-rule pin (R9 protocol)

RULE: FT = (-1)^(sum of propagator powers) x loop_package.

| example          | nprops | factor | status |
|------------------|--------|--------|--------|
| box (1L)         | 4      | +1     | VALIDATED: FT box values verified analytically in the test suite AND end-to-end vs the box1l pin (section 2.1) |
| box_bubble       | 5      | -1     | VALIDATED end-to-end 2026-06-11 at all 4 orders (section 2.2; Docs/FeynmanTrickBoxFamilyStatus.md:19-26) |
| pentagon         | 5      | -1     | PREDICTED ONLY — pentagon pipeline has never passed; pin held PROVISIONAL (R9) |
| box_triangle     | 6      | +1     | PREDICTED ONLY — pipeline fails for unrelated (D2) reasons; re-run protocol before M5 judgment (R9) |
| double_box_planar| 7      | -1     | PREDICTED ONLY — never integrated |

Protocol (what "re-run" means): generate the loop_package value with the
vendored generator, generate the FT-family value for the same integral with
`Scripts/export_pysecdec_family_specs.m` + `Scripts/pysecdec_family_driver.py`
(exact FT normalization, no sign ambiguity), and confirm the ratio.  That is
exactly how the box (+1) case was pinned
(`Tests/refs/generators/pysecdec_box1l_pin.py` docstring).  Mixing the two
routes without the factor produces phantom mismatches.

---

## 2. End-to-end example pins (the M5 ladder)

### 2.1 box (1-loop, 4 props, sign +1) — VALIDATED

Pin (loop_package, Cuhre epsrel 1e-7, epsabs 1e-12, maxeval 5e6; raw record
`Tests/refs/pysecdec/box1l_pin_result.log`):

| order  | pin value            | Cuhre err | old-core record (Tests/refs/oracle_logs/l2_box.log FINAL) | agreement |
|--------|----------------------|-----------|------------------------------------------------------------|-----------|
| eps^-2 | 12 (exact)           | 8.1e-15   | 12.000000000000000000...                                   | exact     |
| eps^-1 | -0.3349142468097521  | 1.9e-13   | -0.33491424680973617928...                                 | 1.6e-14   |
| eps^0  | -41.28416739757452   | 4.7e-10   | -41.28416739754646857143...                                | 2.9e-11 (11 digits, inside ref error) |
| eps^1  | -70.1654426808142    | 8.3e-07   | (not produced by the recorded run — printed through eps^0) | UNTESTED vs old core |

Generator: `Tests/refs/generators/pysecdec_box1l_pin.py`.
Doc record: Docs/FeynmanTrickBoxFamilyStatus.md:15-16,77-79.

### 2.2 box_bubble (2-loop 3-point, 5 props, sign -1) — VALIDATED

Pin: FT = -1 x loop_package.  loop_package values (Cuhre, defaults of
`pysecdec_ft_boxes.py`: epsrel 1e-3, epsabs 1e-7, maxeval 2e5; raw record
`Tests/refs/pysecdec/loop_package_box_results.jsonl`):

| order  | loop_package value    | Cuhre err | FT pin (= -1 x) | FT validated value (Tests/refs/oracle_dumps/m0_oracle_boxbubble.log) | agreement |
|--------|-----------------------|-----------|------------------|------------------------------------------------------------------------|-----------|
| eps^-3 | -0.5 (exact -1/2)     | 3.9e-16   | +0.5             | +0.49999999999999999999...                                             | 1.2e-20   |
| eps^-2 | -0.42278433509626767  | 5.1e-12   | +0.422784335096  | +0.42278433509846713945...                                             | 2.2e-12   |
| eps^-1 | -0.35626900321018007  | 3.0e-06   | +0.356269003210  | +0.35627956058053973479...                                             | 1.1e-05   |
| eps^0  | +4.877038708178345    | 9.6e-06   | -4.877038708178  | -4.87713966245589825427...                                             | 1.0e-04   |

PRECISION HONESTY: the eps^-1/eps^0 agreements exceed the quoted Cuhre
errors by ~3-10x; treat this pin as good to ~1e-5 relative at eps^-1 and
~2e-5 relative at eps^0 until regenerated at higher maxeval (command in
Tests/refs/generators/README.md).  The campaign verdict "validated at the
reference's accuracy" is Docs/FeynmanTrickBoxFamilyStatus.md:19-26.

Budget pins attached to this example (must be reproduced by the plan §3.4
formula, M5 unit): the per-level transports need FULL boundary depth —
eps order 9 at L2 and 11 at L1.  Failure record with the starved budget:
`Tests/refs/oracle_logs/ft_box_bubble_run2_budget_starved.log` lines 29,32
("WARNING level 2 transport capped at eps order 2 (needs 9)", "level 1 ...
(needs 11)").  Fix commit c0b24f3.

### 2.3 pentagon (1-loop 5-point, 5 props, sign -1) — PROVISIONAL (R9)

Pin (loop_package, Cuhre epsrel 1e-8, epsabs 1e-12, maxeval 3e6; raw record
`Tests/refs/pysecdec/pentagon_pin_result.txt`):

| order  | loop_package value       | Cuhre err | FT pin (= -1 x, PROVISIONAL) |
|--------|--------------------------|-----------|-------------------------------|
| eps^-2 | -0.309523809523809978    | 2.6e-16   | +0.309523809524               |
| eps^-1 | +0.0164223473245557472   | 5.4e-14   | -0.0164223473                 |
| eps^0  | +0.0254117766957063611   | 1.5e-09   | -0.0254117770                 |
| eps^1  | -0.279721093096727325    | 3.6e-09   | +0.2797210931                 |

Closed-form anchor: eps^-2 = -13/42 EXACTLY (-0.30952380952380952...),
matched to 4e-16 — use it as the leading-pole unit regardless of sign
protocol.  Generator: `Tests/refs/generators/pysecdec_pentagon.py`
(dot products fixed numerically to match Scripts/FTExamples.m).

PROVISIONAL because: (a) the (-1)^5 factor is predicted, not pinned, for
this example (section 1); (b) the FT pipeline currently FAILS at every
order including the leading pole (FAIL baseline
`Tests/refs/oracle_logs/f_pentagon.log`: FT eps^-2 = 0.25213984... vs
expected +0.30952...; plus the unresolved DeltaPrescriptions configuration
warning — M0 task 16 triages).  Do not pass M5 judgment on pentagon before
the R9 protocol run.

### 2.4 box_triangle (2-loop 4-point, 6 props, sign +1) — pin trusted, FT FAILS

Pin (loop_package, same defaults as box_bubble; raw record
`Tests/refs/pysecdec/loop_package_box_results.jsonl`).  Sign +1 (6 props),
so FT pin = loop_package values directly:

| order  | pin value           | Cuhre err |
|--------|---------------------|-----------|
| eps^-4 | 2.25 (exact 9/4)    | 1.9e-15   |
| eps^-3 | 2.346284806970911   | 9.9e-11   |
| eps^-2 | 4.924311554761632   | 1.0e-05   |
| eps^-1 | 11.944387049203305  | 7.4e-04   |
| eps^0  | 9.470615214142624   | 4.2e-03   |

PRECISION HONESTY: trust eps^-1 to ~6e-5 relative and eps^0 to ~4e-4
relative (quoted Cuhre errors at maxeval 2e5).  Sign re-check before M5
judgments per R9 (plan §9).  Doc record:
Docs/FeynmanTrickBoxFamilyStatus.md:17-18.

FAIL baseline: `Tests/refs/oracle_logs/f_box_triangle.log` (FT eps^-4 =
153.3648... vs pin 2.25; complex contamination from eps^-3) — confirmed
D2-class (Docs/RewritePlan.md §10).  Geometric trigger: interior IBP pole
at xx3 = 12167/12651 ~ 0.9617 near the L3 line's endpoint.

### 2.5 banana (equal-mass 2-loop, FT comparator suite) — VALIDATED

References from the FAMILY-EXPORTER route (exact FT normalization, z0 =
11/23) — hardcoded in `Scripts/compare_stepwise_log.py:18-29`, mirrored
here.  Exact rationals are noted; treat the rest as good to the comparator
TOL 5e-6 relative.

| level, master      | reference Laurent coefficients |
|--------------------|--------------------------------|
| L1 {2,0,1,1}       | eps^0: 5.4025802965 |
| L1 {1,0,1,0}       | eps^-2: -2.0037878788, eps^-1: 2.2439088015, eps^0: -13.5042252469 |
| L1 {1,0,0,1}       | identical to {1,0,1,0} (symmetry — use as cross-check) |
| L1 {1,0,1,1}       | eps^-3: 1/3 (exact), eps^-2: -0.1144868285, eps^-1: 1.3904392399, eps^0: -5.8738356912 |
| L1 {1,0,1,2}       | eps^-2: 1/2 (exact), eps^-1: -0.1717302427, eps^0: -0.9025092631 |
| L1 {1,0,2,1}       | identical to {1,0,1,2} (symmetry) |
| L1 {1,-1,1,1}      | NO direct reference (exporter cannot encode numerator powers); validated indirectly through L0 |
| L0 {1,1,1,1}       | eps^0: 8.26810451329511583109184 (no eps poles) |

L0 provenance: direct numerical Feynman-parameter integration, ~24 digits
(Docs/FeynmanTrickBananaStatus.md:215-218); independently corroborated by
pySecDec's direct L0 run 8.2681 +/- 3e-5 (same doc:336-337).  The L0
generator script (`/tmp/banana_direct_feynman_param_c0.m`) is LOST — see
section 6.

PASS record: `Tests/refs/oracle_logs/l2_banana.log` — comparator failures: 0;
final L0 value 8.2681045358689687... agrees with the reference to 2.7e-9
relative at FT_WORKING_PRECISION=300 (the doc's 2.3e-9 record used the
canonical 500/50/4/8 settings; both limited by accuracy decay at the x=1/2
resonant crossing, Docs/FeynmanTrickBananaStatus.md:482-487).  Deeper-level
boundary values (L3/L2) are recorded as STEPWISE rows in the same log and in
Docs/FeynmanTrickBananaStatus.md:52-57.

Second-point checks at z0 = 2/5 (anti-overfitting pins; underlying
results.json LOST, values preserved at
Docs/FeynmanTrickBananaStatus.md:300-302,330,466-469):
L1 {1,0,0,1} -> -0.694444/eps^2 (leading); {1,-1,1,1} -> 0.166667/eps^3;
{1,0,1,1} eps^-3 -> 1/6 (vs 1/3 at 11/23 — kills the frozen-coefficient
failure mode).  Regenerate via the vendored spec pattern with
FT_FIXED_VALUE=2/5.

### 2.6 bubble and sunrise — VALIDATED (independent refs LOST, values preserved)

Campaign-validated end-to-end finals (Docs/FeynmanTrickBananaStatus.md:478-481):

    bubble  {1,1}   eps^0 = 0.860817881928008
    sunrise {1,1,1} eps^0 = 2.2367927002126465   (no eps poles in either)

Full-precision (50-digit) old-core STEPWISE records, including all
intermediate boundary levels: `Tests/refs/oracle_logs/l2_bubsun.log`.
NOTE these logs are OLD-CORE OUTPUT (parity oracles); the independent
pySecDec results.json files for bubble/sunrise are LOST (section 6).  The
two values above were verified against pySecDec while those files existed;
regeneration commands are in Tests/refs/generators/README.md.

### 2.7 double_box_planar — NO PIN

Never completed (FIRE ~40 min/iteration; needs ReductionCache, plan A3).
Generator support exists (`pysecdec_ft_boxes.py double_box_planar`,
expected sign -1 with 7 props) but no integration was ever run to pin
precision.  M5 stretch goal (R5).

---

## 3. Closed-form pins (exact; no numerics needed)

These are the highest-value pins: check them symbolically to any precision.

1. BOX L2 ENDPOINT TOWER (the campaign's flagship exactness record):
   the L2-line upper-endpoint limit tower equals

       Gamma[1+eps] * Beta[1-eps,1-eps] * (4/23)^(-eps)

   verified EXACT TO 20 DIGITS at all computed eps orders
   (Docs/FeynmanTrickBoxFamilyStatus.md:125-128).  Identification: this is
   eps * M2(t) at t = 11/23, with the L1 analytic structure
   M1(t) = Gamma[eps] B(1-eps,1-eps) t^(-eps),
   M2(t) = Gamma[eps] B(1-eps,1-eps) ((1-t)/3)^(-eps)
   (Docs/FeynmanTrickBoxFamilyStatus.md:43-46).  M2/M3 unit per plan.

2. BOX L1 BOUNDARY TOWERS: the eps-prefactor-normalized L1 boundary rows
   must equal Gamma[1+eps] B(1-eps,1-eps) (11/23)^(-eps) (master {1,0,1,0})
   and Gamma[1+eps] B(1-eps,1-eps) (4/23)^(-eps) (master {1,0,0,1}).
   50-digit recorded values (Tests/refs/oracle_logs/l2_box.log STEPWISE):
   {1,0,1,0}: eps^-1: 1, eps^0: 2.16038327822924628607525340456874...
   {1,0,0,1}: eps^-1: 1, eps^0: 3.17198418990772621136577648878322...
   These were verified EXACT against independent sympy/mpmath quadrature
   (Docs/FeynmanTrickBoxFamilyStatus.md:37-39;
   `Tests/refs/generators/crosschecks/eval_box_specs.py`).

3. PENTAGON LEADING POLE: loop_package eps^-2 = -13/42 exactly
   (FT convention: +13/42 under the predicted -1 factor).

4. BOX_TRIANGLE LEADING POLE: eps^-4 = 9/4 exactly.

5. BOX_BUBBLE LEADING POLE: loop_package eps^-3 = -1/2 exactly
   (FT: +1/2 — reproduced to 1.2e-20 in the validated run).

6. BOX L0 LEADING POLE: eps^-2 = 12 exactly.

7. BANANA EXACT RATIONALS: L1 {1,0,1,1} eps^-3 = 1/3 at z0=11/23 (and 1/6
   at z0=2/5); L1 {1,0,1,2}/{1,0,2,1} eps^-2 = 1/2.

---

## 4. Structural / pointwise pins (M2-M4 units)

INDICIAL PINS:
- Box L2 apparent chart: the L2 transport line crosses an APPARENT
  singularity at t* = 7/11, the root of the quadratic Symanzik coefficient
  beta(t) = -A*B - A/3 + 4B/3 with A = 1 - 11t/23, B = 1 - t
  (Docs/FeynmanTrickBoxFamilyStatus.md:80-88); recorded chart eigenvalue -1
  (plan M2, Docs/RewritePlan.md:364-367).  The function is analytic there:
  all genuine local x^-1/x^-2 content must vanish by analyticity.
  NOTE: the campaign's 20-digit pointwise J1/J2 values at this chart were
  NOT found in any on-disk artifact (section 6) — regenerate at M2.
- Banana L1 endpoints xx1 = 0, 1: NILPOTENT DOUBLE POLE in the numerator-
  master row + 5-fold resonant residues; first and worst in-scope
  higher-order-pole case (Docs/FeynmanTrickBananaStatus.md:222-227,282-285;
  plan §3.2 Indicial, R11).
- Banana L2 upper endpoint sector mix: x^0, x^(-1+eps), x^(2 eps) at one
  integer power — the canonical pseudo-resonance unit
  (Docs/FeynmanTrickBananaStatus.md:343-346; plan I2).
- Box L1 line: IBP-coefficient poles at t = 0, 1/4, 1 — the 1/4 pole is
  INTERIOR to the integration range (first example ever;
  Docs/FeynmanTrickBoxFamilyStatus.md:45-55, commit 9aeb300, tests in
  Tests/test_interior_singular_integration.m).
- Box_triangle L3 line: interior IBP pole at xx3 = 12167/12651 ~ 0.9617
  (exact rational from campaign memory; plan quotes 0.9617 at
  Docs/RewritePlan.md:383,461; M4 unit per plan).

POINTWISE TRANSPORT PINS (box L1 needed-integrand {2,0,1,1}(t), family
route, straddling the t = 1/4 interior pole; vendored
`Tests/refs/pysecdec/results_box_t*.json`):

| t      | eps^-1 (exact-grade)      | eps^0 (noisy)      |
|--------|---------------------------|---------------------|
| 1/20   | -63.1579 (= -1200/19)     | 1357.88 +/- 141.7   |
| 2/5    | -12.5 (= -25/2)           | 280.223 +/- 16.2    |
| 19/20  | -63.1579 (= -1200/19)     | 1190.83 +/- 59.1    |
| 11/23  | -12.0227                  | 311.889 +/- 12.6    |

(11/23 row from `results_box1l.json` entry box_L1_2x0x1x1_needed.)  Use the
eps^-1 column as a hard pin; eps^0 only as a sanity band.  The t=1/20 vs
t=19/20 equality of eps^-1 is an IBP-coefficient symmetry — another
anti-frozen-coefficient check.

BUDGET PINS (plan §3.4 validation unit): box_bubble needs eps depth 9 at L2
and 11 at L1 (section 2.2 above).

---

## 5. Comparator protocol

    # banana / bubble / sunrise (banana refs hardcoded):
    python3 Scripts/compare_stepwise_log.py <stepwise.log>

    # box family vs vendored driver refs:
    python3 Scripts/compare_stepwise_log.py <stepwise.log> \
        --refs Tests/refs/pysecdec/results_box1l.json

    # exit code 0 = all matched within TOL 5e-6 relative; prints per-order
    # PASS/FAIL and flags reference orders missing from the run.

Caveats:
- `--refs` skips entries whose source_name ends in `_needed`
  (compare_stepwise_log.py:65-66), so the t-variant files and the _needed
  rows of results_box1l.json do not auto-compare; read their `summary`
  fields directly (format parsed by the same script's TERM_RE).
- results_box1l.json's box_L0_1x1x1x1 entry is noisy garbage (manifest);
  a box run compared `--refs results_box1l.json` therefore reports ONE
  EXPECTED FAIL at L0 eps^-1 (ref -175.59).  Judge box L0 against the
  box1l pin (section 2.1), not against that entry.  Verified on the
  vendored l2_box.log: 1 failure, exactly this line.

The loop_package pins (sections 2.1-2.4) are FINAL-value pins, not STEPWISE
comparator refs: compare the run's FINAL/L0 row against the table after
applying the sign rule.

---

## 6. NOT FOUND — regeneration required (orchestrator queue)

Searched /tmp, Scripts/, Tests/ on 2026-06-11.  Lost to /tmp cleanup
(values survive in docs; raw reference files do not):

1. `/tmp/pysecdec-ft-family-bubble/results.json`,
   `/tmp/pysecdec-ft-family-sunrise/results.json`,
   `/tmp/pysecdec-ft-family-banana-deep/results.json`,
   `/tmp/pysecdec-ft-family-banana-level1-pos-factored/results.json`
   (the original family-route reference set,
   Docs/FeynmanTrickBananaStatus.md:17-24).  Regenerate: export specs with
   FT_EXAMPLE=bubble/sunrise/banana + run the family driver
   (Tests/refs/generators/README.md; banana specs already vendored as
   Tests/refs/pysecdec/specs_banana.json).
2. `/tmp/banana_direct_feynman_param_c0.m` — generator of the 24-digit
   banana L0 reference 8.26810451329511583109184.  Value preserved
   (compare_stepwise_log.py:28, banana doc:217); generator lost.
3. Banana second-point (z0 = 2/5) results.json — values preserved at
   Docs/FeynmanTrickBananaStatus.md:300-302,330,466-469 only.
4. Box L2 apparent-chart (t* = 7/11) 20-digit pointwise J1/J2 values —
   exist only in campaign session memory; regenerate from the stored
   campaign matrices at M2 (plan M2 gate names this pin).
5. `/tmp/diffexp_banana_l1_dumps_v2final/` — known-good banana L1 Laurent
   dumps (mapping documented at Docs/FeynmanTrickBananaStatus.md:204-209).
   Regenerate with DIFFEXP_DUMP_LAURENT_DIR during a banana run if dump
   replay is wanted for banana (bubble/sunrise dumps ARE vendored:
   Tests/refs/oracle_dumps/laurent_bubsun.tar.gz).
6. Box L0 eps^1 old-core record — the vendored l2_box.log prints through
   eps^0 only; the eps^1 pin value (-70.1654426808142) is from pySecDec and
   has never been reproduced by the old core.
7. `/tmp/l2_boxbubble.log` is truncated (no FINAL row) and was NOT vendored;
   the complete box_bubble PASS record is
   Tests/refs/oracle_dumps/m0_oracle_boxbubble.log (byte-identical to the
   campaign's ft_box_bubble_run3.log — deterministic reproduction).
