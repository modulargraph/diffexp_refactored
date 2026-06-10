# FeynmanTrick box-family campaign status (d = 4-2eps, Euclidean)

Updated: 2026-06-10.  Predecessor: `FeynmanTrickBananaStatus.md` (banana
RESOLVED).  This documents the massless box-family examples defined in
`Scripts/FTExamples.m`: `box`, `box_bubble`, `box_triangle`,
`double_box_planar` (s = -1, t = -1/3) and `pentagon`
({s12,s23,s34,s45,s15} = {-1,-2,-3,-5,-7}).

## Reference values (pySecDec)

loop_package route (`/tmp/pysecdec_ft_boxes.py`, Minkowski propagators):
relative factor (-1)^(number of propagators) vs the FT family
normalization - confirmed +1 for the 1-loop box by the pin below.

- box 1L (pin, epsrel 1e-7):
  12/eps^2 - 0.334914246810/eps - 41.28416739757452 - 70.1654426808142 eps
- box_triangle: 2.25/eps^4 + 2.346284806971/eps^3 + 4.924311554762/eps^2
  + 11.944387049203/eps + 9.470615214143   (sign flip NOT needed: 6 props)
- box_bubble: -0.5/eps^3 - 0.422784335096/eps^2 - 0.356269003210/eps
  + 4.877038708178   (FT convention = MINUS these: 5 props)

Family-exporter route (exact FT normalization, any level):
`Scripts/export_pysecdec_family_specs.m` (FT_EXAMPLE, FT_FIXED_VALUE) +
`Scripts/pysecdec_family_driver.py`.  The exporter now derives U/F
exponents from the configured dimension.

## Box bisection: where the FT recursion stands

FT box stepwise (`FT_EXAMPLES=box FT_EPS_ORDER=1`):

| quantity | status |
|---|---|
| L3, L2, L1 master boundary values | EXACT (match independent sympy/mpmath evaluation of the exported specs to all digits; `/tmp/eval_box_specs.py` pattern) |
| L0 eps^-2, eps^-1 | EXACT vs pin (12 and -0.334914246810, 1e-13) |
| L0 eps^0 | -47.9756 vs pin -41.2842 (deficit +6.69) |

Analytic level-1 structure (from the t-dependence of exported F
polynomials): M1 = Gamma(eps) B(1-eps,1-eps) t^(-eps),
M2 = Gamma(eps) B(1-eps,1-eps) ((1-t)/3)^(-eps).  The L1->L0 integrand
{2,0,1,1}(t) carries IBP-coefficient poles at t = 0, 1/4, 1.

## Fixed this campaign (commit 9aeb300)

1. Complex-log leakage in segments straddling an interior singular
   point (the IBP pole at xx1 = 1/4 sits INSIDE the integration range -
   first example ever; banana-era factors had no roots in (0,1)).
2. Unit-regulator over-regularization of regular powers of meromorphic
   terms (artificial eps-mixing; exact Hadamard finite-part now).
Tests: `Tests/test_interior_singular_integration.m` (11 checks).

## RESOLVED: the endpoint-series corruption (NormalizeLogPower export)

The endpoint-segment series corruption is FOUND and FIXED:
`DiffExp`SeriesOps`NormalizeLogPower` had no usage declaration, so
Integration.m's fully qualified call hit an undefined public symbol and
returned unevaluated; the integer filter in DiffExpIntegrate's entry then
saw no powers and the IMaxLogOrder auto-extension NEVER fired.  Any
solve relying on DiffExpIntegrate with sources deeper than the current
IntReps (the box L1 endpoint blocks: diagonal A = eps diag(-1/x, ...),
per-order A0 = 0, dispatched to SolveSimple = pure integration; sources
gain one Logx power per epsilon order) silently integrated
x^-1 Logx^k (k >= 2) terms AS CONSTANTS (the b*Const fallback).
Breakage exactly from epsilon order 3 (= IntReps depth 1 + 2), matching
the observed ODE residuals {3,4,5,6} and the spurious (r-1)^3 confluent
moments.  Post-fix: endpoint ODE residuals = 0 exactly, towers are clean
single-sector, the salvage warnings are gone, and
Tests/test_integration_log_depth.m locks the behavior.

## Remaining: eps^0 deficit (+6.69) now at the COMBINATION level

Post-solver-fix state: the box L0 totals are bit-identical
({12, -0.334914, -47.9756, -64.1316, -20.1632, ...} all real) because
the eps^0 output only consumes tower offsets <= 2, which were always
correct.  Verified piece by piece: interior segments match NIntegrate;
straddling segment 5 matches its closed form; endpoint segment 1's
machinery integration matches an independent termwise closed-form
reconstruction to 13 digits; the b != 0 endpoint monomial machinery is
exact on synthetic input; the endpoint towers satisfy the level ODE
exactly.  The pin (-41.28417 at eps^0) is confirmed independently by
the analytic all-orders box formula.

RESOLVED FURTHER (forensics complete, fix pending): the dump integrand
IS the true G = c1 M1 + c2 M2 (FIRE coefficients verified: c1 =
(18-6d)/(t-4t^2), c2 = 18(d-3)/(1-5t+4t^2); pointwise eps^0 values
match the dump analytically at t = 0.05, 0.4, 0.95 to all digits; the
exporter's "needed"-spec is a DIFFERENT normalization - same class as
the banana "factor 2" memory note - red herring).  Per-segment
comparison against semi-analytic window integrals (partial fractions +
Beta forms, validated: window sum = pin to Laurent truncation) pins the
ENTIRE deficit on SEGMENT 12 (the xx1 = 1 endpoint segment):

  seg12 true Laurent:      {6, 12.33311519, +4.892941862, -24.49, ...}
  seg12 machinery Laurent: {6, 12.33311519, -1.798476481, -66.23, ...}
  Delta(eps^0) = -6.691418 = exactly the global deficit; eps^-2/-1
  agree to 10 digits; segment 1 (xx1 = 0) is CORRECT (its termwise
  closed-form reconstruction matches the machinery to 13 digits, and
  matches the true window integral).

The machinery's seg12 value equals the termwise closed-form of its own
series under the (+i pi power-phase, +i pi log-value, formal drop at
the singular end) convention - i.e. integration is internally
consistent; the stored seg12 tower carries -i pi w1-relative content
where seg1 carries +i pi w1 (the toUpper transport's branch data for
the (-1+xx1) factor vs the local-side theta resolution).  Pointwise
evaluation of the same data is correct, so the inconsistency is
specifically between the upper-line stored branch structure and the
endpoint-integration convention.  Flipping evaluation-side branches
alone is inconsistent (breaks eps^-1); the fix must align the
convention pair (stored data <-> integration) for upper-anchored
segments.

Probe inventory (all /tmp): per_segment_truth.py (semi-analytic
windows), probe_seg12_exact.m (+ _conj/_split variants),
probe_seg1_exact.m, identity_finite_eps.py, probe_fire_reduction.m.

Previously suspected and now exonerated:
Pointwise (machinery convention: local-side thetas + complex local
log) the dump integrand is real and smooth:
  t = 1/20 : -126.3158/eps - 349.2138
  t = 2/5  : -25/eps      - 43.1342
  t = 19/20: -126.3158/eps - 453.3886
with dump(eps^-1) = -6/(t(1-t)) = 2 x spec(eps^-1) at all three points
(spec eps^-1 = -3/(t(1-t)) exactly, verified both by pySecDec and the
z0-endpoint residue formula 1/(A(A+B))).  A clean global factor 2 with
1x totals at the poles is self-contradictory, so something in this
reading is convention-skewed; the four-way real-eps test
(/tmp/fourway_eps_test.py: 3d quad of the box, int dt of the spec, pin
Laurent, FT Laurent, all at eps = -1/20) discriminates which leg is
broken without any Laurent bookkeeping.

- branch confluent-sectors-wip holds an N-root/confluent Prony
  generalization of FitResidualEndpointSectors built while chasing the
  (now-explained) spurious confluent moments; with the solver fixed the
  towers are plain and the existing 2-root fitter suffices for the box.
  Keep the branch for genuinely multi-sector future families.

## Practical guidance

- POLES of box-family results are trustworthy; finite parts are NOT
  until the endpoint-series defect is fixed.  Each run prints loud
  warnings on the affected segments.
- Compare runs with `Scripts/compare_stepwise_log.py <log> --refs
  <results.json>`; loop_package refs need the (-1)^nprops sign rule.
- Dump-replay loop: DIFFEXP_DUMP_LAURENT_DIR capture +
  `Scripts/eval_dump_generic.m` (seconds per iteration).
- Sector-fit diagnostics: DEBUG_SECTOR_FIT=1 (on the WIP branch).
