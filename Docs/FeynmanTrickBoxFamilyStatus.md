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

## RESOLVED: the eps^0 deficit (apparent-singularity solver fixes)

Box L0 now reproduces the pin to 11 significant digits
(-41.28416739754647 vs -41.28416739757452; eps^-2/-1 exact as before).
Fixed 2026-06-11 in DiffExp/IntegrationStrategies/Recurrence.m.

Root cause: the L2 transport line crosses an APPARENT singularity of the
DE at xx2 = 7/11 (the root of the quadratic Symanzik coefficient
beta(t) = -A*B - A/3 + 4B/3, A = 1-11t/23, B = 1-t) - the first interior
singular chart any example ever exercised (banana had none; the box L1
system is diagonal with poles only at 0 and 1).  The function is
analytic there, so all genuine local x^-1/x^-2 content is zero by
analyticity; what broke was numerics/bookkeeping in the singular
recurrence's particular solution, in compounding layers:

1. Upstream series arithmetic leaves inert one-past-the-end
   SeriesCoefficient[sd, {x, 0, nmax}] tail requests (numerically
   ~1e-25) in the eps>=2 sources; their symbolic presence turns the
   source entry's head into Plus.
2. ComputeSingularParticular probed Head === SeriesData and silently
   fell back to leading power 0, so the genuine x^-1 source coefficient
   was never consumed: the particular lost a CONSTANT (-2.367 at eps^2),
   boundary fixing absorbed the mismatch into the x^-1 homogeneous mode,
   and the wrong-by-a-constant function propagated through the matching
   chain into the L1 anchor of master {1,0,0,1} (eps^2 off by -1.11524),
   rode the L1 transport as an exact constant (order-2 deviations are
   x-independent when eps^0/1 are exact), entered the combined eps^1
   residue, and the x^(-1-eps) endpoint sector integration's 1/eps
   enhancement deposited it in the final eps^0: the -6.6914.
3. The source assembly also leaves a cancellation-residue x^-2
   coefficient (~1e-29 RELATIVE at WP 300, precision-tracked at ~271
   digits - far above both PChop and 10^(-ChopPrecision/2)); reading
   nmin naively then shifts the ansatz onto the eigenvalue (s - lambda
   = 0).
4. The non-resonance guard used IntegerQ on numeric eigenvalue
   differences (always False), and the particular recursion loops had no
   zero-divisor guard - a software-zero division could poison silently.
5. The particular's SD * x^s assembly can come out head-Plus when
   coefficients carry theta/Logx content; downstream only consumes
   SeriesData shapes and silently dropped such particulars.

Fixes: normalize compound sources back to SeriesData (zeroing
out-of-window SeriesCoefficient tails, loudly otherwise); theta-aware,
Logx-probed RELATIVE leading-coefficient skip with threshold
Max[10^(-ChopPrecision/2), 10^-24] (the 1e-24 floor is load-bearing);
numeric-aware resonance guard; zero-divisor guards (defer to the general
solver); particular normalized to SeriesData on return.  Env-gated debug
hooks added: DEBUG_SING_PART=1 (routing/extraction telemetry),
DEBUG_DUMP_DISPATCH_DIR (per-call ctx/bVec dumps),
DEBUG_POWER_COEFF=1 (silent-zero events).

Validation: the L2 upper-endpoint limit tower is exact to 20 digits at
ALL eps orders (Gamma[1+eps]B(1-eps,1-eps)(4/23)^-eps); battery 17/17;
bubble/sunrise/banana comparators 0 failures; box L0 matches the pin.

Exonerated along the way (each by direct measurement): theta/branch
conventions (twice), the combination assembly (verified exact on its
inputs to 13 digits), EvaluateLimitFromTransport (faithful to stored
data), the eps-lookahead budget (warnings gone, value unchanged), the
exported matrices (all slices exact; matrix linear in eps),
PowerCoefficient's silent zeroing (3 benign boundary events), and the
frozen-SC tails' numeric value (~1e-25) - their TYPE was the poison,
not their size.

Latent issues catalogued for follow-up (subagent audit reports):
the general singular solver returns an empty particular (with complex
leaks at higher orders) for resonant-with-source inputs - currently
bypassed, never properly fixed; CheckEpsPoles guard tests the d-form
and can never fire (FeynmanTrickIteration.m:454); negative-index
boundary requests fall through to "direct" and can silently delete
numerator powers; EvaluateLimitFromTransport silently zeroes
endpoint-divergent IBP coefficients; the limit/direct paths lack the
integrate path's incomplete-top-order trimming and
ShiftRawBoundariesToFinite zero-pads unknown orders; hardcoded
expOrd = 30 in the symbolic assembly branch; Normal+re-Series launders
one phantom top order (also in ApplyAnalyticContinuation); the inert
SeriesCoefficient tail emission itself; eps_prefactors.m is written but
never read.

## Practical guidance

- POLES of box-family results are trustworthy; finite parts are NOT
  until the endpoint-series defect is fixed.  Each run prints loud
  warnings on the affected segments.
- Compare runs with `Scripts/compare_stepwise_log.py <log> --refs
  <results.json>`; loop_package refs need the (-1)^nprops sign rule.
- Dump-replay loop: DIFFEXP_DUMP_LAURENT_DIR capture +
  `Scripts/eval_dump_generic.m` (seconds per iteration).
- Sector-fit diagnostics: DEBUG_SECTOR_FIT=1 (on the WIP branch).
