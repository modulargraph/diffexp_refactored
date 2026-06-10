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

## Remaining defect: endpoint-segment SERIES at eps^0+

Pinned by the pointwise reference {2,0,1,1}(t = 1/20) (family exporter
at FT_FIXED_VALUE=1/20 + driver):

- pySecDec: -63.1579/eps + 1357.9 +/- 142 (real)
- dump integrand series at t = 0.05 (inside endpoint segment 1):
  eps^-1 = -126.3158 = 2 x (-63.1579) EXACT (prefactor ratio 2);
  eps^0 complex/wrong (about -349 - 471 i under the local-side theta
  + real-log reading).

So the boundary-anchored endpoint segments (1 and 12) of the L1
transport represent the integrand correctly at eps^-1 but WRONGLY from
eps^0 up.  All integrator-level explanations are exhausted: interior
segments verified against NIntegrate, endpoint towers verified
frame-independently (the q0 = 0 moment diagonal is invariant under the
complex-log <-> real-log rebasing), and the salvage windows are honest.
The defect is in how the endpoint series are CONSTRUCTED (transport ->
boundary handoff / local sector matching), one level upstream.

Diagnostic facts for the continuation:
- endpoint towers show an EXACT confluent structure on the q0 = 0
  moment diagonal: moments {24,0,0,24,72,144,240} = 12(k-1)(k-2), i.e.
  characteristic polynomial (r-1)^3 (triple root) relative to the
  extracted b = -1 - while the true function (M1, M2 above) has two
  DISTINCT sectors t^(-eps) and t^0.  A triple root at +1 on top of
  b = -1 is NOT explicable by the true sector content; it is whatever
  the boundary builder emitted.
- branch confluent-sectors-wip holds a full N-root/confluent Prony
  generalization of FitResidualEndpointSectors (works, explains the
  full tower, but integrating the reassigned sectors loses the eps^-2
  pole - do not merge until the upstream representation is understood).
- The warnings ("omitted N endpoint coefficients... starting at eps^k")
  fire on exactly these segments; trust windows are honest.

## Practical guidance

- POLES of box-family results are trustworthy; finite parts are NOT
  until the endpoint-series defect is fixed.  Each run prints loud
  warnings on the affected segments.
- Compare runs with `Scripts/compare_stepwise_log.py <log> --refs
  <results.json>`; loop_package refs need the (-1)^nprops sign rule.
- Dump-replay loop: DIFFEXP_DUMP_LAURENT_DIR capture +
  `Scripts/eval_dump_generic.m` (seconds per iteration).
- Sector-fit diagnostics: DEBUG_SECTOR_FIT=1 (on the WIP branch).
