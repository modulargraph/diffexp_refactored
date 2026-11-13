# Integral families from the 2026 literature

Feynman integral families from 2026 papers that cite DiffExp
([arXiv:2006.05510](https://arxiv.org/abs/2006.05510)), set up as DiffExp
Feynman-trick calculations. Each family is computed from its propagator list
alone. DiffExp then:

- reduces the integrals with the built-in IBP Solver (no FIRE, Kira or Mathematica);
- builds the Feynman-parameter differential equations;
- fixes the boundary values;
- transports them to the requested point.

No boundary values or differential equations from the papers are used.

This directory has the 17 families that run to completion. Where a paper prints
values at a point DiffExp can reach, they are reproduced (Check column).
Families that do not work yet are listed [at the end](#not-yet-working).

## Layout

| Path | Contents |
| --- | --- |
| `families/<name>.json` | `DiffExp.FeynmanFamily/v1` input: propagators, kinematics, dimension, requested integrals |
| `families/<name>.meta.json` | source paper, the paper's propagators and labels, mapping notes, kinematic point, published reference values and the conversion used to compare with them |
| `families/physical/<name>.json` | the paper's physical point of a family shipped at a Euclidean point, set up to be reached from that Euclidean point (`continuation`); see [below](#physical-points-by-continuation) |
| `tools/add_continuation.py` | adds such a `continuation` object to a physical-point family |
| `results/<name>.json` | DiffExp response and wall time of the run in the table |
| `run_families.py` | runs all (or selected) families with a per-family time limit |
| `make_table.py` | rebuilds the table below from `results/` |
| `compare_references.py` | compares results with the papers' printed numbers (`comparison` in the metadata) |
| `causal_explore.py` | at a physical point, reruns a family with each choice of an i0 side DiffExp cannot derive |
| `compute_one.sh` | minimal shell example for one family |
| `ComputeFamily.wl` | the same from Mathematica through `DiffExp.m` |

## Running

Build DiffExp (see the top-level README). The scripts use `build/diffexp` of
this checkout; set `DIFFEXP` to use another executable.

```sh
cd examples/literature-2026
python3 run_families.py --threads 8                 # all families, one after another
python3 run_families.py hvp_chpt_e5e6_d2            # one family
python3 compare_references.py                       # digits agreeing with the papers
python3 make_table.py --write                       # refresh the table
sh compute_one.sh families/hvp_chpt_e5e6_d2.json    # raw JSON response
```

Exact reductions are cached under `$DIFFEXP_CACHE`
(default `~/.cache/diffexp-literature-2026`), so a second run of a family only
repeats the numerics.

## Conventions

- DiffExp's Feynman-trick integrals use:
  - the measure `d^dk / (i pi^{d/2})` per loop, without `e^{eps gamma_E}`;
  - propagators `1/(m^2 - q^2)`, positive in the Euclidean region.

  A paper with `1/(q^2 - m^2)` propagators differs by `(-1)^{sum of powers}`.
  `mapping_notes` and `normalization` in each metadata file give the relation
  to the paper's labels. The `comparison` object gives the exact coefficients,
  epsilon shift and `exp(c eps)` factor applied before comparing.
- Values in the table are DiffExp's own normalization.
- **Point**: most papers quote numbers only in the physical region. Several of
  those points need an i0 side that DiffExp 2.2 cannot yet derive on
  multi-loop levels (see below). Such families are shipped at a Euclidean point
  (`*_euclid`); the paper's point stays in the metadata as `paper_point` or in
  `mapping_notes`.
  - "paper's": the paper's own point.
  - "Euclidean (chosen here)": all F coefficients positive, real result.
- Cut integrals are computed as their uncut parent families.

## Results

Run with DiffExp 2.2.1 on an AWS c8g.16xlarge (Graviton4, 64 vCPU, 123 GB),
8 IBP threads per family and 8 families side by side. The wall times are
therefore measured on a fully loaded machine, from a cold reduction cache, with
default precision and the full epsilon expansion. "Leading ε power" is the
lowest order with a nonzero coefficient.

<!-- results-table:start -->
| Family | Paper | What | Point | Status | Prep [s] | Numerics [s] | Total [s] | Leading ε power | ε⁰ coefficients | Check |
|---|---|---|---|---|---|---|---|---|---|---|
| `diphoton_heavyloop_pla_euclid` | [2606.30354](https://arxiv.org/abs/2606.30354), [2509.15315](https://arxiv.org/abs/2509.15315) | 2-loop 4-point, massless external legs, closed heavy-quark loop (mass m) + massless lines; family involves … | s,t < 0 < u (complex; i0 sides derived by DiffExp) | ✅ works | 7.1 | 72.8 | 88.4 | -2 | -0.928587618979 + 0.336701088709 i; 1.64451176526 + 0.040332014888 i; 3.81376139557 - 2.47105887201 i |  |
| `dr_mixed_ewqcd_3l_ladder_tbW` | [2603.06571](https://arxiv.org/abs/2603.06571) | 3-loop 2-point ladder (LA) topology: closed fermion loop (top on one side of a W rung, massless b on the … | Euclidean (chosen here) | ✅ works | 18.8 | 82.5 | 112 | -1 | 1.01905166284; -3.56481577782 |  |
| `hh_qqbar_ew_doublebox_PL` | [2601.16924](https://arxiv.org/abs/2601.16924), [2606.25928](https://arxiv.org/abs/2606.25928), [2503.16620](https://arxiv.org/abs/2503.16620) | 2-loop planar 4-point double box, 2 massless + 2 Higgs legs (mH^2), 3 internal lines of mass mV; 45 … | Euclidean (chosen here) | ✅ works | 4.3 | 24.4 | 28.9 | -2 | -0.084329205182; 0.326031225536; 0.919307816855 |  |
| `hvp_chpt_e5e6_d2` | [2603.15252](https://arxiv.org/abs/2603.15252) | 3-loop 2-point, all lines mass M, elliptic kite with bubble on the rung, d=2-2ε | Euclidean (chosen here) | ✅ works | 10.5 | 17.2 | 27.9 | 0 | 0.460352338672; 0.222971397090 | E1 subintegral = shipped banana (38 digits); 55-digit rerun stable |
| `hvp_chpt_e5e6_d4` | [2603.15252](https://arxiv.org/abs/2603.15252) | 3-loop 2-point, all lines mass M, elliptic kite with bubble on the rung, d=4-2ε | Euclidean (chosen here) | ✅ works | 15.6 | 55.0 | 74.2 | -3 | -3.80077909070; -0.040055288205 | leading pole exactly 1/3; paper gives no tables |
| `np_triangle_elliptic_b` | [2606.26691](https://arxiv.org/abs/2606.26691), [2305.13951](https://arxiv.org/abs/2305.13951) | p1^2=p2^2=0, s=(p1+p2)^2; D2,D4 massive (m), others massless; D7=l1^2 ISP | paper's | ✅ works | 6.6 | 10.9 | 19.2 | 0 | 0.206288513718; 0.669763290327; 1.02520066359 | vs paper: M4 24 digits (1 orders); M6 24 digits (1 orders) |
| `np_triangle_elliptic_b_y0p01` | [2606.26691](https://arxiv.org/abs/2606.26691), [2305.13951](https://arxiv.org/abs/2305.13951) | p1^2=p2^2=0, s=(p1+p2)^2; D2,D4 massive (m), others massless; D7=l1^2 ISP | paper's | ✅ works | 6.6 | 13.6 | 21.9 | 0 | 0.031909164539; 0.250959890331; 0.344602371530 | vs paper: M4 24 digits (1 orders); M6 22 digits (1 orders) |
| `pentagon_hp_p1` | [2610.41872](https://arxiv.org/abs/2610.41872) | finite one-loop scalar pentagon with off-shell spacelike legs and massive lines, 'Euclidean' kinematics | paper's | ✅ works | 3.9 | 6.6 | 10.5 | 0 | 1.95853777099e-11 | vs paper: pentagon (LoopTools quad) 16 digits (1 orders) |
| `pentagon_hp_p2` | [2610.41872](https://arxiv.org/abs/2610.41872) | finite one-loop scalar pentagon with off-shell spacelike legs and massive lines, 'Euclidean' kinematics | paper's | ✅ works | 3.7 | 6.6 | 10.4 | 0 | 2.34275324706e-12 | vs paper: pentagon (LoopTools quad) 17 digits (1 orders) |
| `radreturn_pl1_euclid` | [2607.03343](https://arxiv.org/abs/2607.03343) | 2-loop planar 4-point double box with electron mass on 6 lines, off-shell leg q^2; contains the … | Euclidean (chosen here) | ✅ works | 11.9 | 44.4 | 58.9 | -2 | 0.462617192757; 0.020360347190; 1.73159654367 |  |
| `radreturn_pl2_euclid` | [2607.03343](https://arxiv.org/abs/2607.03343) | 2-loop planar 4-point (e+e- -> gamma gamma*), electron mass m on 5 internal lines, one off-shell leg q^2; … | Euclidean (chosen here) | ✅ works | 70.6 | 109 | 185 | 0 | 0.260690998893; 1.73704275634; 0.300032928591 |  |
| `sunrise_two_mass_sr_d2` | [2609.08875](https://arxiv.org/abs/2609.08875) | 2-loop sunrise with masses m,M,m in d=2-2eps; elliptic. (Families D and TT of the paper are not usable: … | Euclidean (chosen here) | ✅ works | 0.1 | 0.2 | 0.3 | -2 | 1.76612007887; 1.87994862507; 2.31128991446 |  |
| `ttgamma_PBE_euclid` | [2610.42062](https://arxiv.org/abs/2610.42062) | 2-loop planar 5-point pentagon-box, 2 on-shell tops + photon + 2 massless; 3 massive internal lines; 154 … | Euclidean (chosen here) | ✅ works | 22.4 | 141 | 181 | -4 | -0.754736846260; -0.686415343127; -1.65367421688 |  |
| `ttgamma_PBF_euclid` | [2610.42062](https://arxiv.org/abs/2610.42062) | 2-loop planar 5-point pentagon-box, 2 on-shell tops + photon + 2 massless; 3 massive internal lines; 142 … | Euclidean (chosen here) | ✅ works | 27.4 | 317 | 355 | -4 | 1.18397137244; 3.28028515273 |  |
| `tth_lightquark_pentabox_T1` | [2611.12661](https://arxiv.org/abs/2611.12661) | 2-loop planar 5-point penta-box, two on-shell tops (m_t) + off-shell boson q^2 + 2 massless legs; massive … | Euclidean (chosen here) | ✅ works | 13.2 | 64.7 | 78.3 | -4 | -0.212887047758; 0.251199310342 |  |
| `ttj_pentabox_B_euclid` | [2603.05336](https://arxiv.org/abs/2603.05336), [2404.12325](https://arxiv.org/abs/2404.12325) | outgoing p1,p2 tops (p^2=mt2), p3,p4,p5 massless; D1=k1^2-mt2, D2=(k1-p2)^2, D3=(k1-p2-p3)^2, … | Euclidean (chosen here) | ✅ works | 26.9 | 117 | 147 | -4 | 6.60508458647; 2.54815575708 |  |
| `vbf_oneloop_pentagon_massive` | [2609.37919](https://arxiv.org/abs/2609.37919) | 1-loop 5-point pentagon, two adjacent internal lines of mass m (W/Z) flanking massive leg q1^2=s1 (Higgs), … | Euclidean (chosen here) | ✅ works | 0.4 | 0.5 | 0.9 | -2 | -0.084122326392 |  |
<!-- results-table:end -->

Notes:

- The ε⁰ column truncates each coefficient to 12 digits; `results/` holds the
  full Arb balls of every order.
- Paper timings are not comparable. The papers evaluate known canonical
  differential equations from boundary values computed beforehand (AMFlow,
  pySecDec). Here DiffExp derives everything, including the boundary values,
  from the propagators.

## Physical points by continuation

DiffExp 2.2 can now reach a physical point P from a Euclidean one. A family's
`continuation` object gives a Euclidean direction V; DiffExp starts at
E = P + tV, where every F coefficient is positive, and transports the masters
along 0 → ih → 1 + ih → 1 in the upper half plane of the line E + x(P − E).
Along that line every F coefficient decreases, so the endpoint is the `F − i0`
value without choosing any i0 side. The IBP Solver now provides the derivative
in a kinematic direction (varying external Gram), and the boundary at E is a
separate ordinary run. See `docs/feynman-families.md`, "Continuation from a
Euclidean point".

Checks:

- These continuations agree to every printed digit with DiffExp's derived i0
  sides (`tests/causal.cmake`): the bubble above threshold, the massless box
  at a physical point, and the sunrise above threshold.
- The sunrise above threshold is also derived directly now, through a new
  statement on the Symanzik F ("dF/dx + c F has one sign").

`families/physical/` holds the paper points of five shipped Euclidean families
with a continuation set up: `radreturn_pl1`, `radreturn_pl2`,
`ttj_pentabox_B`, `ttgamma_PBE` and `ttgamma_PBF`.

- Their directions keep the physical masses fixed, together with on-shell
  p² = m². Otherwise the line crosses m² = 0 next to P.
- For `radreturn_pl1` with m² ≈ 2.6·10⁻⁷, that crossing was 2.6·10⁻⁸ from P
  and the transport could not reach the point.
- They have not been run to completion yet. The boundary at E alone takes the
  85 masters of `radreturn_pl1` about 20 minutes through ε⁴ on a 10-core Mac,
  and the larger families need hours of IBP reduction per line.

Run one with:

```sh
build/diffexp ft examples/literature-2026/families/physical/radreturn_pl1.json --json
```

## Not yet working

These are kept out of this directory for now. Each line gives the paper, the
family, and DiffExp 2.2.1's stopping reason.

**Physical-point variants of shipped Euclidean families** (derived i0 sides;
the first five now have a continuation set-up in `families/physical/`, not yet
run to completion). The i0 side of a real
singular point on a multi-loop level cannot be derived from the level's
Symanzik F. `causal_explore.py` tried both sides of every such point; no
combination completed.

| Family (paper point) | Paper | Stopping reason |
| --- | --- | --- |
| `diphoton_heavyloop_pla_x1`, `_npa_x1`, `_npa_euclid` | [2606.30354](https://arxiv.org/abs/2606.30354) | with either side: "unsupported residue spectrum: eigenvalue is nonlinear in epsilon". Massless 2→2 (s + t + u = 0): no Euclidean region, so no continuation |
| `ttgamma_PBE` | [2610.42062](https://arxiv.org/abs/2610.42062) | with either side at level 1: 3 h timeout |
| `ttgamma_PBF` | [2610.42062](https://arxiv.org/abs/2610.42062) | each choice reveals a further underivable point on level 1 (x = 0.9668, 0.9907, 0.9969, …) |
| `radreturn_pl1` | [2607.03343](https://arxiv.org/abs/2607.03343) | with either side at level 2: 3 h timeout |
| `radreturn_pl2` | [2607.03343](https://arxiv.org/abs/2607.03343) | point at x = 3.4e-12 (from m² ≈ 2.6e-7); with either side: "artifact payload exceeds store size limit" |
| `ttj_pentabox_B` | [2603.05336](https://arxiv.org/abs/2603.05336), [2404.12325](https://arxiv.org/abs/2404.12325) | with either side at level 3: "artifact payload exceeds store size limit" |
| `np_double_box_massless` | [2607.18406](https://arxiv.org/abs/2607.18406) | level 1 at x = 1/4 is fixed; it now stops at a later level whose side is still undecided. Massless 2→2 nonplanar: no Euclidean region, so no continuation either |

**Families that do not complete at any point tried.**

| Family | Paper | Stopping reason |
| --- | --- | --- |
| `pbb_3l5p`, `pbp_3l5p` (+ symmetric points) | [2512.17330](https://arxiv.org/abs/2512.17330), [2607.08477](https://arxiv.org/abs/2607.08477), [2606.26691](https://arxiv.org/abs/2606.26691) | IBP pivot-fill / preparation budget exhausted (3-loop 5-point) |
| `hbbbb_3l_elliptic_k3_uncut` | [2607.03398](https://arxiv.org/abs/2607.03398) | modular preparation time budget exhausted on level 1 |
| `ladder3l_two_mass` | [2609.37536](https://arxiv.org/abs/2609.37536) | modular preparation time budget exhausted on level 1 |
| `np_double_pentagon_hjj_top` | [2610.01154](https://arxiv.org/abs/2610.01154) | modular preparation time budget exhausted on level 1 |
| `sixpt_2l_double_pentagon` (physical and Euclidean, d = 6 − 2ε) | [2603.16831](https://arxiv.org/abs/2603.16831) | artifact store size limit; the Euclidean run exceeded 87 GB |
| `qed_4pt_3l_triple_ladder` | [2602.10807](https://arxiv.org/abs/2602.10807) | negative ε-valuation cycle |
