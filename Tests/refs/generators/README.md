# Pin generators (vendored from /tmp, 2026-06-11)

These scripts mint the pySecDec reference values listed in `Tests/PINS.md`.
They are vendored so a wiped /tmp costs nothing but compute time.  Results
they produced are frozen under `Tests/refs/pysecdec/` — rerun these only to
ADD precision or NEW points, and compare against the frozen copies first.

## Environment rebuild

    python3 -m venv /tmp/pysecdec-venv
    /tmp/pysecdec-venv/bin/pip install pySecDec

Campaign environment (2026-06-10 build): Python 3.14, pySecDec 1.6.6,
pySecDecContrib bundled (FORM 4.3.1, Cuba).  ~5 minutes to install.
Run everything with `/tmp/pysecdec-venv/bin/python3` — the scripts hardcode
that venv path internally (they prepend `pySecDecContrib/bin` to PATH for
`make`; if you relocate the venv, update the two hardcoded
`/tmp/pysecdec-venv` paths in `pysecdec_ft_boxes.py` and
`pysecdec_pentagon.py`).

All scripts cache: if the package's `*_pylink.so` exists they skip
generation/build and go straight to integration, so replays cost seconds to
~1 minute.

## Scripts

### pysecdec_smoke.py
Venv smoke test: generates (does not build/integrate) a 1-loop massive
bubble loop_package in /tmp/pysecdec-smoke.  Runtime: ~1-2 min.  Run this
first after a venv rebuild.

### pysecdec_box1l_pin.py  ->  THE box L0 pin + sign-protocol anchor
1-loop massless on-shell box via loop_package at s=-1, t=-1/3, orders
through eps^1.  Cuhre epsrel 1e-7, epsabs 1e-12, maxeval 5e6.  Prints
`PYSECDEC_BOX1L_PIN <json>`.  Frozen result:
`../pysecdec/box1l_pin_result.log`.  FT box values are independently
verified analytically in the repo test suite, which is what makes this the
anchor that pinned the loop_package-vs-FT sign convention (+1 for 4 props).
Runtime: ~5-10 min first run (build), ~1 min cached.  Workdir:
/tmp/pysecdec-ft-boxes/box1L_pin.

### pysecdec_ft_boxes.py {box_bubble|box_triangle|double_box_planar}
2-loop loop_package route at s=-1, t=-1/3, order eps^0, Minkowski
propagators (REMEMBER: FT = (-1)^nprops x these — PINS.md section 1).
Defaults: Cuhre epsrel 1e-3, epsabs 1e-7, maxeval 2e5; override with
`--maxeval/--epsrel/--epsabs` (DO increase maxeval when re-minting the
box_bubble/box_triangle pins — the eps^0 Cuhre error estimates at the
default were optimistic by ~10x, see PINS.md 2.2).  Prints
`PYSECDEC_JSON <json>`.  Frozen results:
`../pysecdec/loop_package_box_results.jsonl`.  Runtime: ~15-25 min per
example first run (generation + make -j2), minutes cached.
double_box_planar generates but was never integrated to pin quality.
Workdir: /tmp/pysecdec-ft-boxes/<name>.

### pysecdec_pentagon.py  ->  the (PROVISIONAL) pentagon pin
1-loop massless 5-point at {s12,s23,s34,s45,s15} = {-1,-2,-3,-5,-7}
(dot products fixed numerically to match Scripts/FTExamples.m), orders
through eps^1.  Cuhre epsrel 1e-8, epsabs 1e-12, maxeval 3e6.  Frozen
result: `../pysecdec/pentagon_pin_result.txt`.  Measured runtime
2026-06-11: ~21 min end-to-end (generate + make -j3 + integrate).
Workdir: /tmp/pysecdec-ft-pentagon.

## Family route (exact FT normalization — lives in Scripts/, not here)

The second, sign-unambiguous reference route; use it for any per-level
master and for the R9 sign re-verification:

    # 1. export specs (NEEDS THE WOLFRAM KERNEL — orchestrator only):
    env WolframKernel=<kernel> FT_EXAMPLE=box FT_FIXED_VALUE=11/23 \
        PYSECDEC_SPEC_FILE=/tmp/specs.json \
        wolframscript -file Scripts/export_pysecdec_family_specs.m
    # vendored spec files in ../pysecdec/specs_*.json skip this step
    # for banana@11/23 and box@{11/23, 1/20, 2/5, 19/20}.

    # 2. build + integrate + write results.json:
    /tmp/pysecdec-venv/bin/python3 Scripts/pysecdec_family_driver.py \
        ../pysecdec/specs_box.json --output-root /tmp/pysecdec-ft-mybox \
        --factor-monomials

Runtime: ~2-5 min per package (make-dominated); the 11-package box family
~30-45 min fresh.  The driver writes `<output-root>/results.json` in the
format `compare_stepwise_log.py --refs` consumes (minus `_needed` entries —
see ../MANIFEST.md caveat).  The exporter reads $FTConfig
DimensionExpression and cannot encode numerator (negative) indices — it
skips them loudly.

## crosschecks/ — pySecDec-free independent evaluators

sympy/mpmath only (`pip install sympy mpmath`); each runs in seconds to
~2 min.  These broke the box campaign's reference deadlocks and verify the
closed-form pins without any pySecDec build:

| script | what it checks |
|--------|----------------|
| `eval_box_specs.py` | Laurent expansion of exported family specs by direct high-precision quadrature — verified box L1-L3 boundary rows EXACT (Docs/FeynmanTrickBoxFamilyStatus.md:37-39). |
| `eval_2011_exact.py` | eps-expansion of the box needed-integral {2,0,1,1}(t) via z0-endpoint subtraction + log-moment quadrature. |
| `identity_finite_eps.py` | FIRE reduction identity at finite eps = -1/20, t = 1/20: 2-d quad of the spec vs closed-form masters M1/M2. |
| `fourway_eps_test.py` | four-way consistency at eps = -1/20: box integral = Int dt spec_{2011}(t) = pin Laurent = FT dump total — localizes which side of a mismatch is wrong. |
| `spec_hyp.py` | closed-inner-form (2F1) representation of spec_{2011}(t); Laurent towers at three t values "replaces noisy refs". |
