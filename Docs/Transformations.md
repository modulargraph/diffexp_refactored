# Transformations: Mobius, Analytic Continuation, and Frobenius

This document describes three DiffExp subpackages that handle coordinate transformations, branch-cut crossings, and local solution construction near singular points.

---

## Table of Contents

1. [Mobius.m (DiffExp\`Mobius\`)](#mobiusm)
2. [AnalyticContinuation.m (DiffExp\`AnalyticContinuation\`)](#analyticcontinuationm)
3. [Frobenius.m (DiffExp\`Frobenius\`)](#frobeniusm)

---

## Mobius.m

**Package:** `DiffExp`Mobius``
**Dependencies:** `DiffExp`Symbols``, `DiffExp`State``, `DiffExp`Utilities``
**Source:** `DiffExp/Mobius.m` (151 lines)

### Overview

This subpackage provides Mobius (fractional linear) transformations for reparametrizing line segments during numerical transport. The core idea is to map a local coordinate `x` ranging over `(-1, 1)` to a physical interval `(zmin, zmax)` on the integration path, with a designated center point `zmid` mapped to `x = 0`. This conformal mapping concentrates Taylor series coefficients near the center of the interval, improving convergence relative to a naive linear reparametrization.

### Mathematical Background

A Mobius transformation is a map of the form:

```
z(x) = (a*x + b) / (c*x + d)
```

which is uniquely determined (up to an overall scale) by specifying the images of three points. In DiffExp, the three constraints are:

- `z(-1) = zmin` (left singularity bound)
- `z(0) = zmid` (expansion center)
- `z(+1) = zmax` (right singularity bound)

The advantage of Mobius reparametrization over linear reparametrization is that the radius of convergence in the `x`-variable is exactly 1 in all directions (the singularities at `zmin` and `zmax` are mapped to `x = -1` and `x = +1` respectively), making optimal use of the available convergence region.

### Functions

#### GetMobius

```mathematica
GetMobius[{zmin, zmid, zmax}]
```

Returns the Mobius transformation `z(x)` that maps `x in (-1, 1)` to `z in (zmin, zmax)` with `z(0) = zmid`.

**Parameters:**
- `zmin` -- Left boundary (singularity or -Infinity)
- `zmid` -- Center point mapped to `x = 0`
- `zmax` -- Right boundary (singularity or +Infinity)

**Return value:** An expression in `DiffExp`Symbols`x` representing `z(x)`.

**Special cases:**

| Configuration | Formula | Description |
|---|---|---|
| `{-Infinity, zmid, Infinity}` | `x - zmid` | Linear shift (no finite boundaries) |
| `{zmin, zmid, Infinity}` | `(zmid + x*zmid - 2*x*zmin) / (1 - x)` | Right boundary at infinity |
| `{-Infinity, zmid, zmax}` | `(2*x*zmax + zmid - x*zmid) / (1 + x)` | Left boundary at infinity |
| `{zmin, zmid, zmax}` | Full Mobius formula (see below) | Both boundaries finite |
| `{zmin, Infinity, zmax}` | `(zmin*(-1+x) + zmax*(1+x)) / (2*x)` | Center at infinity |

**Full Mobius formula** (both boundaries finite, center finite):

```
z(x) = -((-1+x)*zmid*zmin + zmax*(zmid + x*zmid - 2*x*zmin))
        / ((-1+x)*zmax + zmin + x*(-2*zmid + zmin))
```

The result is passed through `Simplify` for canonical form.

**Verification:** One can verify that `z(-1) = zmin`, `z(0) = zmid`, `z(1) = zmax` for each case.

---

#### GetLineRescaled

```mathematica
GetLineRescaled[line, at, {signsproj, signsim}, nomobius]
```

Produces a rescaled version of the differential equation matrices (represented as an Association of rational functions in `x`) centered at the point `at`, suitable for local series expansion.

**Parameters:**
- `line` -- Association mapping kinematic variables to rational expressions in `x`
- `at` -- The center point for the expansion
- `{signsproj, signsim}` -- Lists of singularity locations projected onto and imaginary-part-sorted along the line
- `nomobius` -- Optional (default `False`). If `True`, forces linear rescaling even when `UseMobius` is enabled.

**Behavior:**

The function operates in two modes depending on the configuration option `UseMobius`:

**Mode 1: Mobius rescaling** (`UseMobius = True`, `nomobius = False`):
1. Finds the nearest singularity to the left (`leftBound`) and right (`rightBound`) of `at` in `signsproj`.
2. Substitutes `x -> GetMobius[{leftBound, at, rightBound}]` into the line.
3. Further rescales by `x -> x / RadiusOfConvergence`.
4. Sets precision to `2 * WorkingPrecision`.
5. Collects numerator and denominator by powers of `x`.

**Mode 2: Linear rescaling** (`UseMobius = False` or `nomobius = True`):
1. Finds `leftBound` and `rightBound` as above.
2. Computes `minDistance = Min[at - leftBound, rightBound - at]`.
3. Substitutes `x -> at + x * minDistance`.
4. Further rescales by `x -> x / RadiusOfConvergence`.
5. Sets precision to `2 * WorkingPrecision`.
6. Expands and reconstructs as an Association.

**Note:** The computation is wrapped in `Block[{$MaxExtraPrecision = 1000}]` to allow high-precision intermediate calculations.

---

#### Center Point Functions

These functions compute center points for the "predivision" segmentation strategy. When a line segment between two singularities is too long for a single Taylor expansion to converge adequately, it is subdivided into smaller segments. The center point functions determine where to place the expansion centers of these sub-segments.

The spacing is controlled by the `DivisionOrder` parameter `k`, which determines how close to a singularity the center of an adjacent segment can be placed. Higher values of `k` place centers closer to singularities (shorter segments, more segments needed).

##### GetMobiusCPL / GetMobiusCPR

```mathematica
GetMobiusCPL[{zmin, zbound, zmax}]
GetMobiusCPR[{zmin, zbound, zmax}]
```

Compute the Mobius-optimal center point to the left (CPL) or right (CPR) of a boundary singularity `zbound`, given the interval `(zmin, zmax)`.

**Parameters:**
- `zmin` -- Left boundary of the interval (or -Infinity)
- `zbound` -- The singularity being stepped away from
- `zmax` -- Right boundary of the interval (or +Infinity)

**Return value:** The coordinate of the next center point.

**Formulas for GetMobiusCPL:**

| Configuration | Formula |
|---|---|
| `{-Infinity, zbound, zmax}` | `(-zbound + k*zbound + 2*zmax) / (1 + k)` |
| `{zmin, zbound, Infinity}` | `(zbound + k*zbound - 2*zmin) / (-1 + k)` |
| `{zmin, zbound, zmax}` | `(zbound*zmax + k*zbound*zmax + zbound*zmin - k*zbound*zmin - 2*zmax*zmin) / (2*zbound - zmax + k*zmax - zmin - k*zmin)` |

**Formulas for GetMobiusCPR:**

| Configuration | Formula |
|---|---|
| `{-Infinity, zbound, zmax}` | `(zbound + k*zbound - 2*zmax) / (-1 + k)` |
| `{zmin, zbound, Infinity}` | `(-zbound + k*zbound + 2*zmin) / (1 + k)` |
| `{zmin, zbound, zmax}` | `(zbound*zmax - k*zbound*zmax + zbound*zmin + k*zbound*zmin - 2*zmax*zmin) / (2*zbound - zmax - k*zmax - zmin + k*zmin)` |

Here `k = DivisionOrder` from the configuration.

##### GetCPL / GetCPR

```mathematica
GetCPL[{zmin, zbound, zmax}, k]
GetCPR[{zmin, zbound, zmax}, k]
```

Compute linear (non-Mobius) center points to the left or right of `zbound`. These use a geometric spacing strategy that accounts for the relative position of `zbound` within `(zmin, zmax)`.

**Parameters:**
- `{zmin, zbound, zmax}` -- Interval endpoints and singularity
- `k` -- Optional division order (defaults to `DivisionOrderVal` from configuration)

**Algorithm:** These functions set up constraint equations involving a step size `s` and solve for the new center point `xnew`. The constraints encode that the ratio of the distance from the new center to the singularity versus the distance to the interval boundary equals `1/k`. Different formula branches handle the cases where `zbound` is in the left half, right half, or at infinity.

##### GetCPLRep

```mathematica
GetCPLRep[MyEq]
```

Helper function that extracts the value of `xnew` from a set of constraint equations by iterated rule application.

##### FindNextCenterPointL / FindNextCenterPointR

```mathematica
FindNextCenterPointL[xbc, singsproj]
FindNextCenterPointR[xbc, singsproj]
```

Top-level center point finding functions that dispatch to the appropriate Mobius or linear variant.

**Parameters:**
- `xbc` -- Current boundary/center position
- `singsproj` -- Sorted list of singularity positions projected onto the line

**Behavior:**
1. Finds `leftBound` and `rightBound` as the nearest singularities to the left and right of `xbc`.
2. If `UseMobius = True`: calls `GetMobiusCPL` or `GetMobiusCPR`.
3. If `UseMobius = False`: calls `GetCPL` or `GetCPR`.

---

### Configuration Parameters

| Parameter | Default | Description |
|---|---|---|
| `UseMobius` | `False` | Whether to use Mobius transformations (vs. linear) |
| `DivisionOrder` | `3` | Controls center point spacing (higher = closer to singularities) |
| `RadiusOfConvergence` | `1` | Additional rescaling factor applied after Mobius/linear transformation |

---

### Example Usage

```mathematica
(* Compute the Mobius transformation for a segment between singularities at 0 and 1,
   centered at 1/3 *)
mobius = GetMobius[{0, 1/3, 1}]
(* Result: a rational function in x that maps x=0 to z=1/3,
   x=-1 to z=0, x=+1 to z=1 *)

(* Find the next center point to the right of x=1/2,
   given singularities at {0, 1, 4} *)
nextCenter = FindNextCenterPointR[1/2, {0, 1, 4}]
```

---

## AnalyticContinuation.m

**Package:** `DiffExp`AnalyticContinuation``
**Dependencies:** `DiffExp`Symbols``, `DiffExp`State``, `DiffExp`Utilities``, `DiffExp`SeriesOps``
**Source:** `DiffExp/AnalyticContinuation.m` (107 lines)

### Overview

This subpackage handles the analytic continuation of solutions around branch points. When the integration path crosses a branch cut, the multi-valued functions (logarithms, fractional powers) must be continued to the correct Riemann sheet. The continuation is encoded via Heaviside-like theta functions (`theta_p`, `theta_m`) that track which sheet the solution lives on, depending on whether the singularity is approached from above or below the real axis (as determined by the `i*delta` prescription).

### Mathematical Background

Consider a branch point at `z = z0` with an `i*delta` prescription indicating that the physical value lies at `z0 + i*delta` (i.e., approached from above). When the line parameter `x` passes through zero (the expansion center mapped to the singularity), multi-valued functions pick up monodromy phases:

- **Logarithm:** `Log(x)` acquires a `2*pi*i` shift when `x` crosses the negative real axis.
- **Square root:** `x^(1/2)` flips sign under continuation around the origin.
- **General fractional power:** `x^b` acquires a phase `exp(-2*pi*i*b)`.

The theta functions encode the result in a form valid on both sides:

```
theta_p = 1  if Im(x) > 0    (approaching from above)
theta_m = 1  if Im(x) < 0    (approaching from below)
```

with the constraint `theta_p + theta_m = 1` (exactly one is active at any point).

### Functions

#### PrepareAnalyticContinuation

```mathematica
PrepareAnalyticContinuation[Line]
```

Prepares the replacement rules that encode how multi-valued functions transform when the integration path crosses a branch cut at the current singularity.

**Parameters:**
- `Line` -- An Association representing the current line segment (mapping kinematic variables to expressions in `x`)

**Side effects:** Sets the following state variables:
- `AnalyticContinuationReplacements` -- List of replacement rules
- `AnalyticContinuationReplacementsAssociation[Line]` -- Cached result for this line
- `CurrentSingularityWasAddedFromSquareRoot` -- Boolean flag
- `CurrentSingularityHasIDeltaPrescription` -- Boolean flag
- `AnalyticContinuationFailed` -- Set to `True` if sign determination fails

**Algorithm:**

1. **Identify vanishing factors:** Determines which entries in `DeltaPrescriptions` vanish on the given line (i.e., which branch points coincide with `x = 0` on this segment).

2. **Check square root origin:** Tests whether the singularity was automatically added from a square root in the differential equations (rather than being user-specified). This is tracked via `SquareRootPrescriptionsAdded[]`.

3. **Determine sign of i*delta:** For each vanishing prescription, computes the sign of the imaginary part of `x` as it approaches the singularity. This determines whether the approach is from above (+1) or below (-1).

4. **Construct replacement rules:** If the sign is -1 (approaching from below), the following replacement rules are generated:

| Original | Replacement | Condition |
|---|---|---|
| `Logx` | `(theta_p + theta_m)*Logx - 2*pi*i*theta_m` | Logarithm |
| `x^b` where `Denominator[b] = 2` | `x^(b-1/2) * (theta_p - theta_m) * Sqrt[x]` | Square root |
| `x^b` where `Denominator[b] > 2` | `(theta_p + Exp[-2*pi*i*b]*theta_m) * x^b` | Higher root |

5. **Cache results:** Stores the replacements in `AnalyticContinuationReplacementsAssociation` keyed by the line, avoiding redundant computation on revisits.

**Error handling:** If the signs are ambiguous (marked "?") or inconsistent (multiple distinct signs found), the function sets `AnalyticContinuationFailed = True` and defaults to sign = +1 (no monodromy).

**Interpretation of replacement rules (sign = -1):**

- **Logarithm:** The replacement `Logx -> (theta_p + theta_m)*Logx - 2*pi*i*theta_m` encodes:
  - For `theta_p = 1` (above): `Logx -> Logx` (no change)
  - For `theta_m = 1` (below): `Logx -> Logx - 2*pi*i` (standard branch cut)

- **Square root:** The replacement `x^(1/2) -> (theta_p - theta_m)*Sqrt[x]` encodes the sign flip:
  - For `theta_p = 1`: `+Sqrt[x]`
  - For `theta_m = 1`: `-Sqrt[x]`

- **General fractional power:** The phase factor `exp(-2*pi*i*b)` is the monodromy of `x^b` around the origin.

---

#### ProjectThetas

```mathematica
Project\[Theta]s[Expr, f]
```

Projects an expression containing multiple occurrences of theta functions back into a canonical linear form in `(theta_p, theta_m)`. After algebraic manipulations, expressions may contain products like `theta_p * theta_m` or higher powers `theta_p^2`, which are unphysical. This function uses the idempotency relations (`theta_p^2 = theta_p`, `theta_m^2 = theta_m`, `theta_p * theta_m = 0`) to simplify.

**Attributes:** `Listable` (automatically maps over lists).

**Parameters:**
- `Expr` -- Expression potentially containing `theta_p` and `theta_m` symbols
- `f` -- Optional simplification function to apply (default: `Expand`)

**Algorithm:**
1. If no theta functions are present, returns `f[Expr]` directly.
2. Extracts the `theta_m` part: sets `theta_p -> 0`, `theta_m -> 1`.
3. Extracts the `theta_p` part: sets `theta_m -> 0`, `theta_p -> 1`.
4. Recombines as:
   ```
   f[theta_p_part + theta_m_part] / 2
   + f[theta_p_part - theta_m_part] / 2 * (theta_p/2 - theta_m/2)
   ```

**Mathematical justification:** Any function `h(theta_p, theta_m)` subject to the constraints `theta_p + theta_m = 1` and `theta_p * theta_m = 0` can be written as `a + b*(theta_p - theta_m)` where:
- `a = (h(1,0) + h(0,1)) / 2` (symmetric part)
- `b = (h(1,0) - h(0,1)) / 2` (antisymmetric part)

The function applies `f` (typically `Expand`) and `PChop` (precision chopping) during recombination for numerical stability.

---

### State Variables

| Variable | Type | Description |
|---|---|---|
| `AnalyticContinuationReplacements` | List | Current replacement rules |
| `AnalyticContinuationReplacementsAssociation` | Association | Cache of replacements keyed by Line |
| `AnalyticContinuationFailed` | Boolean | Flag indicating sign determination failure |
| `CurrentSingularityWasAddedFromSquareRoot` | Boolean | Whether current singularity comes from a square root |
| `CurrentSingularityHasIDeltaPrescription` | Boolean | Whether current singularity has an i*delta prescription |

---

### Example Usage

```mathematica
(* Configure delta prescriptions *)
LoadConfiguration[{
  DeltaPrescriptions -> {{t - 4, 1}, {t - 9, -1}},
  (* ... other options ... *)
}];

(* During transport, when approaching t=4 from below: *)
PrepareAnalyticContinuation[<|t -> 4 + x|>]

(* The resulting AnalyticContinuationReplacements will contain rules
   for continuing Logx and fractional powers across the branch cut *)

(* After applying replacements, project theta functions *)
result = expr /. AnalyticContinuationReplacements;
result = Project\[Theta]s[result];
```

---

## Frobenius.m

**Package:** `DiffExp`Frobenius``
**Dependencies:** `DiffExp`Symbols``, `DiffExp`State``, `DiffExp`Utilities``, `DiffExp`SeriesOps``, `DiffExp`Integration``
**Source:** `DiffExp/Frobenius.m` (141 lines)

### Overview

This subpackage implements the Frobenius method for constructing local series solutions to linear ordinary differential equations (ODEs) near regular singular points. It is used within DiffExp to find a basis of homogeneous solutions (the Wronskian matrix columns) near each singularity on the integration path.

### Mathematical Background

Consider an n-th order linear ODE:

```
a_n(x) * y^(n)(x) + a_{n-1}(x) * y^(n-1)(x) + ... + a_1(x) * y'(x) + a_0(x) * y(x) = 0
```

At a regular singular point (Fuchsian singularity) `x = 0`, the coefficients `a_i(x)` have at most pole singularities of a controlled order. The Frobenius ansatz seeks solutions of the form:

```
y(x) = x^r * sum_{i=0}^{N} c_i * x^i
```

where `r` is determined by the **indicial equation** -- a polynomial equation obtained by substituting `x^r` into the leading-order terms of the ODE.

The indicial equation has `n` roots (counting multiplicity). The largest root `r_max` always yields a valid Frobenius solution. Smaller roots may require logarithmic terms (handled implicitly through the order reduction procedure in `FrobeniusSolutions`).

### Functions

#### Frobenius1

```mathematica
Frobenius1[DEqn]
```

Finds one Frobenius solution corresponding to the largest root of the indicial equation.

**Parameters:**
- `DEqn` -- The differential equation in DiffExp internal form: a list `{a_0, a_1, ..., a_n}` of series coefficients, where the ODE is `sum_i a_i(x) * y^(i)(x) = 0`.

**Return value:** A `SeriesData` object representing the Frobenius solution as a power series in `x`.

**Algorithm:**

1. **Fuchsian check:** Computes the leading powers of all coefficient series. Verifies that `LeadingPower[a_i] - LeadingPower[a_n] >= i - n` for all `i` (the Fuchsian condition). Aborts with an error if violated.

2. **Series expansion:** Expands all coefficients `a_i(x)` to the configured `ExpansionOrder` using `DiffExpSeries`.

3. **Indicial equation:** Constructs the indicial polynomial by substituting the trial solution `x^r` into the leading-order terms:
   ```
   P(r) = sum_i LeadingCoeff[a_i] * (d^i/dx^i x^r) * x^(-r)
   ```
   This yields a polynomial in `r` whose roots are the characteristic exponents.

4. **Largest root:** Solves `P(r) = 0` and takes `r_max = Max[roots]`. Issues a warning if `Denominator[r_max] > 2` (higher-order roots are less thoroughly tested).

5. **Ansatz construction:** Builds the truncated series:
   ```
   y(x) = x^r_max * (1 + c_1*x + c_2*x^2 + ... + c_N*x^N) + O(x^{ExpansionOrder})
   ```
   where `N = Ceiling[ExpansionOrder - r_max]` and `c_0 = 1` (normalization).

6. **Linear system:** Substitutes the ansatz into the ODE, collects coefficients of each power of `x`, and solves the resulting triangular linear system for `{c_1, c_2, ..., c_N}` using `LinearSolve`.

7. **Cross-check:** If numerical instability is detected (via Mathematica's `Check` mechanism), substitutes the computed solution back into the ODE to verify it satisfies the equation to working precision. Reports an error if validation fails.

8. **Return:** The series solution truncated to the appropriate order.

**Configuration dependencies:**
- `ExpansionOrder` -- Determines the number of series terms computed
- `ChopPrecision` -- Controls numerical zero detection
- `WorkingPrecision` -- Precision of the `LinearSolve` computation

---

#### FrobeniusSolutions

```mathematica
FrobeniusSolutions[DEqn]
```

Finds a complete basis of `n` independent Frobenius solutions for an n-th order ODE by recursive order reduction.

**Parameters:**
- `DEqn` -- The differential equation in DiffExp internal form (list of coefficient series)

**Return value:** A list `{y_1, y_2, ..., y_n}` of `SeriesData` objects representing the independent solutions.

**Algorithm:**

1. **First solution:** Calls `Frobenius1[DEqn]` to obtain the solution `y_1` corresponding to `r_max`.

2. **Order reduction (for n > 1):** Uses the multiplicative ansatz `y(x) = y_1(x) * w(x)` to reduce the order. Substituting into the ODE and using Leibniz's rule:
   ```
   DEqnReduced[k] = sum_{i=k}^{n} Binomial(i-1, k-1) * a_i(x) * y_1^(i-k)(x)
   ```
   The first entry `DEqnReduced[1]` should vanish (since `y_1` is a solution). This is cross-checked when the crosscheck flag `"FrobeniusSolutions"` is active.

3. **Recursive solution:** The reduced equation (with the vanishing first entry removed) is a differential equation of order `n-1` in `w'(x)`. Recursively calls `FrobeniusSolutions` on this reduced system.

4. **Integration:** The solutions `w'_i(x)` of the reduced equation are integrated via `DiffExpIntegrate` to obtain `w_i(x)`, then multiplied by `y_1(x)` to get the remaining independent solutions:
   ```
   y_{k+1}(x) = y_1(x) * Integrate[w'_k(x)]
   ```

5. **Final cross-check:** When the crosscheck flag `"FrobeniusSolutions"` is active, substitutes all solutions back into the original ODE to verify they satisfy it. Checks that the residuals are zero (within `CrosscheckChopPrecision`). Reports an error if any solution fails validation.

**Recursive structure:**
```
n-th order ODE
  -> Frobenius1 gives y_1
  -> reduce to (n-1)-th order ODE in w'
     -> Frobenius1 gives w'_1
     -> reduce to (n-2)-th order ODE in u'
        -> ...
        -> base case: 1st order ODE -> Frobenius1 directly
```

**Key series operations used:**
- `SExpand` -- Series expansion with precision management
- `SMultiply` -- Multiplication of two series
- `SD` -- Series differentiation
- `DiffExpIntegrate` -- Series integration (handles `Logx` terms)

---

### Configuration Parameters

| Parameter | Default | Description |
|---|---|---|
| `ExpansionOrder` | `50` | Number of terms in the Frobenius series |
| `ChopPrecision` | `250` | Digits below which coefficients are set to zero |
| `WorkingPrecision` | `500` | Internal computation precision |
| `CrosscheckLevel` | `0` | Level >= 1 activates `"FrobeniusSolutions"` cross-checks |

---

### Example Usage

```mathematica
(* The Euler equation x^2*y'' + x*y' - y = 0
   In DiffExp internal form: {a_0, a_1, a_2} = {-1, x, x^2}
   After series expansion, find solutions: *)

DEqn = DiffExpSeries[{-1 + O[x]^50, x + O[x]^50, x^2 + O[x]^50}];
solutions = FrobeniusSolutions[DEqn]
(* Returns {x^1 + O[x^50], x^(-1) + O[x^50]} *)

(* For the Bessel equation of order 0: x^2*y'' + x*y' + x^2*y = 0
   Indicial roots: r = 0 (double root)
   Frobenius1 gives J_0(x) as a power series
   FrobeniusSolutions gives {J_0(x), Y_0(x)} where Y_0 contains Logx terms *)
```

---

### Interaction with Other Subpackages

The Frobenius solutions are used by the `Wronskian.m` subpackage to construct the Wronskian matrix at each singular point. This matrix is essential for the variation-of-parameters method used in `Transport.m` to propagate boundary conditions along the integration path.

The `DiffExpIntegrate` function from `Integration.m` handles the appearance of `Logx` terms during the order reduction step. When the characteristic exponents differ by an integer, the integration step naturally produces logarithmic solutions (the second-kind solutions).

---

## Interrelation of the Three Subpackages

These three subpackages work together during the transport process:

1. **Mobius.m** handles the geometric setup: how line segments are parametrized and subdivided between singularities.

2. **Frobenius.m** constructs the local solution basis at each singular point where a new series expansion begins.

3. **AnalyticContinuation.m** manages the transition between Riemann sheets when the path crosses a branch cut, ensuring the transported solution remains on the correct physical sheet.

The typical call sequence during `TransportTo` is:
1. The line is segmented (using Mobius center point functions).
2. At each segment center near a singularity, `Frobenius1`/`FrobeniusSolutions` provides the homogeneous solutions.
3. `GetLineRescaled` transforms the matrices to the local coordinate.
4. When a branch point is encountered, `PrepareAnalyticContinuation` sets up the monodromy rules.
5. After transport through the branch point, `ProjectThetas` simplifies the resulting expressions.
