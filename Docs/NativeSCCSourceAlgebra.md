# Native SCC source algebra

This slice is the typed finite-series foundation for moving off-diagonal SCC
propagation out of Wolfram.  It does not yet schedule target-block recurrence
solves; it provides the operations those solves consume without coefficient
JSON round trips.

## Prepared rational multiplier

Exact classification remains on the Wolfram side for the first milestone.
For a rational entry `c(t,eps)`, Wolfram proves that no noncentral pole lies
inside the chart and prepares

```text
c(t,eps) = eps^jmin t^-M sum[j>=0,n>=0] q[j,n] eps^j t^n.
```

For an input local solution with epsilon window `[kmin,kmax]` and Taylor width
`N+1`, C++ retains exactly `kmax-kmin+1` epsilon kernels and `N+1` Taylor
coefficients.  The result has

```text
epsilon window = [kmin+jmin, kmax+jmin]
sector tag     = (a-M, b, p)
```

and coefficients

```text
out[e,n] = sum[j=0..e] sum[m=0..n] q[j,m] in[e-j,n-m].
```

This is the finite triangular recurrence used by
`SectorSeries`MultiplyRational`; coefficients above the complete window are
never filled with assumed zeros.

## Sparse coupling matrices

A prepared sparse matrix stores `(target row, source column, multiplier)`
records.  Application selects only the referenced source component, performs
the prepared rational convolution, embeds it into the target component, and
combines every entry proved nonzero by exact preparation.  Combination uses

```text
min power    = minimum input min power
complete max = minimum input complete max
Taylor max   = minimum input Taylor max
```

with explicit epsilon-power alignment.  Identical exact `(a,b,p)` tags are
merged.  Integer-spaced `a` towers are deliberately not compacted yet:
compaction is optional, and doing it safely requires proving both the exact
integer difference and that a finite Taylor shift discards no nonzero tail.

Source presence is structural provenance.  A prepared multiplier carries an
explicit `proven_zero` fact from exact Wolfram preparation; C++ never infers
absence from an all-zero Acb slab.  A numerically zero term from an exact
nonzero coupling still participates in the honest-window intersection,
matching the Wolfram rule that inexact zeros remain active.

The algebra rejects nonempty error envelopes for now.  Dropping an incoming
certificate would be incorrect; the future native transport certificate must
instead apply the multiplier norm and transfer sensitivity explicitly.

## Next ownership step

The session must retain a composite SCC chart containing:

- the exact graph bound to structural nonzeros of the original coupling
  matrices, not merely a self-consistent caller-supplied certificate;
- one prepared diagonal recurrence operator per component;
- prepared physical coupling matrices and sequential gauge/spectral source
  transforms;
- exact epsilon-halo requirements and target-sector schedules.

It can then walk the topological order, keep predecessor `LocalSolution`
slabs native, form each target source with this algebra, and invoke the target
particular recurrence without returning an intermediate slab to Wolfram.
Pseudo-resonant compensation remains a mandatory part of that transaction;
it cannot be deferred until after a truncated source has crossed the API.
