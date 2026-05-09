# Recursive Fuchsian Local Solver

This note records the intended local algorithm for the next recurrence solver layer.  It is deliberately separated into a symbolic local-connection layer and a finite-width coefficient layer, because those two parts have different failure modes.

## Local Problem

After pulling a multiscale differential equation to one line and choosing a local coordinate `z` at the expansion point, write

```text
theta f = M(z, eps) f,     theta = z d/dz.
```

The desired local model is regular singular.  In a good local basis there is a meromorphic gauge transform

```text
f = T(z, eps) g
```

such that

```text
theta g = B(z, eps) g
```

and `B` is holomorphic at `z = 0`.  If no such finite meromorphic transformation exists, the point is irregular and pure powers/logarithms are not a complete local ansatz.

## Fuchsianization

Fuchsianization can be formulated as stable-lattice saturation.  A local lattice is a free module over the local power-series ring.  If the columns of `T` generate the current lattice, then the transformed connection is

```text
B_T = T^(-1) M T - T^(-1) theta(T).
```

The lattice is stable under the connection exactly when `B_T` has no negative powers of `z`.  Starting from `T = I`, repeatedly:

1. Compute `B_T`.
2. If its minimum `z`-valuation is non-negative, stop.
3. Pick a column with a deepest pole.  This column is the coordinate vector of a connection image not contained in the current lattice.
4. Adjoin that vector to the lattice while keeping a square basis.

Adjoining a vector `u` works as follows.  Let `m = min_i val_z(u_i) < 0` and write `u = z^m(c + O(z))`.  Choose a pivot with nonzero `c_p`, make a constant row operation matrix `E` such that `E c = e_p`, then replace the pivot basis vector by `E u`.  In matrix form, the new basis is `T E^(-1) U`, where `U` is the identity with pivot column replaced by `E u`.

For large systems, the production route should batch the deepest-pole column space, exploit block-triangular sector structure, remove off-diagonal poles with Sylvester recurrences, then run a short global saturation pass.

## Finite-Width Form

Once `B` is holomorphic in `z`, clear only the regular `z`-dependent denominators:

```text
q(z, eps) theta g = C(z, eps) g,
```

where `q(0, eps) != 0`, and `q` and `C` are polynomials in `z` with Laurent-rational coefficients in `eps`.  This finite-width form is the key speed improvement: order `n` depends only on finitely many previous orders `n - j`, not on a length-`n` convolution.

For one exponent sector, use

```text
g = z^alpha Sum[eps^k z^n log(z)^ell c[k,n,ell]].
```

Substitution gives the recurrence

```text
Sum[j,r] q[j,r] ((alpha+n-j)c[k-r,n-j,ell]
                 +(ell+1)c[k-r,n-j,ell+1])
- Sum[j,r] C[j,r] c[k-r,n-j,ell] = 0.
```

Terms with `j = 0` are same-order unknowns and are solved as one small block at fixed `n`.  Terms with `j > 0` are already known.

## Epsilon Poles

If the coefficient matrices contain no negative powers of `eps`, the recurrence is triangular in epsilon order.  If negative epsilon powers occur, the equation for `eps^k` can depend on higher solution coefficients.  The robust algorithm uses a finite work window, enforces equations below the requested output window, and repeats with a larger epsilon buffer until the requested coefficients stabilize.  Non-stabilization is evidence for non-Laurent behavior such as `z^(1/eps)`.

Negative epsilon powers can still be compatible with a Laurent solution if their exponent action is nilpotent.  The finite block solve is designed to cover that case without special-casing the nilpotent normal form.

## Boundary Data

The recurrence determines a local solution space, not a unique solution.  At ordinary points one prescribes `c[k,0,0]`.  At singular points one prescribes leading asymptotic data in each sector.  A production fundamental-solution mode should expose the nullspace of the leading block; the current low-level kernel accepts prescribed coefficient vectors and solves the remaining coefficients.

## Validation

Every local solve should be validated by substituting the truncated result into `q theta g - C g`.  A nonzero residual inside the requested window usually means one of:

- insufficient log cutoff,
- insufficient epsilon work buffer,
- incomplete sector list,
- inconsistent boundary data,
- failed Fuchsianization.

