# Configuration-space oracle for banana integrals

This note fixes the normalization of the independent one-dimensional Bessel
oracle used for Euclidean banana integrals.  It does not depend on FIRE,
DiffExp, or the recurrence solver.

## Exact normalization

Let a banana have `L` loops and `N = L + 1` massive lines.  After Wick
rotation, the package measure and propagator convention become

```text
B_L(nu) = Integral Product[a=1..L, d^d k_a / pi^(d/2)]
          Product[j=1..N, (q_j^2 + m_j^2)^(-nu_j)],
q_N = Q - Sum[a=1..L, k_a],   Q = Sqrt[-p^2] > 0.
```

There is no `Exp[EulerGamma eps]` in this definition.  Such a factor in a
canonical master is a separate basis normalization.

With the Fourier convention

```text
G_nu(x;m) = Integral d^d k/(2 pi)^d Exp[i k.x]/(k^2+m^2)^nu
```

the massive propagator is

```text
G_nu(x;m) = m^(d/2-nu) x^(nu-d/2) K_(d/2-nu)(m x)
            / ((2 pi)^(d/2) 2^(nu-1) Gamma[nu]).
```

The `L` momentum integrals enforce equality of the `N` configuration-space
coordinates.  Thus

```text
B_L(nu) = (4 pi)^(L d/2) Integral d^d x Exp[-i Q.x]
          Product[j=1..N, G_(nu_j)(x;m_j)].
```

Performing the angular Fourier transform gives the exact Hankel
representation

```text
B_L(nu) = C Integral_0^Infinity dr
                    r^(Sum[nu_j] - L d/2)
                    J_(d/2-1)(Q r)
                    Product[j=1..N, K_(d/2-nu_j)(m_j r)],

C = 2^(L d/2 + N - Sum[nu_j]) Q^(1-d/2)
    Product[j=1..N, m_j^(d/2-nu_j)/Gamma[nu_j]].
```

The apparent `Q` singularity at `Q = 0` is removable by taking the small
argument limit of the Bessel `J` factor.

For `d = 2 - 2 eps` and unit propagator powers this simplifies to

```text
B_L(eps) = 2^(L (1-eps)) Q^eps Product[j=1..N, m_j^(-eps)]
           Integral_0^Infinity dr r^(1+L eps) J_(-eps)(Q r)
                                      Product[j=1..N, K_eps(m_j r)].
```

In particular,

```text
B_L(0) = 2^L Integral_0^Infinity dr r J_0(Q r)
             Product[j=1..N, K_0(m_j r)].
```

For the unequal four-loop point used by `banana4_unequal`, the constants in
the propagators are squared masses:

```text
m_j^2 = {2, 3/2, 4/3, 5/4, 1},   Q = 1.
```

Consequently,

```text
B_4(0) = 16 Integral_0^Infinity dr r J_0(r)
         K_0(Sqrt[2] r) K_0(Sqrt[3/2] r) K_0(2 r/Sqrt[3])
         K_0(Sqrt[5] r/2) K_0(r)

       = 28.647167500233570163016705540503359...
```

## Epsilon coefficients and endpoint subtraction

The massive unit-power integral in two dimensions is finite and analytic at
`eps = 0`.  Its small-`r` Frobenius branches have powers

```text
r^(1 + (L - 1 + Sum[j, s_j]) eps),   s_j in {-1,+1},
```

so epsilon derivatives produce only integrable powers of `Log[r]`.  A direct
way to obtain coefficients without differentiating Bessel functions with
respect to their order is Cauchy sampling.  For
`eps_j = rho Exp[2 pi i j/M]`,

```text
[eps^k] B_L(eps) ~= rho^(-k)/M Sum[j=0..M-1,
  Exp[-2 pi i j k/M] B_L(eps_j)].
```

For four loops, `rho = 0.05` is conservative: the direct integral converges
in the strip `-1/L < Re[eps] < 1`.  Double the number of samples and repeat at
`rho/2` to measure aliasing and quadrature error.  The samples are independent
and can be distributed over available workers.  Recovering coefficient `k`
loses roughly `k Log10[1/rho]` decimal digits, which must be included in the
working-precision budget.

For general powers or genuine Laurent cases, expand `J` at the lower endpoint
and retain both noninteger Frobenius branches

```text
K_mu(z) = pi (I_(-mu)(z) - I_mu(z))/(2 Sin[pi mu])
```

until after they have been combined.  This preserves cancellations of fake
`1/eps` poles when the order tends to an integer.  Subtract a finite sum of
terms `c(eps) r^(a+b eps)`, integrate them analytically as

```text
c(eps) R^(a+b eps+1)/(a+b eps+1),
```

and numerically integrate only the smooth remainder.  Derivatives with
respect to `a` handle explicit logarithms.  This procedure exposes genuine
Laurent poles through the denominators.  At the upper endpoint all positive
masses give exponential decay `Exp[-r Sum[m_j]]`; a finite cutoff can be
certified with an incomplete-gamma tail bound.

Positive propagator powers and analytic powers such as
`nu_j = n_j + alpha_j eps` are supported directly by the general formula.
Exactly zero or negative integer powers are distributions/contact terms in
configuration space and should instead be treated as pinched topologies or
carefully regulated limits.  Irreducible numerators require configuration-
space derivatives and tensor angular reduction; they are not represented by
the scalar product of Bessel `K` functions above.

## Standalone finite oracle

The companion script evaluates the finite `d = 2` moment in both `r` and
`r = t^2` coordinates and checks their agreement:

```sh
python3 Scripts/banana_bessel_oracle.py
python3 Scripts/banana_bessel_oracle.py --digits 50
python3 Scripts/banana_bessel_oracle.py \
  --mass-squared 2 3/2 4/3 5/4 1 --q 1
```

It requires only Python and `mpmath`.  `--masses` may be used instead when
the physical masses, rather than their squares, are already known.

The configuration-space treatment of arbitrary-mass sunrise/banana graphs
and their Laurent expansion is described by Groote, Koerner, and Pivovarov,
[arXiv:hep-ph/0403122](https://arxiv.org/abs/hep-ph/0403122).  The same
one-dimensional Bessel moments are used as independent differential-equation
checks in [arXiv:1204.0694](https://arxiv.org/abs/1204.0694).
