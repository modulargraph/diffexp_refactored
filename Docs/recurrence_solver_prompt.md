## Context: Solving Systems of Linear ODEs for Feynman Integrals via Recurrence Relations

### Physical Setup

We compute Feynman integrals in dimensional regularization, where the spacetime dimension is $d = 4 - 2\varepsilon$. The integrals satisfy systems of first-order linear differential equations in a line parameter $x$ that parametrizes a path through the space of kinematic invariants. The key structure is that the differential equation matrix admits a Laurent expansion in $\varepsilon$:

$$\partial_x \vec{f}(x,\varepsilon) = \mathbf{A}(x,\varepsilon)\,\vec{f}(x,\varepsilon), \quad \mathbf{A}(x,\varepsilon) = \sum_{k=0}^{K} \mathbf{A}^{(k)}(x)\,\varepsilon^k$$

where $\varepsilon$ is the dimensional regulator (a small parameter). There are no poles in $\varepsilon$ in the matrix (these can always be removed by rescaling). Expanding the integrals as $\vec{f} = \sum_k \vec{f}^{(k)}(x)\,\varepsilon^k$ and collecting order-by-order gives:

$$\partial_x \vec{f}^{(k)} = \mathbf{A}^{(0)}(x)\,\vec{f}^{(k)} + \vec{b}^{(k)}(x)$$

The leading matrix $\mathbf{A}^{(0)}(x)$ (called $\mathbf{M}(x)$ below) is the **same** at every order in $\varepsilon$. The inhomogeneous term $\vec{b}^{(k)}$ depends on the solutions at lower orders in $\varepsilon$ (which are already known when we solve order $k$) and on contributions from other integrals in the system that have already been solved. So the problem reduces to repeatedly solving:

$$\partial_x \vec{g} = \mathbf{M}(x)\,\vec{g} + \vec{b}(x)$$

for blocks of coupled integrals, order-by-order in $\varepsilon$, where $\mathbf{M}$ is a rational matrix function of $x$ (possibly containing square roots of polynomials in $x$), and $\vec{b}$ is known.

### What We Solve With Recurrence Relations

We seek **series solutions** in $x$ around expansion points (typically $x = 0$). Near a **regular singular point**, $\mathbf{M}(x)$ has a simple pole:

$$\mathbf{M}(x) = \frac{\mathbf{M}_{-1}}{x} + \mathbf{M}_0 + \mathbf{M}_1 x + \mathbf{M}_2 x^2 + \cdots$$

where $\mathbf{M}_i$ are constant $p \times p$ matrices (the coefficients are known numerically to high precision). The eigenvalues $\lambda_1, \ldots, \lambda_p$ of $\mathbf{M}_{-1}$ (the **residue matrix**) determine the leading behavior of solutions.

The general ansatz for solutions near a regular singular point is:

$$g_i(x) = \sum_{\alpha} x^{\lambda_\alpha} \sum_{n=0}^{N} \sum_{k=0}^{K_\alpha} f^{(\alpha)}_{n,k}\, x^n\, (\ln x)^k$$

where:
- $\lambda_\alpha$ are determined by the eigenvalues of $\mathbf{M}_{-1}$ (can be rational numbers, possibly differing by integers — "resonance")
- $K_\alpha$ is the maximum logarithm power (determined by Jordan block structure and resonance)
- The coefficients $f^{(\alpha)}_{n,k}$ are vectors in $\mathbb{C}^p$, computed via a **recurrence relation** in $n$

Substituting the ansatz into $x\partial_x \vec{g} = x\mathbf{M}(x)\vec{g} + x\vec{b}$ and collecting powers of $x^{\lambda+n}(\ln x)^k$ yields:

$$[(\lambda + n)\mathbf{I} - \mathbf{M}_{-1}]\,f_{n,k} = -(k+1)\,f_{n,k+1} + \sum_{i=1}^{n} \mathbf{M}_i\,f_{n-i,k} + \beta_{n,k}$$

where $\beta_{n,k}$ comes from the inhomogeneous term $\vec{b}$.

At **non-singular points** ($\mathbf{M}$ is analytic at $x = 0$), the recurrence simplifies: $\lambda = 0$, no logs, and the recurrence is:

$$(n+1)\,f_{n+1} = \sum_{j=0}^{n} \mathbf{M}_j\,f_{n-j} + \beta_n$$

### The Form of the Inhomogeneous Term $\vec{b}$

The inhomogeneous term $\vec{b}^{(k)}(x)$ is constructed from solutions at lower $\varepsilon$-orders and from other integrals already solved. These solutions have the general **decomposed form**:

$$f_i(x) = \sum_{\text{terms}} x^{a + b\varepsilon} \cdot g(x, \varepsilon)$$

where:
- $a \in \mathbb{Q}$: a rational power (can be negative, fractional, e.g., $-1, -1/2, 1/3$)
- $b \in \mathbb{Q}$: an $\varepsilon$-dependent exponent. The factor $x^{b\varepsilon} = e^{b\varepsilon \ln x}$ generates logarithms when expanded in $\varepsilon$: $x^{b\varepsilon} = 1 + b\varepsilon\ln x + \frac{(b\ln x)^2}{2}\varepsilon^2 + \cdots$
- $g(x,\varepsilon) = \sum_k g^{(k)}(x)\,\varepsilon^k$ where each $g^{(k)}(x)$ is a **power series** in $x$ starting at $x^0$ or higher (no negative powers), with coefficients that may contain powers of $\ln x$

So at any given $\varepsilon$-order, $\vec{b}$ is a finite sum of terms of the form:

$$x^a \cdot (\ln x)^j \cdot [\text{power series in } x \text{ with numerical coefficients}]$$

This means $\vec{b}$ is a **known** series with (possibly fractional) leading power and logarithmic factors. The series coefficients are numerical (known to high working precision, typically 30+ digits).

### Currently Implemented Solvers

We currently have the following recurrence-based solvers working:

1. **Non-singular rational recurrence**: For analytic $\mathbf{M}(x)$, uses the simple recurrence $(n+1)f_{n+1} = \sum_j M_j \cdot f_{n-j} + \beta_n$. Can optionally clear denominators when $\mathbf{M}$ is rational to get a fixed-band recurrence.

2. **Singular recurrence (non-resonant)**: For $\mathbf{M}$ with a simple pole and **diagonalizable** $\mathbf{M}_{-1}$ with **non-resonant** eigenvalues (no two eigenvalues differ by a positive integer). Diagonalizes $\mathbf{M}_{-1}$, works in the eigenbasis where the recurrence decouples at leading order.

3. **General singular recurrence (resonant)**: Handles Jordan blocks and resonant eigenvalues (eigenvalues differing by integers). Uses Jordan decomposition, tracks logarithmic terms $(\ln x)^k$, and resolves the solvability conditions that arise when $(\lambda + n)\mathbf{I} - \mathbf{M}_{-1}$ is singular.

### Cases That Need Better Handling

We want to extend the recurrence approach to handle more difficult cases. Here are the situations where we currently fall back to more expensive methods (Frobenius + Wronskian + variation of parameters):

**A. Higher-order poles.** When $\mathbf{M}(x)$ has a pole of order $\geq 2$ (irregular singular point):
$$\mathbf{M}(x) = \frac{\mathbf{M}_{-r}}{x^r} + \cdots + \frac{\mathbf{M}_{-1}}{x} + \mathbf{M}_0 + \cdots$$

The solution structure is much richer: exponential factors $e^{P(1/x)}$, Stokes phenomena, divergent asymptotic series, etc.

**B. The particular solution in the resonant case.** When $\vec{b} \neq 0$ and the exponent of the particular solution resonates with an eigenvalue (i.e., the leading power $s$ of $\vec{b}$ satisfies $s = \lambda_j$ for some eigenvalue $j$), the simple division by $(s + n - \lambda_j)$ breaks down at $n = 0$. Additional logarithmic terms are needed.

**C. Confluent eigenvalues with non-trivial Jordan structure combined with complicated $\vec{b}$.** When the Jordan blocks are large and the inhomogeneous term has logarithmic structure, the bookkeeping of solvability conditions becomes intricate.

**D. Semi-simple but resonant cases where the resonance shifts are large.** When eigenvalues differ by large integers, the fundamental matrix computation requires careful handling of the singular steps at each resonance.

### What We Want

Given:
- The matrix coefficients $\mathbf{M}_{-1}, \mathbf{M}_0, \mathbf{M}_1, \ldots$ (numerical, known to high precision)
- The inhomogeneous term $\vec{b}$ in the form described above (series with fractional powers and logs, numerical coefficients)
- The desired number of terms $N$ in the expansion

Derive and implement the recurrence relations that produce the series coefficients $f_{n,k}$ for the general solution (homogeneous + particular) in the most general regular-singular-point case, including all resonance and Jordan block configurations. Where possible, also handle irregular singular points or derive criteria for when the recurrence approach is applicable vs. when it fundamentally cannot work.

The output should be a general solution of the form:
$$\vec{g}(x) = \vec{g}_{\text{particular}}(x) + \sum_{i=1}^{p} c_i\,\vec{g}_i^{\text{(hom)}}(x)$$

where $c_i$ are undetermined constants (fixed later by boundary conditions), and each component is a sum of terms $x^{\lambda}(\ln x)^k \cdot [\text{series}]$.

### Key Mathematical Questions

1. **Resonant particular solutions**: When the leading power of $\vec{b}$ resonates with eigenvalues of $\mathbf{M}_{-1}$, what is the correct ansatz for the particular solution? How many additional logarithmic orders are needed?

2. **Irregular singular points**: Under what conditions can we still use a recurrence-based approach (e.g., rank-1 irregular singular points where the formal solution involves $e^{c/x}$ times a regular series)? Can we factor out the exponential and reduce to a regular problem?

3. **Optimal solvability resolution**: In the resonant case, when $(\lambda + n)\mathbf{I} - \mathbf{M}_{-1}$ is singular, the solvability condition determines free parameters from earlier steps. What is the most numerically stable way to resolve this?

4. **Fractional powers in $\vec{b}$**: When $\vec{b}$ contains terms like $x^{a}$ with $a$ not an integer, does the particular solution require a ramified ansatz (series in $x^{1/q}$ for some $q$)?
