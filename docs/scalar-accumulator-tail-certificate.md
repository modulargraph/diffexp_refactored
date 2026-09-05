# Scalar integral tail certificate

The disk certificate recognizes the exact one-way system

\[
y'=a(x,\epsilon)y,\qquad z'=w(x,\epsilon)y.
\]

Recognition requires a two-component matrix whose entire second column is proved zero on the witness disk. Each retained epsilon coefficient of `a` and `w` receives its own outward magnitude bound, `A_j` and `W_j`, on the full enclosing square. Denominator and analytic-sheet exclusions are unchanged. A matrix with feedback or an accumulator self-coupling uses the general matrix-norm certificate.

For disk radius `R`, write `b_j` for the magnitude of the initial source coefficient. Radial integral comparison gives

\[
Y_k=e^{R A_0}\sum_{j=0}^k b_j
 [\epsilon^{k-j}]\exp\left(R\sum_{m=1}^k A_m\epsilon^m\right).
\]

All positive-epsilon convolutions are finite. Since the accumulator does not feed the source,

\[
|z_k(x)-z_k(x_0)|\le V_k
=R\sum_{j=0}^k W_jY_{k-j}
\]

on the disk. Its full circle bound remains `|z_k(x0)|+V_k`. The constant initial accumulator has no omitted positive Taylor coefficients, so it is excluded only from the Taylor-tail bound; its arithmetic/input uncertainty stays in the transported polynomial balls.

For retained degree `N` and `q=|h|/R<1`, a uniform omitted-tail bound for the two components is

\[
\max(Y_k,V_k)\,\frac{q^{N+1}}{1-q}.
\]

Consequently large integral weights enter linearly, rather than artificially contributing to an exponential `exp(R*max(|a|,|w|))`. The finite chart adaptation uses `A_0` as this system's exponential growth rate, so it can reduce `h/R` without unnecessarily shrinking a valid disk because `w` is large.

The initial boundary must already enclose the true solution, including all previous tails. No correlation cancellation or boundary independence is assumed. A small tail does not establish small arithmetic uncertainty; final actual-radius checks remain mandatory.

Regression coverage includes an analytic exponential source with a `10^12*(1+eps^2)` forcing, input accumulator uncertainty, rejection of feedback from this special case, and the former double-box/box-triangle scalar budget failures. Their weighted `J_1` representations are compared with the exact identity `J_2/J_1=(2*eps-3)*U/F`. Both now certify at the original `N=80`, `384` working bits, `28` requested absolute digits, and unchanged chart/proof budgets. The isolated double-box large-weight case changed from exhausting512 charts to completing in about0.7 seconds.
