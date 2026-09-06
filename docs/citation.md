# Citation

When using DiffExp 2.1, cite the software release itself and the method papers
relevant to the workflow used.

For one-dimensional series transport:

> M. Hidding, “DiffExp, a Mathematica package for computing Feynman integrals
> in terms of one-dimensional series expansions,” Computer Physics
> Communications 269 (2021) 108125.
> [doi:10.1016/j.cpc.2021.108125](https://doi.org/10.1016/j.cpc.2021.108125),
> [arXiv:2006.05510](https://arxiv.org/abs/2006.05510).

For the recursive Feynman-trick workflow:

> M. Hidding and J. Usovitsch, “Feynman parameter integration through
> differential equations,” Physical Review D 108 (2023) 036024.
> [doi:10.1103/PhysRevD.108.036024](https://doi.org/10.1103/PhysRevD.108.036024),
> [arXiv:2206.14790](https://arxiv.org/abs/2206.14790).

`CITATION.cff` contains machine-readable metadata for the software repository.

For the sequential-epsilon Chebyshev collocation method used by the spectral backend:

> Yuanche Liu and Yang Zhang, CHESS (2026),
> [arXiv:2606.26691](https://arxiv.org/abs/2606.26691).

The spectral backend is an independent C++ implementation of this method, with
shared-expression preparation, adaptive refinement and DiffExp path/boundary
integration. These implementation changes are not a claim of a new collocation method.
