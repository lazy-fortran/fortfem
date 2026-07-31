# PARALLAX-aligned FCI parallel diffusion

This small fixture exercises FortFEM's field-coordinate-independent support
operator with identity interpolation maps. It applies

\[
L_\parallel u=-W_c^{-1}Q^T W_sK_\parallel Q u
\]

to a smooth cosine profile on an open field-line segment, then checks the
independent telescoping mass oracle and the negative weighted-energy identity.
The zero-gradient end behaviour keeps the boundary flux small in the plot.
The maps are trivial by design; field-line tracing and interpolation geometry
remain separate services, as in the PARALLAX architecture.

Outputs:

- `fci_parallel_profile_1d.png`
- `fci_parallel_dissipation_1d.png`
- `fci_parallel_profile.csv`
- `fci_parallel_benchmark.csv` with the measured matrix-free action time on
  the local runner.

The implementation follows the support-operator construction described by
[Stegmeir et al.](https://doi.org/10.1016/j.cpc.2015.09.016). FortFEM does not
copy PARALLAX source code or benchmark data.
