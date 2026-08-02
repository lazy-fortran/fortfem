# Maxwell open-boundary comparison

This example evaluates the biperiodic planar Maxwell capacity operator on
analytical TE and TM Fourier modes. For tangential wavevector \(\xi\) and
outgoing normal wavenumber
\(\beta=\sqrt{k^2-|\xi|^2}\), the eigenvalue magnitudes are
\(\beta\) for TE and \(k^2/\beta\) for TM.

It also places the exact DtN residual, the measured error from an executable
Nédélec curl-curl PML solve, and a physically larger PML box on one scale. The
two PML solutions are reconstructed at shared interior targets and compared
against the exact complex-stretched plane wave; no hard-coded far-wall value
is used. The comparison follows the transparent boundary formulation of Jiang
et al., arXiv:1811.12449.

The first generated image is the solution itself: the colour field is the
reconstructed Nédélec magnitude on an interior `z` slice, with linear shading
between the displayed physical samples, and the overlaid arrows are the real
in-plane components reconstructed from the solved edge degrees of freedom. The
arrows therefore show the computed vector field rather than a source or
convergence diagnostic. An analytical plane-wave check guards the plotted
sample orientation, so attenuation in `x` cannot silently appear as stripes in
`y`.

The volume-boundary example additionally constructs a tetrahedral box,
reproduces a constant field with first-kind Nédélec elements of orders one
through four, samples its tangential trace on the FFT grid, and pulls the
capacity form back to only the edge and face moments on the selected planar
boundary. A separate tetrahedral box then exercises the complex FortSparse
curl-curl/PML solver and records its solve time and relative edge error.

CI generates `maxwell_pml_field_slice_2d.png`,
`maxwell_domain_comparison_1d.png`, `maxwell_dtn_modes_1d.png`,
`maxwell_reflection_1d.png`, `maxwell_nedelec_dtn_1d.png`, and
`benchmark.txt`; generated media are not committed.
