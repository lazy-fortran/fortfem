# Maxwell open-boundary comparison

This example evaluates the biperiodic planar Maxwell capacity operator on
analytical TE and TM Fourier modes. For tangential wavevector \(\xi\) and
outgoing normal wavenumber
\(\beta=\sqrt{k^2-|\xi|^2}\), the eigenvalue magnitudes are
\(\beta\) for TE and \(k^2/\beta\) for TM.

It also places the exact DtN residual, the analytically predicted reflection
of the quadratic PML used by the scalar and vector tests, and an undamped
far-wall reflection on one scale. The comparison follows the transparent
boundary formulation of Jiang et al., arXiv:1811.12449.

The volume-boundary example additionally constructs a tetrahedral box,
reproduces a constant field with first-kind Nedelec elements of orders one
through four, samples its tangential trace on the FFT grid, and pulls the
capacity form back to only the edge and face moments on the selected planar
boundary. This is the same compact block accepted by the complex FortSparse
curl-curl/PML solver.

CI generates `maxwell_dtn_modes_1d.png`, `maxwell_reflection_1d.png`,
`maxwell_nedelec_dtn_1d.png`, and `benchmark.txt`; generated media are not
committed.
