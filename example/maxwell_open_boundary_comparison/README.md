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

CI generates `maxwell_dtn_modes_1d.png`, `maxwell_reflection_1d.png`, and
`benchmark.txt`; generated media are not committed.
