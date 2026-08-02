# Regularized surface-current layer

This physical-first slab example resolves a tangential surface current
\(\mathbf K=(0,1200,-450)\;\mathrm{A/m}\) as the Gaussian volume profile

\[
\mathbf J_\epsilon(d)=
\frac{\exp[-(d/\epsilon)^2]}{\sqrt{\pi}\epsilon}\,\mathbf K,
\qquad \epsilon=0.015\;\mathrm m.
\]

The first plot shows the two nonzero Cartesian components of the evaluated
volume current against an independently calculated analytical Gaussian. It is
the physical sheet profile, not a convergence graph. The secondary plot shows
that positive trapezoidal normal weights recover the prescribed integrated
surface current as the layer is resolved.

The executable checks the Gaussian profile, monotone integrated-current
convergence, and the finest-grid relative error before producing any files.
It also records the data behind both plots and a small machine-readable timing
record. The colors are from a color-vision-safe palette, while solid lines and
distinct oracle markers retain the distinctions in grayscale.

Outputs, in gallery order:

- `regularized_sheet_profile_1d.png`
- `regularized_sheet_integral_convergence.png`
- `regularized_sheet_profile.csv`
- `regularized_sheet_integral_convergence.csv`
- `benchmark.json`
