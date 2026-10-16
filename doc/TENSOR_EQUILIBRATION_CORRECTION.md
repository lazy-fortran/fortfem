# Opt-in divergence-free tensor correction

`equilibrated_flux_estimate_tensor(..., correction_sweeps=s)` accepts a
nonnegative number of correction sweeps, default zero. Positive counts require
a fully Dirichlet boundary; requests with Neumann edges are rejected. Existing
calls retain their reconstruction and arithmetic. No new quadrature claim is
made. The method is a practical sampled reconstruction, with its measured
remaining divergence defect included in the majorant.

On each vertex patch, let V_a be the RT space with zero divergence and zero
normal trace on artificial patch edges. Domain Dirichlet edges have free
normal trace. For the current defect e = sigma - A grad u_h, solve the local
A^(-1) orthogonal projection delta_a = -P_(V_a) e in the caller's quadrature
inner product. Update sigma by
sum_a delta_a/3. These projections use the existing generated tensor and RT
leaves, with block-scaled saddle systems. No continuous local formula is
redefined in handwritten code.

The flux remains H(div), and its divergence is unchanged in exact arithmetic.
Also `(e,delta_a) = -||delta_a||^2` and
`||sum_a delta_a||^2 <= 3 sum_a ||delta_a||^2`, since exactly three vertex
patches cover a triangle. Expansion gives
`||e_next||^2 <= ||e||^2 - sum_a ||delta_a||^2/3`.
This proves non-increase of the flux defect, not an h- or anisotropy-uniform
contraction. At fixed polynomial degree, bounded patch valence and fixed s,
work/storage are O(N). Cell evaluation, local patch solves and updates are
OpenMP parallel. The local saddle system is assembled/factored anew each pass;
this implementation makes no claim of small cost relative to the primal solve.

Local divergence-free postprocessing also appears in
[Cai and Zhang (2012)](https://www.math.purdue.edu/~zcai/pdf-paper/2011CaZh2.pdf).
Their jump-robustness result uses a scalar, locally quasi-monotone coefficient;
it does not establish robustness for the rotating, highly anisotropic Sovinec
tensor. The KIN6D study finds index 4.5 at nr=8 but 48 at nr=32, ratio 1e9,
after eight sweeps. A field-aware coarse level/localisation remains open.

`fpm test --target test_equilibrated_flux_tensor_2d --profile release
--flag '-fopenmp -O3'` checks exact polynomial recovery, scalar equivalence,
containment and convergence on rotated constant tensors, strict improvement of
the weighted flux defect and an independently measured divergence defect.
