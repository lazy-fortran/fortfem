# Coupled Fourier transport residual majorant

`fortfem_meridian_fourier_equilibration` estimates the full retained real
Fourier-Galerkin energy, using the same uniform phi rule as
`meridian_fourier_diffusion` and `meridian_fourier_load`. The measure is
`R dR dZ` times the phi mean; on a cylinder R is its metric radius.
The reference solution shares the candidate's discrete Dirichlet trace.
This is a practical estimate: neither meridional nor phi quadrature is enclosed.
It does not bound the Fourier complement or boundary approximation error.

Let `j = K grad u_h - f` and let a subscript b denote its real Fourier
coefficients under the phi rule. The meridional flux target is
`tau_b = R (j_R, j_Z)_b`. Retain `(j_phi)_b` as the toroidal flux and set
`s_b = (R f0)_b + (d_phi j_phi)_b`. In particular, for positive n,

```
s_cos = R f0_cos + n j_phi_sin
s_sin = R f0_sin - n j_phi_cos.
```

The RT reconstruction sigma_b has continuous meridional normal trace and
`div_2 sigma_b + Pi_k s_b = d_b`, where the existing patch estimator measures
d_b. Thus the reconstructed three-dimensional flux
`(sigma_R/R, j_phi, sigma_Z/R)` has the required balance mode by mode up to
the explicit data oscillation and measured divergence terms. The RT
coefficients and toroidal samples are returned for inspection. Algebraic
patch defects are retained; a floating-point balance is not claimed exact.

For a test function in the retained Fourier space the residual is

```
r(v) = sum_b w_b [ (s_b - Pi s_b, v_b) + (d_b, v_b)
                  + (sigma_b - tau_b, grad_2 v_b) ],
w_0 = 1, w_cos = w_sin = 1/2.
```

Discrete phi integration by parts follows from orthogonality of retained
modes (`nphi >= 2N+1`), including when the coefficient is not band limited.
If the caller supplies a domain-wide lower bound `K >= chi_lower I`, then
`||v||_energy^2 >= sum_b w_b int R chi_lower |grad_2 v_b|^2`.
Cauchy-Schwarz, the patch estimator's local oscillation bounds and the
supplied Friedrichs inequality give

```
||u - u_h||_energy <= sqrt(sum_b w_b flux_b^2)
                   + sqrt(sum_b w_b oscillation_b^2)
                   + sqrt(sum_b w_b divergence_b^2).
```

This argument does not discard coefficient coupling: K acts on the full
candidate gradient before projection. It also does not require separate
mode solves or a global coupled reconstruction solve. The guaranteed
ellipticity constant is a caller contract, not inferred from samples.
Efficiency can deteriorate with anisotropy because the defect norm uses
`(R chi_lower)^-1`, rather than the inverse coupled tensor.

Sampling uses FortNum complex FFTs with caller-owned per-thread workspaces.
At fixed spatial degree, storage is O(cells * modes), reconstruction is
O(cells * modes), and transforms are O(cells * nphi log nphi). There is no
dense mode-pair matrix. Pointwise tensor residuals, toroidal gradients and
source coupling are emitted by `gen_meridian_fourier_equilibration_products`
using `tools/codegen/fortsym.lock`. The estimator is not differentiated.

`test_meridian_fourier_equilibration` checks P4 exact recovery with a nonzero
strong source, cosine and sine components, and even/odd phi grids. A P2
Galerkin solve for a quartic meridional bubble and phi-dependent scalar K
checks containment, efficiency <= 1.6 and energy convergence >= 1.8.
An independent direct trigonometric sum at twice the phi sample count
computes the error. Optional arguments `mesh_intervals repetitions` run a
larger timing case without expanding the fast test. Solve time includes
matrix/load assembly, Dirichlet elimination and SuperLU; estimator time
includes sampling, reconstruction and all reported terms.

Local gfortran -O3 evidence (shared host, N=2, nphi=8): exact P4 recovery
2.27e-11; P2 efficiency 1.226/1.217/1.213 on 4/8/16 intervals and energy
rates 1.971/1.992. On 64 intervals (83,205 DOFs), three estimator calls
average 2.53 s at one thread and 0.818 s at four threads (3.09x), with an
unchanged 4.2688e-5 bound and 1.210 efficiency. Four-thread assembly plus
SuperLU takes 5.64 s, so estimator/solve = 0.145. An earlier shared-host
repeat gave 6.36/1.88 s (3.39x), estimator/solve = 0.189. These are practical
sampled bounds; the cost is the generic estimator, excluding consumer
phi-dependence detection. The FFT workspace API requires FortNum acc8dce4
or later; consumer validation uses a source override, without changing pins.
