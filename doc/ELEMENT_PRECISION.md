# Batched element kernels and exterior integration

[V] `fortfem_element_batch:p1_element_batch(jac, diffusion, reaction, blocks,
status)` evaluates affine triangular P1 diffusion plus reaction matrices.
`jac(nel,4)` stores the columns of the physical Jacobian as `[a,b,c,d]`,
corresponding to `[[a,c],[b,d]]`. Coefficients are cellwise constants.
`blocks(nel,3,3)` stores element, local row, local column. The element index
is contiguous for CPU vector loads and GPU coalescing. Both windings are
accepted; the integration measure is the absolute determinant.

[V] Input kind selects storage; output kind selects arithmetic. Supported
pairs are real64/real64, real32/real32, and real32/real64. There is no global
precision switch. Status is per cell: 0 success, 1 nonfinite input,
2 singular or numerically singular geometry, 3 nonfinite result. Invalid
blocks are zeroed. Shape errors stop. Callers must handle nonzero statuses.
`csr_assemble_batched` in `fortfem_element_to_csr` consumes the same blocks
and the connectivity used with `csr_build_pattern`; it converts real32
entries in registers and accumulates into the existing fp64 CSR values.
Rows have independent owners, so shared nodes need no atomic additions.
The existing exchange-form diagonal omission is preserved.

[V] `tools/codegen/app/gen_element_precision_products.f90` emits fp64/fp32
P1 matrices, P2 values/gradients/Hessians, affine determinants, and covariant
and contravariant Piola factors. The existing reference generator shares
`triangle_scalar_bases` with this generator. Derivatives come from FortSym.
The runtime has no symbolic dependency. Regeneration uses `fortsym.lock`;
`check_generated.sh` compares both the Fortran and CUDA outputs.
[H] Migration of the remaining element families and fp16/bf16 storage awaits
the shared FortSym/FortGen precision contracts. The pinned emitter exposes
real64, real32, and mixed arithmetic; this change implements no private half
representation. Native half arithmetic and a posteriori precision escalation
are not provided. Estimators and accepted residuals remain consumer-owned fp64.

[V] Build the narrow CMake target with `FORTFEM_ELEMENT_KERNELS_ONLY=ON`;
link `fortfem::element_kernels`. Enable `FORTFEM_BUILD_TESTING` for its CTest.
Enable `FORTFEM_BUILD_BENCHMARKS` for `fortfem_element_precision_scaling`.
Its variants 1/2/3 are fp64/fp32/mixed cell evaluation; 4/5 assemble fp64/fp32
blocks into fp64 CSR. Reported times are medians of five trials on a fixed
524288-cell problem. Run separate processes with OMP_NUM_THREADS=1,2,4,8,16
and OMP_PROC_BIND=close, OMP_PLACES=cores for strong scaling.

[V] CUDA is opt-in: `FORTFEM_ENABLE_CUDA=ON` requires an explicit
`FORTFEM_CUDA_ARCHITECTURES` list. `fortfem_element_cuda.h` provides asynchronous
fp32, fp64, and mixed launches on a caller-supplied stream. Arrays are device
resident, with the same SoA layout; the library allocates and transfers no
memory. This prototype evaluates cell matrices; global GPU scatter/solves
are outside its scope. The scalar device leaves use FortSym emission.

[V] CPU/CUDA tests use independent barycentric gradients and polynomial
reproduction, reversed winding, singular cells, and shared-node scatter.
[H] The test tolerance model is gamma_96(u) * chi(J)^2 * ||K||_max,
where gamma_n(u)=n*u/(1-n*u), u=epsilon/2, and
chi(J)=(sum(abs(J)))^2/abs(det(J)). The operation budget includes determinant,
reciprocal, gradient products, coefficient scaling and storage conversion.
It assumes normal, round-to-nearest arithmetic and well-conditioned cells
with nonnegative diffusion/reaction. This is a conservative regression model,
not an interval certificate or a proved bound for arbitrary cells. It follows
the separation of geometry, arithmetic and assembly errors in
[Croci and Wells](https://arxiv.org/abs/2410.12614); no theorem from that paper
is claimed for this implementation.

[V] FortFEM already provides per-mode exterior closure through
`toroidal_mode_boundary_operators_2d`: `tm_bem_assemble`,
`tm_bem_steklov_symmetric`, and `tm_bem_exterior_neumann`. The normal points
into the exterior and the flux is R*d_n(u); S_n*u is its negative. These
operators require a surface of revolution. They are distinct from the GS
poloidal-vector-potential kernel. General 3D triangulated surfaces have
`assemble_laplace_calderon_p1_p0_3d` and
`evaluate_laplace_representation_triangles_3d`; curved torus implementations,
virtual casing, and harmonic-period machinery also exist in `src/bem`.

[H] For KIN6D F4, first adapt the fixed Gamma_c surface of revolution to the
existing per-mode operators. Verify sign, R measure, n=0 gauge, energy and
n>0 exact harmonic traces through KIN6D's core residual. General chart-free
shells require a separate integration: oriented mesh validation, trace/flux
transfer, Calderon solve, toroidal harmonic periods, coil superposition,
JVP/VJP of the shell coupling and fp64 exterior error estimates. A scalar
Laplace potential alone does not establish the full multiply connected
magnetic exterior. Keep coil physics and shell coupling in KIN6D. Qualify
sphere/cube monopoles, torus harmonic periods and near-boundary quadrature
before hierarchical acceleration. This exceeds the bounded precision change;
no new exterior solver is introduced here.
