# Free-boundary port trace gallery

This small gallery fixture samples a manufactured vector trace on a circular
toroidal boundary and sends it through FortFEM's
`assemble_free_boundary_port_residual` contract.  The first image is the
physical boundary: the torus is shown in Cartesian coordinates, point colour
shows the supplied trace magnitude, and short black segments show its tangent
and normal components.

The residual is the positive weighted mismatch

\[
 r_q = w_q\,(t_q-g_q-k_q),
\]

where `t` is the physical trace, `g` is an external target, and `k` is a
manufactured sheet/jump target.  The program also checks the generated JVP by
a centred difference and the VJP by the real dot-product identity.  CSV files
contain the sampled geometry, traces, residual, and derivative diagnostics.

This is deliberately a **neutral numerical contract example**.  It does not
implement an equilibrium solver, coil model, vacuum/BEM/DTN operator, or
free-boundary physics.  Those choices belong to a caller-owned adapter; this
fixture only demonstrates the trace pairing and its differentiable residual.

The contract is documented in the source module
`src/operators/free_boundary_port_residual.f90` and is exposed through the
canonical `fortfem_boundary` facade.
