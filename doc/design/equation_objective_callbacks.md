# Equation/objective callback dispatch

`fortfem_equation_objective_callbacks` composes caller-owned callbacks with
the deterministic metadata and packing contract in
[`equation_objective_registry`](equation_objective_registry.md).  For each
declared block, the value, JVP, or VJP routine receives a copy of the named
`equation_objective_block_t`; this lets a client dispatch its own profile,
force, diagnostic, or constraint formula without putting that formula in
FortFEM.

The value and JVP paths return one packed real vector in registry order.  The
VJP path scatters packed cotangents through the callbacks and sums their state
cotangents.  Callback failures, shape mismatches, and non-finite data leave an
invalid status and no partially accepted value vector.  The test uses direct
formulas and an independent finite-difference/transpose oracle.

This is a neutral adapter for DESC-like direct-force clients and other
equilibrium applications.  It does not select a profile law, coordinate
normalization, optimizer, plasma model, or file format.
