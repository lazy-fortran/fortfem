---
title: Interface normal-traction balance
---

# Interface normal-traction balance

`fortfem_interface_traction_balance` provides the neutral scalar interface
residual

\[
r = \mathbf n\cdot(\mathbf t^+ - \mathbf t^-) - t_{\rm target},
\qquad \lVert\mathbf n\rVert=1.
\]

The traction vectors may come from the generated CGL tensor, an anisotropic
elastic stress, a Maxwell stress, or any caller-owned constitutive block. The
target is an application-owned scalar with explicitly declared units; FortFEM
does not select a plasma pressure law, wall law, or equilibrium closure.

```fortran
call assemble_normal_traction_jump( &
    plus_traction, minus_traction, normal, target, residual, status)
call assemble_normal_traction_jump_jvp( &
    plus_traction, minus_traction, normal, target, plus_dot, minus_dot, &
    normal_dot, target_dot, residual_dot, status)
call assemble_normal_traction_jump_vjp( &
    plus_traction, minus_traction, normal, target, residual_bar, plus_bar, &
    minus_bar, normal_bar, target_bar, status)
```

The fixed-topology JVP includes traction, normal, and target increments. The
VJP is the exact transpose of that product rule. The normal is validated as a
finite unit vector, and no derivative is taken through a topology or
orientation rebuild.

`test_interface_traction_balance` checks the scalar contraction, product-rule
JVP, real dot-product VJP, and rejection of a non-unit normal. It is an
independent interface oracle rather than a check that merely mirrors the
implementation loops.
