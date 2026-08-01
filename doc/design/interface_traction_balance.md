---
title: Interface traction balance
---

# Interface traction balance

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

## Full vector traction jump

For a law that constrains all traction components, the neutral vector residual
is

\[
\mathbf r=\mathbf t^+-\mathbf t^- -\mathbf t_{\rm target}.
\]

`assemble_traction_jump` and its JVP/VJP companions expose this residual
without assuming whether the traction is CGL pressure, Maxwell stress, or
elastic stress. The target vector can therefore represent a prescribed
surface force, a multiplier contribution, or a caller-owned jump law.

```fortran
call assemble_traction_jump( &
    plus_traction, minus_traction, target, residual, status)
call assemble_traction_jump_jvp( &
    plus_traction, minus_traction, target, plus_dot, minus_dot, target_dot, &
    residual_dot, status)
call assemble_traction_jump_vjp( &
    plus_traction, minus_traction, target, residual_bar, plus_bar, &
    minus_bar, target_bar, status)
```

`test_interface_vector_traction_balance` checks the componentwise analytical
oracle, the product-rule JVP, the real dot-product VJP, and shape rejection.
Normal projection remains a separate scalar contract so callers can choose
normal, tangential, or full-vector interface laws explicitly.
