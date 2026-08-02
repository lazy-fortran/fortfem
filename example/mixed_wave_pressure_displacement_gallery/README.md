# Mixed pressure--displacement wave state

This neutral gallery advances a first-order mixed wave state with the existing
structure-preserving midpoint map.  The coordinate block is intentionally
labelled as one pressure amplitude and two displacement amplitudes; the second
block is labelled as the conjugate momentum and two velocities.  It is not a
second-order displacement-only update.

The first three outputs are physical solutions: a three-dimensional
pressure/displacement trajectory, a two-dimensional pressure contour with
displacement arrows, and one-dimensional pressure, displacement, momentum, and
velocity traces.  Energy, exact-state error, reversibility, and symplectic-map
defects are emitted afterward as diagnostics.  Diagonal positive masses make
the manufactured modal solution independently computable without any plasma or
material-specific model.

The corresponding test derives the modal solution and Cayley map independently
of the example implementation.  This fixture is therefore a reusable
structure-preserving time-integration contract for later FEM/IGA wave,
elasticity, acoustics, or anisotropic systems.
