# Mixed acoustic midpoint wave

This fixture is a small physical wave problem for the common mixed
first-order time-step API. Two independent modal pairs satisfy

\[
\dot q=-C^T v,\qquad \dot v=Cq,
\]

with unit mass matrices and modal frequencies one and 1.7 radians per
second. The implicit-midpoint step is the Cayley map of this
port-Hamiltonian system. It preserves the quadratic energy and is exactly
reversible under a signed time-step reversal.

The first plot shows the propagated modal displacement-like and
velocity-like variables. The phase-space and energy plots are secondary
structure diagnostics. The program compares the numerical trajectory with
the independent closed-form oscillator solution and records the wall time in
`benchmark.txt`.

This is a generic acoustic/wave foundation fixture. Pressure, displacement,
electromagnetic, and elasticity clients can supply their own compatible mass
and interconnection blocks without changing the time integrator.
