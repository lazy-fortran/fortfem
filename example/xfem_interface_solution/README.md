# XFEM interface solution

This physical-first gallery fixture constructs a manufactured scalar and
vector field with a diagonal internal interface using the matrix-level
shifted-Heaviside XFEM contracts. It shows the reconstructed scalar solution,
the vector field with arrows, and the enriched jump contribution before any
diagnostic data.

The example is a space-construction foundation, not a production interface
solver. Cut-cell quadrature, global FEEC numbering, Piola maps, interface
laws, and sparse assembly remain explicit caller-owned layers. The executable
checks an independently written scalar and vector sign oracle and records the
construction time.
