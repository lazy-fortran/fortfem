# FEniCSx interoperability runner

`manufactured_oracles.py` independently solves the same unit-square Poisson and
Ampere cases as FortFEM and FreeFEM. It uses triangular P1 and first-family
lowest-order Nédélec elements, respectively.

DOLFINx requires its compiled C++ core, PETSc, Basix, and FFCx. Consequently the
CI job uses the official versioned `dolfinx/dolfinx` image. A standalone `uv`
environment cannot provide that compiled stack; `uv` is therefore reserved for
lightweight post-processing dependencies.

The runner emits `poisson.csv` and `ampere.csv`, verifies monotone analytical
convergence, and never writes plot images into the repository.
