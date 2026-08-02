# Physical Maxwell boundary parity

`test_maxwell_physical_boundary_parity` supplies a small, license-safe oracle
for the common physical-sample layer used by FEM, BEM, DtN, and PML clients.
It evaluates an analytical complex vector trace on torus/cylinder-labelled
fixed-topology samples, applies independent backend perturbations, and checks
the weighted error and trace-space metadata.  It does not assemble a Maxwell
solver or select a boundary condition; production kernels remain caller-owned.
