# MFEM interoperability runner

This original adapter uses MFEM's triangular `H1_FECollection` and
`ND_FECollection` to solve the shared Poisson and Ampere cases. The Ampere
field and forcing also appear in MFEM's upstream Example 3, providing direct
literature/software provenance for the independent formulation.

The workflow builds MFEM 4.9 commit
`d9d6526cc1749980a2ba1da16e2c1ca1e07d82ec` outside the repository, compiles
the adapter through MFEM's configured example build, and uploads only CSV
records and version metadata. MFEM itself is BSD-3-Clause licensed and is not
vendored.
