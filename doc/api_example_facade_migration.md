# API-05 representative example imports

The API-05 gallery migration is staged by capability rather than by file
count.  The following small set is the fast import gate for the public
facades.  Each migrated program is still a physical or manufactured solution
example; no numerical kernel was copied or rewritten.

| Capability | Representative program | Canonical import | Status |
| --- | --- | --- | --- |
| Scalar Poisson | `example/simple_poisson/simple_poisson.f90` | `fortfem_api` | compatibility smoke path |
| Vector/Nedelec | `example/tetra_nedelec_p_convergence.f90` | `fortfem_feec` | migrated |
| Exterior Laplace BEM | `example/laplace_exterior_bem_circle/laplace_exterior_bem_circle.f90` | `fortfem_boundary` | migrated |
| Circular Helmholtz DtN | `example/circular_dtn_modes.f90` | `fortfem_boundary` | migrated |
| Mixed acoustic wave | `example/mixed_acoustic_wave/mixed_acoustic_wave.f90` | `fortfem_time` | migrated |
| Three-dimensional mixed wave | `example/mixed_wave_3d_structure/mixed_wave_3d_structure.f90` | `fortfem_time` | migrated |

The Poisson program is intentionally the one umbrella-import compatibility
smoke.  It keeps a complete legacy client available while the smaller
facades grow to cover the general FEM form and mesh contracts.  The other
five programs no longer import `fortfem_api`; their canonical modules
re-export the existing implementations directly.

Run the structural gate from the repository root:

```text
python3 tools/check_example_facade_imports.py
```

The gate resolves each expected module from its Fortran `module` declaration
and checks the program's `use` associations.  It deliberately does not treat
that path check as a numerical test.  The independent behavioral oracle is
the short `test_api05_example_facades` target, which checks constant Nedelec
interpolation, an exact logarithmic BEM integral, a circular DtN eigenmode,
and midpoint wave energy/reversibility:

```text
fo test test_api05_example_facades
bash test/test_example_facade_imports.sh
```

Gallery images and build products remain generated outputs and are not part of
this migration commit.
