---
title: Building
---

# Building FortFEM

## fpm

The Fortran library, tests and examples are built with
[fpm](https://github.com/fortran-lang/fpm). It needs a Fortran 2018 compiler,
BLAS/LAPACK and OpenMP, and fetches fortnum, fortsparse and fortplot at the
revisions pinned in `fpm.toml`.

```bash
fpm build
fpm test
fpm run --example simple_poisson
```

## CMake (C/C++ API)

CMake builds the shared C/C++ library `fortfem::capi` and installs a CMake
package. It fetches the pinned fortnum and fortsparse revisions when they are
not already available as targets.

```bash
cmake -S . -B build/cmake -DCMAKE_BUILD_TYPE=Release
cmake --build build/cmake
cmake --install build/cmake --prefix /desired/prefix
```

Downstream projects then use

```cmake
find_package(fortfem CONFIG REQUIRED)
target_link_libraries(my_app PRIVATE fortfem::capi)
```

The installed header [`fortfem.h`](https://github.com/lazy-fortran/fortfem/blob/main/include/fortfem.h)
documents the zero-based array layouts, process-local handles and status
arguments.

To build and run the sparse adapter, derivative and C API tests:

```bash
cmake -S . -B build/cmake -DFORTFEM_BUILD_TESTING=ON
cmake --build build/cmake -j2
ctest --test-dir build/cmake -L sparse --output-on-failure -j2
```

### Options

| Option | Default | Effect |
| --- | --- | --- |
| `FORTFEM_FETCH_DEPENDENCIES` | `ON` | Fetch pinned dependencies when the targets are not already available |
| `FORTFEM_BUILD_TESTING` | `OFF` | Build the C API and sparse adapter tests |
| `FORTFEM_NATIVE_SPARSE_ONLY` | `ON` | Disable SuperLU, UMFPACK, MUMPS and CUDA in the fetched fortsparse |
| `FORTFEM_TRIANGLE_COMPAT_ONLY` | `OFF` | Build only the dependency-free `fortfem::triangle_compat` mesher library |
| `FORTFEM_ELEMENT_KERNELS_ONLY` | `OFF` | Build only the dependency-free batched element and CSR interface |
| `FORTFEM_ENABLE_CUDA` | `OFF` | Build batched CUDA element kernels |
| `FORTFEM_BUILD_BENCHMARKS` | `OFF` | Build the element scaling benchmark |
| `FORTFEM_ENABLE_ENZYME` | `OFF` | Enable optional Flang/Enzyme element-kernel comparisons |
| `FORTFEM_ENZYME_REQUIRED` | `OFF` | Fail if the requested Enzyme toolchain is unavailable |

## Sparse direct solves

Sparse direct solves always select fortsparse's native CPU backend, for real
and complex matrices, including the retained transpose/adjoint factors used by
the derivative paths. There is no fallback to another sparse solver. When AMD or
METIS is available, supernodal LU is used; fortsparse chooses nested dissection
for the default ordering and AMD when a memory cap or low-memory mode is
requested. Builds without an ordering library use the scalar path.

The one-shot and retained factor wrappers accept a fortsparse
`sparse_options_t` as an optional final `options` argument, to select
scalar/natural ordering, AMD or nested dissection, or a memory cap. An
unsupported combination returns the fortsparse status unchanged.
