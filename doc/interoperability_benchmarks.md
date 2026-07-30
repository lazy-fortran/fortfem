---
title: FEM interoperability benchmarks
---

# Interoperability benchmarks

FortFEM's interoperability suite compares independent implementations of the
same mathematical problem. It is deliberately separate from the unit-test and
documentation jobs: a missing third-party package must not hide a FortFEM
failure, and installing several complete FEM stacks must not slow every commit.

## What is compared

The initial cross-code suite covers two manufactured problems on the same
unit-square meshes:

| Case | Space | Exact field | Operator |
|---|---|---|---|
| Poisson | \(H^1_0\) | \(u=\sin(\pi x)\sin(\pi y)\) | \(-\Delta u=f\) |
| Ampere | \(H_0(\mathrm{curl})\) | \(\mathbf E=(\sin(\pi y),\sin(\pi x))\) | \(\mathrm{curl}\,\mathrm{curl}\,\mathbf E+\mathbf E=\mathbf f\) |

Each implementation uses the same geometry, mesh levels, polynomial order,
essential boundary condition, quadrature policy, and manufactured right-hand
side. The analytical fields are the primary oracle. Cross-code agreement is a
secondary oracle that can expose a shared FortFEM implementation error but
never replaces convergence against the exact solution.

The published result records, for every mesh level:

- code name, exact version, compiler, scalar type, and solver configuration;
- cells, degrees of freedom, and polynomial order;
- \(L^2\), energy, and (for Ampere) \(H(\mathrm{curl})\) errors;
- assembly, setup, solve, and total wall-clock time;
- peak resident memory and a description of the runner hardware;
- a fixed probe-grid trace used for solution-difference plots.

Correctness plots show errors and observed convergence rates. Performance plots
are separate and never imply that timings from different hardware are
comparable.

## Repositories and generated data

FortFEM contains only:

- the neutral case definitions and independently written runner programs;
- small validation fixtures;
- the schema and scripts that validate and plot benchmark records;
- documentation describing versions, commands, and provenance.

Third-party source code, containers, meshes, raw logs, solution dumps, and plot
images are not vendored. CI artifacts are the staging area. Accepted benchmark
records belong in a dedicated benchmark-data repository under
`benchmark/fortfem/<suite-version>/<runner-id>/`; the documentation workflow
must check out an explicitly pinned commit from that repository.

There is currently no `lazy-fortran/data` or `lazy-fortran/benchmark-data`
repository. Until one is created, GitHub Actions artifacts are retained for
review and feed the generated Pages gallery, but are not treated as durable
published measurements.

## Runner policy

- **FEniCSx:** use the official DOLFINx container. A plain `uv` environment is
  not sufficient because the Python package requires the compiled DOLFINx C++
  core. `uv` may manage only lightweight reporting dependencies outside that
  container.
- **FreeFEM:** use a pinned FreeFEM release/container or distribution package in
  its own job.
- **MFEM:** build a pinned release in a cached, isolated job, or use an official
  release container when available.
- **FortFEM:** use a release build from the exact commit under test.

All runner versions are immutable in a published record. Scheduled benchmark
jobs may propose a version update, but they do not silently rewrite historical
data.

The runner programs are original adapter code licensed with FortFEM. Merely
executing an external GPL/LGPL/BSD solver does not copy that solver into this
repository. Container images and binary packages retain their upstream
licenses; their names, versions, source links, and licenses are recorded in the
benchmark metadata.

## CI and publication

The ordinary `test` workflow runs FortFEM tests and examples only.
Interoperability runs are manual or scheduled and use one isolated job per FEM
code. A final validation job rejects missing fields, non-finite values,
inconsistent mesh identifiers, failed convergence thresholds, or an unexpected
solver version.

After review, a benchmark-data commit is immutable. GitHub Pages checks out the
pinned commit, validates it again, and generates all tables and Fortplot figures
during the docs build. Plot images remain generated artifacts and are never
checked into FortFEM.

While the durable data repository is unavailable, the interoperability
workflow validates the three raw artifacts, generates four Fortplot figures
directly from those records, and publishes the uncommitted figures to Pages.
The gallery contains separate Poisson and Ampere accuracy plots and separate
runner-local timing diagnostics. Its timing plots are explicitly not a
cross-code performance ranking.

For meaningful performance numbers:

1. use a dedicated runner label with documented CPU, RAM, OS, governor, and
   thread counts;
2. build optimized binaries and disable interactive visualization;
3. perform one warm-up followed by at least five measured repetitions;
4. publish median and dispersion, while retaining every raw repetition;
5. separate mesh generation, assembly, preconditioner/factorization, solve, and
   output;
6. compare codes only on the same runner and numerical contract.

GitHub-hosted runners remain useful for correctness and coarse regression
detection, not for authoritative speed rankings.

## Provenance

Installation choices follow the upstream documentation:

- [DOLFINx installation](https://docs.fenicsproject.org/dolfinx/main/python/installation.html)
- [FreeFEM installation](https://doc.freefem.org/introduction/installation.html)
- [MFEM build documentation](https://mfem.org/building/)

Every published benchmark record additionally links the exact upstream release
and the FortFEM commit that produced it.
