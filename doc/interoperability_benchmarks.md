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

## Versioned oracle manifest

Every external result may be accompanied by the neutral
`fortfem-oracle-manifest-1` text record.  The public
`oracle_manifest_t` contract records the external code name, release and
immutable revision, the applicable license, case revision, coordinate system,
SHA-256 checksums for coordinates and sampled data, declared units and scale
factors, comparison tolerances, runner identity, wall-clock phase timings,
peak memory, repetition counts, and an optional immutable sister-repository
URI.  `write_oracle_manifest` and `read_oracle_manifest` provide a bounded
round-trip format; `validate_oracle_manifest` rejects incomplete provenance,
non-finite scales/tolerances/timings, invalid dimensions, and inconsistent
performance metadata before data are published.

The code-name field is intentionally open rather than an enum.  Thus an
adapter can identify CHEASE, FreeGS, VMEC/PARVMEC, GVEC, DESC, SPEC/SPECTRE,
GPEC, MARS-F, GLISS, STARWALL, JOREK, FreeFEM, MFEM, or FEniCSx without making
FortFEM depend on that project's format or license.  The manifest is metadata
only: sampled fields remain caller-owned artifacts in the separate benchmark
data repository.

### Target-ladder validation fixture

`test/test_oracle_manifest_ladder.f90` is a deterministic, metadata-only
fixture for the external target ladder. It instantiates the same small
manufactured case for CHEASE, FreeGS, VMEC/PARVMEC, GVEC, DESC, SPEC/SIESTA,
GPEC, MARS-F, GLISS, STARWALL, JOREK, FreeFEM, MFEM, and FEniCSx. For every
name it requires a nonempty release, immutable revision, declared license,
case and sample checksums, normalization, tolerances, runner identity, and
phase-level timing and memory fields. It also checks that the declared total
time covers all measured phases and that a manifest cannot be published after
its revision is removed.

The fixture uses synthetic strings such as `fixture-v1` and
`runner-supplies-revision`; they are deliberately not claims about an
upstream release or license. No external executable, reader, source tree,
mesh, or solution data is needed to run it. Real adapters replace those
fields with values obtained from their isolated runner and publish samples in
the separate benchmark-data repository.

The ordinary correctness workflow also runs a lightweight isogeometric oracle
with Nutils 9.2. On the same uniform quadratic tensor-product patch, Nutils
independently assembles the scalar mass and stiffness matrices and verifies:

- 20 global spline degrees of freedom;
- unit constant-field mass;
- unit energy for the reproduced coordinate field; and
- annihilation of constants by the stiffness matrix.

Nutils is installed ephemerally with `uv`; no Nutils source, wheel, or generated
matrix is committed. FortFEM retains separate analytical tests on nonuniform
knot vectors, which Nutils' structured-topology parameterization does not use
as an identical algebraic basis.

The published result records, for every mesh level:

- code name, exact version, compiler, scalar type, and solver configuration;
- cells, degrees of freedom, and polynomial order;
- \(L^2\), energy, and (for Ampere) \(H(\mathrm{curl})\) errors;
- assembly, setup, solve, and total wall-clock time;
- peak resident memory and a description of the runner hardware;
- a fixed probe-grid trace used for solution-difference plots.

Adapters should use FortFEM's neutral `interchange_sample_set_t` contract for
that probe grid. It validates the common physical coordinates, positive
weights, component shape, and provenance before reporting weighted absolute,
relative, and (L^\infty) differences. This keeps comparisons between
different bases or mesh numberings independent of any plasma-specific reader.

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

- [Nutils documentation](https://docs.nutils.org/)
- [DOLFINx installation](https://docs.fenicsproject.org/dolfinx/main/python/installation.html)
- [FreeFEM installation](https://doc.freefem.org/introduction/installation.html)
- [MFEM build documentation](https://mfem.org/building/)

Every published benchmark record additionally links the exact upstream release
and the FortFEM commit that produced it.
