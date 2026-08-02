# FortFEM roadmap

Status: living architecture and verification plan, 2026-08-02

FortFEM is a Fortran library for finite-element, boundary-element, and
compatible discretizations. The long-term goal is to provide the reusable
mathematical ingredients from which three-dimensional equilibrium, linear
stability, and nonlinear resistive-MHD applications can be built. This file
is a roadmap for those ingredients. It is not a promise to reimplement
VMEC, GVEC, SPEC, GPEC, MARS, GLISS, JOREK, CHEASE, FreeGS, or any other
production code in this repository.

The target architecture is:

\[
\boxed{\text{FEEC}+\text{IGA}+\text{Fourier-FEM}+\text{DG/HDG}
      +\text{XFEM/XIGA}+\text{explicit interfaces}
      +\text{distributional sources}+\text{FortSym}
      +\text{differentiable residuals}+\text{structure preservation}}
\]

The public result should be a small, composable library with many toy
problems, analytical solutions, and independent oracles. A future equilibrium
or MHD application can then select the appropriate spaces, geometry, weak
forms, solvers, and time integrator without changing the numerical contracts.

## FortFEM and application boundaries

FortFEM owns numerical foundations. It does not implement plasma applications.
In particular, the repository will not contain GEQDSK or COCOS readers and
writers, equilibrium profile models, coil or vessel physics, Braginskii or
sheath closures, species and reaction networks, neutral or impurity models,
or production CHEASE, FreeGS, VMEC, GVEC, SPEC, GPEC, MARS, GLISS, JOREK,
GRILLIX, BOUT++, or PARALLAX algorithms.

External applications may provide those data and laws through documented
callbacks or adapters. FortFEM supplies the typed geometry, field, trace,
tensor, Fourier, residual, derivative, solver, and structure-preserving
contracts that consume them. License-safe examples use manufactured fields,
small neutral arrays, or data produced by an external adapter. A reference
code name in this roadmap identifies a parity target or oracle, not a planned
FortFEM module.

### Foundation-only delivery rule

Every equilibrium or MHD item below is delivered as a reusable ingredient:
space, geometry map, trace or interface operator, residual block, derivative
action, solver contract, manufactured fixture, and provenance record. A gallery
example may compose these pieces into a small analytical or supplied-data
demonstration. It must not grow into a production equilibrium, stability, or
time-dependent MHD code. Application-specific profiles, closures, readers,
continuation policy, and physical interpretation stay in the external client.
The acceptance question is whether an independent client can assemble and
verify the needed operator, not whether FortFEM can run a named code end to
end.

## How to read this document

The labels below describe the state of a deliverable.

| Label | Meaning |
| --- | --- |
| **done** | Present on `main` and covered by an independent behavioral oracle |
| **active** | Partly present or currently being closed with tests first |
| **planned** | Required architecture or example that has not been implemented |
| **reference** | External code or literature target, not a FortFEM implementation |

The detailed current audit is in [the FEEC, MEPHIT, and open-boundary
roadmap](doc/roadmap_mephit_feec_bem.md). Derivative contracts are specified in
[the differentiation design](doc/design/differentiation.md), and the existing
JOREK-related coverage is tracked in
[the JOREK operator document](doc/jorek_iga_operator_coverage.md).

## 1. Development principles

These principles apply to every phase. They are part of the product design,
not only contributor etiquette.

### 1.1 Mathematical and scientific principles

1. **Structure preservation is a first-class requirement.** Exact sequences,
   orientations, conservation laws, gauge structure, symplectic or Poisson
   structure, and correct dissipation must be represented explicitly. A result
   that looks plausible but violates the relevant invariant is not accepted as
   a reference implementation.
2. **Residuals come first.** Every formulation is an explicit residual
   \(R(u,p)=0\), with declared fields, spaces, traces, constraints, and
   boundary terms. The assembled matrix is a product of this residual, not a
   separate hand-written equation.
3. **FortSym is the source of truth wherever algebra is repetitive.** Weak-form
   reduction, tensor contractions, element kernels, enrichment derivatives,
   manufactured forcing, and analytical comparison data should be derived and
   generated. Generated files are revision-pinned, checked byte for byte, and
   never edited by hand.
4. **Independent oracles are mandatory.** A test must compare against an
   analytical solution, a separately written reference expression, a discrete
   identity, a conservation law, a cross-formulation result, or an external
   code. Checking that a file matches its own generator is a gate, not a
   behavioral test.
5. **Differentiation is part of the API.** New operators expose primal, JVP,
   and VJP actions with documented real and complex inner-product conventions.
   Shape, coefficient, geometry, and interface parameters are differentiated
   when the topology is fixed. Topology changes report piecewise-smooth or
   invalid derivative regions instead of silently returning a misleading
   gradient.
6. **Compatible spaces before stabilizing patches.** Use de Rham complexes,
   commuting projections, Piola maps, and orientation-aware traces before
   adding divergence cleaning or penalty terms. Stabilization must have a
   declared consistency and conservation effect.
7. **The code remains Fortran-first.** Do not move performance-critical work
   to C. Existing interoperability ABIs may remain where they are already
   required, but new numerical kernels are Fortran and generated Fortran.
8. **No Lean dependency is required.** FortSym identities, independent
   numerical oracles, exact-sequence checks, compiler diagnostics, and Enzyme
   checks are the normal proof and verification stack. A formal proof tool can
   be introduced later for a narrowly defined theorem if it provides a clear
   benefit, but it is not a prerequisite for the roadmap.

### 1.2 Engineering principles

- Make the smallest complete change. Add the test and oracle before the
  implementation.
- Keep focused tests short enough for roughly a ten-second feedback loop.
  `fo test` therefore keeps a ten-second per-target wall-clock limit. Tests
  that intentionally exercise large tetrahedral systems, high-order
  convergence, or curved torus operators carry a `_slow` suffix and are run
  explicitly with `fo test --all`; they are not silently allowed to consume the
  fast-suite budget. Run full `fo` suites and expensive examples asynchronously
  in CI or a controlled background job. The build-test workflow uses the same
  `fo test` runner for its per-test timeout and coverage path, rather than an
  unbounded serial `fpm test`. Do not wait for a complete gallery build before
  making an independent, low-risk progress increment.
- Work directly on FortFEM `main` because it is the upstream library. Work on
  PR branches for downstream repositories such as MEPHIT and FortPlot. Never
  duplicate a fix that has already landed upstream. Rebase downstream work on
  the current upstream branch before opening or updating a PR.
- Commit coherent increments, push them promptly, and record the exact source
  revisions, compiler, flags, hardware, and external-code versions for every
  benchmark.
- Keep generated plots and large benchmark output out of the source tree.
  GitHub Actions regenerates the gallery. Large or license-sensitive oracle
  data belongs in the separate benchmark-data repository, with a small
  manifest and immutable revision recorded here.
- Prefer pure, elemental, allocatable-safe Fortran. Avoid hidden global state,
  pointer ownership, aliasing, and allocation inside element kernels.
- Profile before optimizing. Use Fortran and existing `fortnum` or
  `fortsparse` facilities first. Do not add a new dependency merely to make a
  benchmark table look better.
- Treat performance as an acceptance criterion. Generated kernels should be
  allocation-free in their hot path, vectorizable, batchable, and measurable
  with a roofline or operation-count explanation. FortFEM targets competitive
  or better performance for each named workload. A claim such as "faster than
  the competition" is valid only for a declared problem, compiler, hardware,
  thread count, and reference implementation.
- Keep the implementation serial-first and distributed-ready. Global IDs,
  ownership, ghost data, deterministic reductions, local and global index
  maps, and communicator operations are part of the data model before a full
  MPI backend is added. Residual and element kernels remain independent of the
  communicator. This preserves a short single-node feedback loop while
  avoiding a later rewrite of mesh, trace, and sparse interfaces.
- Add a deterministic property-testing path to `fo` and the test helpers.
  Every randomized case records its seed, generated topology, and tolerance.
  Seeded random tests provide fast coverage of orientations, geometry
  degeneracies, mode sets, constitutive tensors, and solver dimensions. They
  complement analytical and manufactured oracles and never replace them.

### 1.3 Public API naming and modularity

FortFEM has no downstream consumers whose source compatibility constrains the
pre-release API. A naming cleanup may therefore be breaking and should happen
in one coordinated change. Do not add compatibility aliases merely to retain
an internal name. Update source, tests, examples, documentation, and generated
API references together, then remove the obsolete spelling.

The canonical procedure vocabulary is:

- `build_*` constructs topology, metadata, maps, or spaces;
- `evaluate_*` evaluates pointwise maps or analytical quantities;
- `assemble_*` forms residuals, matrices, or dense blocks;
- `apply_*` performs a matrix-free or factor action;
- `solve_*` performs a complete solve;
- `advance_*` performs one time step;
- `validate_*` checks metadata or input contracts;
- `compare_*` or `diagnose_*` reports parity, error, or invariant results.

`compute_*` is reserved for an explicitly named reduction or derived-array
operation when `evaluate_*` would obscure that it performs an aggregate
calculation. New names use the verb first and the numerical object second,
with a domain qualifier when it prevents ambiguity. Derivative actions use
the fixed suffixes `_jvp`, `_vjp`, and `_hvp`. Public derived types use the
`*_t` suffix, and metadata contracts use `*_metadata_t`.

Names such as `evaluate_*_parity` are scheduled for conversion to
`compare_*` or `diagnose_*` before the first external release. The same rule
applies to any procedure whose verb does not describe its side effect or
mathematical role. The API review must prefer established numerical terms
(`trace`, `residual`, `quadrature`, `skeleton`, `DtN`, `NtD`, `JVP`, `VJP`)
over project-specific synonyms.

The repository remains a modular monorepo. Public facades are split by
dependency layer so a client can import only what it needs:

1. `fortfem_core`: kinds, status, topology, geometry, spaces, and residual
   contracts;
2. `fortfem_feec`: FEM, FEEC, DG, HDG, XFEM/XIGA, and IGA operators;
3. `fortfem_fourier`: Fourier-FEM, toroidal modes, and special-function
   adapters;
4. `fortfem_boundary`: BEM, DtN, NtD, PML, wall, and free-boundary ports;
5. `fortfem_time`: mixed-wave, symplectic, dissipative, and structure ledgers;
6. `fortfem_interop`: external-code schemas, comparison metrics, and
   provenance;
7. `fortfem_plot`: plotting and gallery helpers.

The umbrella `fortfem_api` may remain as a convenience import while these
facades are introduced, but it is not the canonical dependency boundary.
Each facade must depend only on lower layers, keep generated kernels behind a
stable module interface, and have a focused test target. Splitting build
targets and modules comes before splitting repositories.

### 1.4 API reorganization and rename workstream

The current public surface has grown faster than its namespace: the
3,412-line `src/fortfem_api.f90` umbrella imports the older
`fortfem_api_types`, `fortfem_api_mesh`, `fortfem_api_spaces`,
`fortfem_api_forms`, `fortfem_api_solvers`, and plotting modules as well as
the newer domain modules. This is useful for examples, but it makes ownership,
dependency direction, generated-code visibility, and future semantic versioning
hard to audit. The following workstream is a product task, not cosmetic
renaming. It must be completed before the public API is treated as stable.
The current deterministic public-symbol inventory is documented in
[`doc/api_public_inventory.md`](doc/api_public_inventory.md), its generator is
[`tools/generate_api_public_inventory.py`](tools/generate_api_public_inventory.py),
and the first bounded rename candidates are recorded in
[`doc/api_rename_candidates.md`](doc/api_rename_candidates.md).

The migration is deliberately staged:

1. inventory every exported symbol and its defining module, source file,
   derivative actions, tests, examples, and generator provenance;
2. freeze the dependency layers and create canonical facades
   (`fortfem_core`, `fortfem_feec`, `fortfem_fourier`, `fortfem_boundary`,
   `fortfem_time`, `fortfem_interop`, and `fortfem_plot`);
3. move or re-export symbols behind those facades without changing numerical
   behavior, and make `fortfem_api` a thin convenience import;
4. apply the verb-first names in §1.3 in one pre-release rename, updating all
   source, tests, examples, documentation, and generated API references in
   the same commit; and
5. remove obsolete spellings after the migration gate passes. A temporary
   compatibility alias is allowed only when an external release has already
   depended on the name; it must carry an explicit removal issue and test.

The refactor must not create a second implementation of an operator. A public
facade owns validation and metadata; a generated module owns kernels; a
domain module owns the residual or algebra; and the umbrella only re-exports
the canonical symbol. No generated file may import the umbrella, and no
facade may depend on a higher layer or on a plotting module. All public
procedures retain primal/JVP/VJP/HVP conventions and fixed-topology
differentiability during the move.

#### Agent-ready API task ledger

Each row is independently reviewable and must be done on an agent branch or
worktree with disjoint file ownership. Agents must not mass-rename symbols
until API-00 has produced the inventory. A row closes only with its focused
test, an independent behavioral oracle, and a short migration note.

| ID | Scope and deliverable | Acceptance gate |
| --- | --- | --- |
| **API-00** | Build a deterministic public-symbol inventory from `fortfem_api.f90` and all existing API modules. Record symbol, defining module/file, proposed facade, current spelling, proposed spelling (if any), generator, derivative actions, tests, examples, and status. Detect duplicate exports, orphaned public symbols, cyclic facade dependencies, and generated modules that leak directly into clients. Do not rename code. | Re-running the inventory is byte-stable; every current export has exactly one owner; an independent fixture catches a deliberately duplicated or orphaned entry. |
| **API-01** | Write and enforce the module-layer graph: core → FEEC/forms/operators → boundary/time → solvers/interop/plot, with generated kernels below their owning facade. Add a small cycle/leak check that uses compiler-visible module dependencies rather than string matching alone. | The existing tree passes; a fixture with a reversed edge or umbrella import fails with a useful diagnostic; no numerical source changes. |
| **API-02** | Introduce the canonical facades in small slices, starting with `fortfem_core` and `fortfem_feec`, then the remaining facades. Re-export existing implementations first; keep `fortfem_api` as a thin umbrella. Move implementation ownership only after a facade has a focused test. | Public smoke programs compile using each facade and the umbrella; old and new import paths produce identical values, JVPs, VJPs, and structure ledgers on an analytical case. |
| **API-03** | Perform the coordinated verb-first rename. Prioritize ambiguous `evaluate_*_parity` and other non-semantic `compute_*` names, then types and constructors. Update call sites, examples, docs, generated manifests, and FortSym registration atomically. Remove obsolete spellings unless the release policy explicitly records a compatibility alias. | A compile-time old-spelling scan finds no stale internal use; independent value/JVP/VJP and conservation or symplectic oracles pass for every renamed operator; the generated-code currentness check remains green. |
| **API-04** | Audit generated and handwritten boundaries. Give every generated family a stable internal module and one checked public wrapper; add visibility tests so clients cannot accidentally import an unpinned generated symbol. Record FortSym revision and regeneration command in the inventory. | Regeneration is byte-current; a negative compile fixture cannot use an internal generated symbol; wrapper and direct-kernel outputs agree independently. |
| **API-05** | Migrate examples and documentation in increasing complexity. Every example imports the smallest suitable facade, while one smoke gallery retains the umbrella path. Generate API pages and cross-links from the inventory; do not check in plots or build products. | `fo test`, formatting, example compilation, link checks, and the example documentation coverage test pass; at least one physical solution plot per migrated example remains unchanged within the declared image/data oracle. |
| **API-06** | Add the release gate and deprecation/removal policy. The gate checks public exports, layer cycles, stale names, generated currentness, focused derivative tests, and a clean module dependency graph in the fast ten-second path; slow cross-compiler and gallery jobs remain asynchronous. | A clean tree passes the gate; each failure mode has an independent fixture and actionable diagnostic; the roadmap records the release milestone at which obsolete names disappear. |
| **API-07** | Run a downstream-style consumer audit without importing any application physics. Compile minimal clients for FEM, IGA/Fourier, FEM/BEM/DtN/PML, mixed waves/elasticity, and external-adapter interchange using only canonical facades. | Clients compile with no umbrella dependency, preserve units and metadata, and reproduce the same analytical/oracle results as the pre-refactor baseline. |

The first implementation increment is intentionally small: API-00 and API-01
may land without adding facades, and API-02 may land with only re-exports.
This keeps the product usable while agents parallelize the inventory and
refactor. A rename is not considered complete because `rg` found a new
spelling; it is complete only when the behavioral, derivative, generated, and
documentation gates above pass.

#### API migration status on `main` (2026-08-02)

- **API-00 — complete first slice:** `doc/api_public_inventory.md` and its
  standard-library-only generator inventory 1,761 unique exports, wrapper
  ownership, facade classification, derivative siblings, call-site evidence,
  and generated-code provenance. The inventory now includes the first
  canonical facade and is regenerated byte-current in CI.
- **API-01 — complete first slice:**
  `scripts/check_module_layers.py` audits compiler-visible `module`/`use`
  edges, duplicate providers, cycles, umbrella/facade/generated leaks, and
  has independent negative fixtures in `test/test_module_layer_audit.sh`.
- **API-02 — seven canonical slices complete:** `fortfem_core` directly
  re-exports foundational cell-complex and toroidal-coordinate contracts;
  `fortfem_feec` exposes exact-sequence/commuting-projection and tree--cotree
  gauges; `fortfem_boundary` exposes planar Helmholtz DtN and boundary-port
  metadata; `fortfem_fourier` exposes mode registry/expansion and toroidal
  harmonics; `fortfem_time` exposes mixed-wave/symplectic/dissipative steps;
  and `fortfem_interop` exposes sample, oracle-manifest, and boundary
  comparison contracts; `fortfem_plot` exposes the existing mesh, sampling,
  and solution-plot contracts without generated media. The FEEC facade also
  exposes field-aligned constitutive tensors and 3-D tensor-diffusion
  contractions used by the anisotropic gallery. Each has a direct-import
  facade and focused analytical or structure smoke test. No facade may grow
  into the umbrella.
- **API-03 — six parity renames complete:** the boundary parity family is now
  `compare_boundary_operator_parity{,_jvp,_vjp}` throughout its defining
  module, umbrella exports, tests, and design documentation. The obsolete
  internal spelling is gone; the independent weighted complex value/JVP/VJP
  and metadata oracles pass. The larger-domain family is now
  `compare_larger_domain_solution{,_jvp}` with its report schema preserved and
  independent weighted/JVP oracle. The sheet-current and surface-sheet-current
  representation families are now
  `compare_sheet_current_representations` and
  `compare_sheet_current_surface_representations{,_jvp}`, with independent
  Gaussian/surface-quadrature ledgers and derivative checks. The Beltrami
  families are now `compare_beltrami_two_region_residual` and
  `compare_beltrami_shell_residual`, with independent curl-eigen,
  constraint, and energy ledgers and unchanged report schemas. The
  tree--cotree diagnostic is now `diagnose_tree_cotree_iga_invariance` with
  signed-map and loop-period invariance checks.
- **API-04 — generated visibility complete:**
  `scripts/check_generated_visibility.py` and
  `test/test_generated_visibility.sh` enforce private generated implementation
  modules, explicit stable wrapper exports, and negative client/facade
  fixtures. The gate is standard-library-only and passes on the current tree;
  the temporary analytical-oracle allowlist is explicit and documented by the
  checker.
- **API-05 — first migration slices complete:** nineteen representative examples
  now import the smallest suitable `fortfem_feec`, `fortfem_boundary`, or
  `fortfem_time` facade, with one scalar Poisson umbrella compatibility smoke.
  The added boundary slice covers 2-D Helmholtz/Laplace BEM spectra, CFIE,
  symmetric FEM--BEM transmission, the neutral free-boundary trace gallery,
  the 3-D tree--cotree solution gallery, the field-aligned anisotropic
  3-D solution, the mixed pressure--displacement wave gallery, and the 3-D
  CGL tensor-pressure gallery through canonical capability surfaces.
  `tools/check_example_facade_imports.py` and
  `test/test_api05_example_facades.f90` provide import and behavioral gates;
  the remaining gallery still needs a complete smallest-facade audit and a
  solution-first visual oracle for every page.
  The 2026-08-03 migration wave extended this work without changing numerical
  implementations: scalar/core and anisotropic/FCI consumers landed in
  `999d08b`, six boundary BEM/DtN/PML consumers in `0e8fb7d`, six IGA/Fourier/
  toroidal consumers plus a dedicated consumer oracle in `900f7cd`, and seven
  wave/tensor/elasticity/response consumers in `03b511b`. The representative
  checker now covers 25 cases (one deliberate umbrella smoke), while the
  broader migrated examples are guarded by their focused behavioral tests and
  module-layer audit. The remaining gallery still needs the same smallest-
  facade treatment; these commits do not close API-05 globally.
- **API-06 — complete first gate:**
  `scripts/check_api_release_gate.py` and
  `test/test_api_release_gate.sh` compose byte-current inventory, module-layer,
  generated-visibility, stale-name fixtures, and focused boundary/larger-domain
  derivative tests under a ten-second timeout. Full gallery and cross-compiler
  jobs remain asynchronous; the deprecation/removal release policy is still
  active work.
- **API-07 — first client complete:** the no-umbrella downstream smoke client
  now imports `fortfem_core`, `fortfem_fourier`, `fortfem_boundary`,
  `fortfem_time`, and `fortfem_interop` exclusively (apart from kinds and
  sparse utility modules). It covers core geometry, Fourier modes, planar
  Helmholtz DtN, structure-preserving mixed-wave stepping, and interchange
  metrics with independent analytical, energy, reversibility, and weighted
  error oracles. Additional consumer clients remain useful as the facades
  grow, but the first API-07 gate is closed.

#### Parallel agent queue (API reorganization)

The following queue is part of the roadmap and is intentionally split into
disjoint worktree-owned increments. Agents rebase from the current `main`,
write tests before implementation, run the focused ten-second gate, commit a
coherent slice, and report the hash. No agent edits another slice's files or
waits for the full gallery/CI job before handing off.

| Queue item | Owned files and result | Independent acceptance oracle |
| --- | --- | --- |
| API-03-SHEET | `fortfem_sheet_current_parity` definition, umbrella/facade exports, sheet-current tests/docs; canonical `compare_sheet_current_representations` export and value parity are complete | Fitted surface ledger versus independently integrated regularized layer; old-spelling internal scan |
| API-03-SURFACE | Surface-quadrature definition, exports, tests/docs; canonical `compare_sheet_current_surface_representations{,_jvp}` and interop-facade value/JVP exports are complete | Orientation/measure/toroidal quadrature oracle, central-difference JVP, and invalid-measure rejection |
| API-03-BELTRAMI | **complete**: renamed `evaluate_beltrami_two_region_parity` and `evaluate_beltrami_shell_parity` to the documented `compare_*_residual` names, keeping flux/helicity/energy schemas unchanged; canonical FEEC facade exports both reports and comparators | Independent curl-eigen residual, constraint, and energy value oracles; the underlying residual JVP/VJP closure remains covered; no physics closure or reader changes |
| API-03-TREE | `diagnose_tree_cotree_iga_invariance` is complete in the implementation, umbrella, FEEC facade, tests, and inventory | Signed-map invariance, loop-period, direct-reduction, and fixed-topology oracle |
| API-06-GATE | Release-gate runner, isolated negative fixtures, and fpm-safe fixture layout are complete; compose inventory, layer, generated-visibility, stale-name, and canonical-consumer checks | Clean-tree pass plus one independently failing fixture per gate component |
| API-05-GALLERY | First ten-example import audit and documentation migration are complete, including the 2-D boundary/BEM slice; continue in scalar-to-vector-to-interface order with no media committed | Example compilation, link/coverage checks, and solution-first numerical/visual data oracle |
| GALLERY-BIRO-TEAM | **first slice complete**: `biro_tree_cotree_benchmark` and the 3-D `biro_tree_cotree_3d_gallery` provide provenance-pinned tree--cotree curl--curl reductions, direct reconstruction, face-circulation/energy/residual CSV, and solution-first graph/cube vector plots; `team7_neutral_benchmark` and `team13_neutral_benchmark` provide license-safe TEAM-shaped manufactured fields, probe/energy diagnostics, and solution-first magnitude/vector/geometry plots. The remaining TEAM 3/13/20 ladder stays external-data work; benchmark arrays, readers, nonlinear material laws, and reference curves do not enter FortFEM | Independent tree/cotree restriction/prolongation and reduced-system oracles; finite TEAM field/probe/energy manifests; generated media remain ignored and no paper/TEAM source or reader is copied |
| GALLERY-BIRO-EXACT-DATA | Keep the method-faithful Bíró tree--cotree fixture in FortFEM, but add a provenance-gated adapter slot for the paper’s exact geometry/material/source arrays when the user supplies a redistributable or sister-repository dataset. Never label the manufactured cube as an exact paper figure and never copy copyrighted paper data into the library. | With an approved external manifest, the adapter reproduces the paper’s stated field/energy/error tables and plot metadata; without it, the neutral fixture remains the only checked-in oracle and reports its limitation explicitly |
| GALLERY-TEAM-LADDER | Extend the neutral TEAM gallery one problem at a time (TEAM-3/7/13/20) through external manifests and optional sister-data jobs. Each in-repo fixture remains tiny, analytic, license-safe, and solution-first; exact probes, B-H curves, and reference measurements stay outside FortFEM. | Every enabled problem has a checksum/license/revision manifest, independent field/probe/energy oracle, mesh-complete solution plot, and a bounded optional external runner; absent data is a documented skip, not a fabricated comparison |

| API-05-SOLUTION-FIRST | Migrate the remaining gallery in scalar-to-vector-to-interface order. For each example, assign the smallest canonical facade, make the first generated panel a physical solution (scalar field, vector arrows, surface, or field line), and retain convergence/diagnostic panels after it. Add a bounded shell gate that checks the primary plot and an independent CSV/JSON oracle without committing media. `circular_dtn_modes` now records its physical mode map before diagnostics and has a ten-second PNG/stage-order/discarded-mode oracle. | The example-doc coverage test reports every executable example; each migrated example runs under the ten-second gate, produces a non-empty primary data record, and its independent manufactured/analytical oracle passes |
| API-07-CONSUMER-MATRIX | Add disjoint downstream-style smoke clients for the growing product boundaries: core geometry/mesh, FEEC/IGA, Fourier/special functions, FEM--BEM/DtN/PML, mixed waves/elasticity, solver/gauge, time integration, plotting, and interop. Each client imports only its canonical facade and records the public symbols it intentionally uses. | Each client compiles without `fortfem_api`, reproduces an independent value/JVP/VJP or structure-ledger result, and the inventory has one owner plus one migration note for every imported symbol |
| API-07-MIXED-WAVE | Migrate the mixed pressure/velocity, displacement/momentum, and structure-preserving wave galleries to `fortfem_time` plus the smallest FEEC/core facades. Keep the umbrella path only in the compatibility smoke and record every moved symbol in the public inventory. | A no-umbrella consumer compiles, reproduces the independent energy/reversibility/symplectic-form oracle, and the first generated panel is the physical field/trajectory before diagnostics |
| API-07-TENSOR | Migrate tensor-pressure, stress, gyrotropic, and field-aligned anisotropy galleries to canonical FEEC/core facades. Keep constitutive laws caller-owned and do not duplicate generated kernels in examples. | The client compiles without `fortfem_api`, passes independent contraction/JVP/VJP and power-ledger checks, and emits a solution-first 2-D/3-D field plot |
| API-03-COMPATIBILITY | Audit the post-rename public surface for obsolete spellings, duplicate wrappers, and generated registration drift. Keep aliases only for names already used by an external release, with a removal milestone and a negative test that prevents new internal use. `tools/check_api_compatibility.py` now enforces a byte-stable allowlist for the 16 intentional boundary/interoperability metadata re-exports and has a negative duplicate-facade fixture wired into the release gate. | Release-gate stale-name and duplicate checks are clean; alias inventory is byte-stable; old and canonical entry points agree where an alias is explicitly retained; no generated module imports the umbrella |
| GALLERY-FREE-BOUNDARY | **first neutral slice complete**: `free_boundary_port_gallery` samples a manufactured toroidal trace, renders a connected 3-D surface with vector components first, and verifies the free-boundary-port value/JVP/VJP contract plus 2-D and 1-D diagnostics. It deliberately does not implement equilibrium, coil, BEM, DtN, or plasma physics; those remain caller-owned adapters and sister-data work | Independent trace residual, central-difference JVP, real dot-product VJP, non-empty 3-D surface/vector plot, CSV/JSON provenance, and ten-second execution gate |

The queue is extensible: a new public family gets an inventory row and a
disjoint owner before implementation. A rename may not be merged merely
because `rg` finds the new spelling; the derivative, generated provenance,
consumer, and documentation gates remain part of the same task.

#### Current agent dispatch (2026-08-03)

The following worktrees are the active API-05/API-07 consumer migration
queue. They are deliberately disjoint so that facade refactoring can proceed
without a serial gallery-wide rename or duplicated fixes. Each agent must
rebase its branch on the current `main` before hand-off, add or update a
focused behavioral test first, migrate the smallest suitable canonical facade
(`fortfem_core`, `fortfem_feec`, `fortfem_fourier`, `fortfem_boundary`,
`fortfem_time`, `fortfem_interop`, or `fortfem_plot`), run the ten-second
focused gate, and report the exact files and public symbols moved. Agents may
not edit this roadmap, generated API inventory, or another queue item's files;
the integrator regenerates the inventory and documentation after each verified
slice.

| Worktree task | Owned migration slice | Required gate |
| --- | --- | --- |
| `migrate_scalar_gallery_facades` | Simple Poisson, mesh, scalar/core, and anisotropic/FCI examples and their focused tests; split legacy `fortfem_api_*` imports between `fortfem_core`, `fortfem_feec`, and `fortfem_plot` without adding a second implementation | No-umbrella compile/import audit, independent analytical solution or field/energy oracle, and focused `fo test` |
| `migrate_boundary_gallery_facades` | Adaptive BEM, acoustic DtN/NtD, Helmholtz/PML, Maxwell open-boundary, and boundary solver examples/tests; own the representative facade checker updates | FEM/BEM/DtN/PML physical-solution oracle, derivative/trace parity, checker negative fixture, and focused `fo test` |
| `migrate_iga_fourier_facades` | IGA, Fourier-FEM, toroidal, special-function, and field-aligned examples/tests; expose missing symbols through `fortfem_feec` or `fortfem_fourier` only | Geometry/mode/toroidal analytical oracle, physical-first plot-data gate, no-umbrella compile, and focused `fo test` |
| `migrate_wave_tensor_facades` | Mixed acoustics/waves, elasticity, wall, tensor-pressure, sheet-current, and linear-response consumers; use `fortfem_time`, `fortfem_feec`, `fortfem_boundary`, `fortfem_interop`, and `fortfem_plot` as applicable | Independent energy/symplectic, constitutive, response, or trace oracle; no-umbrella compile, physical-first plot-data gate, and focused `fo test` |

All four 2026-08-03 dispatches have now handed off clean branches. Their
commits are recorded above and the integrator has regenerated the API
inventory at `0fabd83`; future agents must take the next disjoint example
batch rather than reopening these files or duplicating their exports.

This dispatch is a consumer migration, not permission to mass-rename public
procedures. Verb-first renames remain coordinated API-03 work: an agent may
use an already canonical name, but a new rename requires an inventory row,
derivative and generated-code evidence, a negative stale-name test, and one
integrator commit that updates all call sites atomically. A hand-off that
fails its independent oracle is returned to the owning worktree rather than
being papered over in `main`.

## 2. Current baseline

The following capabilities are already on FortFEM `main` or in the verified
documentation baseline. The list is intentionally conservative.

| Area | Current state | Next gate |
| --- | --- | --- |
| Scalar FEM | P1/P2/Q1 and arbitrary-order triangular scalar paths, Poisson and diffusion forms, boundary conditions, plotting | General field and coefficient callbacks in the symbolic form compiler |
| FEEC | Oriented triangular and tetrahedral H1, H(curl), H(div), and DG families, Piola maps, commuting tests, sparse assembly, mixed RT-DG Poisson, neutral real and complex two-field/packed N-field block graph residuals with complete JVP/VJP contracts, real/complex packed-graph CSC adapters, retained real/complex field-split solve JVP/VJP products, and retained coupled Schur value/JVP/VJP composition | FortSym manufactured composed forms, energy/power ledgers, and arbitrary multipatch assembly |
| IGA | Nonuniform B-splines, rational maps, two- and three-dimensional de Rham incidence complexes, cylindrical and toroidal Fourier blocks, initial JOREK magnetic-flux residual/JVP, and arbitrary 2D/3D tensor-patch H1/H(curl)/H(div) packed signed maps with orientation-cycle checks, face swaps, and reversals | General geometry-aware trace/mortar assembly, trimming, enrichment, and the remaining coupled JOREK variables |
| Special functions | FortNum quadrature, ordinary and associated Legendre P/Q, orthonormal complex spherical harmonics with angular derivatives and Gaunt products, Hobson-normalized toroidal P/Q branches, stable high-degree half-integer continuation, Bessel/Hankel paths, and a FortFEM Fourier mode registry with phase/radial derivative contracts; the pinned spherical and toroidal APIs are re-exported through `fortfem_api` and checked by integration oracles | Uniform asymptotics for arbitrarily large degree/order and cross-geometry special-function oracles |
| Sparse algebra | FortSparse CSC assembly, retained factors, real and complex solves, sparse products, tree--cotree CSC direct reductions with fixed-map JVP/VJP, and CG, PCG, GMRES, and BiCGSTAB converged-state derivative contracts; dense and standalone sparse IC(0)/ILU(0) factor/apply paths plus deterministic sparse fixed-factor ILUT, memory-scalable row-oriented ILUT and ICHOL, controlled ICHOL paths, and a solver-gallery timing fixture are public | Production-size measured scaling, flexible Krylov products, and block solver derivatives |
| Linear response interchange | Neutral modal `(m,n)` metadata, provenance, complex equilibrium/inertia/resistive/vacuum/wall blocks, response channels, the (K-omega^2M+iomega R+V+W) operator and forced residual with complex JVP/VJP products, reciprocity/passivity diagnostics, common real and complex weighted physical sample-set contracts with real-part comparison JVP/VJP actions, a weighted inner/outer singular-layer trace-matching block, a bounded versioned text schema round-trip, and deterministic bounded cross-factor actions are public and independently tested | FEM/BEM/DtN/PML assembly-specific response matrices, singular-layer asymptotic models, and external-code sampler fixtures |
| Open boundaries | Planar, circular, and spherical scalar Helmholtz DtN paths, scalar BEM, Maxwell trace and PML components, a mixed RWG/RBC weak Maxwell DtN map assembled from exact-curved torus EFIE/MFIE/mass forms, toroidal Laplace/Helmholtz off-surface reconstruction with optional target gradients and fundamental-solution oracles, fixed-geometry data/target JVP/VJP products, toroidal Laplace and Helmholtz representation geometry JVP/VJP products, assembled toroidal Laplace and Helmholtz DtN single-layer/double-layer/mass geometry JVP/VJP products covering regular and coincident panels, a manufactured toroidal Maxwell FEM--BEM solved-state/field-reconstruction fixture, a neutral complex implicit value/JVP/VJP solve map for concatenated volume/surface states, and a geometry-generated curvilinear tetrahedral Nédélec PML CSC chain with full layer/mesh JVP/VJP products | Assembly-specific curved-object matrix/solver derivatives, larger-domain toroidal Maxwell/PML comparisons, nonzero-scattering vector FEM/BEM/DtN parity, and robust vector field reconstruction fixtures |
| PML | Scalar and curl-curl Cartesian complex-stretching tensors with slab, triangular, and tetrahedral examples, plus full complex 3-by-3 curvilinear scalar and curl-curl coefficient maps, tetrahedral scalar P1 and Nédélec element/CSC/solved-state JVP/VJP products, caller-owned normal-frame curved-layer geometry, and weighted complex reflection/error diagnostics with diagonal-reduction oracles | Automated curved-object layers and derivative coverage for all geometry parameters |
| Differentiation | Analytical FortSym paths, selected Enzyme checks, sparse matrix products, converged CG/PCG/GMRES/BiCGSTAB solves, toroidal coordinate and DtN products | Complete operator inventory, JVP/VJP parity for all public operators, and shape derivatives |
| Parallel readiness | Serial local kernels and deterministic focused tests | Owned/ghost mesh and field data, partition-independent numbering, halo exchange, distributed assembly, checkpointing, and MPI-enabled solver backends |
| Examples | Generated documentation pages for Poisson, Maxwell, Helmholtz, BEM, IGA, torus, PML, solver, and structure-preserving mixed wave/wall examples; the new mixed wave--wall gallery starts with a physical state plot before energy diagnostics | Ordered gallery beginning with simple Poisson and adding 1D, 2D, and 3D application-oriented toy models with manufactured or external oracle data |

The aggregate full test suite can expose resource-sensitive failures when many
independent numerical targets run in parallel. Focused repetitions are the
first diagnostic and remain the merge gate for a local increment. CI remains
the final cross-platform gate. The current fast suite has 430 targets; 25
deliberately long regression targets are classified as `_slow` and retain a
separate explicit run path. A result formatter truncation is never treated as
a pass: batch or named runs are used when a report exceeds the display limit.

## 3. Stable contracts before physics applications

### 3.1 Field, space, and residual contracts

Every public PDE object must declare:

- domain dimension, coordinate chart, orientation, and physical units;
- trial, test, trace, and multiplier spaces;
- scalar, real-vector, complex-vector, or tensor value type;
- primal residual and active constraints;
- boundary and internal-interface terms;
- nullspace or gauge treatment;
- assembled and matrix-free actions;
- JVP and VJP actions, including parameter and geometry derivatives;
- independent manufactured forcing and an analytical or cross-code oracle.

The standard interface is conceptually:

```text
residual(state, parameter) -> R
jvp(state, parameter, state_dot, parameter_dot) -> R_dot
vjp(state, parameter, R_bar) -> state_bar, parameter_bar
solve(residual, constraints) -> converged_state
implicit_jvp(converged_state, ...) -> state_dot
implicit_vjp(converged_state, ...) -> state_bar, parameter_bar
```

The solve derivatives assume a converged state and hold iteration and stopping
branches inactive. This contract must be visible in each solver's documentation.
For complex fields, the VJP convention and conjugation must be stated next to
the API and checked by a dot-product identity.

### 3.2 Topology and geometry contracts

Meshes, spline patches, Fourier charts, and internal manifolds share a common
orientation model. Geometry objects provide values, Jacobians, inverse maps,
surface normals, measures, and derivatives with respect to control parameters.
The same physical problem can be represented by a fitted interface or by an
unfitted level set without changing the interface-law API.

### 3.3 Foundation gap register

The current algebraic slices are reusable building blocks. The following
contracts still have to be composed before FortFEM can serve as the basis for
an arbitrary-topology three-dimensional MHD or edge application.

| Foundation | Contract still missing | Required independent oracle |
| --- | --- | --- |
| Topological complex | A region and cell-complex graph with periodic identifications, orientations, homology, cohomology, harmonic representatives, cuts, and gauges; `cell_complex_homology_cycle_basis` and the dual `cell_complex_cohomology_cocycle_basis` now return deterministic quotient representatives of `ker(boundary_1)/im(boundary_2)` and `ker(boundary_2^T)/im(boundary_1^T)` alongside the metric-harmonic one-form period normalization map and complex edge-period residual JVP/VJP actions. `build_multipatch_signed_dof_map` now provides arbitrary-patch signed local-to-global numbering and rejects inconsistent orientation cycles; geometry, trace maps, and distributed ownership remain separate | Chain-complex identities, Euler characteristic, cycle and flux integrals, filled-loop quotient, normalized periods, and nullspace dimension on slab, cylinder, sphere, and torus cells |
| Sheet-current interface | Neutral open/closed internal-manifold graph, integrated-current junction ledger, fixed-topology loop-current constraints, differentiable geometry-to-edge-flux contraction, topology-only edge-flux balance, scalar and full-vector traction jumps, an independent test/trial tangential surface-current trace residual, and the caller-owned surface-quadrature parity ledger (`compare_sheet_current_surface_representations`) are public. The latter accepts fitted/cut/DG/IGA normals and positive measures, checks fixed row topology and unit orientation, and compares the Ampere trace against a resolved Gaussian layer on slab, cylinder, sphere, and torus-labelled quadratures; constitutive pressure laws and flux/helicity constraints remain | Ampere jump, surface-current conservation, loop current, pressure jump, orientation reversal, invalid-measure rejection, and regularized-layer limits |
| Cut FEEC spaces | Scalar and vector matrix-level shifted-Heaviside enriched-space constructors, a 3D vector-enrichment curl/divergence product-rule diagnostic, the physical vector/metric support Gram contraction, batched 2D/3D covariant/contravariant Piola-enrichment composition, matching 2D/3D affine Piola H(curl)/H(div) differential contracts with value/JVP/VJP actions, a rectangular commuting-projection audit, and signed dense/CSC local-to-global `assemble_glued_feec_sequence` reference compositions with value/JVP/VJP actions are public; the CSC path now also reports sparse `curl(grad)` and `div(curl)` products with product-rule JVP/VJP actions; Piola-aware vector-compatible XFEM/XIGA and DG spaces that preserve or explicitly report the de Rham sequence across cuts remain | Curl-gradient and divergence-curl identities on every generated space, Piola/XIGA constructors, and fitted versus unfitted convergence |
| Coupled field residuals | Neutral real and complex rectangular field-plus-constraint and packed N-field block graph residuals now compose caller-owned vector, tensor, interface, boundary, FEM, BEM, DtN, or PML blocks with complete JVP/VJP actions; real and complex semantic normal/tangential boundary-port residuals cover supplied trace/jump targets and work weights. The neutral tensor volume-work and tensor-weighted diffusion contractions are public. Complex paths use the real-part adjoint convention. The packed graph's scalar value, product-rule JVP, and reverse product are emitted by the pinned FortSym generator and checked against an independent matrix oracle. Real/complex packed-graph CSC adapters provide retained-factor storage without dense assembly, `retained_field_split` composes those fixed factors into real/complex block-diagonal solve JVP/VJP paths, and `assemble_retained_coupled_schur` now supplies the retained elimination/JVP/VJP layer. Plasma state assembly remains in an external client | FortSym manufactured residuals for composed physical forms, energy or power balance, cross-formulation parity, and composition of the retained Schur layer with arbitrary geometry/patch graphs |
| Equilibrium interchange | The neutral external-adapter schema for mapped coordinates, physical samples, named scalar/vector/tensor coefficients and scalar profiles, segmented boundaries, provenance, units, and normalization is public and validated. `build_equilibrium_interchange_sample_set` now projects selected, already-sampled physical coefficient components onto the common weighted sample-set schema, and weighted L2/relative comparison metrics expose fixed-coordinate value/weight JVP/VJP actions. GEQDSK and COCOS parsing remain outside FortFEM | Analytic manufactured data covers toroidal mapped samples, named fields, segmented boundaries, vector/tensor component offsets, copy semantics, rejection, selected-component common-grid projection, and central/adjoint comparison derivatives; license-safe CHEASE and FreeGS outputs still enter through external resamplers |
| Linear response interchange | The neutral external-adapter record now carries modal `(m,n)` metadata, complex frequency, provenance, normalization, response channels, and caller-owned equilibrium, inertia, resistive, vacuum, and wall blocks. The composed operator, forced residual, complex JVP/VJP, normalized reciprocity error, conservative Hermitian passivity lower bound, neutral generalized eigen-residual `K u - lambda M u` with analytic JVP/VJP, common real and complex weighted physical sample-set validation with real-part comparison JVP/VJP actions, independently sized weighted inner/outer singular-layer trace matching, a bounded versioned text schema round-trip, and deterministic bounded cross-factor actions are public; no GPEC/MARS/GLISS/STARWALL reader or closure is included | Assembly-specific FEM/BEM/DtN/PML response matrices, singular-layer asymptotic models, reciprocity/passivity normalization, and external sampler fixtures |
| Fourier and toroidal modes | The fixed-topology mode registry now provides field-period phase, normalization, conjugate packing, retained-triad lookup, an ordered triad map with omitted-sum count, deterministic one-product and bounded repeated-product work registries, radial regularity, complex coordinate derivatives, a same-registry three-factor vector/tensor Fourier product, a distinct-registry bilinear application with JVP/VJP, weighted modal-energy/conjugacy diagnostics with JVP/VJP, and a bounded-memory retained-mode toroidal Green convolution with complex JVP/VJP. The padded one-product registry now grows geometrically instead of reserving four N+N² arrays, reducing OOM risk for sparse/colliding mode windows. FortNum now exposes independently tested ordinary Legendre Q values/derivatives, orthonormal complex spherical harmonics with angular derivatives and Gaunt products, and Hobson-normalized half-integer toroidal P/Q values/derivatives. The neutral nested-surface geometry evaluator now composes caller-owned Fourier/radial coefficients for `(R,Z,lambda)` with `phi=zeta+lambda`, physical coordinates, both Jacobians, metric, volume, and fixed-topology coefficient JVP/VJP; equilibrium equations and model-specific branch data remain external | Nonlinear application composition with physical operators and uniform high-order torus-harmonic envelopes |
| Edge and SOL equations | Equation-as-data fields, generic coefficient and boundary callbacks, conservative sources, FCI events, vector-valued terminal ledgers, a vector-valued volume-source ledger with JVP/VJP actions, and a deterministic composed FCI-map fixture. Species and closures remain client-owned | Manufactured source terms, mass and energy balances, terminal flux tallies, and broader client-owned reproducible FCI maps |
| Mixed waves and elasticity | A common compatible port-Hamiltonian state for pressure, velocity, displacement, momentum, and tensor stress, including boundary power ports, plus separate dissipative Cayley and resistive-wall midpoint blocks; the coupled wave--wall midpoint now preserves port power and reports RL dissipation; a neutral mixed `C sigma - E u`, `D sigma - f` elasticity residual with complete JVP/VJP is public; the generic `assemble_symplectic_map_defect` diagnostic now reports \(S^T\Omega S-\Omega\) with JVP/VJP actions | Discrete energy, symplectic-form or passivity tests, dispersion, reversibility, and mixed versus second-order parity |
| Open boundaries | Curved vector FEM/BEM/DtN/PML coupling on toroidal external surfaces with larger-domain controls, plus a typed boundary-operator metadata/provenance contract shared by FEM, BEM, DtN, PML, NESTOR-like, BIEST-like, and virtual-casing clients | Reciprocity, passivity, far-field, reflection, and interior-field agreement across all four paths |
| Verification and delivery | Seeded random tests, external-code adapters, provenance manifests, mesh-completeness checks, and Pages health checks | Repeated seeds, license and revision records, independent samplers, HTTP link checks, and FortPlot image regression |

Closing a row requires a public API, a focused test, an independent oracle, a
FortSym or provenance entry, and a gallery fixture when the row affects a
visible numerical result. A current primitive does not imply that its coupled
application residual is complete.

The retained-coupled-Schur gap is now closed at the neutral solver-composition
layer: assemble_retained_coupled_schur eliminates caller-owned exterior
couplings against retained real or complex field factors and provides value,
JVP, and VJP actions without dense global assembly. The independent
test_retained_coupled_schur oracle covers diagonal elimination,
finite-difference derivatives, and real/complex adjoint identities. The
topology-only `build_multipatch_signed_dof_map` now supplies arbitrary patch
counts and signed cycle checks for the next composition layer. Geometry-aware
trace construction, distributed ownership, and FortSym-generated physical
residuals remain; no MHD state or closure is selected here.

## 4. Structure-preserving numerics

Structure preservation is measured, tested, and shown in the gallery.

### 4.1 Spatial structure

- Maintain the discrete de Rham chain
  \(H^1\xrightarrow{\nabla}H(\mathrm{curl})
  \xrightarrow{\nabla\times}H(\mathrm{div})
  \xrightarrow{\nabla\cdot}L^2\) whenever the formulation requires it.
- Generate and test incidence matrices independently of metric matrices.
- Check \(\nabla\times\nabla=0\) and
  \(\nabla\cdot\nabla\times=0\) algebraically, including multipatch and
  periodic identifications.
- Preserve normal magnetic flux in H(div), tangential electric or magnetic
  traces in H(curl), and the declared jump when a surface current is present.
- Use exact or compatible divergence control. A cleaning term is never a
  substitute for a space that has the required topology.
- Keep skew, symmetric, positive, and constraint blocks identifiable in block
  systems. These properties drive solver choice and preconditioning.

### 4.2 Time structure

The time integrator is selected from the continuous structure.

- Ideal or nondissipative Hamiltonian components use variational, symplectic,
  Poisson, or multisymplectic updates. Implicit midpoint, Cayley, splitting,
  and composition methods are appropriate building blocks.
- Resistive, viscous, and thermal terms are dissipative. Their discrete energy
  law must be monotone or have the declared balance. They are not described as
  symplectic merely because they are solved implicitly.
- Coupled ideal plus dissipative systems use a structure-aware split, such as
  a symmetric composition of reversible and dissipative substeps, with order
  and invariant tests.
- Time-dependent examples report energy, magnetic helicity, cross-helicity,
  mass, divergence, charge, and constraint defects as appropriate.
- Generic Runge--Kutta stepping is allowed for a deliberately non-geometric
  transport example, but it is not the default for a Hamiltonian or
  structure-preserving path.
- Spacetime FEM is reserved for a later phase. The present roadmap requires a
  sound space-time contract without implementing spacetime FEM now.

### 4.3 Structure tests

Each time or operator increment adds at least one of:

1. an exact algebraic identity;
2. a conservation or dissipation balance;
3. a reversibility or symplectic-form test;
4. a manufactured solution and convergence rate;
5. an independent energy or adjoint identity.

### 4.4 Mixed first-order waves and mechanics

Wave equations are not forced into a second-order pressure-only form when a
first-order mixed form exposes more structure. The reusable state contract is
of the form

\[
\partial_t q + C^T v = 0,\qquad
M_v\partial_t v - Cq = f,
\]

where `q` may be pressure, displacement, electric flux, or another potential,
`v` is its flux or velocity companion, and `C` is a compatible discrete
gradient, divergence, curl, or elasticity-complex map. The semidiscrete
Hamiltonian or port-Hamiltonian form must expose its skew interconnection and
positive mass or constitutive blocks.

The roadmap therefore includes:

- mixed first-order acoustics with pressure and particle velocity;
- mixed elastodynamics with displacement, velocity, and stress or momentum;
- Maxwell and elastic wave systems with compatible electric, magnetic,
  traction, and displacement traces;
- acoustic-elastic and electromagnetic-elastic interface coupling;
- symplectic Euler, implicit midpoint, Cayley, variational, and partitioned
  symplectic updates for the nondissipative part;
- dissipative viscosity, damping, conductivity, and absorbing layers as
  metriplectic or energy-decaying substeps;
- exact semidiscrete energy conservation or a declared discrete power balance,
  plus dispersion and phase-error diagnostics.

This contract applies to acoustics and to the general wave family. A
second-order reduction remains available as a compatibility path, but it must
be demonstrably equivalent to the mixed first-order residual and must not hide
the conserved pairing needed by a geometric time integrator.

### 4.5 Tensor-valued pressure, stress, and constitutive laws

Pressure is a scalar only for an isotropic closure. FortFEM must support a
tensor-valued pressure or stress contribution \(\boldsymbol\Pi(x,u,p)\), with

\[
\boldsymbol\sigma=-\boldsymbol\Pi+\boldsymbol\tau,
\qquad
\mathbf f_{\mathrm{int}}=\nabla\!\cdot\boldsymbol\sigma,
\]

where the symmetric, nonsymmetric, deviatoric, gyrotropic, or field-aligned
parts are declared by the constitutive model. The form compiler must preserve
tensor symmetry and produce stress, traction, JVP, VJP, and shape derivatives.
Required capabilities are:

- anisotropic elasticity and nearly incompressible mixed elasticity;
- Hellinger--Reissner stress-displacement and weak-symmetry formulations;
- tensor constitutive coefficients with field-aligned, gyrotropic, or
  caller-defined anisotropy;
- symmetric and nonsymmetric stress with explicit angular-momentum balance;
- traction and normal-flux traces on fitted, cut, and BEM interfaces;
- constitutive tensors with spatial, parameter, and magnetic-field dependence.

For a field-aligned constitutive oracle, the CGL-shaped tensor is a first test
case,

\[
\mathbf P=p_\perp I+(p_\parallel-p_\perp)\,\mathbf b\mathbf b^T,
\qquad \mathbf b=\mathbf B/|\mathbf B|,
\]

with optional caller-supplied gyrotropic terms. The parallel and perpendicular
coefficients, their force divergence, and their work pairing are separate
residual terms. FortFEM does not provide a plasma closure or a Braginskii
model.

The first constitutive slice is now public on `main`: FortSym generates the
six independent symmetric CGL components and their JVP/VJP, while the
`fortfem_cgl_pressure_tensor` wrapper validates the unit magnetic direction,
packs a full symmetric tensor, and combines full-matrix off-diagonal
cotangents. The generated tensor now feeds a compositional traction
`t=P n` with value/JVP/VJP, including the `P^T` normal cotangent; independent
tests cover the closed-form traction oracle, central differences, adjoint
identities, and invalid directions. Generic force-volume assembly, caller-
supplied correction tensors, and field-aligned assembly remain separate active
work. These blocks are not a claim that an anisotropic MHD application is
implemented.

The elasticity complex is treated as a structure-preserving extension of the
de Rham complex. The mixed weak-symmetry construction of
[Arnold, Falk, and Winther](https://arxiv.org/abs/math/0701506) is the initial
literature oracle. The tangential-displacement, normal-normal-stress (TDNNS)
work from the [JKU Linz group](https://www.numa.uni-linz.ac.at/Talks/abstract/154/)
and [Pechstein and Schöberl](https://doi.org/10.1002/nme.3319) is the
anisotropic and nearly incompressible implementation reference. The first
FortFEM examples are a manufactured anisotropic elasticity patch, a nearly
incompressible block, and a tensor-pressure wave that compares
displacement-only and mixed stress formulations.

## 5. Geometry, topology, and interfaces

### 5.1 Embedded geometry

Implement first-class geometry for:

- level sets \(\phi(x,p)=0\), with gradients, Hessians, Boolean operations,
  periodicity, and parameter derivatives;
- cut triangles, tetrahedra, spline cells, and Fourier cells;
- adaptive, moment-fitting, recursive, and exact polynomial cut-cell
  quadrature;
- closed or open curves and surfaces, intersections, junctions, nested
  interfaces, separatrices, X-points, and rational magnetic surfaces;
- fitted refinement, anisotropic refinement, h, p, and hp refinement,
  hierarchical and truncated hierarchical B-splines;
- shape derivatives with explicit diagnostics when a cut crosses a node or a
  topology changes.

An interface is an object rather than a boundary marker:

```text
interface
    geometry
    plus_domain, minus_domain, orientation
    volume_spaces, trace_spaces, skeleton_spaces
    jump_laws, constraints, surface_unknowns
    quadrature, shape_derivatives
```

### 5.2 Broken and trace spaces

Provide broken versions of H1, H(curl), H(div), and L2. Their trace API must
support scalar, normal, tangential, rotated-tangential, plus, minus, average,
and jump operators with one orientation convention. Skeleton spaces support
scalar traces, normal-flux traces, tangential-field traces, surface current,
pressure-jump multipliers, and interface displacement.

### 5.3 Interface laws and distributional sources

The weak-form compiler must represent both volume and surface terms:

\[
\int_\Omega f v\,dV + \int_\Gamma g\,\operatorname{tr}(v)\,dS.
\]

For electromagnetics, the surface-current law is

\[
\mathbf K=\mathbf n\times[\![\mathbf H]\!],
\qquad
\mathbf J_{\mathrm{sheet}}=\mathbf K\,\delta_\Gamma.
\]

The implementation should support duplicated degrees of freedom, explicit
surface-current unknowns, Lagrange multipliers, mortar coupling, Nitsche
coupling, hybridization, cohesive-like jump laws, and regularized-layer
comparisons. A delta source is not approximated numerically when an explicit
surface integral is available.

### 5.4 Arbitrary topology and sheet-current composition

The interface contract becomes an application-ready foundation only after the
geometry and spaces are connected through a topological region graph. The
graph must represent volume regions, oriented facets, boundary components,
periodic identifications, open edges, closed loops, and junctions. It assigns
stable IDs to cells, facets, traces, cycles, and cohomology representatives so
that a mesh rebuild does not silently change a gauge or a Fourier mode.

The graph owns the following structure-preserving data:

- boundary, coboundary, and metric matrices with separate incidence and Hodge
  factors;
- harmonic one- and two-form bases for multiply connected domains and their
  flux normalization;
- gauge and nullspace constraints for vector potentials, scalar potentials,
  and pressure-like multipliers;
- plus/minus trace ownership for every interface facet, including a single
  owner for cut or overlapping FCI terminal pieces;
- topology-event records when an interface intersects a vertex, an edge, a
  periodic seam, or another interface.

The first quotient-complex slice is now public on `main`: signed vertex,
edge, face, and volume representative maps are composed with integer chain
boundaries, and the output is checked with exact boundary-of-boundary
identities. The interval-to-circle test supplies the independent Betti-number
oracle. This is still metadata only; quotient coordinates, metric Hodge
operators, gauges, and application-owned periodic laws remain later layers.

The cell-complex cycle-space slice is also public: scale-aware real kernels of
the oriented edge boundary and its face-boundary transpose are available for
subsequent cuts, flux normalization, and gauge construction. It intentionally
does not label those cycles or cocycles as harmonic forms or physical fluxes
until metric Hodge data and application constraints are supplied.

The metric-harmonic one-form slice is public as well. Given a positive edge
metric it returns the closed and co-closed kernel. The fixed-topology
`normalize_harmonic_one_forms` slice now maps that basis to caller-selected
cycle periods or flux units with dense-solve JVP/VJP actions, while leaving
cycle labels, physical units, and gauge ownership to the caller. The
fixed-topology tree--cotree gauge slice supplies a spanning-forest split,
cotree restriction/prolongation, direct dense-system reduction, and fixed-map
JVP/VJP actions. Its selector is frozen for differentiation; graph rebuilds
are topology events.
For high-order FEEC and IGA, the tree is built on the lowest-order
mesh/control graph and higher-order moments are retained by an external DOF
map.

The direct-gauge design follows [Manges and Cendes
(1995)](https://doi.org/10.1109/20.376275), [Bíró, Preis, and Richter
(1996)](https://ieeexplore.ieee.org/document/558631), and the high-order
tree-gauge analysis in [the tree-gauge review](https://doi.org/10.3390/j5010004).
The mortared IGA extension is tracked against [Kapidani et al.
(2022)](https://doi.org/10.1016/j.cma.2022.114654) and its
[preprint](https://arxiv.org/abs/2110.15860). These are literature and oracle
references, not imported source code.

An electromagnetic sheet is a coupled unknown, not only a post-processed
delta source. The implementation must support equivalent formulations and
their conversion:

1. a jump in tangential \(\mathbf H\);
2. an explicit tangential \(\mathbf K\) trace space;
3. a distributional volume current \(\mathbf K\delta_\Gamma\);
4. a resolved resistive layer approaching the sheet limit.

Each formulation carries the same orientation, units, edge balance, and
current-loop normalization. On a closed interface the surface divergence and
global current constraints are explicit. On an open interface the edge flux is
owned by a boundary law. For MHD interfaces, the sheet law is coupled to normal
magnetic-flux continuity, tangential field jumps, total-pressure balance,
field-line or loop-flux constraints, and optional helicity constraints.

The cut or enriched H(curl) and H(div) spaces must expose whether an
enrichment preserves the discrete de Rham sequence. A scalar Heaviside basis
cannot be promoted to a compatible vector space by naming convention alone.
The compiler and tests therefore record the commuting diagram, the required
trace regularity, and the exact identity that is intentionally relaxed when a
physical jump is present.

## 6. Approximation families and discontinuities

FortFEM supports two complementary representations of every important
discontinuity.

| Mechanism | Use case | Required ingredients |
| --- | --- | --- |
| Fitted interface | Known rational surface, material boundary, shielding sheet, multi-region equilibrium | Duplicated traces, interface elements, mortar or Nitsche, multipliers, surface current, independent constraints |
| Unfitted XFEM/GFEM | Moving or emerging interface, island separatrix, geometry optimization | Level sets, cut quadrature, partition of unity, shifted enrichment, blending correction, conditioning diagnostics |
| XIGA, immersed, trimmed, CutIGA | Spline geometry with interfaces not aligned to knot lines | NURBS/B-spline enrichments, connected-component activation, trimmed-support stabilization, shape derivatives |
| DG | Shocks, steep layers, nonmatching cells, discontinuous material or pressure | SIPG, NIPG, IIPG, upwind, entropy-stable fluxes, broken FEEC spaces |
| HDG | Local elimination and interface unknowns | Skeleton traces, static condensation, compatible H(curl)/H(div) fluxes |

### 6.1 XFEM and XIGA library

The initial enrichment library includes Heaviside jumps, signed-distance and
absolute-value kinks, logarithmic and algebraic singularities, crack-tip and
corner modes, user-defined analytical enrichments, vector and tensor
enrichments, helical Fourier factors, and resonant radial forms near
\(m-nq=0\). Shifted enrichments
\(N_i(x)(\Psi(x)-\Psi(x_i))\) preserve interpolation where possible.

The corrected blending-element construction from
[Fries (2008)](https://doi.org/10.1002/nme.2259) is required. Enrichment
activation must include support-size thresholds, rank-revealing checks,
orthogonalization or pruning, diagonal scaling, local condensation, and
condition-number diagnostics.

XIGA must support H1, H(curl), H(div), L2, multipatch, and trimmed supports.
Connected components of a spline support receive independent branches when
the physical support is disconnected. Trimmed B-splines must also support the
small-support stabilization described in
[stable trimmed IGA work](https://www.sciencedirect.com/science/article/abs/pii/S0045782516308222).
Physical jump enrichment and trimmed-space stabilization remain separate
composable mechanisms.

### 6.2 DG and hybridization

The planned DG layers are scalar elliptic interior penalty, hyperbolic central
and upwind fluxes, Rusanov/HLL-style fluxes where appropriate, broken vector
forms with tangential and normal penalties, and HDG with global skeleton
unknowns. Mixed CG-DG and XFEM-enriched DG are valid combinations. Every flux
declares whether it preserves conservation, entropy, energy, or only
consistency.

## 7. IGA, Fourier-FEM, and toroidal geometry

### 7.1 IGA and patch graphs

Move from the current two-patch checks to an arbitrary patch graph with
periodic identifications, toroidal closure, orientation reversals,
extraordinary junctions, nonmatching knot vectors, mortar coupling, and local
refinement. Continuity is field-specific: smooth regions may use
\(C^{p-1}\), interfaces may use \(C^0\), repeated knots may encode derivative
jumps, and traces may be fully discontinuous.

Trimmed and cut IGA are needed for coils, conductors, vacuum vessels, islands,
moving rational surfaces, and non-parametric interfaces.

### 7.2 Fourier-FEM contract

Fourier-FEM combines a Fourier or helical expansion in periodic angles with
finite elements or splines in the remaining coordinates. The canonical
representation is documented as

\[
u(s,\theta,\zeta)=
\sum_{m,n}u_{mn}(s)\,e^{i(m\theta-nN_{\mathrm{fp}}\zeta)}.
\]

The implementation must define phase, field-period, real/conjugate packing,
mode ordering, normalization, and de-aliasing. It must support:

- radial FE or B-spline spaces with regular-axis and edge constraints;
- 2D poloidal FE plus toroidal Fourier modes;
- mode-diagonal linear operators;
- nonlinear convolution with explicit retained-triad policy;
- complex block systems and their real equivalent;
- matrix-free mode actions and adjoints;
- parameter derivatives of mode numbers, geometry, and Fourier coefficients;
- Fourier-mode plots and energy per mode.

The neutral nested-surface geometry contract evaluates caller-owned modal
coefficients for `(R,Z,lambda)` at `(rho,theta,zeta)`, maps them with the
explicit `phi=zeta+lambda` convention to Cartesian samples, and returns both
coordinate Jacobians, `J^T J`, and `det(J)`. Its coefficient JVP/VJP is fixed
topology: changing mode labels, radial powers, or field-period metadata is a
registry rebuild. The independent contract and manufactured oracle are
documented in [`nested_surface_geometry.md`](design/nested_surface_geometry.md).
This is geometry only; nested-equilibrium residuals, profiles, fluxes, and
free-boundary laws stay in external clients.

This is the common numerical ingredient for CHEASE-like axisymmetric models,
GLISS and CAS3D-like linear calculations, GPEC/MARS-like perturbations, and
JOREK-like reduced-MHD operators. It is not a claim that all those codes use
the same formulation.

### 7.3 Special-function and geometry scope

FortNum now provides the real ordinary Legendre second-kind contract
`legendre_q`/`legendre_q_derivative` for x > 1, in addition to Ferrers
associated `P_l^m`, and its independently tested Hobson-normalized
half-integer toroidal `P/Q` API. It also provides standard orthonormal complex
`spherical_harmonic` values and analytical theta/phi derivatives on the closed
polar interval (with derivative evaluation intentionally undefined at its
poles), plus a Gaunt product coefficient with an independent Wigner-3j
Racah-sum oracle. The [pinned FortNum API
revision](https://github.com/lazy-fortran/fortnum/blob/3920109/docs/api.md)
records the domains, normalization, derivative convention, DLMF provenance,
and the current continuation scope. FortNum's zero-order toroidal Q
branch now uses the DLMF zero-balanced continuation near x=1 and generates
associated orders from its analytical derivative; half-integer toroidal P uses
an upward degree recurrence and the recessive Q branch uses scaled Miller
continuation at ordinary aspect ratios, with independent degree-80 references
and recurrence checks. FortNum `3920109` additionally exposes analytical
second derivatives for ordinary Ferrers/second-kind Legendre and half-integer
toroidal P/Q branches, checked against independent finite differences and ODE
residuals. The remaining limitation is uniform asymptotics for
arbitrarily large degree/order and cross-geometry special-function oracles.
Analytical solutions and DtN
maps are tested on slabs, cylinders, spheres, and exact toroidal surfaces. The
FortFEM now re-exports the pinned toroidal P/Q values and derivatives through
`fortfem_toroidal_harmonics`, with an independent branch/derivative/domain
oracle. The torus API must state whether a function is a toroidal harmonic, a
Fourier mode in toroidal coordinates, or a numerical continuation of a special
function. These are not interchangeable names.

## 8. PDE and boundary-operator ingredients

The same residual and derivative contracts cover scalar and vector-valued
equations.

### 8.1 Core equations

- Poisson, diffusion-reaction, and generic axisymmetric elliptic/Fourier
  operators that an external Grad-Shafranov client can instantiate;
- scalar Helmholtz with FEM, BEM, DtN, and PML boundaries;
- Ampere, curl-curl, eddy-current, Maxwell, and anisotropic H(curl) forms;
- H(div) flux, mixed Poisson, magnetic induction, and divergence constraints;
- generic linearized field, constitutive, interface, and wall-response blocks
  that an external ideal or resistive MHD client can compose;
- generic skew brackets, energy functionals, and compatible time operators for
  client-supplied state variables.

Pressure-like and stress-like quantities are represented as scalar, vector, or
tensor fields supplied by the caller. A tensor-valued coefficient is not
silently projected to its trace. Its symmetry, divergence, traction, and work
pairing are tested with manufactured data. A CGL-shaped tensor is one generated
constitutive oracle, not a plasma-closure implementation.

### 8.2 Strongly anisotropic and field-aligned operators

Plasma transport and wave models often contain coefficients whose parallel
and perpendicular scales differ by many orders of magnitude. FortFEM must
represent this explicitly rather than hiding it in an isotropic scalar
diffusion coefficient. For a unit magnetic direction
\(\mathbf b=\mathbf B/|\mathbf B|\), the basic decomposition is

\[
\mathbf K
 = k_\parallel\,\mathbf b\mathbf b^T
 + k_\perp\,(I-\mathbf b\mathbf b^T)
 + \mathbf K_{\mathrm{Hall}},
\]

with optional nonsymmetric Hall or gyroviscous terms. The same idea applies to
viscosity, thermal conduction, resistivity, pressure, elasticity, and wave
impedance. Required work includes:

- covariant and contravariant tensor pullbacks on curvilinear meshes;
- field-line and flux-coordinate-independent (FCI) parallel derivatives;
- perpendicular plane operators and parallel line operators with separate
  quadrature and tolerances;
- robust assembly for \(k_\parallel/k_\perp\) ratios spanning many decades;
- anisotropy-aware preconditioners, field-split Krylov methods, and fGMRES;
- exact energy or dissipation tests for symmetric and skew tensor parts;
- B-field, tensor, and geometry JVP/VJP actions;
- limiting tests as \(k_\parallel/k_\perp\to 1\), \(k_\perp\to0\), and
  \(k_\parallel\to\infty\).

The parallel operator may be represented by a compatible H(curl)/H(div)
complex, a Fourier derivative, or an FCI field-line map. These choices share a
residual and oracle interface but are not assumed algebraically identical.
The generic pointwise field-aligned flux

\[
  \mathbf F=k_\perp\mathbf g+(k_\parallel-k_\perp)\mathbf b
  (\mathbf b\cdot\mathbf g)
\]

is now public with FortSym-generated value/JVP/VJP products and independent
unit-direction, finite-difference, and dot-product tests. It is the common
constitutive block for future anisotropic diffusion, conduction, resistivity,
and wave assemblies.
The companion `evaluate_field_aligned_constitutive_tensor` adds the same
parallel/perpendicular projectors and an optional skew Hall cross-product
term. FortSym generates the nonduplicating three-scalar Hall product and its
JVP/VJP while the wrapper owns only the fixed skew layout; strict finite and
unit-direction validation remains independent of any plasma closure.
The neutral `evaluate_tensor_power_split` contract now reports symmetric,
skew, and total pointwise tensor power, with exact real JVP/VJP actions. Its
independent oracle makes the zero instantaneous power of Hall/gyrotropic skew
parts explicit while retaining tensor-valued derivatives; quadrature,
field-aligned pullbacks, and physical closures remain caller-owned.
The fixed 3-D value, JVP, and VJP contractions now delegate to pinned
FortSym-generated kernels; `check_generated.sh` pins their byte-current source,
while the existing deterministic and seeded constitutive tests remain the
independent behavioral oracle.
The first PARALLAX-aligned algebraic slice is now on `main`: a dependency-light
RK4 field-line tracer provides geometry endpoints; mapped upper and lower plane
interpolation matrices assemble a sparse staggered gradient; the support
divergence is constructed as its negative volume-weighted adjoint; and a
matrix-free (P K_\parallel Q) diffusion action is public. The matrix-free
gradient action exposes FortSym-generated value, JVP, and VJP kernels.
Independent tests cover a constant/exponential tracing oracle, a linear-field
map oracle, constants, an explicit flux-balance vector, a central-difference
JVP, the weighted adjoint identity, and the weighted negative-energy identity.
The `fci_parallel_diffusion` gallery example now runs the same matrix-free
action on a manufactured open-line cosine profile and publishes 1D FortPlot
profiles plus CSV values for the mass-rate and dissipation oracles.
The field-only VJP of this diffusion action is also public and is checked by
an independent dot-product identity. The complete fixed-topology diffusion
JVP/VJP now covers interpolation maps, line lengths, parallel coefficients,
canonical and staggered volumes, and the field through pinned FortSym local
contribution kernels; a central-difference and full real dot-product oracle
guard this contract. The public 1D linear interpolation-map builder now checks
partition of unity, affine reproduction, fixed-topology JVP/VJP dot products,
and Cartesian bilinear affine reproduction. A generated quadratic Lagrange
map now accepts explicit three-node stencils and reproduces quadratic fields
on nonuniform slices, with generated fixed-stencil JVP/VJP dot-product and
finite-difference oracles. Higher-order interpolation Jacobians beyond sextic
and degree-six-and-higher curved support-volume measures, and
anisotropy-aware preconditioning
remain active work. A generated fixed-topology quadrilateral area map now supplies
positive unstructured plane-cell measures with value/JVP/VJP actions, independent
shoelace/finite-difference/adjoint oracles, and a FortPlot gallery fixture. A
generic fixed-topology polygon map now accepts boundary-ordered cells with
any number of vertices at least three, with a generated per-edge contribution,
shared-vertex JVP/VJP accumulation, independent pentagon oracle, and gallery
fixture. Its quadratic Bezier-edge extension now supplies arbitrary curved
polygon value/JVP/VJP measures with independent Gauss--Green and gallery
oracles. The cubic Bezier-edge extension now supplies two-control-point curved
polygon value/JVP/VJP measures, generated by FortSym and checked against an
independent degree-five Gauss--Green oracle, central differences, and the real
VJP dot-product identity. The quartic Bezier-edge extension now supplies
three-control-point curved polygon value/JVP/VJP measures with an independent
degree-seven Gauss--Green oracle. Degree-seven-and-higher curved measures and
moving-cell connectivity remain active. The quintic Bezier-edge extension now
supplies four-control-point value/JVP/VJP measures with an independent
degree-nine Gauss--Green oracle. The sextic Bezier-edge extension now supplies
five-control-point value/JVP/VJP measures with an independent degree-eleven
Gauss--Green oracle. A batched Cartesian bilinear adapter now
turns traced forward/backward endpoint arrays into the per-segment FCI map
tensors used by the support operator, with source-grid accumulation in its
fixed-topology JVP/VJP and independent finite-difference/dot-product oracles.
The composed `test_fci_reproducible_map` fixture now rebuilds the same
deterministic endpoints bit-for-bit, checks partition of unity, applies the
sparse gradient to a manufactured affine field-line state, and verifies the
volume-weighted support-divergence adjoint identity.
An unstructured fixed-cell triangle adapter now supplies barycentric maps and
geometry/target JVP/VJP products with affine and dot-product oracles; higher
order triangle and moving-cell connectivity remain active. Its batched
forward/backward endpoint adapter now emits the per-segment tensors consumed
by the support operator and accumulates shared-vertex cotangents.
The positive diagonal of `-W_c^{-1}Q^TW_sK_\parallel Q` is now public as a
FortSym-generated per-stencil Jacobi baseline, with an explicit Q-squared
oracle and a validated matrix-free diagonal apply; plane multigrid and
field-split preconditioners remain active.
The split action now also accepts one validated square CSC perpendicular block
per canonical plane and combines it with the matrix-free FCI parallel term,
with a negative-energy anisotropic oracle. Its combined positive diagonal is
also public: the generated support diagonal is combined with the negative
diagonal of each oriented plane block and checked against an independent
oracle, giving a reproducible scalar preconditioner baseline while plane
multigrid and stronger field splitting remain active work. A convenience
matrix-free anisotropic Jacobi apply now uses the same combined diagonal;
iterative callers may cache the diagonal for repeated solves. A two-level plane
V(1,1) cycle now supplies CSC restriction/prolongation and a replaceable direct
coarse solve; deeper hierarchies and retained coarse factors remain active.

### 8.3 Equilibrium, perturbation, and edge equation data

The named application families share interchangeable data and residual
contracts. FortFEM supplies those contracts and small fixtures. It does not
ship their production solvers.

#### 2D equilibrium data

An external equilibrium adapter may provide a poloidal mesh or spline,
coefficient fields, axis or boundary metadata, coil or wall traces, units, and
normalization. FortFEM consumes this neutral data through a documented schema.
It does not parse or write GEQDSK, interpret COCOS conventions, or implement
equilibrium profiles. A generated manufactured field is the FortFEM baseline.
CHEASE and FreeGS can provide optional numerical oracle data through an
external, license-safe adapter sampled on the same physical points.

#### Linear ideal and resistive response data

The linear schema must distinguish an eigenproblem from a forced response and
must carry complex frequency, Fourier mode set, equilibrium coefficients,
vacuum and wall regions, interface traces, singular-layer matching data,
normalization, and response-matrix conventions. The residual exposes the
energy-principle, inertia, resistive, wall, and vacuum blocks separately. A
response matrix has reciprocity, passivity, and normalization tests. GPEC,
MARS-F, GLISS, and STARWALL remain versioned oracle or interchange targets.

#### Nonlinear MHD composition

The reusable nonlinear state schema must admit generic scalar, vector, and
tensor fields, constitutive coefficients, interfaces, and constraint
multipliers. An external application can map density, velocity, magnetic field,
pressure, current, or other variables into that schema. FortFEM provides the
residual composition, energy or power ledger, surface-current coupling, and
derivative actions. Plasma closures, profile laws, wall physics, and state
selection remain outside the library. Every optional numerical block has an
independent manufactured case before an application combines it.

#### Edge and SOL equation data

An edge application can map its fields and coefficients into generic equation,
source, conservative-flux, boundary-event, and material-ledger callbacks. The
core owns FCI geometry, anisotropic operators, trace spaces, and residual/JVP/
VJP composition. The vector-valued FCI terminal ledger accumulates multiple
conserved channels per canonical cell and also exposes the integrated event
tally, with analytical JVP/VJP actions for weights, fluxes, and cell volumes.
The matching volume ledger integrates caller-supplied source rates against
positive cell measures and exposes cell contributions, global channel totals,
and analytical JVP/VJP actions. Together these ledgers let a client compare
volume sources and terminal fluxes without FortFEM assuming a species model.
Species, Braginskii coefficients, sheath or Bohm laws,
recycling, radiation, neutral or impurity sources, and target material data
remain client-owned. The callback ABI records units and signs so that a
terminal ledger can be compared with the volume balance. This is the
foundation needed by GRILLIX, BOUT++, and PARALLAX-style clients without
copying or implementing their models.

### 8.4 FEM/BEM, DtN, and PML

The open-boundary layer must support scalar and vector equations in 2D and 3D,
including toroidal external surfaces:

- Laplace and Helmholtz single, double, hypersingular, Calderon, CFIE, and
  symmetric transmission blocks;
- Maxwell EFIE, MFIE, CFIE, RWG, BC/RBC, and higher-order trace pairings;
- planar, circular, spherical, and exact-curved-torus DtN and NtD maps;
- FEM/BEM and FEM/DtN coupling with consistent work-conjugate traces;
- Cartesian and curvilinear PML for scalar Helmholtz and curl-curl Maxwell;
- analytical Mie, cylindrical, toroidal, and manufactured field oracles;
- far-field, flux, reciprocity, passivity, and reflection diagnostics.

Every boundary path has an explicit larger-domain comparison: move the
artificial boundary away, compare FEM/BEM, DtN, and PML fields on a common
interior region, and report the expected nonzero but convergent difference.

### 8.5 Free-boundary, vacuum, wall, and response foundations

Free-boundary support is a numerical boundary-operator capability, not an
equilibrium solver. FortFEM must make the following composition possible for
an external VMEC, VMEC++, GVEC, DESC, SPEC/SPECTRE, CHEASE, FreeGS, GPEC,
MARS, GLISS, JOREK, or other client without importing that client's physics,
file formats, or conventions. The client owns profiles, coils, pressure laws,
force balance, state selection, and continuation; FortFEM owns geometry,
traces, operators, residuals, derivatives, and solver contracts.

#### 8.5.1 Boundary and region graph

Implement a neutral `boundary_region_graph_t` with explicitly oriented
interfaces and exterior components. It must record:

- bounded, exterior, and shell regions, including multiple components,
  handles, holes, and nested surfaces;
- a geometry provider (triangles, high-order panels, Fourier surfaces, NURBS,
  B-spline patches, or IGA patch graphs) with a common physical sampler;
- plus/minus ownership, outward normals, tangent frames, periodicity, field
  period, genus, and orientation transformations;
- trace spaces and work pairings on every interface, including scalar normal,
  vector tangential, RWG/BC, H(curl), H(div), and skeleton unknowns;
- fixed-topology identifiers and an explicit event when components split,
  merge, or change genus.

The graph is the only object passed to a free-boundary residual. It must not
assume cylindrical coordinates, nested flux surfaces, stellarator symmetry,
or a particular Fourier convention. Geometry construction and any adapter
from an application-owned representation remain outside this repository.

The public `boundary_region_graph_t` now supplies this neutral graph boundary:
it delegates oriented topology to the existing region graph and attaches
caller-owned interface genus/exterior flags plus contiguous physical
point/normal/weight samples with validated offsets and per-interface accessors.
Mesh generation, shape derivatives, coils, profiles, and application file
formats remain outside the type.

#### 8.5.2 Vacuum and exterior operator backends

Keep several interchangeable backends because no single method is best for
all geometries:

1. **High-order physical-space BEM.** Complete Laplace, Helmholtz, and
   magnetostatic/eddy-current layer potentials, principal-value traces,
   hypersingular regularization, near-neighbour quadrature, and far-field
   evaluation on curved panels. The existing Duffy, wedge, RWG, BC/RBC, and
   torus/sphere machinery is the baseline.
2. **NESTOR-like spectral Green operator.** A periodic surface sampler and
   Fourier convolution layer shall expose regularized toroidal Green kernels,
   modal Neumann solves, zero-mode constraints, and mode truncation/error
   metadata. This is useful for smooth toroidal surfaces and repeated
   free-boundary iterations, but it must share the physical-space trace API.
   The first neutral spectral trace contract is now public: a caller-owned
   P- or Q-branch toroidal modal expansion returns value and outward normal
   derivative on a constant-eta surface, with analytical modal coefficient,
   coordinate, and scale JVP/VJP actions. It is checked against the analytical
   toroidal harmonic/DtN mode and central/adjoint oracles; periodic surface
   sampling is now available through a batched grid wrapper; convolution,
   zero-mode policy, and truncation estimates are now explicit through a
   fixed-topology metadata/energy block: callers can reject an explicit mean,
   report retained/omitted rectangular-window modes, and differentiate supplied
   coefficient energies. Regularized Green convolution, Neumann compatibility,
   and unprovided-tail estimation remain caller-owned. A fixed-surface modal
   Neumann inversion with explicit zero-mode and resonance status plus
   analytical scale/η JVP/VJP actions is also public; global compatibility and
   Green-kernel assembly remain caller-owned. The public
   `apply_toroidal_modal_convolution` action now supplies the bounded-memory
   retained-mode convolution plus complex JVP/VJP products; a client still
   owns its regularized kernel coefficients, normalization, and compatibility
   policy.
3. **BIEST-like generalized-Debye-source operator.** Provide a second-kind
   surface formulation for Beltrami/Taylor/force-free fields, with explicit
   harmonic representatives for each handle and a resonance/status contract.
   This is a reusable vector boundary operator, not a Taylor-state or MRxMHD
   implementation. The first neutral generalized-Debye-source coordinate layer
   is now public: caller-owned scalar/cogradient/harmonic lifts are composed
   with source and period maps, with complex JVP/VJP actions and an independent
   dot-product oracle. The second-kind Green-kernel block, resonance policy,
   and topology-specific harmonic construction remain separate work. A
   matrix-free second-kind composition now applies a caller-owned `K` block to
   that lifted current and returns source/period residuals, with complete
   complex JVP/VJP actions and an independent reassembly/adjoint oracle; dense
   global assembly, Green-kernel construction, and Beltrami closure remain
   caller-owned.
4. **Virtual-casing/Biot--Savart maps.** A caller-owned volume or surface
   current is mapped to the field, normal field, tangential field, and their
   traces on arbitrary target surfaces. The map supports direct, reciprocal,
   and adjoint actions and is valid for non-nested or island-bearing surfaces.
   The affine-triangle RWG off-surface magnetic map now provides analytical
   source-geometry, current-data, target, and wave-number JVP/VJP actions with
   an independent central-reassembly and real-complex-adjoint oracle; the
   curved sphere and torus maps provide the corresponding geometry contracts.
   A batched target-point wrapper reuses that kernel for arbitrary sampled
   target surfaces and accumulates shared-source geometry/current cotangents,
   with independent single-target, central-reassembly, and adjoint tests.
   The target-trace wrapper now supplies normal and tangential magnetic fields
   plus analytical target-normal JVP/VJP actions, so a free-boundary residual
   can consume work-conjugate traces without reimplementing the projection;
   central-reassembly and real-complex-adjoint tests cover the full chain.
5. **Finite-domain/PML backend.** The same vacuum field can be represented in
   a bounded FEM/PML region. PML, DtN, and BEM results are compared on common
   interior targets, including a larger-domain control.

All backends return a typed operator with value, matrix-free action, assembled
matrix (when practical), JVP, VJP, residual contribution, and provenance. The
public `boundary_operator_contract_t` records these capabilities and the
fixed-topology/normalization metadata without owning callbacks or application
formats. A
backend may use dense, H-matrix, ACA, FMM, or low-rank storage internally;
that choice is not part of the physical contract.

#### 8.5.3 Coupling and free-boundary residuals

Add generic coupling forms, without selecting an equilibrium model:

- Johnson--Nédélec, Bielak--MacCamy, Costabel symmetric, and user-defined
  nonsymmetric FEM--BEM couplings for scalar and vector traces;
- FEM--DtN, FEM--NtD, BEM--DtN, mortar, Nitsche, and multiplier couplings for
  nonmatching interior, plasma, vacuum, wall, and coil surfaces;
- a declarative boundary residual for normal-flux continuity, tangential-field
  work pairing, supplied surface-current jumps, supplied total-pressure data,
  and supplied external-field traces;
- optional adjoint variables for nonsymmetric virtual-casing or response
  operators, so a non-energy-derived residual remains differentiable without
  pretending it is a symmetric energy minimization;
- integral constraints and nullspaces for fluxes, loop currents, gauges,
  harmonic fields, mean potentials, and supplied global invariants;
- block assembly and Schur complements that keep the plasma/interior,
  vacuum, wall, and interface blocks separately inspectable.

The public complex rectangular field-plus-constraint residual is the small
frequency-domain composition boundary for these blocks.  It carries no
equilibrium or wall convention: callers supply the rectangular field and
constraint actions and receive value, JVP, and real-part complex VJP products.
This is deliberately separate from the square linear-response interchange
operator and from application-specific FEM--BEM assembly.

The semantic `assemble_boundary_trace_residual` and
`assemble_complex_boundary_trace_residual` APIs now supply the real and
frequency-domain normal/tangential ports for these blocks. They compose
caller-owned trace maps, positive work weights, and supplied trace or jump
targets with complete allocation-free JVP/VJP products. They do not choose a
constitutive law or an equilibrium convention, so clients may use them for
virtual casing, FEM--BEM/DtN/PML coupling, surface currents, total-pressure
data, or wall ports.
Its JVP/VJP path is allocation-free and has a repeated-call regression, which
keeps batched free-boundary and FEM/BEM trace evaluations bounded in memory.

The library must expose the boundary residual in the form

\[
  R(q,\,g,\,\Gamma)=0,
\]

where (q) is a caller-owned field or trace, (g) is caller-owned source or
boundary data, and \(\Gamma\) is the oriented geometry graph. Generated
JVP/VJP products cover all three arguments at fixed topology. A client may
then add its own force-balance or pressure law without changing the operator
or derivative ABI.

#### 8.5.4 Conducting walls and reusable response matrices

Provide a neutral STARWALL-like response layer for thin or volumetric
conductors:

- surface current-potential spaces on triangulated, spline, and IGA wall
  patches, with holes and disconnected conductors;
- ideal-wall algebraic response and resistive-wall inductive/resistive block
  equations;
- response blocks (M_{ee}, M_{ey}, M_{ye}, M_{yy}), explicit units/signs,
  reciprocity or controlled nonsymmetry, and passivity/energy tests;
- reusable response-matrix serialization through a versioned neutral schema,
  never a JOREK/STARWALL file reader;
- dense-to-low-rank compression (SVD/ACA/H-matrix) with error bounds and a
  matrix-free action, while retaining a small dense oracle;
- static condensation and implicit JVP/VJP actions for wall currents;
- structure-preserving time coupling: the wall RL system is dissipative, the
  ideal-wall limit is constraint-like, and the coupled field/wall energy ledger
  is tested independently of a client's time integrator.

The neutral `advance_resistive_wall_midpoint` primitive now supplies the
implicit-midpoint RL step, analytical matrix/state/voltage/step JVP/VJP
actions, and an independent discrete energy/input/dissipation ledger. This is
the reusable structure-preserving wall block; geometry, surface-current basis,
and STARWALL normalization remain caller-owned.

The neutral `condense_wall_response_blocks` primitive now performs the complex
Schur reduction (M_{ee}-M_{ey}M_{yy}^{-1}M_{ye}), with matrix-free analytical
JVP/VJP actions and an independent reassembly/real-complex-adjoint oracle.
Retained factors, low-rank compression, and field-port assembly remain
separate layers.

The versioned `fortfem_linear_response_schema` text interchange now provides
a bounded metadata/block round-trip and small dense oracle for external wall
and field-response adapters. It is intentionally not a reader for any
application-specific format; compressed and production-scale storage remain
caller-owned.

The neutral `fortfem_complex_low_rank_response` contract now adds deterministic
cross approximation, a relative residual certificate, matrix-free `U V^T`
actions, and factor/input JVP/VJP products. It is a bounded small-oracle layer;
production ACA/SVD/H-matrix/FMM implementations remain interchangeable
caller-owned backends.

This response layer is also the reusable foundation for linear ideal/resistive
perturbations, external-kink and resistive-wall problems, and nonlinear
time-dependent clients. It must not contain a plasma closure, an MARS/GPEC
normalization, or a JOREK state vector.

#### 8.5.5 Shape calculus and continuation contracts

At fixed topology, generate analytical derivatives of panel and spline
geometry, normals, surface Jacobians, Green kernels, principal-value
subtractions, Fourier maps, trace transformations, coupling blocks, response
matrices, and solved states. The common operations are geometry JVP/VJP,
data JVP/VJP, implicit solve JVP/VJP, and adjoint objective evaluation.

Continuation is a client policy, but the numerical layer must provide:

- a public pseudo-arclength residual with predictor/tangent value, JVP, and
  VJP actions; predictor construction, corrector solves, and tangent
  construction remain client-owned, while the public Euclidean tangent
  normalization value/JVP/VJP fixes the scale contract;
- line-search/trust-region and pseudo-transient hooks;
- event diagnostics for a cut crossing, separatrix, topology change, or
  resonance;
- frozen-topology derivatives and an explicit invalid-derivative status across
  a discrete topology event;
- independent complex-step/central and real-part adjoint checks for every
  smooth block.

The neutral `assemble_pseudo_arclength_residual` contract now supplies the
fixed-topology augmented residual and complete analytical JVP/VJP products.
It is deliberately free of free-boundary geometry, force balance, or a
nonlinear solver, so the same primitive can be composed by VMEC/GVEC-like,
SPEC-like, MARS/GPEC-like, JOREK-like, elasticity, and wave clients.
The matching `normalize_pseudo_arclength_tangent` primitive supplies the
allocation-free unit tangent blocks and norm with zero-tangent rejection and
complete analytical derivatives.
The `evaluate_residual_merit` primitive supplies the positive-weighted
least-squares merit and analytical JVP/VJP; acceptance rules and nonlinear
solver state remain client-owned.
The `assemble_pseudo_transient_residual` primitive now supplies the separate
positive-step (R+M(u-u_{old})/\Delta t) continuation stabilizer with complete
matrix/state/time JVP/VJP products; it is explicitly not a symplectic time
integrator.
The `evaluate_step_reduction` primitive now supplies smooth actual and
actual/predicted reduction diagnostics with complete JVP/VJP products;
backtracking, trust-region radii, and acceptance thresholds remain explicit
client policy.
The `classify_continuation_event` primitive now reports deterministic
sign-crossing and near-zero margins with the global minimum margin. It is the
explicit topology-event boundary for cuts, separatrices, resonances, and
interface graphs; derivatives across a reported event remain invalid by
contract.

#### 8.5.6 Geometry-independent 2D and 3D benchmark ladder

The following examples are required before an application adapter is trusted.
They use manufactured fields or public analytical solutions, not equilibrium
profiles or production-code inputs:

1. exterior circle and cylinder Laplace/Helmholtz problems, comparing analytic
   Fourier modes, BEM, DtN, FEM--BEM, and PML;
2. an axisymmetric Grad--Shafranov-shaped *scalar* free-boundary analogue with
   supplied source and boundary traces, comparing Johnson--Nédélec and
   symmetric FEM--BEM coupling;
3. a smooth toroidal vacuum Neumann problem comparing NESTOR-like Fourier,
   physical-space BEM, DtN, and larger-domain PML fields;
4. a toroidal-shell Beltrami manufactured field comparing a BIEST-like
   generalized-Debye-source map with compatible H(curl) volume assembly,
   including handle/harmonic modes and resonance rejection;
5. a virtual-casing surface-current test on sphere, cylinder, and torus with
   direct field, normal/tangential trace, reciprocity, and shape-derivative
   oracles;
6. a generic two-region interface with an imposed sheet current and supplied
   total-pressure jump, checking normal-flux continuity, tangential jump, and
   work balance without solving an MHD equilibrium;
7. a thin conducting shell with ideal and resistive limits, response-matrix
   compression, passivity, and a structure-preserving field/wall time step;
8. a non-axisymmetric multi-surface graph with holes and nonmatching IGA,
   Fourier, triangular, and FEM/PML discretizations, including a fixed-
   topology shape derivative;
9. a linear forced-response block that compares vacuum, ideal-wall, and
   resistive-wall limits and exports neutral samples for external GPEC/MARS/
   GLISS/JOREK oracle runs;
10. a free-boundary continuation fixture in which only a manufactured boundary
    trace is varied, reporting residual, force-like supplied functional,
    conditioning, solve time, and derivative error.
11. a nested-surface variational fixture with a toroidal embedding, prescribed
    fluxes, a NESTOR-like vacuum map, and a physical surface-shape plot;
12. an ideal perturbed-response fixture with a rational surface, explicit
    shielding-current trace, penetrated resonant field, and an optional
    supplied resistive inner-layer oracle;
13. an Eulerian relaxation fixture initialized from the nested solution, with
    an island-forming perturbation, force residual, divergence defect, and
    topology-event report;
14. a multi-region relaxed fixture with independent Beltrami blocks, ideal
    interfaces, pressure jumps, flux/helicity constraints, and comparison to
    a compatible H(curl) volume solution;
15. a reduced resistive island/wall continuation fixture that reports
    branch multiplicity, energy input, dissipation, island width, and
    hysteresis without claiming an arbitrary static resistive equilibrium.

Every example's first plot is the physical solution (surface, slice, vector
field, or field line), followed by mesh/geometry, residual/conservation,
convergence, performance, and derivative plots. The gallery stores numerical
CSV/JSON only; generated images stay out of git and are rebuilt by the docs
workflow.

#### 8.5.7 Equilibrium-model selection

FortFEM will expose several equilibrium-model foundations. They share
geometry, traces, compatible spaces, block residuals, continuation, and
derivative contracts, but they do not impose the same magnetic topology.
Nested ideal equilibria, ideal perturbed responses, relaxed non-nested
equilibria, and resistive dynamics are separate model branches.

**Nested ideal branch.** A surface embedding

\[
  \mathbf{x}=\mathbf{x}(\rho,\theta,\zeta)
\]

represents every \(\rho=\mathrm{const.}\) surface as a torus. The magnetic
field representation enforces \(\mathbf B\cdot\nabla\rho=0\) kinematically.
The neutral `evaluate_nested_surface_geometry` contract provides the first
caller-owned Fourier/radial embedding slice, including `(R,Z,lambda)`,
Cartesian samples, metric/volume diagnostics, and fixed-topology coefficient
derivatives. It deliberately does not impose an equilibrium equation or
plasma-specific normalization.
The residual or constrained energy uses supplied fluxes, currents, pressure
profiles, and other physical invariants. Radial IGA or high-order finite
elements and Fourier or spline surface bases provide a coordinate-clean
VMEC/GVEC/DESC-like foundation. A NESTOR-like toroidal vacuum operator is the
preferred exterior backend for smooth periodic surfaces. This branch is
accurate and efficient for nested equilibria. Nestedness is a model
assumption, not a diagnostic that the solver may silently enforce after a
topology change.

**Ideal perturbed-response branch.** A nested equilibrium is followed by a
linear displacement and field-response solve. GPEC/IPEC-like outer response
blocks, vacuum and wall operators, and supplied resonant-layer data share the
same trace interface. At a rational surface, an ideal shielding current is a
surface distribution represented by an explicit \(\mathbf K\) trace or a
field jump. The response layer may match an external resistive inner-layer
model. It must report the sheet current, penetrated resonant field, and the
validity of the fixed-topology derivative. Ideal shielding response is not a
replacement for the equilibrium branch.

**Eulerian relaxed branch.** A SIESTA-like energy relaxation starts from a
nested state and updates volume fields or displacements on a fixed physical
mesh. The representation can form islands and stochastic regions without
requiring flux coordinates everywhere. FortFEM supplies compatible
\(H(\mathrm{curl})\)/\(H(\mathrm{div})\) spaces, force and divergence
residuals, topology-event reporting, pseudo-transient hooks, and physical
preconditioner blocks. A client supplies the pressure, transport, and
relaxation model.

**Multi-region relaxed branch.** A region graph contains independently
represented plasma volumes and ideal interfaces. Beltrami or generalized
relaxed-field blocks, flux and helicity constraints, total-pressure jumps,
and interface shape derivatives provide a SPEC/MRxMHD-like foundation. A
region may contain islands or chaotic field lines while the interfaces remain
explicit. BIEST-like generalized-Debye-source operators and compatible
volume FEEC are interchangeable interior backends.

**Resistive dynamic branch.** Finite resistivity, current drive, flow,
pressure or heat transport, wall response, and boundary forcing define a
time-dependent problem. A steady state is then a result of the specified
dynamics and continuation path, not an arbitrary constrained minimum. This
branch is required for reconnection, finite island saturation, island
overlap, stochastic transport, error-field penetration, wall locking, and
bifurcation studies. Structure-preserving time integration and energy or
helicity ledgers remain mandatory.

Model selection follows the physical question:

| Question | Foundation branch |
| --- | --- |
| Smooth equilibrium, coil design, or reconstruction with nested surfaces | Nested ideal branch with NESTOR-like or physical vacuum BEM |
| Small non-axisymmetric perturbation and ideal shielding | Nested equilibrium plus ideal perturbed-response branch |
| Finite islands or stochastic regions in an equilibrium calculation | Eulerian relaxed or multi-region relaxed branch |
| Reconnection, island saturation, wall locking, or hysteresis | Resistive dynamic branch |

The solver must preserve non-uniqueness as observable structure. Constraints
are accepted only when they represent physical invariants, boundary data, or
gauge conditions. Pseudo-arclength continuation, deflation, multi-start
initialization, topology-event diagnostics, and Hessian or stability spectra
track distinct solution branches instead of selecting one silently. A
current-sheet limit, a relaxed multi-region state, and a finite-resistivity
state are compared as separate solutions with declared assumptions.

The corresponding FortFEM interfaces are:

- `nested_surface_variational` for the fixed nested topology;
- `ideal_perturbed_response` for displacement, shielding-current, vacuum, and
  wall blocks;
- `eulerian_relaxation` for fixed-domain non-nested ideal relaxation;
- `multi_region_relaxed` for explicit region graphs and ideal interfaces;
- `resistive_dynamics` for client-owned time-dependent closures.

These names describe external composition contracts. FortFEM does not contain
VMEC, GVEC, DESC, SIESTA, SPEC, GPEC, MARS, or resistive-MHD application
physics.

#### 8.5.8 Implementation order and method-selection rule

Implement the foundation in this order:

1. oriented region graph, trace work pairings, harmonic/nullspace metadata,
   and neutral boundary sampling;
2. high-order Laplace/Helmholtz and magnetostatic layer potentials with
   singular/near quadrature and analytical data/geometry derivatives;
3. FEM--BEM, FEM--DtN, and PML coupling plus scalar 2D/3D manufactured tests;
4. virtual-casing and target-surface reconstruction, then toroidal periodic
   Fourier/NESTOR-like operators;
5. vector H(curl) BEM with generalized-Debye-source/BIEST-like harmonic
   blocks, tree--cotree gauge reduction, and Beltrami manufactured tests;
6. wall surface-current spaces, response blocks, compression, passivity, and
   structure-preserving wall coupling;
7. shape derivatives, implicit derivatives, continuation hooks, and the full
   sphere/cylinder/torus benchmark ladder;
8. external oracle adapters and neutral interchange files, never application
   readers or production equilibrium/state algorithms.

NESTOR-like Fourier BIE is preferred for smooth periodic toroidal surfaces and
   repeated modal solves; BIEST-like second-kind generalized Debye sources are
   preferred for vector Beltrami fields and nontrivial topology; physical BEM
   and FEM--BEM are preferred for arbitrary geometry or material coefficients;
   DtN/NtD is preferred when the same exterior geometry is reused many times;
   PML is the independent finite-domain control; virtual casing is preferred
   for supplied interior-current equivalence. These are complementary
   backends, not competing replacements. Each choice must be benchmarked
   against at least one independent method and a manufactured solution.

#### 8.5.9 MHD foundation task ledger

This ledger is the authoritative MHD-specific work breakdown. It translates
the generic contracts above into dependencies that an external equilibrium,
linear-response, edge, or time-dependent client can consume. A status applies
to the reusable ingredient only; it never means that FortFEM implements the
named application or its plasma closure.

| ID | Status | Foundation task | Acceptance gate and oracle |
| --- | --- | --- | --- |
| MHD-00 | **done** | Keep the application boundary explicit: no GEQDSK/COCOS readers, equilibrium profiles, coil or vessel models, Braginskii/sheath closures, species, or production VMEC/GVEC/SPEC/GPEC/MARS/GLISS/JOREK algorithms | A neutral-array example and provenance check show that an external client owns the physics and file conversion |
| MHD-01 | **active** | Stabilize one equation-as-data ABI for scalar, vector, tensor, complex-frequency, and multiplier fields, including units, normalization, constraints, traces, residual/JVP/VJP, and converged-state derivatives. The neutral `equation_objective_registry_t` now provides deterministic named equation/objective/constraint blocks, offsets, tags, active/fixed metadata, validation, and pack/JVP/VJP actions; `evaluate_equation_objective_callbacks` dispatches caller-owned value/JVP/VJP callbacks by typed block metadata with deterministic packing and cotangent accumulation; constitutive formulas and callbacks remain external | Manufactured mixed-field residual with independent matrix, finite-difference, real-adjoint, and invalid-shape tests; registry ordering, duplicate rejection, tags, copy semantics, callback failure handling, and pack/callback transpose oracles pass; no plasma state vector is hard-coded |
| MHD-02 | **active** | Finish the compatible spatial core: H1--H(curl)--H(div)--L2, Fourier-FEM, IGA patch graphs, orientation/periodic maps, harmonic forms, and tree--cotree direct gauge reduction for curl--curl systems. The topology-only `build_multipatch_signed_dof_map` and its 2D/3D B-spline compositions now construct arbitrary-patch signed maps, face permutations, and orientation-cycle checks. The geometry-aware `assemble_geometry_mortar_trace_coupling` layer now chains reference quadrature weights with caller-owned physical surface metrics and exposes trace/metric JVP/VJP actions. The neutral `sample_physical_surface_geometry` contract now turns provider-evaluated 3D surface points and tangent columns at fixed quadrature coordinates into oriented normals and positive measures, with strict validation and JVP/VJP actions. The communicator-free `distributed_trace_layout_t` now composes complete `partition_layout_t` rank ledgers, rejects missing global IDs, validates one owner and consistent IDs for each physical trace row, and exposes deterministic duplicate-ghost value/JVP/VJP reductions for FEM/BEM/DtN and mortar rows. The new `physical_trace_ownership_t` carries finite physical coordinates alongside IDs and owner/ghost masks and compares independently ordered partitions by global ID before exchange. The `diagnose_tree_cotree_iga_invariance` contract now checks the high-order/IGA tree selector, direct reduced solve, and oriented loop-period invariance under a global signed-map change. New neutral `collocation_grid_t`/`direct_fourier_plan_t` contracts provide bounded tensor-grid chunking, axis/periodic maps, flattening, and direct transform/adjoint actions; `fourier_zernike_basis_t` provides deterministic scalar \(l,m,n\) enumeration and axis parity metadata. MPI transport, geometry-aware owner selection, vector/tensor parity, and general patch-graph assembly remain open | Slab, cylinder, sphere, and torus complexes satisfy incidence identities, divergence/curl defects, loop periods, gauge reduction, and direct-solve parity; independent cross-product, central-difference, real-adjoint geometry oracles pass; fixed selectors are not differentiated |
| MHD-03 | **active** | Compose tensor-valued pressure/stress and Maxwell traction blocks, including CGL-shaped parallel/perpendicular terms, gyrotropic/nonsymmetric terms, stress work, force divergence, and interface jumps. The neutral `assemble_total_pressure_jump` block now composes caller-owned scalar pressure and vector-field traces into the magnetic total-pressure residual with exact real JVP/VJP actions; the `cgl_pressure_tensor` gallery samples a physical 2D anisotropy map and the `cgl_pressure_3d_gallery` samples a physical 3-D pressure trace, oblique field direction, and traction arrows before their diagnostics; constitutive closures remain external | Independent tensor contraction, traction, force-volume, normal-flux, tangential-jump, power, and real/complex JVP/VJP oracles; angular-momentum and energy defects are reported |
| MHD-04 | **active** | Add a generic force-balance residual composition for supplied magnetic, pressure/stress, inertial, and body-force blocks, with explicit volume, boundary, and sheet-current terms. `assemble_force_balance_residual` now provides the closure-neutral three-term weak composition with complete real JVP/VJP actions, and its scalar volume/boundary/sheet contraction is FortSym-generated with byte-checked value/JVP/VJP kernels; tensor and constitutive blocks remain caller-owned | A manufactured equilibrium-like field on slab/cylinder/torus has zero residual under the supplied data; fitted, cut, DG, and IGA forms agree on common samples |
| MHD-05 | **active** | Add the nested-surface coordinate and variational composition used by VMEC/GVEC clients and consumed by DESC-like direct-force clients: toroidal embedding, radial FE/IGA, Fourier modes, flux/constraint blocks, shape derivatives, and NESTOR-like vacuum trace. The neutral `evaluate_flux_surface_average` reduction now provides fixed-topology scalar/vector surface averages, an explicit positive denominator diagnostic, and value/JVP/VJP actions for samples, quadrature weights, and optional geometry factors; the nested geometry contract now also exposes fixed-topology JVP/VJP actions for sampled ((\rho,\theta,\zeta)) coordinates; the axis-regular radial evaluator now supplies finite complex polynomial value/JVP/VJP actions for caller-selected powers; the `nested_surface_solution` gallery now shows a physical 3D torus before parameter/radial diagnostics; surface labels and physics remain caller-owned | A supplied nested toroidal embedding reproduces manufactured fluxes and force residuals; physical surface plots, mode energies, conditioning, weighted-average, coordinate-adjoint, and shape-derivative checks are published |
| MHD-06 | **active** | Add the multi-region relaxed composition used by SPEC/SPECTRE-like clients: independent region fields, Beltrami/curl-eigen blocks, ideal interfaces, pressure balance, flux/helicity constraints, and harmonic loop data. The neutral `assemble_beltrami_residual` contract now composes caller-supplied `curl(B)-lambda B` blocks with optional divergence, flux, and helicity value/target rows and exact real JVP/VJP actions; `compare_beltrami_two_region_residual` now compares that compatible H(curl) path with an independent curl-eigen oracle and reports constraint norms, while `validate_beltrami_resonance` rejects caller-supplied resonant gauge/eigenvalue margins; the geometry-labelled `compare_beltrami_shell_residual` fixture now independently compares weighted slab/toroidal-shell paths and closes flux, helicity, and quadratic field-energy ledgers; `assemble_total_pressure_jump` supplies the neutral magnetic total-pressure interface row with exact JVP/VJP; field construction, gauge, topology, and equilibrium physics remain external | The current two-region manufactured residual/JVP/VJP/finite-value and zero-sized-constraint tests pass, as does the independent 10-check parity/resonance fixture and the 9-check slab/toroidal-shell ledger fixture. Physical interfaces, pressure balance, and harmonic loop data remain client acceptance gates |
| MHD-07 | **active** | Add the Eulerian non-nested composition used by SIESTA-like clients: fixed-domain compatible fields, supplied force/divergence residuals, pseudo-transient hooks, topology events, and free-boundary traces. The neutral `assemble_eulerian_nonnested_residual` contract now concatenates caller-owned force/divergence residual vectors, optionally adds a precomputed stabilization from `assemble_pseudo_transient_residual`, and reports continuation-event metadata with exact value/JVP/VJP actions; topology-event metadata is held fixed by derivatives and no relaxation closure is selected. The physical-first `eulerian_island_gallery` now visualizes a manufactured figure-eight separatrix with the supplied tangent field before the residual diagnostic, and emits CSV/JSON provenance without a plasma closure | An island-forming manufactured perturbation reports field topology, divergence, force residual, event status, and fixed-topology derivative validity without selecting a relaxation closure; cross-composition with compatible H(curl)/H(div) state samples and free-boundary traces remains |
| MHD-08 | **active** | Add linear ideal/resistive perturbation composition used by GPEC/MARS/GLISS-like clients: inertia, Lorentz, pressure/stress, vacuum, wall, resistive, singular-layer, and complex-frequency blocks. The neutral `assemble_linear_perturbation_operator` now keeps all seven caller-supplied blocks separately inspectable, composes their canonical complex-frequency factors, and exposes exact JVP/VJP actions plus declared reciprocity/passivity/normalization metadata. The manufactured `linear_perturbation_response` gallery now crosses that composition with the forced-response residual, reconstructs a physical complex Fourier state first, and publishes frequency-sweep reciprocity, passivity, JVP, residual, and timing data. The neutral `assemble_linear_response_eigen_cross_residual` now composes a generalized eigen residual, retained complex response/coupling block, and caller-owned ideal-shielding trace under fixed metadata, with independent matrix and exact complex JVP/VJP actions. The neutral `singular_layer_matching_gallery` now exercises independently sized analytical inner/outer trace rows and their matched jump before a deliberate mismatch/JVP diagnostic; plasma closures and readers remain external | Seven-block value, central-difference, real-complex adjoint, metadata, rejection, forced-response, generalized-eigen/retained-response/ideal-shielding cross-composition, and neutral trace-matching oracles pass. Penetrated resonant field, asymptotic layer models, and external response-matrix fixtures remain |
| MHD-09 | **active** | Finish the sheet-current and resonant-layer bridge: tangential \(\mathbf K\) trace, Ampere jump, normal-flux continuity, total-pressure jump, explicit \(\mathbf K\delta_\Gamma\), fitted/cut/DG/IGA representations, and resolved resistive-layer limit. The neutral `compare_sheet_current_representations` slab contract compares the explicit surface ledger \(A\mathbf K\) with the existing normalized Gaussian layer using an independent quadrature oracle; `compare_sheet_current_surface_representations` extends that ledger to caller-owned fixed-topology positive surface quadratures and unit normals labelled slab, cylinder, sphere, and torus, with orientation, measure, and torus-JVP checks; `singular_layer_matching_gallery` supplies a physical-first 1D analytical trace/jump bridge and centered-difference JVP oracle without selecting a singular asymptotic or plasma closure; the new `surface_current_potential_metadata_t` and `assemble_surface_current_potential` map caller-owned tangent gradients plus harmonic loop representatives to oriented \(\mathbf K\), integrated current, period ledger, and exact JVP/VJP actions with an explicit fixed gauge/topology contract; surface geometry and PDE closures remain caller-owned | Slab HKT-style, cylinder, sphere, and torus fixtures compare all four representations, conserve integrated current, and show the declared non-differentiability at topology events |
| MHD-10 | **active** | Compose vacuum/free-boundary and conductor blocks: physical BEM, DtN, NESTOR-like Fourier, virtual casing, BIEST-like vector maps, PML, wall response, and retained Schur/field-split factors. The neutral `assemble_free_boundary_port_residual` now combines fixed-topology physical traces with supplied exterior/vacuum targets, positive work weights, and optional sheet/jump targets with exact JVP/VJP actions; `compare_boundary_operator_parity` compares FEM/BEM/DtN/PML candidate fields on common weighted physical samples while enforcing shared topology/equation-space/units metadata, and `compare_larger_domain_solution` compares two complex fields from the same or different backends on common interior samples with distance-increase diagnostics; `evaluate_weighted_boundary_response_diagnostics` now supplies a work-weighted reciprocity defect and conservative Hermitian/Gershgorin passivity bound for any square complex boundary response; `evaluate_source_trace_map` supplies caller-owned dense source-to-trace value/JVP/VJP products and a weighted reciprocity defect; `boundary_operator_contract_t` now carries explicit scalar/normal/tangential/mixed trace-channel and work-pairing metadata with backward-compatible scalar defaults; the physical Maxwell torus/cylinder trace fixture independently checks weighted FEM/BEM/DtN/PML errors and mixed-space rejection; `test_boundary_backend_registry` now independently exercises all eight backend tags, mixed trace/work metadata, copy semantics, and invalid-tag rejection; boundary kernels and free-boundary physics remain external | Common interior target fields agree across methods within declared discretization error; larger-domain, reciprocity, far-field, passivity, memory, and performance controls are recorded |
| MHD-11 | **active** | Complete strongly anisotropic and field-aligned operators for magnetic diffusion, resistivity, viscosity, heat/conduction, pressure, and wave impedance, including FCI, Fourier, curvilinear, and IGA pullbacks. The physical-first `field_aligned_anisotropic_3d_gallery` now exercises an oblique \(k_\parallel/k_\perp=5000\) tensor on a tetrahedral unit-cube solve with solution/flux arrows, mid-plane slice, CSV diagnostics, and an independent projector/Hessian oracle; constitutive closures and plasma physics remain external | Limits \(k_\parallel/k_\perp\to1,0,\infty\) and Hall/gyrotropic skew parts pass energy/dissipation, convergence, conditioning, and tensor/geometry derivative tests |
| MHD-12 | **active** | Keep ideal and dissipative time structure separate: mixed pressure/velocity, displacement/momentum/stress, Maxwell, wall RL, and nonlinear client-owned state blocks with midpoint, Cayley, symplectic Euler, Strang, and dissipative splits. The mixed-elasticity gallery now adds a physical 2D stress contour and displacement quiver, and `mixed_wave_pressure_displacement_gallery` adds physical 1-D/2-D/3-D pressure/displacement states with an independent modal midpoint oracle | Energy, helicity, cross-helicity, mass, divergence, charge, reversibility, symplectic-form, and wall-dissipation ledgers are checked on 1D/2D/3D manufactured waves; no dissipative step is called symplectic |
| MHD-13 | **active** | Provide generic edge/SOL foundations for GRILLIX/BOUT++/PARALLAX-like clients: FCI tracing, parallel/perpendicular split, conservative fluxes, material traces, terminal events, and vector ledgers. The neutral `evaluate_fci_power_flux_ledger` contract now reports the staggered parallel flux and signed parallel/perpendicular power pairing with exact fixed-topology JVP/VJP actions; FortSym owns the pointwise products while FCI geometry, material traces, and closures remain caller-owned. The physical-first `fci_sol_gallery` now traces a helical torus field-line family and applies the same conservative FCI map to a manufactured SOL potential, with CSV/JSON provenance and a neutral closure label | MMS source/terminal balances, anisotropic negative-energy tests, power/flux ledger and derivative oracles, boundary-event ownership, reproducibility seeds, and 1D/2D/3D solution plots pass without species or sheath physics |
| MHD-14 | **active** | Compose a bounded nonlinear resistive-MHD operator family: Faraday/Ampere-compatible blocks, supplied momentum/pressure/tensor terms, anisotropic transport, wall/free-boundary ports, continuation, and branch/event diagnostics. The neutral `assemble_nonlinear_resistive_mhd_residual` contract now visits caller-owned Faraday, Ampere, momentum, pressure, tensor, anisotropic-transport, wall, and free-boundary callbacks, sums their nonlinear residuals, and reports stored-energy, signed-input, nonnegative-dissipation, and instantaneous power-balance diagnostics with exact JVP/VJP actions; the neutral branch-history contract now validates strictly monotone caller-owned parameter/state/residual/energy samples, reports deterministic branch multiplicity and fixed-grid hysteresis, and supplies a path metric with parameter/state/residual/energy JVP/VJP actions; constitutive laws, state selection, and closure remain external | A reduced manufactured island/wall continuation case reports residual, input, dissipation, branch multiplicity, hysteresis, and derivative validity; the current eight-block nonlinear manufactured oracle and branch-history fixture cover value/JVP/VJP and ledger/diagnostic checks; the seeded 20-case nonlinear composition campaign independently checks caller-owned block sums, central-difference JVPs, stored/input/dissipation balance, and the real residual-plus-ledger adjoint identity; the seeded 20-case branch-history campaign independently checks monotone-path metric and JVP/VJP oracles, deterministic multiplicity/hysteresis, and non-monotone rejection; the physical-first `nonlinear_resistive_mhd_gallery` fixture now renders a neutral island-flux/wall-response state before branch and ledger diagnostics and emits CSV/JSON provenance |
| MHD-15 | **active** | Build the external oracle ladder and neutral data path for CHEASE, FreeGS, VMEC/PARVMEC, GVEC/DESC, SPEC/SIESTA, GPEC/MARS, GLISS, STARWALL, JOREK, FreeFEM, MFEM, and FEniCSx. The public `oracle_manifest_t` and versioned `fortfem-oracle-manifest-2` text round-trip now record code/revision/license, case and coordinate checksums, declared units/normalization, tolerances, phase timings, peak memory, runner identity and required hardware context, runner/repetition metadata, and an optional sister-repository URI; schema-2 round-trip, missing-hardware rejection, manufactured-Poisson, and metadata-only target-ladder fixtures independently validate solution, provenance, and performance fields without launching external codes; code names and data remain caller-owned | Common physical sampling, checksums, licenses, executable revisions, tolerances, performance metadata, runner hardware, and optional sister-repository data; no external source or reader enters FortFEM |
| MHD-16 | **active** | Finish the MHD gallery and benchmark contract in increasing complexity, with physical solution first: 1D slab, 2D circle/cylinder, 3D sphere/torus, vector arrows/surfaces/field lines, mesh completeness, convergence, invariants, and timings. Simple Poisson, tetrahedral H1, tetrahedral Nédélec, regularized sheet-current, CGL tensor, 3-D CGL pressure/traction, linear response, neutral singular-layer matching, nested-torus, FCI/SOL, Eulerian-island, mixed-elasticity, mixed pressure--displacement waves, `biro_tree_cotree_benchmark`, `biro_tree_cotree_3d_gallery`, `team7_neutral_benchmark`, `team13_neutral_benchmark`, and `free_boundary_port_gallery` now lead with physical plots; the tetrahedral previews record physical-before-diagnostics execution order, with readable sampled scalar/vector fields, mesh edges, and quiver arrows, while the toroidal Maxwell radiation preview uses a connected periodic surface rather than an unordered point cloud. The remaining TEAM ladder, exact paper geometries, and external-code parity remain sister-data work. All pages retain ignored generated media plus CSV/JSON provenance | Every page has a solution preview before diagnostics, generated FortPlot media, CSV/JSON provenance, no checked-in images, and a bounded memory/timeout record |
| MHD-17 | **active** | Add the DESC-like direct nested-surface force-balance/optimization foundation without implementing DESC physics: inverse toroidal-flux coordinates \((R,Z,\lambda)\), Fourier--Zernike or equivalent axis-regular radial bases, parity and axis regularity, collocation force components, profile/objective/constraint callback registries, exact JVP/VJP and Hessian-vector hooks, perturbation/continuation/deflation, flux-surface averages/Boozer-transform hooks, near-axis data hooks, and fixed/free-boundary external-field and sheet-current residual ports. The neutral `build_axis_regular_mode_table` validator reports scalar \(\rho\)-power/parity requirements and deterministic conjugate-safe mode ordering; `evaluate_axis_regular_radial_basis` now evaluates caller-selected finite complex radial polynomials with exact coefficient/rho JVP/VJP products and the same minimum-power/parity contract; `evaluate_flux_surface_average` supplies weighted diagnostic reduction; nested geometry supplies coordinate-sample JVP/VJP products; `evaluate_force_balance_objective` now provides a positive-weighted direct-force least-squares value/JVP/VJP contract plus an exact residual/weight Hessian-vector action through `evaluate_force_balance_objective_hvp`; `assemble_free_boundary_port_residual` supplies fixed-topology trace/exterior-target/sheet-jump composition; vector/tensor shifts remain caller-owned | Manufactured axisymmetric and 3D toroidal states show axis regularity, exponential radial/mode convergence, direct force residual closure, objective/constraint derivatives (including weighted averages and coordinate adjoints), and direct-force Hessian-vector finite-difference parity; perturbation and continuation parity, and free-boundary/sheet-current boundary residuals; all profile laws, readers, and production optimization remain external |

The direct-force campaign is now also a physical-first gallery fixture: it renders a
manufactured toroidal state with a visible force-vector field before its parameter
and derivative diagnostics, and emits CSV/JSON timing and finite-difference
provenance. It remains a neutral contract; profiles, readers, and production
optimization stay caller-owned.

The neutral `critical_point_metadata_t` contract now records caller-supplied
two-dimensional critical-point candidates, gradient/Hessian event margins,
regular/O/X/degenerate classification, and limiter/separatrix labels with fixed-
topology JVP/VJP actions. It does not locate nulls or assign plasma physics.

#### Free-boundary completion checklist

The rows above cover the field and residual pieces, but a reusable free-boundary
client also needs the following small, application-neutral contracts.  They are
deliberately phrased in terms of sampled geometry and linear maps: coil models,
profiles, equilibrium selection, and input readers stay outside FortFEM.

1. **Surface-current potential and loop basis.**  Decompose a tangential surface
   current into a single-valued scalar potential plus one harmonic representative
   per handle, with explicit cuts, periods, gauge, orientation, and units.  The
   same map must feed virtual casing, NESTOR-like vacuum maps, BIEST-like
   generalized-Debye sources, and conducting-wall response blocks.
2. **Source-to-trace/inductance maps.**  Provide a neutral map from caller-owned
   source samples (filament, panel, volume, or spline coefficients) to target
   normal/tangential fields and work-conjugate fluxes.  The
   `evaluate_source_trace_map` contract now supplies dense value/JVP/VJP
   products and a weighted reciprocity defect for a fixed caller-owned map;
   self/near-panel regularization status, retained response blocks, and
   shape/parameter actions remain explicit follow-up metadata.  This is a
   numerical map, not a coil or vessel model.
3. **Moving-surface shape calculus.**  Differentiate point positions, tangent
   frames, normals, measures, periodic seam identifications, and singular or
   coincident-panel quadrature consistently for a fixed region topology.  Report
   a topology-event status when a component, handle, cut, or critical point
   changes; never differentiate silently through that event.
4. **Axisymmetric critical-point metadata.**  Add a neutral descriptor for
   magnetic/null points, X/O points, limiter or separatrix contours, and their
   event margins.  It stores sampled geometry and classification only; locating
   points from a supplied scalar field and all equilibrium interpretation remain
   external.  This is the common 2-D hook needed by CHEASE/FreeGS-like clients.
5. **Coupled free-boundary solve bookkeeping.**  Keep interior, vacuum, wall,
   surface-current, and geometric blocks separately inspectable in a residual
   and Schur/Jacobian graph.  Record block scaling, units, nullspaces, retained
   factors, continuation parameter, conditioning, and work/energy ledgers so
   NESTOR/BIEST, FEM--BEM/DtN, PML, and virtual-casing paths can be swapped in
   the same manufactured solve.

The acceptance fixture for this checklist is a supplied circle and toroidal
surface with a manufactured current potential and one handle mode: direct
source-to-trace, spectral vacuum, and compatible volume paths must agree in
weighted flux/work, pass the reciprocal/adjoint and shape-JVP checks, and reject
an intentionally changed topology.  The fixture may be generated in the
gallery, but its numerical data—not any external code or reader—is the oracle.

The equation/objective registry is the neutral packing layer for MHD-01 and
MHD-17 callback composition; all constitutive and profile callbacks remain
external.

The companion `evaluate_equation_objective_callbacks` adapter now dispatches
typed caller-owned value/JVP/VJP callbacks by registry block, packs their
rows deterministically, and accumulates state cotangents.  This is callback
composition only; profile laws and optimization objectives remain external.

The neutral `assemble_deflated_residual` contract now supplies a smooth,
caller-owned root-deflation wrapper with exact state/residual JVP/VJP actions.
It is intended to compose with the existing pseudo-arclength residual for
DESC-like continuation and branch searches; reference-root selection and
nonlinear acceptance remain external.

Dependency order is MHD-00/01 → MHD-02/03 → MHD-04/09 → MHD-05/06/07/17 →
MHD-08/10 → MHD-11/12/13 → MHD-14/15/16. In particular, a generic
coupled Schur layer is an MHD-01/02 solver dependency, not an equilibrium
implementation; its off-diagonal blocks remain caller-owned. A gallery or
oracle may use a name such as “HKT”, “VMEC-like”, or “JOREK-like” only to
identify the mathematical benchmark and provenance, never to imply that the
corresponding production code is being reimplemented here.

#### DESC-specific foundation checklist

DESC is a direct force-balance and optimization target, not merely another
variational VMEC implementation. Its public theory and documentation require
the following reusable ingredients to be explicit in FortFEM before an
external DESC-like client can be built:

1. **Nested inverse representation.** Store a toroidal-flux radial coordinate
   (the DESC reference convention is \(\rho=\sqrt{\psi/\psi_{\mathrm{LCFS}}}\),
   but the normalization must remain metadata), two periodic angles, and the
   surface embedding \((R,Z)\) plus stream/straight-field-line map \(\lambda\). Keep
   the axis, boundary, and field-period conditions as separate constrained
   traces.
2. **Axis-regular spectral spaces.** Provide Fourier modes coupled to a radial
   polynomial/spline basis with parity, regularity, and mode-coupling rules at
   the magnetic axis. Fourier--Zernike is the reference oracle; a compatible
   radial FE/IGA representation is an allowed FortFEM backend.
3. **Direct force residual.** Expose the radial and helical force-balance
   components at caller-owned collocation or quadrature points, including the
   volume Jacobian/metric weights, and provide assembled, matrix-free JVP,
   VJP, and Hessian-vector hooks. This is distinct from minimizing an energy.
4. **Equation and objective registry.** Treat force balance, boundary
   conditions, flux/rotational-transform constraints, profile callbacks,
   quasi-symmetry/omnigenity objectives, Mercier or ballooning diagnostics, and
   user objectives as typed residual/objective blocks. FortFEM supplies the
   derivative and constraint composition, not the plasma formulas.
5. **Continuation and perturbation.** Support implicit derivatives of a
   converged equilibrium, first/second/third-order perturbation blocks,
   pseudo-arclength continuation, trust-region/Newton acceptance, deflation,
   and topology/event diagnostics. Solver iteration paths remain outside the
   derivative contract.
6. **Geometry diagnostics and transforms.** Reserve neutral contracts for
   flux-surface averages, Boozer-like coordinate transforms, near-axis supplied
   data, magnetic-well/Mercier-style scalar diagnostics, and field-line
   samples. These must consume the common physical interchange schema and must
   not become DESC or VMEC file readers.
7. **Free-boundary ports.** Compose the same nested representation with
   external-field, vacuum, conductor, virtual-casing, PML, and explicit
   sheet-current traces. Fixed-boundary and free-boundary residuals must share
   shape/JVP/VJP conventions.

The 2026-08-02 foundation audit found that items 2--7 are currently contracts
or small algebraic fixtures rather than an executable DESC-class numerical
stack. The next implementation slices must therefore add the following
closure-neutral pieces (and keep their own manufactured oracles):

1. **Fourier--Zernike and vector parity basis.** Add deterministic
   \((l,m,n)\) enumeration, orthogonal/Jacobi or Zernike radial evaluation,
   conjugate-safe packing, symmetry reduction, and effective parity rules for
   vector and tensor components. The existing scalar axis-mode validator and
   finite power-sum evaluator are compatibility layers, not substitutes.
2. **Collocation/grid/transform layer.** Provide linear, quadrature, and
   concentric grids with axis and periodic-seam maps, spacing and quadrature
   weights, reshape/flatten and compress/expand operations, direct/FFT
   transforms, and bounded chunking controls. All transforms must have
   independent round-trip and adjoint oracles and must not allocate an
   unbounded full tensor product.
3. **Geometry differential jet.** Extend nested geometry beyond first
   Jacobians to inverse/contravariant bases, mixed \(\rho,\theta,\zeta\) Hessians
   and third derivatives where required, metric/volume derivatives, and
   safe axis-limit (L'Hôpital) evaluations. Keep fixed-topology JVP/VJP
   products explicit; do not hide singular limits in a client closure.
4. **Physical differential operator chain.** Build neutral metric-aware
   gradient, curl, divergence, Hodge, and force-sample composition from the
   geometry jet. The existing force-balance objective HVP is algebraic
   residual/weight differentiation only until a composed \(J^TWJ\) and
   state-residual HVP fixture is available.
5. **Typed objective/constraint derivatives.** Extend callback blocks with
   parameter/geometry tangents, complex packing, numeric targets/bounds/
   weights/scales, constraint Jacobian/KKT/nullspace metadata, and optional
   second-derivative/HVP and implicit-solved-state actions. Keep profile laws,
   readers, and production optimization external.
6. **Neutral diagnostic hooks.** Add multi-surface averages, Boozer-like
   transforms, supplied near-axis coefficient data, magnetic-well/Mercier-like
   reductions, rotational-transform/iota extraction, and field-line coordinate
   samples consuming only common physical fields and geometry. These are
   diagnostics, not DESC physics.
7. **Free-boundary current potential layer.** Add scalar surface-current
   potentials, harmonic loop/cut bases, period constraints, gauge/nullspace
   metadata, and shape derivatives on boundary-region graphs so NESTOR/BIEST-
   like source-to-trace maps can represent both local sheets and global loop
   currents.
8. **A bounded direct-force campaign.** Add axisymmetric and 3-D manufactured
   direct-force galleries with physical solution first, modal/grid convergence,
   geometry and objective derivative checks, and memory/chunking telemetry.
   This campaign must compose the operators above before any external DESC
   client is considered viable.

Since this audit, the neutral implementation now includes:

- `equation_objective_metadata_t` for targets, bounds, weights, scales,
  active/fixed flags, KKT/nullspace IDs, parameter tangents, and HVP capability
  metadata, with deep-copy and validation oracles;
- `evaluate_equation_objective_merit` and its JVP/VJP companions for the
  normalized weighted value algebra of one active metadata block, with
  independent finite-difference and real-adjoint checks. This remains a
  neutral residual/objective composition layer; bounds, profile laws, and
  optimizer decisions stay caller-owned;
- `evaluate_boozer_like_rotational_transform` and
  `evaluate_near_axis_diagnostic_metadata` for supplied-field/geometry samples
  and supplied near-axis coefficient data, with JVP/VJP actions and strict
  axis-regularity checks.
- `fortfem_nested_geometry_differential_jet` for caller-owned first-, second-,
  and third-order map jets, metric/determinant derivatives, inverse-Jacobian
  diagnostics, fixed-topology JVP/VJP actions, and explicit finite axis-limit
  handling. Its polynomial manufactured test passes the independent value,
  finite-difference, adjoint, singular-branch, and invalid-flag checks.
- `evaluate_surface_shape_objective` and its JVP/VJP companions for a
  positive-weighted mismatch of candidate and target physical surface samples.
  Point ordering, quadrature, topology, and geometry maps remain caller-owned;
  the contract is the neutral shape-mismatch port used by fixed/free-boundary
  clients rather than a DESC, VMEC, or plasma objective. Its pointwise weighted
  squared-mismatch value, JVP, and VJP come from the pinned FortSym generator
  `gen_surface_shape_objective_products`, and `check_generated.sh` compares the
  committed generated kernel byte-for-byte when the locked FortSym revision is
  available.
- The seeded `test_surface_shape_objective_properties` fixture samples 1-D,
  2-D, and 3-D fixed-topology surfaces, compares the weighted-loop value and
  product-rule JVP against independent arithmetic and centered differences,
  checks the real VJP identity, and exercises non-positive-weight and topology
  rejection. This adds randomized objective coverage without selecting a
  geometry or free-boundary physics model.
- `evaluate_surface_integral_constraint` and its JVP/VJP companions for a
  caller-owned scalar surface ledger with a target residual. Positive weights,
  samples, units, quadrature, topology, and target selection remain external;
  the primitive is reusable for area, flux, volume, loop, and shape constraints.
  Its pointwise weighted product, JVP, and VJP now come from the pinned FortSym
  generator `gen_surface_integral_products`; `check_generated.sh` compares the
  committed kernel byte-for-byte whenever the locked FortSym revision is
  available.

These are interchange and diagnostic foundations only. They do not implement
DESC profiles, optimization, equilibrium readers, coil models, or plasma
closures. The geometry differential jet, physical operator chain, and composed
direct-force campaign remain open below.

The checklist is grounded in the [DESC documentation](https://desc-docs.readthedocs.io/en/stable/),
including its [basis and collocation contract](https://desc-docs.readthedocs.io/en/stable/notebooks/basis_grid.html),
[transforms](https://desc-docs.readthedocs.io/en/stable/dev_guide/notebooks/transform.html),
[near-axis constraints](https://desc-docs.readthedocs.io/en/stable/notebooks/tutorials/nae_constraint.html),
[continuation](https://desc-docs.readthedocs.io/en/stable/notebooks/tutorials/continuation_step_by_step.html),
[deflation](https://desc-docs.readthedocs.io/en/stable/notebooks/tutorials/deflation.html), and
[free-boundary tutorial](https://desc-docs.readthedocs.io/en/stable/notebooks/tutorials/free_boundary_equilibrium.html),
the original [Dudt--Kolemen DESC formulation](https://doi.org/10.1063/5.0020743),
[Panici et al., Part I force-balance and convergence study](https://doi.org/10.1017/S0022377823000272),
[Conlin et al., Part II perturbation and continuation](https://doi.org/10.1017/S0022377823000399),
and the [high-order free-boundary formulation](https://arxiv.org/abs/2412.05680).
These are provenance and oracle references only; FortFEM still does not ship
DESC inputs, profiles, coil models, or optimization physics.

## 9. Solver, constraints, and differentiation roadmap

### 9.1 Sparse and nonlinear solvers

Complete the common solver product layer for direct solves, CG, PCG, GMRES,
BiCGSTAB, flexible Krylov, block systems, Schur complements, static
condensation, matrix-free actions, field splits, and retained factorizations.
Each converged-state solver has an implicit JVP/VJP. Iteration counts and
stopping decisions are inactive in that derivative contract.

For curl--curl and other gauge-singular systems, direct solves must use an
explicit tree--cotree or equivalent compatible nullspace reduction before
factorization. The selector is topology metadata and is not differentiated
through. Iterative paths must expose compatible nullspace projection and
measured ILU/ILUT, incomplete Cholesky (IC(0)), or incomplete Cholesky with
level/fill control (ICHOL) options where the matrix symmetry and definiteness
permit them. Complex or nonsymmetric blocks use an explicitly declared
replacement (for example ILU rather than ICHOL), with factor breakdown and
loss of positive definiteness reported rather than hidden.

The dense IC(0) factor/apply primitive and the PCG `ichol`/`ic0` option are now
public and independently tested. Standalone sparse IC(0) and ILU(0) paths now
consume FortSparse CSC directly, preserve the input lower/upper patterns
without fill, and expose fixed-factor right-hand-side JVP/VJP actions,
including the transpose solve for ILU. Production-size scaling remains
active solver work. The deterministic
`build_sparse_ilut`/`apply_sparse_ilut` path now supplies drop tolerance,
per-column fill selection, fixed-factor JVP/VJP, and explicit pivot status;
`build_sparse_ichol` supplies the matching SPD drop/fill path on the existing
sparse IC apply/JVP/VJP contract; `solve_sparse` now accepts the
`ichol_controlled` preconditioner with an exact full-fill one-step PCG oracle.
The public `build_sparse_ilut_row` constructor now performs no-pivot ILUT
against retained row factors with O(n + nnz) temporary work storage and an
independent no-fill diagonal oracle; it exports the same sparse CSC apply and
fixed-factor derivative contract as the dense reference builder. Controlled
row-oriented ICHOL now performs the matching symmetric row sweep with a
Cholesky diagonal oracle and the same sparse CSC apply/JVP/VJP contract.
The `solver_benchmark` gallery fixture now records row-ILUT and row-ICHOL
construction times on a shifted symmetric Poisson matrix and emits a log-log
timing plot; production-size memory and accuracy scaling remain. `solve_sparse`
PCG now accepts `sparse_ilut`/`ilut` and routes the legacy `ilu` alias through
the row-oriented factor directly, avoiding the old dense CSC-to-array
conversion. The sparse `ichol`/`ic`/`ic0` aliases likewise use the CSC
IC builder directly. Exact sparse oracles cover both integrations and their
bounded-memory intent. The
`solve_sparse` dispatcher now also has a true sparse GMRES callback path,
with left application of the selected bounded preconditioner, and a
matrix-free sparse BiCGSTAB callback path for ILUT-backed nonsymmetric blocks,
instead of silently falling back to PCG when either method is requested.
Independent exact CSC oracles cover both dispatches. The
neutral `solver_resource_budget_t` API validates caller-owned positive
wall-time, peak-memory, and repetition limits and evaluates measured usage
without introducing timing or allocation side effects; its independent
acceptance/over-budget/invalid-input test is the serial resource-boundary
contract for solver and benchmark runners. The
converged-state PCG JVP/VJP differentiates the exact solve independently of
the inactive preconditioner iteration path; factor rebuilds, breakdowns, and
graph changes are reported events rather than silently differentiated.

Nonlinear infrastructure includes Newton, damped Newton, Newton-Krylov,
pseudo-transient continuation, trust regions, deflation, parameter
continuation, bifurcation diagnostics, and explicit failure reasons.

Constraints include flux, helicity, mass, current, cross-helicity, angular
momentum, gauge, mode normalization, interface volume, and wall-current
conditions. Nullspace support covers constants, electromagnetic gauges, rigid
modes, harmonic representatives, and cohomology bases on multiply connected
domains.

### 9.2 FortSym and Enzyme pipeline

Every new formulation follows this sequence:

1. Define fields, spaces, geometry, and strong equations in FortSym.
2. Integrate by parts and expose volume, boundary, interface, and
   distributional terms.
3. Derive residual, Jacobian, JVP, VJP, and parameter or shape derivatives.
4. Simplify tensor symmetries and generate Fortran kernels.
5. Generate manufactured forcing, an independently simplified oracle, and
   equation documentation.
6. Check generated output byte for byte against the pinned source revision.
7. Use Enzyme as a second derivative route where compiler support is available.
8. Compare analytical derivatives with finite differences away from topology
   events, using step-size sweeps rather than one arbitrary step.

Generated kernels cover basis values, gradients, curls, divergences, Hessians,
Piola maps, traces, jumps, singular limits, cut-cell moments, Fourier triads,
special-function products, and geometry derivatives. Transient generated
files, plots, and compiler output are gitignored.

## 10. External application parity matrix

The following table defines the generic numerical foundation that an external
application could exercise. It does not authorize copying source code, adding
plasma-specific readers or closures, or shipping a replacement for the named
code. A row is complete when the generic contract and an independent oracle
work. Input conversion and application physics remain outside FortFEM.

| Reference target | Physics class | FortFEM ingredient parity target |
| --- | --- | --- |
| [CHEASE](https://crppwww.epfl.ch/~sauter/chease/), [paper](https://doi.org/10.1016/0010-4655(96)00046-X) | 2D fixed-boundary axisymmetric toroidal equilibrium | Generic axisymmetric elliptic/Fourier forms, spline/FEM geometry, axis and boundary trace contracts, and a common sampler. No COCOS or GEQDSK implementation |
| [FreeGS](https://freegs.readthedocs.io/en/stable/creating_equilibria.html) | 2D free-boundary axisymmetric equilibrium | Generic nonlinear residual, external boundary and coil-trace callbacks, X/O-point metadata fields, and manufactured profiles. Coil physics and GEQDSK conversion remain external |
| [VMEC/PARVMEC](https://github.com/ORNL-Fusion/PARVMEC), [VMEC++ numerics](https://arxiv.org/abs/2502.04374) | 3D nested-surface variational ideal equilibrium and free boundary | Kinematically nested toroidal embeddings, Fourier angles, radial FE/IGA, flux and constraint blocks, NESTOR-like vacuum traces, shape JVP/VJP, and an external-data sampler |
| [GVEC](https://gvec.readthedocs.io/develop/index.html), [DESC](https://github.com/PlasmaControl/DESC) | GVEC-like variational and DESC-like direct-force-balance 3D nested-surface equilibrium and optimization | General coordinate maps, radial B-splines or Fourier--Zernike bases, axis regularity, Fourier modes, multiple interfaces, collocation/weak residuals, exact residual derivatives, perturbation/continuation, objective/constraint callbacks, and free-boundary external-field/sheet-current traces. Input and profile models remain external |
| [SIESTA](https://doi.org/10.1063/1.3597155) | Eulerian ideal-MHD relaxation with islands and stochastic regions, including a free-plasma-boundary extension | Fixed-domain compatible volume fields, force and divergence residuals, relaxation/preconditioner contracts, topology events, and free-boundary trace coupling. Pressure and relaxation closures remain external |
| [SPEC](https://princetonuniversity.github.io/SPEC/) | Multi-region relaxed MHD with ideal interfaces | Region graph, independent fields, generic curl-eigenproblem and constraint blocks, total-pressure trace law, and interface shape derivatives. Beltrami and profile selection remain client code |
| [GPEC](https://princetonuniversity.github.io/GPEC/), [references](https://princetonuniversity.github.io/GPEC/references.html) | Linear ideal, kinetic, and resistive perturbed response with ideal outer shielding and optional resistive inner layers | Fourier coupling, explicit resonant sheet-current/penetrated-field traces, singular outer/inner layer contracts, vacuum and wall response, response matrices, normalization, and reciprocity. Equilibrium import remains external |
| [MARS-F response work](https://doi.org/10.1016/j.cpc.2006.09.003) | Linear toroidal ideal/resistive MHD and wall response | Linear block interfaces, complex frequency, generic plasma-vacuum-wall trace coupling, Fourier-FEM assembly, and resistive-layer matching. MARS physics remains external |
| [GLISS](https://github.com/itpplasma/GLISS) | Global linear ideal-MHD stability in 3D toroidal equilibria | Compatible radial spline FE, Fourier mode topology, generic eigenvalue and inertia contracts, and derivative boundaries. GVEC/VMEC input adapters remain external |
| [JOREK](https://www.jorek.eu/), [overview paper](https://arxiv.org/abs/2011.09120) | Nonlinear extended and resistive MHD | 2D FE plus toroidal Fourier blocks, coupled residuals, anisotropic transport, implicit structure-aware stepping, wall and free-boundary traces, operator-level parity tests |
| MEPHIT and STARWALL | Electromagnetic response and resistive-wall coupling | H(curl)/H(div) FEEC, surface traces, FEM/BEM/DtN wall blocks, retained complex factors, interface currents, reciprocity and energy tests |
| [Bíró, Preis, and Richter tree--cotree formulation](https://ieeexplore.ieee.org/document/558631) | Direct-gauge curl--curl finite-element benchmark for multiply connected conducting domains | A license-safe manufactured reproduction of the paper's topology, tree/cotree reduction, direct solve, gauge-invariant magnetic field, and energy/loop-current diagnostics. No paper code or proprietary geometry is copied |
| [TEAM workshops](https://www.osti.gov/servlets/purl/7179128), [TEAM 13 reference reproduction](https://docs.feelpp.org/toolboxes/latest/maxwell/Tws/index.html) | Community electromagnetic magnetostatic and eddy-current benchmark problems | Neutral supplied geometry/material/source arrays, H(curl) residuals, gauge/direct-solver paths, probe/energy/force outputs, and external reference-data manifests. TEAM readers, nonlinear material laws, and benchmark data licensing remain external |
| [BOUT++](https://bout-dev.readthedocs.io/en/stable/user_docs/introduction.html) | General 3D curvilinear fluid framework with model-specific clients | Equation-as-data residuals, curvilinear metric and boundary contracts, field-aligned operators, mixed conservative fluxes, and model-level JVP/VJP. Fluid models remain external |
| [GRILLIX](https://physik.uni-greifswald.de/ag-manz/forschung/codes/grillix/), [FCI paper](https://doi.org/10.1088/1361-6587/aaa373) | 3D edge and scrape-off-layer use of flux-coordinate-independent operators | FCI field-line tracing and interpolation, parallel/perpendicular operator split, immersed boundaries, anisotropic transport, generic material traces, and manufactured MMS fixtures. Braginskii and sheath laws remain external |
| [PARALLAX](https://gitlab.mpcdf.mpg.de/phoenix-public/parallax), [elliptic solver paper](https://arxiv.org/abs/2509.11831) | FCI mesh, magnetic-field handling, 2D elliptic solves, matrix-free 3D actions, and multigrid for GRILLIX and GENE-X | A Fortran-compatible geometry and operator adapter, matrix-free sparse products, anisotropy-aware multigrid contracts, and independent Poisson/Ampere timings. PARALLAX is LGPL-3.0, so no source is copied into FortFEM. |
| Linear elasticity and wave FEM literature | Mixed stress-displacement elasticity and symplectic mixed acoustics | Elasticity-complex spaces, tensor pressure/stress, first-order wave state, port-Hamiltonian pairing, symplectic time step, and cross-physics manufactured tests |

The local `/home/ert/code/GLISS` checkout was inspected for this roadmap. Its
README identifies version 0.0.2 as MIT licensed, with fixed-boundary FEEC
spectra, Mercier diagnostics, symmetric GVEC/VMEC input, radial B-splines,
Fourier angular modes, and Enzyme-generated derivative actions. Its
`PROVENANCE.md` records source-level traceability and a compatibility policy.
The planned FortFEM integration is therefore an oracle and interchange layer,
not a source import.

The name “Brandon Shennahan” did not resolve to a unique plasma-code project
in the local checkouts or the public search results. The relevant public
lineage is nevertheless clear enough to track: BOUT++ provides a general
curvilinear fluid framework, GRILLIX provides flux-coordinate-independent
edge/SOL turbulence, and PARALLAX provides the shared FCI geometry and
elliptic-solver layer. If a specific Shennahan project is intended, add its
repository or paper to this table before claiming parity.

## 11. Ordered example gallery

Examples are generated from one command, run in increasing complexity, and
publish plots plus machine-readable CSV or JSON. Images are build artifacts,
not committed source files. Every example includes at least a field plot, a
mesh or patch plot, an error or convergence plot, and a derivative or
conservation diagnostic when applicable.

### Foundations

1. **1D Poisson and Sturm--Liouville.** P1, high-order, exact polynomial
   reproduction, JVP/VJP, and a simple line plot.
2. **1D Legendre and spherical functions.** Ordinary, associated, and
   half-integer functions with recurrence and differential-equation oracles.
3. **2D simple Poisson.** First gallery entry, with mesh, solution, error,
   and timing plots.
4. **2D anisotropic diffusion.** Tensor coefficients, flux conservation, and
   manufactured solution.
5. **2D FEM/BEM Laplace.** Circle and polygon transmission with a larger-domain
   comparison.

### Interfaces, vectors, and open boundaries

6. **Curved Poisson interface.** Fitted CG, XFEM, XIGA, Nitsche, SIPG, and a
   regularized interface compared against a known jump solution.
7. **Surface delta source.** Explicit interface integral, duplicated traces,
   and narrow regularized source.
8. **H(curl) Ampere jump.** Analytical tangential magnetic jump and explicit
   surface current \(\mathbf K\), with Nitsche and fitted alternatives.
9. **H(div) magnetic flux.** Exact normal-flux continuity and divergence
   diagnostics on a fitted and cut interface.
10. **Slab Helmholtz.** FEM, BEM, DtN, and PML against an analytical plane wave.
11. **Cylinder and sphere Helmholtz.** Bessel/Hankel and Mie oracles, DtN, PML,
    larger-domain comparison, and 1D radial plus 2D/3D plots.
12. **Torus scalar exterior problem.** Exact curved toroidal surface, torus
    harmonics, FEM/BEM, DtN, PML, and far-field error.
13. **Curl-curl torus scattering.** Nédélec FEM, RWG/BC BEM traces, vector DtN,
    PML, reciprocity, and Ampere performance data.

### Application-oriented numerical ingredients

14. **Toroidal Maxwell FEM--BEM manufactured solution.** An affine physical
    H(curl) field, exact edge-integral coefficients, curved RWG trace coupling,
    solved-state oracle, physical cross-section, vector arrows, and curved
    three-dimensional geometry. This is the neutral vector FEM--BEM baseline;
    nonzero scattering remains in the curl-curl fixture above.
15. **Cylindrical axisymmetric elliptic fixture.** Manufactured coefficients,
    axis regularity, Fourier-FEM, and optional CHEASE/FreeGS comparison data
    supplied through an external adapter.
16. **Fourier-FEM slab and cylinder.** Mode diagonal operators, retained
    nonlinear triads, real/conjugate packing, and torus-harmonic diagnostics.
17. **Bíró tree--cotree curl--curl benchmark.** Reproduce the published
    direct-gauge topology on a small multiply connected conducting domain with
    the same tree/cotree reduction, a gauge-invariant \(\mathbf B\) field, loop
    current, energy, and direct-solve residual plots. The fixture is
    manufactured and provenance-pinned; it demonstrates the method without
    shipping paper source or application-specific readers.
18. **TEAM electromagnetic benchmark ladder.** Start with a small
    license-safe TEAM magnetostatic/eddy-current subset (TEAM 3, 7, 13, or 20
    as the external manifest permits), then add supplied probe curves and
    field/eddy-current/energy/force plots. Every case records geometry and
    material-data provenance, direct-gauge versus constrained solve, mesh
    convergence, and a solution-first 2D/3D view. Full benchmark data stays in
    the sister data repository when redistribution is restricted.
19. **Multi-region curl-eigenproblem fixture.** Independent regions, ideal
    interfaces, generic flux and helicity constraints, and a pressure-balance
    residual with manufactured coefficients.
20. **Linear 3D perturbation blocks.** Generic mode response with vacuum, wall,
    singular-layer, and response-matrix toy operators. GPEC or MARS data enter
    only through an external sampler.
21. **Energy-principle toy spectrum.** Radial spline FE, Fourier modes, inertia
    count, eigenpair derivatives, and external-data interchange tests.
22. **Resonant interface sheet.** Ideal current-sheet limit, finite resistive
    layer, XFEM enrichment, fitted sheet, and convergence to a manufactured
    singular solution.
23. **Skew bracket fixture.** Energy-skew nonlinear bracket, Fourier convolution,
    analytical JVP, and long-time structure-preserving integration for a
    caller-defined state.
24. **Resistive diffusion and tearing layer.** Explicit layer, adaptive layer,
    asymptotic enrichment, DG, and ideal-limit comparison.
25. **Generic coupled-field path.** Independently testable magnetic, scalar,
    vector, tensor, interface, and constraint residual blocks. A JOREK-style
    client can map its fields into this path, but FortFEM contains no JOREK
    state or closure implementation.

### Waves, elasticity, and anisotropy

26. **Mixed symplectic acoustics.** Pressure and particle velocity in a first-
    order compatible pair, with energy-preserving symplectic stepping and a
    comparison against the second-order pressure reduction.
27. **General wave family.** The same mixed residual for acoustics, Maxwell,
    elastodynamics, and a scalar wave, with common energy, dispersion, and
    boundary-port plots.
28. **Structure-preserving linear elasticity.** Displacement, velocity, and
    tensor stress with weak symmetry, traction interfaces, and a mixed
    Hellinger--Reissner oracle.
29. **Tensor-pressure wave.** A caller-supplied anisotropic tensor with
    parallel, perpendicular, and gyrotropic parts, including force balance and
    energy diagnostics. The tensor is a generic constitutive fixture.
The current `cgl_pressure_tensor` gallery example provides the first
manufactured constitutive/force-divergence profile and CSV/1D FortPlot
outputs; the coupled wave and higher-dimensional cases remain active.
The `field_aligned_flux` gallery example now provides a generated
parallel/perpendicular profile and \(k_\parallel/k_\perp=100\) flux plot. The
`anisotropic_tensor_diffusion` gallery fixture now adds a physical 2D P1
solution for a \(k_\parallel/k_\perp=1000\) tensor, with a solution contour as
its first plot, a centerline oracle, mesh view, error/energy, and timing data.
Three-dimensional curvilinear and FCI/Fourier parity cases remain active.
28. **Field-aligned diffusion.** A slab, cylinder, and torus with extreme
    \(k_\parallel/k_\perp\), comparing aligned coordinates, FCI maps, Fourier-
    FEM, and an isotropic control case.
29. **Field-aligned edge operator.** A generic anisotropic transport system
    with caller-supplied coefficients, material traces, and a reproducible FCI
    field-line map. This is an operator fixture, not a GRILLIX, BOUT++, or
    Braginskii implementation.

The gallery must show the same case in 1D, 2D, and 3D where a dimensional
reduction exists. Plot scripts use FortPlot and must include mesh edges,
element orientation, patch boundaries, and internal surfaces without dropping
parts of a mesh.

## 12. Verification and benchmark hierarchy

### 12.1 FortFEM-internal levels

1. **Algebraic:** partition of unity, polynomial reproduction, orientation,
   trace signs, exact sequences, rank, nullspaces, and surface measures.
2. **Analytical patch:** constants, affine fields, exact jumps, singular
   enrichments, special functions, torus harmonics, and exact delta integrals.
3. **Manufactured:** FortSym-generated source, boundary data, interface data,
   residual, Jacobian, JVP, VJP, shape derivative, and convergence rates.
4. **Cross-formulation:** fitted CG, XFEM/XIGA, DG/HDG, Nitsche, mortar,
   explicit surface current, regularized layer, FEM/BEM, DtN, and PML on the
   same physical case.
5. **Conservation and structure:** energy, helicity, cross-helicity, flux,
   divergence, charge, reciprocity, passivity, symplectic form, reversibility,
   and dissipative balance.
6. **Performance:** assembly, factorization, solve, matrix-free action, memory,
   iteration count, conditioning, and derivative cost as functions of order,
   mesh size, mode count, and interface complexity.

Report \(L^2\), H1, H(curl), H(div), trace, jump, flux, energy, conservation,
and derivative errors. Report error versus degrees of freedom and versus wall
time. Use a controlled thread count and record compiler, flags, CPU, BLAS,
FortNum/FortSparse revision, and random seed.

### 12.2 External oracle policy

The external comparison matrix covers FreeFEM, MFEM, FEniCSx, deal.II,
GeoPDEs, PetIGA, CHEASE, FreeGS, VMEC/PARVMEC, GVEC, DESC, SPEC, GPEC,
MARS-F, GLISS, STARWALL, and JOREK where a license-safe fixture is available.

- Keep only tiny analytical inputs, adapters, expected tolerances, and a
  provenance manifest in FortFEM.
- Run lightweight FreeFEM, MFEM, and FEniCSx cases in optional jobs. Use `uv`
  for Python environments where practical.
- Do not require heavy, proprietary, or cluster-only packages in the GitHub
  Pages job. Their pre-generated numerical data and performance summaries go
  to a separate `fortfem-benchmarks` data repository, referenced by commit,
  checksum, license, and executable version.
- Never copy external source code or proprietary test data into FortFEM.
- Compare fields and invariants on a common physical sampling set, not by
  assuming that two codes use the same basis or numbering.

## 13. Performance and reproducibility

The performance path is Fortran, FortNum, and FortSparse only. Profile with
the platform tools before changing an algorithm. Track:

- element and cut-cell kernel throughput;
- Fourier transform and mode-convolution cost;
- sparse assembly and matrix-free products;
- direct factorization reuse;
- Krylov iterations and preconditioner setup;
- BEM near-field, far-field, and quadrature cost;
- PML layer overhead;
- JVP/VJP cost relative to primal;
- memory peak and allocation counts.

Small deterministic benchmarks run in focused tests. Larger benchmarks run in
scheduled or manually dispatched Actions, publish JSON and SVG/PNG artifacts,
and update the Pages gallery without committing images. A benchmark result is
not considered comparable unless hardware, compiler, thread count, source
revision, and external data revision are recorded.

Performance gates are workload-specific. The benchmark manifest records the
operator, discretization, order, mesh or mode count, tolerance, thread count,
and reference implementation. It reports time to assemble, time to apply or
solve, peak memory, allocation count, iteration count, derivative overhead,
and error at a common physical sampling set. Reference runs may use MFEM,
Firedrake, NGSolve, deal.II, PETSc, Bempp-cl, or a small analytical kernel,
depending on the operator. The target is a measured improvement on the
declared workload, with accuracy and structure-preservation constraints held
fixed. The project does not use an unqualified speed claim as a substitute for
the benchmark record.

The execution plan is staged:

1. Keep all element and residual kernels local, pure where possible, and
   independent of MPI. Use batched quadrature, generated contractions,
   static condensation, matrix-free actions, and reusable sparse factors on a
   single node.
2. Add the serial ownership and ghost interfaces, deterministic reductions,
   partition-independent numbering, and a no-op communicator backend to the
   core data structures.
3. Add optional MPI halo exchange, distributed assembly, distributed vectors
   and matrices, checkpoint/restart, and solver adapters after the serial
   contracts have independent analytical and cross-code oracles.
4. Add distributed BEM acceleration, domain decomposition, GPU kernels, and
   large-scale strong and weak scaling only when a representative workload
   demonstrates the need.

The single-node path remains a complete supported path at every stage. The
MPI path must reproduce the serial residual, derivative, conservation, and
structure diagnostics within declared reduction tolerances.

## 14. Implementation phases

The phases are ordered by dependency. Each phase ends with a public API,
generated or analytical oracle, focused tests, documentation, and at least one
gallery example.

### Phase 0: Contracts and inventory: **active**

- Complete the primal/JVP/VJP inventory for FEM, BEM, DtN, PML, geometry,
  Fourier, special functions, sparse products, and all iterative solvers.
- GMRES and BiCGSTAB implicit state derivatives now pass finite-difference and
  adjoint-identity tests. Continue the inventory for remaining public
  operators and block solver compositions.
- Publish the complex-adjoint and shape-derivative conventions.
- The `check` helper used by the normal `fo test` path now provides explicit
  per-case xorshift seeds, bounded real/integer generators, callback-based
  property execution, and deterministic integer shrinking with failure-seed
  reporting. The first property test covers reproducibility and bounds;
  the Fourier mode-set property now exercises 20 generated conjugate-packed
  registries and one-product padded closure with the same helper. Geometry,
  tensor, and fixed-topology interface generators can build on this path
  without global random state.
- The public `partition_layout_t` now defines the communicator-free
  partitionable data boundary: unique local-to-global IDs, owner rank,
  explicit owned/ghost masks, deterministic local reductions, and linear JVP/
  VJP actions. It is the serial no-op implementation that future halo and
  all-reduce backends can replace without changing local kernels. MPI
  assembly, checkpointing, and distributed solver adapters remain deferred
  until the serial residual and invariant contracts are stable.
- The neutral `solver_resource_budget_t` contract now validates positive
  caller-owned wall-time, peak-memory, and repetition limits, and reports
  whether measured usage fits those limits without timing, allocation, MPI,
  or solver side effects. A seven-case independent test covers acceptance,
  over-budget measurements, and invalid limits; this is the resource gate
  used by focused solver and gallery runners before larger benchmarks.
- Keep FortSym revision pins and generated-kernel checks green.

### Phase 1: Interface calculus: **active**

- The public orientation-safe scalar/vector trace contract now provides
  plus/minus averages and jumps, normal/tangential projections, and the
  rotated Ampere surface-current jump with an independent sign oracle.
- The explicit `assemble_surface_delta_load` primitive now assembles
  trace-basis transpose times positive surface quadrature/source weights,
  providing the fitted δ_\Gamma weak-load contract.
- `assemble_surface_vector_delta_load` adds the tangential trace/surface-
  current pairing needed for an explicit Ampère sheet. Scalar and vector
  delta-load actions now also expose analytic fixed-topology JVP/VJP products
  for basis, quadrature weights, and source/current values; this keeps explicit
  sheet residuals differentiable for fitted, cut, DG, and IGA callers.
- `assemble_interface_surface_current` now evaluates the oriented Ampere
  jump, its integrated current ledger, and fixed-topology JVP/VJP actions.
  Independent analytical, orientation-reversal, finite-difference, and
  real-adjoint tests cover the generic trace algebra; conservation at
  interface edges and material laws remain higher-level contracts.
- `evaluate_regularized_surface_current_layer` now maps caller-owned signed
  normal distances and tangential sheet currents to the normalized Gaussian
  approximation of \(\mathbf K\delta_\Gamma\). FortSym-generated value/JVP/VJP
  products cover distance, current, and positive thickness, while a separate
  quadrature diagnostic reports profile normalization and integrated current.
  Geometry-specific distance construction and resistive-layer closures remain
  client-owned.
- `assemble_surface_current_trace_residual` now pairs caller-owned,
  independently sized test and trial tangential trace bases with a target
  current. Its full product-rule JVP and real VJP cover basis, quadrature,
  coefficients, and target data, with a direct vector oracle. This is the
  neutral ownership block for fitted duplicated traces, cut/XFEM or XIGA,
  DG/HDG skeletons, and IGA patches; constitutive pressure laws and
  flux/helicity constraints remain client composition.
- `assemble_interface_jump_penalty` assembles the symmetric positive-
  semidefinite plus/minus jump block used by SIPG and Nitsche penalty terms.
- `assemble_symmetric_nitsche_interface` adds the average-flux consistency
  terms under the same orientation convention. Its value/JVP/VJP actions
  differentiate traces, fluxes, surface weights, and penalty parameters;
  scalar SIPG and non-symmetric consistency variants remain.
- `assemble_mortar_trace_coupling` supplies a weighted cross-mass block for
  independently discretized trace spaces.
  `assemble_scalar_sipg_interface` now composes independently sized test and
  trial trace/flux spaces with symmetric, incomplete, or nonsymmetric
  consistency (`theta=1,0,-1`) and value/JVP/VJP actions; vector FEEC
  consistency, HDG trace blocks, local static condensation, and signed-map
  sparse CSC accumulation are now public. The compatible flux elimination
  primitive now provides a differentiable local recovery map and condensed
  state block; global signed sparse flux ownership remains client-owned.
  `assemble_vector_jump_penalty` now supplies the tensor-valued counterpart
  for caller-owned tangential/normal projectors and anisotropic metrics, with
  value/JVP/VJP actions for vector traces, metrics, weights, and penalties.
  `assemble_vector_sipg_interface` now adds independent test/trial vector
  average-flux consistency for `theta=1,0,-1`, including metric-aware
  value/JVP/VJP actions; higher-order HDG trace composition and concrete
  enriched-space construction remain.
  `assemble_hdg_static_condensation` now exposes the local Schur complement
  and condensed load with implicit-solve value/JVP/VJP actions, so global
  skeleton assembly can be layered without differentiating pivot choices.
  `assemble_compatible_flux_elimination` now exposes the corresponding
  flux-specific recovery map, condensed state block, and full real JVP/VJP;
  `assemble_feec_commuting_projection` now audits all three projected
  differentials with complete value/JVP/VJP actions; actual higher-order
  enriched-space construction and client-owned projection maps remain.
  `assemble_scalar_numerical_flux` now provides conservative central, upwind,
  and Lax--Friedrichs choices with a quadratic-entropy diagnostic and complete
  fixed-topology value/JVP/VJP actions; `assemble_vector_numerical_flux` now
  lifts this to metric-weighted vector states, including strongly anisotropic
  constitutive tensors and complete value/JVP/VJP actions. Entropy-stable
  systems remain.
- The neutral `cell_complex_t` contract now stores oriented integer chain
  boundary maps, checks both boundary-of-boundary identities exactly, and
  reports Euler characteristic and compact Betti diagnostics. Independent
  interval, loop, sphere-CW, torus-CW, and malformed-orientation tests cover
  the primitive; quotient boundary maps, cycle/cocycle kernels, and metric
  harmonic one-forms are now public. `normalize_harmonic_one_forms` now
  supplies fixed-topology period/flux normalization with a central-difference
  JVP and real-adjoint oracle; cycle labels, physical units, gauge constraints,
  and cycle ledgers remain higher graph layers.
- The neutral `region_interface_graph_t` contract now adds oriented plus/minus
  region incidence, periodic self-identifications, and compact connectivity
  labels plus a spanning-forest cycle basis satisfying the exact integer
  incidence nullspace. Independent chain, disconnected, reversed, periodic,
  triangle-cycle, and malformed-endpoint tests cover the graph. Surface laws,
  sheet-current balances, and region physics remain application-owned.
- The neutral `cell_identification_t` contract now validates idempotent
  canonical representatives and signed orientations and reports compact
  quotient classes. `identify_boundary_matrix` pushes oriented incidence to
  quotient classes and rejects inconsistent identified columns. Independent
  identity, signed-periodic, interval-to-circle, signed-column,
  inconsistency-rejection, cycle-rejection, canonical-sign, zero-sign, and
  shape-mismatch tests cover the metadata; quotient geometry remains
  higher-level work.
- The neutral tree--cotree gauge contract now selects a spanning forest on an
  oriented graph, exposes cotree restriction/prolongation, extracts the
  fixed-gauge dense direct system, and composes with real/complex CSC direct
  solves through the existing constrained reduction. Independent triangle,
  disconnected-forest, reduction, sparse solve, fixed-map derivative, and
  malformed-incidence tests cover the selector. `build_tree_cotree_dof_map`
  now lifts the frozen control-edge selector to arbitrary high-order or IGA
  global DOF numbering, retaining extra edge/face/cell moments for the caller's
  sparse direct solver; duplicate and out-of-range maps have independent
  rejection tests. The new complex cycle-period constraint composes this
  fixed gauge with caller-owned oriented periods and has independent
  finite-difference and real complex-adjoint tests. Deterministic connected
  component labels and explicit forest/cycle-rank checks are also public, so
  disconnected meshes and future IGA control graphs do not rely on implicit
  component assumptions.
- The neutral `internal_manifold_graph_t` contract now records oriented
  plus/minus region sides, open or boundaryless manifold endpoints, periodic
  self-identifications, junction incidence, closed-manifold flags, and
  manifold connectivity. Independent slab, cylinder, sphere, torus, and
  malformed-endpoint fixtures cover the metadata; geometry, trace spaces,
  surface divergence, sheet-current balance, and pressure laws remain later
  composition layers.
- The integrated surface-current junction ledger now consumes the fixed
  manifold boundary incidence and distributes each manifold current to open
  junctions, with a zero global balance oracle for conservative columns and
  fixed-topology JVP/VJP actions. Geometric edge flux, constitutive current
  laws, and physical pressure/flux laws remain application composition.
- The fixed-topology loop-current constraint now applies an integer cycle basis
  to integrated manifold currents and subtracts caller-owned target values.
  Its residual, JVP, and VJP have an independent cycle oracle. Open-edge
  junction balance and closed-loop constraints are therefore separate linear
  ledgers; geometric surface divergence, constitutive current laws, and
  physical pressure/flux laws remain composition layers.
- The topology-only surface edge-flux ledger now applies an oriented
  vertex-edge boundary map to caller-integrated tangential edge fluxes. It
  exposes vertex divergence, a global conservation scalar, and fixed-topology
  JVP/VJP actions with open-chain and closed-cycle oracles. Conormal geometry,
  edge quadrature, and the physical tangential-current law remain composition
  layers. `assemble_surface_edge_flux` now supplies the missing neutral
  geometry-to-edge contraction from oriented conormal quadrature and surface
  current, with product-rule JVP/VJP actions and independent orientation,
  finite-difference, and adjoint tests. Surface-basis ownership and physical
  current laws remain composition layers; the neutral test/trial trace
  residual is now public for that ownership layer.
- The neutral normal-traction jump now projects caller-supplied plus/minus
  traction vectors onto a validated unit normal and subtracts a caller-owned
  target. Its product-rule JVP/VJP oracle composes generated CGL, elastic, or
  Maxwell-stress blocks without selecting a constitutive or plasma law.
- The neutral full-vector traction jump now subtracts a caller-owned target
  vector from plus/minus tractions, with analytical JVP/VJP products and an
  independent component and real-adjoint oracle. Normal, tangential, and
  full-vector interface laws remain explicit client compositions.
- Oriented triangle surface measures (area plus unit normal) now have a public
  JVP/VJP API with shared-vertex accumulation and independent finite-difference
  and dot-product oracles. A linear 2D triangle level-set cut primitive now
  returns edge intersections, physical segment length, and gradient normal with
  an affine independent oracle. Exact positive/negative subcell areas and
  interface-length consistency are now available for the same linear cut. A
  degree-one clipped-polygon quadrature primitive now adds exact positive and
  negative centroids, oriented normal data, and an affine-integration oracle.
  The fixed-topology 2D level-set interface JVP now differentiates edge
  intersections, segment length, and the normalized physical normal with an
  independent central-difference and topology-event oracle. The matching
  fixed-topology cut-quadrature JVP propagates edge intersections through
  positive/negative areas and centroids and composes the length/normal
  derivative, with central-difference and area/first-moment conservation
  oracles. A 3D tetrahedral level-set interface now returns ordered triangular
  or quadrilateral cut polygons, area, and gradient-oriented normal with
  independent plane/intersection tests. Exact positive/negative tetra cut
  volumes and centroids now close the degree-one volume/first-moment contract
  with analytic and conservation oracles. Fixed-topology tetra cut JVPs now
  propagate moving vertices and level values through clipped-face moments and
  interface area/normal with central-difference and conservation oracles. The
  internal-manifold graph now provides the topology attachment for this
  surface-measure contract; geometry-to-graph construction remains client
  composition work. The existing vector current pairing consumes the
  surface-measure contract.
- The public `broken_space_layout_t` and `skeleton_space_layout_t` now provide
  deterministic cell/facet ownership maps for broken H1, H(curl), H(div), and
  L2 spaces plus scalar, tangential, and normal skeleton traces. Frozen
  edge/face signs and explicit zero-sign boundary sides are validated before
  maps are exported to caller-owned HDG, DG, Nitsche, cut, or IGA assembly.
  Basis evaluation, metric maps, and physical interface laws remain separate.
- Explicit delta-source and surface-current weak terms.
- The public `assemble_fitted_trace_constraint` block now assembles the
  oriented multiplier pairing `\int_\Gamma lambda (u^+ - u^-) dS` for
  independently sized plus/minus and multiplier traces, with complete
  value/JVP/VJP actions and positive-weight validation. Fitted duplicated
  space numbering, constitutive jump targets, and global block constraints
  remain caller composition.
- Build the region and cell-complex graph around the validated chain maps,
  adding periodic identifications, deterministic homology/cohomology
  representatives, harmonic representatives, gauge constraints, and stable
  cycle IDs. The public homology and cohomology quotient bases remove
  face-boundary and exact-cochain components by scale-aware rank-increase
  oracles; metric and period normalization remain separate layers.
- Add closed-loop and open-edge sheet-current constraints, surface divergence,
  pressure balance, and current-ledger oracles on slab, cylinder, sphere, and
  torus fixtures.

### Phase 2: Cut geometry and XFEM/XIGA: **active**

- The shifted Heaviside partition-of-unity primitive is now public. It
  returns the sign-shifted enrichment, has fixed-sign zero JVP/VJP actions,
  and rejects a zero level value as a topology event. Independent sign and
  derivative oracles cover the piecewise-smooth contract. The matching
  scalar product composition `N_i*(H(phi)-H(phi_i))` now has value/JVP/VJP
  actions and a real adjoint oracle. The same activation is now public for
  componentwise vector bases, so H(curl)/H(div) values can be composed while
  retaining the fixed-topology JVP/VJP contract. The corrected-XFEM blending
  ramp `r=sum_i a_i*N_i` and `r*Psi` value/JVP/VJP composition are now public
  with a full-enrichment reproduction oracle. The same ramp and reverse
  contraction are public for vector enrichment arrays, ready for covariant or
  contravariant Piola values and IGA coefficient blocks. Cut-cell geometry,
  support activation, rank/conditioning diagnostics, and Piola-aware vector
  differential operators remain.
- The matrix-level `evaluate_shifted_enriched_space` constructor now builds
  all scalar shifted-Heaviside basis values at once from a base basis matrix,
  quadrature level-set values, and per-basis anchors. Its independent sign,
  finite-difference, adjoint, and topology-event tests make this the first
  actual scalar unfitted-space composition; vector FEEC/XIGA constructors,
  cut integration, and global sparse numbering remain separate layers.
- The matching `evaluate_shifted_vector_enriched_space` constructor now builds
  a component-by-basis-by-point matrix with one shared cut sign per basis
  function. Its independent componentwise sign, finite-difference, adjoint,
  and topology-event tests cover the vector composition before Piola maps;
  exact-sequence policy, vector cut integration, and global numbering remain
  explicit caller layers.
- `evaluate_vector_enrichment_differential_3d` now exposes the explicit
  product-rule terms `curl(psi*b)=psi*curl(b)+grad(psi) x b` and
  `div(psi*b)=psi*div(b)+grad(psi) . b`, with full input JVP/VJP actions and
  independent curl/divergence and adjoint oracles. It is a neutral diagnostic
  for H(curl), H(div), Piola, and IGA composition; exact-sequence preservation
  or intentional jump reporting remains a higher-level space contract.
- The public `evaluate_batched_vector_enrichment_differential_3d` primitive
  now composes all vector basis functions and quadrature points with the exact
  product rules for enriched values, curl, and divergence. Its full JVP/VJP
  actions have independent central-difference and real-adjoint oracles. This
  closes the batched vector differential layer before Piola-aware cut
  integration; global numbering, continuity, and exact-sequence policy remain
  caller-owned.
- Cut-cell classification and high-order quadrature. Exact degree-one triangle
  and tetrahedron rules plus exact degree-two raw-moment tensors with
  fixed-topology JVPs are now public. The 2-D linear level-set path now also
  exposes exact degree-three raw-moment tensors and fixed-topology JVPs for
  positive and negative clipped polygons, with an independent simplex-moment
  and central-difference oracle. The matching 3-D tetrahedral path now exposes
  exact degree-three raw-moment tensors and fixed-topology JVPs through its
  oriented tetrahedral fan, with an independent simplex-moment and
  central-difference oracle. The 2-D triangle path now also exposes an exact
  degree-four raw-moment tensor and fixed-topology JVP through an exact
  Green-theorem binomial edge primitive, with independent simplex-moment,
  conservation, and central-difference oracles. The 3-D tetrahedral path now
  also exposes an exact degree-four raw-moment tensor and fixed-topology JVP
  through its oriented fan and barycentric multinomial primitive, with the
  same independent simplex-moment, conservation, and central-difference
  oracles; curved-cell moment-fitting rules remain. A weighted
  enrichment-support Gram contract now exposes value/JVP/VJP actions with a
  fixed activation mask, and a symmetric-Jacobi rank/condition diagnostic
  reports the active enrichment rank and singularity with an independent
  weighted, finite-difference, and real-adjoint oracle; support activation maps
  now expose CSR sign activation, unique extrema, signed margins, and
  fixed-topology JVP/VJP actions with explicit owner/tie/topology rejection;
  the physical vector/metric support Gram contraction now composes with
  covariant or contravariant Piola and IGA values and exposes value/JVP/VJP
  actions with SPD-metric rejection. The batched 2D/3D covariant and
  contravariant Piola-enrichment value/JVP/VJP composition now uses FortNum’s
  FortSym-derived inverse/determinant products and rejects singular or
  orientation-reversing maps. The affine differential companions now
  reports the covariant H(curl) curl and contravariant H(div) divergence
  product terms in 2D and 3D with geometry/enrichment-gradient value/JVP/VJP
  actions;
  exact-sequence cut spaces, higher-order curved maps, and commuting
  conditioning remain.
- Heaviside, kink, singular, helical, and resonant enrichments.
- Shifted bases, corrected blending elements, pruning, conditioning, and
  connected-component activation.
- Trimmed B-spline stabilization and cut H(curl)/H(div) extensions.
- Verify the commuting diagram for every enriched vector space and document
  the exact sequence identity that a physical jump intentionally changes.

### Phase 3: DG and HDG: **active**

- The public mortar trace cross-mass block now has value/JVP/VJP actions for
  independently discretized skeleton traces. The scalar SIPG interface block
  now has symmetric, incomplete, and nonsymmetric consistency variants with
  value/JVP/VJP actions; vector FEEC consistency, HDG traces, local
  hybridization, signed-map sparse CSC accumulation, and compatible flux
  elimination are now public. The latter returns the local compatible flux
  recovery map plus condensed state system with full real JVP/VJP actions;
  global sparse flux ownership remains client-owned.
  The geometry-aware `assemble_geometry_mortar_trace_coupling` wrapper now
  forms physical quadrature weights from caller-owned reference weights and
  surface metrics, returns that measure ledger, and differentiates traces,
  metric, and weights. NURBS/Fourier/panel/cut geometry samplers remain
  caller-owned; coordinates, normals, and distributed ownership are not
  inferred by this layer.
  A tensor-weighted vector jump penalty now covers component-valued traces and
  anisotropic metric directions; vector consistency fluxes and FEEC commuting
  projections are now public with value/JVP/VJP actions; hybridization and
  static condensation now have a local differentiable Schur primitive;
  `assemble_hdg_global_skeleton` now supplies a signed-map dense reference
  assembler with value/JVP/VJP actions; sparse accumulation remains.
  `assemble_feec_exact_sequence` now exposes the
  metric-independent `curl(grad)` and `div(curl)` defects with independent
  value/JVP/VJP actions for simplicial, IGA, multipatch, and periodic maps.
  `assemble_feec_commuting_projection` now checks projected discrete and
  continuous differentials, including all projection directions in its
  JVP/VJP; concrete generated enriched-space constructors remain a later
  layer.  The neutral `assemble_enriched_feec_sequence` composition now
  accepts caller-owned shifted scalar, H(curl), H(div), and L2 maps, reports
  enriched `curl(grad)` and `div(curl)` defects, and provides complete
  product-rule value/JVP/VJP actions with independent dense matrix oracles;
  concrete map construction, cut integration, and global ownership remain
  client-owned.
  The symmetric jump-penalty block also has value/JVP/VJP actions, including
  penalty and surface-weight directions.
- Conservative/upwind vector flux interfaces now include scalar and
  metric-weighted vector value/JVP/VJP paths; entropy-stable system fluxes
  now have an explicit SPD-metric wrapper with fixed-topology value/JVP/VJP
  actions; nonlinear system entropy variables and characteristic HLL/HLLC
  laws remain application-owned.
- The public `assemble_broken_feec_sequence` contract now embeds caller-owned
  cell-local gradient, curl, and divergence maps in independent cell-major
  blocks, reports both exact-sequence compositions, and differentiates the
  maps and composition diagnostics with complete value/JVP/VJP actions. It
  preserves cell-local `curl(grad)=0` and `div(curl)=0` without introducing
  inter-cell continuity, so DG, HDG, cut/XFEM, and IGA clients can add their
  own trace laws. General mixed CG-DG coupling, higher-order local map
  construction, and global signed sparse ownership remain active work.
- The dense `assemble_glued_feec_sequence` reference path now accumulates
  signed local-to-global gradient, curl, and divergence maps and reports both
  exact-sequence compositions. Its complete local-matrix JVP/VJP actions and
  two-cell shared/reversed-DOF oracle cover the numbering contract needed by
  conforming meshes, DG/HDG, cut cells, and multipatch IGA. Sparse global
  ownership now has the duplicate-compressed `assemble_glued_feec_sequence_csc`
  counterpart with local-matrix JVP/VJP scatter and a dense-action oracle;
  distributed owner/ghost exchange and higher-order map construction remain
  separate work.
  The sparse path also exposes duplicate-compressed `curl(grad)` and
  `div(curl)` products with independent product-rule JVP/VJP diagnostics.
  The topology-only `build_multipatch_signed_dof_map` now supplies signed
  local-to-global numbering for arbitrary patch graphs and detects inconsistent
  orientation cycles. Its 2D and 3D B-spline compositions now lift the same
  contract to H1/H(curl)/H(div) face traces, including axis swaps and
  reversals. Geometry/trace-map construction, physical patch-graph assembly,
  and distributed owner/ghost exchange remain active work.

### Phase 4: Open boundaries and vector operators: **active**

- The neutral `evaluate_surface_vector_trace` contract now decomposes
  caller-owned three-dimensional surface samples into normalized normal and
  tangential traces, with fixed-topology field/normal JVP/VJP actions and an
  independent reconstruction, finite-difference, adjoint, and invalid-normal
  oracle. FEM, BEM, DtN, PML, and surface-plot clients can share this trace
  convention without selecting an equation or boundary physics.
- `compare_boundary_operator_parity` now has fixed-topology complex JVP/VJP
  actions for the common reference, all FEM/BEM/DtN/PML candidates, and
  quadrature weights. The independent multi-backend oracle covers metric
  re-evaluation, central differences, and the real-part complex adjoint while
  rejecting zero norms where the Euclidean derivative is undefined.
- The seeded `test_boundary_operator_parity_properties` fixture now generates
  finite complex reference/candidate traces and positive quadrature weights,
  checks the weighted norm value against an independent loop oracle, verifies
  JVP central differences and the real-part complex VJP identity, and confirms
  that JVP/VJP reject a zero reference norm. This extends randomized coverage
  without selecting a boundary solver or equation-specific normalization.
- Complete scalar and Maxwell FEM/BEM, DtN, and PML parity on slab, circle,
  sphere, cylinder, and torus.
- The public curvilinear PML coefficient contract now accepts a full complex
  three-by-three stretch for scalar Helmholtz and curl--curl Maxwell forms.
  It provides closed-form primal/JVP/VJP products, scale-aware singular-input
  rejection, and a diagonal reduction oracle. The neutral curvilinear
  normal-frame geometry builder now turns caller-owned layer points and unit
  normals into full per-cell stretches and differentiates centroid, normal,
  width, wave-number, and attenuation parameters with generated products. This
  is the metric and local-layer layer needed by curved meshes and IGA; closest-
  point searches, active-cell topology, and complete curved-object assembly
  remain client-owned follow-up fixtures.
- The tetrahedral first-kind Nédélec local PML form now consumes those full
  tensors. Its value/JVP/VJP products cover covariant Piola geometry,
  determinant quadrature, wave number, and all stretch entries; independent
  reassembly, complex-adjoint, and diagonal-Cartesian reduction tests guard
  the contract.
- The tetrahedral scalar P1 Helmholtz path now consumes the same full
  curvilinear scalar tensors. Its local and CSC value/JVP/VJP products compose
  analytical affine geometry derivatives with the generated coefficient map;
  an independent P1 oracle, central-difference check, complex adjoint, and CSC
  pattern test guard the contract. The constrained scalar solved-state wrapper
  now composes this CSC path with the complex implicit direct solve; independent
  re-solves and a state-level complex adjoint cover load, constraints, geometry,
  stretch, and wave-number derivatives.
- Global curvilinear tetrahedral PML CSC value/JVP/VJP assembly is now public.
  It preserves the merged pattern and scatters shared-vertex and per-cell
  stretch adjoints; the two-tetrahedron oracle checks reassembly and the
  complex adjoint. The curvilinear constrained solved-state value/JVP/VJP now
  composes that CSC path with the direct implicit solve and is independently
  checked by manufactured re-solves and a state-level complex adjoint. Curved-
  object geometry chains remain follow-up work; the neutral weighted complex
  reflection/error metrics are now public and independently checked against a
  reflected plane-wave oracle.
- The public geometry-generated tetrahedral Nédélec PML CSC wrapper now
  composes the layer-origin/normal/width builder with the global H(curl)
  assembly. Its JVP/VJP propagate mesh, layer geometry, wave number, and
  attenuation parameters through both chains, with a fixed-active-set
  reassembly and adjoint oracle. Closest-point classification and cell
  activation events remain caller-owned.
- The matching geometry-generated scalar tetrahedral Helmholtz PML CSC
  wrapper now composes the same builder with the scalar assembly, returning
  the full mesh/layer/wave-number/attenuation JVP/VJP chain under the same
  fixed-active-set event rule.
- The neutral `evaluate_surface_vector_trace` contract now decomposes
  caller-owned three-dimensional fields on FEM/BEM/DtN/PML surfaces into
  normal and tangential traces, normalizing supplied normals and exposing
  value/JVP/VJP actions. An independent central-difference and real-adjoint
  oracle covers non-unit normals and rejects degenerate surfaces; it is a
  reusable physical-field plotting and boundary-work layer, not an equation
  or solver selection.
- The public `fortfem_maxwell_curved_dtn` contract now composes weak electric
  and magnetic trace forms into a mixed RWG/RBC map with matrix-free action,
  JVP, and VJP products. `assemble_maxwell_torus_curved_dtn_rwg_3d` wires the
  map to the exact-curved torus EFIE, one-sided MFIE, and RWG mass operators;
  independent algebraic and torus assembly tests cover the fixed-topology
  contract. Caller-owned weak-to-point trace reconstruction now returns
  `B*A^{-1}` from a supplied weak mass `A` and point-evaluation basis `B`, with
  matrix-free, JVP, and VJP products; it is covered by an independent complex
  adjoint and central-difference oracle. The assembled weak map now also has
  form-level JVP/VJP products, so caller-owned EFIE, MFIE, and mass geometry
  derivatives can be composed without differentiating the map implementation.
  Geometry derivatives of those curved forms and resonance-safe regularization
  remain active.
- The curved-torus Laplace and Helmholtz representation APIs now optionally
  return the reconstructed target gradient.  Independent tests sample a
  fundamental or outgoing solution and its normal trace on two meshes, check
  field convergence and the off-surface gradient; Helmholtz BEM/DtN
  geometry-derivative comparisons remain explicit follow-up work.  The Laplace
  representation also exposes fixed-geometry data/target JVP and VJP products
  for both nodal and panel Neumann data; central differences and a real
  adjoint identity independently guard that contract.  Its fixed-topology
  geometry JVP and VJP now propagate torus parameter, major/minor radius, and
  target variations through the FortSym-generated curved-panel map; an
  independent central-difference and real adjoint oracle guards that chain.
  The outgoing Helmholtz representation now has the same fixed-topology
  geometry JVP/VJP contract for complex traces, with the real-part complex
  adjoint identity independently checked against central differences.
- `assemble_laplace_torus_curved_dtn_3d_geometry_jvp/vjp` now differentiates
  the assembled single-layer, double-layer, and trace-mass blocks, including
  the Duffy/wedge coincident-panel quadrature.  Matrix finite differences and
  a full real matrix adjoint identity independently guard the geometry chain.
- `assemble_helmholtz_torus_curved_dtn_3d_geometry_jvp/vjp` now composes the
  FortSym-backed Laplace singular products with analytical derivatives of the
  smooth outgoing-minus-Laplace correction, including coincident quadrature
  points and wave-number derivatives.  Complex matrix finite differences and
  the real-part complex adjoint identity independently guard all three blocks.
- A fast two-surface toroidal Laplace fixture moves the artificial boundary
  while keeping a common interior target.  It solves the curved BEM/DtN map on
  both surfaces, reconstructs the interior field, and checks both against an
  independent fundamental solution and against one another.  Helmholtz,
  Maxwell, PML, and FEM comparisons remain separate follow-up fixtures.
- The matching two-surface Helmholtz fixture uses an outgoing point source,
  solves the exact-curved torus BEM/DtN map on both boundaries, reconstructs a
  common interior target, and records the declared nonzero coarse-mesh
  difference.  Maxwell/PML and FEM coupling remain separate follow-up paths.
- The larger-domain Laplace and Helmholtz torus controls now have independent
  two-surface tests with fundamental/outgoing-field oracles. The Maxwell
  open-boundary gallery now solves a second, physically larger tetrahedral PML
  box, reconstructs the same interior field at shared targets, and reports the
  measured domain difference and solve time; its previous hard-coded far-wall
  placeholder is gone. The common-interior larger-domain parity now also
  exposes a fixed-topology complex JVP for field, weight, and boundary-distance
  perturbations, with an independent weighted metric and central-re-evaluation
  oracle; norm ties and zero differences are rejected as nondifferentiable.
  Exact-curved torus RWG currents now also have a public
  off-surface magnetic-field reconstruction path with a coefficient-linearity
  oracle, so the torus Maxwell gallery can show a computed field slice. That
  reconstruction now also has analytical geometry JVP and VJP products through
  straight-edge lengths, generated torus-panel geometry, radii, target, wave
  number, and coefficients; an independent central-difference and complex
  adjoint fixture covers the complete directional chain. Full assembled
  torus-block derivatives remain active follow-up work; the assembled curved-
  torus RWG mass block now has analytical geometry JVP/VJP products with an
  independent reassembly and real adjoint oracle. The assembled decaying
  torus EFIE now also exposes an exact impedance-parameter JVP/VJP product,
  independently checked by reassembly and the complex real-part adjoint
  identity, and its decaying Green-kernel decay-rate JVP/VJP now propagates
  through regular, coincident, and adaptive panel-pair quadrature with an
  independent central reassembly and complex real-part adjoint oracle;
  the propagating torus EFIE now has the matching fixed-geometry wave-number
  JVP/VJP through the same analytical quadrature chain and independent
  reassembly/adjoint checks; its propagating impedance JVP/VJP now has the
  corresponding reassembly and complex real-part adjoint oracle; geometry
  derivatives remain explicit follow-up work. The curved-sphere propagating
  and decaying EFIE blocks now expose matching fixed-geometry impedance
  JVP/VJP products, independently checked by central reassembly and complex
  real-part adjoint identities. Their wave-number and decay-rate JVP/VJP
  products now also propagate through the regular, coincident, and adaptive
  sphere panel quadrature with independent reassembly and complex adjoint
  checks. The curved-sphere RWG mass block now has analytical vertex/radius
  geometry JVP/VJP products with an independent central reassembly and real
  adjoint oracle. Exact-curved sphere RWG far-field reconstruction now has
  matching analytical geometry, radius, coefficient, direction, wave-number,
  and impedance JVP/VJP products, independently checked by central
  reassembly and the complex real-part adjoint identity. The assembled
  curved-sphere off-surface magnetic reconstruction now also has complete
  geometry, radius, coefficient, target, and wave-number JVP/VJP products;
  independent central reassembly and complex real-part adjoint tests make the
  sphere a second virtual-casing geometry rather than a visualization-only
  fixture. The assembled
  toroidal MFIE offset trace now exposes a
  fixed-geometry relative-offset and wave-number JVP/VJP product, independently
  checked by central reassembly and a complex real-part adjoint identity; its
  complete fixed-topology geometry JVP/VJP chain now also propagates torus
  vertices, angular parameters, radii, wave number, offset, refinement, BC
  traces, normals, and magnetic reconstruction, with independent reassembly
  and adjoint checks. The shared torus-normal reverse map now includes the
  cylindrical-radius contribution of the radial-difference cotangent. Full
  block derivatives remain explicit follow-up work. The exact-curved plane-wave
  RWG load now has matching geometry, direction, polarization, and wave-number
  JVP/VJP products with independent reassembly and complex real-part adjoint
  checks. Exact-curved torus RWG far-field reconstruction now has matching
  geometry, radius, coefficient, direction, wave-number, and impedance JVP/VJP
  products; independent central-difference and complex real-part adjoint
  checks cover the complete radiation-data chain. The
  `maxwell_torus_fem_bem_solution` fixture now closes the neutral vector
  FEM--BEM solved-state gap with exact edge-integral data, a complex coupled
  solve oracle, and a physical H(curl) slice with arrows. Nonzero-scattering
  vector FEM/BEM/DtN/PML parity, assembly-specific matrix/geometry derivatives,
  and larger-domain Maxwell controls remain active follow-up work. The public
  `solve_maxwell_fem_bem_linear_state` value map and its analytical JVP/VJP
  now provide the common implicit solve layer for those future products; an
  independent complex central-difference/adjoint test guards the contract.
  The actual curved-torus block also has an independent nonzero manufactured
  RWG-current state test, so the coupling gate is not limited to the
  zero-current constant-field example. The refined/localized torus RWG basis
  used by the RBC pairing now also exposes analytical geometry JVP/VJP
  products, with an independent reassembly and real adjoint oracle; this is
  the local derivative layer used by the assembled refined-to-coarse pairing.
  The torus barycentric refinement map now
  also exposes fixed-topology analytical JVP/VJP products for refined angles,
  embedded vertices, and major/minor radii, with a central-difference and
  real-adjoint oracle; coarse-to-refined geometry can therefore be composed
  without finite-differencing the refinement routine. The generic
  Buffa--Christiansen/RBC transformation now also exposes analytical geometry
  JVP/VJP products for its vertex-ring and reference-edge inverse-length
  coefficients; an independent reassembly and real adjoint-product test
  guards this BC layer. The exact-curved torus RWG/RBC pairing now composes
  those BC, localized-basis, parameter-map, normal, and parent-RWG products
  into analytical matrix JVP/VJP actions; central reassembly and a full real
  matrix adjoint oracle cover vertices, torus parameters, and both radii.
  The exact-curved sphere path now has the same composable derivative contract:
  radial barycentric refinement, localized and parent RWG basis evaluation,
  inverse parent-panel coordinates, curved normals, and the assembled RWG/RBC
  pairing expose analytical JVP/VJP products.  A central reassembly test and
  full real adjoint oracle cover sphere vertices and radius, so spherical
  geometry is no longer a finite-difference-only exception.
  The torus parent RWG basis now also differentiates moving reference
  coordinates explicitly in both value JVPs and coordinate VJPs; a direct
  geometry-plus-`(xi,eta)` central reassembly and real adjoint oracle guards
  this lower-level contract independently of the assembled RBC pairing.
  The torus Buffa--Christiansen plane-wave trace load now exposes analytical
  JVP/VJP products through barycentric geometry, the BC transformation,
  localized RWG traces, radii, polarization, direction, and wave number; an
  independent complex central reassembly and real-part adjoint oracle covers
  the assembled RHS used by the vector FEM--BEM path.
  The corresponding curved-sphere BC plane-wave load is now public with the
  same analytical JVP/VJP composition through radial refinement, the BC map,
  localized traces, radius, polarization, direction, and wave number; its
  independent central reassembly and real-part adjoint oracle closes the
  spherical vector trace-load path.
- Verify FortPlot mesh and surface rendering for every mesh-bearing plot.
- The free-boundary foundation in section 8.5 is an active cross-domain
  workstream: first close the oriented region graph and scalar exterior
  operators, then add virtual-casing and NESTOR-like Fourier maps, the
  BIEST-like vector/harmonic layer, and finally STARWALL-like response blocks
  with shape derivatives and structure-preserving time coupling. Each stage
  must land with the independent 2D/3D benchmark ladder before an external
  equilibrium or MHD adapter consumes it.

### Phase 5: Fourier-FEM and torus harmonics: **active**

- The public `fourier_mode_registry_t` records field-period phase,
  normalization, real/conjugate packing, retained triads, and caller-selected
  radial regularity. Its analytical mode value, fixed-topology JVP, and
  complex real-part VJP are independently tested. This metadata layer is
  neutral and does not import equilibrium profiles, COCOS/GEQDSK conventions,
  or plasma closures; mode-coupled operators and torus-harmonic radial
  branches remain composition work.
- The neutral `build_axis_regular_mode_table` contract validates the scalar
  axis rule (p\ge |m|, p\equiv |m|\pmod 2), reports required minimum power
  and parity, and canonically orders conjugate-safe `(m,n)` pairs independent
  of input order. `evaluate_axis_regular_radial_basis` now executes a
  caller-selected finite complex polynomial with those scalar rules and exact
  coefficient/rho JVP/VJP actions. Vector/tensor effective-mode shifts,
  alternate radial FE/IGA bases, and DESC physics remain external.
- `assemble_fourier_vector_product` now contracts a caller-owned complex
  coupling tensor over every retained triad for scalar, vector, tensor,
  H(curl), or H(div) coefficient arrays. Its JVP differentiates coupling and
  both factors, while its complex real-part VJP has an independent dot-product
  oracle. `build_fourier_mode_triad_map` now returns every ordered retained
  pair's output index and an omitted-sum count, so padding, projection, or
  rejection is an explicit caller policy rather than a hidden aliasing choice.
- `build_fourier_mode_padded_registry` now keeps the input modes first and
  appends unique ordered pair sums for one bilinear product. It preserves
  field-period/phase, real-conjugate packing, and existing radial
  powers/normalizations; newly added labels receive the neutral `abs(m)` and
  unit-normalization defaults. The constructor is deliberately one-product:
  callers must choose a larger work set for repeated nonlinear products.
- `apply_fourier_bilinear_product` now applies that product with distinct input
  and output registries, so padded/dealiased work sets are executable rather
  than metadata-only. Its analytical JVP and complex real-part VJP share the
  same explicit output-label projection and are independently checked against
  a padded contraction oracle; no model-specific bracket or coefficient
  symmetry is imposed.
- The neutral `apply_fourier_mode_linear_operator` now supplies the
  mode-diagonal block map `y(:,m)=A(:,:,m)c(:,m)` for scalar, vector, tensor,
  H(curl), and H(div) coefficient compositions. Exact matrix/field JVP and
  complex real-part VJP actions are exported through the public API and checked
  against nested-loop, finite-difference, adjoint, shape, and non-finite-input
  oracles; physical closures and basis conventions remain caller-owned.
- The neutral `evaluate_fourier_mode_expansion` now synthesizes a caller-owned
  complex coefficient vector at one `(rho,theta,phi)` point by contracting all
  retained radial/Fourier modes. It exposes analytical coordinate derivatives,
  an exact coefficient/coordinate JVP, and a real-part VJP, each checked against
  independent phase/radial, central-difference, adjoint, shape, and finite-input
  oracles; vector/tensor component composition remains caller-owned.
- `evaluate_fourier_mode_expansion_hessian` now exposes the six symmetric
  coordinate Hessian components `(rho,theta,phi)` and
  `evaluate_fourier_mode_expansion_hvp` applies that Hessian to a fixed-
  coefficient coordinate direction. Independent nested-loop algebra and
  central differences of the first gradient verify the radial, angular, and
  mixed terms; this remains a fixed-topology neutral derivative contract. A
  FortSym generation audit keeps this kernel and its oracle explicit: runtime
  integer radial powers, variable mode counts, and the `rho=0` power branches
  are not a fixed-shape CAS kernel without duplicating the metadata loop. The
  nested-loop oracle therefore remains deliberately independent of generated
  production code, avoiding correlated symbolic and implementation errors.
- The seeded `test_fourier_mode_expansion_properties` fixture now samples
  bounded radial powers, normalizations, phases, coefficients, coordinates,
  and directions across 24 deterministic cases. It compares value, gradient,
  all six Hessian components, and HVP outputs with independent nested-loop
  algebra, extending randomized Fourier coverage without changing the API or
  importing a physics closure.
- `build_fourier_mode_closure_registry` now applies the one-product constructor
  for a caller-selected positive number of rounds. Round one is the padded
  registry; subsequent rounds retain every prior-work-set pair sum. The
  bounded growth policy is checked independently and does not claim unbounded
  nonlinear closure.
- The seeded `test_fourier_mode_properties` fixture now generates bounded
  conjugate-packed mode sets, validates them, and checks every input-pair sum
  by direct label lookup. It records deterministic seeds through the shared
  property helper without importing any plasma-specific mode convention.
- `assemble_fourier_mode_energy` now provides positive caller-weighted energy
  per mode and in total, with analytical JVP/VJP actions for coefficients,
  point weights, and mode weights. `fourier_coefficients_conjugate_symmetric`
  reports the fixed-topology real-packing residual. Direct contraction,
  central-difference, and real-adjoint oracles cover both diagnostics.
- FortNum `legendre_q` and `legendre_q_derivative` are now public on `main`,
  with closed-form Q0--Q3, centered-difference derivative, and invalid-domain
  tests. Its toroidal P/Q branch is likewise covered by independent value and
  ODE checks. Standard orthonormal complex spherical harmonics and angular
  derivatives are now public with closed-form degree-one, conjugacy, pole,
  and invalid-angle tests. Its Gaunt product coefficient is now public and
  independently checked against tensor-product angular quadrature, with
  exact triangle/parity/azimuthal selection rules. Complete the remaining stable high-order and
  near-cut continuation for moderate orders is now covered by a DLMF
  zero-balanced/Euler continuation. Stable degree-80 half-integer continuation
  is now on FortNum `main`: P uses generated upward recurrence and Q uses
  scaled Miller backward recurrence away from the cut, with independent
  50-digit value and three-term-recurrence checks. Uniform asymptotics for
  arbitrarily large degree/order remain before calling this a complete
  production DtN/torus-harmonic envelope.
- The FortFEM public API now re-exports the pinned toroidal P/Q values and
  derivatives through `fortfem_toroidal_harmonics`; an independent adapter
  oracle covers both Hobson branches, derivative values, invalid-domain NaN,
  and the upstream near-cut continuation. Stable degree-80 continuation is
  pinned to the FortNum revision recorded in `fpm.toml`; mode-coupled toroidal operators and uniform
  asymptotic envelopes remain.
- The toroidal P/Q adapter also exposes second radial derivatives from the
  pinned FortNum branch, with independent associated-Legendre ODE residual
  checks for both branches. This closes the local second-order radial hook;
  uniform asymptotics and mode-coupled operators remain separate work.
- The neutral `evaluate_toroidal_spectral_trace` map now composes the P/Q
  branches with complex toroidal Fourier phases and outward normal traces;
  its analytical coefficient/geometry JVP/VJP products are independently
  checked against the existing harmonic/DtN mode and central/adjoint oracles.
- The seeded `test_toroidal_spectral_trace_properties` fixture now samples
  valid degree-at-least-one, order-bounded P/Q mode sets, finite complex
  coefficients, and toroidal angles with the check RNG. It compares scalar
  and grid values with an independent separated modal sum, checks grid and
  pointwise JVPs by central re-evaluation, and exercises negative-index and
  nonpositive-scale/eta rejection without selecting a toroidal convention.
- The seeded `test_toroidal_spectral_neumann_properties` fixture now samples
  the same FortNum-valid nonzero P/Q modes and finite normal data, compares
  modal division with an independent outward-normal factor, checks central
  JVP and real-complex VJP identities, and verifies explicit zero-mode and
  resonance-tolerance rejection.
- The independent `test_toroidal_harmonic_parity` fixture now checks both
  P/Q branches against the separated trace normalization, angular-period and
  conjugate-reflection parity, the outward-normal sign, linear modal
  contraction, and the associated-Legendre radial ODE.  It is an analytical
  fixed-topology contract only; equilibrium profiles and plasma readers
  remain external.
- Define mode normalization, phase, field-period, and real packing.
- Add mode-coupled scalar, H(curl), H(div), and caller-defined nonlinear
  operators.

### Phase 6: Structure-preserving time evolution: **active**

- The neutral resistive-wall RL block now advances by implicit midpoint and
  reports a discrete magnetic-energy/input/dissipation ledger; its analytical
  JVP/VJP actions are independently checked. Coupling this block to a
  caller-owned FEM/BEM or DtN field port is now covered by the neutral
  `advance_mixed_wave_wall_midpoint` composition. The coupled block assembles
  the skew field/wall port terms and dissipative RL term in one midpoint solve,
  reports an independent wave/wall energy ledger, and exposes complete
  matrix/port/state/time JVP/VJP products; geometry and wall normalization
  remain caller-owned.
- The public `advance_mixed_wave_midpoint` step now provides the common
  first-order pressure/velocity, displacement/momentum, and port-Hamiltonian
  Cayley contract. Its analytical JVP/VJP now differentiate the complete
  block map and time-step product; an independent nonidentity-block test
  checks central differences and the real adjoint identity in addition to the
  oscillator map, energy, and signed-step reversibility.
- The public `advance_mixed_wave_symplectic_euler` step now provides a
  partitioned first-order symplectic update for the same mixed state. Its
  independent test checks the two-stage mass-solve oracle and the canonical
  two-state symplectic-form identity; its analytical JVP/VJP now cover the
  complete two-stage map and time-step product, with an independent central
  difference and real adjoint oracle. Dissipative terms remain separate.
- The neutral `assemble_symplectic_map_defect` contract now checks an
  arbitrary caller-owned linear step against an arbitrary state two-form and
  differentiates the map and form. Its independent matrix, central-difference,
  and real-adjoint oracle is reusable by mixed acoustics, vibration,
  elasticity, and electromagnetic clients without selecting their physics.
- The public `advance_dissipative_cayley` step now provides the separate
  `(M+hD/2)^(-1)(M-hD/2)` update for positive-time dissipative blocks. Its
  JVP/VJP differentiate mass, damping, step size, and state, while an
  independent two-by-two oracle checks non-increasing SPD-mass energy. It is
  intentionally not labeled symplectic; splitting and absorbing clients own
  the composition.
- Variational/symplectic and Poisson building blocks for ideal terms.
- Energy-dissipative integrators for resistive and viscous terms.
- The public `advance_mixed_wave_strang` step now composes caller-owned ideal
  mixed-wave Cayley factors as A(Δt/2)-B(Δt)-A(Δt/2). Its analytical JVP/VJP
  propagate state and signed-step products through all three factors, and an
  independent nonidentity-block test checks energy preservation, reversibility,
  central differences, and the real adjoint identity. A separate 2000-step
  campaign now checks the quadratic Hamiltonian against an independent
  mass-matrix oracle and recovers the initial state under signed reversal.
  A seeded 16-case property campaign now repeats the energy and signed-step
  checks for diagonal SPD masses and random coupling blocks, covering both
  midpoint and Strang without global random state.
  The seeded `test_mixed_wave_3d_properties` campaign now adds 16 fixed-size
  three-component cases with dense Cholesky-built SPD mass blocks and
  diagonally dominant port couplings. An independent block-matrix oracle
  checks midpoint values, energy, signed-step reversal, and the mass-weighted
  symplectic form; it also checks the partitioned symplectic-Euler map and
  separates its symplectic defect from its non-energy-preserving behavior.
  The neutral `advance_quadratic_avf` primitive now supplies the
  caller-owned quadratic Hamiltonian average-vector-field/discrete-gradient
  update, with analytical JVP/VJP products in state, Hamiltonian, skew
  interconnection, linear term, and step size.  An independent three-mode
  oracle checks the linear solve, Hamiltonian preservation, signed-step
  reversibility, central-difference JVP, and real-adjoint VJP; dissipative
  terms remain separate.
  The independent `test_mixed_wave_3d_structure_oracle` fixture extends this
  contract to explicit Cartesian (x/y/z) components: an analytical
  three-component solution, energy ledger, canonical six-dimensional
  symplectic-form defect, signed-step reversal, and a positive-time
  dissipative-Cayley energy decrease with a nonzero symplectic defect.  The
  dissipative block is therefore tested as a separate non-symplectic map.
  Dissipative splitting remains separate; discrete-gradient/average-vector-
  field options are now covered by this neutral quadratic slice, while
  broader problem-size campaigns remain active work.
- Mixed first-order pressure-velocity, displacement-velocity-stress, and
  electromagnetic wave states with a common port-Hamiltonian interface.
- Tensor-valued coefficients and anisotropic constitutive blocks with exact power,
  momentum, and stress-work diagnostics.
- The generated CGL-shaped tensor and product-rule force-divergence blocks are
  public and tested as generic constitutive fixtures. The tensor now also
  feeds compositional traction and `P:grad(v)` stress-work value/JVP/VJP
  diagnostics. `assemble_tensor_volume_work` now provides the neutral full
  volume `stress:grad(test)` residual with tensor, geometry, and quadrature
  JVP/VJP actions and an independent contraction oracle; constitutive laws and
  caller-supplied correction tensors remain external. No plasma closure is
  implemented.
- The neutral `evaluate_field_aligned_constitutive_tensor` contract now
  composes the CGL parallel/perpendicular projectors with an optional skew
  Hall/gyrotropic cross-product block. FortSym generates the three independent
  Hall-direction products and their JVP/VJP without duplicating the generated
  symmetric CGL projector; the wrapper owns only their fixed skew layout.
  Independent scalar-product, central-difference, nonsymmetric transpose,
  optional-zero, and invalid-input tests cover both limits. This remains a
  pointwise constitutive ingredient: geometry pullbacks, quadrature, and
  physical closures stay caller-owned.
- The seeded `test_field_aligned_constitutive_properties` fixture now samples
  bounded unit directions and anisotropic/Hall coefficients with the check
  RNG, compares an independent projector-plus-Hall oracle, checks symmetric,
  skew, and total power identities, and exercises near-unit acceptance and
  non-unit rejection bounds. This adds randomized constitutive coverage while
  keeping the closure-neutral tensor contract and the short feedback path.
- `assemble_tensor_diffusion_matrix` now provides the neutral quadrature
  contraction `grad(test)^T K grad(trial)` for strongly anisotropic scalar,
  elastic, resistive, and compatible-field clients, with full tensor,
  geometry, weight, JVP, and real VJP actions. A hand-evaluated anisotropic
  matrix, central-difference, adjoint, and positivity oracle cover the public
  block; global numbering and constitutive laws remain client-owned. The
  compatibility-preserving 2D entry point now shares a strict arbitrary-
  dimension implementation with new 3D and `nd` wrappers, including matching
  JVP/VJP actions and independent 3D contraction/dimension-validation tests.
- The neutral mixed elasticity residual now exposes caller-owned compliance,
  strain, and equilibrium maps as explicit constitutive and force-balance
  blocks with complete value/JVP/VJP actions. Its independent matrix,
  finite-difference, and adjoint oracles make it the reusable algebraic core
  for Hellinger--Reissner, weak-symmetry, TDNNS, and anisotropic clients; the
  elasticity complex, symmetry constraints, and boundary laws remain external.
- The matching neutral weak-symmetry constraint now exposes `W sigma - g`
  with independent matrix, finite-difference, and adjoint oracles. This is
  the reusable multiplier/trace block for Arnold--Falk--Winther and TDNNS
  compositions; stress-space construction and physical symmetry laws remain
  client-owned.
- The neutral `assemble_coupled_field_residual` contract now composes a
  rectangular field operator and a separate rectangular constraint operator,
  with independent matrix, finite-difference, and real-adjoint oracles. It is
  the reusable residual boundary for multi-field FEM/BEM/DtN/PML, tensor, and
  interface clients; global sparse ownership and application closures remain
  external.

### Phase 7: Equilibrium and linear-response ingredients: **active**

- The neutral `equilibrium_interchange_t` record now defines the external
  adapter boundary for mapped/physical samples, named coefficient and profile
  arrays, explicit scalar/vector/tensor component ranks and offsets, segmented
  boundaries, provenance, units, and normalization. It is validated with
  analytical toroidal and 2D manufactured fixtures plus deep-copy and
  rejection checks. FortFEM still does not implement GEQDSK, COCOS, or any
  plasma-specific reader, profile law, coil model, or equilibrium solver.
  The physics-free `build_equilibrium_interchange_sample_set` adapter now
  selects unique coefficient components from those physical samples, retains
  caller-owned weights and provenance, and feeds the common comparison schema;
  interpolation and external-code resampling remain outside FortFEM.
- The neutral `linear_response_interchange_t` record now defines the external
  adapter boundary for modal `(m,n)` metadata, complex frequency, provenance,
  normalization, response channels, and caller-owned equilibrium, inertia,
  resistive, vacuum, and wall blocks. The generated algebraic operator
  (K-\omega^2M+i\omega R+V+W) and its caller-owned forced residual have
  independent complex JVP/VJP products,
  normalized reciprocity error, and a conservative Hermitian passivity lower
  bound with manufactured oracles. GPEC, MARS-F, GLISS, STARWALL readers,
  singular-layer physics, and plasma normalization remain external.
- The neutral linear perturbation composition now separates caller-owned
  Lorentz, pressure/stress, inertia, vacuum, wall, resistive, and singular-layer
  matrices in `L + P + V + W + S - omega^2 M + i omega R`. Exact complex
  JVP/VJP actions cover every block and the frequency, while validated metadata
  records provenance, normalization, and declared reciprocity/passivity
  conventions. Plasma closures, singular asymptotics, and readers remain
  external. The manufactured `linear_perturbation_response` example composes
  the seven blocks with the public forced residual, reconstructs the complex
  Fourier state over physical poloidal angle before diagnostics, and records a
  frequency sweep of reciprocity, passivity, frequency-JVP accuracy, residual,
  and small runner-local timings.
- The neutral generalized complex eigen-residual `K u - lambda M u` now has
  independent matrix-product, finite-difference tangent, and real-complex
  adjoint oracles. It is a solver-independent contract for curl--curl,
  Fourier-FEM, GLISS-like, and other modal clients; eigensolver choice,
  normalization, equilibrium import, and application physics remain external.
- The neutral `interchange_sample_set_t` contract now validates common
  physical coordinates, scalar/vector/tensor value arrays, positive weights,
  and provenance. Its weighted absolute, relative, and componentwise
  L-infinity comparison rejects mismatched point sets, giving license-safe
  CHEASE/FreeGS/GPEC/MARS/GLISS samplers one shared numerical oracle without
  adding any application reader. Fixed-coordinate real weighted L2/relative
  comparison JVP/VJP actions are now public and independently checked; the
  complex counterpart uses the real-part adjoint convention for frequency-
  domain response samples.
- The neutral singular-layer matching block now composes independently sized
  inner and outer complex trace spaces as a weighted oriented jump residual,
  with complete value/JVP/VJP actions and an independent central-difference,
  real-adjoint, and invalid-weight oracle. Asymptotic layer models,
  normalization, and physical jump laws remain external.
- Generic axisymmetric elliptic and Fourier variational fixtures that can be
  sampled by external CHEASE or FreeGS adapters. FortFEM will not implement
  GEQDSK or COCOS readers, profile laws, coil models, or equilibrium solvers.
- Generic multi-region curl-eigenproblem, interface, and constraint blocks that
  external VMEC, GVEC, DESC, or SPEC clients can exercise.
- DESC-specific foundation is tracked separately from the VMEC/GVEC variational
  path: axis-regular Fourier--Zernike-compatible nested coordinates, direct
  radial/helical force residuals, objective/constraint registries, exact
  converged-state derivatives, perturbation/continuation/deflation, neutral
  flux-surface/Boozer/near-axis diagnostics, and fixed/free-boundary sheet and
  vacuum ports. No DESC reader, profile law, or optimization objective is
  implemented here.
- The neutral `evaluate_force_balance_objective` contract now provides the
  positive-weighted quadratic objective and exact residual/weight JVP/VJP
  needed by direct-force optimization clients; it deliberately remains
  independent of profiles, coordinates, optimizers, and plasma closures.
- The metadata-driven `evaluate_equation_objective_merit` value/JVP/VJP
  contract now applies caller-owned targets, positive weights, scales, and
  activation to equation/objective/constraint rows, giving DESC/GVEC-like
  clients a reusable normalized merit without selecting an optimizer or
  physical objective.
- Generic linear ideal and resistive response blocks, singular layers, vacuum,
  conducting-wall traces, and response matrices for external GPEC, MARS, and
  GLISS oracle data.
- Compose a small manufactured multi-field state from independently verified
  field, tensor, surface-current, vacuum, wall, and constraint blocks. Plasma
  state selection and closure remain outside FortFEM.

### Phase 7a: Field-aligned edge and SOL ingredients: **active**

- FCI RK4 field-line tracing, support-operator gradient/divergence algebra, and
  matrix-free (P K_\parallel Q) diffusion (including FortSym-generated
  gradient and full diffusion JVP/VJP products) are public and independently
  tested. The
  `fci_parallel_diffusion` gallery fixture provides a reproducible open-line
  mass-conservation and negative-energy profile. A 1D linear map builder and
  its fixed-topology JVP/VJP provide independent partition/affine and
  dot-product oracles, and a 2D Cartesian bilinear builder now covers a
  genuine poloidal slice. A generated quadratic Lagrange map covers explicit
  three-node nonuniform stencils on 1D slices, including fixed-stencil JVP/VJP
  products. The RK4 tracer also has a tangent callback path with an exponential
  endpoint oracle. The positive staggered
  flux-box volume constructor
  now covers the traced expansion/area/(B_\varphi) product with pinned
  FortSym-generated value/JVP/VJP kernels. Fixed-cell 2D map JVP/VJP products
  now cover source and target coordinate motion; generated quadratic and
  cubic and quartic fixed-stencil maps now expose value/JVP/VJP actions with
  independent polynomial, finite-difference, and real-adjoint oracles. A
  generated quintic fixed-stencil map now exposes the same value/JVP/VJP
  actions with independent quintic, finite-difference, and real-adjoint
  oracles. A generated sextic fixed-stencil map now exposes the same
  FortSym-generated value/JVP/VJP actions with independent sextic,
  finite-difference, and real-adjoint oracles. A generated fixed-topology
  quadrilateral area map now supplies
  positive unstructured plane-cell measures with independent shoelace,
  finite-difference, and real-adjoint oracles plus a gallery fixture. Higher
  interpolation derivatives beyond sextic and curved support-volume measures
  beyond sextic remain separate planned components. The generic polygon map
  covers fixed-topology cells with more than four vertices, and its generated
  quadratic Bezier-edge extension covers arbitrary curved polygon boundaries.
  Generated quadratic, cubic, quartic, quintic, and sextic Bezier-edge area maps now
  supply fixed-topology curved-cell value/JVP/VJP contracts with independent
  Gauss--Green oracles; quadratic through sextic are also sampled in the
  boundary gallery.
- The batched 2D bilinear endpoint-to-map adapter now connects traced
  forward/backward endpoints to the support-operator tensor contract and
  carries fixed-topology source-grid and endpoint JVP/VJP actions. Moving
  connectivity, higher-order curved measures beyond sextic, and stencil
  rebuilds at topology events remain planned; the generic straight polygon map and
  quadratic Bezier-edge quadrilateral map are now the unstructured-cell
  baselines. The deterministic composed map fixture also guards this adapter
  against non-reproducible endpoint ordering or weighted-adjoint regressions.
- The fixed-stencil quadratic 1D map now has a batched segment adapter with
  generated-kernel-backed JVP/VJP accumulation and independent polynomial,
  finite-difference, and adjoint tests. Degenerate local source coordinates
  are rejected before evaluating the Lagrange products. Curved or moving-cell
  connectivity remains a topology-rebuild concern.
- A positive FCI diffusion diagonal is public as the first anisotropy-aware
  preconditioner contract, with a matching Jacobi apply and positivity test;
  PARALLAX-compatible plane multigrid and additive field splitting are now
  public baselines, while coupled stronger blocks remain planned.
- A PARALLAX-style anisotropic split action now combines independent per-plane
  CSC elliptic blocks with the conservative FCI parallel operator; tensor
  coefficient assembly and coupled stronger field splitting remain planned.
  `compute_fci_anisotropic_diffusion_diagonal` now combines the
  positive support diagonal with each plane block's oriented diagonal and
  rejects non-positive results, providing an independently tested scalar
  preconditioning oracle. `apply_fci_anisotropic_jacobi_preconditioner` applies
  that diagonal directly for small matrix-free solves; cached diagonal use and
  `apply_fci_plane_two_level_vcycle` now provides the next plane-solver layer;
  a public recursive multilevel V(1,1) path now covers arbitrary level sizes;
  the retained-factor path and additive field-split composition are tested
  separately. A public recursive W(1,1) variant now repeats coarse visits;
  coupled blocks remain active.
- The split action now also has a field-only VJP that composes the conservative
  FCI transpose with an explicit transpose of every plane CSC block. An
  independent nonsymmetric-plane oracle and real dot-product test guard this
  adjoint contract; coefficient and geometry derivatives remain separate.
- A batched `apply_fci_plane_two_level_vcycles` adapter now applies the
  independently tested two-level plane cycle to a homogeneous plane stack,
  preserving the PARALLAX field-split boundary between per-plane elliptic
  solves and the global FCI line action. A ragged-offset companion now covers
  variable plane sizes without padding. The public additive field-split
  preconditioner combines cached parallel Jacobi and ragged plane cycles with
  explicit nonnegative weights. A public retained coarse-factor path now
  reuses FortSparse factorizations across right-hand sides, and the recursive
  multilevel V-cycle accepts flat level offsets for nonuniform hierarchies, and
  the W-cycle variant repeats coarse corrections; coupled blocks remain active
  work.
- A fixed-cell barycentric triangle interpolation path now covers logically
  unstructured poloidal targets, including geometry and target JVP/VJP actions;
  its batched endpoint-to-map path now feeds the support-operator tensor
  contract; moving-cell connectivity and higher-order stencils remain planned.
- FCI field-line maps, higher-dimensional interpolation Jacobians, and parallel
  derivative JVPs.
- A geometry-only 2D FCI first-hit search now returns the nearest transverse
  oriented wall/target segment, exact hit parameter/point, and facet normal;
  valid no-hit traces and malformed facets have explicit status contracts. Its
  fixed-topology JVP differentiates hit point, parameter, and facet normal
  against central differences. `assemble_fci_terminal_boundary_flux` now
  provides the volume-weighted conservative contribution with fixed-owner JVP
  and VJP dot-product oracles; owner remaps and physical material laws remain
  separate contracts. The matching 3D triangle search and fixed-topology JVP
  now cover triangulated toroidal vessel/divertor surfaces without copying a
  PARALLAX implementation. A traced-polyline wrapper now scans those facets in
  field-line order and returns the first segment/triangle, endpoint, oriented
  normal, connection length, and fixed-event JVP; no-hit paths return their
  final endpoint and total length.
- The generic nonlinear material-surface flux contract is now public: an
  application callback supplies the local tagged law while FortFEM assembles
  the oriented trace residual, integrated per-tag ledger, and fixed-tag
  value/JVP/VJP products with independent central-difference and dot-product
  oracles. Sheath, Bohm, recycling, radiation, and material databases remain
  application-layer clients.
- Strongly anisotropic diffusion, conduction, resistivity, and viscosity.
- Treat material boundaries as explicit geometry rather than only diffuse
  penalization. [#58](https://github.com/lazy-fortran/fortfem/issues/58)
  owns exact first-hit FCI terminal segments on oriented wall/target facets;
  [#59](https://github.com/lazy-fortran/fortfem/issues/59) owns the matching
  conservative terminal boundary-flux contribution and fixed-topology
  derivatives.
- Keep physical boundary laws outside the numerical library while making them
  first-class residual terms. [#60](https://github.com/lazy-fortran/fortfem/issues/60)
  owns the generic nonlinear material-surface flux value/JVP/VJP contract;
  Bohm, sheath, recycling, sputtering, radiation, and heat-transmission models
  are application-layer clients.
- Support an optional fitted boundary-layer patch without forcing the bulk FCI
  mesh to conform to the vessel. [#61](https://github.com/lazy-fortran/fortfem/issues/61)
  now has the neutral `assemble_fci_boundary_patch_mortar` cross-mass,
  constant-preserving transfer, weighted-adjoint, overlap-measure, and
  ownership-multiplicity contract with fixed-topology JVP/VJP products and
  independent matching, reversed, duplicate, zero-measure, rank-deficiency,
  finite-difference, and dot-product tests. Geometry construction, general
  moving Chimera meshes, and application boundary laws remain outside FortFEM.
- GORILLA owns reusable characteristic stepping and material events, not
  FortFEM: [GORILLA #80](https://github.com/itpplasma/GORILLA/issues/80)
  introduces a common characteristic event contract and
  [GORILLA #81](https://github.com/itpplasma/GORILLA/issues/81) adds neutral
  free flight without misusing guiding-centre equations, while
  [GORILLA #82](https://github.com/itpplasma/GORILLA/issues/82) covers only
  those fluid-advection characteristics that genuinely admit marker tracing.
  Marker collisions, weights, and conservative cell/surface tallies belong
  above the pusher in
  [GORILLA_APPLETS #47](https://github.com/itpplasma/GORILLA_APPLETS/issues/47).
- Small MMS and performance cases aligned with PARALLAX, GRILLIX, and BOUT++
  concepts, without copying their implementations. The FCI gallery fixture
  now records the measured matrix-free action cost alongside its conservation
  and dissipation diagnostics; larger problem-size scaling remains active.
- Document the generic equation-as-data callback ABI for caller-owned fields,
  coefficients, sources, boundary laws, and target ledgers. Keep species,
  closures, sheath, and material physics in client code while testing units,
  signs, residuals, and balances through the public callback contract.

### Phase 8: Cross-code oracles and gallery: **active**

- Ordered examples with FortPlot plots, numerical data, convergence, and
  performance.
- The solver benchmark now publishes row-oriented ILUT/ICHOL construction
  timings alongside the physical Poisson solution, direct/PCG timings, and
  residual diagnostics.
- Every example carries the smallest applicable complete contract: a
  manufactured or analytical solution, an independent oracle, primal and
  derivative checks, a conservation or structure diagnostic, a mesh or
  geometry view, and a timing and memory record. Examples progress from a
  one-dimensional Poisson patch test through two-dimensional interfaces and
  waves to three-dimensional torus, sphere, FEM/BEM/DtN, PML, anisotropic,
  and mixed structure-preserving cases.
- The `mixed_acoustic_wave` gallery fixture now exercises the public
  structure-preserving midpoint step against independent two-mode acoustic
  solutions, with physical state, phase-space, energy, error, and timing
  plots. The same neutral block contract remains reusable for elasticity,
  electromagnetic, and tensor-pressure clients.
- The `mixed_wave_3d_structure` gallery fixture now adds the smallest
  explicit Cartesian (x/y/z) physical trajectory: its first FortPlot is a
  connected 3D numerical/analytical solution path, followed by component and
  energy diagnostics. A bounded shell oracle checks nonempty PNG/CSV output
  and `physical_solution` before `diagnostics`; no production physics is
  selected.
- The `mixed_elasticity_wave` gallery fixture now exercises the mixed
  displacement/velocity and constitutive/equilibrium contracts on a
  manufactured one-dimensional elastic bar. Its first plot reconstructs the
  physical displacement and stress fields; analytical modal, reversibility,
  energy, residual, and timing oracles remain secondary diagnostics. This is
  a foundation fixture, not a production elasticity solver.
- Gallery examples use the same generated kernels and public APIs as library
  clients. The gallery does not contain a faster special implementation that
  bypasses the tested residual or derivative path.
- Optional lightweight FEniCSx, FreeFEM, and MFEM runners.
- Sister-repository data for heavy or licensed references.
- GitHub Pages generation, link checks, and periodic deployment monitoring.
- The test workflow now continues through example generation after a failed
  verification stage, uploads its plot artifact, and reports the failure only
  after publication. The Pages workflow consumes the triggering test-run
  artifact even when the test conclusion is red, so known oracle failures do
  not hide the ordered gallery or its plots.
- Seeded property tests and a common sampler compare independent codes without
  assuming matching bases, numbering, or mesh topology.
- The gallery now requires every generated example page to embed a real PNG
  plot; the optional interoperability page retains an explicit SVG provenance
  cover until its external FEniCSx/FreeFEM/MFEM artifact is available. The
  formerly empty circular DtN, circle BEM/CFIE, transmission, mesh, solver,
  and tetrahedral Nédélec pages now emit numerical convergence, response,
  mesh, or performance plots. The test workflow uses the valid
  `fpm run --example --list` compile gate rather than the unsupported
  `fpm build --example` spelling.
- FortPlot mesh-bearing examples have a regression fixture that checks element
  count, boundary edges, patch boundaries, and internal surfaces in the
  rendered output before Pages deployment.
- The paper-magnetic three-dimensional oracle truncates its independently
  converged odd membrane series at mode 399, keeping the reference and curl
  checks while bringing that fast-suite target below the short feedback
  budget. Higher truncation studies remain an explicit benchmark concern.
- `doc/examples/primary_plots.txt` is the gallery's explicit visual contract:
  each example names its physical solution, field, geometry, or mesh artifact.
  The documentation generator copies that artifact to `primary.png` and emits
  it before convergence, timing, conditioning, or other diagnostic plots. The
  generator now fails when a mapped artifact is absent instead of silently
  substituting a generated cover, and the docs test rejects diagnostic names
  such as convergence, error, timing, residual, accuracy, or DOF plots as the
  first view. FortPlot's 3-D scatter path must preserve per-sample colour maps
  so solution samples remain physically interpretable. CI collects nested
  per-example output with an independent media validator and keeps the short
  ten-second limit for focused tests while allowing the more expensive
  physical 3-D gallery products a bounded execution window.
- The executable examples follow the same physical-first order at runtime as
  the Pages gallery.  The circular DtN fixture now emits its filled disk field
  before the boundary-trace and response curves, so running an example
  directly presents the solution before diagnostics as well.
- The Maxwell open-boundary example now uses a filled mid-plane magnitude slice
  of its solved Nédélec/PML edge field as the primary preview; the TE/TM DtN
  spectrum remains a secondary operator diagnostic rather than the gallery
  cover.
- Boundary examples now lead with physical fields rather than abstract spectra:
  circular DtN uses a filled disk mode, circle Laplace/Helmholtz BEM and CFIE
  use incident or harmonic disk fields, symmetric transmission uses a filled
  interior solution, and curved-torus Maxwell scattering uses its 2-D RCS
  surface map. Their trace, spectrum, convergence, and geometry plots remain
  available as secondary diagnostics.
- The public vector-plot API now reconstructs lowest-order Nédélec and
  Raviart--Thomas fields from their edge degrees of freedom, and interpolates
  nodal vector fields barycentrically; it no longer fabricates a synthetic
  vector field for a gallery preview. An independent constant-field oracle
  covers this path. FortPlot `main` was checked at the pinned revision for
  mesh separators, dynamic storage, and per-sample 3-D colours, so no
  downstream copy of an upstream FortPlot fix is maintained here.
- The gallery's physical previews now sample the computed fields rather than
  plotting coefficient vectors or manufactured vertex labels: JOREK's
  B-spline flux is reconstructed on a regular cell-centre grid, the Maxwell
  PML preview reconstructs the complex order-1 Nédélec solution in a physical
  slice, and the tetrahedral H1, mixed RT--DG, and Nédélec previews sample
  their solved finite-element fields in the cells. The torus geometry preview
  draws periodic parameter lines instead of projected triangle diagonals, so
  the curved surface is readable without concealing the coarse solver mesh.
  These are visualization contracts and do not replace convergence,
  conservation, or independent manufactured-solution tests.
- The Maxwell open-boundary gallery additionally plots the analytical plane
  wave and reconstructed fields from both the reference and larger PML boxes;
  an independent test guards the shared-target field and domain-difference
  values before publication.
- The curved-torus Maxwell gallery now leads with a reconstructed scattered
  magnetic-field slice; the RCS surface, trace, weak-DtN response, and exact
  geometry remain secondary diagnostics.
- The toroidal Maxwell FEM--BEM manufactured fixture now leads with its
  reconstructed solved field magnitude, followed by visible H(curl) arrows,
  a log-scale round-off error map, a one-dimensional cross-section oracle, and
  a curved torus geometry view. Its independent test solves the same assembled
  complex block system from a manufactured right-hand side, so the gallery
  image cannot be mistaken for a coefficient or source-field preview.
- The `xfem_interface_solution` fixture now gives the shifted scalar and
  componentwise vector enrichment a physical-first two-dimensional preview:
  the scalar solution is primary, the vector magnitude carries explicit
  normalized arrows, and the signed enrichment contribution is secondary.
  Independent scalar and vector sign reconstructions remain executable
  oracles; this is a space-construction foundation, not an interface PDE
  solver. Cut-cell integration, Piola-compatible FEEC numbering, and sparse
  assembly stay caller-owned as specified above.
- The regularized surface-current gallery now leads with the resolved
  three-component Gaussian sheet profile and independent oracle markers;
  integrated-current convergence and a bounded per-profile timing record are
  secondary. The linear perturbation gallery similarly leads with the
  physical complex Fourier response before reciprocity, passivity, and
  residual diagnostics. Both pages retain ignored generated media and
  machine-readable CSV/JSON provenance.
- The `sheet_current_parity` gallery adds a physical-first slab comparison of
  an explicit surface ledger and a resolved Gaussian layer. Its independent
  profile/quadrature oracle and generated CSV/JSON parity record are reusable
  for later cylinder, sphere, and torus surface parametrizations.
- The nested-surface gallery now leads with a colored, parameter-lined 3-D
  torus solution and follows it with parameter-space and axis-regular radial
  diagnostics. The exterior sphere BEM gallery now leads with the computed
  Laplace field on a 3-D observation shell, using the parametric-surface
  renderer plus per-sample colors and an independent (1/r) oracle; density
  and solver diagnostics remain secondary.

### Phase 9: Future application layer: **reference only**

Use the stabilized ingredients to support future research applications that
may resemble GVEC, VMEC, SPEC, GPEC, MARS, GLISS, or JOREK. Do not start a
production code in this repository until the corresponding ingredient,
analytical example, and external oracle have passed the previous phases.
Steady edge/SOL/divertor transport, Braginskii closures, species and reaction
networks, neutral backends, impurity charge states, sheath and material laws,
GORILLA/EIRENE coupling, and target heat-load accounting belong in a separate
application repository. FortFEM owns only their reusable geometry, trace,
operator, residual/JVP/VJP, field-solver, and preconditioner foundations.

## 15. Definition of done

An ingredient is complete only when all applicable checks below pass:

- public API and units are documented;
- the residual and active constraints are explicit;
- FortSym-generated or independently derived kernel is provenance mapped;
- primal, JVP, VJP, and implicit solve derivatives have dot-product or
  finite-difference checks;
- geometry, trace orientation, and topology events are tested;
- an independent analytical or external oracle agrees within a declared
  tolerance;
- convergence, conditioning, conservation, and performance are reported;
- performance is compared with a named reference on a fixed workload, with
  hardware, compiler, thread count, memory, accuracy, and derivative cost
  recorded;
- structure-preserving properties are tested rather than asserted;
- seeded random properties pass with their seed, generated case, and shrink
  record retained in the test log;
- external-code comparisons carry a license, executable version, source
  revision, sampler, tolerance, and data checksum;
- focused `fo` tests meet the short feedback target;
- full CI passes on supported compilers;
- the example and documentation are generated reproducibly;
- transient output is ignored and no unrelated downstream fix is duplicated.

## 16. Non-goals and explicit boundaries

- No full VMEC, GVEC, SPEC, GPEC, MARS, GLISS, JOREK, CHEASE, or FreeGS
  replacement is planned in FortFEM.
- No GEQDSK or COCOS parser, equilibrium profile library, coil or vessel
  physics, plasma closure, species or reaction model, sheath, neutral,
  impurity, or divertor application is planned in FortFEM. These are external
  clients of the generic contracts.
- No C rewrite of numerical kernels is planned.
- No proprietary source, binary, or license-restricted benchmark data is
  checked into this repository.
- No claim of differentiability is made across a changing interface topology,
  adaptive mesh decision, enrichment activation, or eigenvalue crossing
  without a dedicated event or subgradient contract.
- No generic time integrator is accepted as structure preserving without an
  invariant or geometric test.
- No spacetime FEM implementation is in the current scope.

## 17. Literature and provenance

The following references motivate the architecture. Links are kept to papers,
official documentation, or official repositories where possible.

### Interfaces, enrichment, and IGA

- [Arnold, Falk, and Winther, finite element exterior calculus](https://arxiv.org/abs/0906.4325)
- [Fries, corrected XFEM blending elements](https://doi.org/10.1002/nme.2259)
- [Stable isogeometric analysis of trimmed geometries](https://www.sciencedirect.com/science/article/abs/pii/S0045782516308222)
- [XIGA for multi-material problems](https://doi.org/10.1007/s00466-022-02200-y)
- [Independent-field isogeometric boundary elements](https://arxiv.org/abs/1406.0306)
- [DefElement reference](https://defelement.org/)
- [MFEM features and finite-element reference implementation](https://mfem.org/features/)

### Special functions and toroidal harmonics

- [DLMF 14.3: definitions and hypergeometric forms](https://dlmf.nist.gov/14.3)
- [DLMF 14.10: Legendre recurrences and derivatives](https://dlmf.nist.gov/14.10)
- [DLMF 14.19: toroidal (half-integer) specialization](https://dlmf.nist.gov/14.19)
- [FortNum special-function API at the pinned revision](https://github.com/lazy-fortran/fortnum/blob/d8be030/docs/api.md)

### MHD, equilibria, and linear response

- [CHEASE official page](https://crppwww.epfl.ch/~sauter/chease/)
- [CHEASE toroidal-equilibrium paper](https://doi.org/10.1016/0010-4655(96)00046-X)
- [FreeGS documentation](https://freegs.readthedocs.io/en/stable/creating_equilibria.html)
- [PARVMEC official repository](https://github.com/ORNL-Fusion/PARVMEC)
- [GVEC documentation](https://gvec.readthedocs.io/develop/index.html)
- [DESC official repository](https://github.com/PlasmaControl/DESC)
- [SPEC official documentation](https://princetonuniversity.github.io/SPEC/)
- [SIESTA ideal relaxation with islands and stochastic regions](https://doi.org/10.1063/1.3597155)
- [SIESTA free-plasma-boundary extension](https://doi.org/10.1063/1.4986447)
- [GPEC documentation and examples](https://princetonuniversity.github.io/GPEC/)
- [GPEC software metadata and ideal/resistive layer modes](https://doi.org/10.11578/GPEC/dc.20190207.1)
- [MARS-F response-model paper](https://doi.org/10.1016/j.cpc.2006.09.003)
- [GLISS repository](https://github.com/itpplasma/GLISS)
- [JOREK overview paper](https://arxiv.org/abs/2011.09120)
- [BOUT++ documentation](https://bout-dev.readthedocs.io/en/stable/user_docs/introduction.html)
- [GRILLIX official project page](https://physik.uni-greifswald.de/ag-manz/forschung/codes/grillix/)
- [GRILLIX FCI formulation](https://doi.org/10.1088/1361-6587/aaa373)
- [PARALLAX official repository](https://gitlab.mpcdf.mpg.de/phoenix-public/parallax)
- [PARALLAX/PAccX elliptic solver](https://arxiv.org/abs/2509.11831)
- [Mixed elasticity with weak stress symmetry](https://arxiv.org/abs/math/0701506)
- [JKU Linz TDNNS elasticity formulation](https://www.numa.uni-linz.ac.at/Talks/abstract/154/)
- [Anisotropic mixed finite elements for elasticity](https://doi.org/10.1002/nme.3319)
- [Symplectic mixed finite elements for acoustics](https://doi.org/10.1007/s00211-014-0667-4)
- [Port-Hamiltonian mixed finite elements for linear thermoelasticity](https://doi.org/10.1080/01495739.2021.1917322)
- [Structure-preserving linear elasticity example](https://publiweb.femto-st.fr/tntnet/entries/21775/documents/author/data)
- [HKT resonant current-sheet study](https://collaborate.princeton.edu/en/publications/numerical-study-of-%CE%B4-function-current-sheets-arising-from-resonan/)
- [HKT current sheets, MRxMHD, and the ideal nested limit](https://arxiv.org/abs/2108.09327)
- [Ideal shielding and resonant current sheets](https://conferences.iaea.org/event/98/contributions/11599/)
- [Manges and Cendes, generalized tree--cotree gauge](https://doi.org/10.1109/20.376275)
- [Bíró, Preis, and Richter, magnetic vector potential with tree--cotree gauge](https://ieeexplore.ieee.org/document/558631)
- [TEAM workshops and benchmark problems](https://www.osti.gov/servlets/purl/7179128)
- [TEAM Workshop 13 neutral benchmark description](https://docs.feelpp.org/toolboxes/latest/maxwell/Tws/index.html)
- [Fitzpatrick, bifurcated resistive-wall steady states](https://doi.org/10.1088/0029-5515/33/7/I08)
- [Matched-asymptotic resistive-layer formulation](https://collaborate.princeton.edu/en/publications/computation-of-resistive-instabilities-by-matched-asymptotic-expa/)
- [STARWALL and linear MHD stability](https://arxiv.org/abs/1508.04911)

### Free-boundary, vacuum, and conducting-wall operators

- [Merkel, NESTOR toroidal Neumann solver](https://doi.org/10.1016/0021-9991(86)90055-0)
- [Three-dimensional free-boundary calculations using a spectral Green function](https://doi.org/10.1016/0010-4655(86)90058-5)
- [Malhotra et al., BIEST generalized-Debye-source Taylor states](https://arxiv.org/abs/1902.01205)
- [O'Neil and Cerfon, integral-equation Beltrami/Taylor states](https://arxiv.org/abs/1611.01420)
- [Lazerson, virtual-casing principle for 3D toroidal systems](https://doi.org/10.1088/0741-3335/54/12/122002)
- [Toler, Cerfon, and Malhotra, direct virtual-casing field](https://doi.org/10.1017/S0022377824000527)
- [Conlin et al., high-order free-boundary equilibria in DESC](https://arxiv.org/abs/2412.05680)
- [VMEC++ free-boundary method documentation](https://proximafusion.github.io/vmecpp/api/vmecpp.html)
- [SPEC free-boundary virtual-casing documentation](https://princetonuniversity.github.io/SPEC/group__grp__free-boundary.html)
- [FEM--BEM coupling for axisymmetric free-boundary equilibrium](https://doi.org/10.1016/j.jcp.2017.06.006)
- [Hoelzl et al., JOREK--STARWALL coupling](https://arxiv.org/abs/1206.2748)
- [Cipolletta et al., compressed JOREK wall response matrices](https://arxiv.org/abs/2404.16546)
- [JOREK overview and vacuum response matrices](https://doi.org/10.1088/1741-4326/abf99f)

These references motivate neutral operator contracts only. Their profiles,
coordinate conventions, input formats, plasma closures, and production
equilibrium algorithms remain external application responsibilities.

### Structure-preserving discretization

- [Splitting-based structure-preserving MHD discretization](https://doi.org/10.5802/smai-jcm.34)
- [Structure-preserving and helicity-conserving FEM for MHD](https://doi.org/10.1016/j.jcp.2021.110894)
- [Structure-preserving Hall-MHD finite elements](https://arxiv.org/abs/2202.11586)
- [Variational integrators in plasma physics](https://arxiv.org/abs/1307.5665)
- [Discrete variational integration for ideal MHD](https://www.osti.gov/servlets/purl/1179782/)
- [Structure-preserving transport-stabilized compatible FEM for MHD](https://doi.org/10.1016/j.jcp.2024.112737)

### FortFEM provenance

- [FortFEM differentiation design](doc/design/differentiation.md)
- [FortFEM linear-response interchange contract](doc/design/linear_response_interchange.md)
- [FortFEM FEEC, MEPHIT, and open-boundary audit](doc/roadmap_mephit_feec_bem.md)
- [FortFEM interoperability benchmark protocol](doc/interoperability_benchmarks.md)
- [FortFEM generated-kernel checker](tools/codegen/check_generated.sh)

When a future implementation uses a formula, convention, or compatibility
fixture from another code, add an equation-level entry to the relevant
provenance file and record the external revision. A literature citation alone
does not document a software dependency or justify source reuse.
