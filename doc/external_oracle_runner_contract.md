---
title: Optional external-oracle runner contract
---

# Optional external-oracle runner contract

`benchmark/external_oracles/runner_manifest.json` is the small, machine-readable
contract for running the shared analytical Poisson case in FreeFEM, MFEM, or
FEniCSx.  It is intentionally metadata-only: no third-party source, binary,
container, mesh, solution dump, or plot is part of FortFEM.

Validate the contract without installing any external code:

```sh
python3 tools/validate_external_oracle_runner.py \
  benchmark/external_oracles/runner_manifest.json
bash test/test_external_oracle_runner_contract.sh
```

The contract records, for each optional runner:

- the exact launch argument vector, working directory, timeout, source URI,
  expected license, and immutable revision (when a run actually occurred);
- the physical-coordinate sampling rule and shared nine-point analytical
  Poisson oracle;
- coordinate/sample SHA-256 fields and the sister benchmark-data repository
  URI;
- absolute, relative, coordinate, and residual tolerances;
- mesh, assembly, factorization, solve, output, total time, warmups,
  repetitions, and peak memory;
- an explicit `skip_reason` whenever an external executable was not run.

Skipped entries keep timing, revision, and result checksums null.  This avoids
inventing provenance while making omission visible to a publication job.  A
successful entry must provide all of those fields and an immutable artifact
manifest URI in the separate benchmark-data repository.

FEniCSx/DOLFINx is the only runner that advertises `uv`: `uv` can manage
lightweight CSV checking or plotting after an official DOLFINx container run;
it cannot replace DOLFINx's compiled C++/PETSc stack.  FreeFEM and MFEM remain
native optional executables.  The validator never launches any of them and
completes in milliseconds, so the ordinary FortFEM test path stays lightweight.

The nine-point case is an independent behavioral fixture, not an external
result.  Its exact field is

\[
u(x,y)=\sin(\pi x)\sin(\pi y),\qquad
-\Delta u=2\pi^2\sin(\pi x)\sin(\pi y),
\]

with zero boundary data.  After review, a scheduled runner may replace a
skipped record with real revision, checksum, sampling, and performance values
and publish the resulting manifest and samples at an immutable sister-repo
commit.  The source adapters already present under `benchmark/freefem`,
`benchmark/mfem`, and `benchmark/fenicsx` remain independent implementations;
this contract does not make FortFEM depend on their licenses or formats.
