# Exact Biro-data adapter boundary

The repository contains a method-faithful manufactured tree--cotree gallery
(`example/biro_tree_cotree_3d_gallery`) and its independent orientation oracle.
It does **not** contain the application-specific geometry, material, or source
arrays used by Biro, Preis, and Richter. The exact paper provenance is recorded
at [DOI 10.1109/20.497322](https://doi.org/10.1109/20.497322).

`biro_paper_manifest.json` is a metadata-only contract. It names the expected
array members and points to the sister repository
`lazy-fortran/fortfem-benchmark-data`, but the checked-in record is explicitly
`absent` and therefore carries no artifact checksum. A reviewed external
artifact can be enabled only by publishing a new immutable manifest with its
HTTPS URI and SHA-256 checksum; no paper data or solver source is copied into
FortFEM.

The validator and adapter use only the Python standard library:

```text
python3 tools/validate_biro_external_manifest.py \
  benchmark/external_oracles/biro_paper_manifest.json
python3 benchmark/external_oracles/run_biro_paper_adapter.py
python3 benchmark/external_oracles/run_biro_paper_gallery.py
```

The second command intentionally prints `SKIP` and exits successfully while the
external manifest is absent. With a reviewed sister-repository checkout, pass
`--contract`, `--data-manifest`, and `--data-root`; the adapter verifies the
case/provenance/member contract and archive bytes before reporting `READY`.
It never downloads or invokes an external solver.

The gallery command uses the same gate. When a reviewed artifact supplies a
JSON payload with schema `fortfem-biro-paper-payload-1` (nodes, elements, and
postprocessed solution values), it writes `solution.svg`, `solution.csv`, and
`provenance.json` under `--output-dir`. The SVG deliberately puts the verified
solution first and labels the DOI, checksum, and exact-data status. Vector
solution values are rendered by their magnitude. This is a plotting adapter,
not a claim that FortFEM contains the paper's application data or solver.

The checked-in contract remains `SKIP`; the in-tree
`example/biro_tree_cotree_3d_gallery` is the manufactured, method-faithful
tree--cotree oracle and is labeled separately.

## TEAM exact-data boundary

The four TEAM cases used by the neutral galleries are tracked by
`team_manifest.json`: TEAM-3 (C-core/air gap), TEAM-7 (coil/plate
eddy-current), TEAM-13 (coil/channel/plate), and TEAM-20 (3-D solenoid
static force).  Their provenance points to the public [TEAM catalogue and
workshop report](https://www.osti.gov/biblio/7179128) (and, where useful, the
public TEAM-13 and TEAM-20 descriptions).  No TEAM mesh, B--H curve, source,
probe, or force array is copied into this repository.

The metadata-only contract is intentionally `availability: absent` and its
adapter therefore prints an explicit `SKIP`:

```text
python3 benchmark/external_oracles/run_team_adapter.py
python3 benchmark/external_oracles/run_team_gallery.py
```

An independently licensed checkout of
`lazy-fortran/fortfem-benchmark-data` can supply a pinned
`fortfem-team-data-1` manifest and archive.  The adapter verifies every TEAM
case ID, HTTPS provenance, payload schema, license declaration, repository URI,
and SHA-256 before the gallery consumes it.  It never downloads or invokes an
external solver.  With a verified payload, the gallery writes one
solution-first `solution.svg`, `solution.csv`, and `provenance.json` per case
under `--output-dir/TEAM-*`; the exact arrays remain outside FortFEM.
