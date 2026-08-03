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
solution values are rendered by their magnitude. The output metadata labels the
plot as an `external-benchmark-payload` and keeps the analytical reference
explicitly provenance-only; it is not a manufactured or bundled analytical
solution. Elements with repeated nodes or zero area are rejected before
rendering. This is a plotting adapter, not a claim that FortFEM contains the
paper's application data or solver.

The checked-in contract remains `SKIP`; the in-tree
`example/biro_tree_cotree_3d_gallery` is the manufactured, method-faithful
tree--cotree oracle and is labeled separately.

## Private source registry and public extraction boundary

`source_registry.json` records logical source IDs for the private Nextcloud and
Zotero material that has been inventoried for this work. It is deliberately
not a file index: it contains no absolute local paths, PDF bytes, scans,
source scripts, meshes, reference arrays, or copyrighted figures. A private
attachment remains private even when its bibliographic record or public DOI is
listed.

The Paper Magnetic cases are required acceptance fixtures. The current
inventory includes private references for TEAM-1A, TEAM-1B, TEAM-2, TEAM-6,
TEAM-13, and TEAM-20. Their exact data are not yet enabled in FortFEM. The
public extraction task is limited to a reviewed transcription of geometry,
region/material labels, units, excitations, boundary and interface conditions,
normalization, and probe coordinates. Those fields can be published as a
small case manifest when their source permits it; raw scans, B--H tables,
measurement arrays, private meshes, and solver scripts stay in the private
vault or the separately licensed benchmark-data repository.

The reduced-scalar-potential item in the private Biro inventory is recorded as
a candidate only. It is not silently treated as the published tree--cotree
case. The tree--cotree adapter remains tied to its own public provenance and
will be enabled only after the exact target and rights are confirmed.

For every exported case, the Pages job follows the same sequence:

1. render the system first (complete mesh or spline patches, regions, coils,
   interfaces, and boundary labels);
2. render the computed solution (scalar contours, vector arrows, surface
   fields, or field lines) on the same physical geometry;
3. add analytical or external-reference overlays;
4. append error, convergence, conservation, derivative, and performance
   panels.

The gallery writes SVG/PNG as ignored build artifacts and CSV/JSON containing
the case ID, source ID, sampler, units, normalization, executable revision,
checksums, and skip reason. Exact-data adapters must validate the immutable
manifest and SHA-256 before they can emit a solution plot. The complete policy
and hand-off checklist is in
[`doc/provenance_and_data_policy.md`](../../doc/provenance_and_data_policy.md).

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
under `--output-dir/TEAM-*`. Metadata identifies each output as an
`external-benchmark-payload` with a provenance-only analytical reference, and
degenerate elements are rejected before a misleading plot can be emitted; the
exact arrays remain outside FortFEM.
