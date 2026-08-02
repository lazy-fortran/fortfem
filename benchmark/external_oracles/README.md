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
```

The second command intentionally prints `SKIP` and exits successfully while the
external manifest is absent. With a reviewed sister-repository checkout, pass
`--contract`, `--data-manifest`, and `--data-root`; the adapter verifies the
case/provenance/member contract and archive bytes before reporting `READY`.
It never downloads or invokes an external solver.
