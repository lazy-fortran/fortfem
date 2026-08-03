---
title: Provenance and external benchmark data policy
---

# Provenance and external benchmark data policy

FortFEM provides reusable discretization foundations, not a public mirror of
papers, workshop handouts, private source trees, or proprietary benchmark
archives. The public repository keeps the equations, small manufactured
fixtures, adapters, validators, and provenance needed to reproduce a method.
Exact application data are an optional, separately licensed input.

## What may be published

After a human rights review, a case manifest may contain:

- geometry topology and dimensions;
- region, material, and interface labels;
- units, normalization, excitation, and boundary/interface conditions;
- analytical formulas and FortSym-generated source terms;
- probe coordinates, sampling weights, and a generated mesh when its license
  permits redistribution;
- public DOI, report, or official documentation links.

These fields describe the computational problem. They do not reproduce a
paper's protected expression, figures, tables, or raw reference data.

## What stays private or external

Publisher PDFs, scans, copyrighted plots and tables, private meshes, source
scripts, solver output, unlicensed B--H curves, and raw measurement/reference
arrays are never copied to FortFEM. They may remain in a local Nextcloud or
Zotero library, or be published in the separate benchmark-data repository only
under a license that permits it. The checked-in
[`source_registry.json`](../benchmark/external_oracles/source_registry.json)
contains logical source IDs and rights states, not private file paths.

An exact external run requires all of the following:

1. a reviewed source and redistribution status;
2. a versioned payload schema and immutable HTTPS manifest;
3. a SHA-256 checksum for the payload or archive;
4. code, executable revision, compiler, hardware, units, normalization, and
   sampling metadata;
5. an independent solution/error oracle and a declared skip when any item is
   missing.

Absence is a valid result. The adapter must print `SKIP`, and the gallery must
continue to show the small manufactured or neutral fixture rather than
fabricating an exact comparison.

## Paper Magnetic and TEAM hand-off

The private inventory currently contains TEAM problem references for 1A, 1B,
2, 6, 13, and 20. These are all required future acceptance fixtures. For each
case, the implementation proceeds in four separate artifacts:

1. a public geometry/BC manifest containing only reviewed computational
   metadata;
2. a FortSym-generated manufactured or neutral fixture for the fast test path;
3. an optional exact payload in the sister data repository;
4. a Pages gallery record with system and solution plots before diagnostics.

The exact payload remains disabled until its license and checksum are
recorded. A private Biro item whose title or method does not establish the
tree--cotree target is recorded as a candidate, not used as evidence for the
tree--cotree gallery.

## Solution-first gallery contract

Every example page must make the physical computation inspectable. The first
panel is the complete system: mesh or spline patches, element edges and
orientations, regions, coils, interfaces, and boundary labels. The second
panel is the computed solution on that geometry: scalar contours, vector
arrows, a surface field, or field lines as appropriate. Comparison and
diagnostic plots follow. A generated record carries `solution_first: true`,
the physical coordinate system, units, normalization, and the exact source
and data revisions used.

Generated images remain ignored build products. Pages regenerates them from
the pinned source/data revisions; CSV and JSON records remain the durable
machine-readable artifacts. A plot that shows only an abstract convergence
curve does not satisfy this contract.

## References

The public provenance links used by the current adapters are:

- [TEAM workshop catalogue and report](https://www.osti.gov/biblio/7179128);
- [TEAM Workshop 13 reference documentation](https://docs.feelpp.org/toolboxes/latest/maxwell/Tws/index.html);
- [TEAM-20 public validation description](https://www.simscale.com/docs/validation-cases/team-20-magnetostatics/);
- [Bíró, Preis, and Richter tree--cotree formulation](https://ieeexplore.ieee.org/document/558631);
- [Bíró, Preis, and Richter vector-potential formulation](https://doi.org/10.1109/20.497322);
- [Manges and Cendes generalized tree--cotree gauge](https://doi.org/10.1109/20.376275);
- [Park et al. resonant-layer preprint](https://conferences.iaea.org/event/392/papers/35712/files/13694-ParkJK_FEC2025.pdf);
- [Waltz and Waelbroeck resonant perturbations](https://doi.org/10.1063/1.3692222).

These links support provenance and method selection. They do not grant
redistribution rights for attached files or numerical data.
