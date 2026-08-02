#!/usr/bin/env python3
"""Check the staged canonical imports used by representative gallery examples.

This is deliberately a small, dependency-free structural gate.  It resolves
the expected facade modules from their source declarations, then parses the
example programs' Fortran ``use`` statements.  Behavioral equivalence is
covered separately by ``test_api05_example_facades``; this checker only keeps
the migration map honest and actionable.
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path


USE_RE = re.compile(
    r"^\s*use\s+(?:(?:,\s*intrinsic|,\s*non_intrinsic)\s*::\s*)?"
    r"(?:::\s*)?([a-z][a-z0-9_]*)\b",
    re.IGNORECASE,
)
MODULE_RE = re.compile(r"^\s*module\s+([a-z][a-z0-9_]*)\b", re.IGNORECASE)


@dataclass(frozen=True)
class ExampleCase:
    name: str
    path: str
    facade: str
    legacy_smoke: bool = False
    required_facades: tuple[str, ...] = ()


CASES = (
    ExampleCase(
        "scalar Poisson compatibility smoke",
        "example/simple_poisson/simple_poisson.f90",
        "fortfem_api",
        legacy_smoke=True,
    ),
    ExampleCase(
        "vector Nedelec solution",
        "example/tetra_nedelec_p_convergence.f90",
        "fortfem_feec",
    ),
    ExampleCase(
        "exterior Laplace BEM solution",
        "example/laplace_exterior_bem_circle/laplace_exterior_bem_circle.f90",
        "fortfem_boundary",
    ),
    ExampleCase(
        "circular Helmholtz DtN solution",
        "example/circular_dtn_modes.f90",
        "fortfem_boundary",
    ),
    ExampleCase(
        "mixed acoustic wave solution",
        "example/mixed_acoustic_wave/mixed_acoustic_wave.f90",
        "fortfem_time",
    ),
    ExampleCase(
        "three-dimensional mixed wave solution",
        "example/mixed_wave_3d_structure/mixed_wave_3d_structure.f90",
        "fortfem_time",
    ),
    ExampleCase(
        "mixed pressure-displacement wave solution",
        "example/mixed_wave_pressure_displacement_gallery/"
        "mixed_wave_pressure_displacement_gallery.f90",
        "fortfem_time",
    ),
    ExampleCase(
        "Helmholtz circle BEM spectrum",
        "example/helmholtz_bem_circle_spectrum.f90",
        "fortfem_boundary",
    ),
    ExampleCase(
        "Laplace circle BEM spectrum",
        "example/laplace_bem_circle_spectrum.f90",
        "fortfem_boundary",
    ),
    ExampleCase(
        "Helmholtz circle CFIE solution",
        "example/helmholtz_cfie_circle.f90",
        "fortfem_boundary",
    ),
    ExampleCase(
        "Laplace symmetric FEM-BEM transmission",
        "example/laplace_symmetric_transmission.f90",
        "fortfem_boundary",
    ),
    ExampleCase(
        "neutral free-boundary trace",
        "example/free_boundary_port_gallery/free_boundary_port_gallery.f90",
        "fortfem_boundary",
    ),
    ExampleCase(
        "Bíró 3-D tree-cotree solution",
        "example/biro_tree_cotree_3d_gallery/biro_tree_cotree_3d_gallery.f90",
        "fortfem_feec",
    ),
    ExampleCase(
        "field-aligned anisotropic 3-D solution",
        "example/field_aligned_anisotropic_3d_gallery/field_aligned_anisotropic_3d_gallery.f90",
        "fortfem_feec",
    ),
    ExampleCase(
        "CGL tensor-pressure 3-D solution",
        "example/cgl_pressure_3d_gallery/cgl_pressure_3d_gallery.f90",
        "fortfem_feec",
    ),
    ExampleCase(
        "polar IGA FEEC solution",
        "example/iga_polar_feec/iga_polar_feec.f90",
        "fortfem_feec",
    ),
    ExampleCase(
        "toroidal Maxwell FEM-BEM solution",
        "example/maxwell_torus_fem_bem_solution/maxwell_torus_fem_bem_solution.f90",
        "fortfem_boundary",
        required_facades=("fortfem_boundary", "fortfem_core", "fortfem_feec"),
    ),
    ExampleCase(
        "curved toroidal Maxwell scattering solution",
        "example/maxwell_torus_curved_scattering/maxwell_torus_curved_scattering.f90",
        "fortfem_boundary",
        required_facades=("fortfem_boundary", "fortfem_core"),
    ),
    ExampleCase(
        "three-dimensional sphere BEM solution",
        "example/bem_sphere_3d/bem_sphere_3d.f90",
        "fortfem_boundary",
        required_facades=("fortfem_boundary", "fortfem_core"),
    ),
    ExampleCase(
        "adaptive Laplace BEM solution",
        "example/adaptive_bem_prolate/adaptive_bem_prolate.f90",
        "fortfem_boundary",
        required_facades=("fortfem_boundary", "fortfem_core"),
    ),
    ExampleCase(
        "adaptive Helmholtz BEM solution",
        "example/adaptive_helmholtz_bem_sphere/"
        "adaptive_helmholtz_bem_sphere.f90",
        "fortfem_boundary",
        required_facades=("fortfem_boundary", "fortfem_core"),
    ),
    ExampleCase(
        "planar acoustic FEM-DtN solution",
        "example/acoustic_fem_dtn/acoustic_fem_dtn.f90",
        "fortfem_boundary",
    ),
    ExampleCase(
        "curved acoustic NtD solution",
        "example/curved_acoustic_ntd/curved_acoustic_ntd.f90",
        "fortfem_boundary",
    ),
    ExampleCase(
        "scalar Helmholtz PML comparison",
        "example/helmholtz_open_boundary_comparison/"
        "helmholtz_open_boundary_comparison.f90",
        "fortfem_boundary",
    ),
    ExampleCase(
        "Maxwell DtN/PML comparison",
        "example/maxwell_open_boundary_comparison/"
        "maxwell_open_boundary_comparison.f90",
        "fortfem_boundary",
        required_facades=("fortfem_boundary", "fortfem_core", "fortfem_feec"),
    ),
)


def module_index(root: Path) -> dict[str, Path]:
    """Resolve Fortran module names to their defining source files."""

    index: dict[str, Path] = {}
    for path in sorted((root / "src").rglob("*.f90")):
        for line in path.read_text(encoding="utf-8").splitlines():
            match = MODULE_RE.match(line)
            if not match or line.lower().lstrip().startswith("module procedure"):
                continue
            name = match.group(1).lower()
            previous = index.get(name)
            if previous is not None and previous != path:
                raise ValueError(
                    f"module {name} has multiple providers: {previous} and {path}"
                )
            index[name] = path
    return index


def imported_modules(path: Path) -> set[str]:
    modules: set[str] = set()
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.split("!", 1)[0]
        match = USE_RE.match(line)
        if match:
            modules.add(match.group(1).lower())
    return modules


def check(root: Path) -> list[str]:
    errors: list[str] = []
    try:
        providers = module_index(root)
    except ValueError as error:
        return [str(error)]

    legacy_count = 0
    for case in CASES:
        source = root / case.path
        if not source.is_file():
            errors.append(f"{case.name}: missing source {case.path}")
            continue
        modules = imported_modules(source)
        facade = case.facade.lower()
        required_facades = tuple(
            name.lower() for name in (case.required_facades or (case.facade,))
        )
        for required in required_facades:
            provider = providers.get(required)
            if provider is None:
                errors.append(
                    f"{case.name}: facade {required} has no module provider"
                )
            else:
                expected = root / "src" / f"{required}.f90"
                if provider != expected:
                    errors.append(
                        f"{case.name}: facade {required} resolves to {provider}, "
                        f"expected {expected}"
                    )
            if required not in modules:
                errors.append(f"{case.name}: missing canonical import {required}")
        if case.legacy_smoke:
            legacy_count += 1
            if facade != "fortfem_api":
                errors.append(f"{case.name}: legacy smoke must use fortfem_api")
        elif "fortfem_api" in modules:
            errors.append(
                f"{case.name}: migrated example still imports the umbrella fortfem_api"
            )

    if legacy_count != 1:
        errors.append(
            "exactly one representative example must retain the umbrella smoke path"
        )
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    root = args.root.resolve()
    errors = check(root)
    if errors:
        for error in errors:
            print(f"ERROR: {error}")
        return 1
    print(f"example facade imports: PASS ({len(CASES)} representative cases)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
