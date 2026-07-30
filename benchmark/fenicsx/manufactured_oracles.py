"""Independent FEniCSx runners for the shared Poisson and Ampere cases."""

from __future__ import annotations

import csv
import math
import time
from pathlib import Path

from mpi4py import MPI
from petsc4py import PETSc

import numpy as np
import ufl
from dolfinx import default_scalar_type, fem, mesh
from dolfinx.fem.petsc import LinearProblem


LEVELS = (4, 8, 16, 32, 64)
OUTPUT = Path("benchmark-results/fenicsx")


def global_scalar(local_value: float) -> float:
    return MPI.COMM_WORLD.allreduce(local_value, op=MPI.SUM)


def mesh_counts(domain: mesh.Mesh) -> tuple[int, int]:
    cells = domain.topology.index_map(domain.topology.dim).size_global
    return cells, domain.topology.dim - 1


def boundary_facets(domain: mesh.Mesh) -> np.ndarray:
    facet_dimension = domain.topology.dim - 1
    return mesh.locate_entities_boundary(
        domain,
        facet_dimension,
        lambda coordinates: np.full(coordinates.shape[1], True),
    )


def dof_count(space: fem.FunctionSpace) -> int:
    return (
        space.dofmap.index_map.size_global
        * space.dofmap.index_map_bs
    )


def solve_poisson(n: int) -> dict[str, float | int | str]:
    domain = mesh.create_unit_square(
        MPI.COMM_WORLD, n, n, cell_type=mesh.CellType.triangle
    )
    space = fem.functionspace(domain, ("Lagrange", 1))
    trial = ufl.TrialFunction(space)
    test = ufl.TestFunction(space)
    coordinate = ufl.SpatialCoordinate(domain)
    exact = ufl.sin(math.pi * coordinate[0]) * ufl.sin(
        math.pi * coordinate[1]
    )
    source = 2.0 * math.pi**2 * exact
    facets = boundary_facets(domain)
    boundary_dofs = fem.locate_dofs_topological(
        space, domain.topology.dim - 1, facets
    )
    condition = fem.dirichletbc(
        default_scalar_type(0.0), boundary_dofs, space
    )
    bilinear = ufl.inner(ufl.grad(trial), ufl.grad(test)) * ufl.dx
    linear = source * test * ufl.dx

    total_start = time.perf_counter()
    assembly_start = time.perf_counter()
    problem = LinearProblem(
        bilinear,
        linear,
        bcs=[condition],
        petsc_options_prefix=f"poisson_n{n}_",
        petsc_options={"ksp_type": "preonly", "pc_type": "lu"},
    )
    assembly_seconds = time.perf_counter() - assembly_start
    solve_start = time.perf_counter()
    solution = problem.solve()
    solve_seconds = time.perf_counter() - solve_start
    total_seconds = time.perf_counter() - total_start

    l2_squared = global_scalar(
        fem.assemble_scalar(fem.form((solution - exact) ** 2 * ufl.dx))
    )
    gradient_squared = global_scalar(
        fem.assemble_scalar(
            fem.form(
                ufl.inner(
                    ufl.grad(solution - exact),
                    ufl.grad(solution - exact),
                )
                * ufl.dx
            )
        )
    )
    cells, _ = mesh_counts(domain)
    return {
        "mesh_id": f"unit-square-n{n}",
        "cells": cells,
        "dofs": dof_count(space),
        "h": 1.0 / n,
        "l2_error": math.sqrt(l2_squared),
        "h1_error": math.sqrt(l2_squared + gradient_squared),
        "assembly_seconds": assembly_seconds,
        "solve_seconds": solve_seconds,
        "total_seconds": total_seconds,
    }


def solve_ampere(n: int) -> dict[str, float | int | str]:
    domain = mesh.create_unit_square(
        MPI.COMM_WORLD, n, n, cell_type=mesh.CellType.triangle
    )
    space = fem.functionspace(domain, ("N1curl", 1))
    trial = ufl.TrialFunction(space)
    test = ufl.TestFunction(space)
    coordinate = ufl.SpatialCoordinate(domain)
    exact = ufl.as_vector(
        (
            ufl.sin(math.pi * coordinate[1]),
            ufl.sin(math.pi * coordinate[0]),
        )
    )
    source = (math.pi**2 + 1.0) * exact
    facets = boundary_facets(domain)
    boundary_dofs = fem.locate_dofs_topological(
        space, domain.topology.dim - 1, facets
    )
    zero = fem.Function(space)
    zero.x.array[:] = 0.0
    condition = fem.dirichletbc(zero, boundary_dofs)
    bilinear = (
        ufl.inner(ufl.curl(trial), ufl.curl(test))
        + ufl.inner(trial, test)
    ) * ufl.dx
    linear = ufl.inner(source, test) * ufl.dx

    total_start = time.perf_counter()
    assembly_start = time.perf_counter()
    problem = LinearProblem(
        bilinear,
        linear,
        bcs=[condition],
        petsc_options_prefix=f"ampere_n{n}_",
        petsc_options={"ksp_type": "preonly", "pc_type": "lu"},
    )
    assembly_seconds = time.perf_counter() - assembly_start
    solve_start = time.perf_counter()
    solution = problem.solve()
    solve_seconds = time.perf_counter() - solve_start
    total_seconds = time.perf_counter() - total_start

    difference = solution - exact
    l2_squared = global_scalar(
        fem.assemble_scalar(
            fem.form(ufl.inner(difference, difference) * ufl.dx)
        )
    )
    curl_squared = global_scalar(
        fem.assemble_scalar(
            fem.form(
                (ufl.curl(solution) - ufl.curl(exact)) ** 2 * ufl.dx
            )
        )
    )
    cells, _ = mesh_counts(domain)
    return {
        "mesh_id": f"unit-square-n{n}",
        "cells": cells,
        "dofs": dof_count(space),
        "h": 1.0 / n,
        "l2_error": math.sqrt(l2_squared),
        "hcurl_error": math.sqrt(l2_squared + curl_squared),
        "assembly_seconds": assembly_seconds,
        "solve_seconds": solve_seconds,
        "total_seconds": total_seconds,
    }


def write_records(name: str, records: list[dict[str, object]]) -> None:
    if MPI.COMM_WORLD.rank != 0:
        return
    OUTPUT.mkdir(parents=True, exist_ok=True)
    with (OUTPUT / f"{name}.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def require_decreasing(
    records: list[dict[str, object]], fields: tuple[str, ...]
) -> None:
    for field in fields:
        values = [float(record[field]) for record in records]
        if any(fine >= coarse for coarse, fine in zip(values, values[1:])):
            raise RuntimeError(f"{field} did not decrease: {values}")


def main() -> None:
    poisson = [solve_poisson(n) for n in LEVELS]
    ampere = [solve_ampere(n) for n in LEVELS]
    require_decreasing(poisson, ("l2_error", "h1_error"))
    require_decreasing(ampere, ("l2_error", "hcurl_error"))
    write_records("poisson", poisson)
    write_records("ampere", ampere)


if __name__ == "__main__":
    PETSc.Sys.popErrorHandler()
    main()
