#!/usr/bin/env python3
"""Check the compiler-visible Fortran module dependency graph.

This is intentionally a small, dependency-free release-gate check.  It reads
the same ``module`` declarations and ``use`` associations that determine the
compiler's ``.mod`` dependencies; it does not infer dependencies from public
procedure names or from filenames alone.

The current tree has a transitional set of ``fortfem_api_*`` facades.  The
checker therefore enforces only edges which are unambiguously invalid today:

* the convenience umbrella ``fortfem_api`` is never imported internally;
* generated sources never import an API facade;
* core sources never import application-layer modules or API facades; and
* plotting stays behind the plotting facade.

The full local graph is still built for duplicate-provider and cycle checks,
so this can be tightened as the facade refactor progresses without replacing
the checker with a grep-based heuristic.
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path


FORTRAN_SUFFIXES = {".f", ".f90", ".f95", ".f03", ".f08"}
FACADE_PREFIX = "fortfem_api_"
UMBRELLA = "fortfem_api"
PLOT_MODULES = {"fortplot", "fortplot_figure"}

# This is deliberately a conservative map.  Existing assembly and operator
# implementations are still being split, so only the unambiguous lowest and
# highest edges are rejected at present.
DIRECTORY_LAYERS = {
    "core": 0,
    "utils": 0,
    "topology": 0,
    "geometry": 0,
    "mesh": 0,
    "quadrature": 0,
    "triangulation": 0,
    "special": 0,
    "basis": 1,
    "elements": 1,
    "interpolation": 1,
    "assembly": 1,
    "operators": 1,
    "enrichment": 1,
    "spaces": 1,
    "boundary": 2,
    "time": 2,
    "bem": 3,
    "solvers": 3,
    "interoperability": 3,
    "plot": 3,
    "bindings": 3,
}

MODULE_RE = re.compile(
    r"^\s*module\s+(?!procedure\b)(?!subroutine\b)"
    r"(?!function\b)([a-z_]\w*)\b",
    re.IGNORECASE,
)
# The module name is the first token after USE, whether the statement uses
# ``use m``, ``use :: m`` or ``use, intrinsic :: m``.
USE_RE = re.compile(
    r"^\s*use\b(?:\s*,\s*[a-z_]\w*)*\s*(?:::)?\s*"
    r"([a-z_]\w*)\b",
    re.IGNORECASE,
)


@dataclass(frozen=True)
class UseRecord:
    name: str
    line: int


@dataclass(frozen=True)
class ModuleRecord:
    name: str
    path: Path
    uses: tuple[UseRecord, ...]
    layer: int | None
    kind: str


def strip_comment(line: str) -> str:
    """Remove a Fortran comment while respecting quoted strings."""

    quote: str | None = None
    index = 0
    while index < len(line):
        char = line[index]
        if quote is None and char in "'\"":
            quote = char
        elif quote is not None and char == quote:
            # Fortran escapes a quote by doubling it.
            if index + 1 < len(line) and line[index + 1] == quote:
                index += 1
            else:
                quote = None
        elif quote is None and char == "!":
            return line[:index]
        index += 1
    return line


def classify(path: Path, root: Path, module: str) -> tuple[int | None, str]:
    """Return (layer, kind) for a source module."""

    relative = path.relative_to(root)
    parts = {part.lower() for part in relative.parts[:-1]}
    if "generated" in parts or module.startswith("fortfem_generated_"):
        return -1, "generated"
    if module == UMBRELLA:
        return 4, "umbrella"
    if module.startswith(FACADE_PREFIX):
        return 4, "facade"
    for part in relative.parts[:-1]:
        layer = DIRECTORY_LAYERS.get(part.lower())
        if layer is not None:
            return layer, "source"
    return None, "source"


def discover(root: Path) -> tuple[list[ModuleRecord], list[str]]:
    records: list[ModuleRecord] = []
    diagnostics: list[str] = []
    paths = sorted(
        path
        for path in root.rglob("*")
        if path.is_file() and path.suffix.lower() in FORTRAN_SUFFIXES
    )
    for path in paths:
        modules: list[str] = []
        uses: list[UseRecord] = []
        try:
            lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        except OSError as error:
            diagnostics.append(f"{path}: cannot read source: {error}")
            continue
        for line_number, raw_line in enumerate(lines, start=1):
            line = strip_comment(raw_line)
            module_match = MODULE_RE.match(line)
            if module_match:
                modules.append(module_match.group(1).lower())
            use_match = USE_RE.match(line)
            if use_match:
                uses.append(UseRecord(use_match.group(1).lower(), line_number))
        for module in modules:
            layer, kind = classify(path, root, module)
            records.append(ModuleRecord(module, path, tuple(uses), layer, kind))
    return records, diagnostics


def display_path(path: Path, root: Path) -> str:
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


def check(root: Path) -> list[str]:
    records, diagnostics = discover(root)
    providers: dict[str, ModuleRecord] = {}
    for record in records:
        previous = providers.get(record.name)
        if previous is not None:
            diagnostics.append(
                "duplicate module provider: "
                f"{record.name} in {display_path(previous.path, root)} and "
                f"{display_path(record.path, root)}"
            )
        else:
            providers[record.name] = record

    edges: dict[str, list[str]] = {record.name: [] for record in records}
    for record in records:
        for use in record.uses:
            target = use.name
            target_record = providers.get(target)
            if target_record is not None:
                edges[record.name].append(target)

            if target == UMBRELLA and record.kind != "umbrella":
                diagnostics.append(
                    f"{display_path(record.path, root)}:{use.line}: "
                    "forbidden umbrella import: use fortfem_api; "
                    "internal modules must depend on a narrow facade or core module"
                )
            if target.startswith(FACADE_PREFIX):
                if record.kind == "generated":
                    diagnostics.append(
                        f"{display_path(record.path, root)}:{use.line}: "
                        f"forbidden generated-layer import of facade {target}"
                    )
                elif record.layer == 0:
                    diagnostics.append(
                        f"{display_path(record.path, root)}:{use.line}: "
                        f"forbidden core-layer import of facade {target}"
                    )
            if target in PLOT_MODULES and not record.name.startswith(
                "fortfem_api_plot"
            ):
                diagnostics.append(
                    f"{display_path(record.path, root)}:{use.line}: "
                    f"plotting dependency {target} must stay behind fortfem_api_plot"
                )
            if (
                record.layer == 0
                and target_record is not None
                and target_record.layer == 3
            ):
                diagnostics.append(
                    f"{display_path(record.path, root)}:{use.line}: "
                    f"layer violation: core module {record.name} imports "
                    f"application-layer module {target}"
                )

    # DFS over provider modules catches compiler-visible cycles, including
    # cycles which happen to be hidden behind private/non-public procedures.
    state: dict[str, int] = {}
    stack: list[str] = []
    cycle_keys: set[tuple[str, ...]] = set()

    def visit(module: str) -> None:
        state[module] = 1
        stack.append(module)
        for target in sorted(set(edges[module])):
            if state.get(target, 0) == 0:
                visit(target)
            elif state.get(target) == 1:
                start = stack.index(target)
                cycle = tuple(stack[start:] + [target])
                canonical = min(cycle[i:-1] + cycle[:i] for i in range(len(cycle) - 1))
                if canonical not in cycle_keys:
                    cycle_keys.add(canonical)
                    diagnostics.append(
                        "module dependency cycle: "
                        + " -> ".join(cycle)
                    )
        stack.pop()
        state[module] = 2

    for module in sorted(edges):
        if state.get(module, 0) == 0:
            visit(module)
    return sorted(set(diagnostics))


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Audit Fortran module/use dependencies and API layer leaks."
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("src"),
        help="source tree to scan (default: src)",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    root = args.root.resolve()
    if not root.is_dir():
        print(f"error: source root does not exist: {root}", file=sys.stderr)
        return 2
    diagnostics = check(root)
    if diagnostics:
        for diagnostic in diagnostics:
            print(f"ERROR: {diagnostic}")
        print(f"module-layer audit failed: {len(diagnostics)} diagnostic(s)")
        return 1
    print(f"module-layer audit passed: {root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
