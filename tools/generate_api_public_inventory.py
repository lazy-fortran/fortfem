#!/usr/bin/env python3
"""Emit the pre-release public Fortran API inventory.

This intentionally uses only the Python standard library.  The inventory is
derived from the API wrapper ``public`` and ``use ... only:`` statements; it
does not attempt to infer undocumented numerical semantics.
"""

from __future__ import annotations

import argparse
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
API_FILES = sorted((ROOT / "src").glob("fortfem_api*.f90"))
# Canonical facades are part of the public inventory once they exist.  Keep
# this list explicit so implementation modules are not accidentally treated
# as public just because they happen to contain a ``public ::`` statement.
CANONICAL_FACADE_NAMES = (
    "fortfem_core.f90",
    "fortfem_feec.f90",
    "fortfem_fourier.f90",
    "fortfem_boundary.f90",
    "fortfem_time.f90",
    "fortfem_interop.f90",
    "fortfem_plot.f90",
)
API_FILES.extend(
    path for name in CANONICAL_FACADE_NAMES
    if (path := ROOT / "src" / name).is_file()
)
API_FILES = sorted(set(API_FILES))


@dataclass(frozen=True)
class Export:
    symbol: str
    file: Path
    module: str
    line: int


@dataclass(frozen=True)
class Import:
    module: str
    name: str


def statements(path: Path):
    """Yield (line number, logical free-form statement) pairs."""
    pending = ""
    start = 0
    for number, raw in enumerate(path.read_text().splitlines(), 1):
        line = raw.split("!", 1)[0].rstrip()
        if not pending:
            pending = line
            start = number
        else:
            pending += " " + line.lstrip("& ")
        if pending.rstrip().endswith("&"):
            pending = pending.rstrip()[:-1].rstrip()
        else:
            if pending.strip():
                yield start, pending.strip()
            pending = ""


def module_name(path: Path) -> str:
    for _, statement in statements(path):
        match = re.match(r"module\s+(?!procedure\b)([a-z_]\w*)", statement, re.I)
        if match:
            return match.group(1).lower()
    return path.stem.lower()


def split_names(text: str) -> list[str]:
    return [name.strip() for name in text.split(",") if name.strip()]


def parse_api(path: Path):
    module = module_name(path)
    exports: list[Export] = []
    imports: dict[str, Import] = {}
    wildcards: list[str] = []
    for line, statement in statements(path):
        public = re.match(r"public\s*::\s*(.*)", statement, re.I)
        if public:
            exports.extend(
                Export(name, path, module, line)
                for name in split_names(public.group(1))
            )
            continue
        use = re.match(
            r"use(?:\s*,[^ ]+)?\s+([a-z_]\w*)\s*(?:,\s*only\s*:\s*(.*))?$",
            statement,
            re.I,
        )
        if not use:
            continue
        imported_module = use.group(1).lower()
        if use.group(2) is None:
            wildcards.append(imported_module)
            continue
        for item in split_names(use.group(2)):
            alias = re.split(r"\s*=>\s*", item, maxsplit=1)
            local = alias[0].strip().lower()
            original = (alias[1] if len(alias) == 2 else alias[0]).strip()
            imports[local] = Import(imported_module, original)
    return module, exports, imports, wildcards


def source_modules():
    result: dict[str, Path] = {}
    for path in sorted((ROOT / "src").rglob("*.f90")):
        name = module_name(path)
        result.setdefault(name, path.relative_to(ROOT))
    return result


def resolve(module: str, symbol: str, imports, exports, seen=None):
    """Follow API-wrapper re-exports to the first defining module."""
    seen = set() if seen is None else seen
    key = (module, symbol.lower())
    if key in seen:
        return module, symbol, "cycle"
    seen.add(key)
    imported = imports.get(module, {}).get(symbol.lower())
    if imported is None:
        return module, symbol, "direct"
    return resolve(imported.module, imported.name, imports, exports, seen)


def facade(module: str, api_module: str) -> str:
    """Apply the named §1.4 facade split using stable module-name rules."""
    name = f"{module} {api_module}".lower()
    if "api_plot" in name or "plot" in name:
        return "fortfem_plot"
    if any(word in name for word in (
        "interop", "oracle", "equilibrium", "external", "schema",
        "diagnostic", "comparison", "response_interchange",
    )):
        return "fortfem_interop"
    if any(word in name for word in (
        "fourier", "toroid", "harmonic", "legendre", "special_function",
        "mode_registry", "modal",
    )):
        return "fortfem_fourier"
    if any(word in name for word in (
        "boundary", "bem", "dtn", "ntd", "pml", "wall", "surface_current",
        "free_boundary", "trace_map", "open_boundary",
    )):
        return "fortfem_boundary"
    if any(word in name for word in (
        "time", "symplectic", "wave", "transient", "continuation",
        "arclength", "deflated", "step_reduction", "resistive_mhd", "mhd",
    )):
        return "fortfem_time"
    if any(word in name for word in (
        "api_types", "api_mesh", "api_mesh_boundaries", "kinds", "status",
        "topology", "cell_complex", "geometry", "mesh_",
    )):
        return "fortfem_core"
    if any(word in name for word in ("api_forms", "api_spaces", "api_solvers")):
        return "fortfem_feec"
    return "fortfem_feec"


def spelling(symbol: str) -> str:
    match = re.fullmatch(r"evaluate_(.+)_parity", symbol, re.I)
    if match:
        return f"compare_{match.group(1)} (proposed)"
    return "same"


def derivative_actions(symbols: set[str], symbol: str) -> str:
    actions = []
    for suffix, label in (("_jvp", "JVP"), ("_vjp", "VJP"), ("_hvp", "HVP")):
        if symbol.lower().endswith(suffix):
            actions.append(label)
        elif f"{symbol}{suffix}".lower() in symbols:
            actions.append(label)
    return ", ".join(actions) if actions else "—"


def test_examples(symbol: str, index: dict[str, list[str]]) -> str:
    if not re.match(r"^[A-Za-z_]\w*$", symbol):
        return "— (operator/generic; search manually)"
    matches = index.get(symbol.lower(), [])
    if not matches:
        return "—"
    shown = ", ".join(f"`{path}`" for path in matches[:3])
    return shown if len(matches) <= 3 else f"{shown}; +{len(matches) - 3}"


def provenance(module: str, source: Path | None) -> str:
    if source is None:
        return "unresolved"
    absolute = ROOT / source
    text = absolute.read_text()
    if "src/generated/" in str(source):
        generator = re.search(r"! Generator:\s*(.*)", text)
        revision = re.search(r"! Generator revision:\s*(.*)", text)
        details = generator.group(1).strip() if generator else "unknown generator"
        if revision:
            details += f"; {revision.group(1).strip()}"
        return f"FortSym generated ({details})"
    generated = sorted(set(re.findall(r"use\s+(fortfem_generated_\w+)", text, re.I)))
    if generated:
        return "uses generated kernel(s): " + ", ".join(f"`{name}`" for name in generated)
    return "hand-written / no generator link in defining file"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", type=Path, default=ROOT / "doc/api_public_inventory.md")
    args = parser.parse_args()

    parsed = [parse_api(path) for path in API_FILES]
    exports: list[Export] = [export for _, values, _, _ in parsed for export in values]
    imports = {module: values for module, _, values, _ in parsed}
    api_exports: dict[str, list[Export]] = {
        module: values for module, values, _, _ in parsed
    }
    for module, _, module_imports, wildcards in parsed:
        for imported_module in wildcards:
            for exported in api_exports.get(imported_module, []):
                module_imports.setdefault(exported.symbol.lower(), Import(
                    imported_module, exported.symbol))
    exported_names = {export.symbol.lower() for export in exports}
    module_files = source_modules()

    # Keep first occurrence order (umbrella first), while recording all wrappers.
    records: dict[str, dict] = {}
    for export in exports:
        key = export.symbol.lower()
        record = records.setdefault(key, {
            "symbol": export.symbol,
            "wrappers": [],
            "line": export.line,
        })
        wrapper = f"{export.module} ({export.file.relative_to(ROOT)})"
        if wrapper not in record["wrappers"]:
            record["wrappers"].append(wrapper)
    occurrence_index: dict[str, list[str]] = defaultdict(list)
    for directory in (ROOT / "test", ROOT / "example"):
        if directory.exists():
            for path in sorted(directory.rglob("*.f90")):
                text = path.read_text(errors="replace")
                relative = str(path.relative_to(ROOT))
                for token in set(re.findall(r"[A-Za-z_]\w*", text)):
                    occurrence_index[token.lower()].append(relative)

    rows = []
    for record in records.values():
        symbol = record["symbol"]
        wrapper_module = record["wrappers"][0].split(" ", 1)[0]
        defining_module, _, _ = resolve(wrapper_module, symbol, imports, exports)
        source = module_files.get(defining_module.lower())
        rows.append({
            "symbol": symbol,
            "wrappers": "; ".join(f"`{item}`" for item in record["wrappers"]),
            "defining": f"`{defining_module}` / `{source}`" if source else f"`{defining_module}` / unresolved",
            "facade": facade(defining_module, wrapper_module),
            "spelling": spelling(symbol),
            "derivatives": derivative_actions(exported_names, symbol.lower()),
            "tests": test_examples(symbol, occurrence_index),
            "provenance": provenance(defining_module, source),
        })

    lines = [
        "# Public Fortran API inventory",
        "",
        "This deterministic API-00 inventory is derived from the `public ::` and",
        "`use ... only:` declarations in `src/fortfem_api.f90`, existing",
        "`src/fortfem_api_*.f90` wrappers, and canonical facades. It records the current exported spelling",
        "before any facade migration; it does not rename symbols or infer numerical",
        "semantics.",
        "",
        "Regenerate from the repository root with:",
        "",
        "```text",
        "python3 tools/generate_api_public_inventory.py",
        "```",
        "",
        "## Scope and interpretation",
        "",
        f"The inventory contains **{len(rows)} unique exported symbols** ({len(exports)} declarations, including wrapper duplicates).",
        "The `Defining module/file` column follows wrapper re-exports through `use`",
        "associations where they are explicit. `unresolved` means that the source",
        "does not expose an explicit `only:` mapping; it is deliberately not guessed.",
        "",
        "The `Proposed facade` values are the §1.4 layer names. They are a migration",
        "classification, not new modules: `fortfem_core`, `fortfem_feec`,",
        "`fortfem_fourier`, `fortfem_boundary`, `fortfem_time`,",
        "`fortfem_interop`, and `fortfem_plot`.",
        "",
        "Derivative actions are inferred only from exact exported `_jvp`, `_vjp`,",
        "and `_hvp` companions. Test/example paths are exact textual references in",
        "`test/**/*.f90` and `example/**/*.f90`; `—` is an unresolved mapping, not a",
        "claim that the behavior is untested. Generator provenance is reported from",
        "generated headers or explicit `fortfem_generated_*` imports.",
        "",
        "## Inventory",
        "",
        "| Symbol | Exported by wrapper(s) | Defining module/file | Proposed facade | Current → proposed spelling | Derivative actions | Tests/examples | Generator provenance |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for row in rows:
        lines.append("| " + " | ".join(row[key] for key in (
            "symbol", "wrappers", "defining", "facade", "spelling", "derivatives",
            "tests", "provenance",
        )) + " |")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
