#!/usr/bin/env python3
"""Check that FortSym kernels stay behind stable FortFEM wrappers.

The generated Fortran in ``src/generated`` is an implementation detail.  A
public facade or downstream client must import a wrapper module instead of a
generated module directly.  Implementation wrappers may import generated
kernels, but they must have an explicitly private default and an explicit
public interface.  This keeps generated names out of the public namespace
while allowing the generated code to remain a fast lower-level dependency.

This checker deliberately uses only the Python standard library.  It reads
compiler-visible ``module`` and ``use`` statements; it does not rely on file
names alone or on a build directory containing ``.mod`` files.
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path


FORTRAN_SUFFIXES = {".f", ".f90", ".f95", ".f03", ".f08"}
GENERATED_PREFIX = "fortfem_generated_"
CLIENT_ROOTS = {"example", "app", "benchmark"}
TEST_ROOT = "test"
# This source is an analytical manufactured-solution oracle, not a product
# kernel.  Keep the exception explicit until the example is migrated to a
# generated-oracle facade; all other generated imports from clients fail.
ALLOWED_ORACLE_CLIENT_MODULES = {"fortfem_generated_tetra_h1_oracle"}

MODULE_RE = re.compile(
    r"^\s*module\s+(?!procedure\b)(?!subroutine\b)"
    r"(?!function\b)([a-z_]\w*)\b",
    re.IGNORECASE,
)
USE_RE = re.compile(
    r"^\s*use\b(?:\s*,\s*[a-z_]\w*)*\s*(?:::)?\s*"
    r"([a-z_]\w*)\b",
    re.IGNORECASE,
)
PUBLIC_RE = re.compile(r"^\s*public\s*::\s*(.*)$", re.IGNORECASE)
PRIVATE_RE = re.compile(r"^\s*private(?:\s*::|\s*$)", re.IGNORECASE)


@dataclass(frozen=True)
class UseRecord:
    name: str
    line: int


@dataclass(frozen=True)
class SourceRecord:
    path: Path
    module: str
    kind: str
    uses: tuple[UseRecord, ...]
    has_private_default: bool
    public_names: tuple[str, ...]


def strip_comment(line: str) -> str:
    """Remove a Fortran comment without treating ``!`` in a string as one."""

    quote: str | None = None
    index = 0
    while index < len(line):
        char = line[index]
        if quote is None and char in "'\"":
            quote = char
        elif quote is not None and char == quote:
            if index + 1 < len(line) and line[index + 1] == quote:
                index += 1
            else:
                quote = None
        elif quote is None and char == "!":
            return line[:index]
        index += 1
    return line


def source_kind(path: Path, root: Path, module: str) -> str:
    relative = path.relative_to(root)
    parts = relative.parts
    if "generated" in {part.lower() for part in parts[:-1]}:
        return "generated"
    if parts and parts[0].lower() == "src":
        # Root-level fortfem modules are public facades in the staged API
        # migration.  Nested source modules are implementation wrappers.
        if len(parts) == 2 and module.startswith("fortfem_"):
            return "facade"
        if len(parts) >= 2 and parts[1].lower() == "bindings":
            return "client"
        return "implementation"
    if parts and parts[0].lower() in CLIENT_ROOTS:
        return "client"
    if parts and parts[0].lower() == TEST_ROOT:
        return "test"
    return "implementation"


def iter_sources(root: Path) -> list[Path]:
    """Return source files in the public/implementation trees.

    Fixture directories are skipped when checking the repository itself.  A
    fixture is checked by passing its own directory as ``--root``.
    """

    paths: list[Path] = []
    for directory in ("src", "example", "app", "benchmark", "test"):
        candidate = root / directory
        if not candidate.is_dir():
            continue
        for path in candidate.rglob("*"):
            if not path.is_file() or path.suffix.lower() not in FORTRAN_SUFFIXES:
                continue
            relative = path.relative_to(root)
            if relative.parts[:2] == ("test", "fixtures"):
                continue
            paths.append(path)
    return sorted(paths)


def discover(root: Path) -> tuple[list[SourceRecord], list[str]]:
    records: list[SourceRecord] = []
    diagnostics: list[str] = []
    for path in iter_sources(root):
        try:
            lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        except OSError as error:
            diagnostics.append(f"{path}: cannot read source: {error}")
            continue
        modules: list[str] = []
        uses: list[UseRecord] = []
        public_names: list[str] = []
        has_private_default = False
        for line_number, raw_line in enumerate(lines, start=1):
            line = strip_comment(raw_line)
            module_match = MODULE_RE.match(line)
            if module_match:
                modules.append(module_match.group(1).lower())
            use_match = USE_RE.match(line)
            if use_match:
                uses.append(UseRecord(use_match.group(1).lower(), line_number))
            if PRIVATE_RE.match(line):
                has_private_default = True
            public_match = PUBLIC_RE.match(line)
            if public_match:
                for name in public_match.group(1).split(","):
                    name = name.strip().rstrip("&").strip().lower()
                    if name:
                        public_names.append(name)
        if not modules and not uses:
            continue
        # Programs and other downstream clients have USE statements but no
        # MODULE declaration.  Give those files a private pseudo-provider so
        # their imports are still checked without changing module discovery.
        pseudo_module = "__file__" + "_".join(
            part.lower().replace(".", "_")
            for part in path.relative_to(root).parts
        )
        discovered_modules = modules or [pseudo_module]
        kind = source_kind(path, root, discovered_modules[0])
        for module in discovered_modules:
            records.append(
                SourceRecord(
                    path,
                    module,
                    kind,
                    tuple(uses),
                    has_private_default,
                    tuple(public_names),
                )
            )
    return records, diagnostics


def display_path(path: Path, root: Path) -> str:
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


def check(root: Path) -> list[str]:
    records, diagnostics = discover(root)
    generated_modules = {
        record.module for record in records if record.kind == "generated"
    }
    providers: dict[str, SourceRecord] = {}
    for record in records:
        previous = providers.get(record.module)
        if previous is not None:
            diagnostics.append(
                "duplicate module provider: "
                f"{record.module} in {display_path(previous.path, root)} and "
                f"{display_path(record.path, root)}"
            )
        else:
            providers[record.module] = record

    for record in records:
        for use in record.uses:
            target = use.name
            is_generated = target in generated_modules or target.startswith(
                GENERATED_PREFIX
            )
            if (
                not is_generated
                or record.kind in {"generated", "test"}
                or target in ALLOWED_ORACLE_CLIENT_MODULES
            ):
                continue
            path = display_path(record.path, root)
            if record.kind == "facade":
                diagnostics.append(
                    f"{path}:{use.line}: public facade imports generated module "
                    f"{target}; import a stable wrapper instead"
                )
            elif record.kind == "client":
                diagnostics.append(
                    f"{path}:{use.line}: public client imports generated module "
                    f"{target}; import a stable wrapper instead"
                )
            else:
                if not record.has_private_default:
                    diagnostics.append(
                        f"{path}: generated wrapper {record.module} must declare "
                        "PRIVATE before exposing generated kernels"
                    )
                if not record.public_names:
                    diagnostics.append(
                        f"{path}: generated wrapper {record.module} must expose "
                        "an explicit PUBLIC stable interface"
                    )
                leaked = [
                    name
                    for name in record.public_names
                    if name.startswith("generated_")
                ]
                if leaked:
                    diagnostics.append(
                        f"{path}: generated wrapper {record.module} leaks generated "
                        "names through PUBLIC :: "
                        + ", ".join(sorted(leaked))
                    )
    return sorted(set(diagnostics))


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Audit generated-kernel visibility and stable wrappers."
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("."),
        help="repository or fixture root to scan (default: current directory)",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    root = args.root.resolve()
    if not root.is_dir():
        print(f"error: repository root does not exist: {root}", file=sys.stderr)
        return 2
    diagnostics = check(root)
    if diagnostics:
        for diagnostic in diagnostics:
            print(f"ERROR: {diagnostic}")
        print(f"generated visibility audit failed: {len(diagnostics)} diagnostic(s)")
        return 1
    print(f"generated visibility audit passed: {root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
