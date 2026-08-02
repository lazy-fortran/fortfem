#!/usr/bin/env python3
"""Audit duplicate exports at the canonical-facade boundary.

The compatibility umbrella may re-export every public symbol.  Canonical
facades, however, must have one owner unless an intentional cross-facade
re-export is recorded in ``doc/api_compatibility_allowlist.txt``.  This small
standard-library-only check catches accidental duplicate wrappers and stale
allowlist entries without compiling or importing application physics.
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path


FACADE_NAMES = (
    "fortfem_core",
    "fortfem_feec",
    "fortfem_fourier",
    "fortfem_boundary",
    "fortfem_time",
    "fortfem_interop",
    "fortfem_plot",
)
PUBLIC_RE = re.compile(r"^\s*public\s*::\s*(.*)$", re.IGNORECASE)


def logical_lines(path: Path):
    """Yield free-form statements while joining ampersand continuations."""

    pending = ""
    start = 0
    for number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw.split("!", 1)[0].rstrip()
        if not pending:
            pending = line
            start = number
        else:
            pending += " " + line.lstrip("& ")
        if pending.rstrip().endswith("&"):
            pending = pending.rstrip()[:-1].rstrip()
        elif pending.strip():
            yield start, pending.strip()
            pending = ""


def facade_exports(root: Path) -> dict[str, list[str]]:
    exports: dict[str, list[str]] = {}
    for facade in FACADE_NAMES:
        path = root / "src" / f"{facade}.f90"
        if not path.is_file():
            continue
        names: list[str] = []
        for _, statement in logical_lines(path):
            match = PUBLIC_RE.match(statement)
            if match:
                names.extend(
                    name.strip().lower()
                    for name in match.group(1).split(",")
                    if name.strip()
                )
        exports[facade] = names
    return exports


def allowlist(root: Path) -> dict[str, tuple[str, str]]:
    path = root / "doc" / "api_compatibility_allowlist.txt"
    if not path.is_file():
        return {}
    result: dict[str, tuple[str, str]] = {}
    for number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split("|", 3)
        if len(fields) != 4 or any(not field.strip() for field in fields[:3]):
            raise ValueError(f"{path}:{number}: expected symbol|facade|facade|reason")
        symbol, first, second, _ = (field.strip().lower() for field in fields)
        if first == second:
            raise ValueError(f"{path}:{number}: duplicate facades must differ")
        if symbol in result:
            raise ValueError(f"{path}:{number}: duplicate allowlist symbol {symbol}")
        result[symbol] = tuple(sorted((first, second)))
    return result


def check(root: Path) -> list[str]:
    diagnostics: list[str] = []
    try:
        permitted = allowlist(root)
    except (OSError, ValueError) as error:
        return [str(error)]
    exports = facade_exports(root)
    providers: dict[str, list[str]] = defaultdict(list)
    for facade, names in exports.items():
        seen: set[str] = set()
        for symbol in names:
            if symbol in seen:
                diagnostics.append(
                    f"{facade}: duplicate public declaration {symbol}"
                )
            seen.add(symbol)
            providers[symbol].append(facade)

    actual_duplicates: dict[str, tuple[str, ...]] = {}
    for symbol, facades in providers.items():
        unique = tuple(sorted(set(facades)))
        if len(unique) <= 1:
            continue
        actual_duplicates[symbol] = unique
        expected = permitted.get(symbol)
        if expected != unique:
            diagnostics.append(
                "duplicate canonical export "
                f"{symbol}: {', '.join(unique)}; add the exact intentional "
                "pair to doc/api_compatibility_allowlist.txt or remove a wrapper"
            )

    for symbol, facades in permitted.items():
        actual = actual_duplicates.get(symbol)
        if actual != facades:
            diagnostics.append(
                f"stale compatibility allowlist entry {symbol}: "
                f"declared {', '.join(facades)}, found "
                f"{', '.join(actual) if actual else 'no duplicate export'}"
            )
        for facade in facades:
            if facade not in exports:
                diagnostics.append(
                    f"allowlist entry {symbol}: missing facade {facade}"
                )
    return sorted(set(diagnostics))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    root = args.root.resolve()
    errors = check(root)
    if errors:
        for error in errors:
            print(f"ERROR: {error}")
        print(f"API compatibility check failed: {len(errors)} diagnostic(s)")
        return 1
    # Count documented pairs from the actual source, not the allowlist, so the
    # output is deterministic and useful in release logs.
    exports = facade_exports(root)
    names = defaultdict(set)
    for facade, symbols in exports.items():
        for symbol in symbols:
            names[symbol].add(facade)
    duplicate_count = sum(len(facades) > 1 for facades in names.values())
    print(
        "API compatibility check passed: "
        f"{len(exports)} canonical facades, {duplicate_count} documented duplicate(s)"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
