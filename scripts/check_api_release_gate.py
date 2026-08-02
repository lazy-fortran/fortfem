#!/usr/bin/env python3
"""Check the fast, pre-release API migration invariants.

The release gate intentionally keeps the expensive jobs (full example gallery,
cross-compiler builds and external-code comparisons) out of the ten-second
path.  This checker covers the two deterministic checks which can run without
building Fortran: the checked-in public inventory must be reproducible and
obsolete spellings must not remain in product/client sources.

The module-layer and generated-visibility checks are separate focused scripts;
``test/test_api_release_gate.sh`` composes them with this checker and the
small derivative tests.  ``--stale-only`` is provided for an independent
negative fixture, so the stale-name failure cannot be masked by another gate.
"""

from __future__ import annotations

import argparse
import difflib
import re
import subprocess
import sys
import tempfile
from pathlib import Path


SOURCE_DIRECTORIES = ("src", "example", "app", "benchmark", "test")
SOURCE_SUFFIXES = {".f", ".f90", ".f95", ".f03", ".f08", ".py", ".sh"}

# These names were removed by the API-03 boundary/larger-domain migrations.
# Keep the list explicit until each remaining candidate is migrated; the
# candidate inventory and design notes intentionally continue to mention
# future renames and therefore are not part of this scan.
OBSOLETE_NAMES = (
    "evaluate_boundary_operator_parity",
    "evaluate_boundary_operator_parity_jvp",
    "evaluate_boundary_operator_parity_vjp",
    "evaluate_larger_domain_parity",
    "evaluate_larger_domain_parity_jvp",
    "evaluate_sheet_current_parity",
    "evaluate_sheet_current_surface_parity",
    "evaluate_sheet_current_surface_parity_jvp",
    "evaluate_tree_cotree_iga_parity",
    "evaluate_beltrami_two_region_parity",
    "evaluate_beltrami_shell_parity",
)


def source_files(root: Path):
    """Yield product and downstream source files in deterministic order."""

    for directory in SOURCE_DIRECTORIES:
        base = root / directory
        if not base.is_dir():
            continue
        for path in sorted(base.rglob("*")):
            if path.is_file() and path.suffix.lower() in SOURCE_SUFFIXES:
                # The generated-visibility checker owns its fixture trees;
                # they are deliberately allowed to contain forbidden names.
                if path.relative_to(root).parts[:2] == ("test", "fixtures"):
                    continue
                yield path


def stale_name_diagnostics(root: Path) -> list[str]:
    """Return diagnostics for removed names still used by source/client code."""

    diagnostics: list[str] = []
    patterns = {
        name: re.compile(rf"(?<![A-Za-z0-9_]){re.escape(name)}(?![A-Za-z0-9_])")
        for name in OBSOLETE_NAMES
    }
    for path in source_files(root):
        try:
            lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        except OSError as error:
            diagnostics.append(f"{path}: cannot read source: {error}")
            continue
        for line_number, line in enumerate(lines, start=1):
            for name, pattern in patterns.items():
                if pattern.search(line):
                    diagnostics.append(
                        f"{path.relative_to(root)}:{line_number}: stale internal "
                        f"API spelling {name}; use its canonical compare_* name"
                    )
    return diagnostics


def inventory_diagnostics(root: Path, inventory: Path | None) -> list[str]:
    """Regenerate the public inventory and compare it byte-for-byte."""

    expected = inventory or root / "doc" / "api_public_inventory.md"
    if not expected.is_file():
        return [f"public API inventory is missing: {expected}"]

    generator = root / "tools" / "generate_api_public_inventory.py"
    if not generator.is_file():
        return [f"public API inventory generator is missing: {generator}"]

    with tempfile.TemporaryDirectory(prefix="fortfem-api-inventory-") as directory:
        generated = Path(directory) / "api_public_inventory.md"
        result = subprocess.run(
            [sys.executable, str(generator), "--output", str(generated)],
            cwd=root,
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode != 0:
            detail = (result.stdout + result.stderr).strip()
            return [
                "public API inventory regeneration failed"
                + (f": {detail}" if detail else "")
            ]
        actual = generated.read_text(encoding="utf-8")
        expected_text = expected.read_text(encoding="utf-8")
        if actual == expected_text:
            return []

        diff = list(
            difflib.unified_diff(
                expected_text.splitlines(),
                actual.splitlines(),
                fromfile=str(expected),
                tofile="regenerated inventory",
                n=1,
            )
        )
        # A stale inventory can be very large.  Keep the fast-path diagnostic
        # actionable without flooding CI logs; the complete diff is produced by
        # rerunning the generator locally.
        preview = "\n".join(diff[:24])
        if len(diff) > 24:
            preview += "\n... (diff truncated; regenerate locally for the full diff)"
        return ["public API inventory is stale:\n" + preview]


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Check deterministic API migration release invariants."
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("."),
        help="repository (or fixture) root to scan",
    )
    parser.add_argument(
        "--inventory",
        type=Path,
        help="inventory file to compare (default: ROOT/doc/api_public_inventory.md)",
    )
    parser.add_argument(
        "--stale-only",
        action="store_true",
        help="only scan for removed names (used by the negative fixture)",
    )
    parser.add_argument(
        "--inventory-only",
        action="store_true",
        help="only check inventory regeneration (used by the negative fixture)",
    )
    return parser.parse_args(sys.argv[1:] if argv is None else argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.stale_only and args.inventory_only:
        print("--stale-only and --inventory-only are mutually exclusive", file=sys.stderr)
        return 2
    root = args.root.resolve()
    if not root.is_dir():
        print(f"error: API release-gate root does not exist: {root}", file=sys.stderr)
        return 2

    diagnostics: list[str] = []
    if not args.stale_only:
        inventory = args.inventory.resolve() if args.inventory else None
        diagnostics.extend(inventory_diagnostics(root, inventory))
    if not args.inventory_only:
        diagnostics.extend(stale_name_diagnostics(root))

    if diagnostics:
        for diagnostic in diagnostics:
            print(f"ERROR: {diagnostic}")
        print(f"API release-gate check failed: {len(diagnostics)} diagnostic(s)")
        return 1
    scope = "stale-name" if args.stale_only else "inventory and stale-name"
    if args.inventory_only:
        scope = "inventory"
    print(f"API release-gate {scope} checks passed: {root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
