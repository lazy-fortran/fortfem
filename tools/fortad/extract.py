#!/usr/bin/env python3
"""Extract fortfem's generated primals and record what each product treats as active.

Every kernel under src/generated that ships a primal alongside a `_jvp` and a
`_vjp` is a candidate for the fortad path. The primal is copied out as a
standalone subroutine fortad can read, and the active arguments are read off the
committed `_jvp` signature rather than guessed: an argument is active exactly
when the generator emitted a `_dot` for it.

Writing the list down here rather than in the shell script keeps the shell to
running the generator, and keeps this parse in one place where it can be read.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

GENERATED = Path("src/generated")
KERNELS = Path("tools/fortad/kernels")

# Directives are the consumer's choice, not the primal's meaning, and fortad
# has no reason to see them.
DIRECTIVE = re.compile(r"^\s*!\$(omp|acc)\b")


def flattened(text: str) -> str:
    """Join Fortran continuations so a signature is one line."""
    return re.sub(r"&\s*\n\s*", "", text)


def signature(text: str, name: str) -> list[str] | None:
    """The dummy argument list of a subroutine, in order."""
    match = re.search(r"pure subroutine %s\(([^)]*)\)" % re.escape(name), text)
    if match is None:
        return None
    return [a.strip() for a in match.group(1).split(",") if a.strip()]


def body(text: str, name: str) -> str | None:
    """The whole subroutine, directives removed."""
    start = re.search(r"^ *pure subroutine %s\(" % re.escape(name), text, re.M)
    if start is None:
        return None
    end = re.search(r"^ *end subroutine %s *$" % re.escape(name), text, re.M)
    if end is None:
        return None
    lines = text[start.start():end.end()].split("\n")
    kept = [ln for ln in lines if not DIRECTIVE.match(ln)]
    # Left-align: the source sits inside a module.
    indent = min((len(ln) - len(ln.lstrip()) for ln in kept if ln.strip()),
                 default=0)
    return "\n".join(ln[indent:] if ln.strip() else "" for ln in kept)


def dimensions(source: str, names: list[str]) -> dict[str, str]:
    """The declared shape of each dummy, empty for a scalar.

    An argument's shape is part of the call, not decoration: an output declared
    `point(3)` needs three slots on both sides of the comparison, and packing it
    into one would compile on some operators and not others.
    """
    shapes = {name: "" for name in names}
    for line in flattened(source).split("\n"):
        stripped = line.split("!", 1)[0].strip()
        if "::" not in stripped or not stripped.startswith("real"):
            continue
        for item in stripped.split("::", 1)[1].split(","):
            item = item.strip()
            if "(" not in item:
                continue
            base, _, rest = item.partition("(")
            base = base.strip()
            if base in shapes:
                shapes[base] = rest.rstrip(")").strip()
    return shapes


def declared_out(source: str) -> set[str]:
    """Names the primal declares `intent(out)` or `intent(inout)`."""
    written: set[str] = set()
    for line in flattened(source).split("\n"):
        stripped = line.split("!", 1)[0].strip()
        if "intent(out)" not in stripped and "intent(inout)" not in stripped:
            continue
        if "::" not in stripped:
            continue
        for name in stripped.split("::", 1)[1].split(","):
            name = name.strip()
            if name:
                written.add(name.split("(", 1)[0].strip())
    return written


def main() -> int:
    KERNELS.mkdir(parents=True, exist_ok=True)
    catalogue = []

    for path in sorted(GENERATED.glob("*.f90")):
        text = flattened(path.read_text(encoding="utf-8"))
        names = set(re.findall(r"pure subroutine (\w+)\(", text))
        for name in sorted(names):
            if name.endswith("_jvp") or name.endswith("_vjp"):
                continue
            if f"{name}_jvp" not in names or f"{name}_vjp" not in names:
                continue
            primal_args = signature(text, name)
            jvp_args = signature(text, f"{name}_jvp")
            if primal_args is None or jvp_args is None:
                continue
            source = body(path.read_text(encoding="utf-8"), name)
            if source is None:
                continue
            # Both an active input and a differentiated output carry a `_dot`
            # in the jvp signature, so the two cannot be told apart by name.
            # The intents in the primal's own declarations can.
            written = declared_out(source)
            active = [a for a in primal_args
                      if f"{a}_dot" in jvp_args and a not in written]
            outputs = [a for a in primal_args
                       if f"{a}_dot" in jvp_args and a in written]
            if not active or not outputs:
                continue
            stem = name[len("generated_"):] if name.startswith("generated_") else name
            header = (
                f"! Extracted from src/generated/{path.name} by\n"
                f"! tools/fortad/extract.py. It is the primal fortsym generated\n"
                f"! its products from, with the offload directives removed:\n"
                f"! they are the consumer's choice, not part of what the\n"
                f"! function means.\n\n"
            )
            (KERNELS / f"{stem}.f90").write_text(header + source + "\n",
                                                encoding="utf-8")
            catalogue.append({
                "stem": stem,
                "primal": name,
                "active": active,
                "outputs": outputs,
                "dims": dimensions(source, primal_args),
                "file": path.name,
            })

    Path("tools/fortad/catalogue.json").write_text(
        json.dumps(catalogue, indent=2) + "\n", encoding="utf-8")
    print(f"extracted {len(catalogue)} primals into {KERNELS}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
