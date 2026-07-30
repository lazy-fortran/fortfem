#!/usr/bin/env python3
"""Generate lightweight gallery covers when an example has no plot artifact."""

from html import escape
from pathlib import Path
from textwrap import wrap


ROOT = Path(__file__).resolve().parents[1]
ORDER = ROOT / "doc" / "examples" / "gallery_order.txt"
COVERS = ROOT / "artifacts" / "plots" / "examples"

PALETTES = {
    "Getting started": ("#0b3c5d", "#328cc1", "#d9eef7"),
    "Finite elements": ("#284b63", "#3c6e71", "#d9e7e8"),
    "Open boundaries": ("#3d348b", "#7678ed", "#e7e6fa"),
    "Verification and performance": ("#7f4f24", "#f28e2b", "#f9e7d2"),
    "Boundary elements and coupling": ("#5f0f40", "#9a348e", "#f1deee"),
    "Advanced geometry": ("#183a37", "#4b7f52", "#dcebdd"),
}


def motif(group: str, accent: str) -> str:
    if group == "Getting started":
        return f"""
  <path d="M70 280 L250 105 L365 280 Z M70 280 L365 280
           M250 105 L250 280 M70 280 L250 210 L365 280"
        fill="none" stroke="{accent}" stroke-width="7"/>"""
    if group == "Finite elements":
        return f"""
  <path d="M65 285 L180 100 L360 285 Z M180 100 L230 285
           M65 285 L285 170 L360 285" fill="none"
        stroke="{accent}" stroke-width="6"/>
  <path d="M130 245 L170 190 M170 190 L160 215 M170 190 L145 198
           M240 230 L300 165 M300 165 L286 196 M300 165 L270 176"
        fill="none" stroke="#f28e2b" stroke-width="8"/>"""
    if group == "Open boundaries":
        return f"""
  <circle cx="205" cy="205" r="72" fill="none"
      stroke="{accent}" stroke-width="8"/>
  <path d="M110 105 Q205 25 300 105 M75 75 Q205 -35 335 75
           M110 305 Q205 385 300 305 M75 335 Q205 445 335 335"
      fill="none" stroke="{accent}" stroke-width="7"/>"""
    if group == "Verification and performance":
        return f"""
  <path d="M75 300 V190 H125 V300 M160 300 V135 H210 V300
           M245 300 V80 H295 V300 M55 300 H335"
      fill="none" stroke="{accent}" stroke-width="10"/>
  <path d="M90 145 L155 205 L310 55" fill="none"
      stroke="#f28e2b" stroke-width="12"/>"""
    if group == "Boundary elements and coupling":
        return f"""
  <circle cx="205" cy="205" r="120" fill="none"
      stroke="{accent}" stroke-width="7"/>
  <g fill="#f28e2b">
    <circle cx="205" cy="85" r="10"/><circle cx="309" cy="145" r="10"/>
    <circle cx="309" cy="265" r="10"/><circle cx="205" cy="325" r="10"/>
    <circle cx="101" cy="265" r="10"/><circle cx="101" cy="145" r="10"/>
  </g>
  <path d="M205 205 L205 85 M205 205 L309 145 M205 205 L101 265"
      fill="none" stroke="{accent}" stroke-width="5"/>"""
    return f"""
  <ellipse cx="205" cy="205" rx="145" ry="88" fill="none"
      stroke="{accent}" stroke-width="9"/>
  <ellipse cx="205" cy="205" rx="65" ry="37" fill="none"
      stroke="{accent}" stroke-width="9"/>
  <path d="M60 205 Q205 45 350 205 Q205 365 60 205"
      fill="none" stroke="#f28e2b" stroke-width="6"/>"""


def cover_svg(name: str, group: str) -> str:
    dark, accent, pale = PALETTES[group]
    words = name.replace("_", " ")
    lines = wrap(words, width=18, break_long_words=True)
    title_spans = "\n".join(
        f'<tspan x="455" dy="{0 if index == 0 else 43}">{escape(line)}</tspan>'
        for index, line in enumerate(lines)
    )
    return f"""<svg xmlns="http://www.w3.org/2000/svg" width="800" height="450"
 viewBox="0 0 800 450" role="img" aria-labelledby="title description">
<title id="title">{escape(words)} example cover</title>
<desc id="description">FortFEM {escape(group.lower())} example</desc>
<rect width="800" height="450" rx="24" fill="{pale}"/>
<rect x="420" width="380" height="450" fill="{dark}"/>
{motif(group, accent)}
<text x="455" y="90" fill="#ffffff" font-family="sans-serif"
 font-size="24" font-weight="600">{escape(group.upper())}</text>
<text x="455" y="165" fill="#ffffff" font-family="sans-serif"
 font-size="31" font-weight="700">{title_spans}</text>
<text x="455" y="365" fill="#ffffff" font-family="sans-serif"
 font-size="24">FortFEM executable example</text>
</svg>
"""


def main() -> None:
    group = ""
    entries: list[tuple[str, str]] = []
    for raw in ORDER.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if line.startswith("#"):
            group = line.removeprefix("#").strip()
        elif line:
            entries.append((line, group))

    for name, section in entries:
        destination = COVERS / name
        if any(destination.glob("*.png")):
            continue
        destination.mkdir(parents=True, exist_ok=True)
        (destination / "cover.svg").write_text(
            cover_svg(name, section), encoding="utf-8"
        )

    print(f"ensured gallery covers for {len(entries)} examples")


if __name__ == "__main__":
    main()
