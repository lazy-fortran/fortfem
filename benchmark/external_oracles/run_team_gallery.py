#!/usr/bin/env python3
"""Render verified optional TEAM payloads as solution-first SVG/CSV galleries.

No TEAM arrays are shipped in FortFEM.  This standard-library-only renderer
consumes a JSON payload from the reviewed sister repository only after the
adapter has checked its provenance, schema, archive URI, and SHA-256 digest.
"""

from __future__ import annotations

import argparse
import html
import json
import math
import sys
from pathlib import Path
from typing import Any

from run_team_adapter import (
    CONTRACT_CASES,
    PAYLOAD_SCHEMA,
    ContractError,
    canonical_file_sha256,
    case_map,
    load_json,
    locate_verified_artifact,
    require,
    validate_manifest,
)


def finite_number(value: Any, label: str) -> float:
    require(isinstance(value, (int, float)) and not isinstance(value, bool),
            f"{label} must be numeric")
    number = float(value)
    require(math.isfinite(number), f"{label} must be finite")
    return number


def scalar_solution(value: Any, label: str) -> float:
    if isinstance(value, list):
        require(value, f"{label} vector must not be empty")
        components = [finite_number(entry, f"{label}[{index}]")
                      for index, entry in enumerate(value)]
        return math.sqrt(sum(component * component for component in components))
    return finite_number(value, label)


def require_drawable_element(nodes: list[list[float]], indices: list[int],
                             label: str) -> None:
    """Reject zero-area elements before they become misleading solution plots."""

    require(len(set(indices)) == len(indices),
            f"{label} must not repeat a node")
    origin = nodes[indices[0]]
    scale = max(1.0, *(abs(value) for index in indices for value in nodes[index]))
    largest_cross_squared = 0.0
    for left, right in zip(indices[1:-1], indices[2:]):
        first = [nodes[left][axis] - origin[axis] for axis in range(3)]
        second = [nodes[right][axis] - origin[axis] for axis in range(3)]
        cross = [first[1] * second[2] - first[2] * second[1],
                 first[2] * second[0] - first[0] * second[2],
                 first[0] * second[1] - first[1] * second[0]]
        largest_cross_squared = max(
            largest_cross_squared, sum(component * component for component in cross))
    require(largest_cross_squared > 1.0e-24 * scale**4,
            f"{label} must have non-zero area")


def load_case_payload(artifact: Path, case_id: str, contract: dict[str, Any]) -> tuple[
        list[list[float]], list[list[int]], list[float], str]:
    document = load_json(artifact)
    require(document.get("schema") == PAYLOAD_SCHEMA,
            f"payload schema must be {PAYLOAD_SCHEMA}")
    payload_cases = document.get("cases")
    require(isinstance(payload_cases, dict), "payload.cases must be an object")
    raw_case = payload_cases.get(case_id)
    require(isinstance(raw_case, dict), f"payload has no {case_id} case")
    links = case_map(contract)[case_id]["provenance_links"]
    provenance_uri = raw_case.get("provenance_uri")
    require(provenance_uri in links,
            f"{case_id} payload provenance URI is not a contract link")

    geometry = raw_case.get("geometry")
    require(isinstance(geometry, dict), f"{case_id} geometry is required")
    raw_nodes = geometry.get("nodes")
    raw_elements = geometry.get("elements")
    require(isinstance(raw_nodes, list) and len(raw_nodes) >= 3,
            f"{case_id} geometry.nodes must contain at least three nodes")
    require(isinstance(raw_elements, list) and raw_elements,
            f"{case_id} geometry.elements must not be empty")
    nodes: list[list[float]] = []
    for node_index, node in enumerate(raw_nodes):
        require(isinstance(node, list) and len(node) in (2, 3),
                f"{case_id} geometry.nodes[{node_index}] must have 2 or 3 coordinates")
        coordinates = [finite_number(value, f"{case_id} node[{node_index}]")
                       for value in node]
        if len(coordinates) == 2:
            coordinates.append(0.0)
        nodes.append(coordinates)

    raw_indices: list[list[int]] = []
    for element_index, element in enumerate(raw_elements):
        require(isinstance(element, list) and len(element) >= 3,
                f"{case_id} elements[{element_index}] must have at least 3 nodes")
        indices: list[int] = []
        for local_index, value in enumerate(element):
            require(isinstance(value, int) and not isinstance(value, bool),
                    f"{case_id} elements[{element_index}][{local_index}] must be integer")
            indices.append(value)
        raw_indices.append(indices)
    all_indices = [index for element in raw_indices for index in element]
    offset = 0 if any(index == 0 for index in all_indices) else 1
    elements = [[index - offset for index in element] for element in raw_indices]
    for element_index, element in enumerate(elements):
        require(all(0 <= index < len(nodes) for index in element),
                f"{case_id} elements[{element_index}] references an invalid node")
        require_drawable_element(nodes, element,
                                 f"{case_id} elements[{element_index}]")

    raw_solution = raw_case.get("solution")
    if isinstance(raw_solution, dict):
        raw_solution = raw_solution.get("values")
    require(isinstance(raw_solution, list) and len(raw_solution) == len(nodes),
            f"{case_id} solution.values must match node count")
    solution = [scalar_solution(value, f"{case_id} solution[{index}]")
                for index, value in enumerate(raw_solution)]
    return nodes, elements, solution, provenance_uri


def color(value: float, lower: float, upper: float) -> str:
    span = max(upper - lower, 1.0e-15)
    t = min(1.0, max(0.0, (value - lower) / span))
    red = int(round(255.0 * t))
    green = int(round(100.0 + 120.0 * (1.0 - abs(2.0 * t - 1.0))))
    blue = int(round(255.0 * (1.0 - t)))
    return f"#{red:02x}{green:02x}{blue:02x}"


def write_gallery(output_dir: Path, case_id: str, nodes: list[list[float]],
                  elements: list[list[int]], solution: list[float],
                  provenance_uri: str, digest: str, artifact: Path) -> None:
    projected = [(node[0] - 0.45 * node[1], node[2] + 0.30 * node[1])
                 for node in nodes]
    min_u = min(point[0] for point in projected)
    max_u = max(point[0] for point in projected)
    min_v = min(point[1] for point in projected)
    max_v = max(point[1] for point in projected)
    span = max(max_u - min_u, max_v - min_v, 1.0e-12)

    def screen(point: tuple[float, float]) -> tuple[float, float]:
        return (90.0 + 760.0 * (point[0] - min_u) / span,
                610.0 - 520.0 * (point[1] - min_v) / span)

    lower = min(solution)
    upper = max(solution)
    polygons: list[str] = []
    for element in elements:
        vertices = [screen(projected[index]) for index in element]
        average = sum(solution[index] for index in element) / len(element)
        points = " ".join(f"{x:.3f},{y:.3f}" for x, y in vertices)
        polygons.append(
            f'<polygon points="{points}" fill="{color(average, lower, upper)}" '
            'stroke="#25324a" stroke-width="1.2" fill-opacity="0.88"/>')
    circles: list[str] = []
    for index, point in enumerate(projected):
        x, y = screen(point)
        circles.append(
            f'<circle cx="{x:.3f}" cy="{y:.3f}" r="3.0" '
            f'fill="{color(solution[index], lower, upper)}" stroke="#101820" '
            'stroke-width="0.7"/>')
    title = f"{case_id} exact-data solution (verified)"
    subtitle = "solution-first • external TEAM payload • no solver bundled"
    svg = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="960" height="700" '
        'viewBox="0 0 960 700">',
        '<rect width="960" height="700" fill="#ffffff"/>',
        '<text x="60" y="38" font-family="sans-serif" font-size="22" '
        'font-weight="600" fill="#17233b">' + html.escape(title) + '</text>',
        '<text x="60" y="64" font-family="sans-serif" font-size="13" '
        'fill="#394865">' + html.escape(subtitle) + '</text>',
        *polygons,
        *circles,
        '<line x1="885" y1="130" x2="885" y2="610" stroke="#202b3d" '
        'stroke-width="1"/>',
    ]
    for step in range(11):
        value = lower + (upper - lower) * step / 10.0
        y = 610.0 - 480.0 * step / 10.0
        svg.append(f'<line x1="880" y1="{y:.1f}" x2="890" y2="{y:.1f}" '
                   'stroke="#202b3d" stroke-width="1"/>')
        svg.append(f'<text x="898" y="{y + 4.0:.1f}" font-family="sans-serif" '
                   f'font-size="11" fill="#202b3d">{value:.5g}</text>')
    svg.extend([
        '<text x="850" y="112" font-family="sans-serif" font-size="12" '
        'fill="#202b3d">solution</text>',
        '<text x="60" y="658" font-family="sans-serif" font-size="11" '
        'fill="#52617a">provenance: ' + html.escape(provenance_uri) +
        ' • SHA-256: ' + digest + '</text>',
        '<text x="60" y="678" font-family="sans-serif" font-size="11" '
        'fill="#52617a">artifact: ' + html.escape(artifact.name) +
        ' • exact TEAM arrays remain external</text>',
        '</svg>',
    ])
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "solution.svg").write_text("\n".join(svg) + "\n", encoding="utf-8")
    with (output_dir / "solution.csv").open("w", encoding="utf-8") as stream:
        stream.write("node,x,y,z,solution\n")
        for index, node in enumerate(nodes):
            stream.write(f"{index},{node[0]:.16e},{node[1]:.16e},"
                         f"{node[2]:.16e},{solution[index]:.16e}\n")
    metadata = {
        "status": "ready",
        "case_id": case_id,
        "provenance_uri": provenance_uri,
        "exact_data": True,
        "manufactured": False,
        "comparison_mode": "external-benchmark-payload",
        "analytical_reference": "provenance-only; no analytical arrays bundled",
        "solution_first": True,
        "plot": "solution.svg",
        "artifact_sha256": digest,
        "node_count": len(nodes),
        "element_count": len(elements),
    }
    (output_dir / "provenance.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8")


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--contract", type=Path,
                        default=Path(__file__).with_name("team_manifest.json"))
    result.add_argument("--data-manifest", type=Path)
    result.add_argument("--data-root", type=Path)
    result.add_argument("--case-id", choices=CONTRACT_CASES, action="append",
                        help="render only this case (repeat for several cases)")
    result.add_argument("--output-dir", type=Path,
                        default=Path("output/external_oracle/team_gallery"))
    return result


def main() -> int:
    args = parser().parse_args()
    try:
        contract = load_json(args.contract)
        validate_manifest(contract)
        artifact, _ = locate_verified_artifact(contract, args.data_manifest, args.data_root)
        if artifact is None:
            print("SKIP: exact TEAM payloads are absent; no solution plots generated")
            return 0
        selected = args.case_id or list(CONTRACT_CASES)
        digest = canonical_file_sha256(artifact)
        for case_id in selected:
            nodes, elements, solution, provenance_uri = load_case_payload(
                artifact, case_id, contract)
            write_gallery(args.output_dir / case_id, case_id, nodes, elements,
                          solution, provenance_uri, digest, artifact)
        print(f"SOLUTION_FIRST: wrote verified TEAM gallery to {args.output_dir}")
        print(f"               cases={','.join(selected)}; sha256={digest}")
        return 0
    except (OSError, json.JSONDecodeError, ContractError) as error:
        print(f"TEAM solution gallery: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
