#!/usr/bin/env python3
"""Collect a vmec_scan.sbatch output directory into CSV and markdown tables.

Usage: collect.py SCANDIR [--csv out.csv] [--md out.md]
Columns: code, mpol, ns, ftol, threads, DOFs (ns*mnmax*3), ier, converged, iterations,
native solve s (around vmecpp.run only), total s (process wall incl. imports and sampling),
peak RSS MiB (/usr/bin/time -v, whole process), e_B, e_p, e_iota, e_flux, max.
"""
import argparse
import csv
import json
import re
import sys
from pathlib import Path

COLS = ["code", "mpol", "ns", "ftol", "threads", "dofs", "ier", "converged", "iterations",
        "native_s", "total_s", "rss_mib", "e_B", "e_p", "e_iota", "e_flux", "max", "exit"]


def time_rss(path):
    try:
        m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", Path(path).read_text())
        return int(m.group(1)) / 1024 if m else None
    except OSError:
        return None


def row(d):
    rec = {}
    mj = d / "metrics.json"
    m = json.loads(mj.read_text()) if mj.exists() else {}
    l6 = m.get("l6") or {}
    ex = (d / "exit.txt").read_text().strip() if (d / "exit.txt").exists() else ""
    return dict(code="vmecpp-" + str(m.get("version", "?")), mpol=m.get("mpol"),
                ns=m.get("ns"), ftol=m.get("ftol"), threads=m.get("threads"),
                dofs=m.get("coefficient_dofs"), ier=m.get("ier_flag"),
                converged=m.get("converged"), iterations=m.get("iterations"),
                native_s=m.get("native_solve_wall_s"), total_s=m.get("total_wall_s"),
                rss_mib=time_rss(d / "time.txt"), e_B=l6.get("e_B"), e_p=l6.get("e_p"),
                e_iota=l6.get("e_iota"), e_flux=l6.get("e_flux"), max=l6.get("max"),
                exit=ex), rec


def fmt(v):
    if v is None:
        return ""
    if isinstance(v, float):
        return f"{v:.4g}"
    return str(v)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("scandir")
    ap.add_argument("--csv")
    ap.add_argument("--md")
    a = ap.parse_args()
    rows = []
    for d in sorted(Path(a.scandir).iterdir()):
        if d.is_dir() and (d / "metrics.json").exists() or (d.is_dir() and (d / "exit.txt").exists()):
            r, _ = row(d)
            if r["mpol"] is None:
                r["mpol"] = d.name
            rows.append(r)
    rows.sort(key=lambda r: (str(r["ftol"]), r["ns"] or 0, str(r["mpol"]), r["threads"] or 0))
    out = open(a.csv, "w", newline="") if a.csv else sys.stdout
    w = csv.DictWriter(out, fieldnames=COLS)
    w.writeheader()
    w.writerows(rows)
    if a.csv:
        out.close()
    md = ["| " + " | ".join(COLS) + " |", "|" + "---|" * len(COLS)]
    md += ["| " + " | ".join(fmt(r[c]) for c in COLS) + " |" for r in rows]
    if a.md:
        Path(a.md).write_text("\n".join(md) + "\n")
    elif not a.csv:
        print("\n".join(md))


if __name__ == "__main__":
    main()
