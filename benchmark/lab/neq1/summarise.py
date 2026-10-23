#!/usr/bin/env python3
"""Summarise a NEQ1 sCluster run (lab.sbatch) into tables and a matched-accuracy comparison.

Usage: summarise.py LABRUN_DIR_SCAN [--tag solovev|vmec-solovev-samenode] [--threads N]
LABRUN_DIR_SCAN is $RUN/lab/scan-<job> (case A) or $RUN/lab/caseB-<job> (case B, tags
caseB-r* and vmec-caseB). Lab rows: solve wall (t_solve, no sampling), DOFs, Newton and
Krylov counts, peak RSS, max common-sample error (l6.json). VMEC++ rows: native solve
wall (metrics.json), DOFs, iterations, same error metric. Timings are the minimum over
the repeats (rep tags r1..r3) for the lab; VMEC++ has one run per point.
The matched-accuracy table lists, per error target, the fastest point of each code that
reaches max error <= target at the given thread count.
"""
import argparse
import json
import re
from pathlib import Path


def rss(path):
    m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", path.read_text())
    return int(m.group(1)) / 1024 if m else None


def lab_rows(root, prefix, threads):
    rows = {}
    for rep in sorted(root.glob(prefix + "-r*")):
        for d in sorted(rep.iterdir()):
            m = re.fullmatch(r"n(\d+)-p(\d+)-t(\d+)", d.name)
            if not m or int(m[3]) != threads or not (d / "l6.json").exists():
                continue
            log = (d / "solve.log").read_text()
            f = lambda pat: float(re.search(pat, log)[1])
            key = (int(m[1]), int(m[2]))
            rec = dict(n=key[0], p=key[1], dofs=int(re.search(r"dofs=(\d+)", log)[1]),
                       newton=int(re.search(r"newton=(\d+)", log)[1]),
                       krylov=int(re.search(r"krylov=(\d+)", log)[1]),
                       t_solve=f(r"t_solve=\s*(\S+)"), t_setup=f(r"t_setup=\s*(\S+)"),
                       rss=rss(d / "time.txt"),
                       err=json.loads((d / "l6.json").read_text())["max"])
            if key not in rows or rec["t_solve"] < rows[key]["t_solve"]:
                rows[key] = rec
    return sorted(rows.values(), key=lambda r: r["dofs"])


def vmec_rows(root, tag, threads):
    rows = []
    for d in sorted((root / tag).iterdir()):
        mj = d / "metrics.json"
        if not mj.exists():
            continue
        m = json.loads(mj.read_text())
        if m.get("threads") != threads or not m.get("converged"):
            continue
        l6 = m.get("l6") or {}
        rows.append(dict(mpol=m["mpol"], ns=m["ns"], dofs=m["coefficient_dofs"],
                         iters=m["iterations"], t_solve=m["native_solve_wall_s"],
                         rss=rss(d / "time.txt"), err=l6.get("max")))
    return sorted(rows, key=lambda r: r["dofs"])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("scandir")
    ap.add_argument("--lab-tag", default="solovev")
    ap.add_argument("--vmec-tag", default="vmec-solovev-samenode")
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()
    root = Path(a.scandir)
    lab = lab_rows(root, a.lab_tag, a.threads)
    vm = [r for r in vmec_rows(root, a.vmec_tag, a.threads) if r["err"] is not None]
    print(f"## lab ({a.lab_tag}, {a.threads} threads, min of repeats)")
    print("| n | p | DOFs | Newton | Krylov | t_solve s | RSS MiB | max err |")
    print("|---|---|---|---|---|---|---|---|")
    for r in lab:
        print(f"| {r['n']} | {r['p']} | {r['dofs']} | {r['newton']} | {r['krylov']} | "
              f"{r['t_solve']:.4g} | {r['rss']:.0f} | {r['err']:.3g} |")
    print(f"\n## VMEC++ ({a.vmec_tag}, {a.threads} threads, native solve only)")
    print("| mpol | ns | DOFs | iters | t_solve s | RSS MiB | max err |")
    print("|---|---|---|---|---|---|---|")
    for r in vm:
        print(f"| {r['mpol']} | {r['ns']} | {r['dofs']} | {r['iters']} | "
              f"{r['t_solve']:.4g} | {r['rss']:.0f} | {r['err']:.3g} |")
    print("\n## fastest point reaching each target (max common-sample error)")
    print("| target | lab (n,p) | lab t s | VMEC++ (mpol,ns) | VMEC++ t s | VMEC++/lab |")
    print("|---|---|---|---|---|---|")
    for target in (1e-4, 3e-5, 1e-5, 1e-6, 1e-7, 1e-8):
        bl = min((r for r in lab if r["err"] <= target), key=lambda r: r["t_solve"],
                 default=None)
        bv = min((r for r in vm if r["err"] <= target), key=lambda r: r["t_solve"],
                 default=None)
        ls = f"({bl['n']},{bl['p']})" if bl else "-"
        vs = f"({bv['mpol']},{bv['ns']})" if bv else "-"
        lt = f"{bl['t_solve']:.4g}" if bl else "-"
        vt = f"{bv['t_solve']:.4g}" if bv else "-"
        ratio = f"{bv['t_solve'] / bl['t_solve']:.2f}" if bl and bv else "-"
        print(f"| {target:.0e} | {ls} | {lt} | {vs} | {vt} | {ratio} |")


if __name__ == "__main__":
    main()
