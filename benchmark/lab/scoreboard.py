#!/usr/bin/env python3
"""Lab evidence scoreboard: accuracy vs DOFs vs wall time per discretisation family.

Reads only existing measured rows (benchmark/lab/*.csv, the copies of the round-4
sweeps in benchmark/lab/data/ with their checksums, plus the literal KIN6D N1, MF1
and NEQ1 handoff rows below) and writes SCOREBOARD.md, figures/scoreboard.png
(L1/L3/L6 panels of the BEN1 rows) and figures/scoreboard_r4.png (round-4 evidence).
No solver is run.  Output is deterministic.

    MPLBACKEND=Agg python3 benchmark/lab/scoreboard.py            # regenerate
    MPLBACKEND=Agg python3 benchmark/lab/scoreboard.py --check    # fail if stale
"""
import argparse
import csv
import hashlib
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MD = HERE / "SCOREBOARD.md"
PNG = HERE / "figures" / "scoreboard.png"
PNG4 = HERE / "figures" / "scoreboard_r4.png"
DATA = HERE / "data"

FAMILY = {
    "gll-mf": "GLL matrix-free (p-MG)",
    "gll-polar": "GLL polar chart, assembled (comparison)",
    "tet": "tet P1-P3 assembled",
    "bspline": "tensor B-spline assembled",
    "kin6d": "KIN6D P2/P3 (native LU, N1)",
    "vmecpp": "VMEC++ (reference)",
    "gvec": "GVEC (reference)",
    "gll-fourier-mf": "GLL x Fourier matrix-free",
    "spline-fourier-mf": "B-spline x Fourier matrix-free",
    "spline-mf": "B-spline IGA matrix-free",
    "gll-nk": "GLL Newton-Krylov + p-MG (NEQ1)",
    "fem-bem": "FEM + BEM free boundary (dense BEM/DtN)",
    "gll-island": "GLL island-adapted quad mesh (UNS1)",
    "gll-ogrid": "GLL unaligned O-grid (UNS1)",
    "gll-aligned": "GLL field-aligned + mortar seam (SEAM1)",
}
STYLE = {  # colour (Okabe-Ito), marker, filled
    "gll-mf": ("#0072B2", "o"),
    "gll-polar": ("#56B4E9", "P"),
    "tet": ("#D55E00", "s"),
    "bspline": ("#009E73", "^"),
    "kin6d": ("#CC79A7", "v"),
    "vmecpp": ("#000000", "D"),
}
MISSING = [
    ("Same-host, same-compiler, same-geometry timing of assembled tet/B-spline "
     "against the matrix-free families (assembled rows: faepmac2, affine mesh)",
     "not measured; cross-family time ratios with the assembled rows are [H]"),
    ("Assembled tet/B-spline beyond 1e-5 on L1", "not measured"),
    ("p-robust smoother / weighted quadrature for B-spline MF (spline h-MG "
     "iterations grow 6/14/27/54 for p = 2..5)", "not done (IGA-AX)"),
    ("Chi-robust preconditioner for anisotropy >= 1e7 (p-MG not converged in 500 "
     "iterations)", "no solver yet (UNS1)"),
    ("Non-axisymmetric (3D) equilibrium Newton solve; lab-layer admission of the "
     "N1 L6 cases", "no data yet; NEQ1 is axisymmetric Solov'ev only (lane NEQ1)"),
    ("FEM+BEM beyond 6144 trace DOFs, parallel DtN build, fast BEM",
     "not measured; dense DtN build is serial and ~N^3 (lane FB1)"),
    ("3D multi-helicity / stochastic island case",
     "no data yet; UNS1 is 2D, single helicity"),
    ("Field-aligned mesh at anisotropy 1e6", "PCG fails on both layouts (SEAM1 L3)"),
    ("GVEC integer cases", "incomplete (N1)"),
]

# (resolution, coefficients, B_sample_relative vs N1 oracle, native 1-core phase-sum s,
#  line in N1_HANDOFF.md).  Case L = Landreman integer, KIN6D-main rows, 1 core.
KIN6D_L = [
    ("p2-n1-b8", 369, 0.184988, 0.0396473, 15),
    ("p2-n2-b12", 975, 0.0298218, 0.278583, 16),
    ("p2-n4-b12", 1755, 0.00689125, 1.49858, 17),
    ("p2-n4-b16", 3699, 0.00341811, 3.81389, 18),
    ("p2-n4-b8", 1107, 0.0111185, 0.751717, 19),
    ("p2-n6-b12", 2535, 0.00684697, 3.07832, 20),
    ("p3-n1-b8", 765, 0.184769, 0.174801, 21),
    ("p3-n2-b12", 2040, 0.0285327, 1.50674, 22),
    ("p3-n4-b12", 3672, 0.000816506, 9.92195, 23),
    ("p3-n4-b16", 7965, 0.00112077, 34.3715, 24),
    ("p3-n4-b8", 2295, 0.00161565, 4.16732, 25),
    ("p3-n6-b12", 5304, 0.000637827, 20.8144, 26),
]
KIN6D_SRC = "kin6d_islands/N1_HANDOFF.md"
KIN6D_ENV = ("node7", "GNU12.2 Release", "KIN6D-main (job 3355100)")
MF1_SRC = "fortfem_lab/MF1_HANDOFF.md"

CSVS = [
    "baselines-20261022.csv",
    "additional-20261022.csv",
    "mf-baselines-20261022.csv",
    "mf-extra-20261022.csv",
    "mf-high-20261022.csv",
    "mf-anisotropic-high-20261022.csv",
]


def short_compiler(c):
    return "GCC " + c.split(")")[-1].strip() if "GNU" in c else c


def load_rows():
    rows = []
    for name in CSVS:
        with open(HERE / name, newline="") as f:
            for i, r in enumerate(csv.DictReader(f), start=2):
                if r["status"] != "V":
                    continue
                ratio = float(r["anisotropy_ratio"])
                err = float(r["rel_l2"])
                rows.append(dict(
                    problem=r["problem"], crit=int(r["criterion"]),
                    fam=r["family"], p=int(r["p"]), dofs=int(r["dofs"]),
                    metric="rel L2 vs exact", err=err, wall=float(r["wall_s"]),
                    threads=int(r["threads"]), host=r["host"].split(".")[0],
                    compiler=short_compiler(r["compiler"]),
                    rev=r["revision"][:8], src=f"benchmark/lab/{name}:{i}",
                    note=(f"ratio {ratio:g}" if ratio != 1.0 else "")
                    + ("; err>=0.3, not resolved" if err >= 0.3 else "")))
    for res, c, b, w, ln in KIN6D_L:
        rows.append(dict(
            problem="L6", crit=4, fam="kin6d", p=int(res[1]), dofs=c,
            metric="B sample rel vs N1 oracle", err=b, wall=w, threads=1,
            host=KIN6D_ENV[0], compiler=KIN6D_ENV[1], rev="KIN6D-main",
            src=f"{KIN6D_SRC}:{ln}",
            note=f"{res}; wall = native phase sum (lower bound)"))
    with open(HERE / "l6-n1-reference-20261022.csv", newline="") as f:
        for i, r in enumerate(csv.DictReader(f), start=2):
            if r["code"] not in ("vmecpp", "gvec") or r["case"] == "qa":
                continue
            ok = r["status"] == "sampled" and r["B_sample_relative"]
            if r["threads"] != "1":
                continue
            rows.append(dict(
                problem="L6", crit=4, fam=r["code"], p=0,
                dofs=int(r["coefficient_dofs"]), metric="B sample rel vs N1 oracle",
                err=float(r["B_sample_relative"]) if ok else None,
                wall=float(r["native_solve_wall_s"]) if ok else None,
                threads=1, host=r["host"], compiler="external package",
                rev=r["version"], src=f"benchmark/lab/l6-n1-reference-20261022.csv:{i}",
                note=f"{r['case']} res {r['resolution']}"
                + ("" if ok else "; GVEC incomplete, no error/time")))
    rows.append(dict(
        problem="L1", crit=1, fam="gll-mf", p=0, dofs=117649,
        metric="(solve only; no error recorded)", err=None, wall=2.051032,
        threads=16, host="sCluster node4", compiler="GCC (MF1 sweep)",
        rev="25ec09e2", src=f"{MF1_SRC}:27",
        note="MF1 scaling sweep 3355311, 7 PCG its, 3.536 s at 1 thread; table only"))
    rows.extend(load_r4())
    return rows


# Round-4 evidence.  Measured CSVs are copies in data/ (provenance and checksums in
# data/README.md and data/SHA256SUMS); NEQ1 rows are literal from its handoff.
IGAAX_ENV = ("node7", "GCC 12.2", "3bf52b95", 3356269)
IGAAX_FAM = {"gll": "gll-mf", "gll-fourier": "gll-fourier-mf",
             "spline-fourier": "spline-fourier-mf", "spline": "spline-mf",
             "polar": "gll-polar"}
FB1_ENV = ("node6", "GCC 12.2", "46fefb79", 3356163)
SEAM_FAM = {"gll-unaligned": "gll-mf", "gll-aligned-mortar-seam": "gll-aligned"}
UNS1_ENV = ("3c2fe688", "GCC 12.2")
NEQ1_SRC = "fortfem_lab/NEQ1_HANDOFF.md"
# (matched max error level, lab (n,p), lab DOFs, lab s, VMEC++ (mpol,ns) or None,
#  VMEC++ s or None, line in NEQ1_HANDOFF.md); case A exact Solov'ev, 8 threads, node13,
#  jobs 3355912/3355913.  The level is the matched-accuracy target, not a measured error.
NEQ1_A = [
    (1e-4, "(3,5)", 1156, 0.025, "(8,51)", 0.033, 27),
    (3e-5, "(4,5)", 2041, 0.032, "(10,101)", 0.058, 28),
    (1e-5, "(2,8)", 1313, 0.036, "(16,201)", 0.178, 29),
    (1e-6, "(3,7)", 2248, 0.042, None, None, 30),
    (1e-8, "(6,8)", 11617, 0.127, None, None, 31),
]
# case B (nonlinear profile, self-convergence reference): matched level, lab s, VMEC++ s
NEQ1_B = [(1e-4, 0.027, 0.042), (3e-5, 0.037, 0.072), (1e-5, 0.044, 0.164)]


def read_csv(name):
    with open(DATA / name, newline="") as f:
        return list(csv.DictReader(f))


def fnum(x):
    x = x.strip()
    return float(x) if x else None


def load_r4():
    rows = []
    host, comp, rev, job = IGAAX_ENV
    for name, prob in (("igaax-scan-3356269.csv", "L1"),
                       ("igaax-threads-3356269.csv", "L1"),
                       ("igaax-disk-3356269.csv", "disk")):
        for i, r in enumerate(read_csv(name), start=2):
            rows.append(dict(
                r4=True, problem=prob, crit=1 if prob == "L1" else 2,
                fam=IGAAX_FAM[r["family"]], p=int(r["p"]), dofs=int(r["dofs"]),
                metric="rel L2 vs exact", err=fnum(r["rel_l2"]),
                wall=fnum(r["total_s"]), threads=int(r["threads"]), host=host,
                compiler=comp, rev=rev, src=f"benchmark/lab/data/{name}:{i}",
                note=f"job {job}; "
                + (f"{r['iterations']} PCG its" if int(r["iterations"]) else "direct")
                + (f"; axis L2 {fnum(r['axis_l2']):.2e}" if prob == "disk" else "")))
    host, comp, rev, job = FB1_ENV
    seen = set()
    for i, r in enumerate(read_csv("fb1-sweep-3356163.csv"), start=2):
        if (r["n"], r["threads"]) in seen:
            continue  # the sweep log repeats some thread counts
        seen.add((r["n"], r["threads"]))
        wall = sum(fnum(r[k]) for k in ("t_bem", "t_dtn", "t_solve"))
        rows.append(dict(
            r4=True, problem="FB1", crit=1, fam="fem-bem", p=int(r["p"]),
            dofs=int(r["dofs"]), metric="rel L2 vs transmission oracle",
            err=fnum(r["relL2"]), wall=wall, threads=int(r["threads"]), host=host,
            compiler=comp, rev=rev, src=f"benchmark/lab/data/fb1-sweep-3356163.csv:{i}",
            note=(f"job {job}; trace {r['trace']}; its {r['its']}; wall = t_bem "
                  f"{fnum(r['t_bem']):.3g} + t_dtn {fnum(r['t_dtn']):.3g} + t_solve "
                  f"{fnum(r['t_solve']):.3g}")))
    for i, r in enumerate(read_csv("seam1-l3-3356493.csv"), start=2):
        ok = r["status"] == "V" and r["rel_l2"]
        key = r["geometry"].replace("nfp5-unit-disk-", "")
        rows.append(dict(
            r4=True, problem="L3", crit=3, fam=SEAM_FAM[key], p=int(r["p"]),
            dofs=int(r["dofs"] or 0), metric="rel L2 vs exact",
            err=fnum(r["rel_l2"]) if ok else None, wall=fnum(r["wall_s"]) if ok else None,
            threads=int(r["threads"]), host=r["host"], compiler="GCC 12.2",
            rev=r["revision"][:8], src=f"benchmark/lab/data/seam1-l3-3356493.csv:{i}",
            note=f"job {r['job_id']}; disk; chi ratio {float(r['anisotropy_ratio']):g}; "
            + (f"{r['iterations']} PCG its" if ok else f"{r['status']}")))
    for i, r in enumerate(read_csv("uns1-sweep-3356137.csv"), start=2):
        top = 4 if r["kind"] == "island" else 16
        if int(r["level"]) != top or r["coarse"] not in ("ldlt", "direct"):
            continue
        conv = r["status"].strip() == "0"
        rows.append(dict(
            r4=True, problem="UNS1", crit=3,
            fam="gll-island" if r["kind"] == "island" else "gll-ogrid",
            p=int(r["p"]), dofs=int(r["ndof"]), metric="rel L2 vs island oracle",
            err=fnum(r["rel_l2"]), wall=fnum(r["t_setup"]) + fnum(r["t_solve"]),
            threads=int(r["threads"]), host=r["host"], compiler=UNS1_ENV[1],
            rev=r["revision"], src=f"benchmark/lab/data/uns1-sweep-3356137.csv:{i}",
            note=f"job {r['job']}; 2D disk, chi {fnum(r['ratio']):g}; "
            + ("direct LDLT" if r["coarse"] == "direct" else
               f"p-MG/LDLT {int(r['its'])} its" + ("" if conv else ", not converged"))))
    for lvl, lab, dofs, tl, vm, tv, ln in NEQ1_A:
        rows.append(dict(
            r4=True, problem="NEQ1", crit=4, fam="gll-nk", p=int(lab.strip("()").split(",")[1]), dofs=dofs,
            metric="matched max error level (Solov'ev A)", err=lvl, wall=tl,
            threads=8, host="node13", compiler="GCC 12.2", rev="95ca9ddc",
            src=f"{NEQ1_SRC}:{ln}", note=f"(n,p) {lab}; jobs 3355912/3355913"))
        if vm:
            rows.append(dict(
                r4=True, problem="NEQ1", crit=4, fam="vmecpp", p=0, dofs=0,
                metric="matched max error level (Solov'ev A)", err=lvl, wall=tv,
                threads=8, host="node13", compiler="external package",
                rev="0.9.1", src=f"{NEQ1_SRC}:{ln}",
                note=f"(mpol,ns) {vm}; native solve only; same node"))
    for k, (lvl, tl, tv) in enumerate(NEQ1_B):
        for fam, t in (("gll-nk", tl), ("vmecpp", tv)):
            rows.append(dict(
                r4=True, problem="NEQ1", crit=4, fam=fam, p=0, dofs=0,
                metric="matched max error level (nonlinear profile B)", err=lvl,
                wall=t, threads=8, host="node13", compiler="GCC 12.2",
                rev="95ca9ddc" if fam == "gll-nk" else "0.9.1",
                src=f"{NEQ1_SRC}:33", note="case B, self-convergence reference n8 p12"))
    return rows


def order(rows):
    fo = list(FAMILY)
    return sorted(rows, key=lambda r: (r["problem"], fo.index(r["fam"]), r["dofs"],
                                       r["p"], r["note"]))


def fmt(x, spec):
    return "-" if x is None else format(x, spec)


def headline(rows):
    """Best wall time per family on L1 that reaches a given relative L2 error."""
    out = []
    for thr in (1e-3, 1e-5, 1e-7):
        for fam in ("gll-mf", "tet", "bspline"):
            c = [r for r in rows if r["problem"] == "L1" and r["fam"] == fam
                 and not r.get("r4") and r["err"] is not None and r["err"] <= thr]
            if c:
                r = min(c, key=lambda r: r["wall"])
                out.append((thr, fam, r))
            else:
                out.append((thr, fam, None))
    return out


def headline_r4(rows):
    """Best L1 wall time per matrix-free family at 16 threads, one exclusive node."""
    fams = ("gll-mf", "gll-fourier-mf", "spline-fourier-mf", "spline-mf")
    out = []
    for thr in (1e-4, 1e-6, 1e-8):
        for fam in fams:
            c = [r for r in rows if r.get("r4") and r["problem"] == "L1"
                 and r["fam"] == fam and r["err"] is not None and r["err"] <= thr]
            out.append((thr, fam, min(c, key=lambda r: r["wall"]) if c else None))
    return out


def write_md(rows):
    L = []
    a = L.append
    a("# Lab evidence scoreboard")
    a("")
    a("Generated by `benchmark/lab/scoreboard.py` from existing measured rows only "
      "(no new runs). Do not edit by hand.")
    a("")
    a("## Reading guide")
    a("")
    a("- Error: L1/L3 relative L2 against the exact manufactured solution; "
      "L6 relative sampled B against the N1 exact oracle (a different quantity, "
      "never compare it with an L1/L3 number).")
    a("- Wall time is process wall time (assembly, setup, solve, error integral) for "
      "the lab CSVs. KIN6D rows are native phase-sum seconds on one core, a lower "
      "bound. VMEC++ rows are native solve seconds on one core.")
    a("- **Like-for-like groups.** Assembled tet/B-spline rows ran on faepmac2 with "
      "GCC 16.1 (revision 5a48864d), matrix-free GLL rows on faepmac1 with GCC 16.2, "
      "KIN6D/VMEC++/GVEC on sCluster node7 (GCC 12.2 / external builds). The "
      "round-4 matrix-free rows (IGA-AX, FB1, SEAM1, UNS1, NEQ1) are separate "
      "sCluster sweeps and are only comparable within their own job. Different "
      "hosts, compilers, geometries (affine-tetra vs curved torus) and solvers mean "
      "**every cross-family time ratio is [H], not [V]**. Only DOFs-vs-error within "
      "one family is a clean measurement.")
    a("- Criterion column: 1 FEM+BEM-ready interior, 2 axis regularity, "
      "3 field-alignment optional, 4 time-to-solution vs VMEC-class codes.")
    a("")
    a("## Headline, round-1 rows: best L1 wall time to reach a given error "
      "[H across families]")
    a("")
    a("| target rel L2 | family | p | DOFs | rel L2 | wall s | threads | host | source |")
    a("|---|---|---:|---:|---:|---:|---:|---|---|")
    for thr, fam, r in headline(rows):
        if r is None:
            a(f"| <= {thr:g} | {FAMILY[fam]} | - | - | not reached | - | - | - | "
              "no row reaches this |")
        else:
            a(f"| <= {thr:g} | {FAMILY[fam]} | {r['p']} | {r['dofs']} | "
              f"{r['err']:.3e} | {r['wall']:.3g} | {r['threads']} | {r['host']} | "
              f"{r['src']} |")
    a("")
    a("The matrix-free GLL rows are the only family that reaches 1e-7 within the "
      "measured range; tet P3 and the B-spline rows were not run that far, so "
      "\"not reached\" means \"not measured\", not \"cannot\". Tet P3 reaches 5.4e-6 only at "
      "103860 DOFs and 328 s.")
    a("")
    a("## Figure")
    a("")
    a("![scoreboard](figures/scoreboard.png)")
    a("")
    a("Top row x = DOFs, bottom row x = wall time, y = relative error, log-log. "
      "L3 panel: anisotropy ratio 1e6 only (filled: L3 curved/affine torus; hollow: "
      "L3a 2D square P1). L6 panel: reference VMEC++ (black diamonds) and KIN6D "
      "(Landreman integer case); there is no lab-layer solve of the N1 cases (NEQ1 "
      "solves axisymmetric Solov'ev only, see below).")
    a("")
    a("## Round-4 evidence (2026-10-23): matrix-free families, FEM+BEM, aligned seam, "
      "islands, Newton-Krylov")
    a("")
    a("### Matrix-free families on one exclusive node (L1 shaped nfp=5 torus, "
      "16 threads)")
    a("")
    a("IGA-AX job 3356269, sCluster node7, GCC 12.2, one source revision, all four "
      "families on the same host, compiler and geometry, so these times are a clean "
      "within-lab comparison. Wall time includes setup; best row reaching the target.")
    a("")
    a("| target rel L2 | family | p | DOFs | rel L2 | wall s | source |")
    a("|---|---|---:|---:|---:|---:|---|")
    for thr, fam, r in headline_r4(rows):
        if r is None:
            a(f"| <= {thr:g} | {FAMILY[fam]} | - | - | not reached | - | "
              "no row reaches this |")
        else:
            a(f"| <= {thr:g} | {FAMILY[fam]} | {r['p']} | {r['dofs']} | "
              f"{r['err']:.3e} | {r['wall']:.3g} | {r['src']} |")
    a("")
    a("- [V] Within these rows GLL tensor p-MG is the fastest to every target: at "
      "1e-8 it needs 1.26 s against 2.51 s (GLL x Fourier) and 10.6 s (B-spline x "
      "Fourier); pure B-spline IGA does not reach 1e-8 (best 1.06e-8 is above the "
      "target, see the scan CSV: 354k DOFs, 49 s).")
    a("- [V] Thread scaling 1 to 16 (data/igaax-threads-3356269.csv): GLL 252k DOFs "
      "3.2x, B-spline x Fourier 120k DOFs 6.8x, B-spline 60k DOFs 8.7x, GLL x Fourier "
      "108k DOFs 2.7x. GLL leads in absolute time but scales worst; [H] its setup "
      "is nearly serial (IGA-AX handoff). The spline families scale better and "
      "still lose.")
    a("- [V] Criterion 2 (axis, 2D disk, data/igaax-disk-3356269.csv): the "
      "axis-region L2 error of the through-axis GLL and B-spline O-grids is at most "
      "6.6e-8 on the coarsest meshes and falls to 4e-12 (GLL) and 1.3e-11 "
      "(B-spline) on the finest; the polar chart comparison does the same (at most "
      "2.3e-7 down to 9e-14). Both families are axis-regular on the O-grid. [H] "
      "The causes of the spline cost (quadrature, non-p-robust MG) are discussed in "
      "IGA-AX, not measured here.")
    a("")
    a("### FEM+BEM free boundary (FB1, job 3356163, node6, 16 threads)")
    a("")
    a("- [V] Interior FEM + dense BEM/DtN transmission problem: relative L2 3.4e-3 "
      "to 6.9e-5 for 384 to 6144 trace DOFs. At 6144 trace DOFs the DtN dense build "
      "takes 127.5 s against 5.7 s BEM assembly and 1.7 s for the coupled PCG solve.")
    a("- [V] BEM assembly scales 1 to 16 threads 12.5 s to 1.0 s (n=2) and 28.6 s to "
      "2.6 s (n=3); the DtN dense build does not scale (22.9 s to 21.9 s at n=3) and "
      "grows ~N^3. The dense DtN build, not the FEM solve, is the cost to remove.")
    a("")
    a("### Field-aligned seam against unaligned (SEAM1 L3 disk, job 3356493, node17, "
      "4 threads)")
    a("")
    a("- [V] The field-aligned mesh with a mortar seam is 100 to 360 times slower "
      "than the unaligned GLL disk (p6, ratio 1: 119 s against 0.33 s) and gains "
      "at most a factor 1.4 in error (ratio 1e3, p = 2 and 6); at ratio 1 it is "
      "less accurate for every p (up to 4.4 times at p = 6). At ratio 1e6 PCG "
      "fails on both layouts. Alignment stays an optional accelerator "
      "(criterion 3), and this implementation is not one. [H] The cost of the "
      "mortar path, not the alignment itself, is the likely cause.")
    a("")
    a("### Islands, 2D single helicity (UNS1, jobs 3356137_*, 16 threads, largest "
      "meshes: island L4 against O-grid L16 at 0.8 times the DOFs)")
    a("")
    a("- [V] Up to chi = 1e5 the unaligned O-grid matches or beats the island mesh; "
      "at 1e9 the island mesh is 2.8 to 3.6 times more accurate (direct solve, "
      "p = 2/4/6). p-MG is not chi-robust: p = 4 and 6 do not converge in 500 "
      "iterations from chi = 1e7 on (rows marked not converged); the direct LDLT "
      "reference is slow. No fast stiff-anisotropy solver exists yet.")
    a("")
    a("### Nonlinear equilibrium (NEQ1, axisymmetric Solov'ev, 8 threads, node13)")
    a("")
    a("- [V on the handoff rows] Newton-Krylov with GLL p-MG reaches a matched max "
      "error of 1e-4/3e-5/1e-5 in 0.025/0.032/0.036 s against 0.033/0.058/0.178 s "
      "for VMEC++ (same node, native solve), and reaches 1e-6 and 1e-8 in 0.042 and "
      "0.127 s where the VMEC++ scan (ns <= 201) does not. These are runs of 25 to "
      "180 ms, so scheduling noise matters, and the error levels are matched-accuracy "
      "targets read from the handoff table, not individually measured errors. "
      "Axisymmetric only: no 3D equilibrium and no N1 L6 case has been solved on "
      "the lab layer.")
    a("")
    a("![scoreboard r4](figures/scoreboard_r4.png)")
    a("")
    a("Panels: (a, b) L1 accuracy against DOFs and time, thin markers all rows and "
      "lines the lower-left frontier per family, assembled BEN1 rows hollow (faepmac2, "
      "affine mesh, 4 threads: not the same host); (c) strong scaling; (d) FB1 cost "
      "split; (e) NEQ1 against VMEC++; (f) SEAM1 L3 disk.")
    a("")
    a("## Missing evidence")
    a("")
    a("| family / capability | status |")
    a("|---|---|")
    for k, v in MISSING:
        a(f"| {k} | {v} |")
    a("")
    a("## Full table")
    a("")
    a("| problem | crit | family | p | DOFs | metric | error | wall s | thr | host | "
      "compiler | revision | source | note |")
    a("|---|---:|---|---:|---:|---|---:|---:|---:|---|---|---|---|---|")
    for r in order(rows):
        a(f"| {r['problem']} | {r['crit']} | {FAMILY[r['fam']]} | "
          f"{r['p'] or '-'} | {r['dofs'] or '-'} | {r['metric']} | {fmt(r['err'], '.3e')} | "
          f"{fmt(r['wall'], '.4g')} | {r['threads']} | {r['host']} | {r['compiler']} | "
          f"{r['rev']} | {r['src']} | {r['note'].lstrip('; ')} |")
    a("")
    a("## Limits")
    a("")
    a("- No row is measured on a common host, compiler or geometry across families; "
      "time comparisons are [H].")
    a("- L6 has no lab-layer solution of the N1 cases and no matched-accuracy "
      "admission on the lab layer; the "
      "KIN6D-vs-VMEC++ gap (P3 2295 coefficients, 1.9-4.2 s vs 0.023-0.036 s) is "
      "reported from the N1 study, not re-measured. GVEC integer cases are "
      "incomplete and QA has no exact oracle, so neither is plotted.")
    a("- Source line numbers refer to the CSV row (header = line 1) or to the cited "
      "handoff file under the handoff directory, which is outside this repository. "
      "Rows from benchmark/lab/data/ keep the original job output; "
      "data/README.md records job, host, revision and SHA256.")
    a("- UNS1 and FB1 wall times are phase sums (UNS1: t_setup + t_solve; FB1: t_bem "
      "+ t_dtn + t_solve), not process wall time; FB1 omits mesh and "
      "preconditioner setup columns of the CSV. UNS1 lists only the largest mesh "
      "of each family and the p-MG/LDLT and direct solvers (all 288 rows are in the "
      "CSV). NEQ1 error values are matched-accuracy levels.")
    a("- IGA-AX scan and thread rows overlap (the same configuration can appear in "
      "both); the headline takes the best row of either.")
    a("- L3 at 1000 DOFs: tet P2 (affine mesh) is more accurate than GLL p2 "
      "(curved mesh); the geometries differ, so this is not a method ranking.")
    a("- MF1 117649-DOF row is solve-only with no error recorded; it is listed "
      "but not plotted.")
    a("")
    a("Regenerate: `MPLBACKEND=Agg python3 benchmark/lab/scoreboard.py`; "
      "verify: `... scoreboard.py --check`.")
    return "\n".join(L) + "\n"


def draw(rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.3,
                         "svg.hashsalt": "scoreboard"})
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.6))
    cols = [("L1", "L1 Poisson torus (crit. 1)"),
            ("L3", "L3 anisotropic, ratio 1e6 (crit. 3)"),
            ("L6", "L6 nested equilibrium (crit. 4)")]
    for j, (prob, title) in enumerate(cols):
        for i, xkey in enumerate(("dofs", "wall")):
            ax = axes[i, j]
            for r in order(rows):
                if r["err"] is None or r["wall"] is None or r.get("r4"):
                    continue
                if prob == "L3":
                    if r["problem"] not in ("L3", "L3a") or "ratio 1e+06" not in r["note"]:
                        continue
                elif r["problem"] != prob:
                    continue
                col, mk = STYLE[r["fam"]]
                hollow = r["problem"] == "L3a"
                ax.plot(r[xkey], r["err"], mk, color=col, ms=7,
                        mfc="none" if hollow else col, mew=1.4, alpha=0.9)
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlabel("DOFs" if xkey == "dofs" else "wall time [s]")
            ax.set_ylabel("rel. L2 error" if prob != "L6" else "rel. B sample error")
            if i == 0:
                ax.set_title(title, fontsize=10)
            if prob == "L6":
                ax.text(0.97, 0.96, "no lab-layer solve of the N1 cases",
                        transform=ax.transAxes, fontsize=8, style="italic",
                        ha="right", va="top")
            if prob == "L3" and i == 0:
                ax.text(0.97, 0.04, "filled: L3 3D, hollow: L3a 2D P1",
                        transform=ax.transAxes, fontsize=8, ha="right")
    handles = []
    for k, (col, mk) in STYLE.items():
        if k == "gll-polar":
            continue  # L5-only comparison family, not in these panels
        handles.append(plt.Line2D([], [], marker=mk, color=col, ls="", ms=7,
                                  label=FAMILY[k]))
    handles.append(plt.Line2D([], [], marker="x", color="gray", ls="", ms=0,
                              label="hosts/compilers differ: times are [H]"))
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=8.5)
    fig.suptitle("FortFEM lab scoreboard: error vs DOFs and wall time "
                 "(existing measured rows only)", fontsize=11)
    fig.tight_layout(rect=(0, 0.1, 1, 0.97))
    PNG.parent.mkdir(exist_ok=True)
    fig.savefig(PNG, dpi=110, metadata={"Software": None})
    plt.close(fig)


OI = {"blue": "#0072B2", "sky": "#56B4E9", "green": "#009E73", "orange": "#E69F00",
      "verm": "#D55E00", "purple": "#CC79A7", "black": "#000000"}


def pareto(pts):
    """Lower-left frontier of (cost, error) points."""
    out, best = [], float("inf")
    for c, e in sorted(pts):
        if e < best:
            out.append((c, e))
            best = e
    return out


def draw_r4(rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 9})
    fam_style = {
        "gll-mf": ("GLL tensor, p-MG (MF)", OI["blue"], "o", "-"),
        "gll-fourier-mf": ("GLL x Fourier (MF)", OI["sky"], "D", "--"),
        "spline-fourier-mf": ("B-spline x Fourier (MF)", OI["green"], "^", "-."),
        "spline-mf": ("B-spline IGA (MF)", OI["orange"], "v", ":"),
    }
    ben_style = {"tet": ("tet P1-P3 assembled*", OI["verm"], "s"),
                 "bspline": ("B-spline assembled*", OI["purple"], "P")}
    l1 = [r for r in rows if r.get("r4") and r["problem"] == "L1"]
    fig, ax = plt.subplots(2, 3, figsize=(15, 9.2))
    (a, b, c), (d, e, g) = ax
    for fam, (lab, col, mk, ls) in fam_style.items():
        r = [x for x in l1 if x["fam"] == fam]
        dofs, err, t = ([x[k] for x in r] for k in ("dofs", "err", "wall"))
        a.scatter(dofs, err, c=col, marker=mk, s=22, alpha=0.45, edgecolors="none")
        b.scatter(t, err, c=col, marker=mk, s=22, alpha=0.45, edgecolors="none")
        for axx, xs in ((a, dofs), (b, t)):
            p = pareto(list(zip(xs, err)))
            axx.plot([q[0] for q in p], [q[1] for q in p], color=col, ls=ls,
                     marker=mk, ms=6, lw=1.8, label=lab)
    for fam, (lab, col, mk) in ben_style.items():
        r = [x for x in rows if not x.get("r4") and x["problem"] == "L1"
             and x["fam"] == fam and x["err"] is not None
             and x["host"].startswith("faepmac2")]
        a.scatter([x["dofs"] for x in r], [x["err"] for x in r], facecolors="none",
                  edgecolors=col, marker=mk, s=45, label=lab)
        b.scatter([x["wall"] for x in r], [x["err"] for x in r], facecolors="none",
                  edgecolors=col, marker=mk, s=45, label=lab)
    for axx, xl in ((a, "DOFs"), (b, "wall time incl. setup [s]")):
        axx.set_xscale("log")
        axx.set_yscale("log")
        axx.set_xlabel(xl)
        axx.set_ylabel("relative L2 error")
        axx.grid(True, which="major", alpha=0.3)
        for tgt in (1e-6, 1e-8):
            axx.axhline(tgt, color="0.6", lw=0.8, ls=(0, (2, 2)))
    a.set_title("(a) L1 shaped nfp=5 torus: accuracy vs DOFs")
    b.set_title("(b) L1 torus: accuracy vs time (MF: 16 thr sCluster)")
    b.legend(fontsize=8, loc="lower left", framealpha=0.9)

    # (c) strong scaling: IGA-AX thread series, FB1 phases
    th = read_csv("igaax-threads-3356269.csv")
    fam_of = {"gll": "gll-mf", "gll-fourier": "gll-fourier-mf",
              "spline-fourier": "spline-fourier-mf", "spline": "spline-mf"}
    for fam, (lab, col, mk, ls) in fam_style.items():
        groups = {}
        for x in th:
            if fam_of[x["family"]] == fam:
                groups.setdefault((x["p"], x["n"], x["m"], x["nphi"]), []).append(x)
        for grp in groups.values():
            grp.sort(key=lambda x: int(x["threads"]))
            if int(grp[0]["threads"]) != 1:
                continue
            t1 = fnum(grp[0]["total_s"])
            big = fnum(grp[0]["rel_l2"]) < 1e-7
            c.plot([int(x["threads"]) for x in grp],
                   [t1 / fnum(x["total_s"]) for x in grp], color=col, marker=mk,
                   ls=ls, lw=1.6 if big else 1.0, alpha=1.0 if big else 0.5)
    fb, seen = [], set()
    for x in read_csv("fb1-sweep-3356163.csv"):
        if (x["n"], x["threads"]) not in seen:
            seen.add((x["n"], x["threads"]))
            fb.append(x)
    for n, lw in (("2", 1.0), ("3", 1.8)):
        r = sorted([x for x in fb if x["n"] == n], key=lambda x: int(x["threads"]))
        T = [int(x["threads"]) for x in r]
        t1b, t1d = fnum(r[0]["t_bem"]), fnum(r[0]["t_dtn"])
        c.plot(T, [t1b / fnum(x["t_bem"]) for x in r], color=OI["black"], marker="x",
               lw=lw)
        c.plot(T, [t1d / fnum(x["t_dtn"]) for x in r], color=OI["verm"], marker="s",
               lw=lw, ls=":")
    c.plot([1, 16], [1, 16], color="0.6", lw=0.8)
    c.set_xscale("log", base=2)
    c.set_yscale("log", base=2)
    for setter, tick in ((c.set_xticks, c.set_xticklabels),
                         (c.set_yticks, c.set_yticklabels)):
        setter([1, 2, 4, 8, 16])
        tick(["1", "2", "4", "8", "16"])
    c.set_xlabel("OpenMP threads")
    c.set_ylabel("speedup vs 1 thread")
    c.set_title("(c) strong scaling (thick: largest case)")
    c.text(1.1, 9.5, "BEM assembly", color=OI["black"], fontsize=8)
    c.text(5, 0.82, "DtN dense build (serial)", color=OI["verm"], fontsize=8)
    c.text(1.1, 3.4, "MF families:\ncolours as (a)", fontsize=8, color="0.3")
    c.grid(True, alpha=0.3)

    # (d) FB1 cost split at 16 threads
    byn = {int(x["trace"]): x for x in fb if x["threads"] == "16"}
    N = sorted(byn)
    for key, lab, col, mk, ls in (
            ("t_bem", "BEM assembly (dense)", OI["black"], "x", "-"),
            ("t_dtn", "DtN/Steklov dense build", OI["verm"], "s", ":"),
            ("t_solve", "coupled PCG solve", OI["blue"], "o", "--")):
        d.plot(N, [fnum(byn[k][key]) for k in N], color=col, marker=mk, ls=ls,
               label=lab)
    ref = fnum(byn[N[-1]]["t_dtn"])
    d.plot(N, [ref * (k / N[-1]) ** 3 for k in N], color="0.6", lw=0.8, label="~N^3")
    d.set_xscale("log")
    d.set_yscale("log")
    d.set_xticks(N)
    d.set_xticklabels([str(k) for k in N])
    d.minorticks_off()
    d.set_xlabel("BEM trace DOFs N")
    d.set_ylabel("time [s] (16 thr, node6)")
    d.set_title("(d) FEM+BEM free boundary (FB1): cost split")
    d.legend(fontsize=8, loc="upper left")
    d.grid(True, alpha=0.3)
    for k in N:
        d.annotate(f"relL2 {fnum(byn[k]['relL2']):.0e}", (k, fnum(byn[k]["t_solve"])),
                   textcoords="offset points", xytext=(-10, -14), fontsize=7,
                   color=OI["blue"])

    # (e) NEQ1 against VMEC++ (case A, matched error level)
    lab_pts = [(x[3], x[0]) for x in NEQ1_A if x[3] < 0.1]
    vm_pts = [(x[5], x[0]) for x in NEQ1_A if x[5] is not None]
    e.plot([p[0] for p in lab_pts], [p[1] for p in lab_pts], color=OI["blue"],
           marker="o", label="lab GLL + Newton-Krylov")
    e.plot([p[0] for p in vm_pts], [p[1] for p in vm_pts], color=OI["black"],
           marker="D", ls="--", label="VMEC++ (ns<=201)")
    e.annotate("VMEC++ does not reach 1e-6", (0.05, 1e-6), fontsize=8)
    e.set_xscale("log")
    e.set_yscale("log")
    e.set_xticks([0.025, 0.05, 0.1, 0.2])
    e.set_xticklabels(["0.025", "0.05", "0.1", "0.2"])
    e.minorticks_off()
    e.set_xlabel("wall time [s] (8 thr, node13)")
    e.set_ylabel("matched max error vs Solov'ev oracle")
    e.set_title("(e) NEQ1 axisymmetric equilibrium vs VMEC++")
    e.legend(fontsize=8, loc="upper right")
    e.grid(True, alpha=0.3)

    # (f) SEAM1 L3: aligned mortar seam against unaligned
    sl = [r for r in rows if r.get("r4") and r["problem"] == "L3"
          and r["err"] is not None]
    for fam, col, mk, lab in (("gll-mf", OI["blue"], "o", "unaligned GLL"),
                              ("gll-aligned", OI["verm"], "s",
                               "field-aligned + mortar seam")):
        for ratio, ls in (("1", "-"), ("1000", "--")):
            r = sorted([x for x in sl if x["fam"] == fam
                        and f"chi ratio {ratio};" in x["note"]], key=lambda x: x["p"])
            g.plot([x["wall"] for x in r], [x["err"] for x in r], color=col,
                   marker=mk, ls=ls, label=f"{lab}, chi ratio {ratio}")
            for x in r:
                g.annotate(f"p{x['p']}", (x["wall"], x["err"]), fontsize=7,
                           textcoords="offset points", xytext=(4, 2))
    g.set_xscale("log")
    g.set_yscale("log")
    g.set_xlabel("wall time [s] (4 thr, node17)")
    g.set_ylabel("relative L2 error")
    g.set_title("(f) L3 disk: field-aligned seam vs unaligned")
    g.legend(fontsize=7.5, loc="center", bbox_to_anchor=(0.5, 0.42))
    g.grid(True, alpha=0.3)
    g.text(0.17, 0.97, "ratio 1e6: PCG fails (both)", transform=g.transAxes,
           fontsize=8, va="top")
    fig.suptitle("FortFEM discretisation lab, evidence scoreboard r4 (2026-10-23). "
                 "*BEN1 assembled rows: faepmac2, 4 thr, affine mesh (not same host)",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    PNG4.parent.mkdir(exist_ok=True)
    fig.savefig(PNG4, dpi=110, metadata={"Software": None})
    plt.close(fig)


def check_data():
    """Verify data/ against data/SHA256SUMS; return the list of problems."""
    sums = DATA / "SHA256SUMS"
    if not sums.exists():
        return ["data/SHA256SUMS missing"]
    bad = []
    for line in sums.read_text().splitlines():
        digest, name = line.split(None, 1)
        f = DATA / name.strip().lstrip("*")
        if not f.exists() or hashlib.sha256(f.read_bytes()).hexdigest() != digest:
            bad.append(f"data/{name.strip()} does not match SHA256SUMS")
    return bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="exit 1 if SCOREBOARD.md differs from the regenerated text "
                    "or data/ differs from data/SHA256SUMS")
    args = ap.parse_args()
    rows = load_rows()
    text = write_md(rows)
    if args.check:
        problems = check_data()
        if not MD.exists() or MD.read_text() != text:
            problems.append("SCOREBOARD.md is stale; run scoreboard.py")
        if problems:
            print("\n".join(problems), file=sys.stderr)
            return 1
        print("SCOREBOARD.md up to date")
        return 0
    MD.write_text(text)
    draw(rows)
    draw_r4(rows)
    print(f"wrote {MD.name} ({len(rows)} rows), {PNG.name} and {PNG4.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
