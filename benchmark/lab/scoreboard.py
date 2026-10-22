#!/usr/bin/env python3
"""Lab evidence scoreboard: accuracy vs DOFs vs wall time per discretisation family.

Reads only existing measured rows (benchmark/lab/*.csv plus the literal KIN6D N1 and
MF1 handoff rows below) and writes SCOREBOARD.md and figures/scoreboard.png.
No solver is run.  Output is deterministic.

    MPLBACKEND=Agg python3 benchmark/lab/scoreboard.py            # regenerate
    MPLBACKEND=Agg python3 benchmark/lab/scoreboard.py --check    # fail if stale
"""
import argparse
import csv
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MD = HERE / "SCOREBOARD.md"
PNG = HERE / "figures" / "scoreboard.png"

FAMILY = {
    "gll-mf": "GLL matrix-free (p-MG)",
    "gll-polar": "GLL polar chart, assembled (comparison)",
    "tet": "tet P1-P3 assembled",
    "bspline": "tensor B-spline assembled",
    "kin6d": "KIN6D P2/P3 (native LU, N1)",
    "vmecpp": "VMEC++ (reference)",
    "gvec": "GVEC (reference)",
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
    ("Fourier / mixed (spline x Fourier) matrix-free", "no data yet (lane IGA-AX)"),
    ("B-spline/IGA matrix-free", "no data yet (lane IGA-AX)"),
    ("Nonlinear / Newton equilibrium solve on the lab layer (L4, L6)",
     "no data yet (lane NEQ1)"),
    ("FEM + BEM free boundary", "no data yet (lane FB1)"),
    ("Non-nested / island geometry", "no data yet (lane UNS1)"),
    ("Lab-layer solution on L6 (nested equilibrium)", "no data yet (lane NEQ1)"),
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
                 and r["err"] is not None and r["err"] <= thr]
            if c:
                r = min(c, key=lambda r: r["wall"])
                out.append((thr, fam, r))
            else:
                out.append((thr, fam, None))
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
      "KIN6D/VMEC++/GVEC on sCluster node7 (GCC 12.2 / external builds). Different "
      "hosts, compilers, geometries (affine-tetra vs curved torus) and solvers mean "
      "**every cross-family time ratio is [H], not [V]**. Only DOFs-vs-error within "
      "one family is a clean measurement.")
    a("- Criterion column: 1 FEM+BEM-ready interior, 2 axis regularity, "
      "3 field-alignment optional, 4 time-to-solution vs VMEC-class codes.")
    a("")
    a("## Headline: best L1 wall time to reach a given error [H across families]")
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
      "\"not reached\" means \"not measured\", not \"cannot\". Tet P3 at 1/h=2 and "
      "B-spline p4 at 1/h=2 return errors >= 1 (see full table), a solver/conditioning "
      "outcome of the assembled ichol_row path that is recorded, not explained, here.")
    a("")
    a("## Figure")
    a("")
    a("![scoreboard](figures/scoreboard.png)")
    a("")
    a("Top row x = DOFs, bottom row x = wall time, y = relative error, log-log. "
      "L3 panel: anisotropy ratio 1e6 only (filled: L3 curved/affine torus; hollow: "
      "L3a 2D square P1). L6 panel: reference VMEC++ (black diamonds) and KIN6D "
      "(Landreman integer case); the lab layer has no L6 solve yet.")
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
          f"{r['p'] or '-'} | {r['dofs']} | {r['metric']} | {fmt(r['err'], '.3e')} | "
          f"{fmt(r['wall'], '.4g')} | {r['threads']} | {r['host']} | {r['compiler']} | "
          f"{r['rev']} | {r['src']} | {r['note'].lstrip('; ')} |")
    a("")
    a("## Limits")
    a("")
    a("- No row is measured on a common host, compiler or geometry across families; "
      "time comparisons are [H].")
    a("- L6 has no lab-layer solution and no matched-accuracy admission; the "
      "KIN6D-vs-VMEC++ gap (P3 2295 coefficients, 1.9-4.2 s vs 0.023-0.036 s) is "
      "reported from the N1 study, not re-measured. GVEC integer cases are "
      "incomplete and QA has no exact oracle, so neither is plotted.")
    a("- Source line numbers refer to the CSV row (header = line 1) or to the cited "
      "handoff file under the handoff directory, which is outside this repository.")
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
                if r["err"] is None or r["wall"] is None:
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
                ax.text(0.97, 0.96, "lab layer: no data yet (lane NEQ1)",
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="exit 1 if SCOREBOARD.md differs from regenerated text")
    args = ap.parse_args()
    rows = load_rows()
    text = write_md(rows)
    if args.check:
        if not MD.exists() or MD.read_text() != text:
            print("SCOREBOARD.md is stale; run scoreboard.py", file=sys.stderr)
            return 1
        print("SCOREBOARD.md up to date")
        return 0
    MD.write_text(text)
    draw(rows)
    print(f"wrote {MD.name} ({len(rows)} rows) and {PNG.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
