"""Plot every direct pollution measurement and write the sampled-cost table."""

import argparse
import csv
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

STYLES = {
    ("h1", 1): ("H1 p=1", "#0072B2", "o", "-"),
    ("h1", 2): ("H1 p=2", "#D55E00", "s", "--"),
    ("h1", 3): ("H1 p=3", "#009E73", "^", "-."),
    ("h1", 4): ("H1 p=4", "#CC79A7", "D", ":"),
    ("symmetric", 1): ("Gunter symmetric", "#E69F00", "v", "--"),
    ("fci", 1): ("Bilinear FCI", "#000000", "P", "-"),
}


def group_rows(rows):
    groups = {}
    for row in rows:
        key = (row["benchmark"], row["method"], int(row["p"]),
               float(row["chi_ratio"]))
        groups.setdefault(key, []).append(row)
    for group in groups.values():
        group.sort(key=lambda row: int(row["unknowns"]))
    return groups


def draw(groups, problem, output):
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), layout="constrained",
                             sharey="row")
    for column, ratio in enumerate([1e6, 1e8, 1e10]):
        for (method, degree), (label, color, marker, style) in STYLES.items():
            rows = groups.get((problem, method, degree, ratio), [])
            if not rows:
                continue
            n = [int(row["unknowns"]) for row in rows]
            error = [float(row["signed_pollution"]) for row in rows]
            time = [float(row["setup_s"]) + float(row["solve_s"]) for row in rows]
            axes[0, column].plot(n, error, marker=marker, color=color,
                                 linestyle=style, label=label)
            axes[1, column].loglog(n, time, marker=marker, color=color,
                                   linestyle=style, label=label)
        axes[0, column].set(xscale="log", yscale="symlog",
                             title=rf"$R={ratio:.0e}$")
        axes[0, column].set_yscale("symlog", linthresh=1e-8)
        axes[0, column].axhline(0.1, color="#777777", linestyle=":")
        axes[0, column].axhline(0, color="#777777", linewidth=0.5)
        axes[1, column].set_xlabel("Free unknowns (count)")
        for row in range(2):
            axes[row, column].grid(True, which="major", alpha=0.25)
    axes[0, 0].set_yticks([-1e-2, 0, 1e-6, 1e-2, 0.1, 1e2, 1e6])
    axes[0, 0].set_ylabel("Signed pollution E (dimensionless)")
    axes[1, 0].set_ylabel("Assembly + factorization + solve (s)")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=3, fontsize=9)
    fig.suptitle(f"{problem.capitalize()}: measured meshes; dotted threshold E=0.1")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def draw_circle(groups, output):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout="constrained")
    for method in ["symmetric", "fci"]:
        label, color, marker, style = STYLES[method, 1]
        rows = [groups["circle", method, 1, ratio][0] for ratio in [1e6, 1e8, 1e10]]
        ratio = [float(row["chi_ratio"]) for row in rows]
        axes[0].plot(ratio, [float(row["signed_pollution"]) for row in rows],
                      label=label, color=color, marker=marker, linestyle=style)
        axes[1].loglog(ratio, [float(row["setup_s"]) + float(row["solve_s"])
                               for row in rows], label=label, color=color,
                        marker=marker, linestyle=style)
    axes[0].set_xscale("log")
    axes[0].set_yscale("symlog", linthresh=1e-8)
    axes[0].axhline(0.1, color="#777777", linestyle=":", label="Threshold E=0.1")
    axes[0].axhline(0, color="#777777", linewidth=0.5)
    axes[0].set_ylabel("Signed pollution E (dimensionless)")
    axes[1].set_ylabel("Assembly + factorization + solve (s)")
    for axis in axes:
        axis.set_xlabel("Anisotropy R (dimensionless)")
        axis.grid(True, which="major", alpha=0.25)
    fig.suptitle("Circular problem: both candidates stop at 9 free unknowns")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=3, fontsize=9)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=Path(__file__).with_name(
        "fci_transport_pollution_cost.csv"))
    parser.add_argument("--build-dir", type=Path,
                        help="Collect the registered cost/circle and parallel CSVs")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--circle-output", type=Path, required=True)
    parser.add_argument("--summary-output", type=Path, required=True)
    args = parser.parse_args()
    if args.build_dir:
        paths = [args.build_dir / f"pollution_cost_{method}_p{p}_r{power}.csv"
                 for method, degrees in [("h1", range(1, 5)),
                                         ("symmetric", [1]), ("fci", [1])]
                 for p in degrees for power in [6, 8, 10]]
        paths += [args.build_dir / f"pollution_circle_{method}_r{power}.csv"
                  for method in ["symmetric", "fci"] for power in [6, 8, 10]]
        rows = []
        for path in paths:
            with path.open() as stream:
                reader = csv.DictReader(stream)
                fields = reader.fieldnames
                rows.extend(reader)
        with args.data.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        shutil.copyfile(args.build_dir / "pollution_parallel.csv",
                        args.data.with_name("fci_transport_pollution_parallel.csv"))
    with args.data.open() as stream:
        groups = group_rows(list(csv.DictReader(stream)))
    fields = ["benchmark", "candidate", "chi_ratio", "unknowns", "signed_pollution",
              "setup_s", "solve_s", "total_s", "status", "solver"]
    with args.summary_output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for (problem, method, degree, ratio), rows in sorted(groups.items()):
            reached = [row for row in rows if float(row["signed_pollution"]) <= 0.1]
            row = reached[0] if reached else rows[-1]
            writer.writerow(dict(
                benchmark=problem, candidate=STYLES[method, degree][0],
                chi_ratio=f"{ratio:.0e}", unknowns=row["unknowns"],
                signed_pollution=row["signed_pollution"].strip(),
                setup_s=row["setup_s"].strip(), solve_s=row["solve_s"].strip(),
                total_s=float(row["setup_s"]) + float(row["solve_s"]),
                status="reached" if reached else "not reached at cap",
                solver="FortSparse native SPD"))
    draw(groups, "sovinec", args.output)
    draw_circle(groups, args.circle_output)


if __name__ == "__main__":
    main()
