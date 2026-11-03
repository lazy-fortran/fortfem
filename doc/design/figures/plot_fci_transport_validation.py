"""Collect registered validation CSVs, plot every mesh, and tabulate thresholds."""

import argparse
import csv
from pathlib import Path
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from plot_fci_transport_cost import STYLES, group_rows


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def write(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def plot_manufactured(groups, output):
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), layout="constrained", sharey="row")
    for column, ratio in enumerate([1e6, 1e8, 1e10]):
        for (method, p), (label, color, marker, style) in STYLES.items():
            rows = groups["manufactured", method, p, ratio]
            n = [int(r["unknowns"]) for r in rows]
            for axis, quantity in zip(axes[:, column], ["l2_error", "centre_error"]):
                axis.loglog(n, [float(r[quantity]) for r in rows], label=label,
                            color=color, marker=marker, linestyle=style)
        axes[0, column].set_title(rf"$R={ratio:.0e}$")
        axes[0, column].axhline(0.01, color="0.5", linestyle=":")
        axes[1, column].set_xlabel("Free unknowns (count)")
        for axis in axes[:, column]:
            axis.grid(True, which="major", alpha=0.25)
    axes[0, 0].set_ylabel("Absolute L2 temperature error\n(normalized units)")
    axes[1, 0].set_ylabel("Absolute T(0,0) error\n(normalized units)")
    fig.legend(*axes[0, 0].get_legend_handles_labels(), loc="outside lower center", ncol=3)
    fig.suptitle("Manufactured temperature varies along the field; all sampled meshes")
    fig.savefig(output / "manufactured.png", dpi=180)
    plt.close(fig)


def plot_flux(groups, output):
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), layout="constrained", sharey=True)
    for row, problem in enumerate(["sovinec", "circle"]):
        for column, ratio in enumerate([1e6, 1e8, 1e10]):
            axis = axes[row, column]
            for method in ["fci", "symmetric"]:
                label, color, marker, style = STYLES[method, 1]
                rows = groups[problem, method, 1, ratio]
                axis.plot([int(r["cells"]) for r in rows],
                          [float(r["signed_pollution"]) for r in rows],
                          label=label, color=color, marker=marker, linestyle=style)
            axis.set(xscale="log", yscale="symlog", title=rf"{problem}, $R={ratio:.0e}$")
            axis.set_yscale("symlog", linthresh=1e-8)
            axis.axhline(0.1, color="0.5", linestyle=":")
            axis.axhline(0, color="0.5", linewidth=0.5)
            axis.grid(True, which="major", alpha=0.25)
            axis.set_xlabel("Cells per direction (count)")
        axes[row, 0].set_ylabel("Signed pollution E (dimensionless)")
    fig.legend(*axes[0, 0].get_legend_handles_labels(), loc="outside lower center", ncol=2)
    fig.suptitle("Flux-function tests: every signed measurement, no absolute value or fit")
    fig.savefig(output / "flux-refinement.png", dpi=180)
    plt.close(fig)


def thresholds(groups, cost, path):
    rows = []
    for (method, p), (label, _, _, _) in STYLES.items():
        for ratio in [1e6, 1e8, 1e10]:
            result = dict(candidate=label, chi_ratio=f"{ratio:.0e}")
            for problem in ["sovinec", "circle", "manufactured"]:
                measurements = groups.get((problem, method, p, ratio), [])
                if problem != "manufactured":
                    measurements = cost.get((problem, method, p, ratio), measurements)
                if not measurements:
                    result.update({f"{problem}_unknowns": "not measured",
                                   f"{problem}_error": "not measured",
                                   f"{problem}_time_s": "not measured"})
                    continue
                quantity = "l2_error" if problem == "manufactured" else "signed_pollution"
                limit = 0.01 if problem == "manufactured" else 0.1
                reached = [r for r in measurements if float(r[quantity]) <= limit]
                chosen = reached[0] if reached else measurements[-1]
                count = chosen["unknowns"]
                if not reached:
                    count = f"not reached at {count}"
                result.update({f"{problem}_unknowns": count,
                               f"{problem}_error": chosen[quantity].strip(),
                               f"{problem}_time_s": float(chosen["setup_s"]) +
                               float(chosen["solve_s"])})
            rows.append(result)
    write(path, rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path)
    parser.add_argument("--data", type=Path, default=Path(__file__).with_name(
        "fci_transport_pollution_validation.csv"))
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--summary-output", type=Path, required=True)
    args = parser.parse_args()
    if args.build_dir:
        paths = sorted(args.build_dir.glob("pollution_refine_*.csv"))
        paths += sorted(args.build_dir.glob("pollution_cap_*.csv"))
        if len(paths) != 34:
            raise ValueError(f"Expected 30 refinement and 4 cap CSVs, found {len(paths)}")
        write(args.data, [row for path in paths for row in read(path)])
        shutil.copyfile(args.build_dir / "pollution_candidate_parallel.csv",
                        args.data.with_name("fci_transport_pollution_candidate_parallel.csv"))
    # Independent quadrature establishes that a flux-only profile cannot pass.
    for order in [16, 32, 64, 128]:
        x, w = np.polynomial.legendre.leggauss(order)
        xx, yy = np.meshgrid(x / 2, x / 2, indexing="ij")
        psi = np.cos(np.pi * xx) * np.cos(np.pi * yy)
        field_square = (np.cos(np.pi * xx) * np.sin(np.pi * yy)) ** 2 + (
            np.sin(np.pi * xx) * np.cos(np.pi * yy)) ** 2
        perturbation = 4 * xx * psi * field_square ** 2
        norm = np.sqrt(np.sum(w[:, None] * w[None, :] * perturbation ** 2) / 4)
        print(f"Independent perturbation L2, Gauss order {order}: {norm:.15g}")
    groups = group_rows(read(args.data))
    cost = group_rows(read(Path(__file__).with_name("fci_transport_pollution_cost.csv")))
    thresholds(groups, cost, args.summary_output)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    plot_manufactured(groups, args.output_dir)
    plot_flux(groups, args.output_dir)
    parallel = read(args.data.with_name("fci_transport_pollution_candidate_parallel.csv"))
    fig, axis = plt.subplots(figsize=(7, 4), layout="constrained")
    for method in ["fci", "symmetric"]:
        label, color, marker, style = STYLES[method, 1]
        rows = [r for r in parallel if r["method"] == method]
        axis.loglog([int(r["cells"]) for r in rows],
                    [float(r["l2_parallel_error"]) for r in rows],
                    label=label, color=color, marker=marker, linestyle=style)
    axis.set(xlabel="Cells per direction (count)",
             ylabel="L2 parallel-gradient error (normalized units)",
             title="Nonzero analytic gradient: exp(x+2y)")
    axis.grid(True, which="major", alpha=0.25)
    axis.legend()
    fig.savefig(args.output_dir / "parallel.png", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
