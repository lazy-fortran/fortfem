"""Plot the failed reproduction screen and the independent residual checks."""

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_csv(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    directory = Path(__file__).resolve().parent
    paper = read_csv(directory / "fci_transport_pollution_published.csv")[0]
    checks = read_csv(directory / "fci_transport_pollution_checks.csv")
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout="constrained")

    axes[0].errorbar(0, 6, yerr=1, fmt="s", color="#0072B2", capsize=8,
                     label="Figure reading: [5, 7]")
    axes[0].plot(1, float(paper["effective_perp_error"]), "o", color="#D55E00",
                 label="Prescribed analytic field")
    axes[0].set(xticks=[0, 1], xticklabels=["Sovinec Fig. 9", "H1 reconstruction"],
                ylabel=r"$|1/T(0,0)-1|$ (dimensionless)", ylim=(0, 8),
                title=r"Bilinear, $h=1/8$, $\chi_\parallel/\chi_\perp=10^3$")
    axes[0].legend(fontsize=8, loc="lower left")

    for cells, color, marker, style in [(8, "#0072B2", "s", "-"),
                                        (16, "#D55E00", "o", "--")]:
        rows = [r for r in checks if int(r["cells"]) == cells]
        axes[1].loglog([float(r["chi_ratio"]) for r in rows],
                        [float(r["relative_residual"]) for r in rows],
                        marker=marker, linestyle=style, color=color,
                        label=f"{cells} cells per direction")
    axes[1].axhline(1e-9, color="black", linestyle=":", label="Unchanged check limit")
    axes[1].set(xlabel=r"$\chi_\parallel/\chi_\perp$ (dimensionless)",
                ylabel="Relative original residual (dimensionless)",
                title="Symmetric sine-oracle diagnostics")
    axes[1].legend(fontsize=8)
    for axis in axes:
        axis.grid(True, which="major", alpha=0.25)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
