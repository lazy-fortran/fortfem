"""Plot the figure report, normwise residuals and attempted scaling gates."""

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_csv(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--data-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--fit-output", type=Path)
    args = parser.parse_args()
    directory = args.data_dir
    paper = read_csv(directory / "fci_transport_pollution_published.csv")[0]
    checks = read_csv(directory / "fci_transport_pollution_checks.csv")
    fig, panels = plt.subplots(2, 2, figsize=(11, 8), layout="constrained")
    axes = panels.ravel()

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
                        [float(r["normwise_residual"]) for r in rows],
                        marker=marker, linestyle=style, color=color,
                        label=f"{cells} cells per direction")
    axes[1].axhline(float(checks[0]["residual_limit"]), color="black",
                    linestyle=":", label=r"$10^3$ machine epsilon")
    axes[1].set(xlabel=r"$\chi_\parallel/\chi_\perp$ (dimensionless)",
                ylabel="Normwise residual (dimensionless)",
                title="Symmetric sine-oracle diagnostics")
    axes[1].legend(fontsize=8)

    isotropic = read_csv(directory / "fci_transport_pollution_isotropic.csv")
    palette = ["#0072B2", "#D55E00", "#009E73", "#CC79A7"]
    for p, color, marker in zip(range(1, 5), palette, ["s", "o", "^", "D"]):
        rows = [r for r in isotropic if int(r["p"]) == p]
        h = np.array([1 / int(r["cells"]) for r in rows])
        errors = np.array([float(r["l2_error"]) for r in rows])
        axes[2].loglog(h, errors, marker=marker, color=color, label=f"p={p}")
        axes[2].loglog(h, errors[-1] * (h / h[-1]) ** (p + 1),
                       linestyle=":", color=color)
    axes[2].set(xlabel="Element width h (normalized units)",
                ylabel="L2 temperature error (normalized units)",
                title="Isotropic solution; dotted lines have order p+1")
    axes[2].legend(fontsize=8)

    paths = sorted(directory.glob("fci_transport_pollution_scaling_p*.csv"))
    if paths:
        fits = []
        for path in paths:
            rows = read_csv(path)
            design = np.array([[1, np.log(1 / int(r["cells"])),
                                np.log(float(r["chi_ratio"]))] for r in rows])
            response = np.log([float(r["signed_pollution"]) for r in rows])
            coefficients = np.linalg.lstsq(design, response, rcond=None)[0]
            fits.append(dict(p=int(rows[0]["p"]),
                             a=coefficients[1], b=coefficients[2],
                             C=np.exp(coefficients[0]),
                             log_rms=np.sqrt(np.mean((design @ coefficients - response)**2))))
        if args.fit_output:
            with args.fit_output.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=["p", "a", "b", "C", "log_rms"],
                                        lineterminator="\n")
                writer.writeheader()
                writer.writerows(fits)
        # Display all rows of the last attempted degree, including a failed gate.
        scaling = read_csv(paths[-1])
        p = int(scaling[0]["p"])
        ratios = sorted({float(r["chi_ratio"]) for r in scaling})
        for ratio, color, marker in zip(ratios, palette, ["s", "o", "^", "D"]):
            rows = [r for r in scaling if float(r["chi_ratio"]) == ratio]
            h = np.array([1 / int(r["cells"]) for r in rows])
            errors = np.array([float(r["signed_pollution"]) / ratio for r in rows])
            axes[3].loglog(h, errors, marker=marker, color=color,
                           label=rf"$\chi_\parallel/\chi_\perp={ratio:.0e}$")
        axes[3].loglog(h, errors[-1] * (h / h[-1]) ** (2 * p), color="black",
                       linestyle=":", label=f"Expected order {2*p}")
        axes[3].set(xlabel="Element width h (normalized units)",
                    ylabel=r"$(\chi_{\perp,\rm num}/\chi_\perp)/(\chi_\parallel/\chi_\perp)$"
                           "\n(dimensionless)",
                    title=f"Conforming quadrilaterals, last attempted p={p}")
        axes[3].legend(fontsize=8, loc="lower right")
    else:
        axes[3].set_axis_off()
    for axis in axes:
        axis.grid(True, which="major", alpha=0.25)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
