"""Plot all 24 native nested-surface diagnostic rows without smoothing.

The analytic manufactured solution and all coordinates are dimensionless.
--output selects an untracked local artifact; CSV and source are retained.
"""
import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
with Path(__file__).with_name("fci_transport_nested_accuracy.csv").open() as stream:
    rows = list(csv.DictReader(stream))

plt.rcParams.update({"font.size": 10})
fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.8), layout="constrained")
# Okabe–Ito colors and distinct markers separate refinements in grayscale.
for cells, color, marker in zip([5, 9, 17],
                               ["#0072B2", "#D55E00", "#009E73"],
                               ["o", "s", "^"]):
    values = [row for row in rows if int(row["cells"]) == cells]
    epsilon = [float(row["epsilon"]) for row in values]
    label = f"{cells} × {cells} cells; {values[0]['unknowns']} unknowns"
    for ax, key in zip(axes, ["relative_mass_error", "weak_residual"]):
        ax.plot(epsilon, [float(row[key]) for row in values],
                color=color, marker=marker, markersize=4, linewidth=1.3,
                label=label)
        ax.set_xscale("log")
        ax.set_xlabel("epsilon = chi_perp / chi_parallel (dimensionless)")
        ax.grid(True, which="major", alpha=0.25)

for ax in axes:
    ax.invert_xaxis()
axes[0].set_title("Accuracy gate fails")
axes[0].set_ylabel("Relative mass error (dimensionless)")
axes[0].set_ylim(-0.02, 1.06)
axes[0].axhline(0.1, color="0.35", linestyle="--", linewidth=1)
axes[0].text(3e-10, 0.13, "Gate: error < 0.1", color="0.3")
axes[0].legend(loc="center right", frameon=False, fontsize=9)
axes[1].set_title("Original weak equations are solved")
axes[1].set_ylabel("Relative weak residual (dimensionless)")
axes[1].set_yscale("log")
axes[1].set_ylim(1e-14, 3e-9)
axes[1].axhline(1e-9, color="0.35", linestyle="--", linewidth=1)
axes[1].text(3e-10, 1.3e-9, "Check: residual < 1e-9", color="0.3")
fig.savefig(args.output, dpi=160)
