"""Plot the computed, dimensionless errors reported by the native FCI test.

Run with Python and matplotlib; the CSV is the unfiltered printed epsilon
ladder (alpha = 0). Choose a local output path with --output.
"""
import csv
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, required=True)
parser.add_argument("--preview", type=Path)
args = parser.parse_args()
HERE = Path(__file__).resolve().parent
with (HERE / "fci_transport_reference_accuracy.csv").open() as stream:
    rows = list(csv.DictReader(stream))

# Okabe-Ito colors and distinct markers keep refinement visible in grayscale.
colors = ["#0072B2", "#D55E00", "#009E73"]
markers = ["o", "s", "^"]
plt.rcParams.update({"font.size": 10, "svg.fonttype": "none",
                     "svg.hashsalt": "fci-transport-reference-accuracy"})
fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.8), sharey=True)
for ax, fixture, title in zip(axes, ["periodic", "mixed"],
                            ["Periodic lines", "Closed and anchored open lines"]):
    counts = sorted({int(row["unknowns"]) for row in rows
                     if row["fixture"] == fixture})
    for count, color, marker in zip(counts, colors, markers):
        values = [row for row in rows
                  if row["fixture"] == fixture and int(row["unknowns"]) == count]
        ax.loglog([float(row["epsilon"]) for row in values],
                  [float(row["relative_mass_error"]) for row in values],
                  color=color, marker=marker, markersize=4, linewidth=1.4,
                  label=f"{count} unknowns")
    ax.invert_xaxis()
    ax.grid(True, which="major", alpha=0.25)
    ax.set_title(title, pad=31)
    ax.set_xlabel("epsilon (dimensionless)")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0),
              ncol=3, frameon=False, fontsize=8, columnspacing=0.8)
axes[0].set_ylabel("Relative mass error (dimensionless)")
fig.tight_layout()
fig.savefig(args.output)
if args.preview:
    fig.savefig(args.preview, dpi=160)
