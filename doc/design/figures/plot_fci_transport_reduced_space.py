"""Plot unfiltered reduced-space measurements and coarse SVD kernel vectors.

All coordinates, norms and singular values use dimensionless fixture units.
Kernel vectors have unit Euclidean input norm before orthogonal projection.
"""
import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator
import numpy as np

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
root = Path(__file__).parent
with (root / "fci_transport_reduced_space.csv").open() as stream:
    rows = list(csv.DictReader(stream))

plt.rcParams.update({"font.size": 10})
fig, axes = plt.subplots(2, 2, figsize=(9.2, 7.6), layout="constrained")
for name, color, marker, line in zip(
        ["P0-2h", "Q1-2h"], ["#0072B2", "#D55E00"], ["o", "s"], ["-", "--"]):
    values = [row for row in rows if row["space"] == name]
    cells = [int(row["cells"]) for row in values]
    axes[0, 0].loglog(cells, [float(row["beta_h1"]) for row in values],
                      color=color, marker=marker, linestyle=line, label=name)
    axes[0, 1].loglog(cells, [int(row["kernel_dim"]) for row in values],
                      color=color, marker=marker, linestyle=line, label=name)

axes[0, 0].set_title("H1 / L2 inf-sup constant")
axes[0, 0].set_ylabel("Generalized singular value (dimensionless)")
axes[0, 1].set_title("Kernel and radial resolution")
axes[0, 1].set_ylabel("Number of independent modes (count)")
values = [row for row in rows if row["space"] == "P0-2h"]
axes[0, 1].loglog([int(row["cells"]) for row in values],
                  [int(row["surface_bins"]) for row in values],
                  color="0.3", marker="^", linestyle=":", label="Radial bins")
for ax in axes[0]:
    ax.set_xlabel("Section cells per direction (count)")
    ax.set_xticks([21, 41, 81], labels=["21", "41", "81"])
    ax.xaxis.set_minor_locator(NullLocator())
    ax.grid(True, which="major", alpha=0.25)
    ax.legend(frameon=False)

images = []
for name in ["P0-2h", "Q1-2h"]:
    with (root / f"fci_transport_{name}_kernel.csv").open() as stream:
        data = list(csv.DictReader(stream))
    side = int(np.sqrt(len(data)))
    images.append(np.array([float(row["odd_plane_kernel"]) for row in data])
                  .reshape(side, side))
peak = max(np.max(np.abs(field)) for field in images)
for ax, name, field in zip(axes[1], ["P0-2h", "Q1-2h"], images):
    h = 1 / (field.shape[0] + 1)
    image = ax.imshow(field, origin="lower", extent=[h / 2, 1 - h / 2] * 2,
                      cmap="gray_r", vmin=-peak, vmax=peak, interpolation="nearest")
    value = next(row for row in rows if row["space"] == name)
    ax.set_title(f"{name}: {100 * float(value['checkerboard_retained']):.2f}% "
                 "retained (21 cells)")
    ax.set_xlabel("x (dimensionless)")
    ax.set_ylabel("y (dimensionless)")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.5, -0.21, "Opposite sign on the second plane", ha="center",
            transform=ax.transAxes, fontsize=9)
fig.colorbar(image, ax=list(axes[1]), shrink=0.85,
             label="Projected vector amplitude (dimensionless)")
fig.savefig(args.output, dpi=160)
