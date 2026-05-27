"""
Universal faunus output plotter.

Detects and plots all *.csv.gz (energy-type) and *.dat.gz (RDF-type) files
found in the simulation directory.

Usage:
    python3 plot_faunus.py [--dir PATH] [--window N] [--save]
"""

import argparse
import gzip
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ── helpers ───────────────────────────────────────────────────────────────────

def running_mean(x, w):
    """Centered running mean; returns (trimmed_indices, smoothed_values).
    Uses mode='valid' so no zero-padding artifact at the edges."""
    kernel = np.ones(w) / w
    smoothed = np.convolve(x, kernel, mode="valid")
    # valid output is shorter: skip w//2 points at each end
    trim = w // 2
    indices = np.arange(trim, len(x) - (w - 1 - trim))
    return indices, smoothed


def load_csv_gz(path):
    """Return (header_list, data_array) from a gzipped CSV with a header row."""
    with gzip.open(path, "rt") as f:
        header = f.readline().strip().split(",")
        data = np.loadtxt(f, delimiter=",")
    if data.ndim == 1:
        data = data[np.newaxis, :]
    return header, data


def load_dat_gz(path):
    """Return (col_names, data_array) from a gzipped space-separated file
    whose first line is a comment like '# r g(r)'."""
    with gzip.open(path, "rt") as f:
        first = f.readline().strip()
        if first.startswith("#"):
            col_names = first.lstrip("# ").split()
        else:
            col_names = [f"col{i}" for i in range(len(first.split()))]
        data = np.loadtxt(f)
    if data.ndim == 1:
        data = data[np.newaxis, :]
    return col_names, data


def find_files(directory, pattern):
    return sorted(Path(directory).glob(pattern))


# ── plot helpers ──────────────────────────────────────────────────────────────

COLORS = plt.rcParams["axes.prop_cycle"].by_key()["color"]


def plot_energy(ax, path, window, color_offset=0):
    """Plot all energy columns (except 'step' and 'average') vs step."""
    header, data = load_csv_gz(path)
    if "step" not in header:
        ax.set_title(f"{path.name}\n(no 'step' column)", fontsize=8)
        return

    step_idx = header.index("step")
    steps = data[:, step_idx]
    plotted = 0
    for i, col in enumerate(header):
        if col in ("step", "average"):
            continue
        y = data[:, i]
        c = COLORS[(color_offset + plotted) % len(COLORS)]
        ax.plot(steps, y, lw=0.6, alpha=0.4, color=c)
        if len(y) >= window:
            idx, smoothed = running_mean(y, window)
            ax.plot(steps[idx], smoothed, lw=1.8, color=c, label=col)
        else:
            ax.plot(steps, y, lw=1.8, color=c, label=col)
        plotted += 1

    ax.set_xlabel("MC step")
    ax.set_ylabel("Energy (kJ mol⁻¹)")
    ax.set_title(path.name, fontsize=9)
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)


def plot_rdf(ax, path, color_offset=0):
    """Plot g(r) vs r."""
    col_names, data = load_dat_gz(path)
    if data.size == 0 or data.shape[1] < 2:
        ax.set_title(f"{path.name}\n(empty)", fontsize=8)
        return

    r = data[:, 0]
    gr = data[:, 1]
    c = COLORS[color_offset % len(COLORS)]
    ax.plot(r, gr, lw=1.8, color=c)
    ax.axhline(1.0, color="gray", lw=0.8, ls="--", alpha=0.6)
    ax.set_xlabel("r (Å)")
    ax.set_ylabel("g(r)")
    ax.set_title(path.name, fontsize=9)
    ax.grid(True, alpha=0.3)


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Plot faunus simulation output.")
    parser.add_argument("--dir", default=".", metavar="PATH",
                        help="Simulation directory (default: current)")
    parser.add_argument("--window", type=int, default=50, metavar="N",
                        help="Running-mean window for energy plots (default: 50)")
    parser.add_argument("--save", action="store_true",
                        help="Save figure as PNG instead of showing interactively")
    args = parser.parse_args()

    sim_dir = Path(args.dir).resolve()
    if not sim_dir.is_dir():
        sys.exit(f"Error: '{sim_dir}' is not a directory.")

    csv_files = find_files(sim_dir, "*.csv.gz")
    dat_files = find_files(sim_dir, "*.dat.gz")

    if not csv_files and not dat_files:
        sys.exit(f"No *.csv.gz or *.dat.gz files found in '{sim_dir}'.")

    n_csv = len(csv_files)
    n_dat = len(dat_files)
    n_cols = max(n_csv, n_dat, 1)

    fig = plt.figure(figsize=(5 * n_cols, 4 * (bool(n_csv) + bool(n_dat))))
    gs = gridspec.GridSpec(
        bool(n_csv) + bool(n_dat), n_cols,
        hspace=0.45, wspace=0.35
    )

    # Energy rows
    row = 0
    if n_csv:
        for j, path in enumerate(csv_files):
            ax = fig.add_subplot(gs[row, j])
            plot_energy(ax, path, args.window, color_offset=j)
        row += 1

    # RDF rows
    if n_dat:
        for j, path in enumerate(dat_files):
            ax = fig.add_subplot(gs[row, j])
            plot_rdf(ax, path, color_offset=j)

    fig.suptitle(str(sim_dir), fontsize=10, y=1.01)

    if args.save:
        out = sim_dir / "faunus_observables.png"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"Saved: {out}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
