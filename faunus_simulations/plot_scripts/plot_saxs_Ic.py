#!/usr/bin/env python3
"""Plot SAXS intensity per concentration I(c,s)/c versus s.

Reads a CSV file with the first column as s = 2*sin(theta)/lambda (Å^-1)
and one or more intensity columns I(c,s). Each selected intensity column is
normalized by the concentration c and plotted against s.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot SAXS intensity per concentration I(c,s)/c versus s."
    )
    parser.add_argument(
        "input",
        help="Path to the CSV file containing s and intensity columns.",
    )
    parser.add_argument(
        "-c",
        "--concentration",
        type=float,
        default=1.0,
        help="Concentration value to normalize the intensity (default: 1.0).",
    )
    parser.add_argument(
        "-k",
        "--columns",
        nargs="+",
        help=(
            "List of intensity column names to plot. If omitted, all numeric columns "
            "except the first s column are plotted."
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output filename for the saved figure (PNG, PDF, etc.). If omitted, the plot is shown.",
    )
    parser.add_argument(
        "--delimiter",
        default=",",
        help="Delimiter used in the input CSV file (default: ',').",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    path = Path(args.input)
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")

    df = pd.read_csv(path, delimiter=args.delimiter)
    if df.shape[1] < 2:
        raise ValueError("Input file must contain at least two columns: s and one intensity column.")

    x_col = df.columns[0]
    x = df[x_col].astype(float)

    if args.columns:
        y_columns = args.columns
    else:
        y_columns = [col for col in df.columns[1:] if pd.api.types.is_numeric_dtype(df[col])]
        if not y_columns:
            y_columns = [col for col in df.columns[1:]]

    if not y_columns:
        raise ValueError("No intensity columns found to plot.")

    fig, ax = plt.subplots(figsize=(7, 5))
    for col in y_columns:
        if col not in df.columns:
            raise ValueError(f"Selected column not found in input file: {col}")
        y = df[col].astype(float) / args.concentration
        ax.plot(x, y, marker="o", linestyle="-", label=f"{col} / {args.concentration}")

    ax.set_xlabel(r"$s = 2 \sin(\theta) / \lambda$ (Å$^{-1}$)")
    ax.set_ylabel(r"$I(c,s) / c$")
    ax.set_title(path.name)
    ax.legend(loc="best", frameon=False)
    ax.grid(True, linestyle="--", alpha=0.4)

    fig.tight_layout()
    if args.output:
        fig.savefig(args.output, dpi=300)
        print(f"Saved plot to {args.output}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
