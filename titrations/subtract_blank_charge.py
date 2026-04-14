#!/usr/bin/env python3
"""Subtract blank titration charges from sample titration to get net charge."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import List, Tuple

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_csv_data(csv_path: Path) -> Tuple[List[float], List[float]]:
    """Read pH and net_charge from a CSV file."""
    ph_values: List[float] = []
    net_charge_values: List[float] = []

    with csv_path.open("r", encoding="utf-8", newline="") as fp:
        reader = csv.DictReader(fp)
        if "ph" not in reader.fieldnames or "net_charge" not in reader.fieldnames:
            raise ValueError(
                f"CSV file {csv_path} must contain 'ph' and 'net_charge' columns"
            )

        for row in reader:
            try:
                ph_values.append(float(row["ph"]))
                net_charge_values.append(float(row["net_charge"]))
            except (TypeError, ValueError) as exc:
                raise ValueError(f"Invalid numeric value in {csv_path}: {row}") from exc

    if not ph_values:
        raise ValueError(f"CSV file {csv_path} contains no data rows")

    return ph_values, net_charge_values


def interpolate_charge(
    blank_ph: List[float],
    blank_charge: List[float],
    sample_ph: List[float],
) -> List[float]:
    """Interpolate blank charge values at sample pH points."""
    interpolated = []
    for ph in sample_ph:
        if ph <= blank_ph[0]:
            interpolated.append(blank_charge[0])
        elif ph >= blank_ph[-1]:
            interpolated.append(blank_charge[-1])
        else:
            idx = 0
            while idx < len(blank_ph) - 1 and blank_ph[idx + 1] < ph:
                idx += 1

            ph_low = blank_ph[idx]
            ph_high = blank_ph[idx + 1]
            charge_low = blank_charge[idx]
            charge_high = blank_charge[idx + 1]

            fraction = (ph - ph_low) / (ph_high - ph_low)
            interpolated_charge = charge_low + fraction * (charge_high - charge_low)
            interpolated.append(interpolated_charge)

    return interpolated


def write_net_charge_csv(
    sample_ph: List[float],
    raw_diff: List[float],
    normalized_charge: List[float],
    output_csv: Path,
) -> None:
    """Write the raw and normalized net charge data to CSV."""
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(
            fp,
            fieldnames=["ph", "raw_charge_diff", "net_charge"],
        )
        writer.writeheader()
        for ph, raw, normalized in zip(sample_ph, raw_diff, normalized_charge):
            writer.writerow(
                {
                    "ph": ph,
                    "raw_charge_diff": raw,
                    "net_charge": normalized,
                }
            )


def plot_net_charge(
    sample_ph: List[float],
    net_charges: List[float],
    blank_ph: List[float],
    blank_charges: List[float],
    sample_charges: List[float],
    output_plot: Path,
) -> None:
    """Plot the net charge, blank, and sample data."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Plot all three curves
    ax1.plot(blank_ph, blank_charges, marker="o", linestyle="-", label="Blank (TIT1)", markersize=3, alpha=0.7)
    ax1.plot(sample_ph, sample_charges, marker="s", linestyle="-", label="Sample (TITPO)", markersize=3, alpha=0.7)
    ax1.set_xlabel("pH")
    ax1.set_ylabel("Net Charge (mol)")
    ax1.set_title("Blank and Sample Titration Curves")
    ax1.grid(True, linestyle="--", alpha=0.6)
    ax1.legend()

    # Plot net charge
    ax2.plot(sample_ph, net_charges, marker="D", linestyle="-", color="red", markersize=3)
    ax2.axhline(y=0, color="k", linestyle="--", alpha=0.3)
    ax2.set_xlabel("pH")
    ax2.set_ylabel("Normalized Net Charge")
    ax2.set_title("Normalized Net Charge (Sample - Blank)")
    ax2.grid(True, linestyle="--", alpha=0.6)

    plt.tight_layout()
    output_plot.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_plot, dpi=200)
    plt.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Subtract blank titration charges from sample to get net charge"
    )
    parser.add_argument(
        "--blank",
        required=True,
        help="CSV file or directory containing blank titration data (TIT1)",
    )
    parser.add_argument(
        "--sample",
        required=True,
        help="CSV file or directory containing sample titration data (TITPO)",
    )
    parser.add_argument(
        "--output-csv",
        default="net_charge.csv",
        help="Output CSV filename (default: net_charge.csv)",
    )
    parser.add_argument(
        "--output-plot",
        default="net_charge_plot.png",
        help="Output plot filename (default: net_charge_plot.png)",
    )
    parser.add_argument(
        "--sample-moles",
        type=float,
        default=1.99e-5,
        help="Total number of sample moles used for normalization (default: 1.99e-5)",
    )
    parser.add_argument(
        "--extra-charge",
        type=float,
        default=0.0,
        help="Extra charge term to add after normalization (default: 0.0)",
    )
    return parser.parse_args()


def find_csv_file(path: Path) -> Path:
    """Find a single CSV file for a directory or return a CSV path directly."""
    if path.is_file() and path.suffix.lower() == ".csv":
        return path

    if not path.is_dir():
        raise FileNotFoundError(f"Path does not exist or is not a directory: {path}")

    candidates = sorted(path.glob("*.csv"))
    if not candidates:
        raise FileNotFoundError(f"No CSV files found in directory: {path}")

    if len(candidates) == 1:
        return candidates[0]

    prefix_matches = [candidate for candidate in candidates if candidate.stem.startswith(path.name)]
    if len(prefix_matches) == 1:
        return prefix_matches[0]
    if prefix_matches:
        return sorted(prefix_matches, key=lambda p: p.stat().st_mtime)[-1]

    return sorted(candidates, key=lambda p: p.stat().st_mtime)[-1]


def main() -> None:
    args = parse_args()

    # Find CSV files
    blank_path = Path(args.blank)
    sample_path = Path(args.sample)

    blank_csv = find_csv_file(blank_path)
    sample_csv = find_csv_file(sample_path)

    print(f"Reading blank data from {blank_csv}")
    blank_ph, blank_charge = read_csv_data(blank_csv)

    print(f"Reading sample data from {sample_csv}")
    sample_ph, sample_charge = read_csv_data(sample_csv)

    # Find overlapping pH range
    min_ph = max(min(blank_ph), min(sample_ph))
    max_ph = min(max(blank_ph), max(sample_ph))
    print(f"Overlapping pH range: {min_ph:.2f} - {max_ph:.2f}")

    # Filter data to overlapping range
    blank_ph_filtered = [p for p in blank_ph if min_ph <= p <= max_ph]
    blank_charge_filtered = [c for p, c in zip(blank_ph, blank_charge) if min_ph <= p <= max_ph]
    sample_ph_filtered = [p for p in sample_ph if min_ph <= p <= max_ph]
    sample_charge_filtered = [c for p, c in zip(sample_ph, sample_charge) if min_ph <= p <= max_ph]

    # Sort filtered data by pH (in case order changed)
    blank_sorted = sorted(zip(blank_ph_filtered, blank_charge_filtered), key=lambda x: x[0])
    blank_ph_filtered, blank_charge_filtered = zip(*blank_sorted) if blank_sorted else ([], [])
    sample_sorted = sorted(zip(sample_ph_filtered, sample_charge_filtered), key=lambda x: x[0])
    sample_ph_filtered, sample_charge_filtered = zip(*sample_sorted) if sample_sorted else ([], [])

    if not sample_ph_filtered:
        raise ValueError("No overlapping pH data between blank and sample")

    print(f"Interpolating blank charges to sample pH points...")
    interpolated_blank_charge = interpolate_charge(blank_ph_filtered, blank_charge_filtered, sample_ph_filtered)

    print(f"Calculating raw charge difference: sample - blank")
    raw_charge_diff = [
        sample_c - blank_c for sample_c, blank_c in zip(sample_charge_filtered, interpolated_blank_charge)
    ]

    sample_moles = float(args.sample_moles)
    extra_charge = float(args.extra_charge)
    print(f"Normalizing by sample moles: {sample_moles:.6g}, extra charge term: {extra_charge:.6g}")
    normalized_net_charge = [
        diff / sample_moles + extra_charge for diff in raw_charge_diff
    ]

    output_csv = Path(args.output_csv)
    output_plot = Path(args.output_plot)
    write_net_charge_csv(sample_ph_filtered, raw_charge_diff, normalized_net_charge, output_csv)
    
    plot_net_charge(
        sample_ph_filtered,
        normalized_net_charge,
        blank_ph_filtered,
        blank_charge_filtered,
        sample_charge_filtered,
        output_plot,
    )

    print(f"Wrote net charge data to {output_csv}")
    print(f"Wrote plot to {output_plot}")
    print(f"Processed {len(sample_ph_filtered)} data points")
    print(f"pH range: {min(sample_ph_filtered):.2f} - {max(sample_ph_filtered):.2f}")


if __name__ == "__main__":
    main()
