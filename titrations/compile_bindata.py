#!/usr/bin/env python3
"""Compile titration BinData readings into a CSV plus summary plots."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable, List

import matplotlib

matplotlib.use("Agg")  # ensure plots render without a display
import matplotlib.pyplot as plt


@dataclass
class Measurement:
    """Structured representation of a single reading."""

    dsec: float
    timestamp: datetime
    temperature_c: float
    total_volume: float
    concentration: float
    ph: float


def parse_measurement_file(path: Path) -> Measurement:
    """Parse a `.mxf` file into a Measurement instance."""

    text = path.read_text(encoding="utf-8").strip()
    if not text:
        raise ValueError(f"File {path} is empty")

    fields = {}
    normalized = text.replace("\n", "\t")
    for chunk in normalized.split("\t"):
        chunk = chunk.strip()
        if not chunk:
            continue
        if "=" not in chunk:
            continue
        key, value = chunk.split("=", 1)
        fields[key.upper()] = value

    try:
        dsec = float(fields["DSEC"])
        timestamp = datetime.strptime(fields["TIME"], "%Y-%m-%d %H:%M:%S")
        temperature_c = float(fields["TEMP"])
        total_volume = float(fields["VOL"])
        concentration = float(fields["CONC"])
        ph = float(fields["ELE"])
    except KeyError as exc:
        raise ValueError(f"Missing field {exc} in {path}") from exc

    return Measurement(
        dsec=dsec,
        timestamp=timestamp,
        temperature_c=temperature_c,
        total_volume=total_volume,
        concentration=concentration,
        ph=ph,
    )


def build_rows(measurements: List[Measurement]) -> List[dict]:
    """Convert measurements into ordered rows for the output CSV."""

    rows: List[dict] = []
    prev_time = None
    prev_volume = None

    for measurement in measurements:
        time_seconds = measurement.dsec / 10.0
        if prev_time is None:
            delta_t = 0.0
        else:
            delta_t = (measurement.timestamp - prev_time).total_seconds()

        if prev_volume is None:
            delta_volume = 0.0
        else:
            delta_volume = measurement.total_volume - prev_volume

        rows.append(
            {
                "time_seconds": round(time_seconds, 6),
                "time_between_measurements": round(delta_t, 6),
                "temperature_c": round(measurement.temperature_c, 6),
                "total_volume": round(measurement.total_volume, 6),
                "volume_delta": round(delta_volume, 6),
                "concentration": round(measurement.concentration, 6),
                "ph": round(measurement.ph, 6),
            }
        )

        prev_time = measurement.timestamp
        prev_volume = measurement.total_volume

    return rows


def write_csv(rows: Iterable[dict], csv_path: Path) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "time_seconds",
        "time_between_measurements",
        "temperature_c",
        "total_volume",
        "volume_delta",
        "concentration",
        "ph",
    ]

    with csv_path.open("w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def plot_relationships(rows: List[dict], parent_name: str, output_dir: Path) -> None:
    ph_values = [row["ph"] for row in rows]
    plot_specs = [
        ("time_seconds", "Time Elapsed (s)"),
        ("time_between_measurements", "Delta Time (s)"),
        ("temperature_c", "Temperature (C)"),
        ("total_volume", "Total Volume"),
        ("volume_delta", "Delta Volume"),
        ("concentration", "Concentration"),
    ]
    invert_columns = {"time_seconds", "total_volume", "concentration"}

    output_dir.mkdir(parents=True, exist_ok=True)
    for column, xlabel in plot_specs:
        metric_values = [row[column] for row in rows]
        plt.figure(figsize=(6, 4))
        if column in invert_columns:
            x_values = metric_values
            y_values = ph_values
            plt.plot(x_values, y_values, marker="o", linewidth=1)
            plt.xlabel(xlabel)
            plt.ylabel("pH (ELE)")
            plt.title(f"{parent_name} pH vs {xlabel}")
            plot_filename = f"{parent_name}_pH_vs_{column}.png"
        else:
            x_values = ph_values
            y_values = metric_values
            plt.plot(x_values, y_values, marker="o", linewidth=1)
            plt.xlabel("pH (ELE)")
            plt.ylabel(xlabel)
            plt.title(f"{parent_name} {xlabel} vs pH")
            plot_filename = f"{parent_name}_{column}_vs_pH.png"
        plt.grid(True, linestyle="--", linewidth=0.5, alpha=0.6)
        plt.tight_layout()
        plot_path = output_dir / plot_filename
        plt.savefig(plot_path, dpi=200)
        plt.close()


def determine_output_paths(
    data_dir: Path, first_timestamp: datetime, output_dir: Path | None
) -> tuple[Path, Path, str]:
    parent_dir = data_dir.resolve().parent
    parent_name = parent_dir.name
    date_fragment = first_timestamp.strftime("%Y%m%d")
    base_name = f"{parent_name}_{date_fragment}"

    target_dir = Path(output_dir) if output_dir else parent_dir
    csv_path = target_dir / f"{base_name}.csv"
    plots_dir = target_dir / "plots"
    return csv_path, plots_dir, parent_name


def main() -> None:
    parser = argparse.ArgumentParser(description="Compile BinData .mxf files into CSV + plots")
    parser.add_argument(
        "data_dir",
        nargs="?",
        default="Blank_4/BinData",
        help="Directory containing .mxf files (default: Blank_4/BinData)",
    )
    parser.add_argument(
        "--output-dir",
        dest="output_dir",
        help="Directory where the CSV and plots should be written (default: parent of data_dir)",
    )
    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    if not data_dir.exists():
        raise SystemExit(f"Data directory {data_dir} does not exist")

    files = sorted(data_dir.glob("*.mxf"))
    if not files:
        raise SystemExit(f"No .mxf files found in {data_dir}")

    measurements = [parse_measurement_file(path) for path in files]
    first_file_timestamp = measurements[0].timestamp
    measurements.sort(key=lambda m: m.dsec)

    rows = build_rows(measurements)
    csv_path, plots_dir, parent_name = determine_output_paths(
        data_dir, first_file_timestamp, args.output_dir
    )
    write_csv(rows, csv_path)
    plot_relationships(rows, parent_name, plots_dir)

    print(f"Wrote {csv_path}")
    print(f"Plots saved to {plots_dir}")


if __name__ == "__main__":
    main()
