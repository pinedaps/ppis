#!/usr/bin/env python3
"""Compile titration BinData readings into a CSV plus summary plots."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable, List

from pint import UnitRegistry

import matplotlib

matplotlib.use("Agg")  # ensure plots render without a display
import matplotlib.pyplot as plt

ureg = UnitRegistry()


def microliters_to_liters(value: float) -> float:
    """Convert a microliter-valued float into liters using pint."""

    quantity = value * ureg.microliter
    return quantity.to(ureg.liter).magnitude


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


def parse_td_file(path: Path) -> List[Measurement]:
    """Parse a folder-specific TD.mxt file into a list of Measurement instances."""
    text = path.read_text(encoding="utf-8", errors="replace")
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    header_index = None
    header_columns: list[str] = []
    for index, line in enumerate(lines):
        if line.upper().startswith("X\tY\tDSEC"):
            header_index = index
            header_columns = [column.strip().upper() for column in line.split("\t")]
            break

    if header_index is None:
        raise ValueError(f"No valid data header found in {path}")

    measurements: List[Measurement] = []
    for line in lines[header_index + 1 :]:
        values = [value.strip() for value in line.split("\t")]
        if len(values) < len(header_columns):
            continue

        row = dict(zip(header_columns, values))
        try:
            dsec = float(row["DSEC"])
            temperature_c = float(row["TEMP"])
            total_volume = float(row["VOL"])
            concentration = float(row["CONC"])
            ph = float(row["PH"])
            timestamp = datetime.strptime(row["TIME"], "%Y-%m-%d %H:%M:%S")
        except KeyError as exc:
            raise ValueError(f"Missing field {exc} in {path}") from exc
        except ValueError as exc:
            raise ValueError(f"Invalid row value in {path}: {line}") from exc

        measurements.append(
            Measurement(
                dsec=dsec,
                timestamp=timestamp,
                temperature_c=temperature_c,
                total_volume=total_volume,
                concentration=concentration,
                ph=ph,
            )
        )

    if not measurements:
        raise ValueError(f"No measurement rows found in {path}")
    return measurements


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


def annotate_net_charge(
    rows: List[dict],
    measurements: List[Measurement],
    initial_volume_ul: float,
    hcl_molarity: float,
    extra_charge: float,
    sample_moles: float,
) -> None:
    """Augment each row with the solution net charge in moles."""

    if not measurements:
        return

    initial_ph = measurements[0].ph
    c_naoh = 10 ** (initial_ph - 14.0)
    v_naoh_l = microliters_to_liters(initial_volume_ul)
    n_znaoh = v_naoh_l * c_naoh

    cumulative_vhcl_ul = 0.0
    for row in rows:
        cumulative_vhcl_ul += row["volume_delta"]
        vhcl_l = microliters_to_liters(cumulative_vhcl_ul)
        n_zhcl = vhcl_l * hcl_molarity

        vt_l = microliters_to_liters(row["total_volume"])
        n_zh = (10 ** (-row["ph"])) * vt_l
        n_zoh = (10 ** (row["ph"] - 14.0)) * vt_l

        net_charge = n_zhcl - n_znaoh + n_zh - n_zoh + (extra_charge * sample_moles)
        row["net_charge"] = round(net_charge, 12)
        
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
        "net_charge",
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
        ("net_charge", "Net Charge (mol)"),
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
    parent_name = data_dir.name
    date_fragment = first_timestamp.strftime("%Y%m%d")
    base_name = f"{parent_name}_{date_fragment}"

    target_dir = Path(output_dir) if output_dir else data_dir
    csv_path = target_dir / f"{base_name}.csv"
    plots_dir = target_dir / "plots"
    return csv_path, plots_dir, parent_name


def main() -> None:
    parser = argparse.ArgumentParser(description="Compile TD.mxt measurements into CSV + plots")
    parser.add_argument(
        "data_dir",
        nargs="?",
        default="Blank_4",
        help="Directory containing the folder-specific TD.mxt file (default: Blank_4)",
    )
    parser.add_argument(
        "--output-dir",
        dest="output_dir",
        help="Directory where the CSV and plots should be written (default: parent of data_dir)",
    )
    parser.add_argument(
        "--initial-volume",
        type=float,
        default=None,
        help="Initial solution volume in microliters; defaults to the first measurement total_volume",
    )
    parser.add_argument(
        "--hcl-molarity",
        type=float,
        default=0.5,
        help="Molarity of the titrant HCl in mol/L (default: 0.5)",
    )
    parser.add_argument(
        "--extra-charge",
        type=float,
        default=0.0,
        help="Additional charge per mole carried by the sample (z_ext, default: 0 for blank titration)",
    )
    parser.add_argument(
        "--sample-moles",
        type=float,
        default=0.0,
        help="Number of moles of the sample (n_sample, default: 0 for blank titration)",
    )
    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    if not data_dir.exists():
        raise SystemExit(f"Data directory {data_dir} does not exist")

    if data_dir.is_dir():
        td_filename = f"{data_dir.name}_TD.mxt"
        td_path = data_dir / td_filename
    else:
        td_path = data_dir
        data_dir = data_dir.parent

    if not td_path.exists():
        raise SystemExit(f"TD file {td_path} does not exist")

    measurements = parse_td_file(td_path)
    first_file_timestamp = measurements[0].timestamp
    measurements.sort(key=lambda m: m.dsec)

    rows = build_rows(measurements)
    initial_volume_ul = args.initial_volume if args.initial_volume is not None else measurements[0].total_volume
    annotate_net_charge(
        rows,
        measurements,
        initial_volume_ul=initial_volume_ul,
        hcl_molarity=args.hcl_molarity,
        extra_charge=args.extra_charge,
        sample_moles=args.sample_moles,
    )
    csv_path, plots_dir, parent_name = determine_output_paths(
        data_dir, first_file_timestamp, args.output_dir
    )
    write_csv(rows, csv_path)
    plot_relationships(rows, parent_name, plots_dir)

    print(f"Wrote {csv_path}")
    print(f"Plots saved to {plots_dir}")


if __name__ == "__main__":
    main()
