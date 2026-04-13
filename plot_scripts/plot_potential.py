import argparse
import json
import os
import re

import matplotlib.pyplot as plt
import numpy as np

import common as cmn

#################################### Loading plotting features ##################################

cmn

EXPORT_DPI = 300


def parse_args():
    parser = argparse.ArgumentParser(description="Configuration of data source")
    parser.add_argument(
        "-p",
        "--paths",
        type=str,
        nargs="+",
        required=True,
        help="one or more paths where the computed scan data are stored",
    )
    parser.add_argument(
        "-pe",
        "--path_exp",
        type=str,
        help="path to the experimental data file",
    )
    parser.add_argument(
        "-s",
        "--sigma",
        type=float,
        help="radius of the protein in Angstroms (default = 10 A)",
    )
    return parser.parse_args()


def infer_dataset_label(path):
    abs_path = os.path.abspath(path)
    stripped = abs_path.rstrip(os.sep)
    base = os.path.basename(stripped)
    if base == "scans":
        parent = os.path.dirname(stripped)
        return os.path.basename(parent)
    return base


def infer_plot_dir(path):
    abs_path = os.path.abspath(path)
    stripped = abs_path.rstrip(os.sep)
    if stripped.endswith("scans"):
        return stripped[:- len("scans")] + "plots"
    return os.path.join(stripped, "plots")


def extract_temperature(name):
    match = re.search(r"scan_T([0-9]+(?:\.[0-9]+)?)", name)
    if match:
        try:
            return float(match.group(1))
        except ValueError:
            return None
    return None


def load_scan_columns(file_path):
    try:
        return np.loadtxt(file_path, comments="#", usecols=(0, 1), unpack=True)
    except ValueError:
        data = np.loadtxt(file_path, comments="#")
        if data.ndim == 1:
            data = data.reshape(1, -1)
        if data.shape[1] < 2:
            raise ValueError(
                f"Scan file {file_path} must contain at least two numeric columns"
            )
        r_values = data[:, 0]
        u_column = 2 if data.shape[1] > 2 else data.shape[1] - 1
        u_values = data[:, u_column]
        return r_values, u_values


def determine_label(dataset_label, temp, scans_in_dataset, datasets_count):
    if datasets_count == 1:
        if scans_in_dataset > 1 and temp is not None:
            return f"{temp:g} K"
        if temp is not None:
            return f"{dataset_label} (T={temp:g} K)"
        return dataset_label
    if scans_in_dataset == 1 or temp is None:
        return dataset_label
    return f"{dataset_label} (T={temp:g} K)"


def load_experimental_data(exp_path):
    if not exp_path:
        return None
    exp_file = os.path.abspath(exp_path)
    if not os.path.isfile(exp_file):
        print(f"Experimental data file {exp_file} not found; skipping experimental curve.")
        return None
    
    # Read first line to get label
    label = None
    with open(exp_file, 'r') as f:
        first_line = f.readline()
        if first_line.startswith('#'):
            label = first_line.lstrip('#').strip()
    
    if label is None:
        label = "Experimental data"
    
    print(f"Loading experimental data from {exp_file}")
    data = np.loadtxt(exp_file, unpack=True)
    
    return {
        'label': label,
        'data': data
    }


def plot_potentials(potential_series, plot_dir):
    fig, ax = plt.subplots()
    ax.axhline(y=0, color="black", lw=1, ls="--", alpha=0.5)
    for entry in potential_series:
        ax.plot(
            entry["R"],
            entry["U"],
            ms=2.5,
            marker="o",
            lw=1,
            alpha=0.5,
            label=entry["label"],
        )
    ax.set_ylim(-2, 2)
    ax.set_xlim(23, 60)
    ax.set_xlabel(r"Distance, r [$\AA$]")
    ax.set_ylabel(r"Free Energy, F(r) [$k_BT$]")
    ax.legend(ncol=1)
    ax.set_title("Interaction Free Energy of $\gamma$B-crystallin", pad=15)
    os.makedirs(plot_dir, exist_ok=True)
    plt.savefig(os.path.join(plot_dir, "potential.png"), dpi=EXPORT_DPI)
    plt.close(fig)


def plot_b2(b2_series, exp_data, plot_dir):
    fig, ax = plt.subplots()
    if exp_data is not None:
        label = exp_data['label']
        T_exp, b2_red_exp, b2_min, b2_max = exp_data['data']
        yerr = np.vstack([b2_min, b2_max])
        ax.errorbar(
            T_exp,
            b2_red_exp,
            yerr=yerr,
            fmt="o",
            ms=5,
            lw=1,
            capsize=3,
            alpha=0.5,
            label=label,
        )
    for series in b2_series:
        ax.plot(
            series["T"],
            series["values"],
            ms=5,
            marker="o",
            lw=1,
            alpha=0.5,
            label=series["label"],
        )
    ax.set_title("Reduced second virial coefficient ${b_2}^*$ \n of $\gamma$B-crystallin", pad=15)
    ax.legend()
    ax.set_xlabel(r"Temperature, T [$K$]")
    ax.set_ylabel(r"$b_2^*=b_2/b_2^{HS}$")
    os.makedirs(plot_dir, exist_ok=True)
    plt.savefig(os.path.join(plot_dir, "B2.png"), dpi=EXPORT_DPI)
    plt.close(fig)


def main():
    args = parse_args()
    data_dirs = [os.path.abspath(path) for path in args.paths]
    for directory in data_dirs:
        if not os.path.isdir(directory):
            raise FileNotFoundError(f"Data directory {directory} not found")

    if args.sigma:
        sigma = args.sigma
    else:
        sigma = 10
        print("Protein radius not provided, a default radius of 1 nm was used!")
    B2_HS = 2 / 3 * np.pi * sigma**3
    print("B2_HS = " + str(B2_HS))

    plot_dir = infer_plot_dir(data_dirs[0])
    potential_series = []
    b2_series = []
    datasets_count = len(data_dirs)

    for directory in data_dirs:
        dataset_label = infer_dataset_label(directory)
        entries = sorted([name for name in os.listdir(directory) if name.endswith(".dat")])
        if not entries:
            print(f"Warning: no .dat scan files found in {directory}")
            continue
        scans_in_dataset = len(entries)
        for entry in entries:
            file_path = os.path.join(directory, entry)
            print("Processing file:", file_path)
            R, U = load_scan_columns(file_path)
            temp = extract_temperature(entry)
            label = determine_label(dataset_label, temp, scans_in_dataset, datasets_count)
            potential_series.append({"label": label, "R": R, "U": U})

        json_files = sorted([name for name in os.listdir(directory) if name.endswith(".json")])
        dataset_b2 = {"label": dataset_label, "T": [], "values": []}
        for json_file in json_files:
            json_path = os.path.join(directory, json_file)
            temp = extract_temperature(json_file)
            try:
                with open(json_path, "r") as handle:
                    payload = json.load(handle)
            except (OSError, json.JSONDecodeError) as err:
                print(f"Warning: unable to read {json_path}: {err}")
                continue
            if "B2" not in payload:
                print(f"Warning: missing B2 entry in {json_path}")
                continue
            dataset_b2["T"].append(temp if temp is not None else len(dataset_b2["T"]))
            dataset_b2["values"].append(payload["B2"] / B2_HS)
        if dataset_b2["T"]:
            b2_series.append(dataset_b2)
        elif json_files:
            print(f"Warning: no usable B2 data found in {directory}")
        else:
            print(f"No B2 JSON files found in {directory}; skipping this dataset in the B2 plot.")

    if not potential_series:
        raise RuntimeError("No potential data found. Please check the provided paths.")

    plot_potentials(potential_series, plot_dir)

    exp_data = load_experimental_data(args.path_exp)
    if b2_series:
        plot_b2(b2_series, exp_data, plot_dir)
    else:
        print("No B2 data found across the provided datasets. Skipping B2 plot.")


if __name__ == "__main__":
    main()
