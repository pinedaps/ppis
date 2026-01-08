import numpy as np
from collections import defaultdict
import argparse
import os 

parser = argparse.ArgumentParser(description="Average SASA")
parser.add_argument("-i", nargs="+", type=str, help="input path")
parser.add_argument("-o", type=str, default=None, help="output filename (default: auto-generated from inputs)")
args = parser.parse_args()

files           = args.i

if args.o is None:
    base_names = [os.path.splitext(os.path.basename(f))[0] for f in files]
    output_filename = "sasa_average_" + "_".join(base_names).replace('_SASA','')
else:
    output_filename = args.o

data = defaultdict(list)
for fname in files:
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # skip header / comment lines
            if line.startswith("Computed") or line.startswith("#"):
                continue

            res, value = line.split()
            data[res].append(float(value))

# write output
with open(output_filename, "w") as out:
    out.write("Residue  Mean  StdDev  StdError\n")
    for res in sorted(data):
        values = np.asarray(data[res], dtype=float)
        n = len(values)

        if n == 0:
            mean = np.nan
            stderr = np.nan
            stderr = np.nan

        elif n == 1:
            mean = values[0]
            stddev = 0.0   # no spread with one sample
            stderr = 0.0   # no uncertainty with one sample

        else:
            mean = values.mean()
            stddev = values.std(ddof=1)          # sample SD
            stderr = stddev / np.sqrt(n)         # SE = SD / sqrt(N)

        out.write(f"{res:4s} {mean:8.3f} {stddev:8.3f} {stderr:8.3f}\n")


