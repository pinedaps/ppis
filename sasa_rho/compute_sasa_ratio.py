#!/usr/bin/env python3

import sys
import math

def read_sasa_file(filename):
    """
    Reads a SASA file into a dict:
    { 'ALA': (mean, stdev, stderr), ... }
    """
    data = {}
    with open(filename) as f:
        header = f.readline()
        for line in f:
            if not line.strip():
                continue
            aa, mean, dev, err = line.split()
            data[aa] = (float(mean), float(dev), float(err))
    return data


def main(file1, file2, outfile):
    sasa1 = read_sasa_file(file1)
    sasa2 = read_sasa_file(file2)

    # Union of all amino acids appearing in either file
    amino_acids = sorted(set(sasa1.keys()) | set(sasa2.keys()))

    with open(outfile, "w") as out:
        out.write("AA\tRatio_SASA\tProp_SD\tProp_SE\n")

        for aa in amino_acids:
            # SASA_1: missing → assume 0
            mean1, dev1, err1 = sasa1.get(aa, (0.0, 0.0, 0.0))
            # SASA_2: missing → treated as zero
            mean2, dev2, err2 = sasa2.get(aa, (0.0, 0.0, 0.0))

            if mean2 == 0.0:
                ratio = 0.0
                prop_sd = 0.0
                prop_se = 0.0
            else:
                ratio = mean1 / mean2

                # Propagate StdDev
                term1_sd = (dev1 / mean1) ** 2 if mean1 != 0 else 0.0
                term2_sd = (dev2 / mean2) ** 2

                prop_sd = ratio * math.sqrt(term1_sd + term2_sd)

                # Propagate StdErr
                term1_se = (err1 / mean1) ** 2 if mean1 != 0 else 0.0
                term2_se = (err2 / mean2) ** 2

                prop_se = ratio * math.sqrt(term1_se + term2_se)

            out.write(f"{aa}\t{ratio:.6g}\t{prop_sd:.6g}\t{prop_se:.6g}\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(
            "Usage: python3 compute_sasa_ratio.py SASA_1 SASA_2 SASA_ratio"
        )

    main(sys.argv[1], sys.argv[2], sys.argv[3])
