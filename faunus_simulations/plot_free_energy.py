#!/usr/bin/env python3
"""
Plot free energy comparison between Monte Carlo simulation and partition function scan.

This script generates a figure comparing the potential of mean force (PMF) obtained from
Monte Carlo simulations (from COM distance data) with results from a partition function scan.
"""

import argparse
import matplotlib
import numpy as np
import matplotlib.pyplot as plt


def com_distance_to_pmf(file, bins):
    """Load COM distance as a function of MC steps and calculate potential of mean force (PMF)
    
    Args:
        file (str): Path to the COM distance data file (gzipped or plain text)
        bins (int): Number of bins for the histogram
        
    Returns:
        tuple: (r, pmf) where r is the distance array and pmf is the potential of mean force
    """
    p, r = np.histogram(np.loadtxt(file), bins=bins)
    r = r[1:] / 2 + r[0:-1] / 2
    dr = r[1] - r[0]
    r_max = r[-1]
    gofr = p / p.sum() * r_max / dr  # g(r)
    pmf = -np.log(gofr)              # w(r)/kT
    print(f"r_min = {r[0]}, r_max = {r_max} Å, dr = {dr} Å")
    return r, pmf


def main():
    """Main function to create and save the free energy plot"""
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Plot free energy comparison between MC simulation and partition function scan"
    )
    parser.add_argument(
        "--mc",
        type=str,
        help="Path to Monte Carlo COM distance data file"
    )
    parser.add_argument(
        "--pf",
        type=str,
        help="Path to partition function scan file"
    )
    parser.add_argument(
        "--bins",
        type=int,
        default=100,
        help="Number of bins for histogram (default: 100)"
    )
    parser.add_argument(
        "--output",
        type=str,
        default="free_energy.png",
        help="Output file name for the plot (default: free_energy.png)"
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="DPI for output image (default: 300)"
    )
    args = parser.parse_args()
    
    fig, ax = plt.subplots(nrows=1, figsize=(5, 5))
    
    # Load and plot Monte Carlo results
    r, pmf = com_distance_to_pmf(args.mc, bins=args.bins)
    ax.hlines(0, xmin=r[0]-10, xmax=r[-1]+10, color='k', alpha=0.5, ls='--', lw=2)
    ax.plot(r, pmf - 0.7, 'bo', lw=3, label="Metropolis-Hastings Monte-Carlo", 
            alpha=0.7, ms=4)

    # Load and plot partition function scan results
    r, pmf = np.loadtxt(args.pf, usecols=[0, 1], unpack=True)
    ax.plot(r, pmf, label="From partition function", lw=2, alpha=0.8, color='red')

    # Configure plot
    ax.set_ylabel("Free energy ($k_BT$)")
    ax.set_xlabel("Mass center separation, $R$ (Å)")
    ax.legend(loc=0, frameon=False)
    ax.set_xlim(23, 120)
    ax.set_ylim(-4.2, 14)
    
    # Save and display
    plt.savefig(args.output, bbox_inches="tight", dpi=args.dpi)
    print(f"Plot saved as: {args.output}")
    plt.show()


if __name__ == "__main__":
    main()
