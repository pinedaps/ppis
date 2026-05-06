#!/usr/bin/env python3

"""
Compute S(q) and S_eff(q) from a many-body I(q) simulation and a P(q)/beta(q) file.

Usage:
    python3 compute_sq.py -N 50 -i I_q/intensity.csv -f p_q/ff_4lzt.csv [-o sq_output.png]
"""

import argparse
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter,
                                prog='python3 compute_sq.py')
    p.add_argument("-N", type=int, required=True,
                   help="Number of proteins in the simulation box")
    p.add_argument("-i", "--intensity", default="I_q/intensity.csv",
                   help="I(q) CSV file with columns q,total  (default: I_q/intensity.csv)")
    p.add_argument("-f", "--formfactor", default="p_q/ff_4lzt.csv",
                   help="P(q) CSV file with columns q,total,beta  (default: p_q/ff_4lzt.csv)")
    p.add_argument("-o", "--output", default=None,
                   help="Output plot filename (default: show interactively)")
    return p.parse_args()


def main():
    args = parse_args()

    iq_df = pd.read_csv(args.intensity)
    pq_df = pd.read_csv(args.formfactor)

    q_iq = iq_df["q"].values
    I_q  = iq_df["total"].values

    q_pq   = pq_df["q"].values
    P_q    = pq_df["total"].values
    beta_q = pq_df["beta"].values

    # Cubic-spline interpolation of P(q) and beta(q) onto the I(q) q-grid.
    # P(q) from pripps starts at q=0; skip q=0 for the spline to avoid
    # extrapolation issues, then interpolate only within the sampled range.
    mask = q_pq > 0
    cs_P    = CubicSpline(q_pq[mask], P_q[mask],    extrapolate=False)
    cs_beta = CubicSpline(q_pq[mask], beta_q[mask], extrapolate=False)

    P_interp    = cs_P(q_iq)
    beta_interp = cs_beta(q_iq)

    # Warn if any q value falls outside the P(q) range
    out_of_range = np.isnan(P_interp) | np.isnan(beta_interp)
    if out_of_range.any():
        print(f"Warning: {out_of_range.sum()} q-points are outside the P(q) range "
              f"[{q_pq[mask].min():.4f}, {q_pq[mask].max():.4f}] A-1 and will be skipped.")
        valid = ~out_of_range
        q_iq        = q_iq[valid]
        I_q         = I_q[valid]
        P_interp    = P_interp[valid]
        beta_interp = beta_interp[valid]

    S_q     = I_q / (args.N * P_interp)
    S_eff_q = 1.0 + beta_interp * (S_q - 1.0)

    # Save results
    out_df = pd.DataFrame({
        "q":     q_iq,
        "I_q":   I_q,
        "P_q":   P_interp,
        "beta":  beta_interp,
        "S_q":   S_q,
        "S_eff": S_eff_q,
    })
    csv_out = args.output.replace(".png", ".csv") if args.output else "sq_results.csv"
    out_df.to_csv(csv_out, index=False, float_format="%.6f")
    print(f"Results saved to {csv_out}")

    def _save_or_show(fig, path, label):
        if path:
            fig.savefig(path, dpi=200)
            print(f"Plot saved to {path}")
        else:
            fig.suptitle(label)
            plt.show()
        plt.close(fig)

    def _stem(output, suffix):
        """Replace .ext with _suffix.ext, or return None for interactive mode."""
        if output is None:
            return None
        base, ext = output.rsplit(".", 1) if "." in output else (output, "png")
        return f"{base}_{suffix}.{ext}"

    # I(q) — log-log
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.loglog(q_iq, I_q, "o-", color="darkorange", linewidth=1.5, markersize=5)
    ax.set_xlabel(r"$q$ / Å$^{-1}$", fontsize=13)
    ax.set_ylabel(r"$I(q)$ / a.u.", fontsize=13)
    ax.set_title(rf"Scattering intensity  ($N={args.N}$)", fontsize=13)
    ax.tick_params(which="both", direction="in")
    plt.tight_layout()
    _save_or_show(fig, _stem(args.output, "Iq"), "I(q)")

    # P(q) — log-log
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.loglog(q_iq, P_interp, "o-", color="seagreen", linewidth=1.5, markersize=5)
    ax.set_xlabel(r"$q$ / Å$^{-1}$", fontsize=13)
    ax.set_ylabel(r"$P(q)$ / a.u.", fontsize=13)
    ax.set_title(r"Form factor $P(q)$", fontsize=13)
    ax.tick_params(which="both", direction="in")
    plt.tight_layout()
    _save_or_show(fig, _stem(args.output, "Pq"), "P(q)")

    # beta(q) — semi-log x
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.semilogx(q_iq, beta_interp, "o-", color="mediumpurple", linewidth=1.5, markersize=5)
    ax.set_xlabel(r"$q$ / Å$^{-1}$", fontsize=13)
    ax.set_ylabel(r"$\beta(q)$", fontsize=13)
    ax.set_title(r"Contrast factor $\beta(q)$", fontsize=13)
    ax.tick_params(which="both", direction="in")
    plt.tight_layout()
    _save_or_show(fig, _stem(args.output, "beta"), r"beta(q)")

    # S(q) — semi-log x
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.semilogx(q_iq, S_q, "o-", color="crimson", linewidth=1.5,
                markersize=5, label=r"$S(q)$")
    ax.axhline(1.0, color="gray", linestyle="--", linewidth=0.8, label="S=1")
    ax.set_xlabel(r"$q$ / Å$^{-1}$", fontsize=13)
    ax.set_ylabel(r"$S(q)$", fontsize=13)
    ax.set_title(rf"Structure factor  ($N={args.N}$)", fontsize=13)
    ax.legend(fontsize=11)
    ax.tick_params(which="both", direction="in")
    plt.tight_layout()
    _save_or_show(fig, _stem(args.output, "Sq"), "S(q)")

    # S_eff(q) — semi-log x
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.semilogx(q_iq, S_eff_q, "o-", color="steelblue", linewidth=1.5,
                markersize=5, label=r"$S_\mathrm{eff}(q)$")
    ax.axhline(1.0, color="gray", linestyle="--", linewidth=0.8, label="S=1")
    ax.set_xlabel(r"$q$ / Å$^{-1}$", fontsize=13)
    ax.set_ylabel(r"$S_\mathrm{eff}(q)$", fontsize=13)
    ax.set_title(rf"Effective structure factor  ($N={args.N}$)", fontsize=13)
    ax.legend(fontsize=11)
    ax.tick_params(which="both", direction="in")
    plt.tight_layout()
    _save_or_show(fig, _stem(args.output, "Seff"), "S_eff(q)")


if __name__ == "__main__":
    main()
