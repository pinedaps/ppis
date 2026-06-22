import argparse
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.expanduser('~/projects/ppis/plot_scripts'))
import common_presentation  # noqa: F401 — applies rcParams side-effects


def percus_yevick_Sq(q, sigma, phi):
    """
    Analytical Percus-Yevick hard sphere structure factor.

    Parameters
    ----------
    q     : array, wavevectors in same units as sigma (e.g. Å⁻¹)
    sigma : float, hard sphere diameter (Å)
    phi   : float, volume fraction (0 to ~0.64)

    Returns
    -------
    S(q) : array
    """
    x = q * sigma

    a = (1 + 2*phi)**2 / (1 - phi)**4
    b = -6*phi*(1 + phi/2)**2 / (1 - phi)**4
    g = phi * a / 2

    x = np.where(x < 1e-10, 1e-10, x)

    A = a * (np.sin(x) - x*np.cos(x)) / x**3
    B = b * (2*x*np.sin(x) - (x**2 - 2)*np.cos(x) - 2) / x**4
    G = g * (-x**4*np.cos(x) + 4*((3*x**2 - 6)*np.cos(x)
             + (x**3 - 6*x)*np.sin(x) + 6)) / x**6

    # ĉ(q) = -4πσ³(A+B+G) with units Å³; ρ = 6φ/(πσ³) in Å⁻³
    # ρ·ĉ(q) = (6φ/πσ³)·(-4πσ³)·(A+B+G) = -24φ(A+B+G)  [dimensionless]
    return 1.0 / (1.0 + 24.0 * phi * (A + B + G))


def hard_sphere_Pq(q, sigma):
    """Analytical sphere form factor P(q), normalised to 1 at q=0."""
    qr = q * sigma / 2
    qr = np.where(qr < 1e-10, 1e-10, qr)
    F = 3 * (np.sin(qr) - qr*np.cos(qr)) / qr**3
    return F**2


# ── CLI ───────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(
    description="Hard-sphere Percus-Yevick benchmark vs. pripps simulation.")
parser.add_argument("--sigma", type=float, default=40.0,
                    help="Hard sphere diameter in Å (default: 40.0)")
parser.add_argument("--phi", type=float, default=0.30,
                    help="Volume fraction (default: 0.30)")
parser.add_argument("--N", type=int, default=100,
                    help="Number of spheres in the simulation (default: 100)")
parser.add_argument("--csv", type=str, default=None,
                    help="Path to pripps intensity.csv output (optional)")
args = parser.parse_args()

sigma = args.sigma
phi   = args.phi
N     = args.N

# ── Simulation (pripps, f=1 constant ff) ──────────────────────────────────────
has_sim = False
if args.csv is not None:
    try:
        data  = np.loadtxt(args.csv, delimiter=",", skiprows=1)
        q_sim = data[:, 0]
        S_sim = data[:, 1] / N
        P_sim = hard_sphere_Pq(q_sim, sigma)
        I_full_sim = P_sim * S_sim
        has_sim = True
        print(f"Loaded {len(q_sim)} q-points from '{args.csv}'")
    except FileNotFoundError:
        print(f"Warning: file '{args.csv}' not found — plotting analytical only.")

# ── Analytical PY (clipped to simulation q range if available) ─────────────────
q_max = q_sim.max() if has_sim else 0.8
q    = np.linspace(0.01, q_max, 500)
S_PY = percus_yevick_Sq(q, sigma, phi)
P_HS = hard_sphere_Pq(q, sigma)
I_PY = P_HS * S_PY

stem = f'HS_benchmark_sigma{sigma:.0f}_phi{phi:.2f}'

# ── S(q) ──────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(q, S_PY, color='#EF9F27', lw=2.5, label='Percus-Yevick analytical')
if has_sim:
    ax.plot(q_sim, S_sim, color='#E05252', lw=1.5, ls='--', label='Simulation')
ax.axhline(1, color='gray', lw=0.8, ls=':')
ax.set_xlabel(r"$q$ [Å$^{-1}$]")
ax.set_ylabel(r"Structure factor, $S(q)$")
ax.legend()
ax.tick_params(which='both', direction='in')
ax.yaxis.set_major_locator(plt.MultipleLocator(0.5))
ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.1f'))
fig.savefig(f'{stem}_Sq.png', dpi=150, transparent=True, bbox_inches='tight')
plt.show()
plt.close(fig)

# ── P(q) ──────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7, 5))
ax.semilogy(q, P_HS, color='#1D9E75', lw=2.5)
ax.set_xlabel(r"$q$ [Å$^{-1}$]")
ax.set_ylabel(r"Form factor, $P(q)$")
ax.tick_params(which='both', direction='in')
fig.savefig(f'{stem}_Pq.png', dpi=150, bbox_inches='tight')
plt.show()
plt.close(fig)

# ── I(q) ──────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7, 5))
ax.semilogy(q, I_PY, color='#378ADD', lw=2.5, label='Percus-Yevick analytical')
if has_sim:
    ax.semilogy(q_sim, I_full_sim, color='#E05252', lw=1.5, ls='--',
                label=r'$P(q)\cdot S_\mathrm{sim}(q)$')
ax.set_xlabel(r"$q$ [Å$^{-1}$]")
ax.set_ylabel(r"Intensity, $I(q)$")
ax.legend()
ax.tick_params(which='both', direction='in')
fig.savefig(f'{stem}_Iq.png', dpi=150, bbox_inches='tight')
plt.show()
plt.close(fig)

# ── Key features ──────────────────────────────────────────────────────────────
q_peak = q[np.argmax(S_PY)]
print(f"PY  S(q) peak:  q = {q_peak:.4f} Å⁻¹  →  d = {2*np.pi/q_peak:.1f} Å  (expected ≈ σ = {sigma:.0f} Å)")
print(f"PY  S(q=0) = {S_PY[0]:.4f}  (compressibility limit)")
print(f"Volume fraction φ = {phi:.3f},  N = {N}")
if has_sim:
    q_peak_sim = q_sim[np.argmax(S_sim)]
    print(f"Sim S(q) peak:  q = {q_peak_sim:.4f} Å⁻¹  →  d = {2*np.pi/q_peak_sim:.1f} Å")
    print(f"Sim S(q→0) ≈ {S_sim[0]:.4f}")
