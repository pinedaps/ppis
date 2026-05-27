import argparse
import numpy as np
import matplotlib.pyplot as plt

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

# ── Analytical PY ─────────────────────────────────────────────────────────────
q    = np.linspace(0.01, 0.8, 500)
S_PY = percus_yevick_Sq(q, sigma, phi)
P_HS = hard_sphere_Pq(q, sigma)
I_PY = P_HS * S_PY

# ── Simulation (pripps, f=1 constant ff) ──────────────────────────────────────
# pripps Debye: I(q) = Σᵢ fᵢ² + 2Σᵢ<ⱼ fᵢfⱼ sinc(q·rᵢⱼ)
# with fᵢ=1  →  I_pripps(q) = N·S(q)  ⟹  S_sim = I_pripps / N
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

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

axes[0].plot(q, S_PY, color='#EF9F27', lw=2.5, label='PY analytical')
if has_sim:
    axes[0].plot(q_sim, S_sim, color='#E05252', lw=1.5, ls='--', label='Simulation')
axes[0].axhline(1, color='gray', lw=0.8, ls=':')
axes[0].set_xlabel('q (Å⁻¹)')
axes[0].set_ylabel('S(q)')
axes[0].set_title(f'Structure factor  φ={phi}  σ={sigma}Å')
axes[0].legend()
axes[0].grid(True)

axes[1].semilogy(q, P_HS, color='#1D9E75', lw=2.5)
axes[1].set_xlabel('q (Å⁻¹)')
axes[1].set_ylabel('P(q)')
axes[1].set_title('Hard sphere form factor P(q)')
axes[1].grid(True)

axes[2].semilogy(q, I_PY, color='#378ADD', lw=2.5, label='PY analytical')
if has_sim:
    axes[2].semilogy(q_sim, I_full_sim, color='#E05252', lw=1.5, ls='--',
                     label='P(q)·S_sim(q)')
axes[2].set_xlabel('q (Å⁻¹)')
axes[2].set_ylabel('I(q)')
axes[2].set_title('I(q) = P(q)·S(q)')
axes[2].legend()
axes[2].grid(True)

plt.suptitle('Hard sphere Percus-Yevick benchmark', color='#EF9F27', fontsize=13)
plt.tight_layout()
plt.savefig(f'HS_benchmark_sigma{sigma:.0f}_phi{phi:.2f}.png', dpi=150, bbox_inches='tight')
plt.show()

# ── Key features ──────────────────────────────────────────────────────────────
q_peak = q[np.argmax(S_PY)]
print(f"PY  S(q) peak:  q = {q_peak:.4f} Å⁻¹  →  d = {2*np.pi/q_peak:.1f} Å  (expected ≈ σ = {sigma:.0f} Å)")
print(f"PY  S(q=0) = {S_PY[0]:.4f}  (compressibility limit)")
print(f"Volume fraction φ = {phi:.3f},  N = {N}")
if has_sim:
    q_peak_sim = q_sim[np.argmax(S_sim)]
    print(f"Sim S(q) peak:  q = {q_peak_sim:.4f} Å⁻¹  →  d = {2*np.pi/q_peak_sim:.1f} Å")
    print(f"Sim S(q→0) ≈ {S_sim[0]:.4f}")
