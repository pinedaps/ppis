import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os
import glob

parser = argparse.ArgumentParser(description='Plot data from Faunus simulations.')
parser.add_argument('--dat_dir', required=True, help='Root output directory containing per-epsilon subdirectories')
parser.add_argument('--plots_dir', required=True, help='Directory to save the plot images')
args = parser.parse_args()

dat_dir = args.dat_dir
plots_dir = args.plots_dir

os.makedirs(plots_dir, exist_ok=True)

# Discover per-epsilon dat directories: <dat_dir>/*/results/dat/
eps_dat_dirs = sorted(glob.glob(os.path.join(dat_dir, '*/results/dat')))

if not eps_dat_dirs:
    raise FileNotFoundError(f"No per-epsilon subdirectories found under {dat_dir}")

def eps_label(path):
    # Extract the immediate parent of 'results/dat', e.g. 'eps_ana_0.8368'
    return os.path.basename(os.path.dirname(os.path.dirname(path)))

# Plot energy data
fig_e, ax_e = plt.subplots(figsize=(10, 6))
for d in eps_dat_dirs:
    df = pd.read_csv(f'{d}/energy.csv.gz', compression='gzip').apply(pd.to_numeric)
    cut = int(0.2 * len(df))
    ax_e.plot(df['step'][cut:], df['total'][cut:], label=eps_label(d))
ax_e.set_xlabel('Step')
ax_e.set_ylabel('Energy')
ax_e.set_title('Energy Components Over Steps')
ax_e.legend()
fig_e.savefig(f'{plots_dir}/energy_plot.png')
plt.close(fig_e)

# Plot hydrophobic energy data
fig_h, ax_h = plt.subplots(figsize=(10, 6))
for d in eps_dat_dirs:
    df = pd.read_csv(f'{d}/hydrophobic_energy.csv.gz', compression='gzip').apply(pd.to_numeric)
    cut = int(0.2 * len(df))
    ax_h.plot(df['step'][cut:], df['energy'][cut:], label=eps_label(d))
ax_h.set_xlabel('Step')
ax_h.set_ylabel('Energy')
ax_h.set_title('Hydrophobic Energy Over Steps')
ax_h.legend()
fig_h.savefig(f'{plots_dir}/hydrophobic_energy_plot.png')
plt.close(fig_h)

# Plot RDF data
fig_r, ax_r = plt.subplots(figsize=(10, 6))
for d in eps_dat_dirs:
    df = pd.read_csv(f'{d}/rdf_com.dat.gz', compression='gzip',
                     sep=r'\s+', comment='#', header=None, names=['r', 'g(r)']).apply(pd.to_numeric)
    ax_r.plot(df['r'], df['g(r)'], label=eps_label(d))
ax_r.set_xlabel('r')
ax_r.set_ylabel('g(r)')
ax_r.set_ylim(-0.1, 4)
ax_r.set_title('Radial Distribution Function')
ax_r.legend()
fig_r.savefig(f'{plots_dir}/rdf_plot.png')
plt.close(fig_r)
