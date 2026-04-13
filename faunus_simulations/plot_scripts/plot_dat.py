import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os

parser = argparse.ArgumentParser(description='Plot data from Faunus simulations.')
parser.add_argument('--dat_dir', required=True, help='Directory containing the .gz data files')
parser.add_argument('--plots_dir', required=True, help='Directory to save the plot images')
args = parser.parse_args()

dat_dir = args.dat_dir
plots_dir = args.plots_dir

os.makedirs(plots_dir, exist_ok=True)

# Plot energy data
df_energy = pd.read_csv(f'{dat_dir}/energy.csv.gz', compression='gzip', sep='\s+', header=None, names=['step', 'celloverlap', 'nonbonded', 'intramolecular', 'total'], skiprows=1)
df_energy = df_energy.apply(pd.to_numeric)
plt.figure(figsize=(10, 6))
cut_e = int(0.2 * len(df_energy))
plt.plot(df_energy['step'][cut_e:], df_energy['total'][cut_e:], label='Total Energy')
#plt.plot(df_energy['step'], df_energy['nonbonded'], label='Nonbonded Energy')
plt.xlabel('Step')
plt.ylabel('Energy')
plt.title('Energy Components Over Steps')
plt.legend()
plt.savefig(f'{plots_dir}/energy_plot.png')

# Plot hydrophobic energy data
df_hydro = pd.read_csv(f'{dat_dir}/hydrophobic_energy.csv.gz', compression='gzip', sep='\s+', header=None, names=['step', 'energy', 'average'], skiprows=1)
df_hydro = df_hydro.apply(pd.to_numeric)
plt.figure(figsize=(10, 6))
cut_h = int(0.2 * len(df_hydro))
plt.plot(df_hydro['step'][cut_h:], df_hydro['energy'][cut_h:], label='Hydrophobic Energy')
#plt.plot(df_hydro['step'], df_hydro['average'], label='Average Hydrophobic Energy')
plt.xlabel('Step')
plt.ylabel('Energy')
plt.title('Hydrophobic Energy Over Steps')
plt.legend()
plt.savefig(f'{plots_dir}/hydrophobic_energy_plot.png')

# Plot RDF data
df_rdf = pd.read_csv(f'{dat_dir}/rdf_com.dat.gz', compression='gzip', sep='\s+', header=None, names=['r', 'g(r)'], skiprows=1)
df_rdf = df_rdf.apply(pd.to_numeric)
plt.figure(figsize=(10, 6))
plt.plot(df_rdf['r'], df_rdf['g(r)'])
plt.xlabel('r')
plt.ylabel('g(r)')
plt.ylim(-0.1, 4)
plt.title('Radial Distribution Function')
plt.savefig(f'{plots_dir}/rdf_plot.png')
