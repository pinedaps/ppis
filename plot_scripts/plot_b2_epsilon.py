import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata
from math import ceil
import argparse
import os
import json
import common_presentation as cmn
#import common as cmn

#################################### Loading plotting features ##################################

cmn

######################################### Acquiring path ########################################

parser = argparse.ArgumentParser(description='Configuration of data source')
parser.add_argument('-p', '--path', type=str, 
                    help='path where the computed data is stored')
parser.add_argument('-pe', '--path_exp', type=str,
                    help='path where the experimental data is stored')
parser.add_argument('-s', '--sigma', type=float,
                    help='radius of the protein in Angtroms (default= 10 A)')
args = parser.parse_args()

data_dir = args.path
plot_dir = args.path.replace('scans','plots')
exp_dir  = args.path_exp

# Calculation of the B2 for a hard sphere of radius 'sigma'

if args.sigma:
	sigma = args.sigma
else:
	sigma = 10
	print('Protein radius not provided, a default radius of 1 nm was used!')
B2_HS = 2/3*np.pi*sigma**3
print('\nB2_HS = ' +str(B2_HS)+'\n')

# Extract T values from filename and create the list of files 

epsilon   = []
b2        = []
diff      = []
for name in sorted(os.listdir(data_dir),reverse=True):
    if name.endswith(".json"):
        epsilon.append(float(name.replace('scan_epsilon','').replace('.json','')))
        with open(data_dir+name, 'r') as file:
            b2_val = json.load(file)['B2']/B2_HS
            b2.append(b2_val)
            diff.append(b2_val-B2_HS)
# Extract B2_reduced from experimental data

T_exp, b2_red_exp, b2_min, b2_max = np.loadtxt(exp_dir+os.listdir(exp_dir)[0], unpack=True)

print('Experimental reduced B2 at 293K :')
print(b2_red_exp[0],'\n')
print('Evaluated epsilon :')
print(epsilon,'\n')
print('Computed reduced B2:')
print(b2,'\n')
print('Difference from experiments')
print(b2-b2_red_exp[0])

fig, ax = plt.subplots(nrows=1)

ax.plot(b2, epsilon, ms=5, marker='o', lw=1)
plt.axhline(y=b2_red_exp[0], color='black', lw=1, ls='--', alpha=0.5)
plt.xlabel(r"Epsilon, [$\epsilon$]")
plt.ylabel(r"$b_2^*=b_2/b_2^{HS}$")
plt.savefig(plot_dir+"B2.png")
