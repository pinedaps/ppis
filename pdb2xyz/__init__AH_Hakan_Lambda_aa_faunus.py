#!/usr/bin/env python3

# Copyright 2025 Mikael Lund
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse
import jinja2
import mdtraj as md
import logging

def parse_args():
    """Parse command line arguments for the script."""
    parser = argparse.ArgumentParser(description="Convert PDB files to XYZ format")
    parser.add_argument(
        "-i", "--infile", type=str, required=True, help="Input PDB file path"
    )
    parser.add_argument(
        "-o", "--outfile", type=str, required=True, help="Output XYZ file path"
    )
    parser.add_argument(
        "-t",
        "--top",
        type=str,
        required=False,
        help="Output topology path (default: topology.yaml)",
        default="topology.yaml",
    )

    parser.add_argument(
        "--pH", type=float, required=False, help="pH value (default: 7.0)", default=7.0
    )
    parser.add_argument(
        "--alpha",
        type=float,
        required=False,
        help="Excess polarizability (default: 0.0)",
        default=0.0,
    )
    parser.add_argument(
        "--T",
        type=float,
        required=False,
        help="Temperature (default: 293K)",
        default=293,
    )
    parser.add_argument(
        "--epsilon",
        type=float,
        required=False,
        help="Epsilon at 293 K (default: 0.8368)",
        default=0.8368,
    )
    parser.add_argument(
        "--sc",
        type=float,
        required=False,
        help="Salt concentration in mol/L (default: 0.115)",
        default=0.115,
    )
    parser.add_argument(
        "--sidechains",
        action="store_true",
        help="Off-center ionizable sidechains (default: disabled)",
        default=False,
    )
    parser.add_argument(
        "--sasa",
        help="Calculate the sasa ratio per aminoacid between the current xyz and the IDP reference (default: disabled)",
        default=False,
    )
    # take list of chain IDs to include (list of strings)
    parser.add_argument(
        "--chains",
        type=str,
        nargs="*",
        required=False,
        help="List of chain IDs to include (default: all chains)",
        default=None,
    )
    
    return parser.parse_args()


def render_template(context: dict):
    template_str = calvados_template()
    return jinja2.Template(template_str).render(context)


def ssbonds(traj):
    """return set of cysteine indices participating in SS-bonds"""
    bonds = traj.topology.bonds
    ss_bonds = []
    for bond in bonds:
        atom1, atom2 = bond
        if (
            atom1.name == "SG"
            and atom1.residue.name == "CYS"
            and atom2.name == "SG"
            and atom2.residue.name == "CYS"
        ):
            ss_bonds.append((atom1.residue.index, atom2.residue.index))
    return set(res for pair in ss_bonds for res in pair)


def convert_pdb(pdb_file: str, output_xyz_file: str, use_sidechains: bool, chains=None):
    """Convert PDB to coarse grained XYZ file; one bead per amino acid"""
    traj = md.load_pdb(pdb_file, frame=0)
    SASA, map = md.shrake_rupley(traj, probe_radius=0.15, n_sphere_points=960, mode='residue', get_mapping=True) 
    cys_with_ssbond = ssbonds(traj)
    residues = []
    for index,res in enumerate(traj.topology.residues):
        if not res.is_protein:
            continue

        if chains is not None and res.chain.chain_id not in chains:
            continue

        cm = [0.0, 0.0, 0.0]  # residue mass center
        mw = 0.0  # residue weight
        for a in res.atoms:
            # Add N-terminal
            if res.index == 0 and a.index == 0 and a.name == "N":
                residues.append(dict(name="NTR", cm=traj.xyz[0][a.index] * 10, sasa=0))
                logging.info("Adding N-terminal bead")

            # Add C-terminal
            if a.name == "OXT":
                residues.append(dict(name="CTR", cm=traj.xyz[0][a.index] * 10, sasa=0))
                logging.info("Adding C-terminal bead")

            # Add coarse grained bead
            cm = cm + a.element.mass * traj.xyz[0][a.index]
            mw = mw + a.element.mass

        # rename CYS -> CSS participating in SS-bonds

        if res.name == "CYS" and res.index in cys_with_ssbond:
            name = "CSS"+str(res.resSeq)
            logging.info(f"Renaming SS-bonded CYS{res.resSeq} to {name}")
        else:
            name = str(res.name)+str(res.resSeq)

        residues.append(dict(name=name, cm=cm / mw * 10, sasa=SASA[0][index] * 100))
        if use_sidechains and not name.startswith("CSS"):
            side_chain = add_sidechain(traj, res)
            if side_chain is not None:
                residues.append(side_chain)

    with open(output_xyz_file, "w") as f:
        f.write(f"{len(residues)}\n")
        f.write(
            f"Converted with Duello pdb2xyz.py with {pdb_file} (https://github.com/mlund/pdb2xyz)\n"
        )
        for i in residues:
            f.write(f"{i['name']} {i['cm'][0]:.3f} {i['cm'][1]:.3f} {i['cm'][2]:.3f}\n")
        logging.info(
            f"Converted {pdb_file} -> {output_xyz_file} with {len(residues)} residues."
        )

    with open(output_xyz_file+'_SASA', "w") as f:
        f.write(
            f"Computed using the Shrake and Rupley algorithm implemented in mdtraj using Duello pdb2xyz.py with {pdb_file} (https://github.com/mlund/pdb2xyz)\n"
        )
        for i in residues:
            f.write(f"{i['name']} {i['sasa']:.3f}\n")
        logging.info(
            f"Converted {pdb_file} -> {output_xyz_file+'_SASA'} with {len(residues)} residues."
        )
    return residues

def add_sidechain(traj, res):
    """Add sidechain bead for ionizable amino acids"""
    # Map residue and atom names to sidechain bead names
    sidechain_map = {
        ("ASP", "OD1"): "Dsc",
        ("GLU", "OE1"): "Esc",
        ("ARG", "CZ"): "Rsc",
        ("LYS", "NZ"): "Ksc",
        ("HIS", "NE2"): "Hsc",
        ("CYS", "SG"): "Csc",
    }
    for atom in res.atoms:
        bead_name = sidechain_map.get((res.name, atom.name))
        if bead_name:
            return dict(name=bead_name, cm=traj.xyz[0][atom.index] * 10)

    if res.name in ["ASP", "GLU", "ARG", "LYS", "HIS", "CYS"]:
        logging.warning(f"Missing sidechain bead for {res.name}{res.index}")
    return None


def write_topology(output_path: str, context: dict):
    """Render and write the topology template."""
    template = calvados_template()
    rendered = jinja2.Template(template).render(context)
    with open(output_path, "w") as file:
        file.write(rendered)
        logging.info(f"Topology written to {output_path}")


def main():
    logging.basicConfig(level=logging.INFO)
    args     = parse_args()
    residues = convert_pdb(args.infile, args.outfile, args.sidechains, args.chains)
    
    context = {
        "residues": residues,
        "pH": args.pH,
        "alpha": args.alpha,
        "sidechains": args.sidechains,
  	    "T": args.T,
	    "ec": args.epsilon,
        "sc": args.sc,
        "xyz_path": args.outfile,
        "sasa": args.sasa,
    }
    write_topology(args.top, context)


# Average pKa values from https://doi.org/10.1093/database/baz024

# Temperature dependece of the hydrophobic interactions via λ(T), using as template Eq.15 described in H. Wennerstrom & B. Lindman (https://doi.org/10.1016/j.molliq.2025.128169), as follows:

# λ(T) = λ_c * (T * ( 1/λ_c + C - C * np.log(T/T_c)) - C * T_c)

# because jinja2 does not support log we can use the Series expansion of the exact solution around Tc, up to the second order

# λ(T) = λ_c * (1 / T_c) * (T - (C / 2) * (T - T_c)**2) 

# The reference values for λ_c parameter per aminoacid at Tc = 294.72 ± 5.58 was taken from G. Tesei & K. Lindorff-Larsen (https://doi.org/10.12688/openreseurope.14967.2). 
# For alkanes, c = ΔG/ΔCp ≈ 1.9e-2 K-1 (Eq. 12 in H. Wennerstrom & B. Lindman, https://doi.org/10.1016/j.molliq.2025.128169).
# For amino acid side chains, c ≈ 2.52e-2 K-1 based on the thermodynamic properties reported in G. I. Makhatadze (https://doi.org/10.1016/S0301-4622(98)00095-7) and V. Pliška, et al. (https://doi.org/10.1016/S0021-9673(00)82337-7).
# Using these parameters, λ(T) becomes negative for 492.2 < T < 176.3 K. Therefore, the formula is applicable only within this temperature range.

def calvados_template():
    return """

{%- set Tc = 294.72 -%}
{%- set c  = 2.52e-2 -%}
{%- set lT = (1 / Tc) * (T - (c / 2) * (T - Tc)**2) -%} 
{%- set f = 1.0 - sidechains -%}

{%- macro calc_charge(res_name, pH) -%}
    {%- if res_name in residues_info -%}
        {%- set info = residues_info[res_name] -%}
        {%- if info.type == "acid" -%}
            {{ - 1 / (1 + 10**(info.pKa - pH)) }}
        {%- elif info.type == "base" -%}
            {{ 1 / (1 + 10**(pH - info.pKa)) }}
        {%- else -%}
            0.0
        {%- endif -%}
    {%- else -%}
        0.0
    {%- endif -%}
{%- endmacro -%}

{%- macro calc_sasa_ratio(res, residues_info) -%}
    {%- set base_name = res.name[:3] -%}
    {%- if base_name in residues_info and residues_info[base_name].sasa_mean > 0 -%}
        {{ res.sasa / residues_info[base_name].sasa_mean }}
    {%- else -%}
        0.0
    {%- endif -%}
{%- endmacro -%}

{%- set residues_info = {
    "CTR": {"mass": 0, "sigma": 0, "hydrophobicity": 0,"pKa": 3.16, "type": "acid", "sasa_mean": 0.0},
    "NTR": {"mass": 0, "sigma": 0, "hydrophobicity": 0,"pKa": 7.64, "type": "base", "sasa_mean": 0.0},
    "ALA": {"mass": 71.09,  "sigma": 5.12, "hydrophobicity": 0.3377244362031627, "pKa": None, "type": None, "sasa_mean": 70.311},
    "ARG": {"mass": 156.19, "sigma": 6.56, "hydrophobicity": 0.7407902764839954, "pKa": 12.5, "type": "base", "sasa_mean": 91.373},
    "ASN": {"mass": 114.1,  "sigma": 5.68, "hydrophobicity": 0.3706962163690402, "pKa": None, "type": None, "sasa_mean": 86.573},
    "ASP": {"mass": 115.09, "sigma": 5.58, "hydrophobicity": 0.0925875575361580, "pKa": 3.43, "type": "acid", "sasa_mean": 77.661},
    "CSS": {"mass": 103.14, "sigma": 5.48, "hydrophobicity": 0.5922529084601322, "pKa": None, "type": None, "sasa_mean": 0},
    "CYS": {"mass": 103.14, "sigma": 5.48, "hydrophobicity": 0.5922529084601322, "pKa": 6.25, "type": "acid", "sasa_mean": 17.558},
    "GLN": {"mass": 128.13, "sigma": 6.02, "hydrophobicity": 0.3143449791669133, "pKa": None, "type": None, "sasa_mean": 118.064},
    "GLU": {"mass": 129.11, "sigma": 5.92, "hydrophobicity": 0.0002495905394260, "pKa": 4.14, "type": "acid", "sasa_mean": 118.928},
    "GLY": {"mass": 57.05,  "sigma": 4.50, "hydrophobicity": 0.7538308115197386, "pKa": None, "type": None, "sasa_mean": 42.462},
    "HIS": {"mass": 137.14, "sigma": 6.08, "hydrophobicity": 0.4087176216525476, "pKa": 6.45, "type": "base", "sasa_mean": 172.739},
    "ILE": {"mass": 113.16, "sigma": 6.18, "hydrophobicity": 0.5130398874425708, "pKa": None, "type": None, "sasa_mean": 79.603},
    "LEU": {"mass": 113.16, "sigma": 6.18, "hydrophobicity": 0.5548615312993875, "pKa": None, "type": None, "sasa_mean": 111.672},
    "LYS": {"mass": 128.17, "sigma": 6.36, "hydrophobicity": 0.1380602542039267, "pKa": 10.68, "type": "base", "sasa_mean": 162.046},
    "MET": {"mass": 131.2,  "sigma": 6.18, "hydrophobicity": 0.5170874160398543, "pKa": None, "type": None, "sasa_mean": 124.479},
    "PHE": {"mass": 147.18, "sigma": 6.36, "hydrophobicity": 0.8906449355499866, "pKa": None, "type": None, "sasa_mean": 116.899},
    "PRO": {"mass": 97.12,  "sigma": 5.56, "hydrophobicity": 0.3469777523519372, "pKa": None, "type": None, "sasa_mean": 83.630},
    "SER": {"mass": 87.08,  "sigma": 5.18, "hydrophobicity": 0.4473142572693176, "pKa": None, "type": None, "sasa_mean": 72.653},
    "THR": {"mass": 101.11, "sigma": 5.62, "hydrophobicity": 0.2672387936544146, "pKa": None, "type": None, "sasa_mean": 77.480},
    "TYR": {"mass": 163.18, "sigma": 6.46, "hydrophobicity": 0.9506286873011070, "pKa": None, "type": None, "sasa_mean": 138.495},
    "VAL": {"mass": 99.13,  "sigma": 5.86, "hydrophobicity": 0.2936174211771383, "pKa": None, "type": None, "sasa_mean": 84.858},
    "TRP": {"mass": 186.22, "sigma": 6.78, "hydrophobicity": 1.0334501235745120, "pKa": None, "type": None, "sasa_mean": 0.0} 
} -%}

comment: "Calvados 3 coarse grained amino acid model for use with Faunus"

pH: {{ pH }}
salt_c: {{ sc }}
T:  {{ T }}
epsilon_c: {{ ec }}
sidechains: {{ sidechains }}
xyz_path: {{ xyz_path }}
SASA: {{ sasa }}

version: 0.2.0

atoms:
{% for res in residues %}
  {%- set base_name = res.name[:3] -%}
  {%- set z = calc_charge(base_name, pH) | float -%}
  {%- set sasa_ratio = 1 if not sasa else (calc_sasa_ratio(res, residues_info) | float) -%}
  - {name: {{ res.name }},
   charge: {{ "%.2f" % z }},
   hydrophobicity: !Lambda {{ residues_info[base_name].hydrophobicity * lT * sasa_ratio }},
   mass: {{ residues_info[base_name].mass }}, 
   σ: {{ residues_info[base_name].sigma }},
   ε: {{ "%.4f" % ec }}{% if residues_info[base_name].type is not none %},
   custom: { alpha: {{ f * alpha }} }{% endif %} }
{% endfor %}

molecules:
- name: MOL1
  degrees_of_freedom: Rigid
  has_com: true
  from_structure: {{ xyz_path }}
- name: MOL2
  degrees_of_freedom: Rigid
  has_com: true
  from_structure: {{ xyz_path }}

system:
  cell: !Sphere {radius: 60.0}
  medium:
    permittivity: !Water
    temperature: {{ T }}
    salt: [!NaCl, {{ sc }}]
  blocks:
  - molecule: MOL1
    N: 1
    insert: !RandomCOM { filename: {{ xyz_path }}, rotate: true, directions: none, offset: [0.0, 0.0, -40.0] }
  - molecule: MOL2
    N: 1
    insert: !RandomCOM { filename: {{ xyz_path }}, rotate: true, directions: none, offset: [0.0, 0.0, 40.0] }
  energy:
    nonbonded:
      default:
        - !Coulomb {cutoff: 1000.0}
        - !AshbaughHatch {mixing: arithmetic, cutoff: 20.0}

analysis:
- !MassCenterDistance
  molecules: [MOL1, MOL2]
  file: com_distance_{{ T }}.dat.gz
  frequency: !Every 100
- !Trajectory
  file: traj_{{ T }}.xyz
  frequency: !Every 100
- !VirtualTranslate
  molecule: MOL2
  dL: 0.1
  directions: !z
  file: "vt.dat_{{ T }}.gz"
  temperature: {{ T }}
  frequency: !Every 100


propagate:
  seed: Hardware
  criterion: Metropolis
  repeat: 10000000
  collections:
  - !Stochastic
    moves:
    - !RotateMolecule
      molecule: MOL1
      dp: 0.28
      weight: 1.0
    - !RotateMolecule
      molecule: MOL2
      dp: 0.28
      weight: 1.0
    - !TranslateMolecule
      molecule: MOL2
      dp: 0.1
      weight: 1.0
      directions: !z

"""

if __name__ == "__main__":
    main()
