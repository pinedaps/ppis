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
        "--saltcon",
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
            name = "CSS"
            logging.info(f"Renaming SS-bonded CYS{res.index} to {name}")
        else:
            name = res.name

        residues.append(dict(name=name, cm=cm / mw * 10, sasa=SASA[0][index] * 100))
        if use_sidechains and name != "CSS":
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
    args = parse_args()
    convert_pdb(args.infile, args.outfile, args.sidechains, args.chains)

    context = {
        "pH": args.pH,
        "alpha": args.alpha,
        "sidechains": args.sidechains,
	    "T": args.T,
	    "epsilon": args.epsilon,
        "saltcon": args.saltcon,
        "xyz_path": args.outfile,
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

{%- set f = 1.0 - sidechains -%}
{%- set zCTR = - 10**(pH-3.16) / (1 + 10**(pH-3.16)) -%}
{%- set zASP = - 10**(pH-3.43) / (1 + 10**(pH-3.43)) -%}
{%- set zGLU = - 10**(pH-4.14) / (1 + 10**(pH-4.14)) -%}
{%- set zCYS = 10**(pH-6.25) / (1 + 10**(pH-6.25)) -%}
{%- set zHIS = 1 - 10**(pH-6.45) / (1 + 10**(pH-6.45)) -%}
{%- set zNTR = 1 - 10**(pH-7.64) / (1 + 10**(pH-7.64)) -%}
{%- set zLYS = 1 - 10**(pH-10.68) / (1 + 10**(pH-10.68)) -%}
{%- set zARG = 1 - 10**(pH-12.5) / (1 + 10**(pH-12.5)) -%}

{%- set Tc = 294.72 -%}
{%- set c  = 2.52e-2 -%}
{%- set LambdaT = (1 / Tc) * (T - (c / 2) * (T - Tc)**2) -%} 

comment: "Calvados 3 coarse grained amino acid model for use with Faunus"

T:  {{ T }}
pH: {{ pH }}
salt_c: {{ saltcon }}
epsilon: {{ epsilon }}
sidechains: {{ sidechains }}
xyz_path: {{ xyz_path }}

atoms:
  - {charge: {{ "%.2f" % zCTR }}, hydrophobicity: !Lambda 0, mass: 0, name: CTR, σ: 2.0, ε: {{ "%.4f" % epsilon }}}
  - {charge: {{ "%.2f" % zNTR }}, hydrophobicity: !Lambda 0, mass: 0, name: NTR, σ: 2.0, ε: {{ "%.4f" % epsilon }}}
{%- if sidechains %}
  - {charge: {{ "%.2f" % zGLU }}, hydrophobicity: !Lambda 0, mass: 0, name: Esc, σ: 2.0, ε: {{ "%.4f" % epsilon }}}
  - {charge: {{ "%.2f" % zASP }}, hydrophobicity: !Lambda 0, mass: 0, name: Dsc, σ: 2.0, ε: {{ "%.4f" % epsilon }}}
  - {charge: {{ "%.2f" % zHIS }}, hydrophobicity: !Lambda 0, mass: 0, name: Hsc, σ: 2.0, ε: {{ "%.4f" % epsilon }}}
  - {charge: {{ "%.2f" % zARG }}, hydrophobicity: !Lambda 0, mass: 0, name: Rsc, σ: 2.0, ε: {{ "%.4f" % epsilon }}}
  - {charge: {{ "%.2f" % zLYS }}, hydrophobicity: !Lambda 0, mass: 0, name: Ksc, σ: 2.0, ε: {{ "%.4f" % epsilon }}}
  - {charge: {{ "%.2f" % zCYS }}, hydrophobicity: !Lambda 0, mass: 0, name: Csc, σ: 2.0, ε: {{ "%.4f" % epsilon }}}
{%- endif %}
  - {charge: {{ "%.2f" % (zARG * f) }}, hydrophobicity: !Lambda {{0.7407902764839954 * LambdaT}}, mass: 156.19, name: ARG, σ: 6.56, ε: {{ "%.4f" % epsilon }}, custom: {alpha: {{ f * alpha }}}}
  - {charge: {{ "%.2f" % (zASP * f) }}, hydrophobicity: !Lambda {{0.092587557536158* LambdaT}},  mass: 115.09, name: ASP, σ: 5.58, ε: {{ "%.4f" % epsilon }}, custom: {alpha: {{ f * alpha }}}}
  - {charge: {{ "%.2f" % (zGLU * f) }}, hydrophobicity: !Lambda {{0.000249590539426 * LambdaT}},  mass: 129.11, name: GLU, σ: 5.92, ε: {{ "%.4f" % epsilon }}, custom: {alpha: {{ f * alpha }}}}
  - {charge: {{ "%.2f" % (zLYS * f) }}, hydrophobicity: !Lambda {{0.1380602542039267 * LambdaT}}, mass: 128.17, name: LYS, σ: 6.36, ε: {{ "%.4f" % epsilon }}, custom: {alpha: {{ f * alpha }}}}
  - {charge: {{ "%.2f" % (zHIS * f) }}, hydrophobicity: !Lambda {{0.4087176216525476 * LambdaT}}, mass: 137.14, name: HIS, σ: 6.08, ε: {{ "%.4f" % epsilon }}, custom: {alpha: {{ f * alpha }}}}
  - {charge: {{ "%.2f" % (zCYS * f) }}, hydrophobicity: !Lambda {{0.5922529084601322 * LambdaT}}, mass: 103.14, name: CYS, σ: 5.48, ε: {{ "%.4f" % epsilon }}, custom: {alpha: {{ f * alpha }}}}
  - {charge: 0.0, hydrophobicity: !Lambda {{0.3706962163690402 * LambdaT}}, mass: 114.1,  name: ASN, σ: 5.68, ε: {{ "%.4f" % epsilon }}}
  - {charge: 0.0, hydrophobicity: !Lambda {{0.3143449791669133* LambdaT}}, mass: 128.13, name: GLN, σ: 6.02, ε: {{ "%.4f" % epsilon }}}
  - {charge: 0.0, hydrophobicity: !Lambda {{0.4473142572693176 * LambdaT}}, mass: 87.08,  name: SER, σ: 5.18, ε: {{ "%.4f" % epsilon }}}
  - {charge: 0.0, hydrophobicity: !Lambda {{0.7538308115197386 * LambdaT}}, mass: 57.05,  name: GLY, σ: 4.5,  ε: {{ "%.4f" % epsilon }}}
  - {charge: 0.0, hydrophobicity: !Lambda {{0.2672387936544146 * LambdaT}}, mass: 101.11, name: THR, σ: 5.62, ε: {{ "%.4f" % epsilon }}}
  - {charge: 0.0, hydrophobicity: !Lambda {{0.1897587086430935 * LambdaT}}, mass: 71.09,  name: ALA, σ: 5.12, ε: {{ "%.4f" % epsilon }}}
  - {charge: 0.0, hydrophobicity: !Lambda {{0.5170874160398543 * LambdaT}}, mass: 131.2,  name: MET, σ: 6.18, ε: {{ "%.4f" % epsilon }}}
  - {charge: 0.0, hydrophobicity: !Lambda {{0.950628687301107 * LambdaT}},  mass: 163.18, name: TYR, σ: 6.46, ε: {{ "%.4f" % epsilon }}}
  - {charge: 0.0, hydrophobicity: !Lambda {{0.2936174211771383 * LambdaT}}, mass: 99.13,  name: VAL, σ: 5.86, ε: {{ "%.4f" % epsilon }}}
  - {charge: 0.0, hydrophobicity: !Lambda {{1.033450123574512 * LambdaT}},  mass: 186.22, name: TRP, σ: 6.78, ε: {{ "%.4f" % epsilon }}}
  - {charge: 0.0, hydrophobicity: !Lambda {{0.5548615312993875 * LambdaT}}, mass: 113.16, name: LEU, σ: 6.18, ε: {{ "%.4f" % epsilon }}}
  - {charge: 0.0, hydrophobicity: !Lambda {{0.5130398874425708 * LambdaT}}, mass: 113.16, name: ILE, σ: 6.18, ε: {{ "%.4f" % epsilon }}}
  - {charge: 0.0, hydrophobicity: !Lambda {{0.3469777523519372 * LambdaT}}, mass: 97.12,  name: PRO, σ: 5.56, ε: {{ "%.4f" % epsilon }}}
  - {charge: 0.0, hydrophobicity: !Lambda {{0.8906449355499866 * LambdaT }}, mass: 147.18, name: PHE, σ: 6.36, ε: {{ "%.4f" % epsilon }}}
  - {charge: 0.0, hydrophobicity: !Lambda {{0.5922529084601322 * LambdaT }}, mass: 103.14, name: CSS, σ: 5.48, ε: {{ "%.4f" % epsilon }}}

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
  cell: !Sphere {radius: 120.0}
  medium:
    permittivity: !Water
    temperature: {{ T }}
    salt: [!NaCl, {{ saltcon }}]
  blocks:
  - molecule: MOL1
    N: 1
    insert: !RandomCOM { filename: {{ xyz_path }}, rotate: true, directions: none, offset: [0.0, 0.0, 0.0] }
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
  selections: ["molecule MOL1", "molecule MOL2"]
  file: com_distance_{{ T }}.dat.gz
  frequency: !Every 5
- !Trajectory
  file: traj_{{ T }}.xyz
  frequency: !Every 100
- !VirtualTranslate
  selections: ["molecule MOL1"]
  molecule: MOL2
  dL: 0.05
  directions: !z
  file: "vt.dat_{{ T }}.gz"
  temperature: {{ T }}
  frequency: !Every 5


propagate:
  seed: Hardware
  criterion: Metropolis
  repeat: 500000
  collections:
  - !Stochastic
    moves:
    - !RotateMolecule
      molecule: MOL1
      dp: 1
      weight: 1.0
    - !RotateMolecule
      molecule: MOL2
      dp: 1
      weight: 1.0
    - !TranslateMolecule
      molecule: MOL2
      dp: 10
      weight: 1.0
      directions: !z

"""

if __name__ == "__main__":
    main()
