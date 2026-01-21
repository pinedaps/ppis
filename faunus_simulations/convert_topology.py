#!/usr/bin/env python3
"""
Convert topology from Duello format to Faunus format.
Generates a complete Faunus topology file based on the template structure.
"""

import yaml
import sys
import argparse
from pathlib import Path


class CustomYAMLLoader(yaml.SafeLoader):
    """Custom YAML loader that handles custom tags from Duello format."""
    pass


def custom_tag_constructor(loader, tag_suffix, node):
    """Handle custom tags like !Lambda, !AshbaughHatch, etc."""
    if isinstance(node, yaml.ScalarNode):
        return f"!{tag_suffix.split('!')[-1]} {loader.construct_scalar(node)}"
    elif isinstance(node, yaml.MappingNode):
        return {f"!{tag_suffix.split('!')[-1]}": loader.construct_mapping(node)}
    elif isinstance(node, yaml.SequenceNode):
        return {f"!{tag_suffix.split('!')[-1]}": loader.construct_sequence(node)}
    return None


# Register custom tag constructors
CustomYAMLLoader.add_multi_constructor('!', custom_tag_constructor)


def convert_atoms_section(duello_atoms):
    """
    Convert atoms from Duello format to Faunus format.
    
    Duello format (inline):
    - {charge: -1.00, hydrophobicity: !Lambda 0, mass: 0, name: CTR, σ: 2.0, ε: 0.8368}
    
    Faunus format (expanded):
    - name: CTR
      mass: 0
      charge: -1.0
      sigma: 2.0
      epsilon: 0.8368
      hydrophobicity: 0  # For Ashbaugh-Hatch potential (lambda parameter)
    """
    faunus_atoms = []
    
    for atom in duello_atoms:
        faunus_atom = {}
        
        # Required fields - map from duello to faunus naming
        if 'name' in atom:
            faunus_atom['name'] = atom['name']
        if 'mass' in atom:
            faunus_atom['mass'] = atom['mass']
        if 'charge' in atom:
            faunus_atom['charge'] = atom['charge']
        if 'σ' in atom:  # Greek sigma to 'sigma'
            faunus_atom['sigma'] = atom['σ']
        if 'ε' in atom:  # Greek epsilon to 'epsilon'
            faunus_atom['epsilon'] = atom['ε']
        
        # Extract hydrophobicity (lambda) value from the duello format
        # It comes as a string like "!Lambda 0.5" or as a float
        if 'hydrophobicity' in atom:
            hydro_val = atom['hydrophobicity']
            # If it's a string starting with !Lambda, extract the numeric part
            if isinstance(hydro_val, str) and '!Lambda' in hydro_val:
                try:
                    faunus_atom['hydrophobicity'] = float(hydro_val.split()[-1])
                except (ValueError, IndexError):
                    faunus_atom['hydrophobicity'] = 0.0
            elif isinstance(hydro_val, (int, float)):
                faunus_atom['hydrophobicity'] = float(hydro_val)
            else:
                faunus_atom['hydrophobicity'] = 0.0
        
        if 'custom' in atom:
            faunus_atom['custom'] = atom['custom']
        
        faunus_atoms.append(faunus_atom)
    
    return faunus_atoms


def load_duello_topology(filepath):
    """Load the Duello topology YAML file."""
    with open(filepath, 'r') as f:
        return yaml.load(f, Loader=CustomYAMLLoader)


def extract_energy_config(duello_topology):
    """
    Extract energy configuration from Duello topology.
    Returns properly formatted AshbaughHatch parameters.
    """
    energy_config = duello_topology.get('system', {}).get('energy', {}).get('nonbonded', {}).get('default', [])
    
    ashbaugh_config = None
    for item in energy_config:
        if isinstance(item, dict):
            for key, val in item.items():
                if 'AshbaughHatch' in key or 'ashbaugh' in key.lower():
                    # Extract mixing and cutoff from the duello config
                    if isinstance(val, dict):
                        mixing = val.get('mixing', 'LB')
                        cutoff = val.get('cutoff', 20)
                        # Format as LB (Lorentz-Berthelot mixing rule) if arithmetic is specified
                        if mixing.lower() == 'arithmetic':
                            mixing = 'LB'
                        # Ensure cutoff is an integer if it's a whole number
                        if isinstance(cutoff, float) and cutoff.is_integer():
                            cutoff = int(cutoff)
                        ashbaugh_config = f"!AshbaughHatch {{mixing: {mixing}, cutoff: {cutoff}}}"
    
    return ashbaugh_config if ashbaugh_config else "!AshbaughHatch {mixing: LB, cutoff: 20}"


def create_faunus_topology_full(duello_topology, structure_file, sphere_radius, 
                                salt_concentration, mol2_offset_z=23.0):
    """
    Create a complete Faunus topology from Duello topology with system configuration.
    
    Args:
        duello_topology: Parsed Duello YAML topology
        structure_file: Name of the structure file (.xyz)
        sphere_radius: Radius of the spherical simulation cell
        salt_concentration: Salt concentration in mol/L
        mol2_offset_z: Z-offset for second molecule (default: 23.0)
    """
    # Extract temperature from duello file
    temperature = duello_topology.get('T', 298.15)
    
    # Extract AshbaughHatch configuration from duello
    ashbaugh_config = extract_energy_config(duello_topology)
    
    faunus_topology = {
        'atoms': convert_atoms_section(duello_topology.get('atoms', [])),
        'molecules': [
            {
                'name': 'MOL1',
                'degrees_of_freedom': 'Rigid',
                'has_com': True,
                'from_structure': structure_file
            },
            {
                'name': 'MOL2',
                'degrees_of_freedom': 'Rigid',
                'has_com': True,
                'from_structure': structure_file
            }
        ],
        'system': {
            'cell': f'!Sphere {{ radius: {sphere_radius} }}',
            'medium': {
                'permittivity': '!Water',
                'temperature': temperature,
                'salt': f'[!NaCl, {salt_concentration}]'
            },
            'blocks': [
                {
                    'molecule': 'MOL1',
                    'N': 1,
                    'insert': f'!RandomCOM {{ filename: "{structure_file}", rotate: true, directions: none, offset: [0.0, 0.0, 0.0] }}'
                },
                {
                    'molecule': 'MOL2',
                    'N': 1,
                    'insert': f'!RandomCOM {{ filename: "{structure_file}", rotate: true, directions: none, offset: [0.0, 0.0, {mol2_offset_z}] }}'
                }
            ],
            'energy': {
                'nonbonded': {
                    'default': ashbaugh_config if ashbaugh_config else {},
                    'coulomb': '!Coulomb {cutoff: 1000}'
                }
            }
        },
        'analysis': [
            '!MassCenterDistance: {molecules: [MOL1, MOL2], file: com_distance.dat.gz, frequency: !Every 100}',
            '!Trajectory: {file: traj.xyz, frequency: !Every 100}'
        ],
        'propagate': {
            'seed': 'Hardware',
            'criterion': 'Metropolis',
            'repeat': 500000,
            'collections': [
                {
                    '!Stochastic': {
                        'moves': [
                            '!RotateMolecule: {molecule: MOL1, dp: 1.0, weight: 1.0}',
                            '!RotateMolecule: {molecule: MOL2, dp: 1.0, weight: 1.0}',
                            '!TranslateMolecule: {molecule: MOL2, dp: 1.0, weight: 1.0, directions: !z}'
                        ]
                    }
                }
            ]
        }
    }
    
    return faunus_topology


def save_faunus_topology(topology, output_filepath):
    """Save the Faunus topology to a YAML file matching the exact format."""
    output_lines = []
    
    # Write atoms section
    output_lines.append('atoms:')
    for atom in topology['atoms']:
        output_lines.append('- name: ' + str(atom['name']))
        output_lines.append('  mass: ' + str(atom['mass']))
        output_lines.append('  charge: ' + str(atom['charge']))
        output_lines.append('  sigma: ' + str(atom['sigma']))
        output_lines.append('  epsilon: ' + str(atom['epsilon']))
        output_lines.append('  hydrophobicity: !Lambda ' + str(atom['hydrophobicity']))
        if 'custom' in atom:
            output_lines.append('  custom:')
            for key, val in atom['custom'].items():
                output_lines.append(f'    {key}: {val}')
    
    # Write molecules section
    output_lines.append('')
    output_lines.append('molecules:')
    for mol in topology['molecules']:
        output_lines.append('- name: ' + mol['name'])
        output_lines.append('  degrees_of_freedom: ' + mol['degrees_of_freedom'])
        output_lines.append('  has_com: ' + str(mol['has_com']).lower())
        output_lines.append('  from_structure: "' + mol['from_structure'] + '"')
    
    # Write system section
    output_lines.append('')
    output_lines.append('system:')
    output_lines.append('  cell: ' + topology['system']['cell'])
    output_lines.append('  medium:')
    output_lines.append('    permittivity: ' + topology['system']['medium']['permittivity'])
    output_lines.append('    temperature: ' + str(int(topology['system']['medium']['temperature'])))
    output_lines.append('    salt: ' + topology['system']['medium']['salt'])
    output_lines.append('  blocks:')
    for block in topology['system']['blocks']:
        output_lines.append('  - molecule: ' + block['molecule'])
        output_lines.append('    N: ' + str(block['N']))
        output_lines.append('    insert: ' + block['insert'])
    output_lines.append('  energy:')
    output_lines.append('    nonbonded:')
    output_lines.append('      default:')
    
    # Write Coulomb and AshbaughHatch in the correct order
    coulomb = topology['system']['energy']['nonbonded']['coulomb']
    ashbaugh = topology['system']['energy']['nonbonded']['default']
    
    output_lines.append('      - ' + coulomb)
    output_lines.append('      - ' + ashbaugh)
    
    # Write analysis section
    output_lines.append('')
    output_lines.append('analysis:')
    output_lines.append('- !MassCenterDistance')
    output_lines.append('  molecules: [MOL1, MOL2]')
    output_lines.append('  file: com_distance.dat.gz')
    output_lines.append('  frequency: !Every 100')
    output_lines.append('- !Trajectory')
    output_lines.append('  file: traj.xyz')
    output_lines.append('  frequency: !Every 100')
    
    # Write propagate section
    output_lines.append('')
    output_lines.append('propagate:')
    output_lines.append('  seed: Hardware')
    output_lines.append('  criterion: Metropolis')
    output_lines.append('  repeat: 500000')
    output_lines.append('  collections:')
    output_lines.append('  - !Stochastic')
    output_lines.append('    moves:')
    output_lines.append('    - !RotateMolecule')
    output_lines.append('      molecule: MOL1')
    output_lines.append('      dp: 1.0')
    output_lines.append('      weight: 1.0')
    output_lines.append('    - !RotateMolecule')
    output_lines.append('      molecule: MOL2')
    output_lines.append('      dp: 1.0')
    output_lines.append('      weight: 1.0')
    output_lines.append('    - !TranslateMolecule')
    output_lines.append('      molecule: MOL2')
    output_lines.append('      dp: 1.0')
    output_lines.append('      weight: 1.0')
    output_lines.append('      directions: !z')
    
    with open(output_filepath, 'w') as f:
        f.write('\n'.join(output_lines))


def main():
    parser = argparse.ArgumentParser(
        description='Convert Duello topology to Faunus format with system configuration.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  %(prog)s topology_duello.yaml --xyz structure.xyz --radius 125.0 --salt 0.005
  %(prog)s topology_duello.yaml -x protein.xyz -r 150 -s 0.01 -z 23.0
        '''
    )
    
    parser.add_argument(
        'duello',
        type=str,
        help='Input Duello topology file'
    )
    
    parser.add_argument(
        '-x', '--xyz',
        type=str,
        required=True,
        help='Path or filename of the .xyz structure file'
    )
    
    parser.add_argument(
        '-r', '--radius',
        type=float,
        required=True,
        help='Radius of the spherical simulation cell'
    )
    
    parser.add_argument(
        '-s', '--salt',
        type=float,
        required=True,
        help='Salt concentration in mol/L'
    )
    
    parser.add_argument(
        '-z', '--offset-z',
        type=float,
        default=23.0,
        help='Z-offset for second molecule (default: 23.0)'
    )
    
    args = parser.parse_args()
    
    # File paths
    duello_file = Path(args.duello)
    
    if not duello_file.exists():
        print(f"Error: {duello_file} not found")
        sys.exit(1)
    
    # Construct output filename: remove extension and add _faunus.yaml
    output_file = duello_file.stem + '_faunus.yaml'
    output_file = Path(output_file)
    
    print(f"Reading Duello topology from: {duello_file}")
    duello_topology = load_duello_topology(duello_file)
    
    temperature = duello_topology.get('T', 298.15)
    print(f"Temperature from Duello file: {temperature} K")
    
    print(f"\nConverting to Faunus format with:")
    print(f"  Structure file: {args.xyz}")
    print(f"  Sphere radius: {args.radius}")
    print(f"  Salt concentration: {args.salt} mol/L")
    print(f"  MOL2 Z-offset: {args.offset_z}")
    
    faunus_topology = create_faunus_topology_full(
        duello_topology,
        args.xyz,
        args.radius,
        args.salt,
        args.offset_z
    )
    
    print(f"\nSaving Faunus topology to: {output_file}")
    save_faunus_topology(faunus_topology, output_file)
    
    print("✓ Conversion complete!")
    print(f"✓ Number of atoms converted: {len(faunus_topology['atoms'])}")


if __name__ == '__main__':
    main()
