#!/usr/bin/env python

#Written by Thomas H. Tugwell
#
#This program generates modredundant files for Gaussian calculations over a specified change in bond distances. 
#It takes command-line arguments specifying the input file, atom indices for the bond (atom_index_1 is fixed), extension range(min,max) in angstroms, step size(number of files), 
# output prefix, and a list of other atoms to move along the same vector.

#usage: modred_maker_improved.py <input_file> <atom_index_1> <atom_index_2> <extension_min> <extension_max> <step_size> <output_prefix> <fragment_list>
#example input: modred_maker_improved.py test.gjf 23 19 0 5 10 mod '[15,16,17,18,20,21,22,32]'

import sys
import os
import numpy as np

print("\nNOTE: MAKE SURE YOU HAVE CHANGED BASIS SET, CHARGE, MULT., ETC... IN THIS SCRIPT BEFORE USE.\n")

def parse_gaussian_input(input_file):
    """
    Parses the Gaussian input file and extracts atomic coordinates.

    Args:
        input_file (str): Path to the Gaussian input file.

    Returns:
        atoms (list): List of atomic symbols.
        coords (list): List of atomic coordinates.
    """
    # Read lines from the input file
    with open(input_file, 'r') as f:
        lines = f.readlines()

    atoms = []
    coords = []

    for line in lines:
        line = line.strip()
        if line.startswith('#') or line.startswith('%') or len(line) == 0:
            continue

        split_line = line.split()
        if len(split_line) < 4:
            continue

        try:
            atom = split_line[0]
            x = float(split_line[1])
            y = float(split_line[2])
            z = float(split_line[3])
            atoms.append(atom)
            coords.append([x, y, z])
        except ValueError:
            continue

    return atoms, coords

def generate_gaussian_input(route_card, title_card, charge, multiplicity, atoms, coords):
    """
    Generates the Gaussian input string with specified route card, title card, charge, multiplicity,
    atoms, and coordinates.

    Args:
        route_card (str): Route card for the Gaussian calculation.
        title_card (str): Title card for the Gaussian calculation.
        charge (int): Charge of the system.
        multiplicity (int): Multiplicity of the system.
        atoms (list): List of atomic symbols.
        coords (list): List of atomic coordinates.

    Returns:
        input_string (str): Gaussian input string.
    """
    
    input_string = "%nproc=16\n"
    input_string += "%mem=32GB\n"
    input_string += f"{route_card}\n"
    input_string += "\n"
    input_string += f"{title_card}\n"
    input_string += "\n"
    input_string += f"{charge} {multiplicity}\n"

    for atom, coord in zip(atoms, coords):
        x, y, z = coord
        input_string += f"{atom:2s} {x:12.6f} {y:12.6f} {z:12.6f}\n"
    #input_string += "\nB 16 5 F" #Only use this line if you need to constrain additional bonds    
    input_string += f"\nB {bond_indices[0]} {bond_indices[1]} F\n\n"
    
    #basis1_atoms = "C H F O Br 0"             # Use these lines to add infor if using mixed basis set
    #basis1 = "6-31+G(d,p)"
    #basis2_atoms = "Ni 0"
    #basis2 = "SDD"
    #input_string += f"{basis1_atoms}\n"    
    #input_string += f"{basis1}\n"
    #input_string += f"****\n"
    #input_string += f"{basis2_atoms}\n"
    #input_string += f"{basis2}\n"
    #input_string += f"****\n\n"
    #input_string += f"{basis2_atoms}\n"
    #input_string += f"{basis2}\n\n"

    return input_string

def main(input_file, bond_indices, extension_range, step_size, output_prefix,fragment_list):
    """
    Main function that generates modredundant files for Gaussian with specified bond extensions.

    Args:
        input_file (str): Path to the Gaussian input file.
        bond_indices (list): List of two integers representing the indices of the atoms involved in the bond.
        extension_range (list): List of two floats representing the minimum and maximum bond extensions.
        step_size (int): Number of steps in the bond extension.
        output_prefix (str): Prefix for the output file names.

    Returns:
        output_files (list): List of output file names.
    """
    atoms, coords = parse_gaussian_input(input_file)

    atom1_index = bond_indices[0] - 1
    atom2_index = bond_indices[1] - 1

    fragment_index = [x - 1 for x in fragment_list]

    bond_length = np.linalg.norm(np.array(coords[atom1_index]) - np.array(coords[atom2_index]))

    output_files = []

    for i in range(step_size):
        factor = float(i) / (step_size - 1)
        delta_length = extension_range[0] + factor * (extension_range[1] - extension_range[0])
        new_bond_length = bond_length + delta_length

        new_coords = coords.copy()
        new_coords[atom2_index] = np.array(coords[atom1_index]) + (np.array(coords[atom2_index]) - np.array(coords[atom1_index])) * (new_bond_length / bond_length)

        change_vector = new_coords[atom2_index] - coords[atom2_index] #calculates the vector that atom2 moved from atom1 (used for moving rest of fragment)
        
        for idx in fragment_index:
            new_fragment_coords = [
            coords[idx][0] + change_vector[0],
            coords[idx][1] + change_vector[1],
            coords[idx][2] + change_vector[2]
            ]
            new_coords[idx] = new_fragment_coords
        

        output_file = f"{atoms[atom1_index]}-{atoms[atom2_index]}_{new_bond_length:.1f}.gjf"
        #output_file = f"{output_prefix}_{i}.gjf"
        output_files.append(output_file)

        route_card = "# opt=modredundant freq 6-31g(d) m062x"
        title_card = "# modreds"
        charge = -1
        multiplicity = 1
        
        gaussian_input = generate_gaussian_input(route_card, title_card, charge, multiplicity, atoms, new_coords)

        with open(output_file, 'w') as f:
            f.write(gaussian_input)

    return output_files


if __name__ == '__main__':
    if len(sys.argv) < 8:
        print("Usage: python modred_maker.py <input_file> <atom_index_1> <atom_index_2> <extension_min> <extension_max> <step_size> <output_prefix> [<fragment_list>]")
        print("Example input: modred_maker.py gaussian_input.gjf 23 19 0 5 10 mod '[15,16,17,18,20,21,22,32]")
        sys.exit(1)

    input_file = sys.argv[1]
    bond_indices = [int(sys.argv[2]), int(sys.argv[3])]
    extension_range = [float(sys.argv[4]), float(sys.argv[5])]
    step_size = int(sys.argv[6])
    output_prefix = sys.argv[7]
    
    fragment_index = []
    if len(sys.argv) > 8:
        fragment_list = sys.argv[8].strip('[]').split(',')
        fragment_index = [int(item) for item in fragment_list]

    output_files = main(input_file, bond_indices, extension_range, step_size, output_prefix, fragment_index)

    print("Output files created:")
    for output_file in output_files:
        print(output_file)
