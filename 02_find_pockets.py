#!/usr/bin/python3

########################################################################################
# Authorship
########################################################################################

__author__ = "Mikhail Magnitov"
__email__ = "mikhail.magnitov@phystech.edu"

########################################################################################
# Modules to import
########################################################################################

import argparse
import pandas as pd
import numpy as np
import os
import collections
import warnings
warnings.filterwarnings('ignore')

########################################################################################
# Functions
########################################################################################

def get_protein_and_ligand_coordinates(dataset, pdb_code):
    """
    Get reference protein and ligand coordinates from the superimposed PDB file

    Arguments:
    dataset -- annotated dataset containing information about protein-ligand interactions
    pdb_code -- PDB code of protein-ligand pair in the dataset
    
    Output:
    Two data frames with ligand and protein atoms
    """

    # Read superimposed coordinates
    protein_coords_sup_df = []
    protein_coords_sup = open('./pockets/protein_ligand_' + pdb_code + '.pdb', 'r').readlines()
    for j in range(0, len(protein_coords_sup)):
        protein_coords_sup_df.append([protein_coords_sup[j][12:16].replace(' ', ''), protein_coords_sup[j][17:21].replace(' ', ''),
                                      protein_coords_sup[j][21:22].replace(' ', ''), protein_coords_sup[j][22:27].replace(' ', ''),
                                      protein_coords_sup[j][30:38].replace(' ', ''), protein_coords_sup[j][38:46].replace(' ', ''),
                                      protein_coords_sup[j][46:54].replace(' ', '')])
    protein_coords_sup = pd.DataFrame(protein_coords_sup_df, columns = ['atom', 'residue', 'chain', 'res_number', 'x', 'y', 'z'])
    protein_coords_sup = protein_coords_sup[:-1] # chop last line with END, which is empty after previous command
    protein_coords_sup['x'] = [float(x) for x in protein_coords_sup['x'].values]
    protein_coords_sup['y'] = [float(x) for x in protein_coords_sup['y'].values]
    protein_coords_sup['z'] = [float(x) for x in protein_coords_sup['z'].values]

    # Split protein and ligand
    ligand_chain = dataset[dataset['PDB'] == pdb_code]['Ligand_chain'].values[0]
    ligand_coords_sup = protein_coords_sup[protein_coords_sup['chain'] == ligand_chain]
    protein_coords_sup = protein_coords_sup[protein_coords_sup['chain'] != ligand_chain]
    
    return(ligand_coords_sup, protein_coords_sup)

##################################################################################

def get_reference_pocket(ligand_coordinates, protein_coordinates, max_distance):
    """
    Get pocket residues on the protein structure

    Arguments:
    ligand_coordinates -- data frame with ligand atoms from PDB
    protein_coordinates -- data frame with protein atoms from PDB
    max_distance -- distance to search for pocket residues
    
    Output:
    A list of pocket residues for each position
    """
    
    # Make a dictionary with aminoacids and their synthetic analogues to avoid cases,
    # when HOH and other small ions around the ligand would be mapped as pocket residues
    aminoacids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE',
                  'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
    comp = open('components.cif', 'r').read().split('data_')[1:]
    for i in range(0, len(comp)):
        comp[i] = comp[i].split('#')[1]
        comp[i] = comp[i].split('_chem_comp.id')[1].split('\n')[0].replace(' ', '') + ' '\
                  +  comp[i].split('_chem_comp.mon_nstd_parent_comp_id')[1].split('\n')[0].replace(' ', '')
    codes = []
    for i in range(0, len(comp)):
        if comp[i][3:] != ' ':
            check = comp[i].split()[1].replace('"', '').split(',')
            for k in range(0, len(check)):
                if check[k] in aminoacids:
                    codes.append(comp[i].split()[0])
                    break
    aminoacids = aminoacids + codes

    # Take protein residues that have at least 1 atom within a pre-defined distance from ligand residues
    done_residues = []
    pocket = []
    for res_num in np.unique(ligand_coordinates['res_number'].values):
        if res_num not in done_residues:
            done_residues.append(res_num)
            res_coords = ligand_coordinates[ligand_coordinates['res_number'] == res_num]
            pocket_residues_temp, pocket_residues_numbers_temp = [], []
            # Calculate distances between each atom in ligand and in protein, then take those that are closer than max_distance
            for j in range(0, len(res_coords)):
                for k in range(0, len(protein_coordinates)):
                    if (np.sqrt((res_coords['x'].values[j] - protein_coordinates['x'].values[k])**2 +\
                                (res_coords['y'].values[j] - protein_coordinates['y'].values[k])**2 +\
                                (res_coords['z'].values[j] - protein_coordinates['z'].values[k])**2) < max_distance)\
                                 and protein_coordinates['residue'].values[k] in aminoacids:
                        pocket_residues_temp.append(protein_coordinates['residue'].values[k] + protein_coordinates['res_number'].values[k])
                        pocket_residues_numbers_temp.append(int(protein_coordinates['res_number'].values[k]))
            sort_rule = dict(zip(pocket_residues_temp, pocket_residues_numbers_temp))
            pocket.append(sorted(list(np.unique(pocket_residues_temp)), key=sort_rule.__getitem__))

    return(pocket)

##################################################################################

def align_ligand(target_ligand_coordinates, reference_ligand_coordinates, max_distance, initial_ligand):
    """
    Align target ligand to reference ligand. 
    The algorithm is the following: we take each residue of reference ligand and search for atoms in target 
    ligand that are within a certain distance from reference ligand atoms. We count how many times each 
    residue of target was found around certain residue in reference. Then we take the most frequent 
    reference residue on each position and assign it the same position as in treference.

    Arguments:
    target_ligand_coordinates -- data frame with target ligand atoms from PDB
    reference_ligand_coordinates -- data frame with reference ligand atoms from PDB
    max_distance -- distance to project target reference residues on ligand residues
    initial_ligand -- ligand sequence
    
    Output:
    Data frame with target ligand alignment into S4-S4' pockets
    """

    pocket_names = ["S4", "S3", "S2", "S1", "S1'", "S2'", "S3'", "S4'"]
    done_residues = []
    ligand_residue_names, ligand_residue_counts, pocket_aligned_name = [], [], []
    for res_num in np.unique(reference_ligand_coordinates['res_number'].values):
        if res_num not in done_residues:
            done_residues.append(res_num)
            res_coords = reference_ligand_coordinates[reference_ligand_coordinates['res_number'] == res_num] # reference residue atoms
            close_ligand_residues = []
            for k in range(0, len(res_coords)):
                for m in range(0, len(target_ligand_coordinates)):
                    # Checking that atom from target is closer than cutoff
                    if (np.sqrt((res_coords['x'].values[k] - target_ligand_coordinates['x'].values[m])**2 +\
                                (res_coords['y'].values[k] - target_ligand_coordinates['y'].values[m])**2 +\
                                (res_coords['z'].values[k] - target_ligand_coordinates['z'].values[m])**2) < max_distance):
                        if target_ligand_coordinates['residue'].values[m] not in ['ACE', '0QE', 'CF0', 'NH2', '0HQ', 'POL']:
                            close_ligand_residues.append(target_ligand_coordinates['residue'].values[m] + target_ligand_coordinates['res_number'].values[m])
            # If there are any found target residues within cutoff distance, we take most common residue found in this pocket
            if len(collections.Counter(close_ligand_residues).most_common(1)) > 0:
                ligand_residue_names.append(collections.Counter(close_ligand_residues).most_common(1)[0][0])
                ligand_residue_counts.append(collections.Counter(close_ligand_residues).most_common(1)[0][1])
            # If there are no target residues within the distance, append zero and empty string
            else:
                ligand_residue_names.append('')
                ligand_residue_counts.append(0)

    alignment = pd.DataFrame({'Residues' : ligand_residue_names, 'Counts' : ligand_residue_counts, 'Pocket' : pocket_names}, 
                             columns = ['Residues', 'Counts', 'Pocket'])
    # Sometimes terminal residues are aligned on two positions, the actual one, and the one nearby. We want to only take the position where it has
    # most number of atoms aligned. We look for incorrect position (taking into account that alignment can't have gaps) and drop it.
    for i in range(0, len(alignment)):
        if len(alignment[alignment['Residues'] == alignment['Residues'].values[i]]) > 1:
            alignment_temp = alignment[alignment['Residues'] == alignment['Residues'].values[i]]
            index = int(alignment_temp.loc[alignment_temp['Counts'].idxmin()].name)
            alignment.at[index, 'Residues'] = ''
            alignment.at[index, 'Counts'] = 0

    # Sometimes ligand is partially aligned, e.g. ligand consists of 1-2-3-4-5 residues and 1-2 are aligned to 1'-2' in reference,
    # and then ligand makes a turn and reference position 3' alignes to 5 in target. We want to only take consistent parts of ligand
    # Therefore, we need to find this outlier position in target and drop it. We search for it using position numbers of residues.
    if len(alignment[alignment['Residues'] != '']) <= 2:
        alignment['Counts'] = [0, 0, 0, 0, 0, 0, 0, 0]
        alignment['Residues'] = ['', '', '', '', '', '', '', '']

    if len(alignment[alignment['Residues'] != '']) > 2:
        while len(np.unique(np.array([int(x[3:]) for x in alignment[alignment['Residues'] != '']['Residues'].values][:-1]) -\
              np.array([int(x[3:]) for x in alignment[alignment['Residues'] != '']['Residues'].values][1:]))) != 1:
            alignment_filtered = alignment
            for i in range(0, len(alignment[alignment['Residues'] != ''])-1):
                dist = int(alignment[alignment['Residues'] != '']['Residues'].values[i][3:]) - int(alignment[alignment['Residues'] != '']['Residues'].values[i+1][3:])
                # If the difference between position numbers is +-1, then these residues are actually following each other in original ligand
                if dist == -1 or dist == 1:
                    pass
                # If difference is more than +-1, then these two residues have other residues in between, and we need to drop one of them
                elif dist < 0:
                    drop = alignment[alignment['Residues'] != '']['Residues'].values[i+1]
                    index = alignment[alignment['Residues'] == drop].index[0]
                    alignment_filtered.at[index, 'Residues'] = ''
                    alignment_filtered.at[index, 'Counts'] = 0
                    break
                else:
                    drop = alignment[alignment['Residues'] != '']['Residues'].values[i+1]
                    index = alignment[alignment['Residues'] == drop].index[0]
                    alignment_filtered.at[index, 'Residues'] = ''
                    alignment_filtered.at[index, 'Counts'] = 0
                    break
            alignment = alignment_filtered
            
            # If after the above filtration two residues with gap between left, we drop them
            if len(alignment[alignment['Residues'] != '']) == 2:
                 dist = int(alignment[alignment['Residues'] != '']['Residues'].values[0][3:]) - int(alignment[alignment['Residues'] != '']['Residues'].values[1][3:])
                 if dist != 1 or dist != -1:
                     alignment['Counts'] = [0, 0, 0, 0, 0, 0, 0, 0]
                     alignment['Residues'] = ['', '', '', '', '', '', '', '']
                     break

    alignment = alignment.drop(['Counts'], axis = 1)

    # Sometimes after these operations there are gaps left in the alignment that we need to fill (e.g. residues 1-2-3 are in 
    # S3-S2-S1 pockets and the terminal residue 4 is in pocket S2', but since we can't have gaps, it should be moved to S1').
    # We condider larger fraction of consistent residues to be aligned properly, and move only smaller part.
    if len(alignment[alignment['Residues'] != '']) > 0:
        alignment_indicators, alignment_indexing = [], []

        # Make a list of indicators aligned/not and a list of real indexes for list
        for i in range(0, len(alignment)):
            alignment_indexing.append(i)
            if alignment['Residues'].values[i] != '':
                alignment_indicators.append(1)
            else:
                alignment_indicators.append(0)

        # Slice only part with aligned residues in both indicators list and indexes
        while alignment_indicators[0] != 1:
            alignment_indicators = alignment_indicators[1:]
            alignment_indexing = alignment_indexing[1:]
        while alignment_indicators[-1] != 1:
            alignment_indicators = alignment_indicators[:-1]
            alignment_indexing = alignment_indexing[:-1]

        # Determine gap position and its index in full alignment
        if 0 in alignment_indicators:
            index = alignment_indicators.index(0) + 1
            real_index = alignment_indexing[index] - 1
            # Determine the part of substrate that we should move in annotation
            if index > len(alignment_indicators) / 2:
                # If we are moving the second half, we should move residues backward
                for j in range(real_index, len(alignment)-1):
                    alignment['Residues'].values[j] = alignment['Residues'].values[j+1]
            else:
                # If we are moving the first half, we should move residues forward
                for j in range(real_index, 0, -1):
                    alignment['Residues'].values[j] = alignment['Residues'].values[j-1]
                alignment['Residues'].values[0] = ''

    return(alignment)

##################################################################################

def process_alignment(alignment, initial_ligand):
    """
    Process the aligned ligand to determine whether it has been truncated or reversed

    Arguments:
    alignment -- data frame with ligand alignment into S4-S4' pockets
    initial_ligand -- ligand sequence
    
    Output:
    Truncated ligand sequence, numbering and pockets
    """

    # Change the format of extracted alignment data
    alignment = alignment[alignment['Residues'] != '']
    if '-'.join([x[:3] for x in alignment['Residues'].values]) in  initial_ligand:
        truncated_ligand = '-'.join([x[:3] for x in alignment['Residues'].values])
        truncated_numbering = [int(x[3:]) for x in alignment['Residues'].values]
        pockets = list(alignment['Pocket'].values)

    # Sometimes ligand in dataset is reversed (direction is not the same as it is in active site)
    # Here we detect these cases and reverse ligand, numbering and pockets
    else:
        truncated_ligand = '-'.join([x[:3] for x in alignment['Residues'].values][::-1])
        truncated_numbering = [int(x[3:]) for x in alignment['Residues'].values][::-1]
        pockets = list(alignment['Pocket'].values)[::-1]

    return(truncated_ligand, truncated_numbering, pockets)

########################################################################################
# Main
########################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('--input', default = 'sample_data.csv', type = str, required = True,
                    help = 'Dataset with protease-ligand pairs to annotate')
parser.add_argument('--ref', default = None, type = str, required = True,
                    help = 'Reference protease-ligand structure')
parser.add_argument('--extra_ref', default = False, type = str, required = False,
                    help = 'Extra reference protease-ligand structure for S1 pocket refinement')
parser.add_argument('--radius', default = 4.5, type = float, required = True,
                    help = 'Distance to search for pocket residues')
parser.add_argument('--cutoff', default = 2, type = int, required = True,
                    help = 'Distance to project target reference residues on ligand residues')
parser.add_argument('--output', default = 'annotated_pockets.csv', type = str, required = False, 
                    help='Path to save the dataset with annotated pockets')
args = parser.parse_args()

# Read the dataset
data = pd.read_table(args.input, sep = ',', header = 0)
# Get pocket search parameters
reference_structure = args.ref
radius = args.radius
cutoff = args.cutoff

# Taking reference protein and ligand coordinates
reference_ligand_coordinates, reference_protein_coordinates = get_protein_and_ligand_coordinates(data, reference_structure)
# Drop one reference ligand residue that doesn't fit in S4-S4' pockets
# NOTE: If you use other reference ligand, please modify the next line accordingly
reference_ligand_coordinates = reference_ligand_coordinates[reference_ligand_coordinates['res_number'] != '9']
# Find pocket residues on the reference protease
reference_pocket = get_reference_pocket(reference_ligand_coordinates, reference_protein_coordinates, radius)
# Refine the reference pocket using other structures that have different ligand residues at P1 position by placing them into the reference protein coordinates
if args.extra_ref:
    for pdb in args.extra_ref.split(','): 
        extra_reference_ligand_coordinates, extra_reference_protein_coordinates = get_protein_and_ligand_coordinates(data, pdb)
        refined_residues = get_reference_pocket(extra_reference_ligand_coordinates, reference_protein_coordinates, radius)[4]
        reference_pocket[3] = np.unique(np.hstack((reference_pocket[3], refined_residues)))
# Print obtained reference protease pocket residues
print('Pocket residues annotated:')
for (name, pocket) in zip(["S4", "S3", "S2", "S1", "S1'", "S2'", "S3'", "S4'"], reference_pocket):
    print (name + ':', ', '.join(np.sort(pocket)))

# Mapping ligands with pockets
truncated_ligand_all, truncated_numbering_all, truncated_pockets_all, target_pocket_residues_all = [], [], [], []
for (pdb, ligand) in zip(data['PDB'].values, data['Ligand'].values):
    # Coordinates of target ligand and protein
    target_ligand_coordinates, target_protein_coordinates = get_protein_and_ligand_coordinates(data, pdb)
    # Aligning target ligand to reference ligand          
    alignment = align_ligand(target_ligand_coordinates, reference_ligand_coordinates, cutoff, ligand)
    # Processing alignment data to get truncated ligand data
    truncated_ligand, truncated_numbering, truncated_pockets = process_alignment(alignment, ligand)
    
    truncated_ligand_all.append(truncated_ligand)
    truncated_numbering_all.append(truncated_numbering)
    truncated_pockets_all.append(truncated_pockets)

data['Truncated_ligand'] = truncated_ligand_all
data['Truncated_numbering'] = truncated_numbering_all
data['Truncated_pocket_names'] = truncated_pockets_all

# Save annotated dataset
data.to_csv(args.output, sep = '\t', index = 0)
