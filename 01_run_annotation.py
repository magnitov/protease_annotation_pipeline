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
import re
from urllib import parse, request
import warnings
warnings.filterwarnings('ignore')

########################################################################################
# Functions
########################################################################################

def find_ligand(ligand, rasmol_file):
    """
    Find ligand chain and ligand residue numbers

    Arguments:
    ligand -- a list of ligand residue names
    rasmol_file -- rasmol.dfn file from the PDBsum database
    
    Output:
    Ligand chain and a list of residues as present in the PDB file
    """

    ligand_chain = ""
    residues_numbers = []
    length_ligand = len(ligand)
    for row in rasmol_file:
        if row[0] == "L" and all(ligand_iter in row for ligand_iter in ligand) and ':' in row:
            row = row.replace(",", "").replace("[", "").replace("]", "").split()
            del row[0]
            ligand_str = "".join(ligand)
            ligand_str_row = "".join((iter[0:3] for iter in row))
            if ligand_str_row == ligand_str:
                for i in range(len(row)):
                    if row[i][-3].isalpha():
                        row[i] = row[i][:-3] + row[i][-2:]
                ligand_chain = row[0][-1]
                residues_numbers = [int(row[i][3:-2]) for i in range(length_ligand)]
                break
        else: # the cases when ligand chain is missing from the rasmol file
            ligand_chain = ''
            residues_numbers = np.nan
    return (residues_numbers, ligand_chain)

##################################################################################

def find_prot_chain_bound_to_ligand(ligand, residue_number, ligand_chain, grow_file):
    """
    Find protein chains bound to ligand

    Arguments:
    ligand -- a list of ligand residue names
    residue_number -- a list of ligand residue numbers on its PDB chain
    ligand_chain -- ligand PDB chain
    grow_file -- grow.out file from the PDBsum database
    
    Output:
    A list of protein chains bound to ligand chain
    """

    chains_bound_to_ligand = []
    for row in grow_file:
        if row[14:17] == row[14:17].upper() and row[54] == ligand_chain\
           and (row[44:47], int(row[49:53])) in zip(ligand, residue_number):
            if row[6] not in chains_bound_to_ligand:
                chains_bound_to_ligand.append(row[6])
    return (chains_bound_to_ligand)

##################################################################################

def find_prot_residues_bound_to_ligand(ligand, residue_number, ligand_chain, grow_file):
    """
    Find protein residues bound to ligand

    Arguments:
    ligand -- a list of ligand residue names
    residue_number -- a list of ligand residue numbers on its PDB chain
    ligand_chain -- ligand PDB chain
    grow_file -- grow.out file from the PDBsum database
    
    Output:
    A list of protein chains bound to ligand chain
    """

    residues_bound_to_ligand = []
    for row in grow_file:
        if row[14:17] == row[14:17].upper() and row[54] == ligand_chain\
           and (row[44:47], int(row[49:53])) in zip(ligand, residue_number):
            if (row[5:17].replace(' ', '') not in residues_bound_to_ligand) and len(re.findall('\d+', row[5:17].replace(' ', '')[-3:])) == 0:
                residues_bound_to_ligand.append(row[5:17].replace(' ', ''))
    return (residues_bound_to_ligand)

##################################################################################

def find_chains(data):
    """
    Perform the search for ligand and protein chain and residues as described above

    Arguments:
    data -- a data frame with PDB code and ligand
    
    Output:
    An updated data frame with information about the ligand and protein added
    """

    # Length of the ligand
    data['Ligand_size'] = [len(x.split('-')) for x in data['Ligand'].values]

    ligand_chain, ligand_residues, bound_chains, bound_residues = [], [], [], []
    for j, PDB, ligand in zip(data.index, data['PDB'], [x.split('-') for x in data['Ligand'].values]):
        # Find ligand position
        with open('./annotation/rasmol/' + PDB + '_rasmol.dfn', 'r') as file_ligand:
            ligand_residues_tmp, ligand_chain_tmp = find_ligand(ligand, file_ligand)
        ligand_chain.append(ligand_chain_tmp)
        ligand_residues.append(ligand_residues_tmp)

        # Find protein chains bound to ligand
        with open('./annotation/grow/' + PDB + '_grow.out', 'r') as grow_file:
            bound_chains_tmp = find_prot_chain_bound_to_ligand(ligand, ligand_residues_tmp, ligand_chain_tmp, grow_file)
            if len(bound_chains_tmp) == 0:
                bound_chains_tmp = ''
        bound_chains.append(bound_chains_tmp)

        # Find protein residues bound to ligand
        with open('./annotation/grow/' + PDB + '_grow.out', "r") as grow_file:
            bound_residues_tmp = find_prot_residues_bound_to_ligand(ligand, ligand_residues_tmp, ligand_chain_tmp, grow_file)
            if len(bound_residues_tmp) == 0:
                bound_residues_tmp = ''
            else:
                bound_residues_tmp = list(np.unique(bound_residues_tmp))
        bound_residues.append(bound_residues_tmp)

    data['Ligand_chain'] = ligand_chain
    data['Ligand_residues_number'] = ligand_residues
    data['Protein_chains'] = bound_chains
    data['Protein_residues'] = bound_residues
    return(data)

##################################################################################

def find_ids(data):
    """
    Annotate EC number, UniProt ID and MEROPS ID

    Arguments:
    data -- a data frame with PDB code and information about ligand and protein
    
    Output:
    An updated data frame with IDs added
    """

    # Retrieving EC and UniProt IDs from PDBsum annotation file
    ec = pd.read_table('pdb_chain_sp_ec', header = None, sep = '\s', engine = 'python')
    for i in [2, 4, 6, 8, 10]:
        ec[i] = [x.replace('"', '') for x in ec[i].values]

    ecs, uniprots = [], []
    for i in range(0, len(data)):
        ec_temp, uniprot_temp = [], []
        chains_take = str(data['Protein_chains'].values[i]).replace('[', '')\
                         .replace(']', '').replace("'", "").split(', ')
        for j in range(0, len(ec[ec[2] == data['PDB'].values[i]])):
            if ec[ec[2] == data['PDB'].values[i]][4].values[j] in chains_take:
                ec_temp.append(ec[ec[2] == data['PDB'].values[i]][10].values[j])
                uniprot_temp.append(ec[ec[2] == data['PDB'].values[i]][6].values[j])

        ec_temp = ', '.join(np.unique(ec_temp))
        if ec_temp == '':
            ec_temp = '-'
        if ec_temp[0] == ',':
            ec_temp = ec_temp[2:]
        ecs.append(ec_temp)          

        uniprot_temp = ', '.join(np.unique(uniprot_temp))
        if uniprot_temp == '':
            uniprot_temp = '-'
        if uniprot_temp[0] == ',':
            uniprot_temp = uniprot_temp[2:]
        uniprots.append(uniprot_temp)

    # Retrieving MEROPS IDs using UniProt REST API
    url = 'https://www.uniprot.org/uploadlists/'
    params = {'from': 'ACC+ID', 'to': 'MEROPS_ID', 'format': 'tab', 'query': ', '.join(uniprots)}
    convertion = parse.urlencode(params)
    convertion = convertion.encode('utf-8')
    req = request.Request(url, convertion)
    with request.urlopen(req) as f:
       merops_temp = f.read()
    merops_temp = merops_temp.decode('utf-8').split()[2:]

    merops_part = [x for (i, x) in enumerate(merops_temp) if i in np.arange(1, len(merops_temp), 2)]
    uniprot_part = [x for (i, x) in enumerate(merops_temp) if i in np.arange(0, len(merops_temp), 2)]
    uniprot_merops_convertion = dict(zip(uniprot_part, merops_part))
          
    merops = []
    for uniprot_temp in uniprots:
        if uniprot_temp == '-':
            merops.append('-')
        else:
            merops_temp = []
            for uniprot_id in uniprot_temp.split(', '):
                if uniprot_id in uniprot_part:
                    merops_temp.append(uniprot_merops_convertion[uniprot_id])
            merops_temp = ', '.join(merops_temp)
            if merops_temp == '':
                merops_temp = '-'
            merops.append(merops_temp)

    data['EC_number'] = ecs
    data['UniProt_ID'] = uniprots
    data['MEROPS_ID'] = merops
    return(data)

##################################################################################

def find_catalytic_residues(data):
    """
    Find catalytic residues

    Arguments:
    data -- a data frame with PDB code and information about ligand and protein
    
    Output:
    An updated data frame with catalytic residues added
    """
    with open('homologues_residues.json', 'r') as f:
        residues = f.readlines() 
    residues = residues[0].split('}, {')
    residues = [x.split('[{')[1] if len(x.split('[{')) > 1 else x for x in residues]

    catalytic_residues = []
    for (pdb, chains) in zip(data['PDB'].values, data['Protein_chains'].values):
        cat_res_tmp = []
        for res in residues:
            if pdb in res:
                res = [x.replace('"', '').replace(',', '') for (i, x) in enumerate(res.split()) if i in [1, 5, 9]]
                if res[1] != 'null' and res[2] in chains:
                     cat_res_tmp.append(res[2]+res[1]+res[0].upper())
        cat_res_tmp = ', '.join(np.unique(cat_res_tmp))
        if cat_res_tmp == '':
            cat_res_tmp = '-'
        catalytic_residues.append(cat_res_tmp)

    data['Catalytic_residues'] = catalytic_residues
    return(data)


########################################################################################
# Main
########################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('--input', default = 'sample_data.csv', type = str, required = True,
                    help = 'Dataset with protease-ligand pairs to annotate')
parser.add_argument('--output', default = 'annotated_dataset.csv', type = str, required = False, 
                    help='Path to save the annotated dataset')
args = parser.parse_args()

# Read the dataset
data = pd.read_table(args.input, sep = ',', header = None, names = ['PDB', 'Ligand'])
# Find ligand and protein chains
data = find_chains(data)
# Perform annotation of IDs
data = find_ids(data)
# Annotate catalytic residues
data = find_catalytic_residues(data)
# Save annotated dataset
data.to_csv(args.output, sep = '\t', index = 0)
