# -*- coding: utf-8 -*-

import bagpype
import glob
import numpy as np
import pandas as pd
import os

def generate_mask(size):
    mask = [[1 if i != j else 0 for j in range(size)] for i in range(size)]
    return mask

folder = '/Users/kevinmicha/Downloads/pdb_files'
file_list = [file for file in glob.glob(os.path.join(folder, '*.pdb')) if 'stripped' not in file]
std_amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

for file in file_list:
    pdb = file[-8:-4]
    myprot = bagpype.molecules.Protein()

    parser = bagpype.parsing.PDBParser(file)
    parser.parse(myprot, strip = {'res_name': ['HOH']})
    
    try:
        ggenerator = bagpype.construction.Graph_constructor()
        adjacency_matrix = ggenerator.construct_graph(myprot, exclude_backbone=False, 
            atoms_file_name=f'/Users/kevinmicha/Downloads/atoms/{pdb}.csv', 
            bonds_file_name=f'/Users/kevinmicha/Downloads/bonds/{pdb}.csv')

        # Build the membership matrix 'H'
        atoms_data = pd.read_csv(f'/Users/kevinmicha/Downloads/atoms/{pdb}.csv', delimiter=',')
        atoms_data_ca = atoms_data[atoms_data['name'] == 'CA'] # Residue iff contains CA
        atoms_data_ca = atoms_data_ca[atoms_data_ca['res_name'].isin(std_amino_acids)] # Removing non-standard amino acids
        atoms_data['residue_key'] = atoms_data['res_num'].astype(str) + '_' + atoms_data['chain'].astype(str)
        membership_matrix = np.zeros((len(atoms_data), len(atoms_data_ca)))
        
        for i in range(len(atoms_data_ca)):
            res_key = str(atoms_data_ca['res_num'].values[i]) + '_' + str(atoms_data_ca['chain'].values[i])
            indices = atoms_data[atoms_data['residue_key'] == res_key]['id']
            membership_matrix[indices, i] = 1

        atoms_data.drop(columns=['residue_key'], inplace=True) # We do not need this anymore

        # Compute L_{res}
        adj_res = np.multiply(generate_mask(len(atoms_data_ca)), membership_matrix.T @ adjacency_matrix @ membership_matrix)   
        degree_matrix = np.diag(np.sum(adj_res, axis=1))
        laplacian_res = degree_matrix - adj_res
        np.save(f'/Users/kevinmicha/Downloads/inv_laplacians_bagpype/{pdb}.npy', np.linalg.pinv(laplacian_res))

    except:
        pass
