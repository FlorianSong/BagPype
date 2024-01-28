# -*- coding: utf-8 -*-

import bagpype
import glob
import os

folder = 'pdbs/'
file_list = glob.glob(os.path.join(folder, '*.pdb'))

for file in file_list:
    myprot = bagpype.molecules.Protein()

    parser = bagpype.parsing.PDBParser(file)
    parser.parse(myprot, strip = {'res_name': ['HOH']})

    ggenerator = bagpype.construction.Graph_constructor()
    ggenerator.construct_graph(myprot)