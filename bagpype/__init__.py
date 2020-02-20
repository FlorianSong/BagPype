import bagpype.molecules
import bagpype.parameters

import bagpype.energies
import bagpype.construction
import bagpype.parsing


head_string = u"\
    \n\
 |\u203E|                Welcome to              \n\
 | |__   __ _  __ _ _ __  _   _ _ __   ___  \n\
 | '_ \ / _` |/ _` | '_ \| | | | '_ \ / _ \ \n\
 | |_) | (_| | (_| | |_) | |_| | |_) |  __/ \n\
 |_.__/ \__,_|\__, | .__/ \__, | .__/ \___| \n\
               __/ | |     __/ | |          \n\
              |___/|_|    |___/|_|          \n\
Biochemical, Atomistic Graph construction   \n\
software in PYthon for Proteins, Etc.       \n\
(C) Yaliraki group @ Imperial College London\n\
"

print(head_string)


################################################

def standard_run(pdb_file):

    myprot = bagpype.molecules.Protein()

    parser = bagpype.parsing.PDBParser(pdb_file)
    parser.parse(myprot, strip = {'res_name': ['HOH']})

    ggenerator = bagpype.construction.Graph_constructor()
    ggenerator.construct_graph(myprot)

    return(myprot)