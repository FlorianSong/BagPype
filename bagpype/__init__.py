import bagpype.molecules
import bagpype.settings
import bagpype.parameters

import bagpype.energies
import bagpype.construction
import bagpype.parsing

print(bagpype.settings.head_string)


################################################

def standard_run(pdb_file):

    myprot = bagpype.molecules.Protein()

    parser = bagpype.parsing.PDBParser(pdb_file)
    parser.parse(myprot, strip = {'res_name': ['HOH']})

    ggenerator = bagpype.construction.Graph_constructor()
    ggenerator.construct_graph(myprot)

    return(myprot)