import networkx as nx
import os
import subprocess
from collections import defaultdict
import string
import numpy as np
import csv
import bagpype.parameters
import bagpype.settings
import bagpype.molecules
import warnings



#####################################
#                                   #
# The following functions deal with #
#      parsing of the PDB file      #
#                                   #
#####################################

class PDBParser(object):
    """Class for loading information on atoms and residues 
    from specified PDB file. 

    Parameters
    -----------
    pdb_filename : str 
      PDB ID or file name
    """

    def __init__(self, pdb_filename):
        self.pdb = pdb_filename
        # filename after adding hydrogens/stripping unwanted atoms etc.
        self.pdb_final = self.pdb         
        
    def parse(self, protein, bio=False, model=None, chain='all', 
                strip='default', strip_ANISOU=True, remove_LINK=False,
              add_H=None):
        """ Takes a bagpype Protein object and loads it with data 
        from the given PDB file.

        Parameters
        ----------
        bio : bool
          Indicates whether the file is a biological pdb file.
          If True, multiple models in the same file will be combined.
        model : str
          Model to be loaded.  Should be specified if there
          are multiple models in the same file (e.g. NMR structures)
        chain : tuple of str
          Tuple of chains to load.  
        strip : tuple of str 
          name of objects to be removed from PDB file
        strip_ANISOU : bool
          specifies whether remove ANISOU entries should be
          removed from PDB file
        remove_LINK : bool 
          specifies whether LINK entries should be removed 
          from the PDB file
        """
    
        if strip_ANISOU:
            print("Removing ANISOU entries")
            self.strip_ANISOU()

        # Symmetric subunits are often stored as separate models in bio files         
        if self.check_models():
            print("Combining models")
            self.combine_models()

        if remove_LINK:
            print("Removing LINK entries")
            self.remove_LINK_entries()

        print("Stripping unwanted atom types from the PDB file")
        if not (strip=='default' or isinstance(strip, dict)):
            raise TypeError("'strip' should be a dict")
        self.strip_atoms(strip)

        self.renumber_atoms()

        if add_H or not self.has_hydrogens():
            print("Adding hydrogens to PDB file")
            self.add_hydrogens()
        
        # Parse individual lines of the pdb file
        print("Loading atoms into protein object")
        (protein.atoms,
         protein.residues,
         protein.chains,
         protein.pdbNum_to_atomID) = self.parse_pdb_lines(model, chain)
        
        protein.pdb_id = self.pdb_final


    def strip_ANISOU(self):
        """Strips ANISOU records from the PDB file by calling 
        grep and cat bash commands

        """
        subprocess.call("grep -v ANISOU " + self.pdb_final + " > " + 
                        self.pdb_final[0:-4] + "_temp.pdb",
                        shell = True)
        subprocess.call("cat " + self.pdb_final[0:-4] + "_temp.pdb > " + 
                        self.pdb_final, 
                        shell = True) 
        os.remove(self.pdb_final[0:-4] + "_temp.pdb")

    def remove_LINK_entries(self):
        subprocess.call("grep -v LINK " + self.pdb_final + " > " + 
                        self.pdb_final[0:-4] + "_temp.pdb", shell = True)
        subprocess.call("cat " + self.pdb[0:-4] + "_temp.pdb > " + 
                        self.pdb_final, shell = True)
        os.remove(self.pdb_final[0:-4] + "_temp.pdb")

    def strip_atoms(self, strip):
        """Creates a new PDB file with atoms specified in the dictionary
        'strip' removed.
        
        Parameters
        ----------
        strip : dict
          Dictionary 
        """
        aa_to_eliminate = bagpype.settings.aa_to_eliminate

        if strip == 'default':
            strip = {}
            strip['res_name'] = aa_to_eliminate
        else:
            strip['res_name'] = strip.get('res_name', []) + aa_to_eliminate
            
        old_f = open(self.pdb_final, 'r')
        new_f = open(self.pdb_final[0:-4] + '_stripped.pdb', 'w')
        print(( strip['res_name'] ))
        for line in old_f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom = parse_atom_line(line)
                if (
                        atom.PDBnum not in strip.get('PDBnum', []) and
                        atom.name not in strip.get('name', []) and
                        atom.element not in strip.get('element', []) and
                        atom.chain not in strip.get('chain', []) and
                        atom.res_num not in strip.get('res_num', []) and
                        atom.res_name not in strip.get('res_name', []) and
                        [atom.res_num, atom.chain] not in strip.get('residues', [])
                ):
                    new_f.write(line)
            elif line.startswith('CONECT'):
                continue
            else:
                new_f.write(line)

        self.pdb_final = self.pdb_final[0:-4] + '_stripped.pdb'

    
    def renumber_atoms(self):
        subprocess.call('python ' + bagpype.settings.DEPENDENCIES_ROOT +
                        '/atom_renumbering.py ' + self.pdb_final + ' > ' + 
                        self.pdb_final[0:-4] + '_temp.pdb',
                        shell=True)
        subprocess.call('mv ' + self.pdb_final[0:-4] + '_temp.pdb ' + self.pdb_final, 
                        shell=True)

    def has_hydrogens(self):
        """Search PDB file for hydrogen atoms. If one is found the 
        file is assumed to have had hydrogens previously added.
        """

        try:
            pdbf = open(self.pdb_final,'r')
        except:
            raise IOError("Couldn't open PDB file " + self.pdb_final + 
                            "to check for hydrogens")

        for line in pdbf:
            if line[0:4] == 'ATOM' and line[13] == 'H':
                return True
        return False


    def add_hydrogens(self):
        """ Runs the command-line program Reduce to add hydrogens
        to the specified pdb file
        """

        subprocess.call(bagpype.settings.REDUCE + ' -Quiet -BUILD -DB ' +
                        bagpype.settings.DEPENDENCIES_ROOT + '/het_dict.txt ' 
                        + self.pdb_final + ' > ' + self.pdb_final[0:-4] 
                        + '_H_temp.pdb', shell=True)
        
        subprocess.call('python ' + bagpype.settings.DEPENDENCIES_ROOT +
                        '/atom_renumbering.py ' + self.pdb_final[0:-4] + 
                        '_H_temp.pdb' + ' > ' + self.pdb_final[0:-4] + '_H.pdb',
                        shell=True)
        
        os.remove(self.pdb_final[0:-4] + '_H_temp.pdb')
        self.pdb_final = self.pdb_final[0:-4] + '_H.pdb'
        

    def parse_pdb_lines(self, model, chain):
        """Parses the details of the atoms from a pdb file.

        Returns
        -------
        pdb_data : tuple
          tuple of data for the atoms in the protein containing
            1. AtomList of the atoms in the PDB file 
            2. Dict mapping residues to lists of atoms
            3. Dict mapping chains to atom IDs
            4. Dict mapping PDB numbers to unique atom ids
        """

        # Containers for data from PDB files
        residues = defaultdict(list)
        chains = defaultdict(list)
        pdbNum_to_atomID = defaultdict(list)
                        
        #TO DO : use a generator here
        atoms = load_atoms(self.pdb_final, model=model, chain=chain)
        
        for atom in atoms:
            residues[(atom.res_num, atom.chain)].append(atom.id)
            # Add atom id to the chain -> atom map
            chains[atom.chain].append(atom.id)
            # Map between atom number in PDB file and its atom ID in protein
            pdbNum_to_atomID[atom.PDBnum] = atom.id

        # Turn the defaultdict automatic adding of keys off, so
        # strange things do not happen when the user accesses the dicts
        residues.default_factory = None
        chains.default_factory = None
        pdbNum_to_atomID.default_factory = None
        pdb_data = (atoms, residues, chains, pdbNum_to_atomID)
        
        return pdb_data

    def check_models(self):
        """ Checks if there are multiple models in the PDB file
        """
    
        pdbf = open(self.pdb_final, 'r')
        models = False
        for line in pdbf:
            if line.startswith('MODEL'):
                models = True
                break
        return models


    def combine_models(self):
        """ Combine multiple models into a single model. Chains and atoms
        are renamed so that they do not clash.
        """

        pdbf = open(self.pdb_final, 'r')
        lines = pdbf.readlines()
        lines_new = []
        for i, line in enumerate(lines):
            lines_new.append(list(line.rstrip('\r\n')))

        # Find MODEL and ENDMDL lines in file
        MODEL_lines = [(i, line) for i, line
                    in enumerate(lines)
                    if line.startswith('MODEL')]
        ENDMDL_lines = [(i, line) for i, line
                        in enumerate(lines)
                        if line.startswith('ENDMDL')]

        # Check there are the same number of MODEL and ENDMDL lines
        if not (len(MODEL_lines) == len(ENDMDL_lines)):
            raise Exception('The number of MODEL and ENDMDL lines do not match')

        # Get models
        models = [int(x[1][13]) for i, x in enumerate(MODEL_lines)]
        if not (models == list(range(1, len((models))+1))):
            raise Exception('The model numbers are not consecutive')
        nModels = len(models)

        # Get start and end lines for each model
        lineNumbers = dict((i+1, (MODEL_lines[i][0], ENDMDL_lines[i][0]))
                        for i in range(0, nModels))

        # Check that all entries between the start and end
        # of the MODEL are ATOM/HETATM/TER
        for model in models:
            for line in lines[lineNumbers[model][0]+1:lineNumbers[model][1]-1]:
                if not (line.startswith('ATOM')
                        or line.startswith('HETATM')
                        or line.startswith('TER')):
                    raise Exception('There are non-atom lines inside MODEL ' +
                                    ' environment at line: ' + str(i))

        # Get the chain IDs for each model
        # (these will probably be duplicate for each model)
        chainIDs_old = dict((model, []) for model in models)
        for model in models:
            for line in lines[lineNumbers[model][0]:lineNumbers[model][1]]:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    chainIDs_old[model].append(line[21])
            chainIDs_old[model] = sorted(list(set(chainIDs_old[model])))

        # Create new unique chain IDs
        chainIDs_new = dict((model, []) for model in models)
        chainCount = 0
        for model in models:
            for chain in chainIDs_old[model]:
                chainIDs_new[model].append(string.ascii_uppercase[chainCount])
                chainCount = chainCount+1

        # Map old chain id to new chain id
        chainIDmapping = dict((el, {}) for el in models)
        for model in models:
            for i, chain in enumerate(chainIDs_old[model]):
                chainIDmapping[model][chain] = chainIDs_new[model][i]

        for model in models:
            for i, line in enumerate(lines_new[lineNumbers[model][0]:lineNumbers[model][1]]):
                if ''.join(line[0:4]) == 'ATOM' or "".join(line[0:6]) == 'HETATM':
                    line[21] = chainIDmapping[model][line[21]]
                    lines_new[lineNumbers[model][0]+i] = line
        for i, line in enumerate(lines_new):
            lines_new[i] = "".join(line)

        i = 0
        for model in models:
            lines_new.pop(lineNumbers[model][0]-i)
            i += 1
            lines_new.pop(lineNumbers[model][1]-i)
            i += 1

        f = open(self.pdb_final, 'w')
        for line in lines_new:
            print(line, file = f)
        f.close()

        return True


def load_atoms(pdb, PDBnum='all', name='all', res_name='all', chain = 'all',
               res_num='all', model=None):
    """ Load a list of atoms from a pdb file
    
    Returns
    --------
    atoms : :bagpype:AtomList
      AtomList of atoms with given PDBnums/name/residue/model
    """
    
    atoms = bagpype.molecules.AtomList()
    currentModel = False
    id_counter = 0
    with open(pdb, 'r') as f:
        for line in f:
            if not currentModel:
                if ((line.startswith('MODEL') and int(line[10:14]) == model)
                    or model is None):
                    currentModel = True
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom = parse_atom_line(line)
                if (
                    currentModel and 
                    (PDBnum == 'all' or atom.PDBnum in PDBnum) and
                    (name == 'all' or atom.name in name) and
                    (res_name == 'all' or atom.res_name in res_name) and
                    (chain == 'all' or atom.chain in chain) and
                    (res_num == 'all' or atom.res_num in res_num)
                ):
                    atom.id = id_counter
                    atoms.append(atom)
                    id_counter += 1
                    if line.startswith('ENDMDL'):
                        currentModel = False
    return atoms
        
def parse_atom_line(line):
    """ Extract data from a single ATOM line in a pdb file

    Returns
    -------
    atom : :bagpype:Atom
      atom object representing the atom in the input line
    """

    element = line[76:78].strip()
    name = line[12:16].strip()
    PDBnum = int(line[6:11])
    chainID = line[21]
    res_num = line[22:27].strip()
    res_name = line[17:20].strip()# + line[26].strip()
    try:
        bfactor = float(line[60:66])
    except ValueError:
        bfactor = None

    coordinates = np.array([float(line[30:38]),
                            float(line[38:46]),
                            float(line[46:54])])
                
    # Create Atom object with this information
    atom = bagpype.molecules.Atom(0, PDBnum, element, name,
                                        chainID, res_num, res_name,
                                        bfactor, coordinates)
    return atom


