from __future__ import print_function
import networkx as nx
import os
import warnings
import sys
import array
import shutil
import subprocess
from collections import defaultdict
import string
import numpy as np
import csv

import bagpype.parameters
import bagpype.settings
import bagpype.molecules
# from bagpype.construction import within_cov_bonding_distance
import bagpype.construction

import dependencies.MakeMultimer

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

    def __init__(self, pdb_filename, download = False):
        # self.pdb = pdb_filename
        # filename after adding hydrogens/stripping unwanted atoms etc.
        self.pdb_filename = pdb_filename

        if download:
            import urllib.request
            url = "https://files.rcsb.org/download/" + pdb_filename
            urllib.request.urlretrieve(url, pdb_filename)

        
    def parse(self, protein, model=None, 
                strip='default', strip_ANISOU=True, remove_LINK=False,
                add_H=None, alternate_location=None, MakeMultimer_number=1):
        """ Takes a bagpype Protein object and loads it with data 
        from the given PDB file.

        Parameters
        ----------
        # bio : bool
        #   Indicates whether the file is a biological pdb file.
        #   If True, multiple models in the same file will be combined.
        model : int
          Model to be loaded.  Should be specified if there
          are multiple models in the same file (e.g. NMR structures)
        strip : tuple of str 
          name of objects to be removed from PDB file
        strip_ANISOU : bool
          specifies whether remove ANISOU entries should be
          removed from PDB file
        remove_LINK : bool 
          specifies whether LINK entries should be removed 
          from the PDB file
        """

        # Sanity check whether file actually exists
        if not os.path.isfile(self.pdb_filename):
            raise IOError("Couldn't open PDB file " + self.pdb_filename)


        # MakeMultimer step
        are_biomt_entries_present = False
        with open(self.pdb_filename) as f:
            for line in f:
                if line.startswith("REMARK 350"):
                    are_biomt_entries_present = True
                    break
        if MakeMultimer_number is not None and are_biomt_entries_present:
            print("Applying MakeMultimer")
            self._MakeMultimer_wrapper(MakeMultimer_number)
        

        #######################################
        # The area below uses self.pdb_lines
        #######################################
        # Read in lines from this file, so we don't need so many file read/write events.
        self.pdb_lines = open(self.pdb_filename, "r").readlines()
        print(self.pdb_filename)

        # Detect whether ANISOU entries are present and strip if wanted.
        are_anisou_entries_present = False
        with open(self.pdb_filename) as f:
            for line in f:
                if line.startswith("ANISOU"):
                    are_anisou_entries_present = True
                    break
        if strip_ANISOU and are_anisou_entries_present:
            print("Removing ANISOU entries")
            self._strip_ANISOU()

        # Symmetric subunits are often stored as separate models in bio files         
        # if self._check_models():
        #     print("Combining models")
        #     self._combine_models()

        if model is not None:
            self._strip_models(model)
        
        if remove_LINK:
            print("Removing LINK entries")
            self._remove_LINK_entries()

        if not (strip=='default' or isinstance(strip, dict)):
            raise TypeError("'strip' should be a dict")
        self._strip_atoms(strip)

        if alternate_location is not None:
            self._strip_alternate_location(alternate_location)

        self.pdb_filename = self.pdb_filename[0:-4] + '_stripped.pdb'
        open(self.pdb_filename, "w").writelines(self.pdb_lines)
        #######################################
        # The area above uses self.pdb_lines
        #######################################

        self._renumber_atoms()

        if add_H or not self._has_hydrogens():
            print("Adding hydrogens to PDB file")
            self._add_hydrogens()
            print("Finished adding hydrogens")

        self._renumber_atoms()

        # Parse individual lines of the pdb file
        print("Loading atoms into protein object")
        (protein.atoms,
         protein.residues,
         protein.chains,
         protein.pdbNum_to_atomID) = self._parse_pdb_lines()
        
        protein.pdb_id = self.pdb_filename

        protein.LINKs = self._parse_LINK()


    def _strip_ANISOU(self):
        """Strips ANISOU records from the PDB file by calling 
        grep and cat bash commands

        """
        out_lines = []
        for line in self.pdb_lines:
            if not line.startswith("ANISOU"):
                out_lines.append(line)
        self.pdb_lines = out_lines

    def _remove_LINK_entries(self):
        """Removes LINK entries.
        """
        out_lines = []
        for line in self.pdb_lines:
            if not line.startswith("LINK"):
                out_lines.append(line)
        self.pdb_lines = out_lines

    def _strip_atoms(self, strip):
        """Creates a new PDB file with atoms specified in the dictionary
        'strip' removed.
        # "CONECT" entries are also removed.
        
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
        print("Stripping unwanted atom types from the PDB file", ( strip['res_name'] ))

        out_lines = []
        for line in self.pdb_lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom = _parse_atom_line(line)
                if (
                        atom.PDBnum not in strip.get('PDBnum', []) and
                        atom.name not in strip.get('name', []) and
                        atom.element not in strip.get('element', []) and
                        atom.chain not in strip.get('chain', []) and
                        atom.res_num not in strip.get('res_num', []) and
                        atom.res_name not in strip.get('res_name', []) and
                        [atom.res_num, atom.chain] not in strip.get('residues', [])
                ):
                    out_lines.append(line)
            # elif line.startswith('CONECT'):
            #     continue
            else:
                out_lines.append(line)

        self.pdb_lines = out_lines


    def _strip_alternate_location(self, alternate_location):
        out_lines = []
        for line in self.pdb_lines:
            if line.startswith("ATOM") or line.startswith("HETATM"): 
                if line[16:17] in [" "] + [alternate_location]:
                    out_lines += line
            else:
                out_lines += line
        self.pdb_lines = out_lines

    def _strip_models(self, model):
        out_lines = []
        currentModel = False
        for line in self.pdb_lines:
            if line.startswith("ATOM") or line.startswith("HETATM") or line.startswith("TER "):
                if currentModel:
                    out_lines += line
                else:
                    continue
            elif (line.startswith('MODEL') and int(line[10:14]) == model):
                currentModel = True
            elif line.startswith("MODEL"):
                continue
            elif line.startswith('ENDMDL'):
                currentModel = False
            else:
                out_lines += line

        self.pdb_lines = out_lines


    def _MakeMultimer_wrapper(self, MakeMultimer_number):
        header = []
        with open(self.pdb_filename) as temp_file:
            for line in temp_file:
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    header.append(line)
                else:
                    break

        MakeMultimer_options = dict(
        backbone = False,
        nowater = False,
        nohetatm = False,
        renamechains = 1,
        renumberresidues = 0)
        
        pdblurb = open(self.pdb_filename).read()
        r = dependencies.MakeMultimer.PdbReplicator(pdblurb, MakeMultimer_options)
        outfile_template = self.pdb_filename.split('.')[0] + '_mm%s.pdb'

        for i, bm in enumerate(r.biomolecules):
            outfile = outfile_template % (i+1)
            open(outfile, 'w').write("".join(header) + bm.output(self.pdb_filename))

        self.pdb_filename = outfile_template % (MakeMultimer_number)



    def _renumber_atoms(self):
        # This function renumbers the ATOM and HETATM records following the output from reduce such that they can be treated by FIRST
        with open(self.pdb_filename, 'r') as fin:
            number=0

            final_lines = []

            for line in fin:
                temp = list(line)
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    number = number+1
                    temp[6+(5-len(str(number))):11] =  list(str(number))
                    dummy='     '
                    temp[6:6+(5-len(str(number)))] = list(dummy[0:(5-len(str(number)))])
                final_lines+= ["".join(temp)]
        with open(self.pdb_filename, "w") as fout:
            fout.writelines(final_lines)      


    def _has_hydrogens(self):
        """Search PDB file for hydrogen atoms. If one is found the 
        file is assumed to have had hydrogens previously added.
        """

        with open(self.pdb_filename,'r') as pdbf:
            for line in pdbf:
                if line[0:4] == 'ATOM' and line[13] == 'H':
                    return True
            return False


    def _add_hydrogens(self):
        """ Runs the command-line program Reduce to add hydrogens
        to the specified pdb file
        Runs with subprocess.call to retain compatibility with Python2.7
        """
        if sys.platform.startswith("linux"):
            subprocess.call(bagpype.settings.REDUCE + ' -Quiet -BUILD -DB ' +
                        bagpype.settings.DEPENDENCIES_ROOT + '/reduce_wwPDB_het_dict.txt ' 
                        + self.pdb_filename + ' > ' + self.pdb_filename[0:-4] 
                        + '_H.pdb', shell=True)
        elif sys.platform.startswith("darwin"):
            subprocess.call(bagpype.settings.DEPENDENCIES_ROOT + "/reduce.macosx" + ' -Quiet -BUILD -DB ' +
                        bagpype.settings.DEPENDENCIES_ROOT + '/reduce_wwPDB_het_dict.txt ' 
                        + self.pdb_filename + ' > ' + self.pdb_filename[0:-4] 
                        + '_H.pdb', shell=True)
        else:
            print("Sorry, but adding hydrogens with Reduce is currently only implemented for UNIX based operating systems.")

        self.pdb_filename = self.pdb_filename[0:-4] + '_H.pdb'
        

    def _parse_pdb_lines(self):
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
                        
        atoms = _load_atoms(self.pdb_filename)
        
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

    def _check_models(self):
        """ Checks if there are multiple models in the PDB file
        """

        for line in self.pdb_lines:
            if line.startswith('MODEL'):
                return True
        return False


    def _combine_models(self):
        """ Combine multiple models into a single model. Chains and atoms
        are renamed so that they do not clash.
        """

        lines_new = []
        for i, line in enumerate(self.pdb_lines):
            lines_new.append(list(line.rstrip('\r\n')))

        # Find MODEL and ENDMDL lines in file
        MODEL_lines = [(i, line) for i, line
                    in enumerate(self.pdb_lines)
                    if line.startswith('MODEL')]
        ENDMDL_lines = [(i, line) for i, line
                        in enumerate(self.pdb_lines)
                        if line.startswith('ENDMDL')]

        # Check there are the same number of MODEL and ENDMDL lines
        if not (len(MODEL_lines) == len(ENDMDL_lines)):
            raise Exception('The number of MODEL and ENDMDL lines do not match')

        # Get models
        models = [int(x[1][12:14]) for i, x in enumerate(MODEL_lines)]
        if not (models == list(range(1, len((models))+1))):
            raise Exception('The model numbers are not consecutive')
        nModels = len(models)

        # Get start and end lines for each model
        lineNumbers = dict((i+1, (MODEL_lines[i][0], ENDMDL_lines[i][0]))
                        for i in range(0, nModels))

        # Check that all entries between the start and end
        # of the MODEL are ATOM/HETATM/TER
        for model in models:
            for line in self.pdb_lines[lineNumbers[model][0]+1:lineNumbers[model][1]-1]:
                if not (line.startswith('ATOM')
                        or line.startswith('HETATM')
                        or line.startswith('TER')):
                    raise Exception('There are non-atom lines inside MODEL ' +
                                    ' environment at line: ' + str(i))

        # Get the chain IDs for each model
        # (these will probably be duplicate for each model)
        chainIDs_old = dict((model, []) for model in models)
        for model in models:
            for line in self.pdb_lines[lineNumbers[model][0]:lineNumbers[model][1]]:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    chainIDs_old[model].append(line[21])
            chainIDs_old[model] = sorted(list(set(chainIDs_old[model])))

        new_chain_creator = string.ascii_uppercase + string.ascii_lowercase

        # Create new unique chain IDs
        chainIDs_new = dict((model, []) for model in models)
        chainCount = 0
        for model in models:
            for chain in chainIDs_old[model]:
                chainIDs_new[model].append(new_chain_creator[chainCount])
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

        with open(self.pdb_filename, 'w') as f:
            for line in lines_new:
                print(line, file = f)

        return True



    def _parse_LINK(self):
        """Find the PDB numbers of the atoms in the LINK entries of
        a pdb file and the bond length of the specified interactions
        """

        with open(self.pdb_filename, 'r') as pdb:
            lines = pdb.readlines()
        LINK_bonds = []
        for line in lines:
            if line.startswith('LINK'):
                atom1 = {'name':line[12:16].strip(), 
                         'res_name': line[17:20].strip(), 
                         'res_num': line[22:26].strip(),
                         'chain': line[21]}
                atom2 = {'name': line[42:46].strip(), 
                         'res_name': line[47:50].strip(), 
                         'res_num': line[52:56].strip(),
                         'chain': line[51]}
                distance_between = float(line[74:78])
                LINK_bonds.append((atom1, atom2, distance_between))
        return LINK_bonds
        






def _load_atoms(pdb, PDBnum='all', name='all', res_name='all', chain = 'all',
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
                atom = _parse_atom_line(line)
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
        
def _parse_atom_line(line):
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
    atom = bagpype.molecules.Atom(-1, PDBnum, element, name,
                                        chainID, res_num, res_name,
                                        bfactor, coordinates)
    return atom


