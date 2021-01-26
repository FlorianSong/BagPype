import networkx as nx
import os
import sys
import subprocess
from collections import defaultdict
import string
import numpy as np

import bagpype.molecules

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

    def __init__(self, pdb_filename, download=False):
        # self.pdb = pdb_filename
        # filename after adding hydrogens/stripping unwanted atoms etc.
        self.pdb_filename = pdb_filename

        if download:
            import urllib.request
            
            path_split = pdb_filename.split("/")
            folders = path_split[:-1]
            path = "/".join(folders)
            if len(path) > 0 and not os.path.exists(path):
                os.makedirs(path)
            name = path_split[-1]
            url = "https://files.rcsb.org/download/" + name
            urllib.request.urlretrieve(url, pdb_filename)

    def parse(
        self,
        protein,
        model=1,
        strip="default",
        strip_ANISOU=True,
        strip_LINK=False,
        strip_CONECT=True,
        strip_HETATM=False,
        add_H=None,
        trim_H="",
        alternate_location="A",
        MakeMultimer_number=1,
        strip_weird_H=[],
    ):
        """ Takes a bagpype Protein object and loads it with data 
        from the given PDB file.

        Parameters
        ----------
        model : int
          Model to be loaded.  Should be specified if there
          are multiple models in the same file (e.g. NMR structures)
        strip : dict where keys are from: [PDBnum, name, element, chain, res_num, res_name, residues]
          name of objects to be removed from PDB file
        strip_ANISOU : bool
          specifies whether remove ANISOU entries should be removed from PDB file
        strip_LINK : bool 
          specifies whether LINK entries should be removed from the PDB file
        strip_HETATM : bool
          specifies whether HETATM entries should be removed from the PDB file
        trim_H : str "all" or "notHOH" or ""
          trims all hydrogens from file, except for waters in the case of "notHOH"
        add_H : bool
          specifies whether REDUCE should be run to add H's to the PDB file
        alternate_location : str (1)
          specifies alternate location letter to be considered
        MakeMultimer_number : int or None
          specifies whether MakeMultimer should be run to obtain biological assemblies.
        strip_weird_H : list of int
          specifies which H's should be removed from the structure even if Reduce placed them there
        """

        # Sanity check whether file actually exists
        if not os.path.isfile(self.pdb_filename):
            raise IOError("Couldn't open PDB file " + self.pdb_filename)

        ### Here, pdb_filename = 1xyz.pdb

        #######################################
        # The area below uses self.pdb_lines
        #######################################

        # Read in lines from this file, so we don't need so many file read/write events.
        self.pdb_lines = open(self.pdb_filename, "r").readlines()

        # Detect whether ANISOU entries are present and strip if wanted.
        if strip_ANISOU and any([l.startswith("ANISOU") for l in self.pdb_lines]):
            print("Removing ANISOU entries")
            self._strip_ANISOU()

        # Strip LINK entries if wanted
        if strip_LINK and any([l.startswith("LINK") for l in self.pdb_lines]):
            print("Removing LINK entries")
            self._strip_LINK_entries()

        # Strip all HETATM entries if wanted
        if strip_HETATM and any([l.startswith("HETATM") for l in self.pdb_lines]):
            print("Removing HETATM entries")
            self._strip_HETATM_entries()

        # Strip all HETATM entries if wanted
        if strip_CONECT and any([l.startswith("CONECT") for l in self.pdb_lines]):
            print("Removing CONECT entries")
            self._strip_CONECT_entries()

        # Symmetric subunits are often stored as separate models in bio files
        # if self._check_models():
        #     print("Combining models")
        #     self._combine_models()

        # Strip all models BUT model (int)
        if model is not None and any([l.startswith("MODEL") for l in self.pdb_lines]):
            print("Selecting model number", model, "and deleting all others")
            self._strip_models(model)

        # Strip certain entries, according to residue name, atom number etc.
        if not (strip == "default" or isinstance(strip, dict)):
            raise TypeError("'strip' should be a dict")
        self._strip_atoms(strip)

        if alternate_location is not None:
            self._strip_alternate_location(alternate_location)

        if trim_H:
            if trim_H == True:
                trim_H = "all"
            print(
                {"all": "Deleting ALL H's",
                 "notHOH": "Deleting H's except in waters"
                }[trim_H]
            )
            self._strip_Hs(trim_H)

        self.pdb_filename = self.pdb_filename[0:-4] + "_stripped.pdb"
        open(self.pdb_filename, "w").writelines(self.pdb_lines)

        #######################################
        # The area above uses self.pdb_lines
        #######################################

        ### Here, pdb_filename = 1xyz_stripped.pdb

        # MakeMultimer step
        if MakeMultimer_number is not None and any(
            [l.startswith("REMARK 350") for l in self.pdb_lines]
        ):
            print("Applying MakeMultimer and using Biomolecule number:", MakeMultimer_number)
            self._MakeMultimer_wrapper(MakeMultimer_number)

        ### Here, pdb_filename = 1xyz_stripped_mm#.pdb

        self._renumber_atoms()

        if add_H or not self._has_hydrogens() or trim_H == "notHOH":
            print("Adding hydrogens to PDB file")
            self._add_hydrogens()
            print("Finished adding hydrogens")

        ### Here, pdb_filename = 1xyz_stripped_mm#_H.pdb

        self._renumber_atoms()

        # Sometimes, certain H's added by REDUCE are weird, so here's a way to strip them...
        if len(strip_weird_H) > 0:
            self._strip_weird_H(strip_weird_H)
            self._renumber_atoms()

        # Parse individual lines of the pdb file
        print("Loading atoms into protein object from file:", self.pdb_filename)
        (
            protein.atoms,
            protein.residues,
            protein.chains,
            protein.pdbNum_to_atomID,
        ) = self._parse_pdb_lines()

        protein.pdb_id = self.pdb_filename

        protein.LINKs = self._parse_LINK()

    ####################
    # Helper functions #
    ####################

    def _MakeMultimer_wrapper(
        self,
        MakeMultimer_number,
        MakeMultimer_renamechains=1,
        MakeMultimer_renumberresidues=0,
    ):

        MakeMultimer_options = dict(
            backbone=False,
            nowater=False,
            nohetatm=False,
            renamechains=int(MakeMultimer_renamechains),
            renumberresidues=int(MakeMultimer_renumberresidues),
        )

        pdblurb = open(self.pdb_filename, "r").read()
        r = dependencies.MakeMultimer.PdbReplicator(pdblurb, MakeMultimer_options)
        outfile_template = self.pdb_filename.split(".")[0] + "_mm%s.pdb"

        for i, bm in enumerate(r.biomolecules):
            header = []
            footer = []

            encountered_ATOM_lines = False
            with open(self.pdb_filename) as temp_file:
                for line in temp_file:
                    if not encountered_ATOM_lines and not (
                        line.startswith("ATOM") or line.startswith("HETATM")
                    ):
                        header.append(line)
                    elif encountered_ATOM_lines and not (
                        line.startswith("ATOM") or line.startswith("HETATM")
                    ):
                        footer.append(line)
                    elif not encountered_ATOM_lines and (
                        line.startswith("ATOM") or line.startswith("HETATM")
                    ):
                        encountered_ATOM_lines = True

            outfile = outfile_template % (i + 1)
            MakeMultimer_output = bm.output(self.pdb_filename)

            """ Fix header elements: LINK, HELIX, SHEET """
            new_chain_names = {}
            for chain in bm.collected_chains:
                new_chain_names.setdefault(chain[0], []).append(chain[1])

            new_header = []
            for i, line in enumerate(header):
                if line.startswith("LINK  "):
                    try:
                        chainID1_new_names = new_chain_names[line[21]]
                    except KeyError:
                        # print("WARNING: Link entry between {} and {} was ignores as chain {} is not part of the final biologically relevant molecule.".format())
                        continue
                    try:
                        chainID2_new_names = new_chain_names[line[51]]
                    except KeyError:
                        continue

                    for combo in zip(chainID1_new_names, chainID2_new_names):
                        new_header.append(
                            line[:21] + combo[0] + line[22:51] + combo[1] + line[52:]
                        )

                elif line.startswith("HELIX "):
                    try:
                        chainID1_new_names = new_chain_names[line[19]]
                    except KeyError:
                        # print("WARNING: HELIX entry between {} and {} was ignores as chain {} is not part of the final biologically relevant molecule.".format())
                        continue
                    try:
                        chainID2_new_names = new_chain_names[line[31]]
                    except KeyError:
                        continue

                    for combo in zip(chainID1_new_names, chainID2_new_names):
                        new_header.append(
                            line[:19] + combo[0] + line[20:31] + combo[1] + line[32:]
                        )

                elif line.startswith("SHEET "):
                    try:
                        chainID1_new_names = new_chain_names[line[21]]
                    except KeyError:
                        # print("WARNING: HELIX entry between {} and {} was ignores as chain {} is not part of the final biologically relevant molecule.".format())
                        continue
                    try:
                        chainID2_new_names = new_chain_names[line[32]]
                    except KeyError:
                        continue

                    for combo in zip(chainID1_new_names, chainID2_new_names):
                        new_header.append(
                            line[:21]
                            + combo[0]
                            + line[22:32]
                            + combo[1]
                            + line[33:40]
                            + "\n"
                        )

                else:
                    new_header.append(line)

            open(outfile, "w").write(
                "".join(new_header) + MakeMultimer_output + "".join(footer)
            )

        self.pdb_filename = outfile_template % (MakeMultimer_number)
        # self._MakeMultimer_out = r

    def _strip_ANISOU(self):
        """Strips ANISOU records from the PDB file
        """
        out_lines = []
        for line in self.pdb_lines:
            if not line.startswith("ANISOU"):
                out_lines.append(line)
        self.pdb_lines = out_lines

    def _strip_HETATM_entries(self):
        """Strips ANISOU records from the PDB file
        """
        out_lines = []
        for line in self.pdb_lines:
            if not line.startswith("HETATM"):
                out_lines.append(line)
        self.pdb_lines = out_lines

    def _strip_LINK_entries(self):
        """Removes LINK entries.
        """
        out_lines = []
        for line in self.pdb_lines:
            if not line.startswith("LINK"):
                out_lines.append(line)
        self.pdb_lines = out_lines

    def _strip_CONECT_entries(self):
        """Removes CONECT entries.
        """
        out_lines = []
        for line in self.pdb_lines:
            if not line.startswith("CONECT"):
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
        aa_to_eliminate = ["UNK", "UNL", "UNX"]

        if strip == "default":
            strip = {}
            strip["res_name"] = aa_to_eliminate
        else:
            strip["res_name"] = strip.get("res_name", []) + aa_to_eliminate
        print("Stripping unwanted atom types from the PDB file", (strip["res_name"]))

        out_lines = []
        for line in self.pdb_lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom = _parse_atom_line(line)
                if (
                    atom.PDBnum not in strip.get("PDBnum", [])
                    and atom.name not in strip.get("name", [])
                    and atom.element not in strip.get("element", [])
                    and atom.chain not in strip.get("chain", [])
                    and atom.res_num not in strip.get("res_num", [])
                    and atom.res_name not in strip.get("res_name", [])
                    and (atom.res_num, atom.chain) not in strip.get("residues", [])
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
                    out_lines.append(line)
            else:
                out_lines.append(line)
        self.pdb_lines = out_lines

    def _strip_models(self, model):
        out_lines = []
        currentModel = False
        for line in self.pdb_lines:
            if (
                line.startswith("ATOM")
                or line.startswith("HETATM")
                or line.startswith("TER ")
            ):
                if currentModel:
                    out_lines.append(line)
                else:
                    continue
            elif line.startswith("MODEL") and int(line[10:14]) == model:
                currentModel = True
            elif line.startswith("MODEL"):
                continue
            elif line.startswith("ENDMDL"):
                currentModel = False
            else:
                out_lines.append(line)

        self.pdb_lines = out_lines

    def _strip_Hs(self, option):
        out_lines = []
        for line in self.pdb_lines:
            if (option == "all" and line[76:78].strip() == "H") or (
                option == "notHOH"
                and line[76:78].strip() == "H"
                and line[17:20].strip() != "HOH"
            ):
                continue
            else:
                out_lines.append(line)
        self.pdb_lines = out_lines

    def _strip_weird_H(self, list_of_atom_ids):
        with open(self.pdb_filename, "r") as fin:
            out_lines = []
            for line in fin:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    PDBnum = int(line[6:11])
                    atom_id = PDBnum - 1
                    if atom_id in list_of_atom_ids:
                        continue
                    else:
                        out_lines.append(line)
        with open(self.pdb_filename, "w") as fout:
            fout.writelines(out_lines)

    def _renumber_atoms(self):
        # This function renumbers the ATOM and HETATM records following the output from reduce such that they can be treated by FIRST
        with open(self.pdb_filename, "r") as fin:
            number = 0

            out_lines = []

            for line in fin:
                temp = list(line)
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    number = number + 1
                    if number > 99999:
                        number = 1

                    temp[6 + (5 - len(str(number))) : 11] = list(str(number))
                    dummy = "     "
                    temp[6 : 6 + (5 - len(str(number)))] = list(
                        dummy[0 : (5 - len(str(number)))]
                    )
                out_lines.append("".join(temp))
        with open(self.pdb_filename, "w") as fout:
            fout.writelines(out_lines)

    def _has_hydrogens(self):
        """Search PDB file for hydrogen atoms. If one is found the 
        file is assumed to have had hydrogens previously added.
        """
        with open(self.pdb_filename, "r") as pdbf:
            for line in pdbf:
                if line[0:4] == "ATOM" and line[13] == "H":
                    return True
            return False

    def _add_hydrogens(self):
        """ Runs the command-line program Reduce to add hydrogens
        to the specified pdb file
        Runs with subprocess.call to retain compatibility with Python2.7
        """
        current_directory = os.path.dirname(os.path.dirname(bagpype.__file__))

        if sys.platform.startswith("linux"):

            subprocess.call(
                current_directory
                + "/dependencies/reduce"
                + " -Quiet -BUILD -DB "
                + current_directory
                + "/dependencies/reduce_wwPDB_het_dict.txt "
                + "'"
                + self.pdb_filename
                + "'"
                + " > "
                + "'"
                + self.pdb_filename[0:-4]
                + "_H.pdb"
                + "'",
                shell=True,
            )
            self.pdb_filename = self.pdb_filename[0:-4] + "_H.pdb"

        elif sys.platform.startswith("darwin"):

            subprocess.call(
                current_directory
                + "/dependencies/reduce.macosx"
                + " -Quiet -BUILD -DB "
                + current_directory
                + "/dependencies/reduce_wwPDB_het_dict.txt "
                + "'"
                + self.pdb_filename
                + "'"
                + " > "
                + "'"
                + self.pdb_filename[0:-4]
                + "_H.pdb"
                + "'",
                shell=True,
            )
            self.pdb_filename = self.pdb_filename[0:-4] + "_H.pdb"

        else:
            print(
                "Sorry, but adding hydrogens with Reduce is currently only implemented for UNIX based operating systems."
            )

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

    def _combine_models(self):
        """ Combine multiple models into a single model. Chains and atoms
        are renamed so that they do not clash.
        """

        lines_new = []
        for i, line in enumerate(self.pdb_lines):
            lines_new.append(list(line.rstrip("\r\n")))

        # Find MODEL and ENDMDL lines in file
        MODEL_lines = [
            (i, line) for i, line in enumerate(self.pdb_lines) if line.startswith("MODEL")
        ]
        ENDMDL_lines = [
            (i, line) for i, line in enumerate(self.pdb_lines) if line.startswith("ENDMDL")
        ]

        # Check there are the same number of MODEL and ENDMDL lines
        if not (len(MODEL_lines) == len(ENDMDL_lines)):
            raise Exception("The number of MODEL and ENDMDL lines do not match")

        # Get models
        models = [int(x[1][12:14]) for i, x in enumerate(MODEL_lines)]
        if not (models == list(range(1, len((models)) + 1))):
            raise Exception("The model numbers are not consecutive")
        nModels = len(models)

        # Get start and end lines for each model
        lineNumbers = dict(
            (i + 1, (MODEL_lines[i][0], ENDMDL_lines[i][0])) for i in range(0, nModels)
        )

        # Check that all entries between the start and end
        # of the MODEL are ATOM/HETATM/TER
        for model in models:
            for line in self.pdb_lines[
                lineNumbers[model][0] + 1 : lineNumbers[model][1] - 1
            ]:
                if not (
                    line.startswith("ATOM")
                    or line.startswith("HETATM")
                    or line.startswith("TER")
                ):
                    raise Exception(
                        "There are non-atom lines inside MODEL "
                        + " environment at line: "
                        + str(i)
                    )

        # Get the chain IDs for each model
        # (these will probably be duplicate for each model)
        chainIDs_old = dict((model, []) for model in models)
        for model in models:
            for line in self.pdb_lines[lineNumbers[model][0] : lineNumbers[model][1]]:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    chainIDs_old[model].append(line[21])
            chainIDs_old[model] = sorted(list(set(chainIDs_old[model])))

        new_chain_creator = string.ascii_uppercase + string.ascii_lowercase

        # Create new unique chain IDs
        chainIDs_new = dict((model, []) for model in models)
        chainCount = 0
        for model in models:
            for chain in chainIDs_old[model]:
                chainIDs_new[model].append(new_chain_creator[chainCount])
                chainCount = chainCount + 1

        # Map old chain id to new chain id
        chainIDmapping = dict((el, {}) for el in models)
        for model in models:
            for i, chain in enumerate(chainIDs_old[model]):
                chainIDmapping[model][chain] = chainIDs_new[model][i]

        for model in models:
            for i, line in enumerate(
                lines_new[lineNumbers[model][0] : lineNumbers[model][1]]
            ):
                if "".join(line[0:4]) == "ATOM" or "".join(line[0:6]) == "HETATM":
                    line[21] = chainIDmapping[model][line[21]]
                    lines_new[lineNumbers[model][0] + i] = line
        for i, line in enumerate(lines_new):
            lines_new[i] = "".join(line)

        i = 0
        for model in models:
            lines_new.pop(lineNumbers[model][0] - i)
            i += 1
            lines_new.pop(lineNumbers[model][1] - i)
            i += 1

        with open(self.pdb_filename, "w") as f:
            for line in lines_new:
                print(line, file=f)

        return True

    def _parse_LINK(self):
        """Find the PDB numbers of the atoms in the LINK entries of
        a pdb file and the bond length of the specified interactions
        """

        with open(self.pdb_filename, "r") as pdb:
            lines = pdb.readlines()
        LINK_bonds = []
        for line in lines:
            if line.startswith("LINK  "):
                atom1 = {
                    "name": line[12:16].strip(),
                    "res_name": line[17:20].strip(),
                    "res_num": line[22:27].strip(),
                    "chain": line[21],
                }
                atom2 = {
                    "name": line[42:46].strip(),
                    "res_name": line[47:50].strip(),
                    "res_num": line[52:57].strip(),
                    "chain": line[51],
                }
                distance_between = float(line[74:78])
                LINK_bonds.append((atom1, atom2, distance_between))
        return LINK_bonds


def _load_atoms(
    pdb, PDBnum="all", name="all", res_name="all", chain="all", res_num="all", model=None
):
    """ Load a list of atoms from a pdb file
    Note that the selection process is currently obsolete. This is done by _strip_atoms().
    
    Returns
    --------
    atoms : :bagpype:AtomList
      AtomList of atoms with given PDBnums/name/residue/model
    """

    atoms = bagpype.molecules.AtomList()
    currentModel = False
    id_counter = 0
    with open(pdb, "r") as f:
        for line in f:
            if not currentModel:
                if (line.startswith("MODEL") and int(line[10:14]) == model) or model is None:
                    currentModel = True
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom = _parse_atom_line(line)
                if (
                    currentModel
                    and (PDBnum == "all" or atom.PDBnum in PDBnum)
                    and (name == "all" or atom.name in name)
                    and (res_name == "all" or atom.res_name in res_name)
                    and (chain == "all" or atom.chain in chain)
                    and (res_num == "all" or atom.res_num in res_num)
                ):
                    atom.id = id_counter
                    atoms.append(atom)
                    id_counter += 1
            if line.startswith("ENDMDL"):
                currentModel = False
    return atoms


def _parse_atom_line(line):
    """ Extract data from a single ATOM line in a pdb file

    Returns
    -------
    atom : :bagpype:Atom
      atom object representing the atom in the input line
    """

    element = line[76:78].strip().upper()
    name = line[12:16].strip()
    PDBnum = int(line[6:11])
    chainID = line[21]
    res_num = line[22:26].strip()
    iCode = line[26:27].strip()
    res_num += iCode
    res_name = line[17:20].strip()  # + line[26].strip()
    try:
        bfactor = float(line[60:66])
    except ValueError:
        bfactor = None

    coordinates = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])

    if element == "D":
        element = "H"  # This is for compatibility with heavy waters

    # Create Atom object with this information
    atom = bagpype.molecules.Atom(
        -1, PDBnum, element, name, chainID, res_num, res_name, bfactor, coordinates
    )
    return atom
