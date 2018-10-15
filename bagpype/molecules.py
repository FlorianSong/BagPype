import numpy as np
import warnings

##################################
#                                #
# The following classes are the  #
# key data members, representing #
# atoms, bonds and proteins      #
#                                #
##################################

class Protein(object):
    """ Stores information about atoms/residues/bonds and 
    an associated network of atoms and bonds.
    """

    def __init__(self, **kwargs):
        self.name = []
        self.pdb_id = []
        self.graph = []
        self.atoms = AtomList()
        self.bonds = BondList()
        # dict for mapping chains to atom ids
        self.chains = []
        # dict for mapping residues to atom ids
        self.residues = []
        self.pdbNum_to_atomID = []
    
    def get_atoms(self, selection):
        """ Returns atoms in specified residues.

        :param selection: List of residues with format (res_num, chain)
        :type selection: list of tuples

        :returns: AtomList of atoms in the specified residues  
        """

        atoms = AtomList()
        for residue in selection:
            atoms += AtomList([self.atoms[i] 
                               for i 
                               in self.residues[residue]])
        if atoms:
            return atoms
        else:
            warnings.warn('Residue not found.  AtomList is empty', 
                          RuntimeWarning)
            return atoms

class GNMProtein(object):

    def __init__(self, **kwargs):
        self.name = []
        self.pdb_id = []
        self.graph = []
        self.residues = []
        self.residue_id_map = {}

class Residue(object):

    def __init__(self, id_, res_name, chain, res_num, coordinates, bfactor):
        self.id_ = id_
        self.res_name = res_name
        self.chain = chain
        self.res_num = res_num
        self.bfactor = bfactor
        self.xyz = coordinates

class ResidueList(object):

    def __init__(self):
        self.residues = []

    def append(self, residue):
        self.residues.append(residue)

    def id(self):
        return [residue.id_ for residue in self.residues]
    
    def __getitem__(self, idx):
        return self.residues[idx]

class AtomList(object):
    """ Create a list of atoms. Provides interface to atoms which allows
        accessing members of the list by their PDB number, retrieving 
        the coordinates of the atoms, and getting the id numbers of
        all atoms in the list.
    
    Parameters
    -----------
    atoms : list 
      Initial list of atoms
    """

    def __init__(self, atoms=None):
        if not (atoms is None):
            self.atoms = list(atoms)
        else:
            self.atoms = []

    def coordinates(self):
        """ Returns an N x 3 matrix containing the coordinates of the 
            atoms in the list"""

        nAtoms = len(self.atoms)
        coordinates = np.zeros((nAtoms,3))
        for i in range(nAtoms):
            coordinates[i,:] = self.atoms[i].coordinates()
        return coordinates

    def append(self, atoms):
        self.atoms.append(atoms)

    def extend(self, atoms):
        self.atoms.extend(atoms)

    def id(self):
        return [atom.id for atom in self.atoms]
    
    def find_by_PDBnum(self, idx):
        """Returns an AtomList containing atoms whose PDB ID number
        matches those in idx
        
        Parameters
        ----------
        idx : int or list of int 
          PDB ID numbers of atoms to be found
        """

        atoms = [atom 
                for atom in self.atoms
                if atom.PDBnum == idx]
        if len(atoms) > 1:
            IOError('More than one atom has the same PDB number')
        elif len(atoms) == 0:
            IOError('Atom with given PDB number not found')
        else:
            return atoms[0]

    def index(self, atom):
        return self.atoms.index(atom)

    def __getitem__(self, idx):
        if isinstance(idx, int):
            return self.atoms[idx] 
        elif isinstance(idx, list):
            return AtomList([self.atoms[i] for i in idx])
        elif isinstance(idx, slice):
            return AtomList(self.atoms[idx])

    def __call__(self, selection):
        # Parsing PyMol-like selection
        return None

    def __len__(self):
        return len(self.atoms)

    def __repr__(self):
        return "<AtomList containing {0} atoms>".format(len(self.atoms))

    def __contains__(self, atom):
        return atom in self.atoms

    def __add__(self, x):
        if isinstance(x, AtomList):
            return AtomList(self.atoms + x.atoms)
        elif isinstance(x, Atom):
            new_atoms = list(self.atoms).append(x)
            return AtomList(new_atoms)

    def __iter__(self):
        return iter(self.atoms)


class BondList(object):
    """ Create a list of bonds

    Parameters
    ----------
    bonds : list of bonds
      Bonds to initialise list with.  If None list will initially be empty.
    """

    def __init__(self, bonds=None):
        if bonds:
            self.bonds = bonds
        else:
            self.bonds = []

    def append(self, bonds):
        self.bonds.append(bonds)

    def extend(self, bonds):
        self.bonds.extend(bonds)

    def id(self):
        return [bond.id for bond in self.bonds]

    def return_by_id(self, idx):
        """ Return a list of bonds whose id matches those
        in idx.
        
        Parameters
        ----------
        idx : list of ints 
          List of ids for the bonds you would like returned
        """
        if isinstance(idx, int):
            bond_out = [bond for bond in self.bonds if bond.id == idx]
            if len(bond_out) > 0:
                IOError('More than one bond in the list has the same ID')
            if len(bond_out) == 0:
                IOError('No bond found with that id')
            else: 
                return bond_out[0]

        elif isinstance(idx, list):
            return BondList([bond for bond in self.bonds if bond.id in idx])
        else:
            IOError('Input should be a single integer or a list of integers')
           
    def __getitem__(self, idx):
        if isinstance(idx, int):
            return self.bonds[idx] 
        elif isinstance(idx, list):
            return BondList([self.bonds[i] for i in idx])
        elif isinstance(idx, slice):
            return BondList(self.bonds[idx])
    
    def __len__(self):
        return len(self.bonds)

    def __repr__(self):
        return "<BondList containing %s bonds>" % len(self.bonds)

    def index(self, bond):
        return self.bonds.index(bond)

    def __contains__(self, bond):
        return bond in self.bonds

    def __add__(self, x):
        if isinstance(x, BondList):
            return BondList(self.bonds + x.bonds)
        elif isinstance(x, Bond):
            new_bonds = list(self.bonds).append(x)
            return BondList(new_bonds)

    def __iter__(self):
        return iter(self.bonds)

       
class Atom(object):
    """Create an object representing an atom.

    Parameters
    ----------
    id : int 
      unique id number for atom
    PDBnum : int
      number of the atom in the PDB file
    element : str
      Chemical element of atom
    name : str
      Atom name in the pdb file
    chain : str
      Chain that the atom belongs to.
    res_num : str
      Number of the residue the atom belongs to (may include insertion code)
    res_name : str
      Name of the residue the atom belongs to
    bfactor : float 
      B-factor of the atom (measures uncertainty/disorder)
    xyz : tuple
      co-ordinates of the atom from the PDB file (x, y, z)
    """

    def __init__(self, id, PDBnum, element, name, chain, res_num,
                 res_name, bfactor, coordinates):
        self.id = id
        self.PDBnum = PDBnum
        self.element = element
        self.name = name
        self.chain = chain
        self.res_num = res_num
        self.res_name = res_name
        self.bfactor = bfactor
        self.xyz = coordinates
        
    def coordinates(self):
        return np.atleast_2d(self.xyz)

    def __eq__(self, atom):
        if ((self.PDBnum == atom.PDBnum) and
            (self.name == atom.name) and
            (self.chain == atom.chain) and
            (self.res_num == atom.res_num) and
            (self.res_name == atom.res_name)):
            return True
        else:
            return False

    def __ne__(self, atom):
        return not self.__eq__(atom)

    def __hash__(self):
        return hash((self.PDBnum, self.name, self.chain, 
                     self.res_num, self.res_name))

    def __repr__(self):
        output = ('<Atom ID: {0}. Atom {1} from residue {2} {3} in chain {4}>'
                  .format(self.id, 
                          self.name, 
                          self.res_name, 
                          self.res_num,
                          self.chain,
                          ))
        return output


class Bond(object):
    """Create an object representing a bond

    Parameters
    ----------
    id : int
      unique id number for the bond
    atom1 : :bagpype:Atom
      first atom forming the bond
    atom2 : :bagpype:Atom 
      second atom forming the bond
    weight : float 
      strength of the bond (in J/A^2)
    bond_type : list of str
      Types of interaction contributing to the bond
    """

    def __init__(self, id, atom1, atom2, weight, types):
        self.id = id
        self.atom1 = atom1
        self.atom2 = atom2
        self.weight = weight
        if isinstance(types, list):
            self.bond_type = types
        else:
            self.bond_type = [types] 
        
    def add_interaction(self, weight, types):
        """ Adds a new contribution to the bond from 
            a certain type of interation (covalent/hydrogen/
            hydrophobic/electrostatic)
        """
        
        self.weight += weight
        for bond_type in types:
            #if bond_type not in self.bond_type:
            self.bond_type.append(bond_type)

    def __eq__(self, bond):
        """ Check if two bonds are the same (i.e. they link 
        the same two atoms.
        """

        if (set((self.atom1.id, self.atom2.id)) == 
            set((bond.atom1.id, bond.atom2.id))):
            return True
        else:
            return False

    def __hash__(self):
        return hash(frozenset((self.atom1.id, self.atom2.id)))

    def __ne__(self, bond):
        return not self.__eq__(bond)

    def __repr__(self):
        output = ('<Bond ID: {0}.  {1} bond between atoms {2} '
                  .format(self.id, self.bond_type, self.atom1.PDBnum) + 
                  'and {0} of strength {1}>'
                  .format(self.atom2.PDBnum, self.weight))
        return output
