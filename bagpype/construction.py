import networkx as nx
import bagpype.parameters
import bagpype.settings
import bagpype.molecules
import bagpype.parsing
import scipy 
import numpy as np
import time
import itertools
import pandas as pd


#####################################
#                                   #
# The following functions deal with #
# generation of the protein graph   #
#                                   #
#####################################

class Graph_constructor(object):
    """ This class is responsible for construction of the atomistic
    graph of a protein.  It identifies covalent bonds, hydrogen
    bonds, hydrophobic tethers, stacking interactions in DNA and .
    Atom data should already have been loaded into the protein.
    """

    def __init__(self):

        ### Hydrogen bond related parameters:
        # old values commented at the very right
        self.H_acceptor_cutoff = 4 #2.6 
        self.donor_acceptor_cutoff = 5 #3.6 
        self.donor_H_acceptor_min = 100 #90 
        self.donor_H_acceptor_max = 180
        self.acceptor_H_neighbour_cutoff = 80

        self.H_bond_energy_cutoff = -0.01

        self.hydrophobic_RMST_gamma = 0.05
        self.hydrophobic_neighbourhood = 1

        # Used to initialise all possible bonds, should always be LARGER than the longest possible bond
        self.max_cutoff = 9 

        # Currently obsolete.
        self.k_factor = 6.022/4.184


        # List of all bonds
        self.bonds = bagpype.molecules.BondList()

    def construct_graph(self, protein, atoms_file_name = "atoms.csv", bonds_file_name = "bonds.csv", gexf_file_name = None):
        """ This is the main driver function which calls
        the functions to generate each of the different types
        of bonds.
        """
        # The protein should have already had atom data loaded
        # from the pdb file

        print()
        print("Graph construction started")

        self.protein = protein
        
        time1 = time.time()

        self._initialise_possible_bonds()
        if self.protein.LINKs:
            self._process_LINKs()

        self.find_covalent_bonds()

        # This has to be done so that the other functions can use the covalent bonds too:
        self.covalent_bonds_graph = nx.Graph()
        for bond in self.bonds:
            if bond.bond_type == ["COVALENT"]:
                self.covalent_bonds_graph.add_edge(bond.atom1.id, bond.atom2.id)


        # The following functions require knowledge of covalent bonds!
        self.find_hydrogen_bonds(energy_cutoff = self.H_bond_energy_cutoff) 
        self.find_hydrophobics(gamma = self.hydrophobic_RMST_gamma, neighbourhood_extent = self.hydrophobic_neighbourhood)
        self.find_stacked()
        self.find_DNA_backbone()
        self.find_electrostatics()

        # Make list of bonds unique in the sense that there is only
        # one bond per pair of atoms (that can include multiple
        # contributions from hydrophobic, hydrogen, electrostatics etc.)
        self.bonds = uniquify(self.bonds)

        graph = nx.Graph()
        graph.add_nodes_from(self.protein.atoms.id())

        # Construct graph from list of bonds
        id_counter = 0
        for bond in self.bonds:
            bond.id = id_counter
            graph.add_edge(bond.atom1.id, bond.atom2.id,
                            weight = bond.weight, 
                            id = id_counter,
                            bond_type = bond.bond_type)
            id_counter += 1
            

        # Check graph is connected
        if nx.number_connected_components(graph) > 1:
            print('WARNING: Number of connected components is greater than 1.')
        
        protein.bonds = self.bonds
        protein.graph = graph

        if atoms_file_name is not None:
            self.write_atoms_to_csv_file(atoms_file_name)
        if bonds_file_name is not None:
            self.write_bonds_to_csv_file(bonds_file_name)


        if gexf_file_name is not None:
            graph_modded = graph.copy()
            for _,__,d in graph_modded.edges(data=True):
                d["bond_type"] = "".join(d["bond_type"])
            nx.write_gexf(graph_modded, gexf_file_name)

        time2 = time.time()

        print("Finished constructing the graph!")
        print(("    Time taken = %.2fs" % (time2-time1)))
        print()


    def _initialise_possible_bonds(self):
        """ Initialises all possible bonds (ie all atoms within certain distance of atom), which saves a lot of computing time later 
        Relies on the fact that atom ids run from 0 to # of atoms
        """
        print("Initialising all possible bonds...")
        self.possible_bonds = nx.Graph()

        # Get coordinates of all atoms as an array
        atom_coords = self.protein.atoms.coordinates()

        # Compute a distance matrix between each atom
        dist_matrix = scipy.spatial.distance.cdist(atom_coords, atom_coords)

        # Check that there are as many atoms as rows in the matrix
        if dist_matrix.shape[0] != len(self.protein.atoms) or dist_matrix.shape[1] != len(self.protein.atoms):
            raise Exception("The distance matrix does not have the same number of rows and columns as the number of atoms, so the indices will get messed up.")
        
        # Only keep upper triangular entries (setting lower triangle to max cutoff +1 so they get filtered out in the line below)
        # il1 = np.tril_indices(dist_matrix.shape[0])
        # dist_matrix[il1] = self.max_cutoff + 1

        # self.distance_matrix = dist_matrix

        # Find all entries that are below the cutoff and write all those pairs of atoms to a nx graph
        ind1, ind2 = np.where(dist_matrix <= self.max_cutoff)
        all_possible_bonds =  list(zip(ind1, ind2))
        self.possible_bonds.add_edges_from(all_possible_bonds)


    def write_bonds_to_csv_file(self, name):
        bond_data = []
        for bond in self.bonds:
            bond_data.append( [bond.id, ",".join(bond.bond_type), 
                             bond.weight, distance_between_two_atoms(bond.atom1, bond.atom2),
                             bond.atom1.id, bond.atom1.name, bond.atom1.res_name, bond.atom1.res_num, bond.atom1.chain,
                             bond.atom2.id, bond.atom2.name, bond.atom2.res_name, bond.atom2.res_num, bond.atom2.chain,
                             ] )
        bond_df = pd.DataFrame(data = bond_data, columns = [
                               "bond_id", "bond_type", 
                               "bond_weight", "bond_distance",
                               "atom1_id", "atom1_name", "atom1_res_name", "atom1_res_num", "atom1_chain",
                               "atom2_id", "atom2_name", "atom2_res_name", "atom2_res_num", "atom2_chain"
                               ])
        bond_df.to_csv(name, header = True, index = False)

    def write_atoms_to_csv_file(self, name):
        atom_data = []
        for atom in self.protein.atoms:
            atom_data.append([atom.id, atom.PDBnum, atom.name, atom.element, atom.chain, 
                             atom.res_name, atom.res_num, ", ".join([str(i) for i in list(atom.xyz)])
                              ])
        atom_df = pd.DataFrame(data = atom_data, columns = [
                               "id", "PDBnum", "name", "element", "chain",
                               "res_name", "res_num", "xyz"
                               ])
        atom_df.to_csv(name, header = True, index = False)


    def _process_LINKs(self):
        print("Processing LINK entries...")
        LINK_bonds = self.protein.LINKs
        LINK_list = []

        for LINK_bond in LINK_bonds:
            found_atoms = [None, None]
            for i in range(2):
                found = False
                LINK_atom = LINK_bond[i]
                for atom in self.protein.atoms:
                    if ( LINK_atom['name'] == atom.name and
                        LINK_atom['res_name'] == atom.res_name and
                        LINK_atom['res_num'] == atom.res_num and
                        LINK_atom['chain'] == atom.chain ):
                        found = True
                        found_atoms[i] = atom
                        break
                if not found:
                    print(
                        ('WARNING: LINK atom {0} {1} {2} {3} not found'.format(LINK_atom['name'], 
                                                                      LINK_atom['res_name'],
                                                                      LINK_atom['res_num'],
                                                                      LINK_atom['chain'])
                         + '. Ignoring this LINK entry.')
                    )
                    found = False
            if found_atoms[0] and found_atoms[1]:
                is_covalent = within_cov_bonding_distance(found_atoms[0], found_atoms[1])
                new_LINK_entry = (found_atoms, {"is_covalent":is_covalent, "distance":LINK_bond[2]})
                LINK_list.append(new_LINK_entry)
        self._LINK_list = LINK_list



    ##################
    # COVALENT BONDS #
    ##################
    def find_covalent_bonds(self):
        print("Finding covalent bonds...")
        self.cov_bond_energies = self.initialise_covbond_dictionary()

        for atom1 in self.protein.atoms:
            for atom2 in self.protein.atoms[list(self.possible_bonds.neighbors(atom1.id))]:
                if atom2.id > atom1.id:

                    if in_same_residue(atom1, atom2):

                        self.add_intra_residue_bond(atom1, atom2, self.cov_bond_energies)
                    else:
                        self.add_inter_residue_bond(atom1, atom2)
            
        #################
        # Covalent LINK entries        
        if self.protein.LINKs:
            are_there_LINK_covalent_bonds = bool(sum([x[1]["is_covalent"] for x in self._LINK_list]))
            if are_there_LINK_covalent_bonds:
                print("    Considering covalent LINK entries...")

                # To avoid doubly adding covalent LINK entries where they have already been detected, we need the following list of covalent bonds
                cov_bonds = [(bond.atom1.id, bond.atom2.id) for bond in self.bonds if "COVALENT" in bond.bond_type]
                cov_bonds.extend([bond[::-1] for bond in cov_bonds])

                for item in self._LINK_list:
                    atom1, atom2 = item[0]

                    # Is this LINK entry covalent and is this bond not already included?
                    if item[1]["is_covalent"] and (atom1.id, atom2.id) not in cov_bonds:

                        bond_strength = self.add_covalent_bonds_using_distance_contraints(atom1, atom2, print_warnings = True)
                        if bond_strength is not None:
                            bond_strength *= self.k_factor/6.022
                            self.bonds.append(bagpype.molecules.Bond([], atom1, atom2, bond_strength, 'COVALENT'))
        #
        #################

    def initialise_covbond_dictionary(self):
        list_of_present_residues = []
        for atom in self.protein.atoms:
            list_of_present_residues.append(atom.res_name)
        
        list_of_present_residues = sorted(list(set(list_of_present_residues)))
        dictionary = bagpype.energies.generate_energies_dictionary(list_of_present_residues)
        
        for key in bagpype.parameters.non_standard_residues:
            dictionary[key] = {}
        dictionary.update(bagpype.parameters.non_standard_residues)
        
        return dictionary

    def add_intra_residue_bond(self, atom1, atom2, cov_bond_energies_input):
        residue = atom1.res_name
        bond_dict = cov_bond_energies_input[residue]

        try:
            # Checking only one way is enough, because covalent bond energies are generated both ways round by default.
            bond_strength = bond_dict[atom1.name][atom2.name]
        except KeyError:
            bond_strength = self.add_covalent_bonds_using_distance_contraints(atom1, atom2, print_warnings = True)
            
        if bond_strength is not None:
            bond_strength *= self.k_factor/6.022
            self.bonds.append(bagpype.molecules.Bond([], atom1, atom2, bond_strength, 'COVALENT'))

    def add_inter_residue_bond(self, atom1, atom2):
        # peptide bond, disulfide bridge and nucleic bond
        conditions = \
        (atom1.element == "C" and atom2.element == "N") or \
        (atom1.element == "N" and atom2.element == "C") or \
        (atom1.element == "S" and atom2.element == "S") or \
        (atom1.element == 'O' and atom2.element == 'P') or \
        (atom1.element == 'P' and atom2.element == 'O')

        if conditions:
            bond_strength = self.add_covalent_bonds_using_distance_contraints(atom1, atom2)
            if bond_strength is not None:
                bond_strength *= self.k_factor/6.022
                self.bonds.append(bagpype.molecules.Bond([], atom1, atom2, bond_strength, 'COVALENT'))

    def add_covalent_bonds_using_distance_contraints(self, atom1, atom2, print_warnings = False):
        if within_cov_bonding_distance(atom1, atom2):
            if print_warnings:
                print(("WARNING: Adding bond between " + str(atom1.id) + " " + atom1.name + " in res " + atom1.res_name + atom1.res_num + 
                       " and " + str(atom2.id) + " " + atom2.name + " in res " + atom2.res_name + atom2.res_num + " based on distance constraints using single bond energies!"))
            strength = bagpype.parameters.single_bond_energies[atom1.element][atom2.element]
        else:
            strength = None
        return strength






    ##################
    # HYDROGEN BONDS #
    ##################

    def find_hydrogen_bonds(self, energy_cutoff):
        """ Determine hydrogen bonds 
        """
        print("Finding hydrogen bonds at cutoff = " + str(energy_cutoff) + "kcal/mol...")

        print("    Assigning H-bond status...")
        # Assign H-bond status to all N, O, S atoms in the structure
        # Very efficient, so ok if some might not get used
        status = {}
        for atom in self.protein.atoms:
            if atom.element in ["N", "O", "S"]:
                status[atom.id] = self.assign_Hbond_status(atom)
        self.Hbond_status = status

        print("    Applying constraints and computing bond strengths...")
        # This is where hydrogen bonds are identified
        for hydrogen in [atom for atom in self.protein.atoms if atom.element == "H"]:

            donor = self.find_donor(hydrogen)

            if donor.element not in ["N", "O", "S"]:
                continue
            if self.Hbond_status[donor.id] == None:
                continue

            for acceptor in [atom for atom in self.protein.atoms[list(self.possible_bonds.neighbors(hydrogen.id))] if atom.element in ["N", "O", "S"]]:
                
                if self.Hbond_status[acceptor.id] == None:
                    continue

                if self.is_donor(donor) and self.is_acceptor(acceptor) and not in_third_neighbourhood(self.covalent_bonds_graph, hydrogen, acceptor):
                    # Check criteria:  
                    r = distance_between_two_atoms(acceptor, hydrogen)
                    d = distance_between_two_atoms(acceptor, donor) 
                    
                    if d < self.donor_acceptor_cutoff and r < self.H_acceptor_cutoff:
                        theta = self.get_angle(donor, hydrogen, acceptor)

                        if self.donor_H_acceptor_min <= np.rad2deg(theta) <= self.donor_H_acceptor_max:

                            salt_bridge_indicator = self.is_salt_bridge(hydrogen, donor, acceptor)
                            
                            if salt_bridge_indicator:
                                bond_strength = self.compute_salt_bridge_energy(d)
                            else:
                                bond_strength = self.compute_hydrogen_bond_energy(hydrogen, donor, acceptor, d, theta)

                            if bond_strength is not None and (bond_strength < energy_cutoff or (salt_bridge_indicator and bond_strength < 0)):
                                atom1, atom2 = hydrogen, acceptor
                                bond_strength *= -self.k_factor*4.184/6.022
                                # bond_strength = abs(bond_strength)
                                self.bonds.append(bagpype.molecules.Bond([], atom1, atom2, bond_strength, {True: "SALTBRIDGE", False:"HYDROGEN"}[salt_bridge_indicator]))

    def get_angle(self, atom1, atom2, atom3):
        """ Finds the angle formed by the three input atoms, where atom2 is the vertex of the angle.
        input: three atom objects
        output: angle IN RAD!!
        """
        v21 = atom2.xyz - atom1.xyz
        v23 = atom2.xyz - atom3.xyz

        cosine_angle = np.dot(v21, v23)/(np.linalg.norm(v21) * np.linalg.norm(v23))
        
        # Note: as the arccos is only defined on [0, pi], we do not have to worry about this.
        return np.arccos(cosine_angle)

    def get_out_of_plane_angle(self, donor_list, acceptor_list):
        """ Find the angle between the two planes formed by donor_list and acceptor_list 
        input: two lists of three atom ids
        output: angle IN RAD!!
        """

        def normal_unit_vector(point1, point2, point3):
            """ Function to find normal unit vector for 3 points in 3D space
            """
            v12 = point1 - point2
            v13 = point1 - point3
            normal_vector = np.cross(v12, v13)
            return normal_vector/np.linalg.norm(normal_vector)

        # Find normal unit vector for both planes
        donor_plane_normal = normal_unit_vector(self.protein.atoms[donor_list[0]].xyz,
                                                self.protein.atoms[donor_list[1]].xyz,
                                                self.protein.atoms[donor_list[2]].xyz)
        acceptor_plane_normal = normal_unit_vector(self.protein.atoms[acceptor_list[0]].xyz,
                                                   self.protein.atoms[acceptor_list[1]].xyz,
                                                   self.protein.atoms[acceptor_list[2]].xyz)

        # The angle is now arccos(normal1 *dot* normal2)
        gamma = np.arccos(np.dot(donor_plane_normal, acceptor_plane_normal))

        if gamma < np.pi/2.:
            gamma = np.pi - gamma
        return gamma

    def find_donor(self, atom):
        """" This identifies the covalent bonding partner of a hydrogen
        atom. There should be exactly one. This will raise an error if
        the atom is not a hydrogen or the atom has more than 1 covalent bonding
        partner.
        """
        if atom.element != 'H':
            raise ValueError("You have tried to find the donor of atom {0}, "
                             "but atom {0} is not a hydrogen.".format(atom.id))

        cov_bonding_partners_ids = list(self.covalent_bonds_graph.neighbors(atom.id))
        # cov_bonding_partners_ids = [(bond.atom1.id, bond.atom2.id) for bond in covalent_bonds if bond.atom1.id == atom.id or bond.atom2.id == atom.id]
        if len(cov_bonding_partners_ids) > 1:
            raise ValueError("Atom {0} is a hydrogen, but has more than one covalent "
                             "bonding partner.".format(atom.id))
        elif len(cov_bonding_partners_ids) == 0:
            raise ValueError("Atom {0} is a hydrogen, but has no covalent "
                             "bonding partners.".format(atom.id))
        else:
            return self.protein.atoms[cov_bonding_partners_ids[0]]

    def is_donor(self, atom):
        return True if self.Hbond_status[atom.id]["DonorAcceptor"] in ("donor", "both") else False
            
    def is_acceptor(self, atom):
        return True if self.Hbond_status[atom.id]["DonorAcceptor"] in ("acceptor", "both") else False


    def is_salt_bridge(self, hydrogen, donor, acceptor):
        """ Check if triplet is a salt bridge
        """

        # Donor and acceptor need to be charged
        if self.Hbond_status[donor.id]["Charged"] and self.Hbond_status[acceptor.id]["Charged"]:
            
            for n in self.covalent_bonds_graph.neighbors(acceptor.id):
                neighbor_atom = self.protein.atoms[n]
                phi = self.get_angle(hydrogen, acceptor, neighbor_atom)
                if np.rad2deg(phi) < self.acceptor_H_neighbour_cutoff:
                    return False
            return True
        else:
            return False

    def compute_salt_bridge_energy(self, donor_acceptor_distance):
        R_s = 3.2
        x = .375
        V_0 = 10

        Rnorm = R_s/(donor_acceptor_distance + x)
        E = V_0 * (5 * (Rnorm ** 12) - 6 * (Rnorm ** 10))
        return E

    def compute_hydrogen_bond_energy(self, hydrogen, donor, acceptor, donor_acceptor_distance, theta):
        """ Computes hydrogen bond energy according to formula which is printed in FIRST manual, derived from Mayo potential
        """

        D_SP2 = False
        D_SP3 = False
        A_SP2 = False
        A_SP3 = False

        donor_triplet_list = self.find_all_atom_triples(hydrogen)
        acceptor_triplet_list = self.find_all_atom_triples(acceptor)

        R_s = 2.8
        V_0 = 8.
        Rnorm = R_s/donor_acceptor_distance

        E_distance = V_0 * (5 * (Rnorm ** 12) - 6 * (Rnorm ** 10))

        prefactor = (np.cos(theta) ** 2) * np.exp( - (np.pi - theta) ** 6)

        # Set a few booleans denoting Hybridisation, in particular D_SP2, D_SP3, A_SP2, A_SP3
        if self.Hbond_status[donor.id]["Hybridisation"] == "sp2":
            D_SP2 = True
        elif self.Hbond_status[donor.id]["Hybridisation"] == "sp3":
            D_SP3 = True
        else:
            raise ValueError("Hybridisation state for atom " + str(donor.id) + " is invalid: " + str(self.Hbond_status[donor.id]["Hybridisation"]))
        
        if self.Hbond_status[acceptor.id]["Hybridisation"] == "sp2":
            A_SP2 = True
        elif self.Hbond_status[acceptor.id]["Hybridisation"] == "sp3":
            A_SP3 = True
        else:
            raise ValueError("Hybridisation state for atom " + str(acceptor.id) + " is invalid: " + str(self.Hbond_status[acceptor.id]["Hybridisation"]))


        if D_SP2 and A_SP2:

            angle = 0.
            best_angle = 0.


            for donor_triplet in donor_triplet_list:
                for acceptor_triplet in acceptor_triplet_list:

                    phi = self.get_angle(hydrogen, acceptor, self.protein.atoms[acceptor_triplet[1]])
                    gamma = self.get_out_of_plane_angle(donor_triplet, acceptor_triplet)

                    if np.rad2deg(phi) <= 90:
                        return None

                    if phi > gamma:
                        angle = phi
                    else:
                        angle = gamma

                    if angle > best_angle:
                        best_angle = angle

            E_angular = np.cos(best_angle) ** 2 * prefactor

        elif D_SP3 and A_SP2:

            best_phi = 0.
            for acceptor_triplet in acceptor_triplet_list:
                phi = self.get_angle(hydrogen, acceptor, self.protein.atoms[acceptor_triplet[1]])

                if np.rad2deg(phi) <= 90:
                    return None

                if phi > best_phi:
                    best_phi = phi
            E_angular = np.cos(best_phi) ** 2 * prefactor

        elif D_SP2 and A_SP3:

            E_angular = prefactor ** 2

        elif D_SP3 and A_SP3:
            
            diff_angle = np.deg2rad(109.5)
            best_phi = 100.

            for acceptor_triplet in acceptor_triplet_list:

                phi = self.get_angle(hydrogen, acceptor, self.protein.atoms[acceptor_triplet[1]])

                if phi - diff_angle > np.deg2rad(90.):
                    return None

                if phi < best_phi:
                    best_phi = phi

                E_angular = np.cos(best_phi - diff_angle) ** 2 * prefactor


        E = E_distance * E_angular
        return E

    def find_all_atom_triples(self, atom):
        """ Find all triplets 
        """
        triple_list = []

        site1 = atom.id
        degree = nx.degree(self.covalent_bonds_graph, site1)
        
        if degree < 1:
            raise ValueError("Couldn't find triple for single atom with no bonds!")

        elif degree == 1:
            site2 = list(self.covalent_bonds_graph.neighbors(atom.id))[0]
            for site3 in self.covalent_bonds_graph.neighbors(site2):
                if site3 != site1 and (self.protein.atoms[site1].res_name == "HOH" or self.protein.atoms[site2].element != "H"):
                    triple_list.append([site1, site2, site3])
        else:
             for site2 in self.covalent_bonds_graph.neighbors(site1):
                for site3 in self.covalent_bonds_graph.neighbors(site1):
                    if site2 != site3 and (self.protein.atoms[site1].res_name == "HOH" or self.protein.atoms[site2].element != "H"):
                        triple_list.append([site1, site2, site3])

        return triple_list

    def assign_Hbond_status(self, atom):
        """Determine and assign H-bonding status to atoms
        Input: one atom object
        Writes status into a dictionary with the following possible values: 
        Hybridisation: "sp2", "sp3"
        Charged: True, False
        DonorAcceptor: "donor", "acceptor", "both"
        """        

        status_out = {}


        if atom.element == "N":
            if self.is_mainchain(atom):
                # tests for terminal Nitrogens
                if nx.degree(self.covalent_bonds_graph, atom.id) == 4:
                    status_out.update(DonorAcceptor = "donor", Charged = True, Hybridisation = "sp3")
                else:
                    status_out.update(DonorAcceptor = "donor", Charged = False, Hybridisation = "sp2")
            elif self.get_charged_residue(atom) is not None:
                status_out = self.get_charged_residue(atom)
            else:
                status_out = self.determine_status_N(atom)

        elif atom.element == "O":
            if self.is_mainchain(atom):
                status_out.update(DonorAcceptor = "acceptor", Charged = False, Hybridisation = "sp2")
            elif self.get_charged_residue(atom) is not None:
                status_out = self.get_charged_residue(atom)
            else:
                status_out = self.determine_status_O(atom)

        elif atom.element == "S": 
            status_out = self.determine_status_S(atom)

        else:
            raise ValueError("Function for Hbond status assignment did not receive an atom of element N, O or S.")

        # perform a check:
        if (status_out is not None and
            (status_out["DonorAcceptor"] not in ("donor", "acceptor", "both") or 
            status_out["Charged"] not in (True, False) or 
            status_out["Hybridisation"] not in ("sp2", "sp3"))):
            raise ValueError("Atom " + str(atom.id) + " was not correctly given a Hydrogen status. Something must have gone wrong in the code.")

        return status_out


    def is_mainchain(self, atom):
        """ Check if atom is on peptide mainchain
        """

        neighbour_atom_names = [self.protein.atoms[a].name for a in self.covalent_bonds_graph.neighbors(atom.id)]

        # First 3 criteria (ie CA, N and C) refer to the backbone.

        if atom.name == 'CA':
            # Check for N and C as direct neighbours
            if ('N' in neighbour_atom_names and 'C' in neighbour_atom_names):
                return True

        elif atom.name == 'N':
            # Check for CA as direct and C neighbour with dist 2
            if ('CA' in neighbour_atom_names):
                sec_neighbour_ids = sec_neighborhood(self.covalent_bonds_graph, atom.id)
                sec_neighbour_names = [self.protein.atoms[a].name for a in sec_neighbour_ids]
                if 'C' in sec_neighbour_names:
                    return True

        elif atom.name == 'C':
            # check for CA and O as direct neighbours
            if ('CA' in neighbour_atom_names and 'O' in neighbour_atom_names):
                return True

        elif atom.name == 'O':
            # if C is direct neighbour and CA is second nearest neighbour
            if ('C' in neighbour_atom_names):
                sec_neighbour_ids = sec_neighborhood(self.covalent_bonds_graph, atom.id)
                sec_neighbour_names = [self.protein.atoms[a].name for a in sec_neighbour_ids]
                if 'CA' in sec_neighbour_names:
                    return True

        elif atom.name == 'H':
            # if N direct and CA second nearest neighbout
            if ('N' in neighbour_atom_names):
                sec_neighbour_ids = sec_neighborhood(self.covalent_bonds_graph, atom.id)
                sec_neighbour_names = [self.protein.atoms[a].name for a in sec_neighbour_ids]
                if 'CA' in sec_neighbour_names:
                    return True

        # For DNA, RNA...
        # but this case will never actually come up...
        elif atom.name == 'P':
            # if O5 direct neighbour
            if ('O5' in neighbour_atom_names):
                return True

        # Default
        return False


    def get_charged_residue(self, atom):
        """ Check if atom is in a charged residue and assign orbid and charge status
        """

        if atom.element == 'N':
            if atom.res_name == 'ARG':
                return dict(DonorAcceptor = "donor", Charged = True, Hybridisation = "sp2")

            elif atom.res_name == 'LYS':
                return dict(DonorAcceptor = "donor", Charged = True, Hybridisation = "sp3")

            else:
                neighbour_atom_names = [self.protein.atoms[a].name for a in self.covalent_bonds_graph.neighbors(atom.id)]
                if ('CA' in neighbour_atom_names and 'C' not in neighbour_atom_names):
                    return dict(DonorAcceptor = "donor", Charged = True, Hybridisation = "sp3")

        elif atom.element == 'O':
            if atom.res_name == 'ASP' or atom.res_name == 'GLU':
                return dict(DonorAcceptor = "acceptor", Charged = True, Hybridisation = "sp2")

            else:
                neighbour_atom_names = [self.protein.atoms[a].name for a in self.covalent_bonds_graph.neighbors(atom.id)]
                if ('CB' not in neighbour_atom_names):
                    sec_neighbour_ids = sec_neighborhood(self.covalent_bonds_graph, atom.id)
                    sec_neighbour_names = [self.protein.atoms[a].name for a in sec_neighbour_ids]
                    if 'CA' in sec_neighbour_names and 'N' not in sec_neighbour_names:
                        return dict(DonorAcceptor = "acceptor", Charged = True, Hybridisation = "sp2")
        return None

    def determine_status_N(self, atom):
        """ Determines and assigns hydrogen-bonding status to Nitrogen atoms
        """

        if atom.res_name == 'ASN' or atom.res_name == 'GLN':
            return dict(DonorAcceptor = "donor", Charged = False, Hybridisation = "sp2")

        elif atom.res_name == 'HIS':

                hist_state = self.get_histidine_protonation_state(atom)

                if hist_state == 1 or hist_state == 3:
                    return dict(DonorAcceptor = "donor", Charged = True, Hybridisation = "sp2")

                else:
                    return dict(DonorAcceptor = "acceptor", Charged = True, Hybridisation = "sp2")

        elif atom.res_name == 'TRP':
            return dict(DonorAcceptor = "donor", Charged = False, Hybridisation = "sp2")

        else:
            degree = nx.degree(self.covalent_bonds_graph, atom.id)

            if degree == 2:
                for n in self.covalent_bonds_graph.neighbors(atom.id):

                     if self.is_double_bond(atom, self.protein.atoms[n]):
                        return dict(DonorAcceptor = "acceptor", Charged = False, Hybridisation = "sp2")

                return dict(DonorAcceptor = "acceptor", Charged = True, Hybridisation = "sp2")

            elif degree == 3:
                return dict(DonorAcceptor = "donor", Charged = False, Hybridisation = "sp2")

            elif degree == 4:
                return dict(DonorAcceptor = "donor", Charged = True, Hybridisation = "sp3")

            elif degree > 4:
                print("Something went wrong atom {0} is a Nitrogen".format(atom.id) \
                      + "but has more than 4 covalent bonds")
                return None

        #This is a relict in FIRST, this will never happen as the above covers all possible cases. Included for completeness's sake
        return dict(DonorAcceptor = "donor", Charged = False, Hybridisation = "sp2")


    def get_histidine_protonation_state(self, atom):
        """ Determines Histidine protonation state
        """

        state = 0
        for neighbour_atom in [self.protein.atoms[j] for j in self.covalent_bonds_graph.neighbors(atom.id)]:
            if neighbour_atom.element == 'H':
                state += 1
            if neighbour_atom.name == 'CE1':
                CE_atom = neighbour_atom

                for ce_neighbour_atom in [self.protein.atoms[j] for j in self.covalent_bonds_graph.neighbors(CE_atom.id)]:
                    if ce_neighbour_atom.element == 'N' and ce_neighbour_atom.id != atom.id:
                        other_nitrogen = ce_neighbour_atom
                        for other_N_neighbour_atom in [self.protein.atoms[k] for k in self.covalent_bonds_graph.neighbors(other_nitrogen.id)]:
                            if other_N_neighbour_atom.element == 'H':
                                state += 2
        return state

    def is_double_bond(self, atom1, atom2):
        """Check if bond is actually a double bond -- it has to be within 1.4 A
        """

        return True if ((atom1.element == 'N' and atom2.element == 'C') or 
                        (atom2.element == 'C' and atom1.element == 'N')) and \
                        distance_between_two_atoms(atom1, atom2) <= 1.4  \
                        else False

    def determine_status_O(self, atom):
        """ Determines and assigns hydrogen-bonding status to Oxygen atoms
        """

        if atom.res_name == "SER" or atom.res_name == "THR":
            return dict(DonorAcceptor = "both", Charged = False, Hybridisation = "sp3")

        elif atom.res_name == "TYR":
            return dict(DonorAcceptor = "both", Charged = False, Hybridisation = "sp2")

        elif atom.res_name == "ASN" or atom.res_name == "GLN":
            # kept this case here separately in case we go for donor/acceptor assignment later..
            return dict(DonorAcceptor = "acceptor", Charged = False, Hybridisation = "sp2")

        else:
            degree = nx.degree(self.covalent_bonds_graph, atom.id)

            if degree == 0:
                raise Exception("The Oxygen atom " + str(atom.id) + " has no covalent bonding partners.")

            if degree == 1:
                sec_neighbors_id = sec_neighborhood(self.covalent_bonds_graph, atom.id)
                other_O = None

                # check if there is a second O in neighborhood
                for i in sec_neighbors_id:
                    if self.protein.atoms[i].name == 'O':
                        other_O = i
                # if there is than check if it is also single bonded
                if other_O is not None:

                    if nx.degree(self.covalent_bonds_graph, other_O) == 1:
                        return dict(DonorAcceptor = "acceptor", Charged = True, Hybridisation = "sp2")

                    else:
                        return dict(DonorAcceptor = "acceptor", Charged = False, Hybridisation = "sp2")
                # no other O in surroundings..
                else:
                    return dict(DonorAcceptor = "acceptor", Charged = False, Hybridisation = "sp2")

            elif degree == 2:
                any_neighbour_H = False
                for neighbour_atom in [self.protein.atoms[j] for j in self.covalent_bonds_graph.neighbors(atom.id)]:
                    if neighbour_atom.element == 'H':
                        any_neighbour_H = True

                if any_neighbour_H:
                    return dict(DonorAcceptor = "both", Charged = False, Hybridisation = "sp3")

                else:
                    return dict(DonorAcceptor = "acceptor", Charged = False, Hybridisation = "sp3")

            elif degree > 2:
                print("Something went wrong. Atom {0} is an Oxygen ".format(atom.id) \
                      + "but has more than 2 covalent bonds")
                return None

    def determine_status_S(self, atom):
        """ Determines and assigns hydrogen-bonding status to Sulfur atoms
        """

        if atom.res_name == "MET":
            return dict(DonorAcceptor = "acceptor", Charged = False, Hybridisation = "sp3")

        elif atom.res_name == "CYS" or atom.res_name == 'CYX':
            neighbours = self.covalent_bonds_graph.neighbors(atom.id)

            if 'H' in [self.protein.atoms[n].name for n in neighbours]:
                return dict(DonorAcceptor = "both", Charged = False, Hybridisation = "sp3")

            else:
                return dict(DonorAcceptor = "acceptor", Charged = False, Hybridisation = "sp3")

        else:
            degree = nx.degree(self.covalent_bonds_graph, atom.id)

            if degree == 1:
                return dict(DonorAcceptor = "acceptor", Charged = True, Hybridisation = "sp2")

            elif degree == 2:
                neighbours = self.covalent_bonds_graph.neighbors(atom.id)

                if 'H' in [self.protein.atoms[n].name for n in neighbours]:
                    return dict(DonorAcceptor = "both", Charged = False, Hybridisation = "sp3")

                else:
                    return dict(DonorAcceptor = "acceptor", Charged = False, Hybridisation = "sp3")

            elif degree == 3:
                return dict(DonorAcceptor = "donor", Charged = True, Hybridisation = "sp3")

            else:
                print("Something went wrong. Atom {0} is a Sulfur ".format(atom.id) \
                      + "but has more than 3 covalent bonds. Couldn't determine H-bonding status.")
                return None








    ############################
    # HYDROPHOBIC INTERACTIONS #
    ############################

    def find_hydrophobics(self, gamma, neighbourhood_extent):
        """ 
        """
        print(("Finding hydrophobics..."))

        # Initiate list of hydrophobic interactions
        # self._hphobes_list = []
        hphobic_graph = nx.Graph()
        hphobic_graph.add_nodes_from([atom.id for atom in self.protein.atoms])

        # Loop through all atoms
        for atom1 in self.protein.atoms:

            # Requirements for atom1 to be eligible: only Carbon or Sulfur and only bonded to Carbon, Sulfur or Hydrogen
            if atom1.element in ('C', 'S') and self.only_bonded_to_CSH(atom1, neighbourhood_extent):

                for atom2 in self.protein.atoms[ list(self.possible_bonds.neighbors(atom1.id)) ]:

                    # Remaining requirements for hydrophbic interactions to be feasible
                    # Python's logical and/or are short circuit evaluated, so putting all conditions in one is fine, if the most basic condition is in first place
                    conditions = (atom1.id < atom2.id and 
                                  atom2.element in ('C', 'S') and 
                                  self.only_bonded_to_CSH(atom2, neighbourhood_extent) and
                                  not in_same_residue(atom1, atom2) and 
                                  not in_third_neighbourhood(self.covalent_bonds_graph, atom1, atom2)
                                  )
                    if conditions:
                        distance = distance_between_two_atoms(atom1, atom2)
                        energy = self.hydrophobic_potential(distance, atom1.element, atom2.element)
                    else:
                        energy = None

                    if energy is not None:
                        hphobic_graph.add_edge( atom1.id, atom2.id, weight = -energy, distance = distance, energy = energy)

        matches = self.hydrophobic_selection(hphobic_graph, gamma = gamma)

        # Finally, write all bonds to self.bonds
        for bond in matches:
            atom1, atom2 = self.protein.atoms[bond[0]], self.protein.atoms[bond[1]]
            bond_strength = hphobic_graph[bond[0]][bond[1]]["energy"]
            bond_strength *= -self.k_factor*4.184/6.022
            self.bonds.append(bagpype.molecules.Bond([], atom1, atom2, 
                                                          bond_strength, 'HYDROPHOBIC'))


    def hydrophobic_selection(self, graph, gamma, include_MST = True):
        mst = nx.minimum_spanning_tree(graph, weight='energy')
        # print(len([i for i in list(nx.connected_components(mst)) if len(i) != 1]))
        notmst = nx.difference(graph, mst)

        if include_MST:
            matches = sorted( mst.edges )
        else:
            matches = []
        # matches = sorted( mst.edges )

        accepted = 0
        not_accepted = 0

        # find distance to nearest neighbours
        D = nx.adjacency_matrix(graph, weight="energy")
        d = D.min(axis = 0).toarray().flatten().tolist()
        

        for bond_not_in_mst in notmst.edges:
            i,j = bond_not_in_mst[0], bond_not_in_mst[1]

            path = nx.shortest_path(mst, source = i, target = j)
            weights_along_path = [graph[u][v]["energy"] for u,v in zip(path[:-1], path[1:])]

            mlink = max(weights_along_path)


            # d_i = min( [hphobic_graph[m][n]["energy"] for m,n in hphobic_graph.edges(nbunch=i)] )
            # d_j = min( [hphobic_graph[m][n]["energy"] for m,n in hphobic_graph.edges(nbunch=j)] )

            left_hand_side = mlink + gamma*abs(d[i] + d[j])
            right_hand_side = graph[i][j]["energy"]

            if (left_hand_side > right_hand_side ):
                matches += [(i,j)]

                accepted +=1
            else:
                not_accepted +=1
        
        print("    RMST sparsification used. Accepted: " + str(accepted) + " vs. not accepted: " + str(not_accepted) + "; MST size: " + str(len(mst.edges))  )
        matches.sort()

        return matches


    def hydrophobic_potential(self, r, element1, element2):
        c = np.array([3.81679, 5.46692, 7.11677])
        w = np.array([1.68589, 1.39064, 1.57417])
        h = np.array([-.73080, 0.20016, -.09055])
        presum = h * np.exp(- ( (r-c)/w )**2)

        def Lennard_Jones(r, element1, element2):
            epsilon_dict = {"C": -0.11, "S": -0.45}
            epsilon1 = epsilon_dict[element1]
            epsilon2 = epsilon_dict[element2]

            epsilon = np.sqrt(epsilon1*epsilon2)
            r_min = 2

            temp = r_min/r
            return epsilon* ( (temp)**12 - 2*(temp)**6 )  

        return np.sum(presum) + Lennard_Jones(r, element1, element2)

    def only_bonded_to_CSH(self, atom, neighbourhood_extent):
        """ Check whether atom is only bonded to Carbon, Sulfur or Hydrogen
        """
        if neighbourhood_extent == 1:
            for nhbr in self.covalent_bonds_graph.neighbors(atom.id):
                if self.protein.atoms[nhbr].element not in ['C', 'S', 'H']:
                    return False
            return True

        elif neighbourhood_extent == 2:   
            first_and_second_nbhd = set(sec_neighborhood(self.covalent_bonds_graph, atom.id)).union(set(nhbr for nhbr in self.covalent_bonds_graph.neighbors(atom.id)))
            for item in first_and_second_nbhd:
                if self.protein.atoms[item].element not in ['C', 'S', 'H']:
                    return False
            return True

        else: 
            raise Exception("The extent for the neighbourhood variable needs to be either 1 or 2.")

    def VdW_within_cutoff(self, atom1, atom2, cutoff):
        """ Check whether two atoms are within the cutoff (+ radii of both atoms)
        """

        radius1 = {'C':1.7, 'S':1.8}[atom1.element]
        radius2 = {'C':1.7, 'S':1.8}[atom2.element]
        return distance_between_two_atoms(atom1, atom2) <= cutoff + radius1 + radius2








    ########################
    # STACKED INTERACTIONS #
    ########################

    def find_stacked(self):
        """ Function that implements pi stacking interactions, with inspiration taken from Antoine Delmotte's DNA code.
        No input, function will add bonds to internal bond list

        Florian, May 2017
        """
        print("Finding stacked interactions in DNA...")

        # Parameters for DNA stacking
        energy_threshold =  0.0019872041*300*(self.k_factor*4.184/6.022) # Energy of thermal fluctuations at temperature 300 K
        A = 0.214
        C = 4.7e4
        alpha = 12.35
        epsilon = 4.

        def get_atom_specific_parameters(atom_name):
            """ Function to get atom-specific parameters, in particular K (???) & R (vdW radii)
            input: atom name
            output: tuple of relevant K and R
            """
            K = (1., 1., 1., 1.18, 1.36) # Respectively for H, C, C_arom, N and O 
            R = (1.20, 1.70, 1.77, 1.6, 1.5) # Respectively for H, C, C_arom, N and O 

            atom_name = atom_name.strip()
            if atom_name[0] == "H":
                element_identifier = 0
            elif atom_name[0] == "C":
                if atom_name == "C7":
                    element_identifier = 1
                else:
                    element_identifier = 2
            elif atom_name[0] == "N":
                element_identifier = 3
            elif atom_name[0] == "O":
                element_identifier = 4
            else: 
                raise ValueError("Could not find atom: " + atom_name)
            return K[element_identifier], R[element_identifier]

        def get_pisigma_charges(base, normal, coord, base_atom_name):
            """ Function to get sigma and pi charges from data in energies.py 
            input: 
            base = one of DA, DC, DG, DT
            normal = np vector in 3d, the normal to the ring plane
            coord = np vector in 3d, coordinates of the relevant atom
            base_atom_name = name of relevant atom, eg. "C2"
            """

            # retrieve the relevant dictionary from energies.py
            dictio = {
                "DA": bagpype.parameters.DA_sigma_pi,
                "DC": bagpype.parameters.DC_sigma_pi,
                "DT": bagpype.parameters.DT_sigma_pi,
                "DG": bagpype.parameters.DG_sigma_pi,
            }[base]

            # retrieve charges and build list: [pi, sigma, pi]
            atom_charges = dictio[base_atom_name]
            q = [atom_charges[1]]
            q.extend(atom_charges)

            # list of coordinates for point charges
            above = coord + .47*normal
            below = coord - .47*normal
            c = [above, coord, below]

            return q, c


        # lists of indices for all atoms in the nucleobase (on the basis of which the energies are calculated)
        # as well as indices for the rings only
        dna_atom_inds = {}
        dna_atom_ring = {}
        for atom in self.protein.atoms:
            if atom.res_name in ['DT', 'DC', 'DA', 'DG'] and "'" not in atom.name and "P" not in atom.name:
                dna_atom_inds.setdefault((atom.chain, int(atom.res_num), atom.res_name), []).append(atom.id)
                if atom.name in ['C2', 'C4', 'C5', 'C6', 'N1', 'N3']:
                    dna_atom_ring.setdefault((atom.chain, int(atom.res_num), atom.res_name), {}).update({atom.name: atom.id})



        close_DNA_residues = []

        def determine_if_DNA_residues_are_linked(set1, set2):
            for n1 in set1:
                for n2 in self.possible_bonds.neighbors(n1):
                    if n2 in set2:
                        return True
            return False
                       
        set_nuclear_bases = dna_atom_inds.keys()
        for pair_of_bases in itertools.combinations(set_nuclear_bases, 2):
            base1,base2 = pair_of_bases
            set_of_atoms1, set_of_atoms2 = set(dna_atom_inds[base1]), set(dna_atom_inds[base2])
            
            if determine_if_DNA_residues_are_linked(set_of_atoms1, set_of_atoms2):
                close_DNA_residues.append(pair_of_bases)
        
        # for i in close_DNA_residues: print(i)
        # print(len([i for i in itertools.combinations(set_nuclear_bases,2)]))
        # print("Preselected bases", len(close_DNA_residues))




        # Initialisations
        sum_vdw_electro = {}

        # ENERGY CALCULATIONS:
        # for each pair of nucleobases calculate the energies between them using every single pair of atoms
        # formulae for these are in Antoine's thesis
        for pair_i in close_DNA_residues:

            base_key1, base_key2 = pair_i

            # retrieve list of atoms in particular base, get the name of base and compute the normal to the plane
            id_list1 = dna_atom_inds[base_key1]
            base1 = base_key1[2]
            normal1 = np.cross( self.protein.atoms[id_list1[0]].xyz - self.protein.atoms[id_list1[1]].xyz ,  self.protein.atoms[id_list1[0]].xyz - self.protein.atoms[id_list1[2]].xyz )
            normal1 = normal1/np.linalg.norm(normal1)


            id_list2 = dna_atom_inds[base_key2]
            base2 = base_key2[2]
            normal2 = np.cross( self.protein.atoms[id_list2[0]].xyz - self.protein.atoms[id_list2[1]].xyz ,  self.protein.atoms[id_list2[0]].xyz - self.protein.atoms[id_list2[2]].xyz )
            normal2 = normal2/np.linalg.norm(normal2)

            # compute van der Waals component for every pair of atoms between both nucleobases
            vdw = 0
            for id1 in id_list1:
                K1, R1 = get_atom_specific_parameters(self.protein.atoms[id1].name)

                for id2 in id_list2:
                    K2, R2 = get_atom_specific_parameters(self.protein.atoms[id2].name)

                    z = np.asscalar(  np.linalg.norm(self.protein.atoms[id1].xyz - self.protein.atoms[id2].xyz) / (2.*np.sqrt( R1 * R2))  )
                    vdw += (-self.k_factor*4.184/6.022)*K1*K2*(C*np.exp(-alpha*z) - A/z**6.)

            # do the same with the electrostatic component, including two 3x loops for the 3 point charges
            electro = 0
            for id1 in id_list1:
                q1, c1 = get_pisigma_charges(base1, normal1, self.protein.atoms[id1].xyz, self.protein.atoms[id1].name)
                
                for id2 in id_list2:
                    q2, c2 = get_pisigma_charges(base2, normal2, self.protein.atoms[id2].xyz, self.protein.atoms[id2].name)
                    
                    for m in range(3):
                        for n in range(3):
                            electro += (-self.k_factor*4.184/6.022)*(332/epsilon)*q1[m]*q2[n]/(np.linalg.norm( c1[m]-c2[n] ))
            
            sum_vdw_electro[(base_key1, base_key2)] = vdw + electro



        # sum both energies and only keep those that surpass the threshold
        for key in list(sum_vdw_electro.keys()):
            if not sum_vdw_electro[key] > energy_threshold:
                sum_vdw_electro.pop(key)
            
        # for key in sum_vdw_electro:
        #     print(key, sum_vdw_electro[key])
        # print("Final number of accepted interactions", len(sum_vdw_electro))        

        # spread this energy onto the 6 corresponding atoms of the rings
        for key in sorted(sum_vdw_electro, key = lambda x: (x[0][0], x[0][1], x[1][0], x[1][1])):
            base_key1, base_key2 = key[0], key[1]
            
            edge_dic1, edge_dic2 = dna_atom_ring[base_key1], dna_atom_ring[base_key2]

            for ring_atom in ['C2', 'C4', 'C5', 'C6', 'N1', 'N3']:
                bond_strength = sum_vdw_electro[key]/6.

                atom1_id, atom2_id = edge_dic1[ring_atom], edge_dic2[ring_atom]
                atom1, atom2 = self.protein.atoms[atom1_id], self.protein.atoms[atom2_id]

                self.bonds.append(bagpype.molecules.Bond([], atom1, atom2, bond_strength, 'STACKED'))





    #############################
    # DNA BACKBONE INTERACTIONS #
    #############################

    def find_DNA_backbone(self):
        """ Function that implements interactions between phosphates in DNA as per Antoine Delmotte's DNA code.
        No input, function will add bonds to internal bond list.bond

        Florian, May 2017
        """
        print("Finding DNA backbone interactions...")

        # Various parameters
        delta = .24 #0.1522
        epsilon = 4 #80.
        concentration = 0.1 #Mol/l
        debye = 1/(0.329 * (concentration**0.5))

        # get a list of all relevant atoms - all "OP1" atoms
        OP1_atoms = []
        for atom in self.protein.atoms:
            if atom.res_name in ['DT', 'DC', 'DA', 'DG'] and atom.name == "OP1":
                OP1_atoms.append(atom.id)

        # For each pair of OP1 atoms, compute the interaction energy
        # (This is so inefficient but doesn't make too much of a difference in the grand scheme of things 
        # -> think about more efficient implementation if DNA is huge)
        for i in range(len(OP1_atoms)-1):
            for j in range(i, len(OP1_atoms)):
                atom1_id = OP1_atoms[i]
                atom2_id = OP1_atoms[j]

                # Consider only OP1's that are in the same chain and within one residue of each other
                if self.protein.atoms[atom1_id].chain == self.protein.atoms[atom2_id].chain and \
                 (int(self.protein.atoms[atom1_id].res_num) == int(self.protein.atoms[atom2_id].res_num) + 1 or 
                    int(self.protein.atoms[atom1_id].res_num) == int(self.protein.atoms[atom2_id].res_num) - 1):
                    
                    # compute the bond energy using the formula stated in Antoine's thesis
                    rij = np.linalg.norm( self.protein.atoms[atom1_id].xyz - self.protein.atoms[atom2_id].xyz  ) 
                    energy = abs((-self.k_factor*4.184/6.022)*(332/epsilon)*(delta**2)/rij)*np.exp(-1*rij/debye)
                    
                    # add this to the bond dictionary
                    atom1 = self.protein.atoms[atom1_id]
                    atom2 = self.protein.atoms[atom2_id]
                    self.bonds.append(bagpype.molecules.Bond([], atom1, atom2, energy, 'BACKBONE'))








    ##############################
    # ELECTROSTATIC INTERACTIONS #
    ##############################

    def find_electrostatics(self):
        """ Load electrostatic bonds using LINK entries in the PDB file. Some of these may be covalent bonds though. These are added as well.
        """

        if self.protein.LINKs:
            print(("Finding electrostatic interactions using LINK entries...")) 

            # epsilon is the dielectric constant
            epsilon = 4

            for LINK_entry in self._LINK_list:
                
                # Load the two atoms from the LINK entry parsed in parsing.
                atom1, atom2 = LINK_entry[0]
                
                # The bond distance given in the LINK entry is not always accurate.
                bond_length = distance_between_two_atoms(atom1, atom2)
                if LINK_entry[1]["distance"] != round(bond_length,2):
                    print(("WARNING: The LINK entry between {} and {} has a different distance value than the computed one (rounded to two decimal places):" \
                                    " Computed = {}, in PDB = {}. The computed one will be used.").format(atom1, atom2, bond_length, LINK_entry[1]["distance"]))

                if LINK_entry[1]["is_covalent"]:
                    # For this case, see find_covalent().
                    pass 
                else:
                    try:
                        # Find charges of the two atoms using charge database 
                        # in bagpype.parameters file
                        charge1 = bagpype.parameters.charges[atom1.res_name][atom1.name]
                        charge2 = bagpype.parameters.charges[atom2.res_name][atom2.name]

                        # Calculate bond strength using Coulomb's law
                        bond_strength = (332*charge1*charge2)/(epsilon*bond_length)
                        bond_strength *= (-self.k_factor*4.184/6.022)
                        if bond_strength > 0:
                            bond = bagpype.molecules.Bond([], atom1, atom2, 
                                                                bond_strength, 
                                                                'ELECTROSTATIC')
                            self.bonds.append(bond)
                    except(KeyError):
                        print("WARNING: Charge of atom {0} from residue {1} ".format(atom1.name, 
                                                                            atom1.res_name) + 
                                "or atom {0} from residue {1} ".format(atom2.name, 
                                                                        atom2.res_name) + 
                                "could not be found in the charge database, please update. Ignoring this LINK entry for now.")



# Functions that are used across different bond types and generally a lot throughout the script

def in_same_residue(atom1, atom2):
    return (atom1.res_num == atom2.res_num and
            atom1.chain == atom2.chain and
            atom1.res_name == atom2.res_name)

def distance_between_two_atoms(atom1, atom2):
    return np.asscalar( np.linalg.norm(atom1.xyz - atom2.xyz) )

def sec_neighborhood(G, node):
    """ Returns the second neighbourhood of a node in a Graph G as a set;
        does not include the node itself
    """
    sec_nbhood = set()
    for n in G[node]:
        sec_nbhood.update(G.neighbors(n))
    sec_nbhood.remove(node)
    return list(sec_nbhood)

def within_cov_bonding_distance(atom1, atom2):
    """ Checks distance between two atoms is within covalent
    bonding distance for that pair of elements.  Uses the
    parameters specified in the distance_cutoffs dict.
    """
    try:
        cutoff = bagpype.parameters.distance_cutoffs[(atom1.element, atom2.element)]
    except KeyError:
        try:
            cutoff = bagpype.parameters.distance_cutoffs[(atom2.element, atom1.element)]
        except KeyError:
            raise KeyError("The distance cutoff between " + atom1.element + 
                " and " + atom2.element + " does not exist in the parameters. Please add the relevant value.")
    return True if distance_between_two_atoms(atom1, atom2) < cutoff else False


def in_third_neighbourhood(G,source_atom, target_atom):
    """ check if a node is within distance at most d=3 from another node
    """
    length = nx.single_source_shortest_path_length(G, source_atom.id, cutoff=3)
    return True if target_atom.id in length else False

def uniquify(bond_list):
    """Takes a list of bonds which may contain multiple bonds between
    the same two atoms and returns a list where these have been combined 
    into a single bond with the interaction energies summed.

    Parameters
    ----------
    bond_list : BondList 
      List of non-unique bonds

    Returns
    -------
    bonds_unique : BondList 
      List of bonds where separate Bonds between the same two atoms are
      combined into a single Bond with the interaction energies summed
    """

    bonds_unique = bagpype.molecules.BondList()
    count = 0
    seen = {}
    for bond in bond_list:
        try: 
            bonds_unique[seen[bond]].add_interaction(bond.weight, bond.bond_type)
        except:
            seen[bond] = count*1
            bonds_unique.append(bond)
            count += 1
    return bonds_unique
