import networkx as nx
import numpy as np
import pandas as pd
import bagpype
from scipy.integrate import quad
import time

class constrain(object):
    def __init__(self, generator, FIRST = False, RMST_gamma = -0.2, angle_energy = 0.13, dihedral_energy = 0.04):
        self.generator = generator
        self.sparsify_RMST_gamma = RMST_gamma
        self.sparsify_burn_bridges = False
        self.angle_energy_multiplier = angle_energy
        self.dihedral_energy_multiplier = dihedral_energy
        if FIRST == True:
            self.legacy()
        self.find_total_covalent_energy()

    def legacy(self):
        self.generator.covalent_bonds_graph = nx.Graph()
        for bond in (x for x in self.generator.protein.bonds if x.bond_type == ['COVALENT']):
            self.generator.covalent_bonds_graph.add_edge(bond.atom1.id, bond.atom2.id, weight = bond.weight)

    def find_total_covalent_energy(self):
        self.total_covalent_energy = 0
        for bond in (x for x in self.generator.protein.bonds if x.bond_type == ['COVALENT']):
           self.total_covalent_energy += bond.weight 
        self.total_covalent_energy = self.total_covalent_energy/len([x for x in self.generator.protein.bonds if x.bond_type == ['COVALENT']])

    def uniquify(self, plane_list):
        forward_set = set()
        for plane in plane_list:
            if tuple(plane[::-1]) not in forward_set:
                forward_set.add(tuple(plane))
        return list(forward_set)

    def unpack(self, plane_list):
        if plane_list == []:
            return None

        unpacked = []
        for plane in plane_list:
            unpacked.append([self.generator.protein.atoms[i] for i in plane])
        return unpacked

    def find_angle(self, site1):
        angle_list = []
        lam = lambda site: (x for x in self.generator.covalent_bonds_graph.neighbors(site) if not self.generator.protein.atoms[x].element == 'H')   
        for site2 in lam(site1):
            for site3 in lam(site2):
                if site1 != site3:
                    angle_list.append([site1, site2, site3])

        return self.unpack(angle_list)

    def find_dihedral(self, site1):
        dihedral_list = []
        lam = lambda site: (x for x in self.generator.covalent_bonds_graph.neighbors(site) if not self.generator.protein.atoms[x].element == 'H')
        for site2 in lam(site1):
            for site3 in lam(site2):
                for site4 in lam(site3):
                    if site1 != site3 and site2 != site4:
                        dihedral_list.append([site1, site2, site3, site4])

        return self.unpack(dihedral_list) 

    def sparsify(self):
        print('RMST sparsification started')
        mst = nx.minimum_spanning_tree(self.generator.covalent_bonds_graph, weight = 'weight')
        notmst = self.constraint_graph
        matches = list()

        D = nx.adjacency_matrix(mst, weight = 'weight')
        d = D.max(axis = 0).toarray().flatten().tolist()
        node_list = list(mst.nodes)

        for edge_not_in_mst in notmst.edges:
            i,j = edge_not_in_mst[0], edge_not_in_mst[1]

            path = nx.shortest_path(mst, source = i, target = j)
            weights_along_path = [mst[u][v]["weight"] for u,v in zip(path[:-1], path[1:])]

            mlink = max(weights_along_path)

            left_hand_side = mlink + self.sparsify_RMST_gamma*abs(d[node_list.index(i)] + d[node_list.index(j)])
            right_hand_side = notmst[i][j]["weight"]
            if (left_hand_side > right_hand_side ):
                matches.append([i,j])

        if self.sparsify_burn_bridges == True:
            rmst_graph = nx.Graph(matches)
            bridges = list(nx.bridges(rmst_graph))
            rmst_graph.remove_edges_from(bridges)
            matches = rmst_graph.edges
        
        print("RMST sparsification finished. Accepted: " + str(len(matches)) + ", rejected: " + str(len(self.constraint_graph.edges) - len(matches)) + "; MST size: " + str(len(mst.edges)))
   
        self.record_bonds(matches)
       
    def append_constraints(self, energy, plane):
        if len(plane[0]) == 3:
            TYPE = 'ANGLE'
        else:
            TYPE = 'DIHEDRAL'
        
        for i, E in enumerate(energy):
            if E < 1:
                continue
            else:
                self.constraint_bonds.append(bagpype.molecules.Bond(i, plane[i][0], plane[i][-1], E, TYPE))
                self.constraint_graph.add_edge(plane[i][0].id, plane[i][-1].id, weight = E,
                                           distance = distance_between_two_atoms(plane[i][0], plane[i][-1]), 
                                           energy = E, type = TYPE, plane = plane[i])
        
    def record_bonds(self, matches):
        for bond in matches:
            TYPE = self.constraint_graph[bond[0]][bond[1]]["type"]
            atom1, atom2 = self.generator.protein.atoms[bond[0]], self.generator.protein.atoms[bond[1]]
            bond_strength = self.constraint_graph[bond[0]][bond[1]]["energy"]
            self.generator.bonds.append(bagpype.molecules.Bond([], atom1, atom2, 
                                                          bond_strength, TYPE))
        
        for i, bond in enumerate(self.generator.bonds):
            bond.id = i
            self.generator.protein.graph.add_edge(bond.atom1.id, bond.atom2.id,
                                                  weight = bond.weight, 
                                                  id = i,
                                                  bond_type = bond.bond_type)

        self.generator.protein.bonds = self.generator.bonds

        self.write_bonds_to_csv_file('bonds.csv')

    def write_bonds_to_csv_file(self, name):
        print('Writing bond data to CSV')
        bond_data = []
        for bond in self.generator.protein.bonds:
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
    
    def activate(self, energy_list):
        return np.array(energy_list)/sum(np.array(energy_list))

    def angles(self):
        print('Finding Angles')
        def find_energy(angle, th):
            th0, ka = bagpype.constrain_parameters.angle[angle[0]][angle[1]][angle[2]]
            integrand = lambda th: 0.5*(ka)*((th - th0)**2) 

            return integrand(th)

        angle_list = []
        energy_list = []
        for atom in (node for node in self.generator.covalent_bonds_graph.nodes() if not self.generator.protein.atoms[node].element == 'H'):
            plane = self.find_angle(atom)
            if plane == None:
                continue
            else:
                for i in plane:
                    angle_list.append(i)

        angle_list = self.uniquify(angle_list)

        for i, angle in enumerate(angle_list):
            try:
                ff_names = [bagpype.constrain_parameters.translate[atom.res_name][atom.name] for atom in angle]
                energy_list.append(find_energy(ff_names, angle_between_two_atoms(angle)))
            except:
                del angle_list[i]

        energy_list = self.activate(energy_list)*self.total_covalent_energy*self.angle_energy_multiplier*len(energy_list)
        return energy_list, angle_list

    def dihedrals(self):
        print('Finding Dihedrals')
        def find_energy(dihedral, atoms, om):
            def integrand(om, om0, kt, n):
                om0 = np.radians(om0)
                return kt*(1+np.cos(n*om-om0))

            try:
                param_list = bagpype.constrain_parameters.dihedral[dihedral[0]][dihedral[1]][dihedral[2]][dihedral[3]]
            except KeyError:
                param_list = bagpype.constrain_parameters.dihedral_wildtype[dihedral[1]][dihedral[2]]

            E = 0

            for param in param_list:
                om0, kt, n = param
                E += integrand(om, om0, kt, n)

            return E

        dihedral_list = []
        energy_list = []

        for atom in (node for node in self.generator.covalent_bonds_graph.nodes() if not self.generator.protein.atoms[node].element == 'H'):
            plane = self.find_dihedral(atom)
            if plane == None:
                continue
            else:
                for i in plane:
                    dihedral_list.append(i)

        dihedral_list = self.uniquify(dihedral_list)

        for i, dihedral in enumerate(dihedral_list):
            try:
                ff_names = [bagpype.constrain_parameters.translate[atom.res_name][atom.name] for atom in dihedral]
                energy_list.append(find_energy(ff_names, dihedral, angle_between_two_planes(dihedral)))
            except:
                del dihedral_list[i]

        energy_list = self.activate(energy_list)*self.total_covalent_energy*self.dihedral_energy_multiplier*len(energy_list)
        return energy_list, dihedral_list

    def constrain(self, protein):
        t0 = time.time()
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        print('Adding angle and dihedral edges to graph')
        self.constraint_graph = nx.Graph()
        self.constraint_bonds = bagpype.molecules.BondList()

        self.append_constraints(*self.angles())
        self.append_constraints(*self.dihedrals())
        self.sparsify()

        protein.graph = self.generator.protein.graph
        protein.bonds = self.generator.protein.bonds

        print('Constrain finished. Time taken: ', time.time()-t0)

def angle_between_two_atoms(atom_list):
    v21 = atom_list[1].xyz - atom_list[0].xyz
    v23 = atom_list[1].xyz - atom_list[2].xyz
    return np.rad2deg(np.arccos(np.dot(v21, v23)/(np.linalg.norm(v21) * np.linalg.norm(v23))))

def angle_between_two_planes(atom_list):
    vijk = np.cross((atom_list[2].xyz-atom_list[1].xyz), (atom_list[0].xyz-atom_list[1].xyz))
    vjkl = np.cross((atom_list[1].xyz-atom_list[2].xyz), (atom_list[3].xyz-atom_list[2].xyz))
    return np.arccos(np.dot(vijk, vjkl)/(np.linalg.norm(vijk) * np.linalg.norm(vjkl)))

def distance_between_two_atoms(atom1, atom2):
    return np.asscalar(np.linalg.norm(atom1.xyz - atom2.xyz))        