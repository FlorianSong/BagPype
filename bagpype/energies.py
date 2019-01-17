# Script to generate the amino acid covalent bond energies dictionary in "energies.py".
# Florian Song, Feb 2017
# The needed library can be downloaded at:
# ftp://ftp.ebi.ac.uk/pub/databases/msd/pdbechem/files/mmcif.tar.gz
# After the download, unpack the tar and deposit in the same dir as this script.
# References for bond energies: 
# Huheey - Inorganic Chemistry 
# Yu-Ran Luo - Comprehensive handbook of chemical bond energies

import bagpype.settings
import bagpype.parameters
import shlex
from CifFile import ReadCif

def generate_energies_dictionary(AA):
    """
    Input: LIST of amino acids or any residue identifier
    Output: Dictionary of bond energies as required by parsing.py
    """
    final_final_dict = {}

    print(("    Generating covalent bond energies for " + ", ".join(AA)))

    # Single bond energies in kJ/mol taken from parameters.py
    singlebond_energies = dict()
    single_bond_energies = bagpype.parameters.single_bond_energies
    for key1 in single_bond_energies:
        partners = single_bond_energies[key1]
        for key2 in partners:
            singlebond_energies[(key1, key2)] = single_bond_energies[key1][key2]
    
    # Double bond energies in kJ/mol
    doublebond_energies = {
    ('C', 'C'): 602,
    ('C', 'N'): 615,
    ('C', 'O'): 799,
    ('P', 'O'): 544,
    ('S', 'O'): 517,
    ('N', 'O'): 607,
    ('P', 'S'): 335,
    }

    triplebond_energies = {
    ("C", "N"): 887,
    ("C", "O"): 1071,
    }


    if type(AA) == "str":
        AA = [AA]

    for aa in AA:

        try:
            file = open(bagpype.settings.DEPENDENCIES_ROOT +'/mmcif/' + aa + '.cif', "r")
        except IOError:
            print(("The residue named " + str(aa) + " was not found in the cif file library of residues. Please add it manually in 'parameters.py'."))
            continue


        ################################
        # Parse corresponding CIF file #
        ################################

        # # Read the cif file in & convert information to a dictionary
        # # Written myself, so needs refining, see DNA faulty line
        # cif_parsed = {}
        # for line in file:
        #     if line.startswith("_chem_comp"):
        #         entries =  shlex.split(line)
        #         cif_parsed[entries[0]] = entries[1]
        #     elif line.startswith("#"):
        #         pass

        #     elif line.startswith("loop_"):
        #         loop_list = []
        #         for line in file:
        #             if line.startswith("#"):
        #                 break
        #             else:
        #                 if line.startswith(';'):
        #                     continue
        #                 elif line.startswith("_"):
        #                     cif_parsed[line.strip()] = []
        #                     loop_list.append(line.strip())
        #                 else:
        #                     splitline = shlex.split(line)
        #                     if not len(splitline) == len(loop_list):
        #                         if aa not in ["DT", "DC", "DA", "DG", "A", "T", "C", "G"]:
        #                             print("Something went wrong in the loops! This is the faulty line:")
        #                             print(line)
        #                     for i in range(len(splitline)):
        #                         cif_parsed[loop_list[i]].append(splitline[i])
        parsed_cif = ReadCif(file)
        lower_case_aa = aa.lower()


        # Extract covalent bond information from the dictionary
        atom_list = parsed_cif[lower_case_aa]['_chem_comp_atom.atom_id']

        try:
            bond_list1 = parsed_cif[lower_case_aa]['_chem_comp_bond.atom_id_1']
            bond_list2 = parsed_cif[lower_case_aa]['_chem_comp_bond.atom_id_2']
            bond_types = parsed_cif[lower_case_aa]['_chem_comp_bond.value_order']
        except KeyError:
            print("    Just to check: ", aa, " is a single atom without bonds.")
            continue


        #####################################################
        # Building of the bond and bond-strength dictionary #
        #####################################################

        # Initialisation
        final_dict = {}
        bond_list_for_cycle_detection = []


        # Initialise all keys, ie every atom
        for atom in atom_list:
            final_dict[atom] = {}

        # Quick check of bond_lists, just a cautionary thing so that the next bit of code can work fine
        if len(bond_list1) != len(bond_list2):
            raise Exception("Something went wrong with the bonds!")

        # Loop over bonds in the bond_list
        for i in range(len(bond_list1)):

            # Extract information, ie both atom names, bond type and both atom elements
            atom1 = bond_list1[i]
            atom2 = bond_list2[i]
            bond_type = bond_types[i]
            atom1_element = parsed_cif[lower_case_aa]['_chem_comp_atom.type_symbol'][parsed_cif[lower_case_aa]['_chem_comp_atom.atom_id'].index(atom1)]
            atom2_element = parsed_cif[lower_case_aa]['_chem_comp_atom.type_symbol'][parsed_cif[lower_case_aa]['_chem_comp_atom.atom_id'].index(atom2)]

            # This is a tuple for looking up the bond strength from the dictionaries defined above
            # Obviously the reverse is also needed, since the dictionaries above are only one way
            lookup = tuple([atom1_element, atom2_element])
            lookup_reversed = tuple([atom2_element, atom1_element])

            # Look up bond strength depending on bond type
            if bond_type == "SING":
                try:
                    bond_strength = singlebond_energies[lookup]
                except KeyError:
                    try:
                        bond_strength = singlebond_energies[lookup_reversed]
                    except KeyError:
                        raise Exception("Please add the single bond between " + atom1_element + " and " + atom2_element)
            elif bond_type == "DOUB":
                try:
                    bond_strength = doublebond_energies[lookup]
                except KeyError:
                    try:
                        bond_strength = doublebond_energies[lookup_reversed]
                    except KeyError:
                        raise Exception("Please add the double bond between " + atom1_element + " and " + atom2_element)
            elif bond_type == "TRIP":
                try:
                    bond_strength = triplebond_energies[lookup]
                except KeyError:
                    try:
                        bond_strength = triplebond_energies[lookup_reversed]
                    except KeyError:
                        raise Exception("Please add the triple bond between " + atom1_element + " and " + atom2_element)
            else:
                raise Exception("Something went wrong with the bond types!")


            # Add bond to the dictionary both ways
            final_dict[atom1][atom2] = final_dict[atom2][atom1] = bond_strength


            # Build a list of bonds for cycle (aromatic ring) detection later on
            # In order to limit computational expense, bonds with Hs are excluded, as they will never be part of aromatic rings
            if atom1_element != 'H' and atom2_element != 'H':
                bond_list_for_cycle_detection.append([atom1, atom2])

        # Usually, DNA endpoints have a H atom instead of another phosphate link. This takes that into account:
        if parsed_cif[lower_case_aa]["_chem_comp.type"] == "DNA LINKING":
            final_dict["HO5'"] = {}
            # final_dict["O5'"] = {}
            final_dict["O5'"]["HO5'"] = final_dict["HO5'"]["O5'"] = 459



        ##################
        # Bond averaging #
        ##################
        def average_bond(atom, neighbours, bonds_dictionary):
            """ Input: 
            atom = single atom_name
            neighbours = 2-element list of neighbours
            bonds_dictionary = dictionary where the bonds are to be changed
            Output: returns updated bonds_dictionary
            """
            original_bond_strength1 = bonds_dictionary[atom][neighbours[0]]
            original_bond_strength2 = bonds_dictionary[atom][neighbours[1]]
            new_bond_strength = (original_bond_strength1 + original_bond_strength2) / 2.
            bonds_dictionary[atom][neighbours[0]] = new_bond_strength
            bonds_dictionary[atom][neighbours[1]] = new_bond_strength
            bonds_dictionary[neighbours[0]][atom] = new_bond_strength
            bonds_dictionary[neighbours[1]][atom] = new_bond_strength

            return bonds_dictionary


        ###########################
        # Aromatic Ring Detection # 
        ###########################

        def find_cycles(list_of_nodes):
            """ Function for finding cycles in undirected graphs 
            Taken from networkx's cycle_basis()
            Why not use networkx one may ask? This way, the output can be better manipulated and the source code is right here for reference
            Input: A list of lists, each sublist containing two atom names.
            Output: A list of all cycles (themselves in the form of lists), ie in this case all aromatic rings!
            """
            unique = set([i for sub in list_of_nodes for i in sub])
            cycles = []
            while unique:
                root = unique.pop()
                stack = [root]
                pred = {root:root}
                used = {root:set()}
                while stack:
                    z = stack.pop()
                    zused = used[z]

                    neighbours = [item for sub in [[j for j in i if j is not z] for i in list_of_nodes if z in i] for item in sub]
                    for nbr in neighbours:
                        if nbr not in used:   # new node
                            pred[nbr]=z
                            stack.append(nbr)
                            used[nbr]=set([z])
                        elif nbr == z:        # self loops
                            cycles.append([z])
                        elif nbr not in zused:# found a cycle
                            pn=used[nbr]
                            cycle=[nbr,z]
                            p=pred[z]
                            while p not in pn:
                                cycle.append(p)
                                p=pred[p]
                            cycle.append(p)
                            cycles.append(cycle)
                            used[nbr].add(z)
                unique-=set(pred)
                root=None
            return [i for i in cycles if len(i) > 1]
            # return [ item for sub in [i for i in cycles if len(i) > 1] for item in sub]

        # Run above function and find aromatic rings
        aromatic_rings = find_cycles(bond_list_for_cycle_detection)

        # Change corresponding entries in the dictionary
        for ring in aromatic_rings:
            # check whether it's 6*C
            if len(ring) == 6 and all([i[0] == 'C' for i in ring]):
                for i, atom in enumerate(ring):
                    if i%2 == 0:
                        neighbours = [i for i in list(final_dict[atom].keys()) if i in ring]
                        final_dict = average_bond(atom, neighbours, final_dict)


        ################################################
        # Carboxylic acid and Guanidin group detection #
        ################################################

        # Please put in all the amino acids that should have averaging HERE:
        amino_acids_for_nonring_bond_averaging = ['ARG', 'ASP', 'GLU', 'KCX']

        # Detection of the specific groups that should be averaged
        if parsed_cif[lower_case_aa]['_chem_comp.id'] in amino_acids_for_nonring_bond_averaging:

            # Find all nontrivial Carbons (ie not "C") which have not been flagged as aromatic in the cif file
            nontrivial_C = [i for i in list(final_dict.keys()) if i != 'C' and i[0] == 'C' and parsed_cif[lower_case_aa]["_chem_comp_atom.pdbx_aromatic_flag"][atom_list.index(i)] == 'N']

            for carbon in nontrivial_C:
                neighbours2 = [i for i in list(final_dict[carbon].keys()) if i[0] not in ['C', 'H']] 
                if len(neighbours2) > 1:

                    # Select the two nonC atoms with the same first two letters, not sure whether this is always the case, but works for the three amino acids
                    neighbours2 = [i for i in neighbours2 for j in [x for x in neighbours2 if x != i] if i[0:2] == j[0:2]] # This is wayyy too long...

                    # Finally, average the bonds:
                    final_dict = average_bond(carbon, neighbours2, final_dict)


        ############################
        # Include "old" atom_names #
        ############################
        # Note that old here means the "old" nomenclature not yet included in the above dictionary and "new" is the one with which the above was generated


        # old_atom_names = parsed_cif[lower_case_aa]['_chem_comp_atom.alt_atom_id']

        # difference_new_list = [x for i,x in enumerate(atom_list) if x != old_atom_names[i]]
        # difference_old_list = [old_atom_names[i] for i,x in enumerate(atom_list) if x != old_atom_names[i]]

        # print(lower_case_aa)
        # print(difference_new_list)
        # print(difference_old_list)

        # for i, new_atom in enumerate(difference_new_list): # new nomenclature

        #     old_atom = difference_old_list[i] # old nomenclature

        #     neighbours3 = list(final_dict[new_atom].keys()) # new nomenclature

        #     print(new_atom, old_atom)
        #     print(neighbours3)
        #     print(final_dict)

        #     previous_dict = final_dict[new_atom] # new nomenclature
        #     dict_to_add = {}

        #     for nb in neighbours3: # new nomenclature
        #         if nb in difference_new_list:

        #             which_position = difference_new_list.index(nb) 
        #             old_name = difference_old_list[which_position] # old nomenclature
                    
        #             # add same bond strength
        #             dict_to_add[old_name] = previous_dict[nb]

        #         else: 

        #             # Don't forget to change the corresponding entries going the other way
        #             final_dict[nb][old_atom] = final_dict[nb][new_atom]

        #             # Add the same entry pretty much, since "nb" doesn't have a old nomenclature name
        #             dict_to_add[nb] = previous_dict[nb]

        #     # Finally add the dictionary with old names to supplement the general one
        #     final_dict[old_atom] = dict_to_add

        # Put the resulting dictionary into a dictionary with the id at the front and then add that to the outputs
        to_add_to_final_final_dict = {}
        to_add_to_final_final_dict[parsed_cif[lower_case_aa]['_chem_comp.id']] = final_dict
        final_final_dict.update(to_add_to_final_final_dict)

    ###################
    # Post-processing #
    ###################

    return final_final_dict

