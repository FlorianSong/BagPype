from pymol import cmd
import sys
import csv

def visualise(list_of_types = ["COVALENT", "DISULFIDE", "HYDROGEN", "SALTBRIDGE", "HYDROPHOBIC", "STACKED", "BACKBONE", "ELECTROSTATIC"], 
              file_name = "bonds.csv", specific_residues = [] # e.g. [('1', 'A')]
              ):
    list_of_types = [x.upper() for x in list_of_types] if type(list_of_types) is list else [list_of_types.upper()]
    print("Visualising: {} edges.".format(", ".join(list_of_types).lower()))

    colours = {"covalent": [187,187,187], #BBBBBB
               "disulfide": [153,153,51], #999933
               "hydrogen": [0,119,187], #0077BB
               "saltbridge": [51,187,238], #33BBEE
               "hydrophobic": [238,119,51], #EE7733
               "stacked": [204,51,17], #CC3311
               "backbone": [238,51,119], #EE3377
               "electrostatic": [0,153,136], #009988
               "other": [0,0,0], #FFFFFF
               }

    for type_of_edge in colours: 
        cmd.set_color(type_of_edge.upper(), colours[type_of_edge])

    # prefix = "edge"

    f = open(file_name, 'r')
    reader = csv.reader(f)
    header = reader.next()

    for row in reader:
        bond_data = {}
        for i, item in enumerate(header):
            bond_data[item] = row[i]

        if len(specific_residues)>0:
            current_residue1 = tuple([ str(bond_data["atom1_res_num"]), str(bond_data["atom1_chain"]) ])
            current_residue2 = tuple([ str(bond_data["atom2_res_num"]), str(bond_data["atom2_chain"]) ])
            if not any([current_residue1 in specific_residues, current_residue2 in specific_residues]):
                continue

        # bond_id, bond_type, atom1_id, atom2_id = row[0], row[1], str(int(row[4])+1), str(int(row[9])+1)
        bond_id = bond_data["bond_id"]
        bond_type = bond_data["bond_type"].split(",")
        atom1_id = str(int( bond_data["atom1_id"]) +1 )
        atom2_id = str(int( bond_data["atom2_id"]) +1 )

        if len(bond_type) == 1:

            if bond_type[0] in list_of_types:
                prefix = bond_type[0].lower()
                cmd.distance(prefix + bond_id,
                        'id ' + atom1_id,
                        'id ' + atom2_id)
                cmd.color(bond_type[0], prefix + bond_id)
            

        elif len(bond_type) == 2:

            type1, type2 = tuple(bond_type)
            if type1 in list_of_types or type2 in list_of_types:

                rgb1 = colours[type1.lower()] # eval(type1.lower() + "_colour")
                rgb2 = colours[type2.lower()] # eval(type2.lower() + "_colour")

                rgb_avg = [(x+y)/2 for x,y in zip(rgb1, rgb2)]
                colour = "TEMP_AVERAGE"
                cmd.set_color(colour, rgb_avg)

                prefix = "".join([type1, type2]).lower()
                print(prefix)
                cmd.distance(prefix + bond_id,
                        'id ' + atom1_id,
                        'id ' + atom2_id)
                cmd.color(colour, prefix + bond_id)

        if bond_type[0] in ["COVALENT", "DISULFIDE"] and bond_type[0] in list_of_types:
            cmd.set("dash_gap", 0, prefix + bond_id)
            cmd.set("dash_radius", 0.15, prefix + bond_id)
        elif bond_type[0] in ["HYDROGEN", "SALTBRIDGE"] and bond_type[0] in list_of_types:
            cmd.set("dash_gap", 0, prefix + bond_id)
            cmd.set("dash_radius", 0.1, prefix + bond_id)
        elif bond_type[0] in list_of_types:
            cmd.set("dash_gap", 0.3, prefix + bond_id)
            cmd.set("dash_radius", 0.1, prefix + bond_id)

    cmd.hide('labels')
    # cmd.set('dash_gap',0)

cmd.extend("visualise", visualise)