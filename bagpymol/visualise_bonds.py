from pymol import cmd
import sys
import csv

def visualise(list_of_types = ["COVALENT", "HYDROGEN", "SALTBRIDGE", "HYDROPHOBIC", "STACKED", "BACKBONE", "ELECTROSTATIC"], file_name = "bonds.csv", specific_residue = None):
    list_of_types = [x.upper() for x in list_of_types] if type(list_of_types) is list else [list_of_types.upper()]
    print("Visualising: ", list_of_types, " edges.")

    covalent_colour = [187,187,187] #BBBBBB
    hydrogen_colour = [0,119,187] #0077BB
    saltbridge_colour = [51,187,238] #33BBEE
    hydrophobic_colour = [238,119,51] #EE7733
    stacked_colour = [204,51,17] #CC3311
    backbone_colour = [238,51,119] #EE3377
    electrostatic_colour = [0,153,136] #009988
    other_colour = [0,0,0] #FFFFFF

    cmd.set_color("COVALENT", covalent_colour)
    cmd.set_color("HYDROGEN", hydrogen_colour)
    cmd.set_color("HYDROPHOBIC", hydrophobic_colour)
    cmd.set_color("SALTBRIDGE", saltbridge_colour)
    cmd.set_color("STACKED", stacked_colour)
    cmd.set_color("BACKBONE", backbone_colour)
    cmd.set_color("ELECTROSTATIC", electrostatic_colour)
    cmd.set_color("OTHER", other_colour)

    # file_name = "bonds.csv"
    prefix = "edge"

    f = open(file_name, 'r')
    reader = csv.reader(f)
    """with open(file_name, 'r') as f:
        header = f.readline()
        lines = f.readlines()[1:]"""
    header = reader.next()
    print("HEADER:", header)

    for row in reader:
        #row = row.split(",")
        bond_id, bond_type, atom1_id, atom2_id = row[0], row[1], str(int(row[4])+1), str(int(row[9])+1)

        # print(bond_type.split(","))
        bond_type = bond_type.split(",")

        if len(bond_type) == 1:

            if bond_type[0] in list_of_types:
                cmd.distance(prefix + bond_id,
                        'id ' + atom1_id,
                        'id ' + atom2_id)
                cmd.color(bond_type[0], prefix + bond_id)
            

        elif len(bond_type) == 2:

            type1, type2 = tuple(bond_type)
            if type1 in list_of_types or type2 in list_of_types:

                rgb1 = eval(type1.lower() + "_colour")
                rgb2 = eval(type2.lower() + "_colour")

                rgb_avg = [(x+y)/2 for x,y in zip(rgb1, rgb2)]
                colour = "TEMP_AVERAGE"
                cmd.set_color("TEMP_AVERAGE", rgb_avg)

                cmd.distance(prefix + bond_id,
                        'id ' + atom1_id,
                        'id ' + atom2_id)
                cmd.color(colour, prefix + bond_id)

        if bond_type[0] in ["COVALENT", "ASDF"] and bond_type[0] in list_of_types:
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