from pymol import cmd
import sys
import csv

def visualise(list_of_types = ["COVALENT", "HYDROGEN", "HYDROPHOBIC", "STACKED", "BACKBONE", "ELECTROSTATIC"], file_name = "bonds.csv"):
    list_of_types = [x.upper() for x in list_of_types] if type(list_of_types) is list else [list_of_types.upper()]
    print("Visualising: ", list_of_types, " edges.")

    covalent_colour = [251,101,66]
    hydrogen_colour = [55,94,151]
    saltbridge_colour = [127, 152, 189]
    hydrophobic_colour = [255,187,0]
    stacked_colour = [63,104,28]
    backbone_colour = [63,104,28]
    electrostatic_colour = [255,255,255]
    other_colour = [0,0,0]

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

    # header = reader.next()
    # if header != ['', 'atom1', 'atom2', 'bond_name', 'distance', 'pp', 'pp_adjusted', 'pp_normalised', 'weight', 'qs', 'qs_test_set']:
    #     raise ValueError('header is not correct format')
    header = reader.next()

    for row in reader:
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

    cmd.hide('labels')
    cmd.set('dash_gap',0)

cmd.extend("visualise", visualise)