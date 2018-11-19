def visualise(list_of_types = ["COVALENT", "HYDROGEN", "HYDROPHOBIC", "STACKED", "BACKBONE", "ELECTROSTATIC"], filename = "bonds.csv"):
    import sys
    import csv

    covalent_colour = [251,101,66]
    hydrogen_colour = [55,94,151]
    hydrophobic_colour = [255,187,0]
    dna_colour = [63,104,28]

    cmd.set_color("COVALENT", covalent_colour)
    cmd.set_color("HYDROGEN", hydrogen_colour)
    cmd.set_color("HYDROPHOBIC", hydrophobic_colour)
    cmd.set_color("DNA", dna_colour)


    def determine_colour(bond_type):
        if bond_type in ["COVALENT", "HYDROGEN", "HYDROPHOBIC"]:
            color = bond_type
        elif bond_type in ["STACKED", "BACKBONE"]: 
            color = "DNA"
        else:
            color = "white"
        return(color)

    # filename = "bonds.csv"
    prefix = "bond"

    f = open(filename, 'r')
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
                colour = determine_colour(bond_type[0])
                cmd.distance(prefix + bond_id,
                        'id ' + atom1_id,
                        'id ' + atom2_id)
                cmd.color(colour, prefix + bond_id)

        elif len(bond_type) <3:
            type1, type2 = tuple(bond_type)

            if type1 in list_of_types or type2 in list_of_types:

                colour1, colour2 = determine_colour(type1), determine_colour(type2)

                rgb1 = {"COVALENT":covalent_colour, "HYDROGEN":hydrogen_colour, "HYDROPHOBIC":hydrophobic_colour,
                            "DNA":dna_colour, "white":[0,0,0]}[colour1]
                rgb2 = {"COVALENT":covalent_colour, "HYDROGEN":hydrogen_colour, "HYDROPHOBIC":hydrophobic_colour,
                            "DNA":dna_colour, "white":[0,0,0]}[colour2]

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