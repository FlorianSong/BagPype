import bagpype.molecules
import bagpype.parameters

import bagpype.covalent
import bagpype.construction
import bagpype.parsing

import bagpype.settings  # for compatibility with the server for now
from bagpype.molecules import Protein
from bagpype.parsing import PDBParser
from bagpype.construction import Graph_constructor


head_string = (
    u"\n"
    u" |\u203E|                Welcome to              \n"
    r" | |__   __ _  __ _ _ __  _   _ _ __   ___  " + "\n"
    r" | '_ \ / _` |/ _` | '_ \| | | | '_ \ / _ \ " + "\n"
    r" | |_) | (_| | (_| | |_) | |_| | |_) |  __/ " + "\n"
    r" |_.__/ \__,_|\__, | .__/ \__, | .__/ \___| " + "\n"
    r"               __/ | |     __/ | |          " + "\n"
    r"              |___/|_|    |___/|_|          " + "\n"
    u"Biochemical, Atomistic Graph construction   \n"
    u"software in PYthon for Proteins, Etc.       \n"
    u"(C) Yaliraki group @ Imperial College London\n"
    u""
)

print(head_string)
# Finish importing with a random quote.
import requests
try:
    response = requests.get("https://api.quotable.io/random")
    if response.status_code == 200:
        response = response.json()
        print('"' + response["content"] + '"') 
        print("by " + response["author"] + "\n" )
except:
    pass
################################################


def standard_run(pdb_file):

    myprot = bagpype.molecules.Protein()

    parser = bagpype.parsing.PDBParser(pdb_file, download=True)
    parser.parse(myprot, strip={"res_name": ["HOH"]})

    ggenerator = bagpype.construction.Graph_constructor()
    ggenerator.construct_graph(myprot)

    return myprot
