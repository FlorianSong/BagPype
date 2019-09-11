
class GraphConstructionError(Exception):
    def __init__(self,*args,**kwargs):
        Exception.__init__(self,*args,**kwargs)

class MissingEnergyError(Exception):
    def __init__(self, atom1_el, atom2_el, residue, valence = "SING",
                 *args,**kwargs):
        output = "Sadly, we have no parameters for a covalent {} bond between elements {} and {}. This error was likely caused by one or more residues: {}".format(
                        {"SING": "single",
                         "DOUB": "double",
                         "TRIP": "triple"}[valence],
                        atom1_el, atom2_el, residue)
        Exception.__init__(self, output, *args,**kwargs)

class UnknownResidueError(Exception):
    def __init__(self,*args,**kwargs):
        Exception.__init__(self,*args,**kwargs)

class UnusualHydrogenError(Exception):
    def __init__(self,*args,**kwargs):
        Exception.__init__(self,*args,**kwargs)
        
class UnusualCIFFile(Exception):
    def __init__(self, ciffile,
                 *args,**kwargs):
        output = "Something is unusual about the cif file relating to the following three-letter code: {}".format(ciffile)
        Exception.__init__(self, output, *args,**kwargs)
