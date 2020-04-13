import urllib

class GraphConstructionError(Exception):
    def __init__(self, message, *args, **kwargs):
        message = urllib.parse.quote(message, safe='')
        message = message.replace(".", "%2E")
        Exception.__init__(self, message, *args, **kwargs)


class MissingEnergyError(Exception):
    def __init__(self, atom1_el, atom2_el, residue, valence="SING", *args, **kwargs):
        output = (
            "Sadly, we have no parameters for a covalent {} bond between elements {} and {}."
            " This error was likely caused by one or more residues: {}".format(
                {"SING": "single", "DOUB": "double", "TRIP": "triple"}[valence],
                atom1_el,
                atom2_el,
                residue,
            )
        )
        Exception.__init__(self, output, *args, **kwargs)


class UnknownResidueError(Exception):
    def __init__(self, aa, *args, **kwargs):
        output = (
            "The residue named {} was not found in the cif file library of residues. "
            "Please add it manually in 'parameters.py'."
        ).format(aa)
        Exception.__init__(self, output, *args, **kwargs)


class UnusualHydrogenError(Exception):
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)


class UnusualCIFFile(Exception):
    def __init__(self, ciffile, *args, **kwargs):
        output = (
            "Something is unusual about the cif file relating to "
            "the following three-letter code: {}".format(ciffile)
        )
        Exception.__init__(self, output, *args, **kwargs)
