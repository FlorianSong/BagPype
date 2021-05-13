import bagpype
import os

# This is obsolete, but left in the package because of compatibility with further packages
DEPENDENCIES_ROOT = os.path.dirname(os.path.dirname(bagpype.__file__)) + "/bagpype/dependencies/"

# Capping parameter - number of decimals to cap to
# This is to improve numeric stability across OS's and computers
capping_decimals = 10
