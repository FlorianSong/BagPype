import bagpype
import os

DEPENDENCIES_ROOT = os.path.dirname(os.path.dirname(bagpype.__file__)) + "/dependencies/"

# Capping parameter - number of decimals to cap to
# This is to improve numeric stability across OS's and computers
capping_decimals = 10
