# Directories for command line calls to reduce
DEPENDENCIES_ROOT = '/home/fv611/Documents/code_shared/bagpype/dependencies' # Modify this line
REDUCE = DEPENDENCIES_ROOT + '/reduce'

# Names for amino acids you would like automatically removed
# aa_to_eliminate = ['SOL', 'NA+', 'CL', 'EDO']
aa_to_eliminate = []


head_string = u"\
    \n\
 |\u203E|                Welcome to              \n\
 | |__   __ _  __ _ _ __  _   _ _ __   ___  \n\
 | '_ \ / _` |/ _` | '_ \| | | | '_ \ / _ \ \n\
 | |_) | (_| | (_| | |_) | |_| | |_) |  __/ \n\
 |_.__/ \__,_|\__, | .__/ \__, | .__/ \___| \n\
               __/ | |     __/ | |          \n\
              |___/|_|    |___/|_|          \n\
Biochemical, Atomistic Graph construction   \n\
software in PYthon for Proteins, Etc.       \n\
(C) Yaliraki group @ Imperial College London\n\
"
