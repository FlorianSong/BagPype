#!/bin/bash

wget http://kinemage.biochem.duke.edu/php/downlode-3.php?filename=/../downloads/software/reduce31/reduce.3.23.130521.linuxi386.gz
gunzip *.gz
mv *.linuxi386 dependencies/reduce

#wget http://kinemage.biochem.duke.edu/php/downlode-3.php?filename=/../downloads/software/reduce31/reduce.3.23.130521.macosx.zip
#unzip *.macosx.zip
#mv *.macosx dependencies/reduce.macosx
#rm *.macosx.zip

wget ftp://ftp.ebi.ac.uk/pub/databases/msd/pdbechem/files/mmcif.tar.gz
tar -zxf mmcif.tar.gz
mv mmcif/* dependencies/mmcif/
rm -r mmcif*

# echo # added by proteinClass installer >> ~/.bashrc
# export PYTHONPATH="$PWD:$PYTHONPATH"

chmod +x dependencies/reduce

# This bit changes the settings.py file so that no modification there is necessary
# line="DEPENDENCIES_ROOT = '$PWD/dependencies' # Modify this line"
# sed -i "2s|.*|$line|" $PWD/bagpype/settings.py

#python proteingraph/setup.py install
