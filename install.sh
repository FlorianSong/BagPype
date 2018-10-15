#!/bin/bash

wget http://kinemage.biochem.duke.edu/php/downlode-3.php?filename=/../downloads/software/reduce31/reduce.3.23.130521.linuxi386.gz
gunzip *.gz
mv *.linuxi386 dependencies/reduce

# echo # added by proteinClass installer >> ~/.bashrc
# export PYTHONPATH="$PWD:$PYTHONPATH"

chmod +x dependencies/reduce

# This bit changes the settings.py file so that no modification there is necessary
line="DEPENDENCIES_ROOT = '$PWD/dependencies' # Modify this line"
sed -i "2s|.*|$line|" $PWD/biomolgraph/settings.py

# python proteingraph/setup.py install