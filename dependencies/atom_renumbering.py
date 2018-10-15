# -*- coding: utf-8 -*-
# This function renumbers the ATOM and HETATM records following the output from reduce such that they can be treated by FIRST

import sys
import array

if (len(sys.argv) > 1):
  print sys.argv[1] 
  fin = open(sys.argv[1], 'r')
  number=0;
  for line in fin :
    temp = array.array('c', line)
    if (str(line[0:4]) == 'ATOM' or str(line[0:4]) == 'HETA'):
      if (str(line[16]) == 'A' or str(line[16]) == ' '):
        number = number+1
        temp[6+(5-len(str(number))):11] =  array.array('c', str(number))
        dummy='     ';
        temp[6:6+(5-len(str(number)))] = array.array('c', dummy[0:(5-len(str(number)))])
      print (temp.tostring()).rstrip('\n')
    else:
      print (temp.tostring()).rstrip('\n')
  fin.close();
    
      
