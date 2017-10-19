# -*- coding: utf-8 -*-
"""

The objective of this script is to parse the chem.inp file 
and generate a m4 template for setting non-unity lewis numbers 
for any mechanism.


#use default file names:
python m4XMLLewisNumberInputfile.py chem.inp
#specify file name 
python removeTabs.py chem.inp lewisTemplate.m4

then:
m4 lewisTemplate.m4 > LewisNumbersInput.xml 

@author: eb656
"""
import sys
import os
import numpy as np
import math


def isConvertableToFloat(value):
  try:
    float(value)
    return True
  except:
    return False


#sys argument import 
pathfile=os.getcwd()

#check for input filenames
if len(sys.argv) == 3:
    filename=str(sys.argv[1])
    fout=str(sys.argv[2])
else:
    filename='chem.inp'
    fout='lewisTemplate.m4'
File=os.path.join(pathfile, filename)
Fout=os.path.join(pathfile,fout)

#IO solution
F= open(Fout,'w')
F.write("define(setLewis,1.0)\n")
F.write("<Camflow>\n")
F.write("\t<Lewis>\n")
with open(File,'r') as f:
    linNum=0
    section=0
    for line in f:
         # break up line into words
        wordList=line.split()
        
        if len(wordList)==0:
            continue
        elif (wordList[0] =='!'):
            continue
        # Check Section:
        elif (   wordList[0] == "ELEMENTS" 
            or wordList[0] == "SPECIES" 
            ):
            section +=1
        elif (wordList[0] == "REACTION" 
            or wordList[0] == "REACTIONS"):
            section +=1
            continue # go to next line
        
        #print(line)    
        #Add OR logic here to leave other sections unchanged
        if section == 2:
            if wordList[0] == '!':
                print("found comment")
                # continue # to next line
            
            else:
                for word in wordList:
                    if word == "!":
                        break
                    print(word)
                    if (word != "SPECIES" and word != "END"):
                        F.write("\t\t<species name=\""
                                +word
                                +"\">setLewis</species>\n")
# Clean up end of file:
F.write("\t</Lewis>\n</Camflow>")
F.close()