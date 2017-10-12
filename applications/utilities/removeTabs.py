# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 13:56:53 2017

SPROGC needs Ea in CAL/MOL for chem.inp
Very Rough brute force to change K --> Cal/Mole

reads in a chem.inp file:
    filename='chem.inp'
it will write a new chem.inp file:
    fout='chemCalPerMol.inp'


run with following commands

#use default file names:
python kelvinsToCal-Mol.py
#specify file name 
python kelvinsToCal-Mol.py chem_kelvin.inp chem_calmol.inp

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
    fout='chemCalPerMol.inp'
File=os.path.join(pathfile, filename)
Fout=os.path.join(pathfile,fout)

#IO solution
F= open(Fout,'w')
F.write("! Parsed with removeTabs.py")
with open(File,'r') as f:
    linNum=0
    section=0
    for line in f:
        linNum +=1
        #print(linNum)
        #print(line)
        # indicate new line:
        i=1

        # break up line into words
        wordList=line.split()
        
        #print(wordList)
        if len(wordList)==0:
            #print("EmptyLine")
            continue
        elif (wordList[0] =='!'):
            F.write(line)
            continue
        # Check Section:
        elif (   wordList[0] == "ELEMENTS" 
            or wordList[0] == "SPECIES" 
            ):
            section +=1
        elif (wordList[0] == "REACTION" 
            or wordList[0] == "REACTIONS"):
            section +=1
            F.write("REACTION") #Uses assumed units of MOLES and CAL/MOLE
            F.write('\n')
            continue # go to next line
        if section == 1 or section == 3:
            F.write(line)
    
        # Edit wordList to edit Ea 
        # if section == 3:
            #Find lines that represent reactions
            #Word one (or more) will represent the reaction string
            #wordList[-1] will be Ea (units of K)
            #wordList[-2] will be b (unitless)
            #wordList[-3] will be A (MOLES)
            # if (len(wordList) >4 
            #     and isConvertableToFloat(wordList[-1])
            #     and isConvertableToFloat(wordList[-2])
            #     and isConvertableToFloat(wordList[-3])
            #     ):
            #     #multiply by R and convert to cal
            #     tmp = float(wordList[-1])*8.134/4.184
            #     #round to 3 places
            #     tmp = math.ceil(tmp*1000)/1000
            #     wordList[-1]=str(tmp)
                
        
        else:
            for word in wordList:
                # Write regular line 
                #   write first word
                if i==1:
                    F.write(str(word))
                    i=2
                #   write trailing words
                else:
                    F.write('  '+str(word))
            # write end of line
            F.write('\n') 
F.close()