# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 13:56:53 2017

Convert profile.dat camflow output files to CSV

takes two command line arguments
datFileName csvFileName

it will write csvFileName to cwd. 


run with following command

python datToCSV.py datFileName csvFileName
e.g.
python datToCSV.py profile.dat profile.csv

@author: eb656
"""
import sys
import os
import numpy as np

#sys argument import 


pathfile=os.getcwd()

#check for input filenames
if len(sys.argv) == 3:
    filename=str(sys.argv[1])
    fout=str(sys.argv[2])
else:
    filename='profile.dat'
    fout='profile.csv'
File=os.path.join(pathfile, filename)
Fout=os.path.join(pathfile,fout)

#IO solution
F= open(Fout,'w')
with open(File,'r') as f:
    for line in f:
        i=1
        for word in line.split():
            if i==1:
                F.write(str(word))
                i=2
            else:
                F.write(','+str(word))
        F.write('\r\n')
F.close()


#NumPy Implemenataion: Not working on Vienna due to old NumPy version
#Generate header String
#header = np.genfromtxt(File,dtype=str,max_rows=1)
#headerstr=''
#for head in header:
#    if head==header[0]:
#        headerstr=head
#    else:
#        headerstr=headerstr+','+head
#
#    
#data = np.genfromtxt(File,dtype=float,skip_header=1)
#
##write CSV file
##fname, X, fmt='%.18e', delimiter=' ', +
##newline='\n', header='', footer='', comments='# '
#np.savetxt(Fout,data,delimiter=',',header=headerstr,comments='',fmt='%.5E')
