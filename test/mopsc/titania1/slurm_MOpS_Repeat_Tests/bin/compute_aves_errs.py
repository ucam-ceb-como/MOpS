#!/bin/python 

from shutil import copyfile
import numpy
import csv 
import math
import sys

# Start
nruns = int(sys.argv[1])
stgid = sys.argv[2]

# Get data from folder 1
chembar  = numpy.genfromtxt('1/Network(stage'+stgid+')-chem.csv', dtype=float, delimiter=',')
nsteps   = chembar.shape[0]-1

chembar  = numpy.genfromtxt('1/Network(stage'+stgid+')-chem.csv', dtype=float, delimiter=',', skip_header=1)
partbar  = numpy.genfromtxt('1/Network(stage'+stgid+')-part.csv', dtype=float, delimiter=',', skip_header=1)
ptempbar = numpy.genfromtxt('1/Network(stage'+stgid+')-part-temp.csv', dtype=float, delimiter=',', skip_header=1)
pratebar = numpy.genfromtxt('1/Network(stage'+stgid+')-part-rates.csv', dtype=float, delimiter=',', skip_header=1)

cheads   = numpy.genfromtxt('1/Network(stage'+stgid+')-chem.csv', dtype=str, delimiter=',', skip_footer=nsteps)
pheads   = numpy.genfromtxt('1/Network(stage'+stgid+')-part.csv', dtype=str, delimiter=',', skip_footer=nsteps)
ptheads  = numpy.genfromtxt('1/Network(stage'+stgid+')-part-temp.csv', dtype=str, delimiter=',', skip_footer=nsteps)
prheads  = numpy.genfromtxt('1/Network(stage'+stgid+')-part-rates.csv', dtype=str, delimiter=',', skip_footer=nsteps)

# Data columns
nchem  = chembar.shape[1]
npart  = partbar.shape[1]
nptemp = ptempbar.shape[1]
nprate = pratebar.shape[1]

# Get data from other folders and add to sums 	
for f in range (2,nruns+1):
    chembar  = chembar  + numpy.genfromtxt(str(f)+'/Network(stage'+stgid+')-chem.csv', dtype=float, delimiter=',', skip_header=1)
    partbar  = partbar  + numpy.genfromtxt(str(f)+'/Network(stage'+stgid+')-part.csv', dtype=float, delimiter=',', skip_header=1)
    ptempbar = ptempbar + numpy.genfromtxt(str(f)+'/Network(stage'+stgid+')-part-temp.csv', dtype=float, delimiter=',', skip_header=1)
    pratebar = pratebar + numpy.genfromtxt(str(f)+'/Network(stage'+stgid+')-part-rates.csv', dtype=float, delimiter=',', skip_header=1)

# Divide by number of runs	
chembar  = numpy.divide(chembar,nruns)
partbar  = numpy.divide(partbar,nruns)
ptempbar = numpy.divide(ptempbar,nruns)
pratebar = numpy.divide(pratebar,nruns)

# Get data from folder 1
chem  = numpy.genfromtxt('1/Network(stage'+stgid+')-chem.csv', dtype=float, delimiter=',', skip_header=1)
part  = numpy.genfromtxt('1/Network(stage'+stgid+')-part.csv', dtype=float, delimiter=',', skip_header=1)
ptemp = numpy.genfromtxt('1/Network(stage'+stgid+')-part-temp.csv', dtype=float, delimiter=',', skip_header=1)
prate = numpy.genfromtxt('1/Network(stage'+stgid+')-part-rates.csv', dtype=float, delimiter=',', skip_header=1)

# Start to compile std devs
chemerr  = pow(chem - chembar,2.0)
parterr  = pow(part - partbar,2.0)
ptemperr = pow(ptemp - ptempbar,2.0)
prateerr = pow(prate - pratebar,2.0)

# Get data from other folders and add to std devs
for f in range (1,nruns+1):
    chem     = numpy.genfromtxt(str(f)+'/Network(stage'+stgid+')-chem.csv', dtype=float, delimiter=',', skip_header=1)
    part     = numpy.genfromtxt(str(f)+'/Network(stage'+stgid+')-part.csv', dtype=float, delimiter=',', skip_header=1)
    ptemp    = numpy.genfromtxt(str(f)+'/Network(stage'+stgid+')-part-temp.csv', dtype=float, delimiter=',', skip_header=1)
    prate    = numpy.genfromtxt(str(f)+'/Network(stage'+stgid+')-part-rates.csv', dtype=float, delimiter=',', skip_header=1)
    chemerr  = chemerr + pow(chem - chembar,2.0)
    parterr  = parterr + pow(part - partbar,2.0)
    ptemperr = ptemperr + pow(ptemp - ptempbar,2.0)
    prateerr = prateerr + pow(prate - pratebar,2.0)

# Divide by number of runs
chemerr  = numpy.divide(chemerr,nruns-1)
parterr  = numpy.divide(parterr,nruns-1)
ptemperr = numpy.divide(ptemperr,nruns-1)
prateerr = numpy.divide(prateerr,nruns-1)

# Square roots
chemerr  = pow(chemerr,0.5)
parterr  = pow(parterr,0.5)
ptemperr = pow(ptemperr,0.5)
prateerr = pow(prateerr,0.5)

# Combine data 
chembar[:,3::2]  = chemerr[:,2::2]
partbar[:,3::2]  = parterr[:,2::2]
ptempbar[:,3::2] = ptemperr[:,2::2]
pratebar[:,3::2] = prateerr[:,2::2]
	
# New file names
cfile  = 'Network(stage'+stgid+')-chem-avgd.csv'
pfile  = 'Network(stage'+stgid+')-part-avgd.csv'
ptfile = 'Network(stage'+stgid+')-part-temp-avgd.csv'
prfile = 'Network(stage'+stgid+')-part-rates-avgd.csv'

with open(cfile,'wb+') as csvfile:
	file = csv.writer(csvfile, delimiter=",")
	file.writerow(cheads)
	for j in range (0,nsteps):
		file = csv.writer(csvfile, delimiter=",")
		file.writerow(chembar[j,:])

with open(pfile,'wb+') as csvfile:
	file = csv.writer(csvfile, delimiter=",")
	file.writerow(pheads)
	for j in range (0,nsteps):
		file = csv.writer(csvfile, delimiter=",")
		file.writerow(partbar[j,:])

with open(ptfile,'wb+') as csvfile:
	file = csv.writer(csvfile, delimiter=",")
	file.writerow(ptheads)
	for j in range (0,nsteps):
		file = csv.writer(csvfile, delimiter=",")
		file.writerow(ptempbar[j,:])

with open(prfile,'wb+') as csvfile:
	file = csv.writer(csvfile, delimiter=",")
	file.writerow(prheads)
	for j in range (0,nsteps):
		file = csv.writer(csvfile, delimiter=",")
		file.writerow(pratebar[j,:])
