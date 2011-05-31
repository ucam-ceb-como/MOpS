# Script to simulate constant pressure reactor for a given time it spits out
# temps and species concentrations at all the time steps.

# 17 NOV 2009 Started reactor simulator based on Equilibrium.py and an
# old reactor simulator I wrote for the TEOS work - ras81 


from Cantera import *
from Cantera.Reactor import *
from Cantera.Func import *
from math import *
import random
#from pylab import *
from scipy import *
import numpy


##### USER VARIABLES ##########################################################

## PROCESS PARAMETERS
P1 = 10132.5 #Pressure Pa 
Tini = 2000.0    #Temperature K
time = 0.16     # total time to run s
dt = 0.001   # step size s

Ntimesteps = int(time/dt + 1)   #Num timesteps + 1 for 0s


## DATA FROM CTI FILES
#cti_file = 'chemOHPref.cti' # cti file

# SET INITIAL CONCENTRATIONS
#  the following expects mole fractions:
initialcomposition = 'H2:0.625, O2:0.375'

plotgraphs = False  #plots concentrations against T
LowC = 1E-60   #Lowest concentraion on log plots

molefractions = True    # If False writes mass fractions

##### END OF USER VARIABLES ##########################################################

## Get all species names out of cti file - following method gets from thermo if it is formatted as in Raphs
# polygenerator script
#chemfile = open('./' + cti_file, 'r')
chemfile = open('h2o2.xml', 'r')
speciesname = []
species = []
for line in chemfile:
    
    #if line[0:14] == 'species(name =':
    #    speciesname = line.split()[2]
    #    speciesname = speciesname.strip('",')
    #    species = species + [speciesname]
    #else:
    #    continue
    #print line[0:19]
    
    if line.startswith('    <species name='):
        
        speciesname = line.split('"')[1]
        print speciesname
        species = species + [speciesname]
    else:
        continue
chemfile.close()

print 'species = ', species



###############################

#### DEFINING VARIABLES #####

## Define number of species
num_species = len(species)
print 'there are', num_species, 'species'

## Import cti file with species data
#gas = importPhase(cti_file, 'gas')
gas = importPhase('h2o2.xml')

## Define array to hold species index
ispecies = range(num_species)

## Define array to hold mole-fractions [this is actually a 2D-array]
mf_species = range(Ntimesteps)

################################

#### Get values of species index
for i in range(0, num_species):
    ispecies[i] = gas.speciesIndex(species[i])
    #print ispecies[i]

#### WRITE CSV FILE ####
    
file2 = open('results.csv','w')
speciesandtemp = range(num_species+1)
speciesandtemp = ['Time (s)'] + ['Temp (K)'] + ['Pressure (Pa)'] + species 
writeCSV(file2, speciesandtemp )

## Define a two-dimensional array to hold mole-fractions and temps of each species at all time steps
#Changed to + 3 to add pressure into the mix ...
num_speciesplusthree = num_species + 3   # + 1 for time and  + 1 for temp and + 1 for pressure
for i in range(Ntimesteps):
    mf_species[i] = range(0, num_speciesplusthree)
    

## set initial conditions!
gas.set(T = Tini, P = P1, X = "%s" %(initialcomposition) )


######################### SET UP THE REACTOR
r   = Reactor(gas, energy ="off")

gas_b = Air()
gas_b.set(P = P1)
env = Reservoir(gas_b)


# Define a wall between the reactor and the environment, and
# make it flexible, so that the pressure in the reactor is held
# at the environment pressure.
w = Wall(r,env)
w.set(K = 1.0e6)   # set expansion parameter. dV/dt = KA(P_1 - P_2)
w.set(A = 1.0)
######################### SET UP THE REACTOR

sim = ReactorNet([r])
time = 0.0
tim = zeros(Ntimesteps,'d')
#Temp = range(Ntimesteps)


for run in range(Ntimesteps):
    tim[run] = time
    mf_species[run][0] = time
    mf_species[run][1] = gas.temperature()
    #New bit
    mf_species[run][2] = gas.pressure()
    #Changed from two to three ...
    for i in range(3, num_speciesplusthree):
        if molefractions == True:
            mf_species[run][i] = gas.moleFraction(ispecies[i-3])  # mole fractions
        else:
            mf_species[run][i] = gas.massFraction(ispecies[i-3])# mass fractions
    writeCSV(file2,   mf_species[run] )
    time += dt
    sim.advance(time)

 
file2.close()

###############################

##### PLOT MOLE_FRACTION Vs T #######

## Define an array with line_types for the graph [worst job ever!!!]
line_type = ('-b', '-g', '-r', '-c', '-m', '-y', '-k',
             '--b', '--g', '--r', '--c', '--m', '--y', '--k',
             '-.b', '-.g', '-.r', '-.c', '-.m', '-.y', '-.k',
             ':b', ':g', ':r', ':c', ':m', ':y', ':k',
             'ob', 'og', 'or', 'oc', 'om', 'oy', 'ok',
             '+b', '+g', '+r', '+c', '+m', '+y', '+k',
             'xb', 'xg', 'xr', 'xc', 'xm', 'xy', 'xk',
             'vb', 'vg', 'vr', 'vc', 'vm', 'vy', 'vk',
             'db', 'dg', 'dr', 'dc', 'dm', 'dy', 'dk',
             '1b', '1g', '1r', '1c', '1m', '1y', '1k',
             '2b', '2g', '2r', '2c', '2m', '2y', '2k',
             '3b', '3g', '3r', '3c', '3m', '3y', '3k',
             '4b', '4g', '4r', '4c', '4m', '4y', '4k',
             'pb', 'pg', 'pr', 'pc', 'pm', 'py', 'pk',
             'hb', 'hg', 'hr', 'hc', 'hm', 'hy', 'hk',
             '.b', '.g', '.r', '.c', '.m', '.y', '.k',
             '^b', '^g', '^r', '^c', '^m', '^y', '^k',
             '<b', '<g', '<r', '<c', '<m', '<y', '<k',
             '>b', '>g', '>r', '>c', '>m', '>y', '>k',
             'sb', 'sg', 'sr', 'sc', 'sm', 'sy', 'sk',
             'Db', 'Dg', 'Dr', 'Dc', 'Dm', 'Dy', 'Dk',
             'Hb', 'Hg', 'Hr', 'Hc', 'Hm', 'Hy', 'Hk',
             ',b', ',g', ',r', ',c', ',m', ',y', ',k',
             '_b', '_g', '_r', '_c', '_m', '_y', '_k',
             '|b', '|g', '|r', '|c', '|m', '|y', '|k'
             )


if plotgraphs == True:
    timeslist = range(Ntimesteps)
    for run in range(Ntimesteps):
        timeslist[run] = mf_species[run][0]

    speciesconcentrations = range(Ntimesteps)
    for i in range(num_species):
        for run in range(Ntimesteps):
            speciesconcentrations[run] = mf_species[run][i+1]
        semilogy(timeslist, speciesconcentrations, line_type[i], label=species[i])
    axis([0, time, LowC, 1])
    xlabel('time (s)')
    ylabel('mol fraction')

    figure()
    
    for i in range(num_species):
        for run in range(Ntimesteps):
            speciesconcentrations[run] = mf_species[run][i+1]
        semilogy(timeslist, speciesconcentrations, line_type[i], label=species[i])
    axis([0, time, LowC, 1])
    legend(loc='best')
    xlabel('time (s)')
    ylabel('mol fraction')

    show()



##    -     # solid line
##    --    # dashed line
##    -.    # dash-dot line
##    :     # dotted line

##    b  # blue
##    g  # green
##    r  # red
##    c  # cyan
##    m  # magenta
##    y  # yellow
##    k  # black
##    w  # white

##    .     # points
##    ,     # pixels
##    o     # circle symbols
##    ^     # triangle up symbols
##    v     # triangle down symbols
##    <     # triangle left symbols
##    >     # triangle right symbols
##    s     # square symbols
##    +     # plus symbols
##    x     # cross symbols
##    D     # diamond symbols
##    d     # thin diamond symbols
##    1     # tripod down symbols
##    2     # tripod up symbols
##    3     # tripod left symbols
##    4     # tripod right symbols
##    h     # hexagon symbols
##    H     # rotated hexagon symbols
##    p     # pentagon symbols
##    |     # vertical line symbols
##    _     # horizontal line symbols
##    steps # use gnuplot style 'steps' # kwarg only
## 
##The following color abbreviations are supported::
## 



