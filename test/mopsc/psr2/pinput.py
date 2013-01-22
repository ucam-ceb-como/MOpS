# PLACE CONSTANTS HERE FOR THE BATCH/PSR SOURCE FILE.
# CURRENTLY SET-UP FOR CASE D

######################################################################
# SET NUMERICAL & SYSTEM PARAMETERS HERE
######################################################################
# NMAX is the number of particle equations to solve
NMAX = 64
# TAU is the process time to integrate over
TAU  = 10
# STEPS is the number of steps to take over the process time
# (Assumes constant step spacing)
STEPS = 100
# OUT is the output filename
OUT = "results.csv"
# RES is the number of particles written to MOPS per entry.
# (initialised distributions only)
RES = 32
# RESTIME is the residence time of the reactor, assumed to be V/Q
# (activates a PSR if specified.)
RESTIME = 1.0

######################################################################
# SET PROCESS PARAMETERS HERE
######################################################################
# K is the coagulation kernel value (assumed a function of time)
K = lambda t: t*0.0 + 0.0
# I is the inception rate of type 1 parts (assumed a function of time)
I = lambda t: t*0.0 + 0.0
# S is the surface reaction rate (assumed a function of time)
S = lambda t: t*0.0 + 0.0
# OX is the oxidation reaction rate (assumed a function of time)
OX = lambda t: t*0.0 + 0.0
# NINIT is the initial conditions, as a normal Python list
# Length of NINIT must be the same as NMAX.
NINIT = [0.0]
# or generate a lognormal distribution with gmean and gstdev (>1) given.
GMEAN = 0.0
GSTDEV = 0.0
# NIN is the inflow conditions, doesn't have to be the same as NINIT
# (needed for PSR)
NIN = [1.0]
