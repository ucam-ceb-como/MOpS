# PYTHON SOURCE FOR COMPARING A MOPS CALCULATION WITH CANTERA.
# THIS FILE AIMS TO CHECK THE CALCULATION OF LOW-PRESSURE LIMITS FOR
# THE MECHANISM OF HO ET AL., 1994

# Copyright 2012 William Menz (wjm34)

import numpy as np
import Cantera as can
import Cantera.Reactor as reac
import matplotlib.pyplot as plt

def getConcs():
    # Multiply by P/RT factor to convert yi into Ci
    mult = reactor.pressure()/(1.0e6*8.315*reactor.temperature())
    conc_sih4.append(mult* reactor.moleFraction("SIH4"))
    conc_sih2.append(mult* reactor.moleFraction("SIH2"))
    conc_h2.append(mult* reactor.moleFraction("H2"))
    press.append(reactor.pressure())
    temp.append(reactor.temperature())
    dens.append(1.0e-3*gas.molarDensity())
    rates=gas.netRatesOfProgress()
    brates=[]
    for i in range(0, len(rates)):
        #Convert to mol/cm3s
        brates.append(rates[i]*1.0e3)
    rxnrates.append(brates)

def writeLine(ostream, items, delimiter):
    outstr = ""
    L = len(items)-1
    j = 0
    for i in items:
        outstr += "{0:4e}".format(i)
        if j!=L: outstr += delimiter
        j += 1
    outstr += "\n"
    ostream.write(outstr)

wkdir="."
ifile=wkdir+"/"+"chem.troe.xml"

gas = can.importPhase(ifile)
gas.set(\
    T = 1400, \
    P = 2500, \
    X = "SIH4:0.1, AR:0.90" \
    )

# Create the atmosphere
atm = can.importPhase(ifile)
atm.set(T=1400, P=101325, X="AR:1.0")
res = reac.Reservoir(atm)

# Create the batch reactor
reactor = reac.Reactor(gas)

# Create the wall
w1 = reac.Wall(reactor, res);
w1.setHeatTransferCoeff(1.0e10)
w1.setArea(1.0)

# Create the network and begin
network = reac.ReactorNet([reactor])
t = 0.0
tau = 100.0
steps = 2000
mesh = []
for a in range(0, steps):
    mesh.append((tau/steps)*a)

conc_sih4 = []
conc_sih2 = []
conc_h2   = []
press     = []
temp      = []
rxnrates  = []
dens      = []
getConcs()
for t in mesh[1:]:
    network.advance(t) 
    getConcs()

# Write the results
ostream = open(wkdir+"/"+"results.csv", "w")
ostream.write("time,sih4,sih2,h2,P,T,dens\n")
for t, a, b, c, p, T, d in zip(mesh, conc_sih4, conc_sih2, conc_h2, press, temp, dens):
    items = [t, a, b, c, p, T, d]
    writeLine(ostream, items, ",")
ostream.close()

ostream = open(wkdir+"/"+"results-rates.csv", "w")
ohead="time"
for i in rxnrates[0]:
    ohead+=",rxn"
ohead+="\n"
ostream.write(ohead)
for t, line in zip(mesh, rxnrates):
    writeLine(ostream, [t]+line, ",")
ostream.close()
