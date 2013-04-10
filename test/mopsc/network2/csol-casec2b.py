# CASE C2b example
# i.e. Recycle network, adiabatic constant volume.

# chem.xml obtained by converting the chem.inp using ck2cti.

import math
import imp
psr = imp.load_source('', 'cwrapper.py')

import Cantera as can
import Cantera.Reactor as reac

##########################################################
# Set up mechanism and mixtures
# Note: copying phases doesn't seem to work, always importPhase.
imech = "chem.xml"
inert = can.importPhase(imech)
inert.set(T=1000, P=101325, X="AR:1.0")
inert2 = can.importPhase(imech)
inert2.set(T=1000, P=101325, X="AR:1.0")

inGas = can.importPhase(imech)
inGas.set(
    T=1000,
    P=101325,
    X="SIH4:0.1, AR:0.9"
    )
resGas = can.importPhase(imech)
resGas.set(T=1000, P=101325, X="AR:1.0")
exGas  = can.importPhase(imech)
exGas.set(T=1000, P=101325, X="AR:1.0")

##########################################################
# Create reactors

inflow = reac.Reservoir(inGas)
exhaust = reac.Reservoir(exGas)
res = reac.Reservoir(resGas)

r1 = psr.Reactor(inert, 1.0, energy="on")
r1.set_constv(res)
r1.add_inflow(inflow, 0.6)

r2 = psr.Reactor(inert, 1.0, energy="on")
r2.set_constv(res)
r1.add_outflow(r2.reac, 1.0)
r2.add_outflow(exhaust, 0.6)
r1.add_inflow(r2.reac, 0.4)

rlist = [r1, r2]
net = psr.Network(rlist)

##########################################################
# Create timesteps in a logarithmic style mesh
times = []
fTime = 3.0     # final time to end at
steps = 1000      # steps per stop interval

tstop = [fTime/100, fTime/10, fTime]
for i, ts in enumerate(tstop):
    if i == 0: tstart = 0.0
    else: tstart = tstop[i-1]
    dt = (ts - tstart)/steps
    for j in range(1, steps+1):
        times.append(tstart + j*dt)

##########################################################
# Load simulator and solve
sim = psr.Simulator(times)
#sim.set_track_list(["AR"])
sim.solve(net)
sim.write_output("cantera")
