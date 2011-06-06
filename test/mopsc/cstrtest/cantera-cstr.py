# CSTR modelling the thermal decomposition of silane

from Cantera import *
from Cantera.Reactor import *
from Cantera.Phase import *

# KEY INPUT PARAMETERS
p_1 = 2500	# pressure, Pa 
t_1 = 1373	# temperature, K
v_1 = 1		# volume, m3
restime = 10	# residence time, s
time = 20    	# total time to run, s
dt = 0.01   	# step size, s
comp_i = 'AR:0.96, SIH4:0.04' #in mol fractions
gas_inp = importPhase('cantera-chem.xml', 'gas')
output = open('cantera-results.csv','w')

# Other initialisation...
flow_in = gas_inp.density() * v_1 / restime
Ntimesteps = int(time/dt + 1)   #Num timesteps + 1 for 0s
gas_inp.setMoleFractions(comp_i)
gas_inp.set(T = t_1, P = p_1)
res_inp = Reservoir(gas_inp)
exhaust = Reservoir(gas_inp)
bathgas = Reservoir(gas_inp)
r1 = Reactor(gas_inp, v_1)
m1 = MassFlowController(upstream = res_inp, downstream = r1, mdot = flow_in)
exhaust_valve = Valve(r1,exhaust)
exhaust_valve.setValveCoeff(1.0)
w1 = Wall(r1,bathgas)
w1.setHeatTransferCoeff(10000)

# Define file & arrays for output data
num_species = gas_inp.nSpecies()
nam_species = gas_inp.speciesNames()
ind_species = range(num_species)
mfr_species = range(Ntimesteps)
print 'There are', num_species, 'species:', nam_species
speciesandtemp = ['Time (s)'] + ['Temp (K)'] + ['Pressure (Pa)'] + ['Vol (m3)'] + ['Dens (mol/cm3)'] + nam_species 
out_length = len(speciesandtemp)
writeCSV(output, speciesandtemp)
for i in range(Ntimesteps):
    mfr_species[i] = range(0, out_length)

# Define environment gas
# gas_env = Air()
sim = ReactorNet([r1])
time = 0.0
tim = zeros(Ntimesteps,'d')

# Run simulation
for run in range(Ntimesteps):
    tim[run] = time
    mfr_species[run][0] = time
    mfr_species[run][1] = r1.temperature()
    mfr_species[run][2] = r1.pressure()
    mfr_species[run][3] = r1.volume()
    mol_density = r1.pressure() / (r1.temperature() * 8.314) # in mol/m3
    mol_density /= 1000000 # convert to cm3
    mfr_species[run][4] = mol_density
    for i in range(5, out_length):
            mfr_species[run][i] = mol_density * gas_inp.moleFraction(ind_species[i-5])
    writeCSV(output,   mfr_species[run] )
    time += dt
    sim.advance(time)

output.close()

print 'Done!'
