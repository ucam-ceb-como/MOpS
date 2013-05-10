# Networked PSR wrapper around Cantera
# Provides a convenient way of adding networked reactors and
# postprocessing the results. 
#
# Note that the 'Simulator' object can also be used with a normal
# Cantera ReactorNet object.
#
# William Menz (wjm34) 2013

import time
import csv
import sys

# Load cantera libraries
import Cantera as can
import Cantera.Reactor as reac

class Reactor:
    # Wrapper class for a cantera reactor. Provides simplified way of
    # specifying inflows and outflows.
    
    def __init__(self, mix, rt=1.0, **kwargs):
        # Initialise a Cantera reactor object
        self.reac = reac.Reactor(mix, **kwargs)
        # The residence time
        self.rt = rt
        
        # A list of references to in/outflows (reacs or reservoirs)
        self.ilist = []
        self.olist = []
        
        self.ifracs = []
        self.ofracs = []
        
        # A list of references to mass flow controllers corresponding
        # to the order in the lists above
        self.icontlist = []
        self.ocontlist = []
    
    def density(self): return self.reac.density()
    
    def moldensity(self): 
        return self.reac.pressure() / (8.314 * self.reac.temperature())
    
    def temperature(self): return self.reac.temperature()
    
    def pressure(self): return self.reac.pressure()
    
    def volume(self): return self.reac.volume()
    
    def contents(self): return self.reac.contents()
    
    def name(self): return self.reac.name()
    
    def mass_flow(self):
        return self.reac.density() * self.reac.volume() / self.rt
    
    def set_constv(self, res):
        reac.Wall(self.reac, res, K=0.0)
    
    def set_constp(self, res):
        reac.Wall(self.reac, res, K=1.0e11)
    
    def add_inflow(self, r, frac=1.0):
        # Adds an inflow (r should be a Res or Reactor)
        fr = r.density() * self.reac.volume() / self.rt
        cont = reac.MassFlowController(upstream=r, downstream=self.reac, mdot=fr)
        self.ilist.append(r)
        self.ifracs.append(frac)
        self.icontlist.append(cont)
    
    def add_outflow(self, r, frac=1.0):
        # Adds an outflow
        cont = reac.MassFlowController(upstream=self.reac, downstream=r, mdot=frac*self.mass_flow())
        self.olist.append(r)
        self.ofracs.append(frac)
        self.ocontlist.append(cont)
    
    def update_flows(self):
        # Updates the mass flow controllers
        
        # Inflows
        assert(len(self.icontlist) == len(self.ilist))
        for c, r, f in zip(self.icontlist, self.ilist, self.ifracs):
            fr = f * r.density() * self.reac.volume() / self.rt
            c.set(mdot = fr)
        
        # Outflows
        assert(len(self.ocontlist) == len(self.olist))
        for c, r, f in zip(self.ocontlist, self.olist, self.ofracs):
            c.set(mdot = f * self.mass_flow())

class Network:

    def __init__(self, rlist):
        
        self.rlist = rlist
        
        self.raw_rlist = []
        for r in self.rlist: self.raw_rlist.append(r.reac)
        self.net   = reac.ReactorNet(self.raw_rlist)
    
    def reactors(self):
        return self.net.reactors()
    
    def update_flowrates(self):
        # Updates the flowrates in and out of a reactor to ensure that
        # the mdot is correct
        
        for r in self.rlist: r.update_flows()
    
    def advance(self, t):
        self.update_flowrates()
        try:
            self.net.advance(t)
        except Exception:
            print("Failed solving reactor.")
            raise
            sys.exit(4)

class Simulator:
    
    def __init__(self, times):
        
        # The times over which to solve.
        self.times = times
        
        # List of species to track
        self.slist = []
        # List of reactor names for writing later
        self.reacs = []
        # Output data storage dictionary
        self.odata = {}
    
    def __get_output_headers(self, r):
        # Gets the headers for writing
        heads = ["Time(s)", "T(K)", "P(Pa)", "V(m3)", "rho(mol/cm3"]
        
        if self.slist:
            for s in self.slist: heads.append(str(s)+"(mol/cm3)")
        else:
            c = r.contents()
            for s in c.speciesNames(): heads.append(str(s)+"(mol/cm3)")
        
        return heads
    
    def __get_data_from_reactor(self, r):
        # Gets a list of data from the reactor.
        data = []
        
        # Multiplier to mol/cm3
        mult = r.pressure()/(1.0e6*8.315*r.temperature())
        
        # Add T, P, V,  rho
        data.append(r.temperature())
        data.append(r.pressure())
        data.append(r.volume())
        data.append(mult)
        
        c = r.contents()
        
        if self.slist:
            for s in self.slist: data.append(mult * c.moleFraction(s))
        else:
            for s in c.speciesNames():
                data.append(mult * c.moleFraction(s))
        
        return data
    
    def __store_data(self, t, net):
        # Stores the data from each reactor for writing later
        
        # Initialise the output data array if haven't done already
        if not self.odata:
            for r in net.reactors(): self.odata[r.name()] = []
        
        for r in net.reactors():
            self.odata[r.name()].append(
                [t] + self.__get_data_from_reactor(r))
            
    
    def set_track_list(self, slist):
        # Sets a list of species to track through time
        self.slist = slist
    
    def solve(self, net):
        # Solve a network and record its data
        self.reacs = []
        self.odata = {}
        
        # Store the reactor references for use later
        for r in net.reactors():
            self.reacs.append(r)
        self.__store_data(0.0, net)
        
        # Step through the network by advancing each time
        tstart = time.time()
        for t in self.times:
            t0 = time.time()
            net.advance(t)
            self.__store_data(t, net)
            
            print("Advance to {0:.2e}s, took {1:.2e}s".format(t, time.time()-t0))
        
        print("Simulation finished, taking {:.2f}s".format(time.time()-tstart))
    
    def write_output(self, ofile="output"):
        
        # Loop over reactors
        for r in self.reacs:
            print("Writing output data for " + r.name() + ".")
            
            ostr = open(ofile + "-" + r.name() + ".csv", 'wb')
            cwriter = csv.writer(ostr, delimiter=",")
            
            # Get the data for writing
            lines = []
            lines.append(self.__get_output_headers(r))
            
            for l in self.odata[r.name()]: lines.append(l)
            
            # Write the data
            for l in lines: cwriter.writerow(l)
            
            # Clost the file
            ostr.close()
