#!/usr/bin/python

# PYTHON IMPLEMENTED ODE SOLUTION OF SIMPLE POPULATION EQUATIONS
# BASED ON THE "DISCRETE MODEL" OF FRENKLACH & HARRIS 1987
# (C) WILLIAM MENZ WJM34 2013

# PSR OR BATCH REACTOR WITH COAGULATION, SURFACE REACTION, OXIDATION
# OR INCEPTION. REQUIRES pinput.py FILE FOR SETTING PARAMETERS.

# Imports
import random
import math
import csv
import numpy as np
import time
from scipy import integrate
"""
# Import matplotlib for optional plotting interface
try:
    import matplotlib.pyplot as plt
except:
    print("Could not find matplotlib dependencies.")

# Call to mview for PSD calculation
# Obtain at CoMoCommon/Codes/Utilities/mview-py-git
import imp
try:
    mstats = imp.load_source('', '/home/wjm34/programs/mview-py-git/structure/stats.py')
except:
    print("Could not find mview stats source. PSD calculation disabled.")
"""

######################################################################
# IMPORT INPUTS AS GLOBAL VARIABLES
######################################################################
from pinput import *

######################################################################
# OBJECT AND FUNCTION DECLARATIONS
######################################################################
def dN_dt(N, t=0.0):
    # The system of equations function to be passed to scipy.integrate
    # Relies on declaration of (icky) global variables and functions.
    
    # Initialise output storage
    derivs = []
    
    # Loop over each of the discrete particle sizes
    # Here i = 0 corresponds to N1 (type 1)
    for i in range(0, NMAX):
        # Initialise value of derivative
        val = 0.0
        # Add inception for type 1 (i=9) particles
        if i == 0: val += I(t)
        
        # Add sr term
        s1 = -N[i]
        if i > 0: s1 += N[i-1]
        val += S(t) * s1
        
        # Add oxidation term
        o1 = -N[i]
        if i != NMAX - 1: o1 += N[i+1]
        val += OX(t) * o1
        
        # Do the first coagulation term (formation of type due to coag)
        c1 = 0.0
        for j in range(0, i):
            c1 += N[i-j-1] * N[j]
        
        # Do the second coag term
        c2 = 0.0
        for j in range(0, NMAX):
            c2 += N[i] * N[j]
        
        # Add the coagulation contribution
        val += K(t) * (0.5 * c1 - c2)
        
        # Add the transport term for a PSR
        if RESTIME: val += (NIN[i] - N[i]) / RESTIME
        
        derivs.append(val)
    
    return np.array(derivs)

def initialise_array(iarray, nmax):
    # Ensures the array is the same length as the number of equations.
    if nmax > len(iarray):
        print("Warning: More equations than initial conditions.")
        o = np.hstack((np.array(iarray), np.zeros(nmax - len(iarray))))
    elif nmax < len(iarray):
        print("Warning: Fewer equations than initial conditions.")
        o = np.array(iarray[:nmax])
    else:
        # Initial conditions as an array
        o = np.array(iarray)
    assert len(o) == nmax
    return o

def solve(ninit):
    # Generate some time mesh
    mesh = np.linspace(0.0, TAU, STEPS)
    n0 = initialise_array(ninit, NMAX)
    
    # Solve the system!
    t1 = time.time()
    N, info = integrate.odeint(dN_dt, n0, mesh, full_output=True)
    print(info["message"])
    print("Took " + str(time.time()-t1) + "s.")
    
    return N, mesh

"""
def plot_results(times, cols, labels=[]):
    # Plot the results from a integrated calculation.
    #print cols
    #print times
    
    plt.figure()
    for i, c in enumerate(cols):
        try: l = labels[i]
        except: l = "?"
        plt.plot(times, c, label=l)
    plt.legend(loc=0)
    plt.xlabel("Time, s")
    plt.ylabel("Population, ?")
    plt.show()
"""

def pdf_lognorm(val, gmean, gstdev):
    # Return the PDF for a lognormal distribution
    # http://en.wikipedia.org/wiki/Log-normal_distribution
    mu = math.log(gmean)
    sig = math.log(gstdev)
    ans = 1.0/(val*math.sqrt(2*math.pi*sig*sig))
    ans *= math.exp(- pow(math.log(val)-mu, 2.0)/(2.0*sig*sig))
    return ans

def generate_distribution(npts, mean, stdev, pdfn):
    # Generate an initial starting condition for n points with a given
    # mean and standard deviation. Use the PDF to evaluate the 
    # probability at the requested point.
    
    vals = []
    for i in range(npts):
        vals.append(pdfn(i+1, mean, stdev))
    
    # Normalise the distribution.
    vsum = sum(vals)
    for i, v in enumerate(vals): vals[i] /= vsum
    
    return vals

def write_mops_input(distrib, res=1, fname="input.xml"):
    # Shows an input for a MOPS.inx file to ensure the same starting
    # distribution is used (all hardcoded).
    cid = "silicon"
    stream = open(fname, "w")
    for i, d in enumerate(distrib):
        fstr = ""
        fstr = '<particle count="' + str(res) + '">\n'
        fstr += '  <weight>' + str(d) + '</weight>\n'
        fstr += '  <component id="' + cid + '" dx="' + str(i+1) + '"/>\n'
        fstr += '</particle>\n'
        stream.write(fstr)
    stream.close()

def generate_reference(times, num):
    # Using the formaula in Patterson et al 2011 (coagulation an mono-
    # disperse initial conditions), generate the analytical solutions.
    cols = np.zeros((num, len(times)))
    # Loop over times
    for t in range(len(times)):
        tt = times[t]
        # Loop over bin index
        for i in range(num):
            # note just i instead of i-1 due to zero indexing
            cols[i][t] = pow(1+0.5*tt,-2.0) * pow(tt/(tt+2.0), i)
    
    return cols

def get_stats(cols, fn=lambda x: x*1.0):
    # From the discrete trajectories in each column, get some basic
    # stats through time. Use fn to scale the size descriptor.
    m0 = np.zeros(len(cols[0])).T
    m1 = np.zeros(len(cols[0])).T
    m2 = np.zeros(len(cols[0])).T
    m3 = np.zeros(len(cols[0])).T
    avgs = np.zeros(len(cols[0])).T
    
    # Loop over column data
    for i, c in enumerate(cols.T):
        avgs[i] = get_avg(c)
        # Loop over size indicies in 
        for sz, conc in enumerate(c):
            m0[i] += conc
            v = fn(sz+1)
            m1[i] += (sz+1) * conc
            # Only convert M2 and M3 to be consistent with MOPS
            # (outputs M2 and M3 as mass instead of volume moments)
            m2[i] += v * v * conc
            m3[i] += v * v * v * conc
    return m0, avgs, m1, m2, m3

def get_avg(num_concs):
    # From a state of the system (given by an array of the number concs)
    # get the average of that quantity
    
    # Loop over num concs, only consider entries > 0.
    total = 0.0     # Sum of the 'weights'
    num   = 0       # Sum of (particle size index * weight)
    for i, c in enumerate(num_concs):
        if c > 0.0:
            total += c
            num += (i + 1) * c
    if total > 0: return num/total
    else: return 0.0

"""
def write_psd(state, fn):
    # From a state of the system (given by an array of the number concs)
    # get a number PSD. Use the fn argument to convert the size values
    # to another format (e.g. volume -> diameter)
    
    # Quit if stats weren't found.
    if globals().get(mstats, False):
        print("Failed finding mstats.")
        return
    
    vals = []
    wts = []
    for i, s in enumerate(state):
        vals.append(fn(i+1))
        wts.append(s)
    
    s = mstats.LDensity("", wts, vals, "ODE Sol")
    #s.setBounds(1.0e8, 1.4e10)
    print("Got bandwidth of " + str(s.stats["bandwidth"]))
    s.generate()
    
    out = []
    for m, v in zip(s.stats["mesh"], s.stats["pdf"]):
        out.append([m, v])
    print("Writing PSD.")
    write_csv(out, "psd.csv", ["Mesh Value", "Kernel Density"])
"""

def get_header(num):
    header = []
    for i in range(num):
        header.append("N"+str(i+1))
    return header

def write_csv(adata, fname, header=[]):
    # Writes a csv from a np.array
    cwriter = csv.writer(open(fname, "wb"), delimiter=",")
    if header: cwriter.writerow(header)
    for row in adata:
        cwriter.writerow(row)

######################################################################
# MAIN PROGRAM BODY
######################################################################
if __name__ == "__main__":
    # Write the initial conditions for MOPS
    if GMEAN and GSTDEV:
        NINIT = generate_distribution(NMAX, GMEAN, GSTDEV, pdf_lognorm)
        write_mops_input(NINIT, RES)
    
    # Check the NIN array if a PSR reactor.
    if RESTIME: NIN = initialise_array(NIN, NMAX)
    
    # Solve the system.
    res, t = solve(NINIT)
    
    # View the solutions to all
    #plot_results(t, res.T, get_header(NMAX))
    #write_csv(np.vstack((t, res.T)).T, OUT, ["Time"]+get_header(NMAX))
    
    # Look at analytical solution (coag only with monodisperse parts)
    #ref = generate_reference(t, NMAX)
    #plot_results(t, ref, get_header(NMAX))
    
    # Calculate some 'moments', using specified conversion function.
    sums, avgs, m1, m2, m3 = get_stats(res.T, lambda x: x*1000)
    #plot_results(t, np.array([sums]), ["Sum"])
    #plot_results(t, np.array([avgs]), ["Avg"])
    write_csv(
        np.vstack((t, avgs.T, sums.T, m1.T, m2.T, m3.T)).T,
        OUT,
        ["Time","Avg","M0","M1(Vol)","M2(Mass)","M3(Mass)"])
    
    # Calculate a PSD
    #to_diameter = lambda x: 1.0e9 * pow(6.0 * x / math.pi, 1.0/3.0)
    #write_psd(res[len(res)-1], to_diameter)
