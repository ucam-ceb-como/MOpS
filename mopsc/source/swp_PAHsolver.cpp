/*
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2009 Markus Sander.

  File purpose:
    Implementation of the PAHSolver class declared in the
    swp_PAHsolver.h header file.

  Licence:
    This file is part of "mops".

    mops is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

  Contact:
    Dr Markus Kraft
    Dept of Chemical Engineering
    University of Cambridge
    New Museums Site
    Pembroke Street
    Cambridge
    CB2 3RA
    UK

    Email:       mk306@cam.ac.uk
    Website:     http://como.cheng.cam.ac.uk
*/

#include "swp_PAHsolver.h"
#include "swp_gas_profile.h"
#include "mops_timeinterval.h"
#include "mops_reactor.h"
#include "sweep.h"
#include "string_functions.h"
#include "csv_io.h"
#include "swp_PAH_primary.h"
#include <fstream>
#include <stdexcept>
#include <string>
#include "swp_PAH_trajectory.h"

using namespace Sweep;
using namespace Sweep::ActSites;
using namespace std;
using namespace Strings;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
PAHSolver::PAHSolver()
{
}

// Default destructor.
PAHSolver::~PAHSolver()
{
}


 
// Reads a flame gas-phase profile from a TAB formatted file.
void PAHSolver::LoadGasProfile(const std::string &file, Mops::Mechanism &mech)
{
    map<real,real> alpha_prof;

    // Clear the current gas-phase profile.
    m_gasprof.clear();

    // Open the file to read.
    ifstream fin(file.c_str(), ios::in);
    if (!fin.good()) {
        throw runtime_error("Unable to open gas profile input "
                            "file (Mops, Sweep::FlameSolver::LoadGasProfile).");
    }

    // Variables read from file.
    vector<string> subs;
    string delim = ",\t \r"; // Possible file delimiters (comma, tab and space).
    string line;

    // Get the first line (should be the header defining the columns).
    if (!getline(fin, line).eof()) {

        // Split the line to get the column headings.
        split(line, subs, delim);

        // Get important column indices (time, temperature and pressure).
        int tcol=-1, Tcol=-1, Pcol=-1, Acol = -1,Rcol=-1;
        tcol = findinlist(string("Time"), subs);
        Tcol = findinlist(string("T"), subs);
        Pcol = findinlist(string("P"), subs);
        Acol = findinlist(string("Alpha"), subs);
        Rcol = findinlist(string("wdotA4"), subs);

        // Check that the file contains required columns.
        if (tcol < 0) {
            fin.close();
            throw runtime_error("Gas-phase profile contains no Time "
                                "column (Mops, Sweep::FlameSolver::LoadGasProfile).");
        }
        if (Tcol < 0) {
            fin.close();
            throw runtime_error("Gas-phase profile contains no temperature "
                                "column (Mops, Sweep::FlameSolver::LoadGasProfile).");
        }
        if (Pcol < 0) {
            fin.close();
            throw runtime_error("Gas-phase profile contains no pressure "
                                "column (Mops, Sweep::FlameSolver::LoadGasProfile).");
        }


        // All other columns are chemical species.  Add them, and note
        // their columns.
        map<unsigned int,int> spcols;
        for (int i=0; (unsigned)i!=subs.size(); ++i) {
            if ((i!=tcol) && (i!=Tcol) && (i!=Pcol) && (i!=Acol) && (i!=Rcol)) {
                Sprog::Species *sp = mech.AddSpecies();
                sp->SetName(subs[i]);
                spcols[i] = mech.FindSpecies(subs[i]);
            }
        }

        // Now we can read the profile.
        while(!getline(fin, line).eof()) {
            // Set t=0 and create a new IdealGas object.
            real t = 0.0;
            real T = 0.0;
            real P = 0.0;
            real alpha = 0.0;
            real PAHRate = 0.0;
            GasPoint gpoint(mech.Species());
//            Sprog::Thermo::IdealGas gas(mech.Species());

            // Split the line by columns.
            split(line, subs, delim);

            // Loop over all the elements in the line and save them
            // to the correct gas-phase variable.
            for (int i=0; (unsigned)i!=subs.size(); ++i) {
                if (i==tcol) {
                    // This is the Time column.
                    t = cdble(subs[i]);
                    gpoint.Time = t;
                } else if (i==Tcol) {
                    // This is the temperature column.
                    T = cdble(subs[i]);
                } else if (i==Pcol) {
                    // This is the pressure column.
                    P = cdble(subs[i]);
                } else if (i==Acol) {
                    alpha = cdble(subs[i]);
                } else if (i==Rcol) {
                    PAHRate = cdble(subs[i]);
                } else {
                    // This is a gas-phase species column.
                    map<unsigned int,int>::iterator isp = spcols.find(i);
                    if (isp != spcols.end()) {
                        gpoint.Gas.RawData()[isp->second] = cdble(subs[i]);
                    }
                }
            }

            // Set up the gas-phase by setting temperature, pressure and
            // normalising the mixture fractions.
            // TODO:  This will give the wrong component densities
            //        unless all species are specified!
            gpoint.Gas.SetTemperature(T);
            gpoint.Gas.SetPressure(P*1.0e5);
            gpoint.Gas.Normalise();
            gpoint.Gas.SetPAHFormationRate(PAHRate*1E6);

            // Add the profile point.
            alpha_prof[t] = alpha;
            m_gasprof.push_back(gpoint);
        }


        // Close the input file.
        fin.close();

        // Sort the profile by time.
        SortGasProfile(m_gasprof);
    } else {
        // There was no data in the file.
        fin.close();
        throw runtime_error("Input file contains no data "
                            "(Mops, Sweep::FlameSolver::LoadGasProfile).");
    }
}


void PAHSolver::Solve(Mops::Reactor &r, real tstop, int nsteps, int niter, 
                        Mops::Solver::OutFnPtr out, void *data)
{

    real tsplit, dtg, dt, jrate;
    const Sweep::Mechanism &mech = r.Mech()->ParticleMech();
    fvector rates(mech.TermCount(), 0.0);
 

	// Global maximum time step.
    real t  = r.Time();
    dtg     = tstop - t;
    m_maxdt = dtg / 3.0;
    m_tstop = tstop;

    // Loop over time until we reach the stop time.

    while (t < tstop)
    {
        tsplit = tstop;

        //save the old reactors density
        double old_dens=r.Mixture()->MassDensity();
        // Calculate the chemical conditions.
        linInterpGas(t, m_gasprof, *r.Mixture());

        // Scale M0 according to gas-phase expansion.
        real m0 = r.Mixture()->ParticleCount()/r.Mixture()->SampleVolume();
        r.Mixture()->SetM0(r.Mixture()->MassDensity() *  m0 / old_dens);

        // Perform stochastic jump processes.
        while (t < tsplit) {
            // Calculate jump rates.
            jrate = mech.CalcJumpRateTerms(t, *r.Mixture(), Geometry::LocalGeometry1d(), rates);

            // Perform time step.
            dt = timeStep(t, *r.Mixture(), mech, rates, jrate);

            if (dt >= 0.0) {
                t += dt; t = min(t, tstop);
            } else {
                return;
            }
        }

        // Update the PAH molecules inside the particles
        linInterpGas(t, m_gasprof, *r.Mixture());
        mech.LPDA(t, *r.Mixture());

    }


    // Call the output function.
    if (out) 
    {
        out(nsteps, niter, r, *this, data);
    }
	
	
    r.SetTime(t);
	
    return;
}


void PAHSolver::linInterpGas(Sweep::real t, 
                               const GasProfile &gasphase, 
                               Sprog::Thermo::IdealGas &gas) const
{
    // Get the time point after the required time.
    GasProfile::const_iterator j = LocateGasPoint(gasphase, t); //gasphase.upper_bound(t);
    
    if (j == gasphase.begin()) {
        // This time is before the beginning of the profile.  Return
        // the first time point.
        gas = j->Gas;
    } else if (j == gasphase.end()) {
        // This time is after the profile.  Return the last time
        // point
        --j;
        gas = j->Gas;
    } else {
        // Get the time point before the required time.
        GasProfile::const_iterator i = j; --i;

        // Assign the conditions to this point.
        gas = i->Gas;
        
        // Calculate time interval between points i and j.
        real dt_pro = j->Time - i->Time;

        // Calculate time interval between point i and current time.
        real dt = t - i->Time;

        // Calculate the intermediate gas-phase mole fractions by linear
        // interpolation of the molar concentrations.
        real dens = 0.0;
        for (unsigned int k=0; k<gas.Species()->size(); ++k) {
            real dc = (j->Gas.MolarConc(k) - i->Gas.MolarConc(k)) * dt / dt_pro;
            gas.RawData()[k] = gas.MolarConc(k) + dc;
            dens += gas.RawData()[k];
        }
        gas.Normalise();

        // Now use linear interpolation to calculate the temperature.
        real dT = (j->Gas.Temperature() - i->Gas.Temperature()) * dt / dt_pro;
        gas.SetTemperature(gas.Temperature()+dT);

        // Now set the gas density, calculated using the values above.
        gas.SetDensity(dens);
    }
}
