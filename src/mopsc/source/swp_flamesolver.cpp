/*
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the FlameSolver class declared in the
    swp_flamesolver.h header file.

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

#include "swp_flamesolver.h"
#include "swp_gas_profile.h"
#include "mops_timeinterval.h"
#include "mops_reactor.h"
#include "sweep.h"
#include "string_functions.h"
#include "csv_io.h"
#include <fstream>
#include <stdexcept>
#include <string>
#include <iomanip>

using namespace Sweep;
using namespace Sweep::ActSites;
using namespace std;
using namespace Strings;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
FlameSolver::FlameSolver()
{
}

// Default destructor.
FlameSolver::~FlameSolver()
{
}

// PROFILE INPUT.

// Reads a flame gas-phase profile from a TAB formatted file.
// Species concentrations should be expressed in mole fractions with molar
// concentrations being inferred from the temperature and pressure.  This
// rescaling means that quantities proportional to the mole fractions, such
// as molar concentrations may also be used, but will still be rescaled to satisfy
// PV = nRT.  This rescaling means that in all cases the input data must cover
// all gas phase species, if species are omitted the molar concentrations of the
// remaining species will be overstated.
void FlameSolver::LoadGasProfile(const std::string &file, Mops::Mechanism &mech)
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
        int tcol=-1, Tcol=-1, Pcol=-1, Acol = -1, Rcol=-1;
        tcol = findinlist(string("Time"), subs);

        Tcol = findinlist(string("T"), subs);
        if(Tcol < 0)
            Tcol = findinlist(string("T[K]"), subs);

        Pcol = findinlist(string("P"), subs);
        Acol = findinlist(string("Alpha"), subs);
        Rcol = findinlist(string("wdotA4"), subs);

        // Columns to ignore, but which are useful to have in files for brush compatibility
        int Xcol = findinlist(string("X[cm]"), subs);
        int Dcol = findinlist(string("RHO[g/cm3]"), subs);
        int Vcol = findinlist(string("V[cm/s]"), subs);
        int Gcol = findinlist(string("GradT"), subs);


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

        if (Acol < 0 && Rcol <0) {
            fin.close();
            throw runtime_error("Gas-phase profile contains no alpha "
                                "column (Mops, Sweep::FlameSolver::LoadGasProfile).");
        }

        // All other columns are chemical species.  Add them, and note
        // their columns.
        map<unsigned int,int> spcols;
        for (int i=0; (unsigned)i!=subs.size(); ++i) {
            if ((i!=tcol) && (i!=Tcol) && (i!=Pcol) && (i!=Acol) && (i!=Rcol) &&
                (i!=Xcol) && (i!=Dcol) && (i!=Vcol) && (i!=Gcol)) {
                // Try to find this species in the mechanism
                const int speciesMechIndex = mech.FindSpecies(subs[i]);

                if(speciesMechIndex < 0) {
                    std::ostringstream msg("Failed to find species ");
                    msg << subs[i] << " in mechanism (Mops, Sweep::FlameSolver::LoadGasProfile).";
                    throw std::runtime_error(msg.str());
                }
                // Found species
                spcols[i] = speciesMechIndex;
            }
        }

        // riap For logging input in way that should match PSDF_input.dat files
        /*unsigned int logIndices[7];
        logIndices[0] = mech.FindSpecies("C2H2");
        logIndices[1] = mech.FindSpecies("H2");
        logIndices[2] = mech.FindSpecies("H");
        logIndices[3] = mech.FindSpecies("O2");
        logIndices[4] = mech.FindSpecies("OH");
        logIndices[5] = mech.FindSpecies("H2O");
        logIndices[6] = mech.FindSpecies("A4");
        unsigned int logIndices[14];
        logIndices[0]  = mech.FindSpecies("H2");
        logIndices[1]  = mech.FindSpecies("H");
        logIndices[2]  = mech.FindSpecies("O");
        logIndices[3]  = mech.FindSpecies("O2");
        logIndices[4]  = mech.FindSpecies("OH");
        logIndices[5]  = mech.FindSpecies("H2O");
        logIndices[6]  = mech.FindSpecies("CH4");
        logIndices[7]  = mech.FindSpecies("CO");
        logIndices[8]  = mech.FindSpecies("CO2");
        logIndices[9]  = mech.FindSpecies("C2H2");
        logIndices[10] = mech.FindSpecies("C2H4");
        logIndices[11] = mech.FindSpecies("AR");
        logIndices[12] = mech.FindSpecies("N2");
        logIndices[13] = mech.FindSpecies("A4");*/

        // Now we can read the profile.
        while(!getline(fin, line).eof()) {
            // Set t=0 and create a new IdealGas object.
            real t = 0.0;
            real T = 0.0;
            real P = 0.0;
            real alpha = 0.0;
            real PAHRate = 0.0;
            GasPoint gpoint(mech.Species());

            // Split the line by columns.
            split(line, subs, delim);

            // Check that the mole fractions sum to 1
            real checkSum = 0.0;

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
                        const real frac = cdble(subs[i]);
                        assert(isp->second >= 0);
                        gpoint.Gas.RawData()[isp->second] = frac;
                        checkSum += frac;
                    }
                }
            }

            if((checkSum < 0.997) || checkSum > 1.003) {
                std::ostringstream msg;
                msg << "Mole fractions sum to " << checkSum
                    << ", but should sum to 1.000 (FlameSolver::LoadGasProfile)";
                throw std::runtime_error(msg.str());
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

            // Output in PSDF_input.dat format
            //std::cout << t << '\t' << std::scientific << std::setprecision(6) << T << '\t';
            //for(unsigned int j = 0; j != 14; ++j) {
            //for(unsigned int j = 0; j != 5; ++j) {
            //    std::cout << gpoint.Gas.MolarConc(logIndices[j]) * 1e-6 << '\t';
            //}
            //std::cout << alpha << '\t' << gpoint.Gas.Pressure() << '\t'
            //          << gpoint.Gas.Density() * 1e-6 << '\n';
        }

        // Set up ABF model to use alpha profile.
        ABFModel::Instance().Initialise(mech.ParticleMech());
        ABFModel::Instance().SetAlphaProfile(alpha_prof);

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


// SOLUTION AND POST-PROCESSING.

// Performs stochastic stepping algorithm up to specified stop time using
// the given mechanism to define the stochastic processes.  Updates given
// system accordingly.  On error returns <0, otherwise returns 0.  In this
// flavour the gas-phase chemistry is interpolated from a vector of
// IdealGas objects rather than being taken from the given system object.
// However, the particles in the system object are updated accordingly.
void FlameSolver::Solve(Mops::Reactor &r, real tstop, int nsteps, int niter,
                        int (*rand_int)(int, int), Sweep::real (*rand_u01)(),
                        Mops::Solver::OutFnPtr out, void *data)
{
//    int err = 0;
    real tsplit, dtg, dt, jrate;
    const Sweep::Mechanism &mech = r.Mech()->ParticleMech();
    fvector rates(mech.TermCount(), 0.0);
    
    // Save the initial chemical conditions in sys so that we
    // can restore them at the end of the run.
    const Sprog::Thermo::IdealGas chem = *r.Mixture();

    // Store if chemical conditions are fixed at present, because we
    // shall set them to be fixed during this run, to be restored afterwards.
    bool fixedchem = r.Mixture()->FixedChem();
    r.Mixture()->SetFixedChem();

    // Global maximum time step.
    real t  = r.Time();
    dtg     = tstop - t;
    setMaxTimeStep(dtg / 3.0);

    // Set the chemical conditions.
    linInterpGas(t, m_gasprof, *r.Mixture());

    // Loop over time until we reach the stop time.
    while (t < tstop)
    {

        //save the old gas phase mass density
        double old_dens = r.Mixture()->MassDensity();

        // Update the chemical conditions.
        linInterpGas(t, m_gasprof, *r.Mixture());

        // Scale particle M0 according to gas-phase expansion.
        r.Mixture()->AdjustSampleVolume(old_dens / r.Mixture()->MassDensity());

        // Get the process jump rates (and the total rate).
        jrate = mech.CalcJumpRateTerms(t, *r.Mixture(), Geometry::LocalGeometry1d(), rates);

        // Calculate the splitting end time.
        tsplit = calcSplitTime(t, tstop, jrate, r.Mixture()->ParticleCount(), dtg);

        // Perform stochastic jump processes.
        while (t < tsplit) {
            // Calculate jump rates.
            jrate = mech.CalcJumpRateTerms(t, *r.Mixture(), Geometry::LocalGeometry1d(), rates);

            // Perform time step.
            dt = timeStep(t, tsplit, *r.Mixture(), mech, rates, jrate, rand_int, rand_u01);
            if (dt >= 0.0) {
                t += dt; 

                // The following line is suspicious, it seems to be a mathematically
                // inadmissible solution to the problem of stochastic time steps
                // going past the end of the period under consideration.  The
                // correct solution is to simulate the random time step, but
                // only to perform the event if it falls within the interval
                // being considered, otherwise one should just move the time 
                // counter to the end of the period. (riap 25/05/2009)
                t = std::min(t, tstop);

            } else {
                return;
            }
        }

        // Perform Linear Process Deferment Algorithm to
        // update all deferred processes.
        // Perhaps better to use t - 0.5 * dtg not just t
        mech.LPDA(t, *r.Mixture(), rand_int, rand_u01);

        r.SetTime(t);
    }
    
    // Restore initial chemical conditions to sys.
    r.Mixture()->SetFixedChem(fixedchem);

    out(nsteps, niter, r, *this, data);

    return;
}

// HELPER ROUTINES.

void FlameSolver::linInterpGas(Sweep::real t, 
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
