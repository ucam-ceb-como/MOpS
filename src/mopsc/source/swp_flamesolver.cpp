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
using namespace std;
using namespace Strings;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
FlameSolver::FlameSolver()
{
	m_stagnation = false;
	m_endconditions = false;
}

//! Copy constructor
FlameSolver::FlameSolver(const FlameSolver &sol)
: ParticleSolver(sol),
  Sweep::Solver(sol),
  m_gas_prof(sol.m_gas_prof) {}

//! Clone the object
FlameSolver *const FlameSolver::Clone() const {
    return new FlameSolver(*this);
}

// Default destructor.
FlameSolver::~FlameSolver()
{
}

// Set the end conditions flag
void FlameSolver::SetEndConditions(bool val)
{
	m_endconditions = val;
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
    // Clear the current gas-phase profile.
    m_gas_prof.clear();
	
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
        if(tcol < 0)
            tcol = findinlist(string("Time[s]"),subs);

        Tcol = findinlist(string("T"), subs);
        if(Tcol < 0)
            Tcol = findinlist(string("T[K]"), subs);

        Pcol = findinlist(string("P"), subs);
        Acol = findinlist(string("Alpha"), subs);
        Rcol = findinlist(string("wdotA4"), subs);

		//! Columns necessary for post-process with stagnation flame transport correction
		int ucol=-1,vcol=-1,Diffcol= -1;
		// Convective velocity
		ucol = findinlist(string("ConvectiveVelocity[m/s]"), subs);
		// Thermophoretic velocity
		vcol = findinlist(string("ThermophoreticVelocity[m/s]"), subs);
		// Diffusion term
		Diffcol = findinlist(string("DiffusionTerm"), subs);
		// if velocity columns are supplied than turn on stagnation flame correction
		if (ucol > 0 && vcol > 0){
			m_stagnation = true;
			cout << "Stagnation flame correction turned on. \n";
			if(Diffcol < 0) cout << "Diffusion correction not supplied. \n";
		}

        // Columns to ignore, but which are useful to have in files for brush compatibility
        int Xcol = findinlist(string("X[cm]"), subs);
        int Dcol = findinlist(string("RHO[g/cm3]"), subs);
        int Vcol = findinlist(string("V[cm/s]"), subs);

        int Gcol = findinlist(string("GradT"), subs);
        if(Gcol < 0)
            Gcol = findinlist(string("GradT[K/cm]"), subs);


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
                (i!=Xcol) && (i!=Dcol) && (i!=Vcol) && (i!=Gcol) &&
				(i!=ucol) && (i!=vcol) && (i!=Diffcol)) {
                // Try to find this species in the mechanism
                const int speciesMechIndex = mech.GasMech().FindSpecies(subs[i]);

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
            double t = 0.0;
            double T = 0.0;
            double P = 0.0;
            double alpha = 0.0;
            double PAHRate = 0.0;
			double u_conv = 0.0;
			double v_thermo = 0.0;
			double diffusion_term = 0.0;
            GasPoint gpoint(mech.GasMech().Species());

            // Split the line by columns.
            split(line, subs, delim);

            // Check that the mole fractions sum to 1
            double checkSum = 0.0;

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
				} else if (i==ucol) {
					//this is the convective velocity column
					u_conv = cdble(subs[i]);
				} else if (i==vcol) {
					v_thermo = cdble(subs[i]);
				} else if (i==Diffcol) {
					diffusion_term = cdble(subs[i]);
                } else {
                    // This is a gas-phase species column.
                    map<unsigned int,int>::iterator isp = spcols.find(i);
                    if (isp != spcols.end()) {
                        const double frac = cdble(subs[i]);
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
            gpoint.Gas.SetPressure(P*1.0e5);//also set the molar density of gas mixture
            gpoint.Gas.Normalise();

			//! If using the sample volume correction (for a stagnation flame) 
			//! then set the convective and  thermophoretic velocities, and diffusion term
			if(m_stagnation == true){
				gpoint.Gas.SetConvectiveVelocity(u_conv);
				gpoint.Gas.SetThermophoreticVelocity(v_thermo);
				if(Diffcol > 0) {
					gpoint.Gas.SetDiffusionTerm(diffusion_term);
				}else{
					// If diffusion term not supplied then set to 0
					gpoint.Gas.SetDiffusionTerm(0.0);
				}
			}else{
				// set terms to 0
				gpoint.Gas.SetConvectiveVelocity(0.0);
				gpoint.Gas.SetThermophoreticVelocity(0.0);
				gpoint.Gas.SetDiffusionTerm(0.0);
			}

            //! If postprocessing based on the molar rate of production by
            //! chemical reaction of the inception species per unit volume
            //! wdotA4 (mol/cm3/s), PAHRate determines the rate at which the
            //! inception species is inserted into the ensemble. Otherwise it
            //! is not required; the number of stochastic particles made up of
            //! a single primary corresponding to the inception species is
            //! adjusted to match the inception species concentration in
            //! FlameSolver::Solve.
            if (mech.ParticleMech().Postprocessing() == ParticleModel::wdotA4) {
                gpoint.Gas.SetPAHFormationRate(PAHRate*1E6);    //!< Convert from mol/(cm3*s) to mol/(m3*s).
            } else {
                gpoint.Gas.SetPAHFormationRate(PAHRate*0);      //!< Explicitly set to 0 in case the wdotA4 column in the gasphase.inp file is non-zero.
            }

            gpoint.Gas.SetAlpha(alpha);

            // Add the profile point.
            m_gas_prof.push_back(gpoint);

            // Output in PSDF_input.dat format
            //std::cout << t << '\t' << std::scientific << std::setprecision(6) << T << '\t';
            //for(unsigned int j = 0; j != 14; ++j) {
            //for(unsigned int j = 0; j != 5; ++j) {
            //    std::cout << gpoint.Gas.MolarConc(logIndices[j]) * 1e-6 << '\t';
            //}
            //std::cout << alpha << '\t' << gpoint.Gas.Pressure() << '\t'
            //          << gpoint.Gas.Density() * 1e-6 << '\n';
        }

        // Close the input file.
        fin.close();

        // Sort the profile by time.
        SortGasProfile(m_gas_prof);

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
void FlameSolver::Solve(Mops::Reactor &r, double tstop, int nsteps, int niter,
                        rng_type &rng, Mops::Solver::OutFnPtr out, void *data)
{
    double tsplit, dtg, jrate;

    // construct kmcsimulater and initialize gasphase info for kmcsimulater 
    if (r.Mixture()->Particles().Simulator()==NULL)
    {
        r.Mixture()->Particles().SetSimulator(*(Gasphase()));
    }

    const Sweep::Mechanism &mech = r.Mech()->ParticleMech();
    fvector rates(mech.TermCount(), 0.0);

    // Save the initial chemical conditions in sys so that we
    // can restore them at the end of the run.
    const Sprog::Thermo::IdealGas chem = r.Mixture()->GasPhase();

    // Store if chemical conditions are fixed at present, because we
    // shall set them to be fixed during this run, to be restored afterwards.
    bool fixedchem = r.Mixture()->FixedChem();
    r.Mixture()->SetFixedChem();

    // Global maximum time step.
    double t  = r.Time();
    dtg     = tstop - t;

	//declare variables for sample volume adjustment
	double old_dens(0.0);
	//! sample volume adjustment
	if(m_stagnation == false){
    //! If the initial composition was not specified, linearly interpolate the
    //! gas-phrase profile to obtain properties at the initial time step which
    //! may not necessarily be zero.
    if (!(r.Mixture()->GasPhase().MassDensity() >= 0))
        linInterpGas(t, r.Mixture()->GasPhase());

    //! Save density from previous time step for sample volume adjustment.
		old_dens = r.Mixture()->GasPhase().MassDensity();

    //! Update the chemical conditions.
    linInterpGas(t, r.Mixture()->GasPhase());

    //! Adjust sample volume using the change in the density. Note that the
    //! code has to go through the loop below twice for the sample volume to be
    //! adjusted. However, it was found that the loop is only performed once;
    //! therefore, the sample volume adjustment has to be performed here.
		 r.Mixture()->AdjustSampleVolume(old_dens / r.Mixture()->GasPhase().MassDensity() );
	}

    // Loop over time until we reach the stop time.
    while (t < tstop)
    {
		//if (t == 0.0){
		//	if (r.Mixture()->ParticleModel()->Components(0)->WeightedPAHs()){
		//		Particle *sp = NULL;
		//		sp = mech.CreateParticle(t);
		//		sp->setStatisticalWeight(r.Mixture()->Particles().Capacity() / 4.0);
		//		sp->UpdateCache();
		//		(r.Mixture()->Particles()).Add(*sp, rng);
		//	}
		//	else{
		//		for (int iii = 0; iii < r.Mixture()->Particles().Capacity()/4.0; iii++){
		//			Particle *sp = NULL;
		//			sp = mech.CreateParticle(t);
		//			//sp->setStatisticalWeight(1.0);
		//			sp->UpdateCache();
		//			(r.Mixture()->Particles()).Add(*sp, rng);
		//		}
		//	}

		//}

        //save the old gas phase mass density
		old_dens = r.Mixture()->GasPhase().MassDensity();

        // Update the chemical conditions to time t.
        const double gasTimeStep = linInterpGas(t, r.Mixture()->GasPhase());

		//save conditions for stagnation flame correction
		double u_old = r.Mixture()->GasPhase().GetConvectiveVelocity();
		double v_old = r.Mixture()->GasPhase().GetThermophoreticVelocity();
		double rho_old = r.Mixture()->GasPhase().MassDensity();
		double diffusion_term = r.Mixture()->GasPhase().GetDiffusionTerm();
		double t_old = t;
		
		//!sample volume adjustment
		if(m_stagnation == false){
        // Scale particle M0 according to gas-phase expansion.
        // (considering mass const, V'smpvol*massdens' = Vsmpvol*massdens)
        r.Mixture()->AdjustSampleVolume(old_dens / r.Mixture()->GasPhase().MassDensity());
		}

		//! impose from the end of the time step
		if (m_endconditions == true){
			//update gas-phase conditions to right end point
			linInterpGas(tstop, r.Mixture()->GasPhase());
			//save new conditions
			double u_new = r.Mixture()->GasPhase().GetConvectiveVelocity();
			double v_new = r.Mixture()->GasPhase().GetThermophoreticVelocity();
			double rho_new = r.Mixture()->GasPhase().MassDensity();
			//apply sample volume adjustment now
			if (m_stagnation == true){
				//calculate correction
				double scale_factor = 1.0;
				if (t > 0.0) {
					scale_factor = 1.0 + (v_new - v_old - u_new*(rho_new - rho_old) / rho_new) / (u_new + v_new); //no diffusion term
				}
				//adjust sample volume
				r.Mixture()->AdjustSampleVolume(scale_factor);
			}
		}

        //! Tried and tested only for the PAH-PP/KMC-ARS model, binary tree and
        //! the spherical particle model. Only relevant if postprocessing based
        //! on the inception species concentration.
        if ((mech.AggModel() == AggModels::PAH_KMC_ID || mech.AggModel() == AggModels::BinTree_ID || mech.AggModel() == AggModels::Spherical_ID) && 
            r.Mixture()->ParticleModel()->Postprocessing() == ParticleModel::XA4) {
            int index;
            double numCarbons;

            switch (mech.InceptedPAH()){
                case ParticleModel::A1:
                    index=r.Mech()->GasMech().FindSpecies("A1");
                    numCarbons = 6;
                    break;
                case ParticleModel::A2:
                    index=r.Mech()->GasMech().FindSpecies("A2");
                    numCarbons = 10;
                    break;
                case ParticleModel::A4:
                    index=r.Mech()->GasMech().FindSpecies("A4");
                    numCarbons = 16;
                    break;
                case ParticleModel::A5:
                    index=r.Mech()->GasMech().FindSpecies("A5");
                    numCarbons = 20;
                    break;
                default:
                    throw std::runtime_error("no information about the incepted PAH is available, only A1 A2, A4 and A5 are supported now (Sweep::FlameSolver::Solve())");
            }
            if (mech.Inceptions(0) != NULL) {
                if (mech.AggModel() == AggModels::PAH_KMC_ID || mech.AggModel() == AggModels::BinTree_ID || mech.AggModel() == AggModels::Spherical_ID &&
                    mech.Inceptions(0)->ParticleComp(0) == numCarbons) {
					// calculate the amount of stochastic pyrene particles in the ensemble
					unsigned int Pamount=r.Mixture()->NumOfStartingSpecies(index);
					// if Pmount exceeds the capacity of the ensemble at the begining of the simulation,
					// the process should be terminated since further running is meaningless.
					if (t < 1.0e-20 && Pamount >= r.Mixture()->Particles().Capacity())
						throw std::runtime_error("increase the M0 in mops.inx please, current choice is too small (Sweep::FlameSolver::Solve)");
							mech.MassTransfer(Pamount, t, *r.Mixture(), rng, Geometry::LocalGeometry1d());
				}
            }
        }

        // Get the process jump rates (and the total rate).
        jrate = mech.CalcJumpRateTerms(t, *r.Mixture(), Geometry::LocalGeometry1d(), rates);

        // Calculate the splitting end time.
		if (m_endconditions == true){
			//enforce single split
			tsplit = tstop;
		} else{
			tsplit = calcSplitTime(t, std::min(t + std::min(dtg, gasTimeStep), tstop), jrate, r.Mixture()->ParticleWeightSum());
		}

        //std::cout << "At time " << t << " split time is " << tsplit << ", spacing of gas data is " << gasTimeStep << '\n';

        // Perform stochastic jump processes.
        while (t < tsplit) {
            // Calculate jump rates.
            jrate = mech.CalcJumpRateTerms(t, *r.Mixture(), Geometry::LocalGeometry1d(), rates);

            // Perform time step.
            timeStep(t, std::min(t + dtg / 3.0, tsplit), *r.Mixture(), Geometry::LocalGeometry1d(),
                     mech, rates, jrate, rng);

			if (r.Mixture()->ParticleCount() < r.Mixture()->Particles().DoubleLimit() && 
				r.Mixture()->Particles().IsDoublingOn() && 
				r.Mixture()->ParticleModel()->Components(0)->WeightedPAHs()){
				break;
        }

        }
        // Perform Linear Process Deferment Algorithm to
        // update all deferred processes.
        // Perhaps better to use t - 0.5 * dtg not just t
        mech.LPDA(t, *r.Mixture(), rng);

        r.SetTime(t);

		//! sample volume adjustment with stagnation flame correction
		if(m_stagnation == true && m_endconditions == false){
			//interpolate gas-phase conditions to t
			linInterpGas(t, r.Mixture()->GasPhase());
			double u_new = r.Mixture()->GasPhase().GetConvectiveVelocity();
			double v_new = r.Mixture()->GasPhase().GetThermophoreticVelocity();
			double rho_new = r.Mixture()->GasPhase().MassDensity();
			//calculate correction
			double scale_factor = 1.0;
			if( t > 0.0) {
				scale_factor = 1.0 + (v_new-v_old - u_old*(rho_new-rho_old)/rho_old)/(u_old+v_old) + diffusion_term * (t-t_old);
			}
			//adjust sample volume
			r.Mixture()->AdjustSampleVolume(scale_factor);
		}
    
    }

    // Restore initial chemical conditions to sys.
    r.Mixture()->SetFixedChem(fixedchem);

    out(nsteps, niter, r, *this, data);

    return;
}

/*!
 * @param[in]   t       Time at which to interpolate gas profile
 *
 * @param[out]  gas     Mixture into which to insert the newly interpolated gas properties
 *
 * @return      An estimate of the time over which the data remains approximately constant,
 *              that is for which it should be reasonable to use the gas object.
 */

double FlameSolver::linInterpGas(double t,
                               Sprog::Thermo::IdealGas &gas) const
{
    // Get the time point after the required time.
    GasProfile::const_iterator j = Sweep::LocateGasPoint(m_gas_prof, t); //gasphase.upper_bound(t);
    
    if (j == m_gas_prof.begin()) {
        // This time is before the beginning of the profile.  Return
        // the first time point.
        gas = j->Gas;
    } else if (j == m_gas_prof.end()) {
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
        double dt_pro = j->Time - i->Time;

        // Calculate time interval between point i and current time.
        double dt = t - i->Time;

        // Calculate the intermediate gas-phase mole fractions by linear
        // interpolation of the molar concentrations.
        double dens = 0.0;
        for (unsigned int k=0; k<gas.Species()->size(); ++k) {
            double dc = (j->Gas.MolarConc(k) - i->Gas.MolarConc(k)) * dt / dt_pro;
            gas.RawData()[k] = gas.MolarConc(k) + dc;
            dens += gas.RawData()[k];
        }
        gas.Normalise();

        // Now use linear interpolation to calculate the temperature.
        double dT = (j->Gas.Temperature() - i->Gas.Temperature()) * dt / dt_pro;
        gas.SetTemperature(gas.Temperature()+dT);

		// Interpolate the convective and thermophoretic velocities, and diffusion term
		double du =  (j->Gas.GetConvectiveVelocity() - i->Gas.GetConvectiveVelocity()) * dt / dt_pro;
		gas.SetConvectiveVelocity(gas.GetConvectiveVelocity() + du);
		double dv =  (j->Gas.GetThermophoreticVelocity() - i->Gas.GetThermophoreticVelocity()) * dt / dt_pro;
		gas.SetThermophoreticVelocity(gas.GetThermophoreticVelocity() + dv);
		double dD =  (j->Gas.GetDiffusionTerm() - i->Gas.GetDiffusionTerm()) * dt / dt_pro;
		gas.SetDiffusionTerm(gas.GetDiffusionTerm() + dD);

		//! Interpolate A4 rate of production
		double dwdotA4 = (j->Gas.PAHFormationRate() - i->Gas.PAHFormationRate()) * dt / dt_pro;
		gas.SetPAHFormationRate(gas.PAHFormationRate() + dwdotA4);

        // Now set the gas density, calculated using the values above.
        gas.SetDensity(dens);
    }

    // Give some indication of the data spacing
    if((j != m_gas_prof.end()) && ((j+1) != m_gas_prof.end()))
        return ((j+1)->Time - j->Time);
    else
    // Past the end of the data there is no spacing
        return std::numeric_limits<double>::max();
}

GasProfile* FlameSolver::Gasphase(void){
	return &m_gas_prof;
}
