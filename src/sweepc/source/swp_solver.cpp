/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the Solver class declared in the
    swp_solver.h header file.

  Licence:
    This file is part of "sweepc".

    sweepc is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
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

#include "swp_solver.h"
#include "local_geometry1d.h"

#include "choose_index.hpp"

#include <stdlib.h>
#include <cmath>
#include "string_functions.h"
#include <stdexcept>
#include <ctime>
#include <limits>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/lognormal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

using namespace Sweep;
using namespace std;
using namespace Strings;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Solver::Solver(void)
: m_splitratio(1.0e9)
{
//  srnd(time(0));			//added by ms785
//    srnd(getpid());
 #ifdef USE_MPI
  Sweep::init_genrand(getpid());			//added by ms785
 #endif


}

//! Copy constructor
Solver::Solver(const Solver &sol)
: m_splitratio(sol.m_splitratio) {}

// Default destructor.
Solver::~Solver(void)
{
    // Nothing special to destruct.
}


// RUNNING THE SOLVER.

// Performs stochastic stepping algorithm up to specified stop time using
// the given mechanism to define the stochastic processes.  Updates given
// system accordingly.  On error returns <0, otherwise returns 0.
int Solver::Run(double &t, double tstop, Cell &sys, const Mechanism &mech,
                rng_type &rng)
{
    int err = 0;
    double tsplit, dtg, jrate, tflow(t);
    fvector rates(mech.TermCount(), 0.0);
    // Global maximum time step.
    dtg     = tstop - t;
	
    double tin = t; //store start time 

	// aab64 Variables used to shift incepting weight over time
	double nmax = sys.Particles().Capacity();
	double nnew = sys.ParticleCount();
	double nmin = 1.0; // minimum particles for which to activate scaling
	double wmax = 1.0; // maximum incepting weight
	double wmin = 1.0; // minimum incepting weight
	double wnew = 1.0; // incepting weight for next step
	double a = 1.0;    // constant in scaling
	double b = 1.0;    // constant in scaling
	double c = 1.0;    // constant in scaling
	char *wtfn = "L";  // inception weight function

    // Loop over time until we reach the stop time.
    while (t < tstop)
    {
		// aab64 Shift incepting particle weight over time, as the ensemble fills up
		/* Check if particle weight should be updated if SWA is in play
		Only step up weights after ensemble capacity hits a certain point nmin
		(possibly find a way to do this incrementally at a few points, rather than 
		continuously every time this process occurs). */
		if (mech.IsWeightedCoag() && mech.IsVariableWeightedInception()) {
			nnew = sys.ParticleCount();
			wmax = mech.GetMaxInceptionWeight();
			wmin = mech.GetMinInceptionWeight();
			nmin = mech.GetMinSPForAIWOnset();
			mech.GetWeightScalingFn(wtfn);
			wnew = wmin;
			// If the number of particles is large enough, update the incepting weight
			if (nnew > nmin) {
				// Compute and set new weighting (must occur before rate is calculated) 
				if (wtfn == "E") {
					// Exponential scaling
					b = log(wmax / wmin);
					b *= (1.0 / (nmax - nmin));
					a = wmin * exp(-1.0 * b * nmin);
					c = 0.0;
					wnew = (a * exp(b * nnew)) + c;
				}
				else if (wtfn == "Q") {
					// Quadratic scaling
					a = (wmax - wmin) / ((nmax * nmax) - (2.0 * nmax * nmin) + (nmin * nmin));
					b = -2.0 * a * nmin;
					c = wmin - (a * nmin * nmin) - (b * nmin);
					wnew = (a * nnew * nnew) + (b * nnew) + c;
				}
				else {
					// Linear scaling
					a = 0.0;
					b = ((wmax - wmin) / (nmax - nmin));
					c = wmin - (b * nmin);
					wnew = (b * nnew) + c;
				}				
			}
			// Set new incepting weight
			sys.SetInceptingWeight(wnew);
		}
		
		// aab64 If the ensemble is more than half full already / has passed a minimum average size, 
		// adjust inception factor by sampling from lognormal distribution
		// and constrain newIncFactor to [1,100]. 
		// Note that this should really rather check if a certain minimum size has been reached. 
		bool heavyAllowed = mech.GetIsHeavy();
		double dcol_lim = mech.GetHeavyValue();
		double dcol_ave = sys.Particles().GetSum(Sweep::iDW) / sys.Particles().GetSum(Sweep::iW);
		if (heavyAllowed)
		{
			double meanIFdist = 0.0;
			double stdIFdist = 1.0;
			//if (sys.ParticleCount() > ceil(1.0 * sys.Particles().Capacity() / 2.0))
			if (dcol_ave > dcol_lim)
			{
				double newIncFactor = boost::random::lognormal_distribution<double>(meanIFdist, stdIFdist)(rng);
				if (newIncFactor < 1.0)
					newIncFactor = 1.0;
				else if (newIncFactor > 100.0)
					newIncFactor = 100.0;
				sys.SetInceptionFactor((int)newIncFactor);
			}
		}

        if (mech.AnyDeferred() && (sys.ParticleCount() > 0))  {
            // Get the process jump rates (and the total rate).
            jrate = mech.CalcJumpRateTerms(t, sys, Geometry::LocalGeometry1d(), rates);

            // Calculate split end time.
            tsplit = calcSplitTime(t, std::min(t+dtg, tstop), jrate, sys.ParticleCount());
        } else {
            // There are no deferred processes, therefore there
            // is no need to perform LPDA splitting steps.
            tsplit = tstop;
        }

        // Perform stochastic jump processes.
        while (t < tsplit) {

			// aab64 Could also shift incepting weights here for intense inception rate escalation. 
			
            // Sweep does not do transport
            jrate = mech.CalcJumpRateTerms(t, sys, Geometry::LocalGeometry1d(), rates);
            timeStep(t, std::min(t + dtg / 3.0, tsplit), sys, Geometry::LocalGeometry1d(),
                     mech, rates, jrate, rng);

            // Do particle transport
            if (sys.OutflowCount() > 0 || sys.InflowCount() > 0)
                mech.DoParticleFlow(t, t - tflow, sys, Geometry::LocalGeometry1d(), rng);
            tflow = t;
        }

        sys.SetCurrentProcessTau(t - tin); // aab64 store time passed in current loop for heat transfer

        // Perform Linear Process Deferment Algorithm to
        // update all deferred processes.
	    mech.LPDA(t, sys, rng);

	    if (sys.ParticleCount() > 0 && sys.GetIsAdiabaticFlag())
	        mech.DoProcess(1000000, dtg, sys, Geometry::LocalGeometry1d(), rng); // aab64 Call DoProcess with special tag to do just heat transfer
    }

    return err;
}


// TIME STEPPING ROUTINES.

/*!
 * Calculates the splitting end time after which all particles
 * shall be updated using LPDA.
 *
 *@param[in]	t		Current time
 *@param[in]	tstop	Latest possible end for the splitting
 *@param[in]	jrate	Computational jump rate
 *@param[in]	n		Number of computational particles
 *
 *@return	End time for next splitting step
 */
double Solver::calcSplitTime(double t, double tstop, double jrate,
                           unsigned int n) const
{
    // Calculate the splitting time step, ensuring that it is
    // not longer than the maximum allowable time.
    double tsplit = (n + 1) * m_splitratio / (jrate + 1.0);

    // Now put the split end time into tsplit, again
    // checking that it is not beyond the stop time.
    return min(tsplit+t, tstop);
}

/*!
 * Performs a single stochastic event on the ensemble from the given
 * mechanism or move to t_stop, whichever comes first.
 *
 *@param[in,out]    t           Current time, which will be updated
 *@param[in]        t_stop      Time past which step may not go
 *@param[in,out]    sys         System in which jump will take place
 *@param[in]        geom        Specify size and neighbours of cell (use a default constructed object which will apply unit scaling when no geometry information present)
 *@param[in]        mech        Mechanism specifying the jump
 *@param[in]        rates       Vector of computational jump rates, one for each jump process
 *@param[in]        jrate       Sum of entries in rates (total jump rate)
 *@param[in,out]    rng         Random number generator
 *
 *@pre      t <= t_stop
 *@post     t <= t_stop
 */
void Solver::timeStep(double &t, double t_stop, Cell &sys, const Geometry::LocalGeometry1d &geom,
                      const Mechanism &mech, const fvector &rates, double jrate,
                      rng_type &rng)
{
    // The purpose of this routine is to perform a single stochastic jump process.  This
    // involves summing the total rate of all processes, generating a waiting time,
    // selecting a process and performing that process.
    double dt;

    //std::cout << "Solver::timeStep in cell " << &sys << " from " << t << " with rate " << jrate;

    // Calculate exponentially distributed time step size.
    if (jrate > 0.0) {
        boost::exponential_distribution<double> waitDistrib(jrate);
        boost::variate_generator<Sweep::rng_type&, boost::exponential_distribution<double> > waitGenerator(rng, waitDistrib);
        dt = waitGenerator();
        //std::cout << ' ' << dt;
    } else {
        // Avoid divide by zero.
        dt = std::numeric_limits<double>::max();
    }

    // Truncate if step is too long or select a process
    // to perform.
    if (t+dt <= t_stop) {
        boost::uniform_01<rng_type &> uniformGenerator(rng);
        const int i = chooseIndex(rates, uniformGenerator);
        mech.DoProcess(i, t+dt, sys, geom, rng);
        t += dt;
    } else {
        t = t_stop;
        //std:cout << " step truncated to end at " << t_stop;
    }
    //std::cout << std::endl;
}

// Selects a process using a DIV algorithm and the process rates
// as weights.
int Solver::chooseProcess(const fvector &rates, double (*rand_u01)())
{
    // This routine implements a DIV algorithm to select a process from a
    // discrete list of processes when each process's rate is given.

    // Add together all process rates.
    fvector::const_iterator i;
    double sum = 0.0;
    for (i=rates.begin(); i!=rates.end(); ++i) {
        sum += *i;
    }

    // Generate a uniform random number.
    double r = rand_u01() * sum ;

    // Select a process using DIV algorithm.
    int j=-1;
    i = rates.begin();
    while((r > 0.0) && (i != rates.end())) {
        r -= *i;
        ++i; ++j;
    }
    return j;
};
