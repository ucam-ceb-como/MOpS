/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the DeathProcess class declared in the
    swp_death_process.h header file.

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
#include <boost/random/bernoulli_distribution.hpp>
#include "swp_death_process.h"
#include "swp_mechanism.h"
#include <stdexcept>

using namespace Sweep;
using namespace Sweep::Processes;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (protected).
DeathProcess::DeathProcess(void)
: Process(),
  m_dtype(DeathProcess::iDeathDelete),
  m_adaptive(false),
  m_toggled(false),
  m_cell(NULL)
{
    m_name = "Death Process";
}

// Initialising constructor.
DeathProcess::DeathProcess(const Sweep::Mechanism &mech)
: Process(mech),
  m_dtype(DeathProcess::iDeathDelete),
  m_adaptive(false),
  m_toggled(false),
  m_cell(NULL)
{
    m_name = "Death Process";
}

// Copy constructor.
DeathProcess::DeathProcess(const DeathProcess &copy)
{
    *this = copy;
}

// Stream-reading constructor.
DeathProcess::DeathProcess(std::istream &in, const Sweep::Mechanism &mech)
{
    Deserialize(in, mech);
}

// Default destructor.
DeathProcess::~DeathProcess(void)
{
}

// OPERATOR OVERLOADS.

// Assignment operator.
DeathProcess &DeathProcess::operator =(const DeathProcess &rhs)
{
    if (this != &rhs) {
        Process::operator =(rhs);

        m_dtype = rhs.m_dtype;
        m_adaptive = rhs.m_adaptive;
        m_toggled = rhs.m_toggled;
        m_cell = rhs.m_cell;
    }
    return *this;
}

// TOTAL RATE CALCULATIONS.

// Returns rate of the process for the given system.
double DeathProcess::Rate(double t, const Cell &sys,
                        const Geometry::LocalGeometry1d &local_geom) const
{
    throw std::runtime_error("DeathProcesses are no longer jump processes."
            " (Sweep::DeathProcess::Perform");
    return 0.0;
}

// RATE TERM CALCULATIONS.

// Returns the number of rate terms for this process (one).
unsigned int DeathProcess::TermCount(void) const {return 1;}

// Calculates the rate terms given an iterator to a double vector. The 
// iterator is advanced to the position after the last term for this
// process.  Returns the sum of all terms.
double DeathProcess::RateTerms(const double t, const Cell &sys,
                             const Geometry::LocalGeometry1d &local_geom,
                             fvector::iterator &iterm) const
{
    *iterm = Rate(t, sys, local_geom);
    return *(iterm++);
}

// PERFORMING THE PROCESS.

/*!
 * Deprecated death process.
 *
 * \param[in]       t           Time
 * \param[in,out]   sys         System to update
 * \param[in]       local_geom  Details of local phsyical layout
 * \param[in]       iterm       Process term responsible for this event
 * \param[in,out]   rng         Random number generator
 *
 * \return      0 on success, otherwise negative.
 */
int DeathProcess::Perform(double t, Sweep::Cell &sys, 
                          const Geometry::LocalGeometry1d& local_geom,
                          unsigned int iterm,
                          rng_type &rng) const
{
    throw std::runtime_error("DeathProcesses are no longer jump processes."
            " (Sweep::DeathProcess::Rate");
    return 0;
}

/*!
 * Remove particles over time dt.
 *
 * @param t     Current time of the system
 * @param dt    Time to remove particles over
 * @param sys   The system to do transport for
 * @param rng   Random number generator
 */
void DeathProcess::PerformDT (
        const double t,
        const double dt,
        Sweep::Cell &sys,
        rng_type &rng) const {

    Process::PerformDT(t, dt, sys, rng);

    if (m_dtype == DeathProcess::iDeathRescale) {
        // Don't delete anything, just rescale the sample volume
        sys.AdjustSampleVolume(1.0/(1.0 - (dt) * A()));
    } else {
        // Set up some variables
        double weightToRemove = sys.Particles().GetSum(iW) * dt * A();
        double wt(0.0);
        int i(0);

        // Remove particles up to the maximum weight
        while (weightToRemove > 0.0) {
            i = sys.Particles().Select(rng);
            wt = sys.Particles().At(i)->getStatisticalWeight();

            if (wt < weightToRemove) {
                DoParticleDeath(t, i, sys, rng);
            } else {
                boost::random::bernoulli_distribution<double> decider(weightToRemove/wt);
                if (decider(rng)) {
                    DoParticleDeath(t, i, sys, rng);
                }
            }

            weightToRemove -= wt;
        }
    }
}

/*!
 * Carry out the death of a particle.
 *
 * @param t     System time
 * @param isp   The index of the particle to remove
 * @param sys   System to remove the particle from
 * @param rng   Random number generator
 */
void DeathProcess::DoParticleDeath(
        const double t,
        const int isp,
        Sweep::Cell &sys,
        rng_type &rng) const {
    Sweep::Particle *sp = sys.Particles().At(isp);

    if (m_dtype == DeathProcess::iDeathDelete) {
        // Just delete the particle
        sys.Particles().Remove(isp, true);

    } else if (m_dtype == DeathProcess::iDeathMove) {
        // Move it to a downstream cell
        if (m_cell == NULL)
            throw std::runtime_error("No cell to move the particle to!"
                    " (Sweep::DeathProcess::DoParticleDeath).");

        double wtFactor = (double)m_cell->Particles().Capacity() * sys.SampleVolume()
                    / ((double)sys.Particles().Capacity() * m_cell->SampleVolume());
        sp->setStatisticalWeight(sp->getStatisticalWeight() / wtFactor);
        sp->SetTime(t);

        // Remove the particle from the ensemble but don't delete
        sys.Particles().Remove(isp, false);

        // Now add to the new ensemble
        m_cell->Particles().Add(*sp, rng);
    }
}

/*!
 * Changes the nature of the death process of the current cell and the inflow
 * processes of any downstream cells. If N(t) = Nmax, we can reduce the
 * number of ensemble contractions by moving outflow particles directly
 * downstream.
 *
 * @param sys
 */
void DeathProcess::Adapt(Sweep::Cell &sys) {

    if (m_adaptive && !m_toggled) {
        if (sys.Particles().Count() == sys.Particles().Capacity()) {
            // We have a full ensemble.

            std::cout << "sweep: Enabling move-on-death and turning-off inflow." << std::endl;

            // Change this Death process to move particles
            m_dtype = DeathProcess::iDeathMove;

            // Turn-off the downstream birth processes
            Sweep::Processes::BirthPtrVector bps = m_cell->Inflows();
            Sweep::Processes::BirthPtrVector::iterator i;
            for (i = bps.begin(); i != bps.end(); ++i) {
                (*i)->SetProcessSwitch(false);
            }

            // Indicate that this process has been toggled-on
            m_toggled = true;
        }
    }

    if (m_adaptive && m_toggled) {
        if (sys.Particles().Count() <
                (unsigned int) (0.8 * (double) sys.Particles().Capacity())) {
            // Uh oh, our ensemble has depleted to 80% capacity (arbitrary)

            // Change to ensemble rescaling to refill the ensemble
            m_dtype = DeathProcess::iDeathRescale;

            // Now re-enable downstream birth processes
            Sweep::Processes::BirthPtrVector bps = m_cell->Inflows();
            Sweep::Processes::BirthPtrVector::iterator i;
            for (i = bps.begin(); i != bps.end(); ++i) {
                (*i)->SetProcessSwitch(true);
            }

            // Toggle off
            m_toggled = false;
        }
    }

}

// READ/WRITE/COPY.

// Creates a copy of the inception.
DeathProcess *const DeathProcess::Clone(void) const {return new DeathProcess(*this);}

// Returns the process type.  Used to identify different
// processes and for serialisation.
ProcessType DeathProcess::ID(void) const {return Death_ID;}
