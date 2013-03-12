 /*!
  * @file   swp_birth_process.cpp
  * @author Matthew Celnik, William Menz
  * @brief  Implementation of a birth process
  *
  *   Licence:
  *      mops is free software; you can redistribute it and/or
  *      modify it under the terms of the GNU Lesser General Public License
  *      as published by the Free Software Foundation; either version 2
  *      of the License, or (at your option) any later version.
  *
  *      This program is distributed in the hope that it will be useful,
  *      but WITHOUT ANY WARRANTY; without even the implied warranty of
  *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *      GNU Lesser General Public License for more details.
  *
  *      You should have received a copy of the GNU Lesser General Public
  *      License along with this program; if not, write to the Free Software
  *      Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  *      02111-1307, USA.
  *
  *   Contact:
  *      Prof Markus Kraft
  *      Dept of Chemical Engineering
  *      University of Cambridge
  *      New Museums Site
  *      Pembroke Street
  *      Cambridge
  *      CB2 3RA, UK
  *
  *      Email:       mk306@cam.ac.uk
  *      Website:     http://como.cheng.cam.ac.uk
  */

#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
#include "swp_birth_process.h"
#include "swp_mechanism.h"
#include <stdexcept>

using namespace Sweep;
using namespace Sweep::Processes;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (protected).
BirthProcess::BirthProcess(void)
: m_cell(NULL),
  m_btype(BirthProcess::iStochastic),
  m_on(true),
  m_ptype(Processes::Weighted_Transition_Coagulation_ID)
{
    m_name = "Birth Process";
}

/*!
 * Initialising constructor
 *
 * @param mech  The mechanism defining the process
 * @return      A new BirthProcess
 */
BirthProcess::BirthProcess(const Sweep::Mechanism &mech)
: Process(mech),
  m_cell(NULL),
  m_btype(BirthProcess::iStochastic),
  m_on(true),
  m_ptype(Processes::Weighted_Transition_Coagulation_ID)
{
    m_name = "Birth Process";

    // Get the coagulation process type from the mechanism.
    const Processes::CoagPtrVector &coags = mech.Coagulations();
    if (coags.size() > 0) m_ptype = coags[0]->ID();
}

/*!
 * Copy constructor
 *
 * @param copy  Process to copy
 * @return      A new BirthProcess, copy of the original
 */
BirthProcess::BirthProcess(const BirthProcess &copy)
: m_cell(copy.m_cell)
{
    *this = copy;
}

// OPERATOR OVERLOADS.

// Assignment operator.
BirthProcess &BirthProcess::operator =(const BirthProcess &rhs)
{
    if (this != &rhs) {
        Process::operator =(rhs);
        // Copy pointer to sampling ensemble.
        m_cell = rhs.m_cell;
        m_btype = rhs.m_btype;
        m_on = rhs.m_on;
        m_ptype = rhs.m_ptype;
    }
    return *this;
}

//! Set the Cell from which this process samples.
void BirthProcess::SetCell(Cell* c) {
    m_cell = c;
}

//! Sets the birth process type
void BirthProcess::SetBirthType(const BirthType t) {
    m_btype = t;
}

//! Turn the process on or off
void BirthProcess::SetProcessSwitch(const bool s) {
    m_on = s;
}


// INFORMATION FOR THE SOLVER
// Does the Cell inflow have particles present?
bool BirthProcess::HasParticlesInCell() const {
    if (m_cell->Particles().Count()>0u) return true;
    else return false;
}

// TOTAL RATE CALCULATIONS.

/*!
 * Weights are scaled by this quantity when particles move from cell to cell.
 *
 * @param sys   System to evaluate scaling factor for
 * @return      The cell-transfer scaling factor
 */
double BirthProcess::F(const Cell &sys) const {
    return sys.SampleVolume() / m_cell->SampleVolume();
}

/*!
 *@param[in]            t           Time at which rate is being calculated
 *@param[in]            sys         System for which rate is to be calculated
 *@param[in]            local_geom  Spatial configuration information (ignored)
 *
 *@return   Process rate
 */
double BirthProcess::Rate(double t, const Cell &sys,
                        const Geometry::LocalGeometry1d &local_geom) const

{
    if (m_btype == BirthProcess::iStochastic && m_on) {
        return InternalRate(t, sys, local_geom);
    }
    return 0.0;
}

/*!
 *@param[in]            t           Time at which rate is being calculated
 *@param[in]            sys         System for which rate is to be calculated
 *@param[in]            local_geom  Spatial configuration information (ignored)
 *
 *@return   Process rate
 */
double BirthProcess::InternalRate(
        double t,
        const Cell &sys,
        const Geometry::LocalGeometry1d &local_geom) const {
    if (m_cell == NULL) throw runtime_error("No cell specified for sampling."
                " (Sweep, BirthProcess::InternalRate)");

    return A() * (double) m_cell->Particles().Count();
}

// RATE TERM CALCULATIONS.

// Returns the number of rate terms for this process (one).
unsigned int BirthProcess::TermCount(void) const {return 1;}

/*!
 *
 * @param t             Time at which rate is being calculated
 * @param sys           System for which rate is to be calculated
 * @param local_geom    Spatial configuration information (ignored)
 * @param iterm         Iterator on rates vector
 * @return              Rate of birth process
 */
double BirthProcess::RateTerms(const double t, const Cell &sys,
                             const Geometry::LocalGeometry1d &local_geom,
                             fvector::iterator &iterm) const
{
    *iterm = Rate(t, sys, local_geom);
    return *(iterm++);
}

// PERFORMING THE PROCESS.

/*!
 * \param[in]       t           Time
 * \param[in,out]   sys         System to update
 * \param[in]       local_geom  Details of local phsyical layout
 * \param[in]       iterm       Process term responsible for this event
 * \param[in,out]   rng         Random number generator
 *
 * \return      0 on success, otherwise negative.
 */
int BirthProcess::Perform(double t, Sweep::Cell &sys, 
                          const Geometry::LocalGeometry1d& local_geom,
                          unsigned int iterm,
                          rng_type &rng) const
{
    if (m_cell == NULL)
        throw runtime_error("No cell specified for sampling."
            " (Sweep, BirthProcess::Perform)");

    int i = m_cell->Particles().Select(rng);

    DoParticleBirth(t, i, sys,
        m_cell->Particles().At(i)->getStatisticalWeight() * F(sys),
        rng);

    return 0;
}

/*!
 * Create particles over time dt.
 *
 * @param t     Current time of the system
 * @param dt    Time to remove particles over
 * @param sys   The system to do transport for
 * @param local_geom    Geometry of the system
 * @param rng   Random number generator
 */
void BirthProcess::PerformDT (
        double t,
        double dt,
        Sweep::Cell &sys,
        const Geometry::LocalGeometry1d& local_geom,
        rng_type &rng) const {

    if (m_btype == BirthProcess::iContinuous) {

        Process::PerformDT(dt, t, sys, local_geom, rng);

        // Get the rate, and do n times like a LPDA process
        double rate = InternalRate(t, sys, local_geom) * dt;
        if (rate > 0.0) {
            boost::random::poisson_distribution<unsigned, double> rpt(rate);
            unsigned num = rpt(rng);
            while (num > 0) {
                // Do the process to the particle.
                Perform(t, sys, local_geom, 0, rng);
                num--;
            }
        }

        // Initialise some variables
        /*
        double weightToAdd = m_cell->Particles().GetSum(iW) * dt * A();
        const double f = F(sys);
        const double div = std::max(0.001 / (dt * A()),  f);
        double wt(0.0);
        int i(0);

        // Add particles up to the maximum weight
        while (weightToAdd > 0.0) {
            i = m_cell->Particles().Select(rng);
            wt = m_cell->Particles().At(i)->getStatisticalWeight() / div;

            if (wt < weightToAdd) {
                DoParticleBirth(t, i, sys, wt * f, rng);

            } else {
                // Use a bernoulli distribution to decide whether to add the
                // particle
                boost::random::bernoulli_distribution<double> decider(weightToAdd / wt);
                if (decider(rng)) DoParticleBirth(t, i, sys, wt * f, rng);

            }
            weightToAdd -= wt;
        }*/
    }
    // Don't do anything for the iStochastic case, as the Perform() will be called
    // in the usual manner for jump processes.

}

/*!
 * Create the particle in this system's cell.
 *
 * @param t     Time to create particle at
 * @param isp   Index of particle to clone
 * @param sys   System to put particle into
 * @param wt    Weight of the particle
 * @param rng   Random number generator
 */
void BirthProcess::DoParticleBirth(
        const double t,
        const int isp,
        Sweep::Cell &sys,
        const double wt,
        rng_type &rng) const {

    // Make a copy of the sampled particle
    Sweep::Particle *sp = m_cell->Particles().At(isp)->Clone();

    // Adjust its weight and add
    if (IsWeighted(m_ptype)) sp->setStatisticalWeight(wt);
    sp->resetCoagCount();
    sp->SetTime(t);     // Set LPDA update time.

    if (IsWeighted(m_ptype)) {
        // If it's a weighted process, just add the particle right away.
        sys.Particles().Add(*sp, rng);
    } else {
        // Otherwise, add some copies of the particle
        double repeats = F(sys);
        while (repeats > 0.0) {
            if (repeats >= 1.0) sys.Particles().Add(*(sp->Clone()), rng);
            else {
                boost::random::bernoulli_distribution<double> decider(repeats);
                if (decider(rng)) sys.Particles().Add(*sp, rng);
            }
            repeats -= 1.0;
        }
    }


}

// READ/WRITE/COPY.

// Creates a copy of the inception.
BirthProcess *const BirthProcess::Clone(void) const {return new BirthProcess(*this);}

// Returns the process type.  Used to identify different
// processes and for serialisation.
ProcessType BirthProcess::ID(void) const {return Birth_ID;}

