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
  m_on(true)
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
  m_on(true)
{
    m_name = "Birth Process";
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
    if ((m_cell->Particles().Count() + m_cell->Particles().GetTotalParticleNumber())>0u) return true;
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
    unsigned int n_total = m_cell->Particles().Count() + m_cell->Particles().GetTotalParticleNumber();
    return A() * (double)n_total;
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

    int i = 0;
    // If all particles in ensemble, select a particle at random
    if (!(m_mech->IsHybrid()))
        i = m_cell->Particles().Select(rng);
    else
    {
        // Select a particle from the ensemble or particle-number list
        // ===========================================================
        // Get totals
        double ntotal_pn = (double)(m_cell->Particles().GetTotalParticleNumber());
        double ntotal_ens = (double)(m_cell->ParticleCount());

        // Here there should be a check that the index chosen is smaller than the threshold size
        // because nothing stops the stream having a larger threshold size than current system.
        // In that instance, particles could be added to the ensemble like with surface growth. 
        // However this cannot be done here easily because it requires construction of a new particles. 
        // Print warning message.
        if (m_cell->Particles().GetHybridThreshold() > sys.Particles().GetHybridThreshold())
            printf("sweep: Mixture PN threshold > reactor PN threshold; "
	           "could inflow particle that cannot be stored\n");

        // Select the particle
        boost::uniform_01<rng_type&, double> unifDistrib(rng);
        double test = unifDistrib() * (ntotal_pn + ntotal_ens);
        if (ntotal_pn >= test)
        {
            // Particle is chosen from the number list
            double repeats = F(sys);
            if (repeats != floor(repeats)) 
            {
                boost::random::bernoulli_distribution<double> decider(repeats);
                repeats = floor(repeats);
                if (decider(rng))
                    repeats += 1.0;
            }
            if (repeats > 0.0)
            {
                unsigned int index = m_mech->SetRandomParticle(m_cell->Particles(), t, test, iUniform, rng);
		// Check we found a valid index
                if (index > 0)
		{
		    sys.Particles().UpdateTotalsWithIndex(index, repeats);
		    sys.Particles().UpdateNumberAtIndex(index, (int)repeats);
		    sys.Particles().UpdateTotalParticleNumber((int)repeats);
		}
            }
            i = -1;
        }
        else
        {
            // Particle is chosen from the ensemble
            i = m_cell->Particles().Select_usingGivenRand(iUniform, test - ntotal_pn, rng);
        }

    }
    if (i >= 0)
    {
        DoParticleBirth(t, i, sys,
        m_cell->Particles().At(i)->getStatisticalWeight() * F(sys),
        rng);
    }
    // ===========================================================

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

    // Get reference to a particle.
    Sweep::Particle *sp = m_cell->Particles().At(isp)->Clone();

    // Reset some properties
    sp->resetCoagCount();
    sp->resetFragCount();
    sp->SetTime(t);     // Set LPDA update time.

    // Note that we could just adjust the weight of isp and add it (SWAs only),
    // but this makes the ensemble more sensitive to depletion for very
    // downstream reactors.

    double repeats = F(sys);
    int counter(-1);
    while (repeats > 0.0) {
        if (repeats >= 1.0) sys.Particles().Add(*(sp->Clone()), rng);
        else {
            boost::random::bernoulli_distribution<double> decider(repeats);
            if (decider(rng))
                counter = sys.Particles().Add(*sp, rng);
        }
        repeats -= 1.0;
    }

    // Delete the particle copy if it wasn't added by the decider
    if (counter < 0) delete sp;

}

// READ/WRITE/COPY.

// Creates a copy of the inception.
BirthProcess *const BirthProcess::Clone(void) const {return new BirthProcess(*this);}

// Returns the process type.  Used to identify different
// processes and for serialisation.
ProcessType BirthProcess::ID(void) const {return Birth_ID;}

