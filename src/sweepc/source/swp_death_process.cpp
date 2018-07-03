 /*!
  * @file   swp_death_process.h
  * @author Matthew Celnik, William Menz
  * @brief  Implementation of a death process
  *
  *   Licence:
  *      sweepc is free software; you can redistribute it and/or
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
  m_dtype(DeathProcess::iContRescale),
  m_toggled(false),
  m_cell(NULL)
{
    m_name = "Death Process";
}

/*!
 * Usual way of creating a death process.
 *
 * @param mech      Mechanism defining system
 * @return          The new death process
 */
DeathProcess::DeathProcess(const Sweep::Mechanism &mech)
: Process(mech),
  m_dtype(DeathProcess::iContRescale),
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

// OPERATOR OVERLOADS.

// Assignment operator.
DeathProcess &DeathProcess::operator =(const DeathProcess &rhs)
{
    if (this != &rhs) {
        Process::operator =(rhs);

        m_dtype = rhs.m_dtype;
        m_toggled = rhs.m_toggled;
        m_cell = rhs.m_cell;
    }
    return *this;
}

/*!
 * @param t     Death process type
 */
void DeathProcess::SetDeathType(const DeathType t)
{
    m_dtype = t;
}

/*!
 * @return      The type of death process
 */
DeathProcess::DeathType DeathProcess::GetDeathType() const {
    return m_dtype;
}

/*!
 * @param c     The downstream cell to move particles to
 */
void DeathProcess::SetCell(Cell* c) {
    m_cell = c;
}

/*!
 * @return      Whether this is a jump-style death process
 */
bool DeathProcess::IsStochastic() const {
    bool ans(false);
    if (m_dtype == DeathProcess::iStochDelete ||
            m_dtype == DeathProcess::iStochMove) ans = true;
    return ans;
}


// TOTAL RATE CALCULATIONS.

/*!
 * Returns the death rate. Will return the actual rate if stochastic is ON,
 * otherwise will return zero so that the process is not selected (the
 * 'continuous' case).
 *
 * @param t             Time
 * @param sys           System to evaluate the rate for
 * @param local_geom    System geometry
 * @return              Outflow rate
 */
double DeathProcess::Rate(double t, const Cell &sys,
                        const Geometry::LocalGeometry1d &local_geom) const
{
    if (IsStochastic()) {
        return InternalRate(t, sys, local_geom);
    }
    // If not stochastic, just return zero.
    return 0.0;
}

/*!
* @param t             Time
* @param sys           System to evaluate the rate for
* @param local_geom    System geometry
* @return              Outflow rate
*/
double DeathProcess::InternalRate(
        double t,
        const Cell &sys,
        const Geometry::LocalGeometry1d &local_geom) const {
	// aab64 if weighted particles, use sum of weights in place of N
	if (m_dtype == DeathProcess::iWtdDelete)
		return m_a * sys.Particles().GetSum(iW);
	else
	{
		unsigned int n_total = sys.Particles().Count();

		// aab64 for hybrid particle model
		if (m_mech->IsHybrid() && sys.GetIncepted() > 0)
		{
			n_total += (unsigned int)sys.GetIncepted();               // Account for particles in the incepting class
		}
		return m_a * n_total;
	}
}

// RATE TERM CALCULATIONS.

// Returns the number of rate terms for this process (one).
unsigned int DeathProcess::TermCount(void) const {return 1;}

/*!
 *
 * @param t             Time at which rate is being calculated
 * @param sys           System for which rate is to be calculated
 * @param local_geom    Spatial configuration information (ignored)
 * @param iterm         Iterator on rates vector
 * @return              Rate of death process
 */
double DeathProcess::RateTerms(const double t, const Cell &sys,
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
 * \param[in]       local_geom  Details of local physical layout
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
    // Get particle index
	int i = 0;
	// aab64 for hybrid particle model
	if (!(m_mech->IsHybrid()))
	    i = sys.Particles().Select(rng);
	else
	{
		// Check if should remove from ensemble or bin
		boost::uniform_01<rng_type&, double> unifDistrib(rng);
		double test = unifDistrib();
		if (test * (sys.ParticleCount() + sys.GetIncepted()) <= sys.GetIncepted())
		{
			m_mech->SetRandomParticle(true, sys, t, (test / (1.0 - test)) * sys.ParticleCount(), true, sys.GetMuLN(), sys.GetSigmaLN(),
				sys.Particles().GetInceptedSP_youngest().CollDiameter(), sys.Particles().GetInceptedSP_oldest().CollDiameter(), 0.0, rng);
			sys.AdjustIncepted(-1.0);                                                        // Reduce the incepting class count
			sys.AdjustInceptingCoagulations();                                               // Increment number of times particles have left the incepting class
			sys.AdjustInceptingCoagulations_tmp();
			sys.SetCoagulationSums(sys.Particles().GetInceptedSP_tmp_d_1().CollDiameter());  // Store the change in total diameter due to losing this particle from the class
			i = -1;
		}
		else
		{
			i = sys.Particles().Select_usingGivenRand(iUniform, test * (sys.ParticleCount() + sys.GetIncepted()) - sys.GetIncepted(), rng);
		}
	}

    if (i >= 0) DoParticleDeath(t, i, sys, rng);

    return 0;
}

// aab64 Weighted perform
// If weighted particles, use weights to randomly choose particles and
// scale/remove particle depending on delta w required for the interval. 
// Use when inception rate is high to free ensemble space for new particles.
/*!
* \param[in]       t           Time
* \param[in,out]   sys         System to update
* \param[in]       local_geom  Details of local physical layout
* \param[in]       iterm       Process term responsible for this event
* \param[in,out]   rng         Random number generator
* \param[in,out]   num         Amount of weight change still required
*
* \return      0 on success, otherwise negative.
*/
int DeathProcess::Perform_wtd(double t, Sweep::Cell &sys,
	const Geometry::LocalGeometry1d& local_geom,
	unsigned int iterm,
	rng_type &rng, 
	double &num) const
{
	// Get particle index
	int i = sys.Particles().Select(iW, rng);

	// aab64 for hybrid particle model
	if (i >= 0)
	{
		Sweep::Particle *sp = sys.Particles().At(i);
		double sp_wt = sp->getStatisticalWeight();
		if (sp_wt <= num)
		{
			// we can keep this particle and 
			// just reduce its weight
			if (sp_wt > 1.0)
			{
				sp->setStatisticalWeight(sp_wt - 1.0);
				sys.Particles().Update(i);
				if (m_mech->IsHybrid() && i == 0)
					sys.AdjustIncepted(-1.0);                  // If incepting class, track the change
				num = num - 1.0;
			}
			// remove the particle to free up 
			// ensemble space for new particles
			else
			{
				if (m_mech->IsHybrid() &&                      // If this is particle hybrid model,
					i == 0 &&                                  // and the incepting class is selected,
					sys.ParticleCount() > 1 &&                 // and N>1 => can't delete this particle,
					sys.GetIncepted() > 0)                     // and class not empty so need to do something...
				{
					sys.AdjustIncepted(-sp_wt);
					sp->setStatisticalWeight(sp_wt - sp_wt);   // Empty the class 
					sys.Particles().Update(i);
				}
				else
					DoParticleDeath(t, i, sys, rng);           // Just remove the particle
				num = num - sp_wt;
			}
		}
		// particle weight is larger than 
		// required change, just rescale
		else
		{
			sp->setStatisticalWeight(sp_wt - num);
			sys.Particles().Update(i);
			if (m_mech->IsHybrid() && i == 0)
				sys.AdjustIncepted(-num);                      // If incepting class, track the change
			num = 0;
		}
	}
	return 0;
}

/*!
 * Remove particles over time dt.
 *
 * @param t     Current time of the system
 * @param dt    Time to remove particles over
 * @param sys   The system to do transport for
 * @param local_geom    Geometry of the system
 * @param rng   Random number generator
 */
void DeathProcess::PerformDT (
        const double t,
        const double dt,
        Sweep::Cell &sys,
        const Geometry::LocalGeometry1d& local_geom,
        rng_type &rng) const {

    // Only do if set to 'continuous' mode.
    if (!IsStochastic()) 
	{
        Process::PerformDT(t, dt, sys, local_geom, rng);

        if (m_dtype == DeathProcess::iContRescale
                || (m_dtype == DeathProcess::iContAdaptive && !m_toggled)) 
		{
            // Don't delete anything, just rescale the sample volume
            sys.AdjustSampleVolume(1.0/(1.0 - (dt) * A()));
        } 
		else 
		{
            // Get the rate of the process
            double rate = InternalRate(t, sys, local_geom) * dt;
            if (rate > 0.0) {
				boost::random::poisson_distribution<unsigned, double> rpt(rate);
				
				if (m_dtype == DeathProcess::iWtdDelete)
				{
                    // aab64 replacement for cdelete for weighted particles 
				    // should be adjusted to allow multiple outflows
					// we know apriori how much the weight needs to change
					// from the mass balance and can effect this change
					// by scaling the weights
				    double num = sys.Particles().GetSum(iW) * A() * dt;
				    while (num > 0 && sys.ParticleCount() > 0) 
					    Perform_wtd(t, sys, local_geom, 0, rng, num);
				}
				else
				{
				    unsigned num = rpt(rng);
				    while (num > 0) 
					{
					    // Do the process to the particle.
					    Perform(t, sys, local_geom, 0, rng);
					    num--;
						if (m_mech->IsHybrid() && sys.Particles().IsFirstSP())
						{
							m_mech->MomentUpdate(t, 0, sys, rng);
						}
				    }
				}
            }
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

    if (m_dtype == DeathProcess::iContDelete
            || m_dtype == DeathProcess::iStochDelete
			|| m_dtype == DeathProcess::iWtdDelete) {
		sys.Particles().Remove(isp, true);                   // Just delete the particle
    } else if (m_dtype == DeathProcess::iContMove
            || m_dtype == DeathProcess::iStochMove
            || (m_dtype == DeathProcess::iContAdaptive && m_toggled)) {
        // Move it to a downstream cell
        if (m_cell == NULL)
            throw std::runtime_error("No cell to move the particle to!"
                    " (Sweep::DeathProcess::DoParticleDeath).");

        // Only adjust the particle weight if SWA coagulation.
        const double F = (double)m_cell->Particles().Capacity() * sys.SampleVolume()
                    / ((double)sys.Particles().Capacity() * m_cell->SampleVolume());

        sp->SetTime(t);

        // Add some copies of the particle
        double repeats = 1.0/F;
        while (repeats > 0.0) {
            if (repeats >= 1.0) m_cell->Particles().Add(*(sp->Clone()), rng);
            else {
                boost::random::bernoulli_distribution<double> decider(repeats);
                if (decider(rng)) m_cell->Particles().Add(*(sp->Clone()), rng);
            }
            repeats -= 1.0;
        }
        // Remove the particle from the ensemble and delete the original
        sys.Particles().Remove(isp, true);

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

    if (m_dtype == DeathProcess::iContAdaptive && !m_toggled) {
        if (sys.Particles().Count() == sys.Particles().Capacity()) {
            // We have a full ensemble.

            std::cout << "sweep: Enabling move-on-death and turning-off downstream inflow."
                    << std::endl;

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

    if (m_dtype == DeathProcess::iContAdaptive && m_toggled) {
        if (sys.Particles().Count() <
                (unsigned int) (0.8 * (double) sys.Particles().Capacity())) {
            // Uh oh, our ensemble has depleted to 80% capacity (arbitrary)

            std::cout << "sweep: Re-enabling downstream inflow and rescaling ensemble."
                    << std::endl;

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

// Returns the process type
ProcessType DeathProcess::ID(void) const {return Death_ID;}
