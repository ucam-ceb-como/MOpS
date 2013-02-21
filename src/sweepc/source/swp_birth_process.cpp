/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the BirthProcess class declared in the
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
}

// Initialising constructor.
BirthProcess::BirthProcess(const Sweep::Mechanism &mech)
: Process(mech),
  m_cell(NULL),
  m_btype(BirthProcess::iStochastic),
  m_on(true)
{}

// Copy constructor.
BirthProcess::BirthProcess(const BirthProcess &copy)
: m_cell(copy.m_cell)
{
    *this = copy;
}

// Stream-reading constructor.
BirthProcess::BirthProcess(std::istream &in, const Sweep::Mechanism &mech)
{
    Deserialize(in, mech);
}

// Default destructor.
BirthProcess::~BirthProcess(void)
{}

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


// INFORMATION FOR THE SOLVER
// Does the Cell inflow have particles present?
bool BirthProcess::HasParticlesInCell() const {
    if (m_cell->Particles().Count()>0u) return true;
    else return false;
}

// TOTAL RATE CALCULATIONS.

// Get the cell-transfer scaling factor.
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
        if (m_cell == NULL)
            throw runtime_error("No cell specified for sampling."
                    " (Sweep, BirthProcess::Rate)");

        if (m_cell->ParticleCount() > 0u)
            return A() * m_cell->Particles().Count();
    }
    return 0.0;
}

// RATE TERM CALCULATIONS.

// Returns the number of rate terms for this process (one).
unsigned int BirthProcess::TermCount(void) const {return 1;}

// Calculates the rate terms given an iterator to a double vector. The
// iterator is advanced to the position after the last term for this
// process.  Returns the sum of all terms.
double BirthProcess::RateTerms(const double t, const Cell &sys,
                             const Geometry::LocalGeometry1d &local_geom,
                             fvector::iterator &iterm) const
{
    *iterm = Rate(t, sys, local_geom);
    return *(iterm++);
}

// PERFORMING THE PROCESS.

/*!
 * Deprecated Perform process.
 *
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

    // Only do if the process is turned-on and stochastic.
    if (m_btype == BirthProcess::iStochastic && m_on) {
        if (m_cell->ParticleCount() > 0u) {
            int i = m_cell->Particles().Select(rng);

            DoParticleBirth(t, i, sys,
                    m_cell->Particles().At(i)->getStatisticalWeight() * F(sys),
                    rng);
        }
    } else if (m_btype == BirthProcess::iContinuous)
        throw runtime_error("Perform should not be called when birth is continuous."
                " (Sweep, BirthProcess::Perform)");

    return 0;
}

/*!
 * Create particles over time dt.
 *
 * @param t     Current time of the system
 * @param dt    Time to remove particles over
 * @param sys   The system to do transport for
 * @param rng   Random number generator
 */
void BirthProcess::PerformDT (
        double t,
        double dt,
        Sweep::Cell &sys,
        rng_type &rng) const {

    if (m_btype == BirthProcess::iContinuous) {

        Process::PerformDT(dt, t, sys, rng);

        // Initialise some variables
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
        }
    }

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
    sp->setStatisticalWeight(wt);
    sp->resetCoagCount();
    sp->SetTime(t);     // Set LPDA update time.
    sys.Particles().Add(*sp, rng);

}

// READ/WRITE/COPY.

// Creates a copy of the inception.
BirthProcess *const BirthProcess::Clone(void) const {return new BirthProcess(*this);}

// Returns the process type.  Used to identify different
// processes and for serialisation.
ProcessType BirthProcess::ID(void) const {return Birth_ID;}

// Writes the object to a binary stream.
void BirthProcess::Serialize(std::ostream &out) const
{
    if (out.good()) {

        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Serialize base class.
        Process::Serialize(out);

    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, BirthProcess::Serialize).");
    }
}

// Reads the object from a binary stream.
void BirthProcess::Deserialize(std::istream &in, const Sweep::Mechanism &mech)
{
    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        switch (version) {
            case 0:
                // Deserialize base class.
                Process::Deserialize(in, mech);

            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, BirthProcess::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, BirthProcess::Deserialize).");
    }
}
