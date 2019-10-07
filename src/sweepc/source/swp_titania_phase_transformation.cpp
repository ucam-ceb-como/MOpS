/*!
 * @file    swp_titania_phase_transformation.cpp
 * @author  Casper Lindberg
 * @brief   Titanium Dioxide crystal phase transformation
 *
 *   Author(s):      Casper Lindberg
 *   Project:        sweepc (population balance solver)
 *   Copyright (C) 2017 Casper Lindberg
 *
 *   File purpose:
 *      
 *
 *   Licence:
 *      This file is part of "sweepc".
 *
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
#include "swp_titania_phase_transformation.h"
#include "swp_mechanism.h"
#include "swp_process_type.h"
#include "swp_particle_model.h"

#include <cmath>

using namespace std;

namespace Sweep {

namespace Processes {

const double TitaniaPhaseTransformation::m_majfactor = 2.0;

// Default constructor
TitaniaPhaseTransformation::TitaniaPhaseTransformation()
:ParticleProcess(), m_arr(0.0,0.0,0.0) 
{
	m_defer = false;
    m_name = "PhaseTransformation";
}

/*!
 * Mechanism constructor
 *
 * @param mech  Mechanism to create Process with
 * @return      Initialised Process
 */
TitaniaPhaseTransformation::TitaniaPhaseTransformation(
        const Sweep::Mechanism &mech)
: ParticleProcess(mech), m_arr(0.0,0.0,0.0)
{
	m_defer = false;
    m_name = "PhaseTransformation";
}


/*!
 * Copy constructor
 *
 * @param copy  The object to be copied
 * @return      A copy of the object
 */
TitaniaPhaseTransformation::TitaniaPhaseTransformation(const TitaniaPhaseTransformation &copy)
{
    *this = copy;
}

/*!
 * Stream-reading constructor
 *
 * @param in    Input binary stream
 * @param mech  Mechanism to create Process with
 * @return      Initialised Process
 */
TitaniaPhaseTransformation::TitaniaPhaseTransformation(
        std::istream &in,
        const Sweep::Mechanism &mech
        )
{
    Deserialize(in, mech);
}

// Assignment operator
TitaniaPhaseTransformation &TitaniaPhaseTransformation::operator =(const TitaniaPhaseTransformation &rhs)
{
    if (this != &rhs) {
        ParticleProcess::operator =(rhs);
		m_arr    = rhs.m_arr;
		m_pid    = rhs.m_pid;
    }
    return *this;
}

// Creates a copy of the particle process
TitaniaPhaseTransformation *const TitaniaPhaseTransformation::Clone(void) const
{
    return new TitaniaPhaseTransformation(*this);
}

// PARTICLE PROPERTY ID.

// Returns the ID number of the particle property to which
// the rate of this process is proportional.
PropID TitaniaPhaseTransformation::PropertyID(void) const {return m_pid;}

// Sets the ID number of the particle property to which
// the rate of this process is proportional.
void TitaniaPhaseTransformation::SetPropertyID(PropID pid)
{
    m_pid = pid;
}

// RATE CONSTANT AND PARAMETERS.

//! Returns the Arrhenius parameter.
Sprog::Kinetics::ARRHENIUS &TitaniaPhaseTransformation::Arrhenius() {return m_arr;}
const Sprog::Kinetics::ARRHENIUS &TitaniaPhaseTransformation::Arrhenius() const {return m_arr;}

//! Sets the fixed rate constant.
void TitaniaPhaseTransformation::SetArrhenius(Sprog::Kinetics::ARRHENIUS &arr) {m_arr = arr;}


// RATE TERM CALCULATIONS.
//   These routines return the individual rate terms for a
//   process, which may have multiple terms

//! Returns the number of rate terms for this process.
unsigned int TitaniaPhaseTransformation::TermCount(void) const {return 1;}

/*!
 * @brief       Returns the majorant rate.
 * 
 * @param[in]   t       Time at which process occurs
 * @param[in]   sys     System for rate calculation
 * @param[in]   sp      Particle for majorant rate calculation
 * 
 * @return      Rate of process
 */
double TitaniaPhaseTransformation::MajorantRate(double t, const Cell &sys, const Particle &sp) const
{
    // Return the single particle rate multiplied by the
    // majorant factor.
    return Rate(t, sys, sp) * m_majfactor;
}

/*!
 * @brief       Passes the system rate to an iterator
 * 
 * Calculates the rate terms given an iterator to a double vector. The
 * iterator is advanced to the position after the last term for this
 * process.
 * 
 * @return      Rate of process
 */
double TitaniaPhaseTransformation::RateTerms(double t, const Cell &sys, const Geometry::LocalGeometry1d &local_geom,
                             fvector::iterator &iterm) const
{
    return *(iterm++) = Rate(t, sys, local_geom);
}

/*!
 * Returns the rate for the whole system
 *
 * @param t             Time
 * @param sys           The cell to calculate the rate for
 * @param local_geom    Location information
 * @return              Value of the rate
 */
double TitaniaPhaseTransformation::Rate(
        double t,
        const Cell &sys,
        const Geometry::LocalGeometry1d &local_geom
        ) const
{
	//rate constant
	double rate = 3.0*m_arr.A;

	//temperature dependence
	double T = sys.GasPhase().Temperature();
	rate *= exp(-m_arr.E / (R * T));

	//particle dependence
	rate *= sys.Particles().GetSum(m_pid);

	return rate;
}

/*!
 * Calculates the single-particle rate
 *
 * @param t             Time
 * @param sys           The cell to calculate the rate for
 * @param sp            Particle for which to calculate rate
 * @return              Value of the rate
 */
double TitaniaPhaseTransformation::Rate(
        double t,
        const Cell &sys,
        const Particle &sp
        ) const
{
	//rate expressed as the change in number of components
	//rate constant
	double rate = 3.0*m_arr.A;

	//temperature dependence
	double T = sys.GasPhase().Temperature();
	rate *= exp(-m_arr.E / (R * T));

	//particle dependence
	rate *= sp.Property(m_pid);

	return rate;
}

// PERFORMING THE PROCESS.

/*!
 * @brief       Performs the process system-wide
 * 
 * Selects a particle, calculates the rate, and updates the system
 * accordingly.
 *
 * @param[in]       t           Time
 * @param[in,out]   sys         System to update
 * @param[in]       local_geom  Details of local physical layout
 * @param[in]       iterm       Process term responsible for this event
 * @param[in,out]   rng         Random number generator
 *
 * @return      0 on success, otherwise negative.
 */
int TitaniaPhaseTransformation::Perform(double t, Sweep::Cell &sys,
                             const Geometry::LocalGeometry1d& local_geom,
                             unsigned int iterm,
                             rng_type &rng) const
{
	int i = sys.Particles().Select(static_cast<Sweep::PropID>(m_pid), rng);

	PartPtrVector dummy;

	unsigned int times;

    if (i >= 0) {
        Particle *sp = sys.Particles().At(i);


        // Update particle with deferred processes.
        if (m_mech->AnyDeferred()) {
            // Calculate majorant rate then update the particle.
            double majr = MajorantRate(t, sys, *sp);
			m_mech->UpdateParticle(*sp, sys, t, i, rng, dummy);

            // Check that the particle is still valid.
            if (sp->IsValid()) {
                double truer = Rate(t, sys, *sp);

                if (!Fictitious(majr, truer, rng)) {
                    // Adjust particle.
                    times = sp->AdjustPhase(m_dcomp, m_dvals, rng, 1);
                    sys.Particles().Update(i);
                }
            } else {
                // If not valid then remove the particle.
                sys.Particles().Remove(i);
            }
        } else {
            // No particle update required, just perform the process
            times = sp->AdjustPhase(m_dcomp, m_dvals, rng, 1);

            if (sp->IsValid()) {
                // Tell the binary tree to recalculate
                sys.Particles().Update(i);
            }
            else {
                // Particle has been removed due to oxidation
                sys.Particles().Remove(i);
            }

        }
    } else {
        // Failed to select a particle.
        return -1;
    }

    return 0;
}

/*!
 * @brief       Performs process on a given particle n times
 * 
 * @param[in]   t   Time for process
 * @param[in]   sys System in which to act
 * @param[in]   sp  Particle to adjust
 * @param[in]   rng Random number generator
 * @param[in]   n   Number of times to do process
 */
int TitaniaPhaseTransformation::Perform(double t, Cell &sys, Particle &sp, rng_type &rng,
                          unsigned int n) const
{
		unsigned int m = sp.AdjustPhase(m_dcomp, m_dvals, rng, n);
		return m;
}


//! Returns the process type.
ProcessType TitaniaPhaseTransformation::ID(void) const {return TitaniaPhase_ID;}

/*!
 * Writes the object to a binary stream.
 *
 * @param out   Output binary stream
 */
void TitaniaPhaseTransformation::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Serialize base class.
        ParticleProcess::Serialize(out);

        // Write arrhenius coefficients.
        double A  = (double)m_arr.A;
        double nn = (double)m_arr.n;
        double E  = (double)m_arr.E;
        out.write((char*)&A, sizeof(A));
        out.write((char*)&nn, sizeof(nn));
        out.write((char*)&E, sizeof(E));

        // Write particle property ID.
        out.write((char*)&m_pid, sizeof(m_pid));
    } else {
        throw invalid_argument("Output stream not ready "
                               "in TitaniaPhaseTransformation::Serialize");
    }
}

/*!
 * Reads object from a binary stream
 *
 * @param in    Input binary stream
 * @param mech  Mechanism for process
 */
void TitaniaPhaseTransformation::Deserialize(
        std::istream &in,
        const Sweep::Mechanism &mech
        )
{
    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        double A = 0.0, nn = 0.0, E = 0.0;

        switch (version) {
            case 0:
                // Deserialize base class.
                ParticleProcess::Deserialize(in, mech);

                // Read arrhenius coefficients.
                in.read(reinterpret_cast<char*>(&A), sizeof(A));
                in.read(reinterpret_cast<char*>(&nn), sizeof(nn));
                in.read(reinterpret_cast<char*>(&E), sizeof(E));
                m_arr.A = (double)A;
                m_arr.n = (double)nn;
                m_arr.E = (double)E;

                // Read particle property ID.
                in.read(reinterpret_cast<char*>(&m_pid), sizeof(m_pid));
                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, TitaniaPhaseTransformation::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "in TitaniaPhaseTransformation::Deserialize");
    }
}

} // Processes

} // Sweep