/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the Condensation class declared in the
    swp_condensation.h header file.

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

#include "swp_condensation.h"
#include "swp_mechanism.h"
#include "swp_primary.h"

#include <cmath>
#include <stdexcept>
#include <cassert>

using namespace Sweep;
using namespace Sweep::Processes;
using namespace std;

const unsigned int Condensation::TERM_COUNT = 3;
const double Condensation::m_majfactor       = 2.0;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (protected).
Condensation::Condensation(void)
: ParticleProcess(), m_efm(2.2), m_kfm1(0.0), m_kfm2(0.0), m_kfm3(0.0)
{
    m_defer = true;
    m_name = "Condensation";
}

// Initialising constructor.
Condensation::Condensation(const Sweep::Mechanism &mech)
: ParticleProcess(mech), m_efm(mech.GetEnhancementFM()),
 m_kfm1(0.0), m_kfm2(0.0), m_kfm3(0.0)
{
    // Assume the condensation is simulated as a deferred process (LPDA).
    m_defer = true;
    m_name = "Condensation";
}

// Copy constructor.
Condensation::Condensation(const Condensation &copy)
: m_efm(copy.m_efm)
{
    *this = copy;
}

// Stream-reading constructor.
Condensation::Condensation(std::istream &in, const Sweep::Mechanism &mech)
: m_efm(mech.GetEnhancementFM())
{
    Deserialize(in, mech);
}

// Default destructor.
Condensation::~Condensation(void)
{
    // Nothing special to destruct.
}


// OPERATOR OVERLOADS.

// Assignment operator.
Condensation &Condensation::operator=(const Condensation &rhs)
{
    if (this != &rhs) {
        ParticleProcess::operator =(rhs);
        m_a    = rhs.m_a;
        m_kfm1 = rhs.m_kfm1;
        m_kfm2 = rhs.m_kfm2;
        m_kfm3 = rhs.m_kfm3;
    }
    return *this;
}


// RATE CONSTANT AND PARAMETERS.

// Sets the coagulation kernel parameters given the mass and
// collision diameter of the condensing species.
void Condensation::SetCondensingSpecies(const double m, const double d)
{
    // Calculate the free-mol terms for condensation.  This must be done
    // before the condensation process is used.
    m_kfm3 = m_efm * CFM / sqrt(m);
    m_kfm2 = d * m_kfm3 * 2.0;
    m_kfm1 = d * m_kfm2 / 2.0;
}


// TOTAL RATE CALCULATIONS (ALL PARTICLES IN A SYSTEM).

/*!
 *@param[in]            t           Time at which rate is being calculated
 *@param[in]            sys         System for which rate is to be calculated
 *@param[in]            local_geom  Spatial configuration information (ignored)
 *
 *@return   Process rate
 */
double Condensation::Rate(double t, const Cell &sys,
                        const Geometry::LocalGeometry1d &local_geom) const
{
	double cterm = 0.0;
	if(sys.ParticleCount() > 0) {
		// Calculate temperature terms.
		cterm = m_a * sqrt(sys.GasPhase().Temperature()) * NA;

		 // Chemical species concentration dependence.
		cterm *= chemRatePart(sys.GasPhase());

		// Free molecular terms.
		cterm *= (m_kfm1 * sys.ParticleCount()) +
				 (m_kfm2 * sys.Particles().GetSum(Sweep::iDcol)) +
				 (m_kfm3 * sys.Particles().GetSum(Sweep::iD2));
	}
	//else cterm was initialised to 0, so nothing to do

    // If the mechanism contains any deferred processes then we must use the
    // majorant form of the rate, in order to account for any changes to
    // particles during deferred process updates.
    if (m_mech->AnyDeferred()) {
        return cterm * m_majfactor;
    } else {
        return cterm;
    }
}


// SINGLE PARTICLE RATE CALCULATIONS.

// Returns the rate of the process for the given particle in
// the system. Process must be linear in particle number.
double Condensation::Rate(double t, const Cell &sys, const Particle &sp) const
{
    // Calculate temperature terms.
    double cterm = m_a * sqrt(sys.GasPhase().Temperature()) * NA;
//    double trm[3];

    // Chemical species concentration dependence.
    cterm *= chemRatePart(sys.GasPhase());

    // Get particle property
    const double d = sp.CollDiameter();
    // Free molecular terms.
//    trm[0] = cterm * m_kfm1;
//    trm[1] = cterm * (m_kfm2 * sp.CollDiameter());
//    trm[2] = cterm * (m_kfm3 * sp.CoagModelCache->CollDiamSquared());
    cterm *= m_kfm1 + 
             (m_kfm2 * d) +
             (m_kfm3 * d * d);
    return cterm; //trm[0] + trm[1] + trm[2];
}

// Returns majorant rate of the process for the given system.
double Condensation::MajorantRate(double t, const Cell &sys, const Particle &sp) const
{
    // Return the single particle rate multiplied by the 
    // condensation majorant factor.
    return Rate(t, sys, sp) * m_majfactor;
}


// RATE TERM CALCULATIONS.
//   These routines return the individual rate terms for a 
//   process, which may have multiple terms (e.g. condensation).

// Returns the number of rate terms for this process.
unsigned int Condensation::TermCount(void) const {return TERM_COUNT;}

// Calculates the rate terms given an iterator to a double vector. The 
// iterator is advanced to the position after the last term for this
// process.
double Condensation::RateTerms(double t, const Cell &sys,
                             const Geometry::LocalGeometry1d &local_geom,
                             fvector::iterator &iterm) const
{
    // Calculate temperature terms.
    double cterm = m_a * sqrt(sys.GasPhase().Temperature()) * NA;

     // Chemical species concentration dependence.
    cterm *= chemRatePart(sys.GasPhase());

    // If the mechanism contains any deferred processes then we must use the
    // majorant form of the rate, in order to account for any changes to
    // particles during deferred process updates.
    if (m_mech->AnyDeferred()) cterm *= m_majfactor;

    // Free molecular terms.
    double sum = 0.0;
    sum += *(iterm++) = m_kfm1 * cterm * sys.ParticleCount();
    sum += *(iterm++) = m_kfm2 * cterm * sys.Particles().GetSum(Sweep::iDcol);
    sum += *(iterm++) = m_kfm3 * cterm * 
                        sys.Particles().GetSum(Sweep::iD2);
    return sum;
}


// PERFORMING THE PROCESS.

/*!
 * 
 *
 * \param[in]       t           Time
 * \param[in,out]   sys         System to update
 * \param[in]       local_geom  Details of local phsyical layout
 * \param[in]       iterm       Process term responsible for this event
 * \param[in,out]   rng         Random number generator
 *
 * \return      0 on success, otherwise negative.
 */
int Condensation::Perform(double t, Sweep::Cell &sys, 
                          const Geometry::LocalGeometry1d& local_geom,
                          unsigned int iterm,
                          rng_type &rng) const
{
	PartPtrVector dummy;
    // Select particle based on which term was called.
    int i  = -1;
    Sweep::PropID id = Sweep::iUniform;
    switch(iterm) {
        case 1:
            id = Sweep::iDcol;
            i  = sys.Particles().Select(id, rng);
            break;
        case 2:
            id = Sweep::iD2;
            i  = sys.Particles().Select(id, rng);
            break;
        case 0:
        default:
            id = Sweep::iUniform;
            i  = sys.Particles().Select(rng);;
            break;
    }

    // Check for a valid particle (i>=0).
    if (i >= 0) {
        Particle *sp = sys.Particles().At(i);

        double majr = MajorantRate(t, sys, *sp);

        if (m_mech->AnyDeferred()) {
            // Update particle with deferred processes.
            m_mech->UpdateParticle(*sp, sys, t, i, rng, dummy);
        }

        // Check that the particle is still valid.
        if (sp->IsValid()) {
            // Get the true process rate (after updates).
            double truer = Rate(t, sys, *sp);

            // Check that the event is not ficticious by comparing the
            // majorant rate with the true rate.

            // The reason for the second part of the || expression is that,
            // in the absence of any deferred processes, the pyrene jump rate
            // is calculated without majorant factor.  However, in this method,
            // the fictitious jump test always uses the majorant factor and so
            // some events would wrongly be treated as fictitious.  The ugly
            // solution used here is to ensure that there are no fictitious
            // events when there are no deferred processes.
            if (!Fictitious(majr, truer, rng) || !m_mech->AnyDeferred()) {
                // Adjust particle.
                sp->Adjust(m_dcomp, m_dvals, rng, 1);
                sys.Particles().Update(i);

                // Apply changes to gas-phase chemistry.
                adjustGas(sys, sp->getStatisticalWeight());
            }
        } else {
            // If not valid then remove the particle.
            sys.Particles().Remove(i);
        }
    } else {
        // Failed to select a particle.
        return -1;
    }
    return 0;
}

// Performs the process on a given particle in the system.  Particle
// is given by index.  The process is performed n times.
int Condensation::Perform(double t, Cell &sys, Particle &sp, rng_type &rng,
                          unsigned int n) const
{
    unsigned int m = sp.Adjust(m_dcomp, m_dvals, rng, n);
    adjustGas(sys, sp.getStatisticalWeight(), m);
    return 0;
}

// Adjusts a primary particle according to the rules of the condensation.
unsigned int Condensation::adjustPri(Sweep::AggModels::Primary &pri, rng_type &rng, unsigned int n) const
{
    return pri.Adjust(m_dcomp, m_dvals, rng, n);
}


// READ/WRITE/COPY.

// Creates a copy of the particle process.
Condensation *const Condensation::Clone(void) const
{
    return new Condensation(*this);
}

// Returns the process type.  Used to identify different
// processes and for serialisation.
ProcessType Condensation::ID(void) const {return Condensation_ID;}

// Writes the object to a binary stream.
void Condensation::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Serialize base class.
        ParticleProcess::Serialize(out);

        // Write free-mol parameters.
        double v(0.0);
        v = (double)m_kfm1;
        out.write((char*)&v, sizeof(v));
        v = (double)m_kfm2;
        out.write((char*)&v, sizeof(v));
        v = (double)m_kfm3;
        out.write((char*)&v, sizeof(v));
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, Condensation::Serialize).");
    }
}

// Reads the object from a binary stream.
void Condensation::Deserialize(std::istream &in, const Sweep::Mechanism &mech)
{
    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        double val = 0.0;

        switch (version) {
            case 0:
                // Deserialize base class.
                ParticleProcess::Deserialize(in, mech);

                // Read free-mol parameter.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_kfm1 = (double)val;
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_kfm2 = (double)val;
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_kfm3 = (double)val;

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, Condensation::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, Condensation::Deserialize).");
    }
}
