/*
  Author(s):	  Shraddha Shekar (ss663)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2010 Shraddha Shekar.

  File purpose:
    Implementation of the InterParticle class declared in the
    swp_InterParticle.h header file.

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

#include "swp_silica_interparticle.h"
#include "swp_silica_primary.h"
#include "swp_sintering_model.h"
#include "swp_particle_model.h"
#include "swp_mechanism.h"
#include "swp_process_type.h"

#include <cmath>
#include <stdexcept>
#include <cassert>

using namespace Sweep;
using namespace Sweep::Processes;
using namespace std;

const real InterParticle::m_majfactor = 2.0;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (protected).
InterParticle::InterParticle(void)
: ParticleProcess(), m_arr(0.0,0.0,0.0)
{
    m_defer = true;
    m_name = "InterParticle";
}

// Initialising constructor.
InterParticle::InterParticle(const Sweep::Mechanism &mech)
: ParticleProcess(mech), m_arr(0.0,0.0,0.0)
{
    // Assume the InterParticle is simulated as a deferred process (LPDA).
    m_defer = true;
    m_name = "InterParticle";
}

// Copy constructor.
InterParticle::InterParticle(const InterParticle &copy)
{
    *this = copy;
}

// Stream-reading constructor.
InterParticle::InterParticle(std::istream &in, const Sweep::Mechanism &mech)
{
    Deserialize(in, mech);
}

// Default destructor.
InterParticle::~InterParticle(void)
{
    // Nothing to destruct.
}


// OPERATOR OVERLOADS.

// Assignment operator.
InterParticle &InterParticle::operator=(const InterParticle &rhs)
{
    if (this != &rhs) {
        ParticleProcess::operator =(rhs);
        m_arr    = rhs.m_arr;
		m_pid     = rhs.m_pid;
    }
    return *this;
}


// RATE CONSTANT AND PARAMETERS.

// Returns the Arrhenius parameter.
Sprog::Kinetics::ARRHENIUS &InterParticle::Arrhenius() {return m_arr;}
const Sprog::Kinetics::ARRHENIUS &InterParticle::Arrhenius() const {return m_arr;}

// Sets the fixed rate constant.
void InterParticle::SetArrhenius(Sprog::Kinetics::ARRHENIUS &arr) {m_arr = arr;}


// PARTICLE PROPERTY ID.

// Returns the ID number of the particle property to which
// the rate of this process is proportional.
unsigned int InterParticle::PropertyID(void) const {return m_pid;}


// Sets the ID number of the particle property to which
// the rate of this process is proportional.
void InterParticle::SetPropertyID(PropID pid)
{
    m_pid = pid;
}


// TOTAL RATE CALCULATIONS (ALL PARTICLES IN A SYSTEM).

// Returns rate of the process for the given system.
real InterParticle::Rate(real t, const Cell &sys) const
{

	// From theoretical calculations: R_int = R_surf - sintering_contribution + R_inc
    real rate = 0.0;

	// First calculate surface reaction contribution:
	real T = sys.Temperature();

	// Get the total number of OH sites from cache
	int numOH = sys.Particles().GetSum(static_cast<Sweep::PropID>(m_pid));

	// Rate of surface reaction
	real R_surf = m_arr.A*chemRatePart(sys.MoleFractions(), sys.Density())*pow(T, m_arr.n)
			* exp(-m_arr.E / (R * T)) * numOH;

	// Forward-declare the total sintering rate
	real total_sint_rate = 0;

	// Get total surface area from cache
	double surface = sys.Particles().GetSum(static_cast<Sweep::PropID>(Sweep::iS));

	// Check if particle surface area exists. If not, set the sintering rate to zero.
	if (surface == 0) {
		surface = 1;
		total_sint_rate = 0;
	} else {
		// Get total sintering rate from cache.
		total_sint_rate = sys.Particles().GetSum(static_cast<Sweep::PropID>(Sweep::iSintRate));
	}

	// Calculate the total surface density of sites
	real rho_s =  numOH/surface;
	//Rate of sintering
	real R_sint = (rho_s*total_sint_rate)/(2.0);

	rate += (R_surf - R_sint);

	if(rate < 0)
		rate = 0;

	if (m_mech->AnyDeferred())
	{
        return rate * m_majfactor;
    }
	else
	{
        return rate;
    }
}



// SINGLE PARTICLE RATE CALCULATIONS.

// Returns the rate of the process for the given particle in
// the system. Process must be linear in particle number.
real InterParticle::Rate(real t, const Cell &sys, const Particle &sp) const
{

	// Rate constant.
    real rate = m_arr.A;

    // Chemical species concentration dependence.
    rate *= chemRatePart(sys.MoleFractions(), sys.Density());

    // Temperature dependance.
    real T = sys.Temperature();
    rate *= pow(T, m_arr.n) * exp(-m_arr.E / (R * T));

    // Forward-declare some parameters
    real sint_rate = 0;
    double rho_s = 0;

	// Get the number of OH sites from cache
	int numOH = sp.Property(static_cast<Sweep::PropID>(m_pid));

	// Get surface area from cache
	double surface = sp.Property(static_cast<Sweep::PropID>(Sweep::iS));

	if(surface == 0) {
		rho_s = 1;
		sint_rate = 0;
	} else {
		rho_s = numOH / surface;
		//Calculate sintering contribution
		sint_rate = (sp.Property(static_cast<Sweep::PropID>(Sweep::iSintRate)));
	}

	real sint_comp = (sint_rate*rho_s)/(2.0);

	// Sintering rate dependency
	rate -= sint_comp;

	if(rate<0)
		rate = 0;

    return rate; //trm[0] + trm[1] + trm[2];
}

// Returns majorant rate of the process for the given system.
real InterParticle::MajorantRate(real t, const Cell &sys, const Particle &sp) const
{
    // Return the single particle rate multiplied by the
    // majorant factor.
    return Rate(t, sys, sp) * m_majfactor;
}


// RATE TERM CALCULATIONS.
//   These routines return the individual rate terms for a
//   process, which may have multiple terms (e.g. InterParticle).

// Returns the number of rate terms for this process.
unsigned int InterParticle::TermCount(void) const {return 1;}

// Calculates the rate terms given an iterator to a real vector. The
// iterator is advanced to the position after the last term for this
// process.
real InterParticle::RateTerms(real t, const Cell &sys,
                             fvector::iterator &iterm) const
{
	 return *(iterm++) = Rate(t, sys);
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
int InterParticle::Perform(Sweep::real t, Sweep::Cell &sys,
                             const Geometry::LocalGeometry1d& local_geom,
                             unsigned int iterm,
                             rng_type &rng) const
{
    int i = sys.Particles().Select(static_cast<Sweep::PropID>(m_pid), rng);

    if (i >= 0) {
        Particle *sp = sys.Particles().At(i);


        // Update particle with deferred processes.
        if (m_mech->AnyDeferred()) {
            // Calculate majorant rate then update the particle.
            real majr = MajorantRate(t, sys, *sp);
            m_mech->UpdateParticle(*sp, sys, t, rng);

            // Check that the particle is still valid.
            if (sp->IsValid()) {
                real truer = Rate(t, sys, *sp);

                if (!Fictitious(majr, truer, rng)) {
                    // Adjust particle.
                    sp->AdjustIntPar(m_dcomp, m_dvals, rng, 1);
                    sys.Particles().Update(i);

                    // Apply changes to gas-phase chemistry.
                    adjustGas(sys);
                }
            } else {
                // If not valid then remove the particle.
                sys.Particles().Remove(i);
            }
        } else {
            // No particle update required, just perform the surface
            // reaction.
            sp->AdjustIntPar(m_dcomp, m_dvals, rng, 1);

            if (sp->IsValid()) {
                // Tell the binary tree to recalculate
                sys.Particles().Update(i);
            }
            else {
                // Particle has been removed due to oxidation
                sys.Particles().Remove(i);
            }

            // Apply changes to gas-phase chemistry.
            adjustGas(sys);
        }
    } else {
        // Failed to select a particle.
        return -1;
    }

    return 0;
}

// Performs the process on a given particle in the system.  Particle
// is given by index.  The process is performed n times.

int InterParticle::Perform(real t, Cell &sys, Particle &sp, rng_type &rng,
                          unsigned int n) const
{
    unsigned int m = sp.AdjustIntPar(m_dcomp, m_dvals, rng, n);
    adjustGas(sys, m);
    return m;
}

// Returns the process type.  Used to identify different
// processes and for serialisation.
ProcessType InterParticle::ID(void) const {return InterParticle_ID;}

// Creates a copy of the particle process.
InterParticle *const InterParticle::Clone(void) const
{
    return new InterParticle(*this);
}

// Writes the object to a binary stream.
void InterParticle::Serialize(std::ostream &out) const
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
        unsigned int n = (unsigned int)m_pid;
        out.write((char*)&n, sizeof(n));

    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, InterParticle::Serialize).");
    }
}

// Reads the object from a binary stream.
void InterParticle::Deserialize(std::istream &in, const Sweep::Mechanism &mech)
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
                m_arr.A = (real)A;
                m_arr.n = (real)nn;
                m_arr.E = (real)E;

                // Read particle property ID.
                in.read(reinterpret_cast<char*>(&m_pid), sizeof(m_pid));
                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, InterParticle::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, InterParticle::Deserialize).");
    }
}
