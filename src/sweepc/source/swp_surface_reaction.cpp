/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the SurfaceReaction class declared in the
    swp_surface_reaction.h header file.

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

#include "swp_surface_reaction.h"

#include "swp_particle_process.h"
#include "swp_mechanism.h"
#include "swp_process_type.h"
#include "swp_primary.h"

#include <cmath>
#include <stdexcept>
#include <iostream>

using namespace Sweep;
using namespace Sweep::Processes;
using namespace std;

const double Processes::SurfaceReaction::m_majfactor = 2.0;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
SurfaceReaction::SurfaceReaction(void)
: ParticleProcess(), m_arr(0.0,0.0,0.0)
{
    m_defer = true;
    m_name = "Surface Reaction";
}

// Initialising constructor.
SurfaceReaction::SurfaceReaction(const Sweep::Mechanism &mech)
: ParticleProcess(mech), m_arr(0.0,0.0,0.0)
{
    m_defer = true;
    m_name = "Surface Reaction";
}

// Copy constructor.
SurfaceReaction::SurfaceReaction(const SurfaceReaction &copy)
{
    *this = copy;
}

// Stream-reading constructor.
SurfaceReaction::SurfaceReaction(std::istream &in, const Sweep::Mechanism &mech)
{
    Deserialize(in, mech);
}

// OPERATOR OVERLOADS.

// Assignment operator.
SurfaceReaction &SurfaceReaction::operator =(const SurfaceReaction &rhs)
{
    if (this != &rhs) {
        ParticleProcess::operator =(rhs);
        m_arr     = rhs.m_arr;
        m_pid     = rhs.m_pid;
    }
    return *this;
}


// ARRHENIUS COEFFICIENTS.

// Returns the Arrhenius parameter.
Sprog::Kinetics::ARRHENIUS &SurfaceReaction::Arrhenius() {return m_arr;}
const Sprog::Kinetics::ARRHENIUS &SurfaceReaction::Arrhenius() const {return m_arr;}

// Sets the Arrhenius parameters.
void SurfaceReaction::SetArrhenius(Sprog::Kinetics::ARRHENIUS &arr) {m_arr = arr;}




// PARTICLE PROPERTY ID.

// Returns the ID number of the particle property to which
// the rate of this process is proportional.
PropID SurfaceReaction::PropertyID(void) const {return m_pid;}

// Sets the ID number of the particle property to which
// the rate of this process is proportional.
void SurfaceReaction::SetPropertyID(PropID pid)
{
    m_pid = pid;
}


// TOTAL RATE CALCULATIONS (ALL PARTICLES IN A SYSTEM).

// Returns rate of the process for the given system.
double SurfaceReaction::Rate(double t, const Cell &sys,
                           const Geometry::LocalGeometry1d &local_geom) const
{
    // Rate constant.
    double rate = m_arr.A;

    // Chemical species concentration dependence.
    rate *= chemRatePart(sys.GasPhase());

    // Temperature dependance.
    double T = sys.GasPhase().Temperature();
    rate *= pow(T, m_arr.n) * exp(-m_arr.E / (R * T));

    // Particle dependence.
    // Note that this assumes that the property m_pid is surface area - for the particle-number model
    rate *= (sys.Particles().GetSum(m_pid) + (PI * sys.Particles().GetTotalDiameter2()));

    if (m_mech->AnyDeferred()) {

        return rate * m_majfactor;
    } else {

        return rate;
    }
}


// SINGLE PARTICLE RATE CALCULATIONS.

// Returns the rate of the process for the given particle in
// the system. Process must be linear in particle number.
double SurfaceReaction::Rate(double t, const Cell &sys, const Particle &sp) const
{
    //! For the purpose of checking consistency with the spherical soot model
    //! solved using the method of moments with interpolative closure which
    //! assumes that only pyrene (A4) is able to incept and condense.
    //if(sp.Primary()->InceptedPAH()){
    //    double rate = 0.0;
    //    return rate;
    //}

    // Rate constant.
    double rate = m_arr.A;

    // Chemical species concentration dependence.
    rate *= chemRatePart(sys.GasPhase());

    // Temperature dependance.
    double T = sys.GasPhase().Temperature();
    rate *= pow(T, m_arr.n) * exp(-m_arr.E / (R * T));

    // Paticle dependence.
    rate *= sp.Property(m_pid);

    return rate;
}

// Return rate constant and chemistry part for hybrid method
double SurfaceReaction::Rate(double t, const Cell &sys) const
{
	// Rate constant.
	double rate = m_arr.A;

	// Chemical species concentration dependence.
	rate *= chemRatePart(sys.GasPhase());

	// Temperature dependance.
	double T = sys.GasPhase().Temperature();
	rate *= pow(T, m_arr.n) * exp(-m_arr.E / (R * T));
	
	return rate;
}

// Returns majorant rate of the process for the given system.
double SurfaceReaction::MajorantRate(double t, const Cell &sys,
                                   const Particle &sp) const
{
    return Rate(t, sys, sp) * m_majfactor;
}


// RATE TERM CALCULATIONS.
//   These routines return the individual rate terms for a
//   process, which may have multiple terms (e.g. condensation).

// Returns the number of rate terms for this process.
unsigned int SurfaceReaction::TermCount(void) const {return 1;}

// Calculates the rate terms given an iterator to a double vector. The
// iterator is advanced to the position after the last term for this
// process.
double SurfaceReaction::RateTerms(double t, const Cell &sys,
                                const Geometry::LocalGeometry1d &local_geom,
                                fvector::iterator &iterm) const
{
    return *(iterm++) = Rate(t, sys, local_geom);
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
int SurfaceReaction::Perform(double t, Sweep::Cell &sys, 
                             const Geometry::LocalGeometry1d& local_geom,
                             unsigned int iterm,
                             rng_type &rng) const
{
	PartPtrVector dummy;

    int i = sys.Particles().Select(static_cast<Sweep::PropID>(m_pid), rng);
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
                    times = sp->Adjust(m_dcomp, m_dvals, rng, 1);
                    sys.Particles().Update(i);

                    // Apply changes to gas-phase chemistry.
		    if (times > 0) {
                        if (!sys.GetIsAdiabaticFlag())
                            adjustGas(sys, sp->getStatisticalWeight());
                        else
                            adjustParticleTemperature(sys, sp->getStatisticalWeight(), 1, m_dcomp[0], 2);
                    }
                }
            } else {
                // If not valid then remove the particle.
                sys.Particles().Remove(i);
            }
        } else {
            // No particle update required, just perform the surface
            // reaction.
            times = sp->Adjust(m_dcomp, m_dvals, rng, 1);

            if (sp->IsValid()) {
                // Tell the binary tree to recalculate
                sys.Particles().Update(i);
            }
            else {
                // Particle has been removed due to oxidation
                sys.Particles().Remove(i);
            }

            // Apply changes to gas-phase chemistry.
            if (times > 0) {
                if (!sys.GetIsAdiabaticFlag())
                    adjustGas(sys, sp->getStatisticalWeight());
                else
                    adjustParticleTemperature(sys, sp->getStatisticalWeight(), 1, m_dcomp[0], 2);
            }
        }
    } else {
        // Failed to select a particle.
        return -1;
    }

    return 0;
}

// Performs the process on a given particle in the system.  Particle
// is given by index.  The process is performed n times.
int SurfaceReaction::Perform(double t, Cell &sys, Particle &sp, rng_type &rng,
                             unsigned int n) const
{
    unsigned int m = sp.Adjust(m_dcomp, m_dvals, rng, n);

    // Do normal update
    if (m > 0)
    {
        if (!sys.GetIsAdiabaticFlag())
            adjustGas(sys, sp.getStatisticalWeight(), m);
        else
            adjustParticleTemperature(sys, sp.getStatisticalWeight(), m, m_dcomp[0], 2);
    }
    return m;
}

// Surface growth updates for the hybrid particle model (particle-number updates)
// ==============================================================================

// Just do gas-phase adjustment for surface growth
int SurfaceReaction::Perform(double t, Cell &sys, rng_type &rng, unsigned int n) const
{
	if (!sys.GetIsAdiabaticFlag())
		adjustGas(sys, 1, n);
	else
		adjustParticleTemperature(sys, 1, n, m_dcomp[0], 2);
	return n;
}

// Performs the process on a given particle in the system.  Particle
// is given by index.  The process is performed n times.
int SurfaceReaction::Perform(double t, Cell &sys, Particle &sp, rng_type &rng, unsigned int n, bool isParticleNumberUpdate) const
{
	// Just update the particle and leave the gas-phase update
	unsigned int m = sp.Adjust(m_dcomp, m_dvals, rng, n);
	return m;
}

// ==============================================================================

// Adjusts a primary particle according to the rules of the reaction.
unsigned int SurfaceReaction::adjustPri(Sweep::AggModels::Primary &pri, rng_type &rng, unsigned int n) const
{
    return pri.Adjust(m_dcomp, m_dvals, rng, n);
}


// READ/WRITE/COPY.

// Creates a copy of the particle process.
SurfaceReaction *const SurfaceReaction::Clone(void) const
{
    return new SurfaceReaction(*this);
}

// Returns the process type.  Used to identify different
// processes and for serialisation.
ProcessType SurfaceReaction::ID(void) const {return SurfaceReaction_ID;}

// Writes the object to a binary stream.
void SurfaceReaction::Serialize(std::ostream &out) const
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
                               "(Sweep, SurfaceReaction::Serialize).");
    }
}

// Reads the object from a binary stream.
void SurfaceReaction::Deserialize(std::istream &in, const Sweep::Mechanism &mech)
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
                                    "(Sweep, SurfaceReaction::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, SurfaceReaction::Deserialize).");
    }
}
