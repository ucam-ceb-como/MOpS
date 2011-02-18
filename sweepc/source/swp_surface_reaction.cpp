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

#include <cmath>
#include <stdexcept>
#include <iostream>

using namespace Sweep;
using namespace Sweep::Processes;
using namespace std;

const real Processes::SurfaceReaction::m_majfactor = 2.0;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
SurfaceReaction::SurfaceReaction(void)
: ParticleProcess(), m_arr(0.0,0.0,0.0), m_pid(0)
{
    m_defer = true;
    m_name = "Surface Reaction";
}

// Initialising constructor.
SurfaceReaction::SurfaceReaction(const Sweep::Mechanism &mech)
: ParticleProcess(mech), m_arr(0.0,0.0,0.0), m_pid(0)
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
unsigned int SurfaceReaction::PropertyID(void) const {return m_pid;}

// Sets the ID number of the particle property to which
// the rate of this process is proportional.
void SurfaceReaction::SetPropertyID(unsigned int i)
{
    m_pid = i;
}


// TOTAL RATE CALCULATIONS (ALL PARTICLES IN A SYSTEM).

// Returns rate of the process for the given system.
real SurfaceReaction::Rate(real t, const Cell &sys) const
{
    // Rate constant.
    real rate = m_arr.A;

    // Chemical species concentration dependence.
    rate *= chemRatePart(sys.MoleFractions(), sys.Density());

    // Temperature dependance.
    real T = sys.Temperature();
    rate *= pow(T, m_arr.n) * exp(-m_arr.E / (R * T));

    // Particle dependence.
    rate *= sys.Particles().GetSum(static_cast<TreeCache::PropID>(m_pid));

    if (m_mech->AnyDeferred()) {
        return rate * m_majfactor;
    } else {
        return rate;
    }
}


// SINGLE PARTICLE RATE CALCULATIONS.

// Returns the rate of the process for the given particle in
// the system. Process must be linear in particle number.
real SurfaceReaction::Rate(real t, const Cell &sys, const Particle &sp) const
{
    // Rate constant.
    real rate = m_arr.A;

    // Chemical species concentration dependence.
    rate *= chemRatePart(sys.MoleFractions(), sys.Density());

    // Temperature dependance.
    real T = sys.Temperature();
    rate *= pow(T, m_arr.n) * exp(-m_arr.E / (R * T));

    // Paticle dependence.
    rate *= sp.Property(static_cast<ParticleCache::PropID>(m_pid));

    return rate;
}

// Returns majorant rate of the process for the given system.
real SurfaceReaction::MajorantRate(real t, const Cell &sys,
                                   const Particle &sp) const
{
    return Rate(t, sys, sp) * m_majfactor;
}


// RATE TERM CALCULATIONS.
//   These routines return the individual rate terms for a
//   process, which may have multiple terms (e.g. condensation).

// Returns the number of rate terms for this process.
unsigned int SurfaceReaction::TermCount(void) const {return 1;}

// Calculates the rate terms given an iterator to a real vector. The
// iterator is advanced to the position after the last term for this
// process.
real SurfaceReaction::RateTerms(real t, const Cell &sys,
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
 * \param[in,out]   rand_int    Pointer to function that generates uniform integers on a range
 * \param[in,out]   rand_u01    Pointer to function that generates U[0,1] deviates
 * \param[out]      out         Details of any particle being transported out of system
 *
 * \return      0 on success, otherwise negative.
 */
int SurfaceReaction::Perform(Sweep::real t, Sweep::Cell &sys, 
                             const Geometry::LocalGeometry1d& local_geom,
                             unsigned int iterm,
                             int (*rand_int)(int, int), 
                             Sweep::real(*rand_u01)(), 
                             Sweep::Transport::TransportOutflow *out) const
{
    int i = sys.Particles().Select(static_cast<TreeCache::PropID>(m_pid), rand_int, rand_u01);

    if (i >= 0) {
        Particle *sp = sys.Particles().At(i);


        // Update particle with deferred processes.
        if (m_mech->AnyDeferred()) {
            // Calculate majorant rate then update the particle.
            real majr = MajorantRate(t, sys, *sp);
            m_mech->UpdateParticle(*sp, sys, t, rand_u01);

            // Check that the particle is still valid.
            if (sp->IsValid()) {
                real truer = Rate(t, sys, *sp);

                if (!Fictitious(majr, truer, rand_u01)) {
                    // Adjust particle.
                    sp->Adjust(m_dcomp, m_dvals, 1);
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
            sp->Adjust(m_dcomp, m_dvals, 1);

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
int SurfaceReaction::Perform(real t, Cell &sys, Particle &sp,
                             unsigned int n) const
{
    unsigned int m = sp.Adjust(m_dcomp, m_dvals, n);
    adjustGas(sys, m);
    return m;
}

// Adjusts a primary particle according to the rules of the reaction.
unsigned int SurfaceReaction::adjustPri(Sweep::Primary &pri, unsigned int n) const
{
    return pri.Adjust(m_dcomp, m_dvals, n);
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
        unsigned int n = (unsigned int)m_pid;
        out.write((char*)&n, sizeof(n));
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
        unsigned int n = 0;

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
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_pid = n;
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
