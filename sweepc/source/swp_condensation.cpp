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
#include <cmath>
#include <stdexcept>

using namespace Sweep;
using namespace Sweep::Processes;
using namespace std;

const unsigned int Condensation::TERM_COUNT = 3;
const real Condensation::m_majfactor       = 2.0;
const real Condensation::m_efm             = 2.2;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (protected).
Condensation::Condensation(void)
: ParticleProcess(), m_a(1.0), m_kfm1(0.0), m_kfm2(0.0), m_kfm3(0.0)
{
    m_defer = true;
    m_name = "Condensation";
}

// Initialising constructor.
Condensation::Condensation(const Sweep::Mechanism &mech)
: ParticleProcess(mech), m_a(1.0), m_kfm1(0.0), m_kfm2(0.0), m_kfm3(0.0)
{
    // Assume the condensation is simulated as a deferred process (LPDA).
    m_defer = true;
    m_name = "Condensation";
}

// Copy constructor.
Condensation::Condensation(const Condensation &copy)
{
    *this = copy;
}

// Stream-reading constructor.
Condensation::Condensation(std::istream &in, const Sweep::Mechanism &mech)
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

// Returns the fixed rate constant.
real Condensation::A() const {return m_a;}

// Sets the fixed rate constant.
void Condensation::SetA(real a) {m_a = a;}

// Sets the coagulation kernel parameters given the mass and
// collision diameter of the condensing species.
void Condensation::SetCondensingSpecies(const real m, const real d)
{
    // Calculate the free-mol terms for condensation.  This must be done
    // before the condensation process is used.
    m_kfm3 = m_efm * CFM / sqrt(m);
    m_kfm2 = d * m_kfm3 * 2.0;
    m_kfm1 = d * m_kfm2 / 2.0;
}


// TOTAL RATE CALCULATIONS (ALL PARTICLES IN A SYSTEM).

// Returns rate of the process for the given system.
real Condensation::Rate(real t, const Cell &sys) const
{
    // Calculate temperature terms.
    real cterm = m_a * sqrt(sys.Temperature()) * NA;

     // Chemical species concentration dependence.
    cterm *= chemRatePart(sys.MoleFractions(), sys.Density());

    // Free molecular terms.
    cterm *= (m_kfm1 * sys.ParticleCount()) + 
             (m_kfm2 * sys.Particles().GetSum(ParticleCache::iDcol)) +
             (m_kfm3 * sys.Particles().GetSum(ParticleCache::iD2));

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
real Condensation::Rate(real t, const Cell &sys, const Particle &sp) const
{
    // Calculate temperature terms.
    real cterm = m_a * sqrt(sys.Temperature()) * NA;
//    real trm[3];

    // Chemical species concentration dependence.
    cterm *= chemRatePart(sys.MoleFractions(), sys.Density());

    // Free molecular terms.
//    trm[0] = cterm * m_kfm1;
//    trm[1] = cterm * (m_kfm2 * sp.CollDiameter());
//    trm[2] = cterm * (m_kfm3 * sp.CoagModelCache->CollDiamSquared());
    cterm *= m_kfm1 + 
             (m_kfm2 * sp.CollDiameter()) +
             (m_kfm3 * sp.CollDiamSquared());
    return cterm; //trm[0] + trm[1] + trm[2];
}

// Returns majorant rate of the process for the given system.
real Condensation::MajorantRate(real t, const Cell &sys, const Particle &sp) const
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

// Calculates the rate terms given an iterator to a real vector. The 
// iterator is advanced to the position after the last term for this
// process.
real Condensation::RateTerms(real t, const Cell &sys, 
                             fvector::iterator &iterm) const
{
    // Calculate temperature terms.
    real cterm = m_a * sqrt(sys.Temperature()) * NA;

     // Chemical species concentration dependence.
    cterm *= chemRatePart(sys.MoleFractions(), sys.Density());

    // If the mechanism contains any deferred processes then we must use the
    // majorant form of the rate, in order to account for any changes to
    // particles during deferred process updates.
    if (m_mech->AnyDeferred()) cterm *= m_majfactor;

    // Free molecular terms.
    real sum = 0.0;
    sum += *(iterm++) = m_kfm1 * cterm * sys.ParticleCount();
    sum += *(iterm++) = m_kfm2 * cterm * sys.Particles().GetSum(ParticleCache::iDcol);
    sum += *(iterm++) = m_kfm3 * cterm * 
                        sys.Particles().GetSum(ParticleCache::iD2);
    return sum;
}


// PERFORMING THE PROCESS.

// Performs the process on the given system.  The responsible rate term is given
// by index.  Returns 0 on success, otherwise negative.
int Condensation::Perform(const real t, Cell &sys, const unsigned int iterm) const
{
    // Select particle based on which term was called.
    int i  = -1;
    ParticleCache::PropID id = ParticleCache::iUniform;
    switch(iterm) {
        case 1:
            id = ParticleCache::iDcol;
            i  = sys.Particles().Select(id);
            break;
        case 2:
            id = ParticleCache::iD2;
            i  = sys.Particles().Select(id);
            break;
        case 0:
        default:
            id = ParticleCache::iUniform;
            i  = sys.Particles().Select();
            break;
    }

    // Check for a valid particle (i>=0).
    if (i >= 0) {
        Particle *sp = sys.Particles().At(i);

        real majr = MajorantRate(t, sys, *sp);

        if (m_mech->AnyDeferred()) {
            // Update particle with deferred processes.
            m_mech->UpdateParticle(*sp, sys, t);
        }

        // Check that the particle is still valid.
        if (sp->IsValid()) {
            // Get the true process rate (after updates).
            real truer = Rate(t, sys, *sp);

            // Check that the event is not ficticious by comparing the
            // majorant rate with the true rate.
            if (!Ficticious(majr, truer) || !m_mech->AnyDeferred()) {
                // Adjust particle.
                sp->Adjust(m_dcomp, m_dvals,SubModels::BasicModel_ID, 
                           ParticleCache::iS, 1); // Use surface area to weight in first instance.
                sys.Particles().Update(i);

                // Apply changes to gas-phase chemistry.
                adjustGas(sys);
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
int Condensation::Perform(real t, Cell &sys, Particle &sp,
                          unsigned int n) const
{
    unsigned int m = sp.Adjust(m_dcomp, m_dvals, SubModels::BasicModel_ID, 
                               ParticleCache::iS, n); // Use surface-area to weight.
    adjustGas(sys, m);
    return 0;
}

// Adjusts a primary particle according to the rules of the condensation.
unsigned int Condensation::adjustPri(Sweep::Primary &pri, unsigned int n) const
{
    return pri.Adjust(m_dcomp, m_dvals, n);
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

        // Write rate constant.
        double v = (double)m_a;
        out.write((char*)&v, sizeof(v));

        // Write free-mol parameters.
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

                // Read rate constant.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_a = (real)val;

                // Read free-mol parameter.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_kfm1 = (real)val;
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_kfm2 = (real)val;
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_kfm3 = (real)val;

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
