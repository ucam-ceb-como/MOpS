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

#include "swp_birth_process.h"
#include "swp_mechanism.h"
#include <stdexcept>

using namespace Sweep;
using namespace Sweep::Processes;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (protected).
BirthProcess::BirthProcess(void)
: m_particle(NULL), m_cell(NULL)
{
}

// Initialising constructor.
BirthProcess::BirthProcess(const Sweep::Mechanism &mech)
: Process(mech), m_particle(NULL), m_cell(NULL)
{
    // Create a default inflow particle.
    m_particle = new Particle(0.0, mech);
}

// Copy constructor.
BirthProcess::BirthProcess(const BirthProcess &copy)
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
{
    delete m_particle;
}

// OPERATOR OVERLOADS.

// Assignment operator.
BirthProcess &BirthProcess::operator =(const BirthProcess &rhs)
{
    if (this != &rhs) {
        Process::operator =(rhs);
        // Copy default particle.
        delete m_particle;
        if (rhs.m_particle) {
            m_particle = rhs.m_particle->Clone();
        } else {
            m_particle = NULL;
        }
        // Copy pointer to sampling ensemble.
        m_cell = rhs.m_cell;
    }
    return *this;
}


// TOTAL RATE CALCULATIONS.

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
    if (m_cell) {
        return m_a * m_cell->ParticleCount(); // Rate depends on birth cell.
    } else {
        return m_a; // Constant rate.
    }
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
    if (m_cell) {
        *iterm = m_a * m_cell->ParticleCount(); // Rate depends on birth cell.
    } else {
        *iterm = m_a; // Constant rate.
    }
    return *(iterm++);
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
int BirthProcess::Perform(double t, Sweep::Cell &sys, 
                          const Geometry::LocalGeometry1d& local_geom,
                          unsigned int iterm,
                          rng_type &rng) const
{
    Particle *p = NULL;

    if (m_cell || (m_cell->ParticleCount()==0)) {
        // Uniformly select a particle from the sampling
        // cell.
        int i = m_cell->Particles().Select(rng);
        if (i >= 0) {
            p = m_cell->Particles().At(i)->Clone();
        } else {
            // If sampling failed, then just use the default
            // particle.  This is a serious error, and really
            // should never happen.
            p = m_particle->Clone();
            printf("sweep: Ensemble sampling failed.  Using default particle "
                   "(Sweep, BirthProcess::Perform).");
        }
    } else {
        // Create a copy of the default particle, if there
        // is no sampling cell.
        p = m_particle->Clone();
    }

    // Add the new particle to the ensemble.
    sys.Particles().Add(*p, rng);
    return 0;
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
        const unsigned int trueval  = 1;
        const unsigned int falseval = 0;

        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Serialize base class.
        Process::Serialize(out);

        // Write default particle.
        if (m_particle) {
            out.write((char*)&trueval, sizeof(trueval));
            m_particle->Serialize(out);
        } else {
            out.write((char*)&falseval, sizeof(falseval));
        }
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, BirthProcess::Serialize).");
    }
}

// Reads the object from a binary stream.
void BirthProcess::Deserialize(std::istream &in, const Sweep::Mechanism &mech)
{
    delete m_particle; m_particle = NULL;

    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        unsigned int n = 0;

        switch (version) {
            case 0:
                // Deserialize base class.
                Process::Deserialize(in, mech);

                // Read default particle.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                if (n==1) {
                    m_particle = new Particle(in, mech);
                }
                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, BirthProcess::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, BirthProcess::Deserialize).");
    }
}
