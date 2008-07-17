/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the DeathProcess class declared in the
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

#include "swp_death_process.h"
#include "swp_mechanism.h"
#include <stdexcept>

using namespace Sweep;
using namespace Sweep::Processes;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (protected).
DeathProcess::DeathProcess(void)
: Process(), m_a(0.0)
{
    m_name = "Death Process";
}

// Initialising constructor.
DeathProcess::DeathProcess(const Sweep::Mechanism &mech)
: Process(mech), m_a(0.0)
{
    m_name = "Death Process";
}

// Copy constructor.
DeathProcess::DeathProcess(const DeathProcess &copy)
{
    *this = copy;
}

// Stream-reading constructor.
DeathProcess::DeathProcess(std::istream &in, const Sweep::Mechanism &mech)
{
    Deserialize(in, mech);
}

// Default destructor.
DeathProcess::~DeathProcess(void)
{
}

// OPERATOR OVERLOADS.

// Assignment operator.
DeathProcess &DeathProcess::operator =(const DeathProcess &rhs)
{
    if (this != &rhs) {
        Process::operator =(rhs);
        m_a = rhs.m_a;
    }
    return *this;
}


// RATE CONSTANT.

// Returns the rate constant.
real DeathProcess::A(void) const {return m_a;}

// Sets the rate constant.
void DeathProcess::SetA(real a) {m_a = a;}


// TOTAL RATE CALCULATIONS.

// Returns rate of the process for the given system.
real DeathProcess::Rate(real t, const Cell &sys) const 
{
    return m_a * sys.ParticleCount();
}

/*
// Calculates the process rate using the given 
// chemical conditions, rather than those conditions in the
// given system.
real DeathProcess::Rate(const real t, const Sprog::Thermo::IdealGas &gas, 
                        const Cell &sys) const
{
    return m_a * sys.ParticleCount();
}
*/

// RATE TERM CALCULATIONS.

// Returns the number of rate terms for this process (one).
unsigned int DeathProcess::TermCount(void) const {return 1;}

// Calculates the rate terms given an iterator to a real vector. The 
// iterator is advanced to the position after the last term for this
// process.  Returns the sum of all terms.
real DeathProcess::RateTerms(const real t, const Cell &sys, 
                             fvector::iterator &iterm) const
{
    *iterm = m_a * sys.ParticleCount();
    return *(iterm++);
}

/*
// Calculates the rate terms given an iterator to a real vector. The 
// iterator is advanced to the position after the last term for this
// process.  The given chemical conditions are used instead of those
// in the given system object.
real DeathProcess::RateTerms(const real t, const Sprog::Thermo::IdealGas &gas,
                             const Cell &sys, fvector::iterator &iterm) const
{
    *iterm = m_a * sys.ParticleCount();
    return *(iterm++);
}
*/

// PERFORMING THE PROCESS.

// Performs the process on the given system.  The responsible rate term is given
// by index.  Returns 0 on success, otherwise negative.
int DeathProcess::Perform(real t, Cell &sys, unsigned int iterm) const 
{
    // Select a particle for deletion.
    int i = sys.Particles().Select();
    // Delete the particle.
    if (i>=0) sys.Particles().Remove(i);
    return 0;
}


// READ/WRITE/COPY.

// Creates a copy of the inception.
DeathProcess *const DeathProcess::Clone(void) const {return new DeathProcess(*this);}

// Returns the process type.  Used to identify different
// processes and for serialisation.
ProcessType DeathProcess::ID(void) const {return Death_ID;}

// Writes the object to a binary stream.
void DeathProcess::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Serialize base class.
        Process::Serialize(out);

        // Write rate constant.
        double v = (double)m_a;
        out.write((char*)&v, sizeof(v));
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, DeathProcess::Serialize).");
    }
}

// Reads the object from a binary stream.
void DeathProcess::Deserialize(std::istream &in, const Sweep::Mechanism &mech)
{
    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        unsigned int n = 0;
        double val = 0.0;

        switch (version) {
            case 0:
                // Deserialize base class.
                Process::Deserialize(in, mech);

                // Read rate constant.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_a = (real)val;
                
                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, DeathProcess::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, DeathProcess::Deserialize).");
    }
}
