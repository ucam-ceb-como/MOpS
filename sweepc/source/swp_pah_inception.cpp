/*
  Author(s):      Robert Patterson and Markus Sander
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the PAHInception class declared in the
    swp_PAHInception.h header file.

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
    Prof Markus Kraft
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

#include "swp_pah_inception.h"
#include "swp_mechanism.h"

#include "local_geometry1d.h"
using namespace Sweep;
using namespace Sweep::Processes;
using namespace std;


// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (protected).
PAHInception::PAHInception(void)
: Inception()
{
    m_name = "PAHInception";
}

// Initialising constructor.
PAHInception::PAHInception(const Sweep::Mechanism &mech)
: Inception(mech)
{
    m_name = "PAHInception";
}

// Copy constructor.
PAHInception::PAHInception(const PAHInception &copy)
{
    *this = copy;
}

// Stream-reading constructor.
PAHInception::PAHInception(std::istream &in, const Sweep::Mechanism &mech)
{
    Deserialize(in, mech);
}

// Default destructor.
PAHInception::~PAHInception(void)
{
}

// OPERATOR OVERLOADS.

// Assignment operator.
PAHInception &PAHInception::operator =(const PAHInception &rhs)
{
    if (this != &rhs) {
        Inception::operator =(rhs);

    }
    return *this;
}


/*!
 * Create a new particle and add it to the ensemble with position uniformly
 * distributed over the grid cell, if it is of positive size.
 *
 * The iterm parameter is included because it will be needed for many process
 * types and this function is meant to have a general signature.
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
int PAHInception::Perform(const real t, Cell &sys,
                          const Geometry::LocalGeometry1d &local_geom,
                          const unsigned int iterm,
                          int (*rand_int)(int, int), 
                          Sweep::real(*rand_u01)(), 
                          Transport::TransportOutflow * out) const {

    Particle *sp = NULL;

    // Get the cell vertices
    fvector vertices = local_geom.cellVertices();

    // Sample a uniformly distributed position, note that this method
    // works whether the vertices come in increasing or decreasing order,
    // but 1d is assumed for now.
    real posn = vertices.front();

    const real width = vertices.back() - posn;

    if(width > 0) {
        // There is some real spatial detail
        posn += width * rand_u01();
        sp = m_mech->CreateParticle(t, posn, rand_int);
    }
    else {
        // Ignore all questions of position
        sp = m_mech->CreateParticle(t, rand_int);
    }

    sp->UpdateCache();

    if(UseSecondary()) {
        if(m_mech->isSecondary(*sp)) {
            // Particle can go in the secondary population
            sys.Particles().AddSecondaryParticle(*sp, rand_int);
        }
        else {
            // Add to primary population, which require adjusting the statistical weight
            sys.AddParticle(sp, 1.0/sys.SecondarySampleVolume(), rand_int, rand_u01);
        }
    }
    else {
        // Add particle to main ensemble.
        sys.Particles().Add(*sp, rand_int);
    }

    // Update gas-phase chemistry of system.
    adjustGas(sys);

    return 0;
}

// TOTAL RATE CALCULATIONS.

// Returns rate of the process for the given system.
real PAHInception::Rate(real t, const Cell &sys) const
{
    const real rate = NA*sys.PAHFormationRate()*A();

    // PAHFormation rate may be negative which means no inception
    if(rate < 0.0)
        return 0.0;

    return rate * (UseSecondary() ? sys.SecondarySampleVolume() : sys.SampleVolume());
}


// RATE TERM CALCULATIONS.

// Returns the number of rate terms for this process (one).
unsigned int PAHInception::TermCount(void) const {return 1;}

// Calculates the rate terms given an iterator to a real vector. The
// iterator is advanced to the position after the last term for this
// process.  Returns the sum of all terms.
real PAHInception::RateTerms(const real t, const Cell &sys,
                          fvector::iterator &iterm) const
{
    // Calculate the single rate term and advance iterator.
    *iterm = Rate(t, sys);
    return *(iterm++);
}


// READ/WRITE/COPY.

// Creates a copy of the inception.
PAHInception *const PAHInception::Clone(void) const {return new PAHInception(*this);}

// Returns the process type.  Used to identify different
// processes and for serialisation.
ProcessType PAHInception::ID(void) const {return PAH_Inception_ID;}

// Writes the object to a binary stream.
void PAHInception::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Serialize base class.
        Inception::Serialize(out);

    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, PAHInception::Serialize).");
    }
}

// Reads the object from a binary stream.
void PAHInception::Deserialize(std::istream &in, const Sweep::Mechanism &mech)
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
                Inception::Deserialize(in, mech);

                 break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, PAHInception::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, PAHInception::Deserialize).");
    }
}
