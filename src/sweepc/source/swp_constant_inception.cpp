/*!
 * \file   swp_constant_inception.cpp
 * \author Robert I A Patterson
 *
 * \brief  Class for inception at a constant rate
 *
 *  Copyright (C) 2010 Robert I A Patterson.

  Licence:
    This file is part of "sweepc".

    sweepc is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
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

#include "swp_constant_inception.h"
#include "swp_mechanism.h"

#include "local_geometry1d.h"

#include <boost/random/uniform_01.hpp>


// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (private).
Sweep::Processes::ConstantInception::ConstantInception()
: Inception()
{
    m_name = "ConstantInception";
}

// Initialising constructor.
Sweep::Processes::ConstantInception::ConstantInception(const Sweep::Mechanism &mech)
: Inception(mech)
{
    m_name = "ConstantInception";
}

// Copy constructor.
Sweep::Processes::ConstantInception::ConstantInception(const ConstantInception &copy)
{
    *this = copy;
}

// Stream-reading constructor.
Sweep::Processes::ConstantInception::ConstantInception(std::istream &in, const Sweep::Mechanism &mech)
{
    Deserialize(in, mech);
}

// Default destructor.
Sweep::Processes::ConstantInception::~ConstantInception(void)
{
}

// OPERATOR OVERLOADS.

// Assignment operator.
Sweep::Processes::ConstantInception &Sweep::Processes::ConstantInception::operator =(const ConstantInception &rhs)
{
    if (this != &rhs) {
        Inception::operator =(rhs);

    }

    mRate = rhs.mRate;

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
 * \param[in,out]   rng         Random number generator
 *
 * \return      0 on success, otherwise negative.
 */
int Sweep::Processes::ConstantInception::Perform(const real t, Cell &sys,
                          const Geometry::LocalGeometry1d &local_geom,
                          const unsigned int iterm,
                          rng_type &rng) const {

    // Create a new particle of the type specified
    // by the system ensemble.
    Particle * const sp = m_mech->CreateParticle(t);

    // Get the cell vertices
    fvector vertices = local_geom.cellVertices();

    // Sample a uniformly distributed position, note that this method
    // works whether the vertices come in increasing or decreasing order,
    // but 1d is assumed for now.
    real posn = vertices.front();

    const real width = vertices.back() - posn;
    boost::uniform_01<rng_type&, real> unifDistrib(rng);
    posn += width * unifDistrib();

    sp->setPositionAndTime(posn, t);


    // Initialise the new particle.
    sp->Primary()->SetComposition(ParticleComp());
    sp->Primary()->SetValues(ParticleTrackers());
    sp->UpdateCache();

    // Add particle to system's ensemble.
    sys.Particles().Add(*sp, rng);

    // Update gas-phase chemistry of system.
    adjustGas(sys);

    return 0;
}

// TOTAL RATE CALCULATIONS.

/*!
 *@param[in]            t           Time at which rate is being calculated
 *@param[in]            sys         System for which rate is to be calculated
 *@param[in]            local_geom  Spatial configuration information (ignored)
 *
 *@return   Process rate
 */
Sweep::real Sweep::Processes::ConstantInception::Rate(real t, const Cell &sys,
                                                      const Geometry::LocalGeometry1d &local_geom) const
{
    return mRate * A() * sys.SampleVolume();
}


// RATE TERM CALCULATIONS.

// Returns the number of rate terms for this process (one).
unsigned int Sweep::Processes::ConstantInception::TermCount() const {return 1;}

// Calculates the rate terms given an iterator to a real vector. The
// iterator is advanced to the position after the last term for this
// process.  Returns the sum of all terms.
Sweep::real Sweep::Processes::ConstantInception::RateTerms(const real t, const Cell &sys,
                                                           const Geometry::LocalGeometry1d &local_geom,
                                                           fvector::iterator &iterm) const
{
    // Calculate the single rate term and advance iterator.
    *iterm = Rate(t, sys, local_geom);
    return *(iterm++);
}


// READ/WRITE/COPY.

// Creates a copy of the inception.
Sweep::Processes::ConstantInception *const Sweep::Processes::ConstantInception::Clone() const
{
	return new ConstantInception(*this);
}

// Returns the process type.  Used to identify different
// processes and for serialisation.
Sweep::Processes::ProcessType Sweep::Processes::ConstantInception::ID() const
{
	return Constant_Inception_ID;
}

// Writes the object to a binary stream.
void Sweep::Processes::ConstantInception::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write(reinterpret_cast<const char*>(&version), sizeof(version));

        // Serialize base class.
        Inception::Serialize(out);

        // Constant rate value
        out.write(reinterpret_cast<const char*>(&mRate), sizeof(mRate));

    } else {
        throw std::invalid_argument("Output stream not ready (Sweep, Sweep::Processes::ConstantInception::Serialize).");
    }
}

// Reads the object from a binary stream.
void Sweep::Processes::ConstantInception::Deserialize(std::istream &in, const Sweep::Mechanism &mech)
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

                // Constant rate
                in.read(reinterpret_cast<char*>(&mRate), sizeof(mRate));

                 break;
            default:
                throw std::runtime_error("Serialized version number is invalid (Sweep, Sweep::Processes::ConstantInception::Deserialize).");
        }
    } else {
        throw std::invalid_argument("Input stream not ready (Sweep, Sweep::Processes::ConstantInception::Deserialize).");
    }
}
