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
#include "swp_primary.h"

#include "local_geometry1d.h"

#include <boost/random/uniform_01.hpp>
#include <boost/random/lognormal_distribution.hpp>


// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (private).
Sweep::Processes::ConstantInception::ConstantInception()
: Inception()
, mUseFixedPosition(false)
, mFixedPosition(0.0)
{
    m_name = "ConstantInception";
}

/*!
 * Initialising constructor
 *
 * @param[in]       mech        Mechanism of which process will be part
 * @param[in]       rate        Rate of inception jumps
 * @param[in]       locations   Location parameters for lognormal component distributions, ordered to match the components in the mechanism
 * @param[in]       scales      Scale parameters for lognormal component distributions, ordered to match the components in the mechanism
 *
 * @exception       std::runtime_error   Distribution parameters not specified for all components
 */
Sweep::Processes::ConstantInception::ConstantInception(const Sweep::Mechanism &mech, const double rate,
                                                       const std::vector<double>& locations,
                                                       const std::vector<double>& scales)
: Inception(mech)
, mUseFixedPosition(false)
, mFixedPosition(0.0)
, mComponentDistributions(mech.ComponentCount())
{
    if((mech.ComponentCount() != locations.size()) || (mech.ComponentCount() != scales.size()))
            throw std::runtime_error("Details must be supplied for all components in ConstantInception::ConstantInception");

    SetA(rate);
    m_name = "ConstantInception";

    // It is important that the ordering corresponding to the ordering of the components in the mechanism,
    // just having the same number of entries is necessary, but not sufficient.
    for(unsigned i = 0; i < mech.ComponentCount(); ++i) {
            mComponentDistributions[i] = std::make_pair(locations[i], scales[i]);
            SetParticleComp(i, exp(locations[i] + 0.5 * scales[i] * scales[i]));
    }
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

    mUseFixedPosition = rhs.mUseFixedPosition;
    mFixedPosition = rhs.mFixedPosition;
    mComponentDistributions  = rhs.mComponentDistributions;

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
int Sweep::Processes::ConstantInception::Perform(const double t, Cell &sys,
                          const Geometry::LocalGeometry1d &local_geom,
                          const unsigned int iterm,
                          rng_type &rng) const {

    // Position of newly incepted particle
    double posn;

    if(mUseFixedPosition) {
        // If there is a fixed position the rate should only be positive for the cell containing the fixed position
        if(local_geom.isInCell(mFixedPosition))
            posn = mFixedPosition;
        else
            return 0;
    }
    else {
        // Get the cell vertices
        fvector vertices = local_geom.cellVertices();

        // Sample a uniformly distributed position, note that this method
        // works whether the vertices come in increasing or decreasing order,
        // but 1d is assumed for now.
        posn = vertices.front();

        const double width = vertices.back() - posn;
        boost::uniform_01<rng_type&, double> unifDistrib(rng);
        posn += width * unifDistrib();
    }


    // Initialise the new particle.
    std::vector<double> newComposition(mComponentDistributions.size());
    for(unsigned i = 0; i < mComponentDistributions.size(); ++i) {
        // Construct the distribution on the fly.  If this becomes a performance
        // bottleneck some optimisations might be possible, but caching the distribution
        // object in the mechanism is a bad idea, because the mechanism is potentially
        // shared between threads, which is very dangerous for cached data!
        newComposition[i] = boost::random::lognormal_distribution<double>(mComponentDistributions[i].first,
                                                                  mComponentDistributions[i].second)(rng);
    }

    // Check if using hybrid particle-number/particle model
    // If hybrid_flag is active, track the number of particles in a list
    // Otherwise add particle to system's ensemble.
    if (!m_mech->IsHybrid())
    {
        // Create a new particle of the type specified
        // by the system ensemble.
        Particle * sp = m_mech->CreateParticle(t);

        sp->setPositionAndTime(posn, t);

        sp->Primary()->SetComposition(newComposition);

        sp->Primary()->SetValues(ParticleTrackers());
        sp->UpdateCache();

        sys.Particles().Add(*sp, rng);

        // Update gas-phase chemistry of system.
        if (!sys.GetIsAdiabaticFlag())
            adjustGas(sys, sp->getStatisticalWeight());
        // Update gas-phase chemistry and temperature of system.
        else
           adjustParticleTemperature(sys, sp->getStatisticalWeight(), 1, ParticleComp()[0], 1);
    }
    else
    {
        // Adjust particle number properties
        sys.Particles().UpdateNumberAtIndex(ParticleComp()[0], 1);
        sys.Particles().UpdateTotalParticleNumber(1);
        sys.Particles().UpdateTotalsWithIndex(ParticleComp()[0], 1.0);
        // Update gas-phase chemistry of system.
        if (!sys.GetIsAdiabaticFlag())
            adjustGas(sys, 1);
        // Update gas-phase chemistry and temperature of system.
        else
            adjustParticleTemperature(sys, 1, 1, ParticleComp()[0], 1);
    }

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
double Sweep::Processes::ConstantInception::Rate(double t, const Cell &sys,
                                                      const Geometry::LocalGeometry1d &local_geom) const
{
    double rate =  A() * sys.SampleVolume();

    if(mUseFixedPosition) {
        if(local_geom.isInCell(mFixedPosition)) {
            // Dividing by cell volume (which is distinct from the sample volume)
            // seems to be the right thing to do, because the concentration
            // is averaged over the whole cell and thus is proportional to the
            // amount of reactor covered by the cell
            rate /= local_geom.cellVolume();
        }
        else
            rate = 0.0;
    }

    return rate;
}


// RATE TERM CALCULATIONS.

// Returns the number of rate terms for this process (one).
unsigned int Sweep::Processes::ConstantInception::TermCount() const {return 1;}

// Calculates the rate terms given an iterator to a double vector. The
// iterator is advanced to the position after the last term for this
// process.  Returns the sum of all terms.
double Sweep::Processes::ConstantInception::RateTerms(const double t, const Cell &sys,
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

        // Fixed position support
        out.write(reinterpret_cast<const char*>(&mUseFixedPosition), sizeof(mUseFixedPosition));
        out.write(reinterpret_cast<const char*>(&mFixedPosition), sizeof(mFixedPosition));

        // Component distributions
        const unsigned int n = mComponentDistributions.size();
        out.write(reinterpret_cast<const char*>(&n), sizeof(n));
        for(unsigned int i = 0; i < n; ++i) {
            out.write(reinterpret_cast<const char*>(&mComponentDistributions[i].first),  sizeof(mComponentDistributions[i].first));
            out.write(reinterpret_cast<const char*>(&mComponentDistributions[i].second), sizeof(mComponentDistributions[i].second));
        }
    }
    else {
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
            {
                // Deserialize base class.
                Inception::Deserialize(in, mech);

                // Fixed position support
                in.read(reinterpret_cast<char*>(&mUseFixedPosition), sizeof(mUseFixedPosition));
                in.read(reinterpret_cast<char*>(&mFixedPosition), sizeof(mFixedPosition));

                // Component distributions
                unsigned int n = 0;
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                mComponentDistributions.resize(n);
                for(unsigned int i = 0; i < n; ++i) {
                    in.read(reinterpret_cast<char*>(&mComponentDistributions[i].first),  sizeof(mComponentDistributions[i].first));
                    in.read(reinterpret_cast<char*>(&mComponentDistributions[i].second), sizeof(mComponentDistributions[i].second));
                }
                break;
            }
            default:
                throw std::runtime_error("Serialized version number is invalid (Sweep, Sweep::Processes::ConstantInception::Deserialize).");
        }
    } else {
        throw std::invalid_argument("Input stream not ready (Sweep, Sweep::Processes::ConstantInception::Deserialize).");
    }
}
