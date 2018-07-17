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

#include <boost/random/uniform_01.hpp>

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
PAHInception::PAHInception(const Sweep::Mechanism &mech,
                          const EnvironmentInterface::PropertyIndex pah_index)
: Inception(mech)
, mPAHFormationIndex(pah_index)
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
 * \param[in,out]   rng         Random number generator
 *
 * \return      0 on success, otherwise negative.
 */
int PAHInception::Perform(const double t, Cell &sys,
                          const Geometry::LocalGeometry1d &local_geom,
                          const unsigned int iterm,
                          rng_type &rng) const {

	Particle *sp = NULL;

	//If using weighted PAHs...
	if (sys.ParticleModel()->Components(0)->WeightedPAHs()){
	//Check to see if a particle that matches that of the incepted PAHs already exists in the emsemble
		int j = sys.Particles().NumOfInceptedPAH(m_mech->AggModel());
		if (j > 0) {
			//There is already an inception PAH in the ensemble (should only be 1). 
			//Just update it's statistical weight
			int Pindex = sys.Particles().IndexOfInceptedPAH(m_mech->AggModel());
			sp = sys.Particles().At(Pindex);
			int StatWeight = sp->getStatisticalWeight();
			sp->setStatisticalWeight(StatWeight + 1.0);
			sys.Particles().Update(Pindex);
			return 0;
		}
	}

	// Get the cell vertices
	fvector vertices = local_geom.cellVertices();

	// Sample a uniformly distributed position, note that this method
	// works whether the vertices come in increasing or decreasing order,
	// but 1d is assumed for now.
	double posn = vertices.front();

	const double width = vertices.back() - posn;

	if (width > 0) {
		boost::uniform_01<rng_type&, double> uniformGenerator(rng);
		// There is some double spatial detail
		posn += width * uniformGenerator();
		sp = m_mech->CreateParticle(t, posn);
	}
	else {
		// Ignore all questions of position
		sp = m_mech->CreateParticle(t);
	}

	sp->UpdateCache();

	// Add particle to main ensemble.
	if (sys.ParticleCount() < sys.Particles().Capacity()){
		sys.Particles().Add(*sp, rng);
	}
	else
	{
		std::cout << "No room in ensemble after PAH inception" << std::endl;
		sys.Particles().Add(*sp, rng);
	}

	// Update gas-phase chemistry of system.
	adjustGas(sys, sp->getStatisticalWeight());


    return 0;
}

/*!
 * Ensure the number of InceptedPAHs in the particle ensemble is consistent with that in the gas phase
 * If the number of InceptedPAHs (i.e. A1, A2 or A4) in the particle ensemble is less than that in the gas phase, new InceptedPAHs will be created (transfer from gas phase) and added to the ensemble. In contrast, if the number of InceptedPAHs in the particle ensemble are more than that in the gas phase, the redundant InceptedPAH will be removed accordingly.
 * This function is only used for PAH-PP model currently
 * The mass form gas phase will be transfered to particle ensemble by using this function
 *
 * \param[in]       i           number of pyrene supposed in the ensemble
 * \param[in]       t           Time
 * \param[in,out]   sys         System to update
 * \param[in]       iterm       Process term responsible for this event
 * \param[in,out]   rng         Random number generator
 *
 * \return      0 on success, otherwise negative.
 */
int PAHInception::AddInceptedPAH(const int i, const double t, Cell &sys,rng_type &rng) const {

    Particle *sp = NULL;

    // return current number of pyrene in the emsemble
    int j = sys.Particles().NumOfInceptedPAH(m_mech->AggModel());

	if (i > j)
	{
		if (sys.ParticleModel()->Components(0)->WeightedPAHs()){
			if (j > 0) {
				//There is already an inception PAH in the ensemble (should only be 1). 
				//Just update it's statistical weight
				int Pindex = sys.Particles().IndexOfInceptedPAH(m_mech->AggModel());
				sp = sys.Particles().At(Pindex);
				int StatWeight = sp->getStatisticalWeight();
				sp->setStatisticalWeight(StatWeight + i - j);
				sys.Particles().Update(Pindex);
			}
			else{
				sp = m_mech->CreateParticle(t);
				sp->setStatisticalWeight(i);
				sp->UpdateCache();
				// Add particle to main ensemble.
				sys.Particles().Add(*sp, rng);
			}
		}
		else{
			while (i > j)
			{

				// Ignore all questions of position
				sp = m_mech->CreateParticle(t);
				sp->UpdateCache();
				// Add particle to main ensemble.
				sys.Particles().Add(*sp, rng);

				j++;
			}
		}
    }
    else if (i<j){
		if (sys.ParticleModel()->Components(0)->WeightedPAHs()){
			int Pindex = sys.Particles().IndexOfInceptedPAH(m_mech->AggModel());
			sp = sys.Particles().At(Pindex);
			int StatWeight = sp->getStatisticalWeight();
			sp->setStatisticalWeight(StatWeight + i - j);
			sys.Particles().Update(Pindex);
		}
		else{
			while (i < j)
			{
				int Pindex = sys.Particles().IndexOfInceptedPAH(m_mech->AggModel());
				if (Pindex < 0)
					throw runtime_error("There are no InceptedPAH in the ensemble, and all the InceptedPAH molecules are consumed due to unknown reason(Mops, Sweep::PAHInception::Perform).");
				sys.Particles().Remove(Pindex);
				//std::cout<<"j-i is "<<j-i<<std::endl;
				j--;
			}
		}
    }
    //i==j do nothong
    else return 0;
    return 0;
}

// TOTAL RATE CALCULATIONS.

// Returns rate of the process for the given system.
double PAHInception::Rate(double t, const Cell &sys,
                        const Geometry::LocalGeometry1d &local_geom) const
{
    const double rate = NA * sys.GasPhase().PropertyValue(mPAHFormationIndex) * A();

    // PAHFormation rate may be negative which means no inception
    if(rate < 0.0)
        return 0.0;

    return rate *  sys.SampleVolume();
}


// RATE TERM CALCULATIONS.

// Returns the number of rate terms for this process (one).
unsigned int PAHInception::TermCount(void) const {return 1;}

// Calculates the rate terms given an iterator to a double vector. The
// iterator is advanced to the position after the last term for this
// process.  Returns the sum of all terms.
double PAHInception::RateTerms(const double t, const Cell &sys,
                             const Geometry::LocalGeometry1d &local_geom,
                             fvector::iterator &iterm) const
{
    // Calculate the single rate term and advance iterator.
    *iterm = Rate(t, sys, local_geom);
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

        out.write(reinterpret_cast<const char*>(&mPAHFormationIndex), sizeof(mPAHFormationIndex));

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

                in.read(reinterpret_cast<char*>(&mPAHFormationIndex), sizeof(mPAHFormationIndex));
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
