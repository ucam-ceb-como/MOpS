/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the Coagulation class declared in the
    swp_coagulation.h header file.

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

#include "swp_PAH_primary.h"
#include "swp_coagulation.h"
#include "swp_mechanism.h"
#include <stdexcept>
#include <iostream>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_01.hpp>

using namespace Sweep;
using namespace Sweep::Processes;
using namespace std;

// Base coagulation class.

/*!
 * @param[in]   mech    Mechanism to which this process should look for services like LPDA
 *
 * Default rate scaling to 1 for backwards compatibility
 */
 Coagulation::Coagulation(const Sweep::Mechanism &mech)
 : Process(mech)
 , mPositionChoice(NoPositionChoice)
 {}

/*!
 * @param[in]       t       Time at which rates are to be calculated
 * @param[in]       sys     System for which rates are to be calculated
 * @param[in]       local_geom  Spatial configuration information
 * @param[in]       coags   Coagulation processes defining the rates
 * @param[in,out]   rates   Vector to which to add the rates of the individual coagulations
 * @param[in]       start   Position in rates at which to start inserting rates
 *
 * @return      Total rate of all the coagulation processes
 */
double Coagulation::CalcRates(double t, const Cell &sys, const Geometry::LocalGeometry1d &local_geom,
                            const CoagPtrVector &coags, fvector &rates, unsigned int start)
{
    // Iterators for the coagulation processes
    CoagPtrVector::const_iterator itCoag = coags.begin();
    const CoagPtrVector::const_iterator itCoagEnd = coags.end();

    // Iterator for the rate vector
    fvector::iterator it = (rates.begin()+start);

    // Use this variable to accumulate the overall sum of the rates
    double sum = 0.0;
    while(itCoag != itCoagEnd) {
        // Store the rate and move on to the next coagulation process
        *it = (*itCoag++)->Rate(t, sys, local_geom);

        // Add the rate to the sum and move the rates vector iterator onto the next position
        sum += *it++;
    }
    return sum;
}

/*!
 * @param[in]       t       Time at which rates are to be calculated
 * @param[in]       sys     System for which rates are to be calculated
 * @param[in]       local_geom  Spatial configuration information
 * @param[in]       coags   Coagulation processes defining the rates
 * @param[in,out]   iterm   Iterator to point at which to put rate terms of the individual coagulations
 *
 * @return      Total rate of all the coagulation processes
 */
double Coagulation::CalcRateTerms(double t, const Cell &sys, const Geometry::LocalGeometry1d &local_geom,
                                const CoagPtrVector &coags, fvector::iterator &iterm) {
    // Iterators for the coagulation processes
    CoagPtrVector::const_iterator itCoag = coags.begin();
    const CoagPtrVector::const_iterator itCoagEnd = coags.end();

    // Use this variable to accumulate the overall sum of the rates
    double sum = 0.0;
    while(itCoag != itCoagEnd) {
        // The next line does three things, effectively in the following order
        // i) Calls RateTerms on *itCoag
        // ii) Advances itCoag
        // iii) Adds the return value of RateTerms to sum
        // Note that the order of 2 and 3 could be reversed without having any effect
        sum += (*itCoag++)->RateTerms(t, sys, local_geom, iterm);
    }
    return sum;
 }

/*!
 *@param[in]        t           Time at which coagulation is being performed
 *@param[in]        ip1         Index of first particle in ensemble
 *@param[in,out]    sp1         Pointer to first particle
 *@param[in]        ip2         Index of second particle in ensemble
 *@param[in,out]    sp2         Pointer to second particle
 *@param[in]        sys         Cell containing particles that are coagulating
 *@param[in,out]    rng         Random number generator
 *
 *@return       Index of new, larger particle
 */
int Coagulation::JoinParticles(const double t, const int ip1, Particle *sp1,
                               const int ip2, Particle *sp2,
                               Cell &sys, rng_type &rng) const {

    // Position for particle after coagulation, default is to take whatever happens
    // to be in sp1
    double newPos = sp1->getPosition();
    double newPosTime = sp1->getPositionTime();
    if(mPositionChoice == UniformPositionChoice) {
        // Change to position of sp2 with prob 0.5 (bernoulli distribution defaults to prob 0.5)
        boost::bernoulli_distribution<> bernoulliDistrib;
        boost::variate_generator<rng_type &, boost::bernoulli_distribution<> > positionChooser(rng, bernoulliDistrib);
        if(positionChooser()) {
            newPos = sp2->getPosition();
            newPosTime = sp2->getPositionTime();
        }
    }
    else if (mPositionChoice == MassPositionChoice) {
        // Change to position of sp2 with prob sp2->Mass()/(sp1->Mass() + sp2->Mass())
        boost::bernoulli_distribution<double> bernoulliDistrib(sp2->Mass()/(sp1->Mass() + sp2->Mass()));
        boost::variate_generator<rng_type &, boost::bernoulli_distribution<double> > positionChooser(rng, bernoulliDistrib);

        if(positionChooser()) {
            newPos = sp2->getPosition();
            newPosTime = sp2->getPositionTime();
        }
    }
    else if (mPositionChoice == LargestMassPositionChoice) {
        // Change to position of particle with largest mass (default to first particle if masses equal)
        if(sp1->Mass() < sp2->Mass()) {
            newPos = sp2->getPosition();
            newPosTime = sp2->getPositionTime();
        }
        else {
            newPos = sp1->getPosition();
            newPosTime = sp1->getPositionTime();
        }
    }
    else if (mPositionChoice == MidpointPositionChoice) {
        // Pick the half way point
        newPos = (sp1->getPosition() + sp2->getPosition()) / 2;
        newPosTime = (sp1->getPositionTime() + sp2->getPositionTime()) / 2;
    }
    else if (mPositionChoice == CentreOfMassPositionChoice) {
        const double m1 = sp1->Mass();
        const double m2 = sp2->Mass();
        newPos = (m1 * sp1->getPosition() + m2 * sp2->getPosition()) / (m1 + m2);
        newPosTime = (m1 * sp1->getPositionTime() + m2 * sp2->getPositionTime()) / (m1 + m2);
    }


    //sys.Particles().SetNumOfInceptedPAH(-1,sp1->Primary());

    // Add contents of particle 2 onto particle 1
    sp1->Coagulate(*sp2, rng);
    sp1->setPositionAndTime(newPos, newPosTime);
    sp1->SetTime(t);
    sp1->incrementCoagCount();

	//if one of the coagulating particles is tracked then we wish to keep tracking it 
	if (sys.Particles().TrackedParticleNumber() > 0){
		sys.Particles().UpdateTracking(ip2, ip1);
	}

    // Tell the ensemble that particle 1 has changed
    sys.Particles().Update(ip1);
    // Particle 2 is now part of particle 1

    // For hybrid particle model
    // if this particle was introduced from the particle-number model, 
    // it does not exist in the ensemble
    if (!(m_mech->IsHybrid() && ip2 == -2)) 
        sys.Particles().Remove(ip2, true);
    return ip1;
}

/*!
 *@param[in]        t           Time at which coagulation is being performed
 *@param[in]        prop1       Rule for choosing first particle
 *@param[in]        prop2       Rule for choosing second particle
 *@param[in]        weight_rule Specify how to combine particle weights
 *@param[in,out]    sys         Cell containing particles that are coagulating
 *@param[in,out]    rng         Random number generator
 *@param[in]        maj         Specify which majorant to use
 *
 *@return       Negative on failure, 0 on success
 *
 * Weighted coagulation is not symmetric, nothing happens to the second particle,
 * it simply defined a size increment for the first particle.
 */
int Coagulation::WeightedPerform(const double t, const Sweep::PropID prop1,
                                 const Sweep::PropID prop2,
                                 const Sweep::Processes::CoagWeightRule weight_rule,
                                 Cell &sys, rng_type &rng,
                                 MajorantType maj) const {

	PartPtrVector dummy;

    int ip1 = sys.Particles().Select(prop1, rng);
    int ip2 = sys.Particles().Select(prop2, rng);

    // Choose and get first particle, then update it.
    Particle *sp1=NULL;
    if (ip1 >= 0) {
        sp1 = sys.Particles().At(ip1);
    } else {
        // Failed to choose a particle.
        return -1;
    }

    // Choose and get unique second particle, then update it.  Note, we are allowed to do
    // this even if the first particle was invalidated.
    unsigned int guard = 0;
    while ((ip2 == ip1) && (++guard<1000))
            ip2 = sys.Particles().Select(prop2, rng);

    Particle *sp2=NULL;
    if ((ip2>=0) && (ip2!=ip1)) {
        sp2 = sys.Particles().At(ip2);
    } else {
        // Failed to select a unique particle.
        return -1;
    }

    //Calculate the majorant rate before updating the particles
    const double majk = MajorantKernel(*sp1, *sp2, sys, maj);

    //Update the particles
	m_mech->UpdateParticle(*sp1, sys, t, ip1, rng, dummy);
    // Check that particle is still valid.  If not,
    // remove it and cease coagulating.
    if (!sp1->IsValid()) {
        // Must remove first particle now.
        sys.Particles().Remove(ip1);

        // Invalidating the index tells this routine not to perform coagulation.
        ip1 = -1;
        return 0;
    }

	m_mech->UpdateParticle(*sp2, sys, t, ip2, rng, dummy);
    // Check validity of particles after update.
    if (!sp2->IsValid()) {
        // Tell the ensemble to update particle one before we confuse things
        // by removing particle 2
        sys.Particles().Update(ip1);

        // Must remove second particle now.
        sys.Particles().Remove(ip2);

        // Invalidating the index tells this routine not to perform coagulation.
        ip2 = -1;

        return 0;
    }

    // Check that both the particles are still valid.
    if ((ip1>=0) && (ip2>=0)) {
        // Must check for ficticious event now by comparing the original
        // majorant rate and the current (after updates) true rate.

        double truek = CoagKernel(*sp1, *sp2, sys);
        double ceff=0;
        if (majk<truek)
            std::cout << "maj< true"<< std::endl;

		//added by ms785 to include the collision efficiency in the calculation of the rate
		if (sys.ParticleModel()->AggModel()==AggModels::PAH_KMC_ID)
		{
			 ceff=sys.ParticleModel()->CollisionEff(sp1,sp2);
			 truek*=ceff;
		}

        if (!Fictitious(majk, truek, rng)) {
            //Adjust the statistical weight
            switch(weight_rule) {
            case Sweep::Processes::CoagWeightHarmonic :
                sp1->setStatisticalWeight(1.0 / (1.0 / sp1->getStatisticalWeight() +
                                                 1.0 / sp2->getStatisticalWeight()));
                break;
            case Sweep::Processes::CoagWeightHalf :
                sp1->setStatisticalWeight(0.5 * sp1->getStatisticalWeight());
                break;
            case Sweep::Processes::CoagWeightMass :
                sp1->setStatisticalWeight(sp1->getStatisticalWeight() * sp1->Mass() /
                                          (sp1->Mass() + sp2->Mass()));
                break;
            case Sweep::Processes::CoagWeightRule4 : {
                // This is an arbitrary weighting for illustrative purposes
                const double x1 = sp1->Mass() / std::sqrt(sp1->getStatisticalWeight());
                const double x2 = sp1->Mass() / std::sqrt(sp2->getStatisticalWeight());
                sp1->setStatisticalWeight(sp1->getStatisticalWeight() * x1 /
                                          (x1 + x2));
                break;
                }
            default:
                throw std::logic_error("Unrecognised weight rule (Coagulation::WeightedPerform)");
            }

            // In the weighted particle method the contents of particle 2
            // is added to particle 1 while particle 2 is left unchanged.
            // If particle 1 is a single primary particle made up one PAH:
            // the incepted PAH, after the coagulate event there is then
            // one less incepted PAH. This is what the check does. If
            // particle 1 is a PAH other than the incepted PAH, there is
            // not a need to made an adjustment to the number of incepted PAHs.
			//sys.Particles().SetNumOfInceptedPAH(-1,sp1->Primary());

            // Add contents of particle 2 onto particle 1
            sp1->Coagulate(*sp2, rng);
            sp1->SetTime(t);
            sp1->incrementCoagCount();

            assert(sp1->IsValid());
            // Tell the ensemble that particles 1 and 2 have changed
            sys.Particles().Update(ip1);
            sys.Particles().Update(ip2);
        } else {
            sys.Particles().Update(ip1);
            sys.Particles().Update(ip2);
            return 1; // Ficticious event.
        }
    } else {
        // One or both particles were invalidated on update,
        // but that's not a problem.  Information on the update
        // of valid particles must be propagated into the binary
        // tree
        if(ip1 >= 0)
            sys.Particles().Update(ip1);

        if(ip2 >= 0)
            sys.Particles().Update(ip2);
    }

    return 0;
}

// Writes the object to a binary stream.
void Coagulation::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Output position choice rule
        out.write(reinterpret_cast<const char*>(&mPositionChoice), sizeof(mPositionChoice));

        // Serialize base class.
        Process::Serialize(out);

    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, Coagulation::Serialize).");
    }
}

// Reads the object from a binary stream.
void Coagulation::Deserialize(std::istream &in, const Sweep::Mechanism &mech)
{
    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        switch (version) {
            case 0:
                in.read(reinterpret_cast<char*>(&mPositionChoice), sizeof(mPositionChoice));

                // Deserialize base class.
                Process::Deserialize(in, mech);

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, Coagulation::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, Coagulation::Deserialize).");
    }
}


//==============================================================
