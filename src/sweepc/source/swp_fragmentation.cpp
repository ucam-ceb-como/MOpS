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
#include "swp_fragmentation.h"
#include "swp_mechanism.h"
#include <stdexcept>
#include <iostream>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/variate_generator.hpp>

using namespace Sweep;
using namespace Sweep::Processes;
using namespace std;

// Base coagulation class.

/*!
 * @param[in]   mech    Mechanism to which this process should look for services like LPDA
 *
 * Default rate scaling to 1 for backwards compatibility
 */
 Fragmentation::Fragmentation(const Sweep::Mechanism &mech)
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
double Fragmentation::CalcRates(double t, const Cell &sys, const Geometry::LocalGeometry1d &local_geom,
                            const FragPtrVector &frags, fvector &rates, unsigned int start)
{
    // Iterators for the coagulation processes
    FragPtrVector::const_iterator itFrag = frags.begin();
    const FragPtrVector::const_iterator itFragEnd = frags.end();

    // Iterator for the rate vector
    fvector::iterator it = (rates.begin()+start);

    // Use this variable to accumulate the overall sum of the rates
    double sum = 0.0;
    while(itFrag != itFragEnd) {
        // Store the rate and move on to the next coagulation process
        *it = (*itFrag++)->Rate(t, sys, local_geom);

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
double Fragmentation::CalcRateTerms(double t, const Cell &sys, const Geometry::LocalGeometry1d &local_geom,
                                const FragPtrVector &frags, fvector::iterator &iterm) {
    // Iterators for the coagulation processes
    FragPtrVector::const_iterator itFrag = frags.begin();
    const FragPtrVector::const_iterator itFragEnd = frags.end();

    // Use this variable to accumulate the overall sum of the rates
    double sum = 0.0;
    while(itFrag != itFragEnd) {
        // The next line does three things, effectively in the following order
        // i) Calls RateTerms on *itCoag
        // ii) Advances itCoag
        // iii) Adds the return value of RateTerms to sum
        // Note that the order of 2 and 3 could be reversed without having any effect
        sum += (*itFrag++)->RateTerms(t, sys, local_geom, iterm);
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
int Fragmentation::JoinParticles(const double t, const int ip1, Particle *sp1,
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
    sp1->Fragment(*sp2, rng);
    sp1->setPositionAndTime(newPos, newPosTime);
    sp1->SetTime(t);
    sp1->incrementFragCount();

    // Tell the ensemble that particle 1 has changed
    sys.Particles().Update(ip1);
    // Particle 2 is now part of particle 1
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
int Fragmentation::WeightedPerform(const double t, const Sweep::PropID prop,
                                 const Sweep::Processes::FragWeightRule weight_rule,
                                 Cell &sys, rng_type &rng,
                                 MajorantType maj) const {
	PartPtrVector dummy;

    int i = -1;
    i = sys.Particles().Select(m_pid, rng);
    if (i >= 0) {
        Particle *sp = sys.Particles().At(i);
        std::vector<double> m_dcomp(1);
        std::vector<double> m_dvals(1);
        m_dvals[0] = 0.0;
        if (m_mech->AnyDeferred()) {
            m_mech->UpdateParticle(*sp, sys, t, i, rng, dummy);
            if (sp->IsValid()) {
                double majk = MajorantKernel(*sp, sys, Default);
                double truek = FragKernel(*sp, sys);
                if (!Fictitious(majk, truek, rng)) {
                    switch(weight_rule) {
                        case Sweep::Processes::FragWeightSymmetric :
                            sp->setStatisticalWeight(2 * sp->getStatisticalWeight());
                            m_dcomp[0] = - sp->NumCarbon() / 2;
                            break;
                        case Sweep::Processes::FragWeightNumber:
                            sp->setStatisticalWeight(sp->getStatisticalWeight() + 1);
                            m_dcomp[0] = - 1.0;
                            break;
                        case Sweep::Processes::FragWeightMass:
                            sp->setStatisticalWeight(sp->getStatisticalWeight() * sp->NumCarbon() / ( sp->NumCarbon() - 1));
                            m_dcomp[0] = - 1.0;
                            break;
                        default:
                            throw std::logic_error("Unrecognised weight rule (Fragmentation::WeightedPerform)");
                    }
                    sp->Adjust(m_dcomp, m_dvals, rng, 1);
                    sys.Particles().Update(i);
                    sp->SetTime(t);
                    sp->incrementFragCount();
                }
            } else {
                sys.Particles().Remove(i);
            }
        } else {
            switch(weight_rule) {
                case Sweep::Processes::FragWeightSymmetric :
                    sp->setStatisticalWeight(2 * sp->getStatisticalWeight());
                    m_dcomp[0] = - sp->NumCarbon() / 2;
                    break;
                case Sweep::Processes::FragWeightNumber:
                    sp->setStatisticalWeight(sp->getStatisticalWeight() + 1);
                    m_dcomp[0] = - 1.0;
                    break;
                case Sweep::Processes::FragWeightMass:
                    sp->setStatisticalWeight(sp->getStatisticalWeight() * sp->NumCarbon() / ( sp->NumCarbon() - 1));
                    m_dcomp[0] = - 1.0;
                    break;
                default:
                    throw std::logic_error("Unrecognised weight rule (Fragmentation::WeightedPerform)");
            }
            sp->Adjust(m_dcomp, m_dvals, rng, 1);
            sys.Particles().Update(i);
            sp->SetTime(t);
            sp->incrementFragCount();
        }
    } else {
        return -1;
    }
    return 0;
}

// Writes the object to a binary stream.
void Fragmentation::Serialize(std::ostream &out) const
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
                               "(Sweep, Fragmentation::Serialize).");
    }
}

// Reads the object from a binary stream.
void Fragmentation::Deserialize(std::istream &in, const Sweep::Mechanism &mech)
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
                                    "(Sweep, Fragmentation::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, Fragmentation::Deserialize).");
    }
}


//==============================================================
