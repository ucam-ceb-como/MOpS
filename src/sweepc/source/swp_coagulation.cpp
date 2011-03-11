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

#include "swp_coagulation.h"
#include "swp_mechanism.h"
#include <stdexcept>
#include <iostream>

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
 , m_a(1.0)
 , mPositionChoice(NoPositionChoice)
 {}

/*!
 * @param[in]       t       Time at which rates are to be calculated
 * @param[in]       sys     System for which rates are to be calculated
 * @param[in]       coags   Coagulation processes defining the rates
 * @param[in,out]   rates   Vector to which to add the rates of the individual coagulations
 * @param[in]       start   Position in rates at which to start inserting rates
 *
 * @return      Total rate of all the coagulation processes
 */
real Coagulation::CalcRates(real t, const Cell &sys, const CoagPtrVector &coags,
                            fvector &rates, unsigned int start)
{
    // Iterators for the coagulation processes
    CoagPtrVector::const_iterator itCoag = coags.begin();
    const CoagPtrVector::const_iterator itCoagEnd = coags.end();

    // Iterator for the rate vector
    fvector::iterator it = (rates.begin()+start);

    // Use this variable to accumulate the overall sum of the rates
    real sum = 0.0;
    while(itCoag != itCoagEnd) {
        // Store the rate and move on to the next coagulation process
        *it = (*itCoag++)->Rate(t, sys);

        // Add the rate to the sum and move the rates vector iterator onto the next position
        sum += *it++;
    }
    return sum;
}

/*!
 * @param[in]       t       Time at which rates are to be calculated
 * @param[in]       sys     System for which rates are to be calculated
 * @param[in]       coags   Coagulation processes defining the rates
 * @param[in,out]   iterm   Iterator to point at which to put rate terms of the individual coagulations
 *
 * @return      Total rate of all the coagulation processes
 */
real Coagulation::CalcRateTerms(real t, const Cell &sys, const CoagPtrVector &coags,
                                fvector::iterator &iterm) {
    // Iterators for the coagulation processes
    CoagPtrVector::const_iterator itCoag = coags.begin();
    const CoagPtrVector::const_iterator itCoagEnd = coags.end();

    // Use this variable to accumulate the overall sum of the rates
    real sum = 0.0;
    while(itCoag != itCoagEnd) {
        // The next line does three things, effectively in the following order
        // i) Calls RateTerms on *itCoag
        // ii) Advances itCoag
        // iii) Adds the return value of RateTerms to sum
        // Note that the order of 2 and 3 could be reversed without having any effect
        sum += (*itCoag++)->RateTerms(t, sys, iterm);
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
 *@param[in,out]    rand_int    Pointer to function that generates integers uniformly on an interval
 *@param[in,out]    rand_u01    Pointer to function that generates U[0,1] deviates
 *
 *@return       Index of new, larger particle
 */
int Coagulation::JoinParticles(const real t, const int ip1, Particle *sp1,
                               const int ip2, Particle *sp2,
                               Cell &sys, int (*rand_int)(int, int),
                               real(*rand_u01)()) const {

    // Position for particle after coagulation, default is to take whatever happens
    // to be in sp1
    real newPos = sp1->getPosition();
    real newPosTime = sp1->getPositionTime();
    if(mPositionChoice == UniformPositionChoice) {
        // Change to position of sp2 with prob 0.5
        if(rand_u01() < 0.5) {
            newPos = sp2->getPosition();
            newPosTime = sp2->getPositionTime();
        }
    }
    else if (mPositionChoice == MassPositionChoice) {
        // Change to position of sp2 with prob sp2->Mass()/(sp1->Mass() + sp2->Mass())
        if(rand_u01() * (sp1->Mass() + sp2->Mass()) < sp2->Mass()) {
            newPos = sp2->getPosition();
            newPosTime = sp2->getPositionTime();
        }
    }

    // Add contents of particle 2 onto particle 1
    sp1->Coagulate(*sp2, rand_int, rand_u01);
    sp1->setPositionAndTime(newPos, newPosTime);
    sp1->SetTime(t);
    sp1->incrementCoagCount();

    // Tell the ensemble that particle 1 has changed
    sys.Particles().Update(ip1);
    // Particle 2 is now part of particle 1
    sys.Particles().Remove(ip2, !m_mech->UseSubPartTree());
    return ip1;
}

/*!
 *@param[in]        t           Time at which coagulation is being performed
 *@param[in]        prop1       Rule for choosing first particle
 *@param[in]        prop2       Rule for choosing second particle
 *@param[in]        sys         Cell containing particles that are coagulating
 *@param[in,out]    rand_int    Pointer to function that generates integers uniformly on an interval
 *@param[in,out]    rand_u01    Pointer to function that generates U[0,1] deviates
 *
 *@return       Negative on failure, 0 on success
 *
 * Weighted coagulation is not symmetric, nothing happens to the second particle,
 * it simply defined a size increment for the first particle.
 */
int Coagulation::WeightedPerform(const real t, const Sweep::PropID prop1,
                                 const Sweep::PropID prop2,
                                 const Sweep::Processes::CoagWeightRule weight_rule,
                                 Cell &sys, int (*rand_int)(int, int),
                                 real(*rand_u01)()) const {
    int ip1 = sys.Particles().Select(prop1, rand_int, rand_u01);
    int ip2 = sys.Particles().Select(prop2, rand_int, rand_u01);

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
            ip2 = sys.Particles().Select(prop2, rand_int, rand_u01);

    Particle *sp2=NULL;
    if ((ip2>=0) && (ip2!=ip1)) {
        sp2 = sys.Particles().At(ip2);
    } else {
        // Failed to select a unique particle.
        return -1;
    }

    //Calculate the majorant rate before updating the particles
    const real majk = MajorantKernel(*sp1, *sp2, sys, Default);

    //Update the particles
    m_mech->UpdateParticle(*sp1, sys, t, rand_u01);
    // Check that particle is still valid.  If not,
    // remove it and cease coagulating.
    if (!sp1->IsValid()) {
        // Must remove first particle now.
        sys.Particles().Remove(ip1);

        // Invalidating the index tells this routine not to perform coagulation.
        ip1 = -1;
        return 0;
    }

    m_mech->UpdateParticle(*sp2, sys, t, rand_u01);
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

        real truek = CoagKernel(*sp1, *sp2, sys);

        if (!Fictitious(majk, truek, rand_u01)) {
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
                const real x1 = sp1->Mass() / std::sqrt(sp1->getStatisticalWeight());
                const real x2 = sp1->Mass() / std::sqrt(sp2->getStatisticalWeight());
                sp1->setStatisticalWeight(sp1->getStatisticalWeight() * x1 /
                                          (x1 + x2));
                break;
                }
            default:
                throw std::logic_error("Unrecognised weight rule (Coagulation::WeightedPerform)");
            }

            // Add contents of particle 2 onto particle 1
            sp1->Coagulate(*sp2, rand_int, rand_u01);
            sp1->SetTime(t);
            sp1->incrementCoagCount();

            assert(sp1->IsValid());
            // Tell the ensemble that particles 1 and 2 have changed
            sys.Particles().Update(ip1);
            sys.Particles().Update(ip1);
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

        // Output rate scaling factor
        out.write(reinterpret_cast<const char*>(&m_a), sizeof(m_a));

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
                real a;
                in.read(reinterpret_cast<char*>(&a), sizeof(a));
                SetA(a);

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
