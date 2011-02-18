/*!
 * \file   swp_secondary_primary_coag.cpp
 * \author Robert I A Patterson
 *  Copyright (C) 2010 Robert I A Patterson.
 *
 *  Project:        sweepc (population balance solver)
 *  Sourceforge:    http://sourceforge.net/projects/mopssuite
 *
 * \brief  Process for free molecular coagulation of small (secondary) particles with particles from the main population

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
#include "swp_secondary_primary_coag.h"

#include "swp_cell.h"
#include "swp_mechanism.h"

#include <sstream>


const Sweep::real Sweep::Processes::SecondaryPrimaryCoag::m_efm = 2.2;

/*!
 * @param[in]   mech    Mechanism to which this process should look for services like LPDA
 */
Sweep::Processes::SecondaryPrimaryCoag::SecondaryPrimaryCoag(const Sweep::Mechanism &mech)
: Coagulation(mech)
{
    m_name = "SecondaryPrimaryCoag";
}

/*!
 * @param[in,out]   in      Binary stream from which to read details of this coagulation process
 * @param[in]       mech    Mechanism to which this process should look for services like LPDA
 */
Sweep::Processes::SecondaryPrimaryCoag::SecondaryPrimaryCoag(std::istream &in, const Sweep::Mechanism &mech)
: Coagulation(mech)
{
    m_name = "SecondaryPrimaryCoag";

    // No data in SecondaryPrimaryCoag so call Coagulation::Deserialise
    // to deal with the data there.
    Deserialize(in, mech);
}

/*!
 * \return      Enum identifying this process as coagulation of secondary particles with particles from the main population
 */
Sweep::Processes::ProcessType Sweep::Processes::SecondaryPrimaryCoag::ID() const {
    return Secondary_Primary_Coagulation_ID;
}

/*!
 * The (majorant) rate of this coagulation process is expressed as
 * a sum of several terms that are of the form \f$a_i b_j\f$ where
 * \f$a_i\f$ is a property of the first particle and \f$b_j\f$ is
 * a property of the second particle.
 *
 * @return      Number of contributions to the process rate
 */
unsigned int Sweep::Processes::SecondaryPrimaryCoag::TermCount() const {
    return mTermCount;
}

/*!
 * Rate of the process for the given system, scaled up by the majorant factor
 * to allow for fictitious events.
 *
 *\param[in]        t                   Time at which rate should be calculated
 *\param[in]        sys                 System for which to calculate the process rate
 *\param[in]        local_geom          Details of geometry around current location
 *
 *\return           Rate of this process
 */
Sweep::real Sweep::Processes::SecondaryPrimaryCoag::Rate(real t, const Cell &sys) const {
    // Create a vector so that we can call RateTerms
    fvector terms(mTermCount);
    fvector::iterator itTerms = terms.begin();
    return RateTerms(t, sys, itTerms);
}

/*!
 * The free molecular coagulation of secondary particles is implemented with
 * three rate terms, but (apart from a factor of 2) the same majorant is used
 * for all three terms.
 *
 *\param[in]        t                   Time at which rate should be calculated
 *\param[in]        sys                 System for which to calculate the rate
 *\param[in]        local_geom          Details of geometry around current location
 *\param[in,out]    iterm               Position in vector at which to insert the rate
 *
 *\return       Total rate of process
 */
 Sweep::real Sweep::Processes::SecondaryPrimaryCoag::RateTerms(real t, const Cell &sys,
                                                            //const Geometry::LocalGeometry1d& local_geom,
                                                            fvector::iterator &iterm) const {

     // Find the larger of the statistical weight of particles in the secondary population and of
     // the particles in the main population.
     const real maxWeight = 1.0 / std::min(sys.SecondarySampleVolume(), sys.SampleVolume());

     // Details of the secondary population
     const unsigned int numSecondaryParticles = sys.Particles().SecondaryCount();
     const real maxCollDiam2 = m_mech->maxSecondaryCollDiam() * m_mech->maxSecondaryCollDiam();
     const real invRootMinMass = 1.0 / std::sqrt(m_mech->minSecondaryMass());

     // Common factor in all terms of the rate expression
     const real r = 2 * m_efm * CFM * std::sqrt(sys.Temperature()) * A() * maxWeight;

     // Use this variable to accumulate the total rate
     real sum = 0;

     // Calculate the first rate term
     *iterm = r * numSecondaryParticles * numSecondaryParticles * maxCollDiam2 * invRootMinMass;
     // Add it to the total and move the iterator on ready to receive the next rate term value
     sum += *iterm++;

     // Calculate the second rate term
     *iterm = r * numSecondaryParticles * invRootMinMass * sys.Particles().GetSum(TreeCache::iD2);
     // Add it to the total and move the iterator on ready to receive the next rate term value
     sum += *iterm++;

     // Calculate the third rate term
     *iterm = r * numSecondaryParticles * maxCollDiam2 * sys.Particles().GetSum(TreeCache::iM_1_2);
     // Add it to the total and move the iterator on ready to receive the next rate term value
     sum += *iterm++;

     // Calculate the fourth rate term
     *iterm = r * numSecondaryParticles *  sys.Particles().GetSum(TreeCache::iD2_M_1_2);
     // Add it to the total and move the iterator on ready to receive the next rate term value
     sum += *iterm++;

     return sum;
 }


/*!
 * \param[in]       t           Time
 * \param[in,out]   sys         System to update
 * \param[in]       local_geom  Details of local phsyical layout
 * \param[in]       iterm       Index of rate term responsible for this event
 * \param[in,out]   rand_int    Pointer to function that generates uniform integers on a range
 * \param[in,out]   rand_u01    Pointer to function that generates U[0,1] deviates
 * \param[out]      out         Details of any particle being transported out of system
 *
 * \return      0 on success, otherwise negative.
 *
 * @exception   std::runtime_error  Invalid rate term index
 */
int  Sweep::Processes::SecondaryPrimaryCoag::Perform(Sweep::real t,
                              Sweep::Cell &sys, 
                              const Geometry::LocalGeometry1d& local_geom,
                              unsigned int iterm,
                              int (*rand_int)(int, int), 
                              Sweep::real(*rand_u01)(), 
                              Sweep::Transport::TransportOutflow *out) const
{
    // Choose first coagulation partner from main population
    int index1 = -1;
    switch (iterm) {
        case 0:
            // Select uniformly
            index1 = sys.Particles().Select(rand_int);
            break;
        case 1:
            // Select by collision diameter squared
            index1 = sys.Particles().Select(TreeCache::iD2, rand_int , rand_u01);
            break;
        case 2:
            // Select by 1 over square root of mass
            index1 = sys.Particles().Select(TreeCache::iM_1_2, rand_int , rand_u01);
            break;
        case 3:
            // Select by 1 over square root of mass times collision diameter squared
            index1 = sys.Particles().Select(TreeCache::iD2_M_1_2, rand_int , rand_u01);
            break;
        default:
            std::ostringstream msg;
            msg << "Invalid rate term index (" ; //<< iterm << ") in SecondaryPrimaryCoag::Perform()\n";
            throw std::runtime_error(msg.str());
    }
    Particle* sp1 = sys.Particles().At(static_cast<unsigned int>(index1));

    // Perform any deferred events on the first particle.  This changes the particle
    // in the ensemble, because we have a pointer to the particle, not a local copy.
    m_mech->UpdateParticle(*sp1, sys, t, rand_u01);
    // Check that particle is still valid.  If not,
    // remove it and cease coagulating.
    if (!sp1->IsValid()) {
        // Must remove first particle now.
        sys.Particles().Remove(index1);

        return 0;
    }

    // Choose second coagulation partner uniformly from the secondary population
    const int index2 = sys.Particles().SelectSecondaryParticle(rand_int);

    // Perform any deferred events on the second particle.  This changes the particle
    // in the ensemble, because we have a pointer to the particle, not a local copy.
    Particle* sp2 = sys.Particles().SecondaryParticleAt(static_cast<unsigned int>(index2));
    m_mech->UpdateParticle(*sp2, sys, t, rand_u01);
    // Check that particle is still valid.  If not,
    // remove it and cease coagulating.
    if (!sp2->IsValid()) {
        // Must remove particle now.
        sys.Particles().RemoveSecondaryParticle(index2);

        // Tell the main ensemble that the first particle has been updated
        sys.Particles().Update(index1);

        return 0;
    }

    // Now test if the event is fictitious
    const real majK = MajorantKernel(sys, *sp1);
    const real trueK = FreeMolKernel(*sp1, *sp2, sys.Temperature());
    if(Fictitious(majK, trueK, rand_u01)) {
        //fictitious event

        // Tell the main ensemble that the first particle has been updated
        sys.Particles().Update(index1);

        return 0;
    }
    else {
        // See which population has the highest statistical weight per particle
        const real primaryWeight = 1.0 / sys.SampleVolume();
        const real secondaryWeight = 1.0 / sys.SecondarySampleVolume();
        const bool secondaryWeightGreater = (secondaryWeight > primaryWeight);

        if(secondaryWeightGreater) {
            // Stick the secondary particle onto the particle that is in the main population
            sp1->Coagulate(*sp2, rand_int, rand_u01);

            // The secondary particle is removed with the following probability.  There is a
            // probability that it will not be removed because it represents more physical
            // particles than the particle from the main population and so to get mass conservation
            // in mean requires that the secondary particle sometimes be kept.
            if(rand_u01() < primaryWeight / secondaryWeight)
                sys.Particles().RemoveSecondaryParticle(index2);
        }
        else {
            if(rand_u01() < secondaryWeight / primaryWeight)
                sp1->Coagulate(*sp2, rand_int, rand_u01);

            sys.Particles().RemoveSecondaryParticle(index2);
        }

        // Tell ensemble to catch up with changes to the particle in the main population
        sys.Particles().Update(index1);

        return 0;
    }

    // Should never get here
    return -1;
}

/*!
 * Evaluate the free molecular collision kernel between two particles
 *
 * @param[in]   sp1             First  particle for which to calculate kernel
 * @param[in]   sp2             Second particle for which to calculate kernel
 * @param[in]   temperature     Temperature of system in which coagulation rate is to be calculated
 *
 * @return      Coagulation kernel
 */
Sweep::real Sweep::Processes::SecondaryPrimaryCoag::FreeMolKernel(const Particle &sp1, const Particle &sp2,
                                                        real temperature) const {
    const real dterm = sp1.CollDiameter()+sp2.CollDiameter();
    return m_efm * CFM * A() *
           sqrt(temperature * ((1.0/sp1.Mass())+(1.0/sp2.Mass()))) *
           dterm * dterm;
}

/*!
 * Evaluate the free molecular collision kernel between two particles
 *
 * @param[in]   sp1         Particle from main population for which to calculate majorant kernel
 * @param[in]   sys         System temperature and extremal secondary particle properties
 *
 * @return      Coagulation kernel
 */
Sweep::real Sweep::Processes::SecondaryPrimaryCoag::MajorantKernel(const Cell& sys, const Particle &sp1) const {
    return 2 * m_efm * CFM * std::sqrt(sys.Temperature()) * A()
           * (sp1.InvSqrtMass() + 1.0 / std::sqrt(sys.ParticleModel()->minSecondaryMass()))
           * (sp1.CollDiamSquared() + sys.ParticleModel()->maxSecondaryCollDiam()
                                      * sys.ParticleModel()->maxSecondaryCollDiam());
}
