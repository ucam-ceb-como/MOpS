/*!
 * \file   swp_secondary_freecoag.cpp
 * \author Robert I A Patterson
 *  Copyright (C) 2010 Robert I A Patterson.
 *
 *  Project:        sweepc (population balance solver)
 *  Sourceforge:    http://sourceforge.net/projects/mopssuite
 *
 * \brief  Process for free molecular coagulation of small (secondary) particles

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
#include "swp_secondary_freecoag.h"

#include "swp_cell.h"
#include "swp_mechanism.h"


const Sweep::real Sweep::Processes::SecondaryFreeCoag::m_efm = 2.2;

/*!
 * @param[in]   mech    Mechanism to which this process should look for services like LPDA
 */
Sweep::Processes::SecondaryFreeCoag::SecondaryFreeCoag(const Sweep::Mechanism &mech)
: Coagulation(mech)
{
    m_name = "SecondaryFreeCoag";
}

/*!
 * @param[in,out]   in      Binary stream from which to read details of this coagulation process
 * @param[in]       mech    Mechanism to which this process should look for services like LPDA
 */
Sweep::Processes::SecondaryFreeCoag::SecondaryFreeCoag(std::istream &in, const Sweep::Mechanism &mech)
: Coagulation(mech)
{
    m_name = "SecondaryFreeCoag";

    // No data in SecondaryFreeCoag so call Coagulation::Deserialise
    // to deal with the data there.
    Deserialize(in, mech);
}

/*!
 *
 * \return      Enum identifying this process as free molecular regime coagulation of secondary particles
 */
Sweep::Processes::ProcessType Sweep::Processes::SecondaryFreeCoag::ID() const {
    return Secondary_FreeCoagulation_ID;
}

/*!
 * For this process the majorant can be expressed as one term, using maximum
 * and minimum values of particle properties within the secondary population.
 *
 * @return      Number of contributions to the process rate
 */
unsigned int Sweep::Processes::SecondaryFreeCoag::TermCount() const {
    return mTermCount;
}

/*!
 * Rate of the process for the given system, scaled up by the majorant factor
 * to allow for fictitious events.
 *
 *\param[in]        t                   Time at which rate should be calculated
 *\param[in]        sys                 System for which to calculate the process rate
 *
 *\return           Rate of this process
 */
Sweep::real Sweep::Processes::SecondaryFreeCoag::Rate(real t, const Cell &sys) const {
    return 0.5 * MajorantKernel(sys) / sys.SecondarySampleVolume()
           * sys.Particles().SecondaryCount() * sys.Particles().SecondaryCount();
}

/*!
 *\param[in]        t                   Time at which rate should be calculated
 *\param[in]        sys                 System for which to calculate the rate
 *\param[in,out]    iterm               Position in vector at which to insert the rate
 *
 *\return       Total rate of process
 */
 Sweep::real Sweep::Processes::SecondaryFreeCoag::RateTerms(real t, const Cell &sys,
                                                            fvector::iterator &iterm) const {
     return (*iterm++ = Rate(t, sys));
}


/*!
 * \param[in]       t           Time
 * \param[in,out]   sys         System to update
 * \param[in]       local_geom  Details of local physical layout
 * \param[in]       iterm       Process term responsible for this event
 * \param[in,out]   rand_int    Pointer to function that generates uniform integers on a range
 * \param[in,out]   rand_u01    Pointer to function that generates U[0,1] deviates
 * \param[out]      out         Details of any particle being transported out of system
 *
 * \return      0 on success, otherwise negative.
 */
int  Sweep::Processes::SecondaryFreeCoag::Perform(Sweep::real t,
                              Sweep::Cell &sys, 
                              const Geometry::LocalGeometry1d& local_geom,
                              unsigned int iterm,
                              int (*rand_int)(int, int), 
                              Sweep::real(*rand_u01)(), 
                              Sweep::Transport::TransportOutflow *out) const
{
    // Choose first coagulation partner
    const int index1 = sys.Particles().SelectSecondaryParticle(rand_int);
    Particle* sp1 = sys.Particles().SecondaryParticleAt(static_cast<unsigned int>(index1));

    // Perform any deferred events on the first particle.  This changes the particle
    // in the ensemble, because we have a pointer to the particle, not a local copy.
    m_mech->UpdateParticle(*sp1, sys, t, rand_u01);
    // Check that particle is still valid.  If not,
    // remove it and cease coagulating.
    if (!sp1->IsValid()) {
        // Must remove first particle now.
        sys.Particles().RemoveSecondaryParticle(index1);

        return 0;
    }

    // Choose second coagulation partner
    const int index2 = sys.Particles().SelectSecondaryParticle(rand_int);

    // Nothing to do if the indices are equal since the diagonal terms were
    // included in the rate calculations.
    if(index1 == index2)
        return 0;



    // Perform any deferred events on the second particle.  This changes the particle
    // in the ensemble, because we have a pointer to the particle, not a local copy.
    Particle* sp2 = sys.Particles().SecondaryParticleAt(static_cast<unsigned int>(index2));
    m_mech->UpdateParticle(*sp2, sys, t, rand_u01);
    // Check that particle is still valid.  If not,
    // remove it and cease coagulating.
    if (!sp2->IsValid()) {
        // Must remove particle now.
        sys.Particles().RemoveSecondaryParticle(index2);

        return 0;
    }

    // Now test if the event is fictitious, note that the majorant kernel in fact
    // does not depend on the particles, if it did, it would be necessary to store
    // its value before performing deferred process updates, see for example
    // AdditiveCoagulation::Perform.
    const real majK = MajorantKernel(*sp1, *sp2, sys, Default);
    const real trueK = CoagKernel(*sp1, *sp2, sys);
    if(Fictitious(majK, trueK, rand_u01)) {
        //fictitious event
        return 0;
    }
    else {
        // Perform the coagulation
        sp1->Coagulate(*sp2, rand_int, rand_u01);

        // Move the particle to the main population with probability proportional to
        // the statistical weight ratio
        MoveToMainPopulation(sp1, sys, rand_int, rand_u01);

        // Remove both particles from the secondary population
        sys.Particles().RemoveTwoSecondaryParticles(index1, index2);

        // Memory for second particle no longer required as its contents have now
        // been copied into *sp1.
        delete sp2;

        return 0;
    }

    // Should never get here
    return -1;
}

/*!
 * Move particle from secondary population to the main population.  This method deals
 * with the difference in statistical weights of particles in the two populations by
 * cloning and or only performing the insertion with a probability less than one.
 *
 * Ownership of the particle is taken by this method, which will either delete it
 * or pass it on to the main particle population.
 *
 * \param[in,out]   sp          Particle to move into main population
 * \param[in,out]   sys         Cell inside which the particle is to be moved
 * \param[in,out]   rand_int    Pointer to function that generates uniform integers on a range
 * \param[in,out]   rand_u01    Pointer to function that generates U[0,1] deviates
 */
void Sweep::Processes::SecondaryFreeCoag::MoveToMainPopulation(Particle* sp, Cell &sys,
                                                               int (*rand_int)(int, int), real(*rand_u01)()) const {
    real secondaryWeight = 1.0 / sys.SecondarySampleVolume();

    sys.AddParticle(sp, secondaryWeight, rand_int, rand_u01);
}

/*!
 * Evaluate the free molecular collision kernel between two particles
 *
 * @param[in]   sp1             First  particle for which to calculate kernel
 * @param[in]   sp2             Second particle for which to calculate kernel
 * @param[in]   sys     System temperature and extremal secondary particle properties
 *
 * @return      Coagulation kernel
 */
Sweep::real Sweep::Processes::SecondaryFreeCoag::CoagKernel(const Particle &sp1,
                                                            const Particle &sp2,
                                                            const Cell &sys) const {
    const real dterm = sp1.CollDiameter()+sp2.CollDiameter();
    return m_efm * CFM * A() *
           sqrt(sys.Temperature() * ((1.0/sp1.Mass())+(1.0/sp2.Mass()))) *
           dterm * dterm;
}

/*!
 * Evaluate the majorant kernel between two secondary particles (details
 * of the particles are not actually used, but rather a general upper bound.
 *
 * @param[in]   sp1         First  particle for which to calculate kernel
 * @param[in]   sp2         Second particle for which to calculate kernel
 * @param[in]   sys         System temperature and extremal secondary particle properties
 * @param[in]   maj         Unused flag to indicate which majorant kernel is required
 *
 * @return      Majorant kernel
 */
Sweep::real Sweep::Processes::SecondaryFreeCoag::MajorantKernel(const Particle &sp1,
                                                                const Particle &sp2,
                                                                const Cell& sys,
                                                                const MajorantType maj) const {
    return MajorantKernel(sys);
}

/*!
 * Evaluate the majorant kernel between two secondary particles (details
 * of the particles are not actually used, but rather a general upper bound.
 *
 * @param[in]   sys         System temperature and extremal secondary particle properties
 *
 * @return      Majorant kernel
 */
Sweep::real Sweep::Processes::SecondaryFreeCoag::MajorantKernel(const Cell& sys) const {
    // 4 * sqrt(2) == 5.656854249
    return 5.656854249 * m_efm * CFM * std::sqrt(sys.Temperature()) * A()
           * sys.ParticleModel()->maxSecondaryCollDiam()
           * sys.ParticleModel()->maxSecondaryCollDiam()
           / std::sqrt(sys.ParticleModel()->minSecondaryMass());
}
