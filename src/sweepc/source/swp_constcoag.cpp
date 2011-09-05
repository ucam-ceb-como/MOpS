/*
  Author(s):      Robert I A Patterson
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2011 Robert I A Patterson.

  File purpose:
    Implementation of constant additive coagulation kernel

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

#include "swp_constcoag.h"

#include "swp_params.h"
#include "swp_cell.h"
#include "swp_mechanism.h"

using namespace Sweep::Processes;

const Sweep::real Sweep::Processes::ConstantCoagulation::s_MajorantFactor = 1.5;

/**
 * Main way of building a new coagulation object
 * @param{in} mech      Mechanism to which coagulation will belong
 *
 */
Sweep::Processes::ConstantCoagulation::ConstantCoagulation(const Sweep::Mechanism &mech)
: Coagulation(mech)
{
    m_name = "ConstantCoagulation";
}

// Stream-reading constructor.
Sweep::Processes::ConstantCoagulation::ConstantCoagulation(std::istream &in, const Sweep::Mechanism &mech)
: Coagulation(mech)
{
    m_name = "ConstantCoagulation";
    Deserialize(in, mech);
}

// Returns the rate of the process for the given system.
Sweep::real Sweep::Processes::ConstantCoagulation::Rate(real t, const Cell &sys) const
{
    // Get the number of particles in the system.
    const unsigned int n = sys.ParticleCount();

    return A() * n * (n - 1) * s_MajorantFactor / sys.SampleVolume() / 2;
}

/**
 * Number of terms in the expression for the sum of the majorant
 * kernel over all particle pairs.
 */
unsigned int Sweep::Processes::ConstantCoagulation::TermCount() const {return TYPE_COUNT;}


/**
 * Calculate the terms in the sum of the majorant kernel over all particle
 * pairs, placing each term in successive positions of the sequence
 * beginning at iterm and return the sum of the terms added to that
 * vector.
 *
 * @param[in] t         Time for which rates are requested
 * @param[in] sys       Details of the particle population and environment
 * @param[inout] iterm  Pointer to start of sequence to hold the rate terms, returned as one past the end.
 */
Sweep::real Sweep::Processes::ConstantCoagulation::RateTerms(real t, const Cell &sys,
                            fvector::iterator &iterm) const
{
    return *(iterm++) = Rate(t, sys);
}

/*!
 * 
 *
 * \param[in]       t           Time
 * \param[in,out]   sys         System to update
 * \param[in]       local_geom  Details of local phsyical layout
 * \param[in]       iterm       Process term responsible for this event
 * \param[in,out]   rand_int    Pointer to function that generates uniform integers on a range
 * \param[in,out]   rand_u01    Pointer to function that generates U[0,1] deviates
 *
 * \return      0 on success, otherwise negative.
 */
int ConstantCoagulation::Perform(Sweep::real t, Sweep::Cell &sys, 
                             const Geometry::LocalGeometry1d& local_geom,
                             unsigned int iterm,
                             int (*rand_int)(int, int), 
                             Sweep::real(*rand_u01)()) const
{
    // Select properties by which to choose particles.
    // Note we need to choose 2 particles.  One particle must be chosen
    // uniformly and one with probability proportional
    // to particle mass.

    if (sys.ParticleCount() < 2) {
        return 1;
    }

    int ip1=-1, ip2=-1;

    ip1 = sys.Particles().Select(rand_int);
    ip2 = sys.Particles().Select(rand_int);

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
            ip2 = sys.Particles().Select(rand_int);

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
            JoinParticles(t, ip1, sp1, ip2, sp2, sys, rand_int, rand_u01);
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

/**
 * Calculate the coagulation kernel between two particles in a given environment.
 * Note that the details of the environment are not currently used.
 *
 *@param[in]    sp1         First particle
 *@param[in]    sp2         Second particle
 *@param[in]    sys         Details of the environment
 *
 *@return       Value of kernel
 */
Sweep::real Sweep::Processes::ConstantCoagulation::CoagKernel(const Particle &sp1,
                                                              const Particle &sp2,
                                                              const Cell& sys) const
{
    return A();
}


/**
 * Calculate the majorant kernel between two particles in a given environment.
 * Note that the details of the environment are not currently used.
 *
 *@param[in]    sp1         First particle
 *@param[in]    sp2         Second particle
 *@param[in]    sys         Details of the environment
 *@param[in]    maj         Unused flag to indicate which majorant kernel is required
 *
 *@return       Value of majorant kernel
 */
Sweep::real Sweep::Processes::ConstantCoagulation::MajorantKernel(const Particle &sp1,
                                                                  const Particle &sp2,
                                                                  const Cell& sys,
                                                                  const MajorantType maj) const
{
    return CoagKernel(sp1, sp2, sys) * s_MajorantFactor;
}

