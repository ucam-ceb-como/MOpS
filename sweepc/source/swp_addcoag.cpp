/*
  Author(s):      Robert I A Patterson
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2009 Robert I A Patterson.

  File purpose:
    Implementation of additive coagulation kernel

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

#include "swp_addcoag.h"

#include "swp_params.h"
#include "swp_cell.h"
#include "swp_mechanism.h"


const Sweep::real Sweep::Processes::AdditiveCoagulation::s_MajorantFactor = 1.5;

/**
 * Main way of building a new coagulation object
 * @param{in} mech      Mechanism to which coagulation will belong
 *
 */
Sweep::Processes::AdditiveCoagulation::AdditiveCoagulation(const Sweep::Mechanism &mech)
: Coagulation(mech)
{
    m_name = "AdditiveCoagulation";
}

// Stream-reading constructor.
Sweep::Processes::AdditiveCoagulation::AdditiveCoagulation(std::istream &in, const Sweep::Mechanism &mech)
: Coagulation(mech)
{
    m_name = "AdditiveCoagulation";
    Deserialize(in, mech);
}

// Returns the rate of the process for the given system.
Sweep::real Sweep::Processes::AdditiveCoagulation::Rate(real t, const Cell &sys) const
{
    // Get the number of particles in the system.
    unsigned int n = sys.ParticleCount();

    // Check that there are at least 2 particles before calculating rate.
    if (n > 1) {
        // Get system properties required to calculate coagulation rate.
        real massSum = sys.Particles().GetSums().Mass();

        //real vol = sys.SampleVolume());
        return A() * (n - 1) * massSum * s_MajorantFactor / sys.SampleVolume();
    } else {
        return 0.0;
    }
}

/**
 * Number of terms in the expression for the sum of the majorant
 * kernel over all particle pairs.
 */
unsigned int Sweep::Processes::AdditiveCoagulation::TermCount() const {return TYPE_COUNT;}


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
Sweep::real Sweep::Processes::AdditiveCoagulation::RateTerms(real t, const Cell &sys,
                            fvector::iterator &iterm) const
{
    return *(iterm++) = Rate(t, sys);
}

// Performs the process on the given system. Must return 0
// on success, otherwise negative.
int Sweep::Processes::AdditiveCoagulation::Perform(real t, Cell &sys, unsigned int iterm, Transport::TransportOutflow*) const
{
    // debugging variables
    //int numParticles = sys.ParticleCount();
    //double massParticles = sys.Particles().GetSums().Mass();

    // Select properties by which to choose particles.
    // Note we need to choose 2 particles.  One particle must be chosen
    // uniformly and one with probability proportional
    // to particle mass.

    if (sys.ParticleCount() < 2) {
        return 1;
    }

    int ip1=-1, ip2=-1;

    ip1 = sys.Particles().Select(ParticleCache::iM, Sweep::irnd, Sweep::rnd);
    ip2 = sys.Particles().Select(Sweep::irnd);

    // Choose and get first particle, then update it.
    Particle *sp1=NULL;
    Particle *sp1old=NULL;
    if (ip1 >= 0) {
        sp1 = sys.Particles().At(ip1);

        // Create a copy of the particle before updating.
        sp1old = sp1->Clone();

        m_mech->UpdateParticle(*sp1, sys, t);
    } else {
        // Failed to choose a particle.
        return -1;
    }

    // Check that particle is still valid.  If not,
    // remove it and cease coagulating.
    if (!sp1->IsValid()) {
        // Must remove first particle now.
        sys.Particles().Remove(ip1);

        // Invalidating the index tells this routine not to perform coagulation.
        ip1 = -1;
    }

    // Choose and get unique second particle, then update it.  Note, we are allowed to do
    // this even if the first particle was invalidated.
    unsigned int guard = 0;
    while ((ip2 == ip1) && (++guard<1000))
            ip2 = sys.Particles().Select(Sweep::irnd);

    Particle *sp2=NULL, *sp2old=NULL;
    if ((ip2>=0) && (ip2!=ip1)) {
        sp2 = sys.Particles().At(ip2);

        // Create a copy of the particle before updating.
        sp2old = sp2->Clone();

        m_mech->UpdateParticle(*sp2, sys, t);
    } else {
        // Failed to select a unique particle.
        delete sp1old;
        return -1;
    }

    // Check validity of particles after update.
    if (!sp2->IsValid()) {
        // Must remove second particle now.
        sys.Particles().Remove(ip2);

        // Invalidating the index tells this routine not to perform coagulation.
        ip2 = -1;
    }

    // Check that both the particles are still valid.
    if ((ip1>=0) && (ip2>=0)) {
        // Must check for ficticious event now by calculate the original
        // majorant rate and the current (after updates) true rate.
        real majk = CoagKernel(*sp1old, *sp2old, sys, FiftyPercentExtra);
        real truek = CoagKernel(*sp1, *sp2, sys, None);

        //added by ms785 to include the collision efficiency in the calculation of the rate
        if (sys.ParticleModel()->AggModel()==AggModels::PAH_ID)
        {
            double ceff=sys.ParticleModel()->CollisionEff(sp1,sp2);
            truek*=ceff;
        }

        if (!Ficticious(majk, truek)) {
            assert(sys.Particles().Count() <= sys.Particles().Capacity());
            //std::cout << "Coag #" << ++sCoagCount << ' ' << t << ' ' << numParticles
            //          << ' ' << massParticles;
            // We can now coagulate the particles, remember to
            // remove second particle afterwards.
            if (ip1 < ip2) {
                *sp1 += *sp2;
                sp1->SetTime(t);
                //massParticles = sys.Particles().GetSums().Mass();
                sys.Particles().Update(ip1);
                //massParticles = sys.Particles().GetSums().Mass();
                sys.Particles().Remove(ip2, !m_mech->UseSubPartTree());
                //massParticles = sys.Particles().GetSums().Mass();
            } else {
                *sp2 += *sp1;
                sp2->SetTime(t);
                sys.Particles().Update(ip2);
                sys.Particles().Remove(ip1, !m_mech->UseSubPartTree());
            }

            // riap debugging
            //numParticles = sys.ParticleCount();
            //massParticles = sys.Particles().GetSums().Mass();
            //std::cout << " after coag " << numParticles<< ' ' << massParticles << ' '
            //          << sp1old->Mass() << ' ' << sp1->Mass() << ' '
            //          << sp2old->Mass() << ' ' << sp2->Mass() << ' '
            //          << '\n';

        } else {
            sys.Particles().Update(ip1);
            sys.Particles().Update(ip2);
            delete sp1old;
            delete sp2old;
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

    assert(sys.Particles().Count() <= sys.Particles().Capacity());
    delete sp1old; delete sp2old;
    return 0;
}

/**
 * Calculate the coagulation kernel between two particles in a given environment
 * and include support for majorant kernels.  Note that the details of the
 * environment are not currently used.
 *
 *@param[in]    sp1         First particle
 *@param[in]    sp2         Second particle
 *@param[in]    sys         Details of the environment
 *@param[in]    maj         Flag to show which, if any, majorant to use
 *
 *@return       Value of kernel
 */
Sweep::real Sweep::Processes::AdditiveCoagulation::CoagKernel(const Particle &sp1, const Particle &sp2,
                                     const Cell& sys, MajorantType maj) const
{
    real kernelValue = (sp1.Mass() + sp2.Mass()) * A();
    switch (maj) {
        case None:
            return kernelValue;
        case FiftyPercentExtra:
            return kernelValue * s_MajorantFactor;
    }

    // Invalid majorant, return zero.
    return 0.0;
}

