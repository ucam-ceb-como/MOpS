/*
  Author(s):      Robert I A Patterson
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2011 Robert I A Patterson.

  File purpose:
    Implementation of additive coagulation kernel for weighted particles

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

#include "swp_weighted_addcoag.h"

#include "swp_params.h"
#include "swp_cell.h"
#include "swp_mechanism.h"
#include "swp_property_indices.h"


const Sweep::real Sweep::Processes::WeightedAdditiveCoagulation::s_MajorantFactor = 1.5;

/**
 * Main way of building a new coagulation object
 * @param[in]   mech            Mechanism to which coagulation will belong
 * @param[in]   weight_rule     Specify how to calculate statistical weight of newly coagulation particles
 *
 */
Sweep::Processes::WeightedAdditiveCoagulation::WeightedAdditiveCoagulation(
        const Sweep::Mechanism &mech,
        const CoagWeightRule weight_rule)
: Coagulation(mech)
, m_CoagWeightRule(weight_rule)
{
    m_name = "WeightedAdditiveCoagulation";
}

/**
 * Load an instance of this process from a binary stream
 *
 * @param[in,out]   in      Input stream
 * @param[in]       mech    Mechanism to which process will belong
 *
 * @exception   runtime_error   Input stream not ready
 */
Sweep::Processes::WeightedAdditiveCoagulation::WeightedAdditiveCoagulation(std::istream &in, const Sweep::Mechanism &mech)
: Coagulation(mech)
{
    m_name = "WeightedAdditiveCoagulation";
    Deserialize(in, mech);

    if(in.good())
        in.read(reinterpret_cast<char*>(&m_CoagWeightRule), sizeof(m_CoagWeightRule));
    else
        throw std::runtime_error("Input stream not ready (WeightedAdditiveCoagulation::WeightedAdditiveCoagulation)");
}

// Returns the rate of the process for the given system.
Sweep::real Sweep::Processes::WeightedAdditiveCoagulation::Rate(real t, const Cell &sys) const
{
    // Create a vector so we can call through to RateTerms
    Sweep::fvector vec(TYPE_COUNT);
    fvector::iterator it = vec.begin();
    return RateTerms(t, sys, it);
}

/**
 * Number of terms in the expression for the sum of the majorant
 * kernel over all particle pairs.
 */
unsigned int Sweep::Processes::WeightedAdditiveCoagulation::TermCount() const {return TYPE_COUNT;}


/**
 * Calculate the terms in the sum of the majorant kernel over all particle
 * pairs, placing each term in successive positions of the sequence
 * beginning at iterm and return the sum of the terms added to that
 * vector.
 *
 * @param[in] t         Time for which rates are requested
 * @param[in] sys       Details of the particle population and environment
 * @param[inout] iterm  Pointer to start of sequence to hold the rate terms, returned as one past the end.
 *
 * @return      Sum of all rate terms for this process
 */
Sweep::real Sweep::Processes::WeightedAdditiveCoagulation::RateTerms(real t, const Cell &sys,
                            fvector::iterator &iterm) const
{
    const unsigned int n = sys.ParticleCount();

    if(n > 1) {
        // this value is used twice below
        const real sum_iWM = sys.Particles().GetSums().Property(iWM);
        // Calculate the two contributions to the rate
        const real r1 = A() * s_MajorantFactor
                       * (sys.Particles().GetSums().Property(Sweep::iM)
                          * sys.Particles().GetSums().Property(Sweep::iW)
                          - sum_iWM) / sys.SampleVolume();
        const real r2 = A() * (n - 1) * s_MajorantFactor
                            * sum_iWM / sys.SampleVolume();

        // Now deal with the output
        *iterm++ = r1;
        *iterm++ = r2;
        return r1 + r2;
    }
    else {
        *iterm++ = 0.0;
        *iterm++ = 0.0;
        return 0.0;
    }
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
 * \param[out]      out         Details of any particle being transported out of system
 *
 * \return      0 on success, otherwise negative.
 * \exception   logic_error     Unrecognised weight rule
 * \exception   logic_error     Unrecognised rate term index (ie iterm value)
 */
int Sweep::Processes::WeightedAdditiveCoagulation::Perform(
        Sweep::real t, Sweep::Cell &sys,
        const Geometry::LocalGeometry1d& local_geom,
        unsigned int iterm,
        int (*rand_int)(int, int),
        Sweep::real(*rand_u01)(),
        Sweep::Transport::TransportOutflow *out) const
{
    assert(iterm < TYPE_COUNT);

    // Select properties by which to choose particles.
    // Note we need to choose 2 particles.  One particle must be chosen
    // uniformly and one with probability proportional
    // to particle mass.

    if (sys.ParticleCount() < 1) {
        return 1;
    }

    int ip1=-1, ip2=-1;

    // Properties to which the probabilities of particle selection will be proportional
    Sweep::PropID prop1, prop2;
    switch(static_cast<TermType>(iterm)) {
        case FirstByMassSecondByWeight:
            prop1 = Sweep::iM;
            prop2 = Sweep::iW;
            break;
        case FirstUniformlySecondByMassTimesWeight:
            prop1 = iUniform;
            prop2 = iWM;
            break;
        default:
            // This could be removed for performance reasons
            throw std::logic_error("Unrecognised term, (Sweep, WeightedAdditiveCoagulation::Perform)");
    }

    ip1 = sys.Particles().Select(prop1, rand_int, rand_u01);
    ip2 = sys.Particles().Select(prop2, rand_int, rand_u01);

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
    const real majk = CoagKernel(*sp1, *sp2, sys, FiftyPercentExtra);

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

        real truek = CoagKernel(*sp1, *sp2, sys, None);

        if (!Fictitious(majk, truek, rand_u01)) {
            //Adjust the statistical weight
            switch(m_CoagWeightRule) {
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
                throw std::logic_error("Unrecognised weight rule (WeightedAdditiveCoagulation::Perform)");
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
Sweep::real Sweep::Processes::WeightedAdditiveCoagulation::CoagKernel(const Particle &sp1, const Particle &sp2,
                                     const Cell& sys, MajorantType maj) const
{
    real kernelValue = (sp1.Mass() + sp2.Mass()) * A() * sp2.getStatisticalWeight();
    switch (maj) {
        case None:
            return kernelValue;
        case FiftyPercentExtra:
            return kernelValue * s_MajorantFactor;
    }

    // Invalid majorant, return zero.
    return 0.0;
}

/*!
 * @param[in,out]   out     Binary output stream
 *
 * @exception   runtime_error   Output stream not ready
 */
void Sweep::Processes::WeightedAdditiveCoagulation::Serialize(std::ostream &out) const
{
    // Serialise the parent class
    Coagulation::Serialize(out);

    // Now the data in this class
    if(out.good())
        out.write(reinterpret_cast<const char*>(&m_CoagWeightRule), sizeof(m_CoagWeightRule));
    else
        throw std::runtime_error("Output stream not ready (WeightedAdditiveCoagulation::Serialize).");
}


