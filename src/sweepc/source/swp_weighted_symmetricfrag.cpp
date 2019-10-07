/*
  Author(s):      Robert I A Patterson
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2011 Robert I A Patterson.

  File purpose:
    Implementation of constant coagulation kernel for weighted particles

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

#include "swp_weighted_symmetricfrag.h"

#include "swp_params.h"
#include "swp_cell.h"
#include "swp_mechanism.h"
#include "swp_property_indices.h"


const double Sweep::Processes::WeightedSymmetricFragmentation::s_MajorantFactor = 2.0;

/**
 * Main way of building a new coagulation object
 * @param[in]   mech            Mechanism to which coagulation will belong
 * @param[in]   weight_rule     Specify how to calculate statistical weight of newly coagulation particles
 *
 */
Sweep::Processes::WeightedSymmetricFragmentation::WeightedSymmetricFragmentation(
        const Sweep::Mechanism &mech,
        const FragWeightRule weight_rule)
: Fragmentation(mech)
, m_FragWeightRule(weight_rule)
{
    m_name = "WeightedSymmetricFragmentation";
}

/**
 * Load an instance of this process from a binary stream
 *
 * @param[in,out]   in      Input stream
 * @param[in]       mech    Mechanism to which process will belong
 *
 * @exception   runtime_error   Input stream not ready
 */
Sweep::Processes::WeightedSymmetricFragmentation::WeightedSymmetricFragmentation(std::istream &in, const Sweep::Mechanism &mech)
: Fragmentation(mech)
{
    m_name = "WeightedSymmetricFragmentation";
    Deserialize(in, mech);

    if(in.good())
        in.read(reinterpret_cast<char*>(&m_FragWeightRule), sizeof(m_FragWeightRule));
    else
        throw std::runtime_error("Input stream not ready (WeightedSymmetricFragmentation::WeightedSymmetricFragmentation)");
}

/**
 * Calculates sum of the majorant kernel over all particle
 * pairs.
 *
 * @param[in] t         Time for which rates are requested
 * @param[in] sys       Details of the particle population and environment
 * @param[in] local_geom Details of spatial position and boundaries
 *
 * @return      Sum of all rate terms for this process
 */
double Sweep::Processes::WeightedSymmetricFragmentation::Rate(double t, const Cell &sys,
                            const Geometry::LocalGeometry1d &local_geom) const {
    // Create a vector so we can call through to RateTerms
    Sweep::fvector vec(TYPE_COUNT);
    fvector::iterator it = vec.begin();
    return RateTerms(t, sys, local_geom, it);
}

/**
 * Number of terms in the expression for the sum of the majorant
 * kernel over all particle pairs.
 */
unsigned int Sweep::Processes::WeightedSymmetricFragmentation::TermCount() const {return TYPE_COUNT;}

/**
 * Calculate the terms in the sum of the majorant kernel over all particle
 * pairs, placing each term in successive positions of the sequence
 * beginning at iterm and return the sum of the terms added to that
 * vector.
 *
 * @param[in]     t          Time for which rates are requested
 * @param[in]     sys        Details of the particle population and environment
 * @param[in]     local_geom Details of spatial position and boundaries
 * @param[in,out] iterm      Pointer to start of sequence to hold the rate terms, returned as one past the end.
 *
 * @return      Sum of all rate terms for this process
 */
double Sweep::Processes::WeightedSymmetricFragmentation::RateTerms(double t, const Cell &sys,
                            const Geometry::LocalGeometry1d &local_geom,
                            fvector::iterator &iterm) const
{
    // Get system properties required to calculate coagulation rate.
    double rate = sys.Particles().GetSum(m_pid);

    rate *= A();

    if (m_mech->AnyDeferred()) {
        rate *= s_MajorantFactor;
    }
    *iterm++ = rate;
    return rate;
}

/*!
 * Select particles and simulate a (possibly fictitious) coagulation event
 *
 * \param[in]       t           Time
 * \param[in,out]   sys         System to update
 * \param[in]       local_geom  Details of local phsyical layout
 * \param[in]       iterm       Process term responsible for this event
 * \param[in,out]   rng         Random number generator
 *
 * \return      0 on success, otherwise negative.
 * \exception   logic_error     Unrecognised weight rule
 * \exception   logic_error     Unrecognised rate term index (ie iterm value)
 */
int Sweep::Processes::WeightedSymmetricFragmentation::Perform(
        double t, Sweep::Cell &sys,
        const Geometry::LocalGeometry1d& local_geom,
        unsigned int iterm,
        rng_type &rng) const
{
    // Properties to which the probabilities of particle selection will be proportional
    Sweep::PropID prop;
    switch(static_cast<TermType>(iterm)) {
        case First:
            prop = m_pid;
            break;
        default:
            // This could be removed for performance reasons
            throw std::logic_error("Unrecognised term, (Sweep, WeightedSymmetricFragmentation::Perform)");
    }

    return WeightedPerform(t, prop, m_FragWeightRule, sys, rng, Default);
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
double Sweep::Processes::WeightedSymmetricFragmentation::FragKernel(const Particle &sp, const Cell& sys) const
{
    return A() * sp.getStatisticalWeight();
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
double Sweep::Processes::WeightedSymmetricFragmentation::MajorantKernel(const Particle &sp, const Cell& sys, const MajorantType maj) const
{
    return FragKernel(sp, sys) * s_MajorantFactor;
}

/*!
 * @param[in,out]   out     Binary output stream
 *
 * @exception   runtime_error   Output stream not ready
 */
void Sweep::Processes::WeightedSymmetricFragmentation::Serialize(std::ostream &out) const
{
    // Serialise the parent class
    Fragmentation::Serialize(out);

    // Now the data in this class
    if(out.good())
        out.write(reinterpret_cast<const char*>(&m_FragWeightRule), sizeof(m_FragWeightRule));
    else
        throw std::runtime_error("Output stream not ready (WeightedSymmetricFragmentation::Serialize).");
}


