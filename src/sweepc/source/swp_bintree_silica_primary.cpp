/*!
 * @file    swp_bintree_silica_primary.cpp
 * @author  William J Menz
 * @brief   Implementation of the silica particle model
 *
 *   Author(s):      William J Menz
 *   Project:        sweepc (population balance solver)
 *   Copyright (C) 2012 William J Menz
 *
 *   File purpose:
 *      Implementation of the silica particle model
 *
 *   Licence:
 *      This file is part of "sweepc".
 *
 *      sweepc is free software; you can redistribute it and/or
 *      modify it under the terms of the GNU Lesser General Public License
 *      as published by the Free Software Foundation; either version 2
 *      of the License, or (at your option) any later version.
 *
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU Lesser General Public License for more details.
 *
 *      You should have received a copy of the GNU Lesser General Public
 *      License along with this program; if not, write to the Free Software
 *      Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 *      02111-1307, USA.
 *
 *   Contact:
 *      Prof Markus Kraft
 *      Dept of Chemical Engineering
 *      University of Cambridge
 *      New Museums Site
 *      Pembroke Street
 *      Cambridge
 *      CB2 3RA, UK
 *
 *      Email:       mk306@cam.ac.uk
 *      Website:     http://como.cheng.cam.ac.uk
*/

#include "swp_bintree_silica_primary.h"
#include "swp_bintree_primary.h"
#include "swp_aggmodel_type.h"
#include "swp_primary.h"
#include "swp_cell.h"
#include "swp_sprog_idealgas_wrapper.h"


using namespace Sweep::AggModels;

// Default constructor (private)
BinTreeSilicaPrimary::BinTreeSilicaPrimary() {}

/*!
 * Initialising constructor
 *
 * @param time      Time at which particle is created
 * @param model     Particle model describing particle
 * @return          Initialised particle
 */
BinTreeSilicaPrimary::BinTreeSilicaPrimary(
        const double time,
        const Sweep::ParticleModel &model
        )
: BinTreePrimary(time, model) {
    if (m_comp.size() < size_t(3u))
        throw(std::runtime_error("BinTreeSilicaPrimary shouldn't be used with"
                " fewer than three components (nSi, nO, nOH)."));
}

// Copy constructor
BinTreeSilicaPrimary::BinTreeSilicaPrimary(const BinTreeSilicaPrimary &copy) {
    *this = copy;

    if (copy.m_leftchild != NULL && copy.m_rightchild !=NULL) {
        BinTreePrimary::CopyTree(&copy);
    }
}

// Stream-reading constructor
BinTreeSilicaPrimary::BinTreeSilicaPrimary(
    std::istream &in,
    const Sweep::ParticleModel &model
    ) {
    BinTreePrimary::Deserialize(in, model);
}

// Destructor
BinTreeSilicaPrimary::~BinTreeSilicaPrimary() {
    Primary::releaseMem();
}

// Returns a copy of the primary.
BinTreeSilicaPrimary *const BinTreeSilicaPrimary::Clone(void) const
{
   return new BinTreeSilicaPrimary(*this);
}

// Primary assignment operator
BinTreeSilicaPrimary &BinTreeSilicaPrimary::operator=(const Primary &rhs) {
    operator=(dynamic_cast<const BinTreeSilicaPrimary&>(rhs));

    return *this;
}

// Returns the aggregation model which this primary describes.
Sweep::AggModels::AggModelType BinTreeSilicaPrimary::AggID(void) const
    {return Sweep::AggModels::BinTreeSilica_ID;}

/*!
 * @brief       Sinters particles for time dt
 *
 * This function only operates on non-leaf nodes. It begins at the root
 * node, which sinters for time dt. It then descends the tree to sinter
 * nodes below the root. If the sintering level rises above 95%, Merge
 * is called and the particles are combined.
 *
 * Class functions are explicitly specified to avoid confusion with
 * parent class functions. Once particle is fully sintered, the gas-phase
 * is adjusted according to the change in hydroxide sites.
 *
 * @param[in]   dt      Time for which to sinter
 * @param[in]   sys     Environment for particles
 * @param[in]   model   Sintering model to apply
 * @param[in]   rng     Random number generator
 * @param[in]   wt      Statistical weight
 *
 * @exception   std::runtime_error  Could not cast gas phase to SprogIdealGasWrapper
 */
void BinTreeSilicaPrimary::Sinter(double dt, Sweep::Cell &sys,
                            const Sweep::Processes::SinteringModel &model,
                            Sweep::rng_type &rng,
                            double wt)
{
    // Only update the time on the root node
    if (m_parent == NULL) {
        m_sint_time += dt;
        SetSinteringTime(m_sint_time);
    }

    // For use later, only on root node
    const double n_OH_i = GetComponent("hydroxide");

    // Do only if there is a particle to sinter
    if (m_leftparticle!=NULL && m_rightparticle!=NULL) {

        BinTreeSilicaPrimary::SinterNode(dt, sys, model, rng, wt);

        // Check if the sintering level is above the threshold, and merge
        if(m_children_sintering > 0.95) CheckSintering();

        // Now sinter any children this node has (if not merged after CheckSintering)
        if (m_leftparticle!=NULL && m_rightparticle!=NULL) {

            // Need to cast left and right children as BinTreeSilicaPrimary* types
            // so that the correct Sinter() function is accessed.

            BinTreeSilicaPrimary* lc =
                static_cast<BinTreeSilicaPrimary*>(m_leftchild);
            lc->BinTreeSilicaPrimary::Sinter(dt, sys, model, rng, wt);

            BinTreeSilicaPrimary* rc =
                static_cast<BinTreeSilicaPrimary*>(m_rightchild);
            rc->BinTreeSilicaPrimary::Sinter(dt, sys, model, rng, wt);
        }

        UpdateCache();

        m_children_sintering = SinteringLevel();
    }

    // Update gas-phase once sintering of the whole particle is finished
    if (m_parent == NULL) {
        // Attempt to cast the gas-phase
        SprogIdealGasWrapper *gasWrapper = dynamic_cast<SprogIdealGasWrapper*>(&sys.GasPhase());
        if(gasWrapper == NULL)
            throw std::runtime_error("Could not cast gas phase to SprogIdealGasWrapper"
                    " in BinTreeSilicaPrimary::Sinter");

        // If excecution reaches here, the cast must have been successful
        Sprog::Thermo::IdealGas *gas = gasWrapper->Implementation();

        // Get the existing concentrations
        Sweep::fvector newConcs;
        gas->GetConcs(newConcs);

        // Adjust the gas-phase
        // Use n_OH (old) - n_OH (new) for better mass conservation
        double d_water =
                0.5 * (n_OH_i - GetComponent("hydroxide")) * wt / (NA * sys.SampleVolume());

        newConcs[Sprog::Species::Find(std::string("H2O"),*(gas->Species()))]
                 += d_water;
        gas->SetConcs(newConcs);
    }

}


/*!
 * Sinters the silica node for time dt.
 *
 * A sintering event on a node makes the following change:
 *  + 1 O unit
 *  - 2 OH units
 *
 * The number of units is calculated assuming a constant site density
 * of hydroxide sites, and split evenly across the L/R particles.
 *
 * @param[in]   dt      Time for which to sinter
 * @param[in]   sys     Environment for particles
 * @param[in]   model   Sintering model to apply
 * @param[in]   rng     Random number generator
 * @param[in]   wt      Statistical weight
 *
 * @exception   std::runtime_error  Sintering gave initial surface < final surface
 */
void BinTreeSilicaPrimary::SinterNode(
        double dt,
        Sweep::Cell &sys,
        const Sweep::Processes::SinteringModel &model,
        Sweep::rng_type &rng,
        double wt
        ) {

    if (m_leftparticle!=NULL && m_rightparticle!=NULL) {

    // Store some useful variables
    double initial_surface = m_children_surf;
    double n_OH = GetComponent("hydroxide");

    // Call to parent class to sinter this node
    BinTreePrimary::SinterNode(dt, sys, model, rng, wt);

    // Estimate the number of chemical units to be changed
    // NOTE: different to SurfVolSilicaPrimary as d_comp is split across
    // the two leaf node particles
    double d_comp = 0.5 * n_OH * (initial_surface - m_children_surf) / m_surf;

    if (d_comp < 0.0) {
        // Somehow we got a unphysical sintering process
        throw std::runtime_error("Sintering particle gave initial surface < final surface"
                " in BinTreeSilicaPrimary::Sinter");

    } else if (d_comp > 0.0) {
        // Here there is some component to adjust

        // Adjust composition of particles
        m_leftparticle->SetComponent("hydroxide",
                std::max(m_leftparticle->GetComponent("hydroxide") - d_comp, 0.0));
        m_leftparticle->SetComponent("oxygen",
                (m_leftparticle->GetComponent("oxygen")  + (0.5 * d_comp)));

        m_rightparticle->SetComponent("hydroxide",
                std::max(m_rightparticle->GetComponent("hydroxide") - d_comp, 0.0));
        m_rightparticle->SetComponent("oxygen",
                (m_rightparticle->GetComponent("oxygen") + (0.5 * d_comp)));

        // Now update the cache storing particle information
        UpdateCache();
    }
    }
}

