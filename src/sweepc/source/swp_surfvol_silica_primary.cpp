/*!
 * @file    swp_surfvol_silica_primary.cpp
 * @author  William J Menz
 * @brief   Implementation of SurfVolSilicaPrimary
 *
 *   Author(s):      William J Menz
 *   Project:        sweepc (population balance solver)
 *   Copyright (C) 2012 William J Menz
 *
 *   File purpose:
 *      Implementation of SurfVolSilicaPrimary
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

#include "swp_aggmodel_type.h"
#include "swp_primary.h"
#include "swp_surfvol_primary.h"
#include "swp_surfvol_silica_primary.h"
#include "swp_cell.h"
#include "swp_sprog_idealgas_wrapper.h"

using namespace Sweep::AggModels;

// Private default constructor
SurfVolSilicaPrimary::SurfVolSilicaPrimary() {}

/*!
 * Initialising constructor
 *
 * @param time      Time at which particle is created
 * @param model     Particle model describing particle
 * @return          Initialised particle
 */
SurfVolSilicaPrimary::SurfVolSilicaPrimary(
        double time,
        const Sweep::ParticleModel &model
        )
: SurfVolPrimary(time, model),
  m_sinter_rate(0.0)
{
    if (m_comp.size() < size_t(3u))
        throw(std::runtime_error("SurfVolSilicaPrimary shouldn't be used with"
                " fewer than three components (nSi, nO, nOH)."));
}

/*!
 * Copy constructor
 *
 * @param copy  Particle to be copied
 * @return      Initialised particle
 */
SurfVolSilicaPrimary::SurfVolSilicaPrimary(
        const SurfVolSilicaPrimary &copy
        )
{
    *this = copy;
}

/*!
 * Stream-reading constructor
 *
 * @param in        Input binary stream
 * @param model     Particle model describing particle
 * @return          Initialised particle
 */
SurfVolSilicaPrimary::SurfVolSilicaPrimary(
        std::istream &in,
        const Sweep::ParticleModel &model
        )
{
    SurfVolSilicaPrimary::Deserialize(in, model);
}

// Default destructor
SurfVolSilicaPrimary::~SurfVolSilicaPrimary() {
    Primary::releaseMem();
}

// Assignment operator
SurfVolSilicaPrimary &SurfVolSilicaPrimary::operator=(const SurfVolSilicaPrimary &rhs)
{
    // First copy everything for spherical primaries.
    SurfVolPrimary::operator=(rhs);

    m_sinter_rate = rhs.m_sinter_rate;

    return *this;
}

// Returns a copy of the model data.
SurfVolSilicaPrimary* const SurfVolSilicaPrimary::Clone() const
{
    return new SurfVolSilicaPrimary(*this);
}

// Returns the aggregation model which this primary describes.
Sweep::AggModels::AggModelType
    SurfVolSilicaPrimary::AggID(void) const {return AggModels::SurfVolSilica_ID;}

/*!
 * Sinters a particle for time dt. In a similar vein to the silica model,
 * the gas-phase is adjusted according to the amount of chemical units lost
 * from the particle.
 *
 * @param[in]   dt      Time for which to sinter
 * @param[in]   sys     Environment for particles
 * @param[in]   model   Sintering model to apply
 * @param[in]   rng     Random number generator
 * @param[in]   wt      Statistical weight
 *
 * @exception   std::runtime_error  Sintering gave initial surface < final surface
 * @exception   std::runtime_error  Could not cast gas phase to SprogIdealGasWrapper
 */
void SurfVolSilicaPrimary::Sinter(
        double dt,
        Sweep::Cell &sys,
        const Sweep::Processes::SinteringModel &model,
        Sweep::rng_type &rng,
        double wt
        )
{
    // Store some values needed after sintering
    double initial_surface = m_surf;

    // Pass through to SurfVolPrimary to sinter the particle
    SurfVolPrimary::Sinter(dt, sys, model, rng, wt);
    m_sinter_rate = model.Rate(0.0, sys, *this);

    // Estimate the number of chemical units to be changed
    double n_OH = GetComponent("hydroxide");
    double n_O  = GetComponent("oxygen");
    double d_comp = n_OH * (initial_surface - m_surf) / m_surf;

    if (d_comp < 0.0)
        throw std::runtime_error("Sintering particle gave initial surface < final surface"
                " in SurfVolSilicaPrimary::Sinter");

    // Adjust composition of particles
    SetComponent("hydroxide", std::max(n_OH - d_comp, 0.0));
    SetComponent("oxygen", (n_O + (0.5 * d_comp)));

    // Now update the cache storing particle information
    SurfVolPrimary::UpdateCache();

    // Attempt to cast the gas-phase
    SprogIdealGasWrapper *gasWrapper = dynamic_cast<SprogIdealGasWrapper*>(&sys.GasPhase());
    if(gasWrapper == NULL)
        throw std::runtime_error("Could not cast gas phase to SprogIdealGasWrapper"
                "in SurfVolSilicaPrimary::Sinter");

    // If excecution reaches here, the cast must have been successful
    Sprog::Thermo::IdealGas *gas = gasWrapper->Implementation();

    // Get the existing concentrations
    Sweep::fvector newConcs;
    gas->GetConcs(newConcs);

    // Adjust the gas-phase
    newConcs[Sprog::Species::Find(std::string("H2O"),*(gas->Species()))]
             += 0.5 * (n_OH - GetComponent("hydroxide")) * wt / (NA * sys.SampleVolume());
    gas->SetConcs(newConcs);
}

// Updates properties of the particle
void SurfVolSilicaPrimary::UpdateCache()
{
    // Use parent class to update particle properties
    SurfVolPrimary::UpdateCache();

    // Estimate diameter using the same formula as SilicaPrimary
    m_dcol = double(PP_Diameter()) * pow(PP_Count(), 1.0/m_pmodel->GetFractDim());
    m_dmob = m_dcol;
}

/*!
 * Get the number of units of a component
 *
 * @param name  Name of component
 * @return      Number of units of a component
 */
double SurfVolSilicaPrimary::GetComponent(std::string name) const
{
    try {
        return m_comp[m_pmodel->ComponentIndex(name)];
    } catch(std::exception& e) {
        throw e.what();
    }
}

/*!
 * Set the number of units of a component
 *
 * @param name  Name of component
 * @param val   Value to set
 */
void SurfVolSilicaPrimary::SetComponent(std::string name, double val)
{
    try {
        m_comp[m_pmodel->ComponentIndex(name)] = val;
    } catch(std::exception& e) {
        throw e.what();
    }
}

void SurfVolSilicaPrimary::Serialize(std::ostream &out) const
{
    if (out.good()) {

        // Output base class.
        SurfVolPrimary::Serialize(out);

        out.write((char*)&m_sinter_rate, sizeof(m_sinter_rate));
    } else {
        throw std::invalid_argument("Output stream not ready "
                               "(Sweep, SurfVolSilicaPrimary::Serialize).");
    }
}

void SurfVolSilicaPrimary::Deserialize(std::istream &in, const Sweep::ParticleModel &model)
{
    if (in.good()) {
        // Read base class.
        SurfVolPrimary::Deserialize(in, model);

        double val(0.0);
        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_sinter_rate = (double)val;
    } else {
        throw std::invalid_argument("Input stream not ready "
                               "(Sweep, SurfVolSilicaPrimary::Deserialize).");
    }
}

