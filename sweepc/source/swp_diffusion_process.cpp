/*!
 * \file   swp_diffusion_process.cpp
 * \author Robert I A Patterson
 *  Copyright (C) 2009 Robert I A Patterson.
 *
 *  Project:        sweepc (population balance solver)
 *  Sourceforge:    http://sourceforge.net/projects/mopssuite
 *
 * \brief  Process for stochastic realisation of particle diffusion

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
#include "swp_diffusion_process.h"

#include "swp_cell.h"
#include "swp_transport_outflow.h"

#include <cmath>

const Sweep::real Sweep::Processes::DiffusionProcess::s_MajorantFactor = 1.5;

/*!
 * Set the ID number of the particle property to which
 * the rate of this process is proportional.
 * Copied from the surface reaction.
 *
 *\param[in]            id          Particle property to which rate is proportional
 */
void Sweep::Processes::DiffusionProcess::SetPropertyID(unsigned int id) {
    m_pid     = id;
}

/*!
 *
 * \return      Enum identifying this process as particle diffusion
 */
Sweep::Processes::ProcessType Sweep::Processes::DiffusionProcess::ID() const {
    return Diffusion_ID;
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
Sweep::real Sweep::Processes::DiffusionProcess::Rate(real t, const Cell &sys,
                                                     const Geometry::LocalGeometry1d& local_geom) const {
    real rate = 0;

    // Rate of leftwards jumps
    if(!local_geom.zeroGradient(Geometry::left)) {
        const real dx_1 = 1.0 / local_geom.calcSpacing(Geometry::left);
        rate += 0.5 * dx_1 * dx_1;
    }

    // Rate of rightwards jumps
    if(!local_geom.zeroGradient(Geometry::right)) {
        const real dx_1 = 1.0 / local_geom.calcSpacing(Geometry::right);
        rate += 0.5 * dx_1 * dx_1;
    }

    rate *= sys.Particles().GetSum(static_cast<TreeCache::PropID>(m_pid));
    rate *= std::pow(sys.Temperature(), m_TemperatureExponent);
    return  rate * A() * s_MajorantFactor;

}

/*!
 * Calculate the particle dependent part of the process rate
 * for a single particle.  Scaling factors depending on the
 * chemical environment such as temperature, pressure and
 * species concentrations will not be included.
 *
 *@param[in]        t       Time at which rate is to be calculated
 *@param[in]        sys     System containing particle
 *@param[in]        sp      Particle for which to calculate rate
 *
 *@return       Particle dependent part of process rate
 */
Sweep::real Sweep::Processes::DiffusionProcess::Rate(real t, const Cell &sys,
                                                     const Particle &sp) const {
    // Static cast is copied from Sweep::Ensemble::Select
    return sp.Property(static_cast<ParticleCache::PropID>(m_pid));
}


/*!
 * Append the rate to a vector
 *
 *\param[in]        t                   Time at which rate should be calculated
 *\param[in]        sys                 System for which to calculate the rate
 *\param[in]        local_geom          Details of geometry around current location
 *\param[in,out]    iterm               Position in vector at which to insert the rate
 *
 *\return       Total rate of process
 */
 Sweep::real Sweep::Processes::DiffusionProcess::RateTerms(real t, const Cell &sys,
                                                           const Geometry::LocalGeometry1d& local_geom,
                                                           fvector::iterator &iterm) const {
     return *iterm++ = Rate(t, sys, local_geom);
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
 */
int  Sweep::Processes::DiffusionProcess::Perform(Sweep::real t, 
                              Sweep::Cell &sys, 
                              const Geometry::LocalGeometry1d& local_geom,
                              unsigned int iterm,
                              int (*rand_int)(int, int), 
                              Sweep::real(*rand_u01)(), 
                              Sweep::Transport::TransportOutflow *out) const
{
    // Rate of leftwards jumps (factors shared with rightwards jumps are ignored)
    real leftRate = 0;
    if(!local_geom.zeroGradient(Geometry::left)) {
        leftRate = 1.0 / local_geom.calcSpacing(Geometry::left);
    }

    // Rate of rightwards jumps (factors shared with leftwards jumps are ignored)
    real rightRate = 0;
    if(!local_geom.zeroGradient(Geometry::right)) {
        rightRate = 1.0 / local_geom.calcSpacing(Geometry::right);
    }

    // Choose direction weighted by rate
    Geometry::Direction direction;
    if(rand_u01() * (leftRate + rightRate) < leftRate) {
        direction = Geometry::left;
    }
    else {
        direction = Geometry::right;
    }

    // Get the particle that will be transported
    const int particleIndex = sys.Particles().Select(static_cast<TreeCache::PropID>(m_pid), rand_int, rand_u01);

    // LPDA updates and fictitious events are dealt with in the function below
    return Outflow(t, sys, local_geom, particleIndex, direction, out);
}
