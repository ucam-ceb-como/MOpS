/*!
 * \file   swp_advection_process.cpp
 * \author Robert I A Patterson
 *  Copyright (C) 2009 Robert I A Patterson.
 *
 *  Project:        sweepc (population balance solver)
 *  Sourceforge:    http://sourceforge.net/projects/mopssuite
 *
 * \brief  Process for stochastic realisation of advective transport

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

#include "swp_advection_process.h"

#include "swp_cell.h"
#include "swp_transport_outflow.h"

#include <cmath>

using namespace Sweep::Processes;

const Sweep::real Sweep::Processes::AdvectionProcess::s_MajorantFactor = 1.5;

/*!
 *
 * \return      Enum identifying this process as particle advection
 */
Sweep::Processes::ProcessType Sweep::Processes::AdvectionProcess::ID() const {
    return Advection_ID;
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
Sweep::real Sweep::Processes::AdvectionProcess::Rate(real t, const Cell &sys,
                                                     const Geometry::LocalGeometry1d& local_geom) const {
    // Bulk gas velocity
    const real velocity = sys.Velocity();

    // Negative velocity means leftwards movement, positive means rightward
    Geometry::Direction direction = (velocity > 0 ? Geometry::right : Geometry::left);

    const real rate = std::abs(velocity) / local_geom.calcSpacing(direction);
    
    return rate * sys.Particles().Count() * A() * s_MajorantFactor;
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
 Sweep::real Sweep::Processes::AdvectionProcess::RateTerms(real t, const Cell &sys,
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
int AdvectionProcess::Perform(Sweep::real t, Sweep::Cell &sys, 
                             const Geometry::LocalGeometry1d& local_geom,
                             unsigned int iterm,
                             int (*rand_int)(int, int), 
                             Sweep::real(*rand_u01)(), 
                             Sweep::Transport::TransportOutflow *out) const
{
    // Negative velocity means leftwards movement, positive means rightward
    Geometry::Direction direction = (sys.Velocity() > 0 ? Geometry::right : Geometry::left);

    // Choose a particle uniformly
    const int particleIndex = sys.Particles().Select(rand_int);;

    return Outflow(t, sys, local_geom, particleIndex, direction, rand_u01, out);
 }
