/*!
 * \file   swp_advection_process.h
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
#ifndef SWP_ADVECTION_PROCESS_H
#define	SWP_ADVECTION_PROCESS_H

#include "swp_transport_process.h"

#include "geometry1d.h"

#include <vector>

namespace Sweep {

namespace Transport {
    struct TransportOutflow;
}

namespace Processes {

//! Transport of particles by jumps
class AdvectionProcess : public TransportProcess {
public:
    //! Create a copy of the transport process.
    virtual AdvectionProcess *const Clone(void) const {return new AdvectionProcess(*this);}

    //! Return the process type for identification during serialisation
    virtual ProcessType ID() const;

    //! Rate of the process for the given system.
    virtual real Rate(
        real t,
        const Cell &sys,
        //const Sweep::Transport::DirectionMask &transport_mask
        const Geometry::LocalGeometry1d& local_geom
        ) const;

    // SINGLE PARTICLE RATE CALCULATIONS.

    /*!
     * For advection the rate is independent of the particles (it only
     * depends on the bulk fluid velocity, so always return the same
     * positive value.
     *
     *@param[in]        t       Time at which rate is to be calculated
     *@param[in]        sys     System containing particle
     *@param[in]        sp      Particle for which to calculate rate
     *
     *@return       Particle dependent part of process rate
     */
    virtual real Rate(
        real t,
        const Cell &sys,
        const Particle &sp
        ) const {return 1.0;}


    /*!
     * For advection the rate is independent of the particles (it only
     * depends on the bulk fluid velocity, so always return the same
     * positive value, which is the majorant factor multiplied by the
     * non-majorant rate.
     *
     *@param[in]        t       Time at which rate is to be calculated
     *@param[in]        sys     System containing particle
     *@param[in]        sp      Particle for which to calculate rate
     *
     *@return       Particle dependent part of process majorant rate
     */
    virtual real MajorantRate(
        real t,
        const Cell &sys,
        const Particle &sp
        ) const {return Rate(t, sys, sp) * s_MajorantFactor;}


    // RATE TERM CALCULATIONS.
    //   These routines return the individual rate terms for a
    //   process, which may have multiple terms (e.g. condensation).

    //! Number of rate terms for this process.
    virtual unsigned int TermCount() const {return 1u;}


    //! Puts rate terms into a vector for transport in permitted directions
    real RateTerms(
        real t,
        const Cell &sys,
        //const Sweep::Transport::DirectionMask &transport_mask,
        const Geometry::LocalGeometry1d& local_geom,
        fvector::iterator &iterm
        ) const;

    // PERFORMING THE PROCESS.

    //! Performs the process on the given system.
    int Perform(
        real t,
        Cell &sys,
        const Geometry::LocalGeometry1d& local_geom,
        unsigned int iterm = 0,
        Transport::TransportOutflow *out = 0
        ) const;

private:
    //! Ratio of majorant rate to true rate
    static const real s_MajorantFactor;

};

} // namespace Processes
} // namespace Sweep

#endif	/* SWP_ADVECTION_PROCESS_H */

