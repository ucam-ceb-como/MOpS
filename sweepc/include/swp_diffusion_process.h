/*!
 * \file   swp_diffusion_process.h
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
#ifndef SWP_DIFFUSION_PROCESS_H
#define	SWP_DIFFUSION_PROCESS_H

#include "swp_transport_process.h"

#include "geometry1d.h"

#include <vector>

namespace Sweep {

namespace Transport
{
    // Forward declare structure to hold details of particle for onward transport
    struct TransportOutflow;
}

namespace Processes {

//! Transport of particles by jumps
class DiffusionProcess : public TransportProcess {
public:
    //! Create a copy of the transport process.
    virtual DiffusionProcess  *const Clone(void) const {return new DiffusionProcess (*this);}

    //! Return the process type for identification during serialisation
    virtual ProcessType ID() const;

    //! Set the ID number of the particle property used for rate calculation
    void SetPropertyID(unsigned int id);

    //! Set the power of temperature to use in rate calculations
    void SetTemperatureExponent(const real e) {m_TemperatureExponent = e;}

    //! Rate of the process for the given system.
    virtual real Rate(
        real t,
        const Cell &sys,
        const Geometry::LocalGeometry1d& local_geom
        ) const;

        // SINGLE PARTICLE RATE CALCULATIONS.

    //* Particle dependent part of process rate for use in ficitious event tests
    virtual real Rate(
        real t,
        const Cell &sys,
        const Particle &sp
        ) const;


    /*!
     * The majorant rate is simply a multiple of the true rate,
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
    virtual unsigned int TermCount() const {return 1;}


    //! Puts rate terms into a vector for transport in permitted directions
    real RateTerms(
        real t,
        const Cell &sys,
        const Geometry::LocalGeometry1d& local_geom,
        fvector::iterator &iterm
        ) const;


    // PERFORMING THE PROCESS.

    //! Performs the process on the given system.
    virtual int Perform(
        real t,
        Cell &sys,
        const Geometry::LocalGeometry1d& local_geom,
        unsigned int iterm = 0,
        Transport::TransportOutflow *out = 0
        ) const;

private:
    //! Particle property to which the rate of the process is proportional.
    unsigned int m_pid;

    //! Exponent of temperature in rate expression
    real m_TemperatureExponent;

    //! Ratio of majorant rate to true rate
    static const real s_MajorantFactor;
};

} // namespace Processes
} // namespace Sweep

#endif	/* SWP_DIFFUSION_PROCESS_H */

