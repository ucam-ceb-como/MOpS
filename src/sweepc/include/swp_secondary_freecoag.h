/*!
 * \file   swp_secondary_freecoag.h
 * \author Robert I A Patterson
 *  Copyright (C) 2010 Robert I A Patterson.
 *
 *  Project:        sweepc (population balance solver)
 *  Sourceforge:    http://sourceforge.net/projects/mopssuite
 *
 * \brief  Process for free molecular coagulation of small (secondary) particles

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
    Prof Markus Kraft
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
#ifndef SWP_SECONDARY_FREECOAG_H
#define	SWP_SECONDARY_FREECOAG_H

#include "swp_coagulation.h"

namespace Sweep {

namespace Transport
{
    // Forward declare structure to hold details of particle for onward transport
    struct TransportOutflow;
}

namespace Processes {

//! Free molecular coagulation of small (secondary) particles
class SecondaryFreeCoag : public Coagulation {
public:
    //! Create an instance that is part of a mechanism
    SecondaryFreeCoag(const Sweep::Mechanism &mech);

    //! Construct from a binary stream
    SecondaryFreeCoag(std::istream &in, const Sweep::Mechanism &mech);

    //! Create a copy
    SecondaryFreeCoag* const Clone() const {return new SecondaryFreeCoag(*(this));}

    //! Number of rate terms for secondary free molecular coagulation
    virtual unsigned int TermCount() const;

    //! Enum identifying the process
    virtual ProcessType ID() const;

    //! Returns rate of the process for the given system.
    virtual real Rate(
        real t,
        const Cell &sys
        ) const;

    // Calculates the rate terms given an iterator to a real vector. The
    // iterator is advanced to the position after the last term for this
    // process.  Returns the sum of all rate terms.
    virtual real RateTerms(
        real t,                  // Time.
        const Cell &sys,         // System for which to calculate rate terms.
        fvector::iterator &iterm // Iterator to the first term.
        ) const;

    //! Perform a coagulation of two secondary particle chosen according to the free molecular regime kernel
    virtual int Perform(
        real t,
        Cell &sys,
        const Geometry::LocalGeometry1d& local_geom,
        unsigned int iterm,
        int (*rand_int)(int, int),
        real(*rand_u01)(),
        Transport::TransportOutflow *out = 0
        ) const;

protected:
    //! Free molecular kernel for two particles
    virtual real CoagKernel(const Particle &sp1, const Particle &sp2,
                       const Cell &sys) const;

    //! Majorant kernel for any two secondary particles
    virtual real MajorantKernel(const Particle &sp1, const Particle &sp2,
                                const Cell &sys, const MajorantType maj) const;

private:
    //! Majorant kernel for any two secondary particles
    real MajorantKernel(const Cell &sys) const;

    //! Move a particle into the main population
    void MoveToMainPopulation(Particle* sp, Cell &sys, int (*rand_int)(int, int), real(*rand_u01)()) const;

    //! Number of rate terms for secondary free molecular coagulation
    static const unsigned int mTermCount = 1;

    //! Free-molecular enhancement factor, specific to soot
    static const real m_efm;
};

} // namespace Processes
} // namespace Sweep

#endif	/* SWP_SECONDARY_FREECOAG_H */

