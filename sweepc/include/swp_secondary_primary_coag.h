/*!
 * \file   swp_secondary_primary_coag.h
 * \author Robert I A Patterson
 *  Copyright (C) 2010 Robert I A Patterson.
 *
 *  Project:        sweepc (population balance solver)
 *  Sourceforge:    http://sourceforge.net/projects/mopssuite
 *
 * \brief  Process for free molecular coagulation of small (secondary) particles with particles from the main population

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
#ifndef SWP_SECONDARY_PRIMARY_COAG_H
#define	SWP_SECONDARY_PRIMARY_COAG_H

#include "swp_coagulation.h"

namespace Sweep {

namespace Transport
{
    // Forward declare structure to hold details of particle for onward transport
    struct TransportOutflow;
}

namespace Processes {

//! Free molecular coagulation of small (secondary) particles with particles from the main population
class SecondaryPrimaryCoag : public Coagulation {
public:
    //! Create an instance that is part of a mechanism
    SecondaryPrimaryCoag(const Sweep::Mechanism &mech);

    //! Construct from a binary stream
    SecondaryPrimaryCoag(std::istream &in, const Sweep::Mechanism &mech);

    //! Create a copy
    SecondaryPrimaryCoag* const Clone() const {return new SecondaryPrimaryCoag(*(this));}

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

private:
    //! Free molecular kernel for two particles
    real FreeMolKernel(const Particle &sp1, const Particle &sp2,
                       real temperature) const;

    //! Majorant kernel for one secondary particle with a specified particle from the main population
    real MajorantKernel(const Cell &sys, const Particle &sp1) const;

    //! Number of rate terms for secondary free molecular coagulation
    static const unsigned int mTermCount = 4;

    //! Free-molecular enhancement factor, specific to soot
    static const real m_efm;
};

} // namespace Processes
} // namespace Sweep

#endif	/* SWP_SECONDARY_PRIMARY_COAG_H */

