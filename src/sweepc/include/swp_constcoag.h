/** 
 * File:   swp_constcoag.h
 * Author: riap2
 * Project:        sweep (population balance solver)
 * Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2011 Robert I A Patterson.

  File purpose:
    Declaration of the coagulation process for the constant kernel.

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

#ifndef SWP_CONSTCOAG_H
#define	SWP_CONSTCOAG_H

#include "swp_coagulation.h"

namespace Sweep
{
// Forward declare Mechanism class.
class Mechanism;

namespace Transport
{
struct TransportOutflow;
} // namespace Transport

namespace Processes
{
class ConstantCoagulation : public Coagulation
{
public:
    // Constructors.
    ConstantCoagulation(const Sweep::Mechanism &mech); // Default constructor.

    ConstantCoagulation(                     // Stream-reading constructor.
        std::istream &in,            //  - Input stream.
        const Sweep::Mechanism &mech //  - Parent mechanism.
        );
    
    //* Virtual destructor
    virtual ~ConstantCoagulation() {};
    
    virtual ConstantCoagulation* const Clone() const {return new ConstantCoagulation(*this);};

    //* Returns the process type for identification during serialisation
    virtual ProcessType ID(void) const {return Constant_Coagulation_ID;};


    // TOTAL RATE CALCULATION.

    // Returns the rate of the process for the given system.
    virtual real Rate(real t,         // Time.
                      const Cell &sys // System for which to calculate rate.
                      ) const;


    // RATE TERM CALCULATION.

    // Returns the number of rate terms for this process.
    virtual unsigned int TermCount(void) const;

    // Calculates the rate terms given an iterator to a real vector. The 
    // iterator is advanced to the position after the last term for this
    // process.  Returns the sum of all rate terms.
    virtual real RateTerms(
        real t,                  // Time.
        const Cell &sys,       // Indicates true kernel (not majorant).
        fvector::iterator &iterm // Iterator to the first term.
        ) const;

    //! Perform a coagulation with particles chosen according to the additive kernel
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
    // Coagulation rate types.  These define how the rate is 
    // calculated and how the particles are chosen.
    static const unsigned int TYPE_COUNT = 1;
    enum TermType {
        ConstantTimesConstant,
    };

    //* Calculate kernel between two particles
    virtual real CoagKernel(const Particle &sp1, const Particle &sp2,
                            const Cell& sys) const;

    //* Calculate majorant kernel between two particles
    virtual real MajorantKernel(const Particle &sp1, const Particle &sp2,
                                const Cell& sys, const MajorantType maj) const;

    //* Arbitrary factor to give some headroom for LPDA
    static const real s_MajorantFactor;
};

} // namespace Processes

} // namespace Sweep

#endif	/* SWP_CONSTCOAG_H */

