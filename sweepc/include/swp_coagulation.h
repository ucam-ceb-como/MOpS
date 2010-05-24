/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Definition of the coagulation process.  This coagulation process uses a transition
    kernel (REF).  It calculates the rates for the free molecular and slip-flow regimes.

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

#ifndef SWEEP_COAGULATION_H
#define SWEEP_COAGULATION_H

#include "swp_process.h"
#include <iostream>
    
namespace Sweep
{
// Forward declare Mechanism class.
class Mechanism;

namespace Processes
{
    //Forward declaration to make typedef possible
    class Coagulation;

    //! Vector of polymorphic coagulation processes
    typedef std::vector<Coagulation*> CoagPtrVector;

//! Processes that stick two particles together to form one
class Coagulation : public Process
{
public:

    //! Ordinary method of construction just passes argument through.
    Coagulation(const Sweep::Mechanism &mech);

    //! Scaling factor for rate.
    real A(void) const {return m_a;}

    //! Sets the rate constant.
    void SetA(real a) {m_a = a;}

    //! Returns a copy of the coagulation process
    virtual Coagulation *const Clone(void) const = 0;

    //! Calculate rates for a sequence of coagulation processes
    static real CalcRates(real t, const Cell &sys, const CoagPtrVector &coags,
                          fvector &rates, unsigned int start = 0);

    //! Calculate rate terms for a sequence of coagulation processes
    static real CalcRateTerms(real t, const Cell &sys, const CoagPtrVector &coags,
                              fvector::iterator &iterm);

    //! Writes the object to a binary stream.
    virtual void Serialize(std::ostream &out) const;

    //! Reads the object from a binary stream.
    virtual void Deserialize(
        std::istream &in,            // Input stream.
        const Sweep::Mechanism &mech // Parent mechanism.
        );

protected:

    /*!
     * \@brief  Default constructor is protected to prevent coagulations being
     *          defined without knowledge of the parent mechanism.
     */
    Coagulation(void) {};
    
    //! Actually stick two particles together
    int JoinParticles(const real t, const int ip1, Particle *sp1,
                      const int ip2, Particle *sp2,
                      Cell &sys, real(*rand_u01)()) const;

private:
    //! Scaling factor for rate
    real m_a;
};

} //namespace Processes
} //namespace Sweep

#endif
