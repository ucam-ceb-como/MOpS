/*!
 * \file   swp_coag_weight_rules.h
 * \author Robert I A Patterson
 *  Copyright (C) 2011 Robert I A Patterson.
 *
 *  Project:        sweepc (population balance solver)
 *  Sourceforge:    http://sourceforge.net/projects/mopssuite
 *
 * \brief  Symbolic indices for the properties of particles
 *
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

#ifndef SWEEP_COAG_WEIGHT_RULES_H
#define SWEEP_COAG_WEIGHT_RULES_H


namespace Sweep
{
namespace Processes
{
    //! Symbolic names for the different ways of choosing the weight of the new particle in weight coagulation events.
    enum CoagWeightRule {
        //! \f$ \left(\frac{1}{u} + \frac{1}{v} \right)^{-1}\f$
        CoagWeightHarmonic,

        //! Divide by 2
        CoagWeightHalf,

        //! Mass conserving
        CoagWeightMass,

        //! \f$ u \frac{f(x,u)}{f(x,u) + f(y,v)}\f$ where \f$ f(x,u) = m(x) / \sqrt{u} \f$
        CoagWeightRule4,
    };
}
}

#endif
