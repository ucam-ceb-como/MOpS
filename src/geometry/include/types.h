/*!
 * \file   types.h
 * \author Robert I A Patterson
 *
 * \brief  Types useful for geometries and their manipulation
 
 Copyright (C) 2009 Robert I A Patterson.

 Licence:
    This file is part of "brush".

    brush is free software; you can redistribute it and/or
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

#ifndef GEOMETRY_TYPES_H
#define GEOMETRY_TYPES_H

#include <vector>
#include <cstddef>

namespace Geometry {

    //! Indicate direction of transport
    /*!
     * Enum is intended to be extended or even replaced with
     * some kind of vector when more complex geometries are
     * implemented.
     *
     * It is current implemented as a kind of bit field
     * so that one can do things like myDirection & left
     * to find out if the move is a leftwards move.  This is
     * useful for indicating transport directions that are prohibited,
     * for example, no possible transport would be indicated by
     * 0x03 == 0x01 & 0x02.  Note that such combinations are not
     * enum values, but integers for which a typedef is provided
     * below.
     */
    enum Direction {
        none  = 0x00,
        left  = 0x01,
        right = 0x02,
    };

    //! Floating point type to use throughout
    typedef double real;

    //! Vector of floating point numbers
    typedef std::vector<real> fvector;

    //! Specify how solution behaves at edges of domain
    enum BoundaryConditionType {
        //! Solution has no particles at boundary
        dirichlet,

        //! Particle distribution is constant along lines normal to boundary
        neumann,
    };


} //namespace Geometry

#endif // GEOMETRY_TYPES_H
