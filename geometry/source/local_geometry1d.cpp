/*!
 * \file   local_geometry1d.cpp
 * \author Robert I A Patterson
 *
 * \brief Local spatial grid layout for 1d system

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

#include "local_geometry1d.h"

#include "geometry1d.h"
#include "types.h"


/*! 
 * Build an object with geometry information around the specified cell
 *
 *@param[in]	geom	     Global geometry information, of which this object will contain a local subset
 *@param[in]	cell_index   Index of cell in geom to which this object is to be localised
 */
Geometry::LocalGeometry1d::LocalGeometry1d(const Geometry1d &geom, const size_t cell_index)
    : mGeom(&geom)
    , mIndex(cell_index)
{}

/*!
 * Work out the index of destination cell for a particle.  All cells will have
 * index >= 0, but this method may return < 0 if a particle was transported
 * out of the domain of this geometry.
 *
 *\param[in]    direction     Direction in which particle is to move
 *
 *\return   Index of destination cell or a negative number of no destination
 */
int Geometry::LocalGeometry1d::calcDestination(const Direction direction) const {
    return mGeom->calcDestination(mIndex, direction);
}

/*!
 * Is transport in the indicated direction permitted by the geometry.  The
 * most likely (currently only) reason for a negative answer is a Neumann
 * boundary condition.
 *
 *\param[in]	direction	Proposed transport direction
 *
 *\return       True if transport is permitted in direction indicated from the local cell
 */
bool Geometry::LocalGeometry1d::permittedDirection(const Direction direction) const {
    // If the direction is included in the blocked directions mask will return non-zero (=true) when
    // the bitwise and is carried out.  This is negated to get the correct return value.
    return !(mGeom->blockedDirections(mIndex) & direction);
}

/*!
 * Calculate the local grid spacing in one direction
 *
 *@param[in]    direction    Direction for which spacing is required
 *
 *@return       Size of grid in specified direction
 */
Geometry::real Geometry::LocalGeometry1d::calcSpacing(const Direction direction) const {
    return mGeom->calcSpacing(mIndex, direction);
}

