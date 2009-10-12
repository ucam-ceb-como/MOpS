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

/*
 * This method is mainly provided for use in diffusion process jump rate
 * calculations.  It checks whether a zero spatial gradient condition is imposed
 * on the solution at a particular point and direction, which would ensure no
 * diffusion.
 *
 *\param[in]    direction       Direction in which to look at gradient conditions
 *
 *\return       True iff there is a zero gradient condition
 */
bool Geometry::LocalGeometry1d::zeroGradient(const Direction direction) const {
    return mGeom->zeroGradient(mIndex, direction);
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

/*!
 * Get all the cell vertices, the order in which the vertices is not guaranteed.
 *
 * If no geometry data is present (mGeom == NULL) then the vertices will be
 * returned as 0.
 *
 *\return       Vector of vertex positions
 */
Geometry::fvector Geometry::LocalGeometry1d::cellVertices() const {
    if(mGeom != NULL) {
        return mGeom->cellVertices(mIndex);
    }

    // No information available so return both end points as 0.
    return fvector(2, 0.0);
}

/*!
 * Check if x is contained in the local cell.  If there is no local cell, then
 * x is not inside a local cell.
 *
 *@param[in]    x       Position to check for inclusion in local cell
 *
 *@return       True iff x is inside the local cell
 */
bool Geometry::LocalGeometry1d::isInCell(const real x) const {
    if(mGeom != NULL) {
        // Have some information to check against
        return mGeom->isInCell(mIndex, x);
    }

    // There is no cell that could possibly contain x
    return false;
}