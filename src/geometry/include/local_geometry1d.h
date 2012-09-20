/*!
 * \file   local_geometry1d.h
 * \author Robert I A Patterson
 *
 * \brief  Local spatial grid layout for 1d system
 *

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
#ifndef GEOMETRY_LOCAL_GEOMETRY1D_H
#define GEOMETRY_LOCAL_GEOMETRY1D_H

#include "types.h"

namespace Geometry {
    //forward declaration
    class Geometry1d;

/*!
 * The idea of this class is to hold the information about the geometry around
 * a particular cell so that it can be passed into the routines that handle
 * processes for a single cell, without those routines having to have any
 * awareness of the grid.  It is itended that instances of this class should have
 * short lifetimes, being created when needed and not kept in existence after use.
 *
 * The current implementation via a pointer back to the main geometry object
 * is easy to code, but may need reviewing.  In particular instances of this class
 * do not own the object pointed to by mGeom and will not delete it.
 *
 * No methods should be called on default constructed class instances, because
 * all operations require a pointer to a Geometry1d instance.
 *
 *@brief    Localised geometry information to avoid global data sharing
 */
class LocalGeometry1d {
public:
    //! Empty instance to act as placeholder in 0d problems
    LocalGeometry1d() : mGeom(NULL), mIndex(0) {}

    //! Build an object with geometry information around the specified cell
    LocalGeometry1d(const Geometry1d &geom, const size_t cell_index);

    //! True is solution is assumed to be spatially homogeneous in direction.
    //bool zeroGradient(const Direction direction) const;

    //! Centre of current cell
    real cellCentre() const;

    //! Find index of destination cell
    int calcDestination(const Direction direction) const;

    //! Distance to next cell centre in specified direction
    real calcSpacing(const Direction direction) const;

    //! Unordered list of cell vertices
    fvector cellVertices() const;

    //! Check if a position is within the local cell
    bool isInCell(const real x) const;

    //! Volume of cell
    real cellVolume() const;

    //! Volume of cell
    real cellVolume(const Direction direction) const;

private:
    //! Original geometry object which holds the layout information
    const Geometry1d *mGeom;

    //! Index of cell to which the object is localised
    size_t mIndex;
};

} //namespace Geometry

#endif // GEOMETRY_LOCAL_GEOMETRY1D_H

