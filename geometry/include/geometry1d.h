/*!
 * \file   geometry1d.h
 * \author Robert I A Patterson
 *
 * \brief Spatial grid layout for 1d system
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
#ifndef GEOMETRY_GEOMETRY1D_H
#define GEOMETRY_GEOMETRY1D_H

#include "types.h"

#include <string>
#include <vector>

// Forward declarations

namespace CamXML {
    class Element;
}

//! Layout and boundaries of grids/meshes
namespace Geometry {

/*!
 *\brief Position information for a linear array of touching reactors
 *
 * The geometry is defined by a list of vertices (cell ends), which
 * define the cells of the spatial discretisation.  The number of cells
 * will be one less than the number of vertices and any indices will
 * be zero based.
 */
class Geometry1d {
public:
    //! Read geometry from xml
    Geometry1d(const CamXML::Element &xml);

    //! Construct a space separated string listing the cell boundaries
    std::string printMesh() const;

    //! Number of cells
    size_t numCells() const;

    //! Position of cell centre
    real cellCentre(const size_t cell_index) const;

    //! Index of cell containing specified position
    int containingCell(const real x) const;

    //! Directions in which transport will not occur from a cell
    DirectionMask blockedDirections(const size_t cell_index) const;

    //! Find destination of a particle being transported
    int calcDestination(const size_t origin_index, const Direction direction) const;

    //! Calculate grid spacing in specified direction
    real calcSpacing(const size_t cell_index, const Direction direction) const;

private:
    //! Cells ends (should be 1 longer than array of reactors)
    fvector mCellEnds;

    //! Boundary condition at cell 0
    BoundaryConditionType mLeftBoundary;

    //! Boundary condition at final cell
    BoundaryConditionType mRightBoundary;
};

} //namespace Geometry

#endif // GEOMETRY_GEOMETRY1D_H

