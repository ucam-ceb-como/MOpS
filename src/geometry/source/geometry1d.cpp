/*!
 * \file   geometry1d.cpp
 * \author Robert I A Patterson
 *
 * \brief Spatial grid layout for 1d system
 *
 * \mainpage Library to represent spatial structure of a one dimensional system.

 Copyright (C) 2009 Robert I A Patterson.

 \section Licence
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

  \section Contact
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

#include "geometry1d.h"

#include "camxml.h"

#include <string>
#include <cstdlib>
#include <stdexcept>
#include <algorithm>
#include <sstream>
#include <limits>

const Geometry::real Geometry::Geometry1d::sCrossSectionalArea = 1.0;

/*!
 * Read the geometry from an xml tree which should have a root node of
 * type "geometry" which contains a list of at least two elements of 
 * type "vertex", which will be interpreted as the cell boundaries.  
 * 
 * Boundary conditions are also read in since they have to be specified
 * on the same grid.
 *
 *\param[in]   xml        XML node of type "geometry" containing grid data
 *
 *\exception   std::runtime_error    Less than two vertices specified
 *\exception   std::runtime_error    Provided XML is not a geometry element
 */
Geometry::Geometry1d::Geometry1d(const CamXML::Element &xml) {
    if(xml.Tag() == "geometry") {
        // Now get the positions of the cell boundaries
        std::vector<CamXML::Element*> verticesXML;
        xml.GetChildren("vertex", verticesXML);

        // Must have at least two cell boundaries
        if(verticesXML.size() < 2) {
            throw std::runtime_error("Geometry XML must contain at least two <vertex> elements\n");
        }

        // Get enough memory to store all the cell boundaries
        mCellEnds.reserve(verticesXML.size());

        std::vector<CamXML::Element*>::const_iterator it = verticesXML.begin();
        const std::vector<CamXML::Element*>::const_iterator itEnd = verticesXML.end();
        while(it != itEnd) {
            // Read the string representing the boundary position and increment the iterator
            std::string textData = (*it++)->Data();

            // Convert the string into a number and store it
            real boundary = atof(textData.c_str());
            mCellEnds.push_back(boundary);
        }

        // Sort the vertices into ascending order
        std::sort(mCellEnds.begin(), mCellEnds.end());

        // Finally read in the boundary conditions
        std::vector<CamXML::Element*> boundariesXML;
        xml.GetChildren("boundary", boundariesXML);

        // Must have a boundary condition for each end
        if(boundariesXML.size() != 2) {
            throw std::runtime_error("Geometry XML must contain exactly two <boundary> elements\n");
        }
        else {
            // reuse the iterator 'it' from above
            for(it = boundariesXML.begin(); it != boundariesXML.end(); ++it) {
                // Only two types of boundary are supported
                BoundaryConditionType boundaryType;
                CamXML::Element *boundaryTypeNode = (*it)->GetFirstChild("dirichlet");
                if(boundaryTypeNode != NULL) {
                    boundaryType = dirichlet;
                }
                // Note deliberate assignment before comparison
                else if ((boundaryTypeNode = (*it)->GetFirstChild("neumann")) != NULL) {
                    boundaryType = neumann;
                }
                else {
                    throw std::runtime_error("Failed to recognise boundary condition type\n");
                }

                // Work out to which end of the domain this boundary condition applies
                CamXML::Element *positionNode = (*it)->GetFirstChild("position");
                if(positionNode == NULL) {
                    throw std::runtime_error("Each <boundary> element must contain exactly one <position> element\n");
                }
                else {
                    // Position should be described by a string saying "left" or "right"
                    std::string positionText = positionNode->Data();
  
                    // Initialise the appropriate data member of this class
                    if(positionText == "left") {
                        mLeftBoundary = boundaryType;
                    }
                    else if (positionText == "right") {
                        mRightBoundary = boundaryType;
                    }
                    else {
                        throw std::runtime_error("The <position> of a <boundary> element must be left or right\n");
                    }
                }
            }
        }
    }
    else {
        throw std::runtime_error("XML passed to Geometry1d::Geometry1d() was not a <geometry> element");
    }
}

/*!
 * Construct geometry according to a detailed specification.
 * 
 *@param[in]    vertices        List of cell boundaries (not necessarily sorted)
 *@param[in]    left_boundary   Boundary condition at cell 0
 *@param[in]    right_boundary  Boundary condition at final cell
 * 
 *@exception   std::invalid_argument    Less than two vertices specified
 */
Geometry::Geometry1d::Geometry1d(const fvector& vertices,
                       const BoundaryConditionType left_boundary,
                       const BoundaryConditionType right_boundary)
    : mCellEnds(vertices)
    , mLeftBoundary(left_boundary)
    , mRightBoundary(right_boundary)
{
    // Sort the vertices into ascending order
    std::sort(mCellEnds.begin(), mCellEnds.end());

    if(vertices.size() < 2u)
        throw std::invalid_argument("Geometry must have at least two vertices (Geometry1d::Geometry1d)");
}

/*!
 * Calculate the number of cells specified by the geometry.
 * This will generally be one less than the number of vertices,
 * except that 0 vertices means 0 cells.
 *
 *\return       Number of cells in the geometry
 */
size_t Geometry::Geometry1d::numCells() const {
    if(mCellEnds.size() > 1)
        return mCellEnds.size() - 1;
    else
        return 0;
}

/*!
 * Calculate centre of cell.  No index checking.
 *
 *\param[in]    cell_index      Index of cell for which position is requested
 *
 *\return       Position of cell centre
 */
Geometry::real Geometry::Geometry1d::cellCentre(const size_t cell_index) const {
    return (mCellEnds[cell_index] + mCellEnds[cell_index + 1]) / 2;
}

/*!
 * Get all the cell vertices, the order in which the vertices is not guaranteed,
 * while it is easy to see the behaviour at the moment, future implementations 
 * my alter this.
 *
 * No index checking.
 *
 *\param[in]    cell_index      Index of cell for which vertices are requested
 *
 *\return       Vector of vertex positions
 */
Geometry::fvector Geometry::Geometry1d::cellVertices(const size_t cell_index) const {

    // In 1d case there will be exactly two vertices, the left and right hand ends
    fvector vertices(2);

    // left end
    vertices[0] = mCellEnds[cell_index];

    // right end
    vertices[1] =  mCellEnds[cell_index + 1];

    return vertices;
}

/*!
 * In the 1d case the cross section area is assumed to be 1, so the volume is
 * proportional to the length of the cell.
 *
 *@param[in]    cell_index      Index of cell for which vertices are requested
 *
 *@return       Physical volume covered by the cell
 */
Geometry::real Geometry::Geometry1d::cellVolume(const size_t cell_index) const {
    return (mCellEnds[cell_index + 1] - mCellEnds[cell_index]) * sCrossSectionalArea;
}

/*!
 * Find the index of the cell containing the specified position or return
 * a negative number to show that the position is outside the area covered
 * by the geometry.  Cells are regarded as closed at the starting (left) end and
 * open at the other (right) end.  The reverse effect would be achieved by using
 * lower_bound in place of upper_bound.  The current convention is chosen on the
 * so that particles can enter the positive half line at 0.
 *
 *@param[in]    x   Position for which cell index requested
 *
 *@return       Cell index or -1 if x is outside the range of the geometry
 */
int Geometry::Geometry1d::containingCell(const real x) const {
    const fvector::const_iterator it = std::upper_bound(mCellEnds.begin(), mCellEnds.end(), x);

   if(it == mCellEnds.begin()) {
       // Position is strictly before start of first cell
       return -1;
   }
   else if(it == mCellEnds.end()) {
       // Position is at or after end of last cell
       return -1;
   }

   return std::distance(mCellEnds.begin(), it) - 1;
}

/*!
 * Work out the index of destination cell for a particle.  All cells will have
 * index >= 0, but this method may return < 0 if a particle was transported
 * out of the domain of this geometry.
 *
 *\param[in]    origin_index    Index of cell from which the transport is starting
 *\param[in]    direction       Direction in which particle is to move
 *
 *\return   Index of destination cell or a negative number if no destination
 *
 *@exception    std::logic_error    Unhandled enum value
 */
int Geometry::Geometry1d::calcDestination(const size_t origin_index, const Direction direction) const {
    int newIndex;
    switch(direction) {
        case left:
            // Moving one cell left will reduce the index by 1.  If origin_index
            // is already 0, the leftmost cell, then the particle will leave the
            // domain of the geometry as shown by the destination of -1.
            newIndex = origin_index - 1;
            break;
        case right:
            if((origin_index + 1) == numCells()) {
                // Moving right from the rightmost cell involves leaving the
                // domain so a negative value must be returned.  Here the most
                // negative possible number is chosen, to give the effect of
                // wrapping round from the maximum value of int.
                newIndex = std::numeric_limits<int>::min();
            }
            else {
                // Moving right is simply a case of going to the next cell
                newIndex = origin_index + 1;
            }
            
            break;
        case none:
            newIndex = -1;
            break;
        default:
            throw std::logic_error("Unhandled enum value");
    }
    return newIndex;
}

/*!
 * Work out the distance to the next cell centre in the specified direction.  If
 * the grid contains no more cells in the specified direction or only one cell
 * return the length of the cell
 *
 *\param[in]    cell_index      Index of cell from which the transport is starting
 *\param[in]    direction       Direction in which particle is to move
 *
 *\return   Distance to next cell centre
 */
Geometry::real Geometry::Geometry1d::calcSpacing(const size_t cell_index, const Direction direction) const {
    real space;

    // Only one cell so use its length
    if(numCells() == 0) {
        space = mCellEnds[1] - mCellEnds[0];
    }
    else if((cell_index == 0               && direction == left) ||
           ((cell_index + 1) == numCells() && direction == right)) {
            // No cells in specified direction
         space = mCellEnds[cell_index + 1] - mCellEnds[cell_index];
    }
    else {
        // Finally get to the simple case
        switch(direction) {
            case left:
                space = cellCentre(cell_index) - cellCentre(cell_index - 1);
                break;
            case right:
                space = cellCentre(cell_index + 1) - cellCentre(cell_index);
                break;
            case none:
                space = mCellEnds[cell_index + 1] - mCellEnds[cell_index];
                break;
            default:
                throw std::logic_error("Unhandled enum value");
        }
    }
    return space;
}

/*!
 * Note that cells are taken to be closed at their left or lower end and open
 * at their right or upper end.  This method need to be consistent with the
 * \ref containingCell method also in this class.
 *
 *@param[in]    cell_index      The cell which may contain the position
 *@param[in]    x               Position to check for inclusion in cell
 *
 *@return       True iff x in contained in cell with index cell_index
 */
bool Geometry::Geometry1d::isInCell(const size_t cell_index, const real x) const {
    return (mCellEnds[cell_index] <= x && x < mCellEnds[cell_index + 1]);
}

/*!
 *@return The cell boundaries in a human readable form useful for logging.
 */
std::string Geometry::Geometry1d::printMesh() const {
    std::ostringstream log;
    log << "Geometry1d: ";

    fvector::const_iterator it = mCellEnds.begin();
    const fvector::const_iterator itEnd = mCellEnds.end();
    while(it != itEnd) {
        log << *it++ << ' ';
    }

    return log.str();
}

