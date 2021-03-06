/*!
 * \file   reactor1d.h
 * \author Robert I A Patterson
 *
 * \brief  1d system of touching mops reactors
 *
 *  Copyright (C) 2009 Robert I A Patterson.
 *

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
#ifndef BRUSH_REACTOR1D_H
#define BRUSH_REACTOR1D_H

#include "brush_params.h"
#include "geometry1d.h"

#include "swp_cell.h"
#include "swp_mechanism.h"
#include "swp_particle.h"

#include "gpc_mech.h"

#include "linear_interpolator.hpp"

#include <vector>
#include <set>
#include <sstream>

namespace Brush {

// Forward declaration
class ResetChemistry;

/*!
 *@brief Struct for passing around particle populations and their positions.
 *
 *The particle list will not delete any particles for which it contains pointers,
 *all responsibility for instance ownership remains with calling code.  This is
 *a dumb container only.
 */
struct ParticlePopulationPoint1d {
    //! Position to which this population applied
    double position;

    //! Number density of this population
    double m0;

    //! Collection of pointers to particles
    Sweep::PartPtrList particleList;
};

/*!
 *\brief Array of mops reactors
 *
 * This 1-dimensional reactor contains 'cells', each of which is bounded by
 * two grid points.  Each cell contains a (sub-)reactor which is responsible
 * for the chemical and particle population.
 *
 * Note that the order of member declaration is important - initialisation
 * takes place in declaration order and I rely on this (standard required)
 * behaviour in one of the constructors.
 *
 * Cells are indexed with a zero based scheme.
 *
 * @todo copy constructor and assignment operator need to be implemented
 * so that all reactors refer to the Mechanism instance mMech.
 */
class Reactor1d {
public:
    //! Type of reactor in each grid cell
    typedef Sweep::Cell ReactorType;

    //! Initialise all cells with the given mech and particle capacity
    Reactor1d(const Geometry::Geometry1d &geom, const Sprog::Mechanism &g_mech,
              const bool constant_gas, const Sweep::Mechanism &p_mech,
              const Utils::LinearInterpolator<double, double> &max_particle_counts,
              const Utils::LinearInterpolator<double, double> &max_m0s);

    Reactor1d(const Reactor1d &src);

    //! Position of centre of cell
    double getCellCentre(const size_t i) const {return mGeometry.cellCentre(i);}

    //! Number of cells
    size_t getNumCells() const {return mGeometry.numCells();}

    //! Read only access to the particle mechanism
    const Sweep::Mechanism& getParticleMechanism() const {return mParticleMech;}

    //! Read only access to the gas mechanism
    const Sprog::Mechanism& getGasMechanism() const {return mGasMech;}

    //! Read only access to the geometry
    const Geometry::Geometry1d& getGeometry() const {return mGeometry;}

    //! Look at the contents of cell i
    const ReactorType& getCell(const size_t i) const {return mReactors[i];}

    //! Write access to the contents of cell i
    ReactorType& getCell(const size_t i) {return mReactors[i];}

    //! Replace chemical contents using supplied object
    void ReplaceChemistry(const ResetChemistry& chem_source, const bool fixed_chem);

    //! Replace particle populations
    template<typename IteratorType> void ReplaceParticles(const IteratorType& itBegin, const IteratorType& itEnd);

    //! Set the time on the constituent reactors
    void setTime(const double t);

private:
    //! Type in which to hold the reactors for the cells
    typedef std::vector<ReactorType> ReactorVector;

    //! Mechanism for the gas phase shared between the reactors in all the cells
    Sprog::Mechanism mGasMech;

    //! Mechanism for particle reactions shared between the reactors in all the cells
    Sweep::Mechanism mParticleMech;

    //! Spatial positions and relationships of the cells
    Geometry::Geometry1d mGeometry;

    //! The reactors defining the contents of the cells
    ReactorVector mReactors;
};

}

/*!
 * Provide some basic support for setting the particles in a cell.  Ownership of the
 * particles and responsibility for deleting them is taken by the reactor.  The method
 * looks up the containing cell for the position of each supplied particle population,
 * it does not perform any kind of interpolation of spreading of particle populations.
 *
 *@tparam      IteratorType          Iterator with value type ParticlePopulationPoint1d
 *@param[in]   itBegin               Start of a sequence of iterators to pointers to particles
 *@param[in]   itEnd                 One past end of a sequence of iterators to pointers to particles
 *
 *
 *
 *@exception   std::runtime_error     Attempt to set two populations in one cell
 */
template<typename IteratorType> void Brush::Reactor1d::ReplaceParticles(const IteratorType& itBegin,
                                                                        const IteratorType& itEnd) {
    // Make a copy of the iterator so that it can be incremented
    IteratorType it(itBegin);

    // Keep track of which cells have had their population replaced so that no cell has its
    // population replaced twice (destroying the information from the first replace with
    // the second).
    std::set<int> visitedCells;

    while(it != itEnd) {
        // Find the cell into which the current load of particles should be put
        const int cellIndex = mGeometry.containingCell(it->position);

        if(visitedCells.count(cellIndex) > 0) {
            // Attempt to replace population twice in one call to this method
            std::ostringstream msg;
            msg << "Attempt to replace particle population in cell " << cellIndex;
            msg << " twice during call to Reactor1d::ReplaceParticles";
            throw std::runtime_error(msg.str());
        }
        else {
            // Note that this cell has now been visited for particle population replacement
            visitedCells.insert(cellIndex);

            // Put the new particles in the cell
            if(it->particleList.size() > 0)
                mReactors[cellIndex].SetParticles(it->particleList.begin(), it->particleList.end(),
                                                             it->m0 / it->particleList.size());
            else
                mReactors[cellIndex].Reset(mReactors[cellIndex].SampleVolume());
        }
        ++it;
    }
}


#endif // BRUSH_REACTOR1D_H

