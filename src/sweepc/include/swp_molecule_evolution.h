/*!
 * \file   swp_molecule_evolution.h
 * \author Robert I A Patterson
 *  Copyright (C) 2010 Robert I A Patterson.
 *
 *  Project:        sweepc (population balance solver)
 *  Sourceforge:    http://sourceforge.net/projects/mopssuite
 *
 * \brief  Store details of offline calculations of molecular growth

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

#ifndef SWP_MOLECULE_EVOLUTION_H
#define SWP_MOLECULE_EVOLUTION_H

#include <vector>
#include <string>

#include "swp_params.h"


namespace Sweep {
    class ParticleModel;

//! Details of offline calculations of molecular growth, shrinkage and other development
namespace MoleculeEvolution {

/*!
 * \brief  Collection of life stories for molecules that begin at different points in time and space.
 *
 * \invariant mMoleculeStories.size() == mFormationTimes.size()
 * \invariant mMoleculeStories.size() == mFormationPositions.size()
 *
 * Because this class contains only some stl containers and no pointers it
 * does not require specially implemented constructors, nor is a handwritten
 * destructor or assignment operator required.  However, this means that
 * making multiples copies of an instance will start to take up large amounts
 * of memory.  Currently this is not a problem because this class is used as
 * a member of ParticleModel, which (in the form of a Mechanism) is largely
 * shared via pointers, rather than copied inside Mops and Brush.
 */
class Database {
public:
    //! Load molecule stories from the specified file
    void addStoriesFromFile(const std::string& file_name, const real formation_time,
                                                          const real formation_position);

    //! Remove all information
    void clear();

    //! State of molecule at a particular moment in time @todo make a template parameter
    typedef int molecule_state;

private:
    //! To hold the life story of a molecule as a sequence of age, state pairs
    typedef std::vector<std::pair<real, molecule_state> > molecule_story;

    //! Type used to hold the main data in this particular implementation
    typedef std::vector<molecule_story> database_type;

public:
    //! Arbitrary unique identifier for molecule in the database
    typedef database_type::size_type molecule_id;

    //! Number of molecules stored in the database
    database_type::size_type size() const {return mMoleculeStories.size();}

    //! A hint to help find the point in the particle story appropriate for a supplied age
    typedef molecule_story::size_type state_lookup_hint;

    //! Select a molecule formed near a specified time
    molecule_id selectMoleculeNearTime(const real t,
                                       int (*rng)(int, int)) const;

    //! Select a molecule formed near a particular position
    molecule_id selectMoleculeNearPosition(const real x,
                                           int (*rng)(int, int)) const;

    //! Find out the state of a molecule given its age
    molecule_state getMoleculeState(const molecule_id id,
                                    const real age,
                                    state_lookup_hint &hint) const;

private:

    //! The collection of molecular development stories
    database_type mMoleculeStories;

    //! Use for caching molecule formation time/position against id
    typedef std::vector<std::pair<real, molecule_id> > formation_lookup;

    //! Formation time of each molecule
    formation_lookup mFormationTimes;

    //! Formation position of each molecule
    formation_lookup mFormationPositions;
};

} //namespace MoleculeEvolution
} //namespace Sweep

#endif /* SWP_MOLECULE_EVOLUTION_H */
