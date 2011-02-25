/*!
 * @file     swp_molecule_evolution.cpp
 * @author    Robert Patterson

  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Robert I A Patterson.

  @brief Store details of offline calculations of molecular growth

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

#include "swp_molecule_evolution.h"

#include <cstdlib>
#include <stdexcept>
#include <algorithm>
#include <cassert>
#include <sstream>

#include "csv_io.h"

#include "first_of_pair_comparator.hpp"

using namespace Sweep::MoleculeEvolution;

/*!
 * @param[in]   file_name           Path of file from which to read data
 * @param[in]   formation_time      Time at which all molecules in file assumed to have formed
 * @param[in]   formation_position  Position at which all molecules in file assumed to have formed
 *
 * File is assumed to consist of columns:
 * 1 Age of molecules as described on the relevant row
 * 2 Ignored
 * 3 Integer size of molecule
 * 4 Additional ignored data regarding molecule
 * Columns 3 and 4 may be repeated as required to record the development of multiple
 * molecules within one file, however at least one molecule must be specified.
 * All rows must be the same length.
 *
 * @exception   std::runtime_error  No molecule information found
 * @exception   std::runtime_error  Molecule stories do not start at age 0
 */
void Database::addStoriesFromFile(const std::string& file_name,
                                  const Sweep::real formation_time,
                                  const Sweep::real formation_position) {
    // Open the file for reading
    CSV_IO *csvinput = new CSV_IO(file_name,false);

    std::vector<std::string> dataLine;

    // Read first line to initialise loop and to find number of molecules in the file
    csvinput->Read(dataLine);

    // Each molecule has data spread over two columns and the first two
    // columns of the file do not contain molecule data (see main
    // comment).
    const molecule_id numStories = (dataLine.size() - 1) / 2;

    if(dataLine.size() < 4) {
        throw std::runtime_error("No molecule information found in " + file_name +
                                 " (Database::addStoriesFromFile)");
    }

    // Older database files need rebasing so that the molecule age is used and not
    // the flow time.
    const real startAge = std::atof(dataLine[0].c_str());
    if(startAge > 1.0e-9) {
        std::ostringstream msg;
        msg << "Molecule stories should start at 0.0 in " << file_name
            << " but actual start at " << startAge << " ("
            << dataLine[0] << ") ,"
            << " (Database::addStoriesFromFile)";
        throw std::runtime_error(msg.str());
    }

    // Collect molecule stories into this temporary object so that the database
    // is not changed until one can be sure that there will be no errors or
    // exceptions while reading the file.
    database_type newStories(numStories);

    // Loop over the lines in the file
    while(!dataLine.empty()) {

        const real age = std::atof(dataLine[0].c_str());

        for(molecule_id i = 0; i != numStories; ++i) {
            // Skip first two columns, which are not molecule state data, then take
            // the first of every two following columns
            const molecule_state data = std::atoi(dataLine[2 * i + 2].c_str());

            // Put the data into the story
            const std::pair<real, molecule_state> dataPoint(age, data);
            newStories[i].push_back(dataPoint);
        }

        // Move on to next line of the input file
        csvinput->Read(dataLine);
    }

    delete csvinput;

    // Set up the formation time and position caches, no exceptions are
    // expected after this point.
    // The initialisation of these caches depends on the way in which ids
    // are assigned to the molecules.  The current implementation assumes
    // the id is just the index into mMoleculeStories.

    const molecule_id nextId = mMoleculeStories.size();
    try {
        for(molecule_id i = 0; i != numStories; ++i) {
            mFormationTimes.push_back(std::make_pair(formation_time, nextId + i));
        }

        for(molecule_id i = 0; i != numStories; ++i) {
            mFormationPositions.push_back(std::make_pair(formation_position, nextId + i));
        }

        // Finally insert the new molecule stories into the main database
        mMoleculeStories.insert(mMoleculeStories.end(), newStories.begin(), newStories.end());
    }
    catch(...) {
        // Something went wrong updating the internal data so try to go back to the
        // initial state, by returning to the initial size and throwing away all data
        // added to the end of the arrays
        mFormationTimes.resize(nextId);
        mFormationPositions.resize(nextId);
        mMoleculeStories.resize(nextId);

        // pass the exception on to the caller
        throw;
    }


    // Check that the newly added molecules respect the time and position
    // ordering and resort if they do not.
    Utils::FirstOfPairComparator<real, molecule_id> comp;

    if(nextId > 0) {
        // Cannot break an ordering if there was nothing there to start with
        if (!comp(mFormationTimes[nextId - 1], formation_time)) {
            // Time ordering has been broken
            std::stable_sort(mFormationTimes.begin(), mFormationTimes.end(), comp);
        }

        if (!comp(mFormationPositions[nextId - 1], formation_time)) {
            // Position ordering has been broken
            std::stable_sort(mFormationPositions.begin(), mFormationPositions.end(), comp);
        }
    }


    // debug info
    /*{
        formation_lookup::const_iterator it = mFormationTimes.begin();
        std::cout << "just before end of addStoriesFromFile\nt,  id\n";
        while(it != mFormationTimes.end()) {
            std::cout << it->first << ',' << it->second <<'\n';
            ++it;
        }
        std::cout << std::endl;
    }*/
}

/*!
 * Empty all of the constituent data containers
 */
void Database::clear() {
    mMoleculeStories.clear();
    mFormationTimes.clear();
    mFormationPositions.clear();
}

/*!
 * @param[in]       molecule_id         Identifier of molecule for which state is required
 * @param[in]       age                 Age of molecule for which state is required
 * @param[in,out]   hint                Initialise to 0 and use value set from previous call for relevant molecule
 *
 * @return      state of molecule at specified age or last data point before
 *
 * @exception   std::runtime_error      Invalid molecule id
 */
Database::molecule_state Database::getMoleculeState(const molecule_id id, const Sweep::real age,
                                                    state_lookup_hint &hint) const {
    // Use a reference to make the following code shorter
    if(id >= mMoleculeStories.size()) {
        //std::cerr << "Bad molecule id: " << id << ", should be < " << mMoleculeStories.size() << std::endl;
        throw std::runtime_error("Invalid id passed to Database::getMoleculeState");
    }
    const molecule_story &story = mMoleculeStories[id];

    // This only works while molecule_state is a built in scalar type
    molecule_state result = std::numeric_limits<molecule_state>::max();

    if((hint + 1 < story.size()) && (age >= story[hint].first) && (age < story[hint + 1].first)) {
        // hint was correct
        //std::cout << "Found size of " << story[hint].second << " at first attempt\n";
        result = story[hint].second;
    }
    else if((hint + 2 < story.size()) && (age >= story[hint + 1].first) && (age < story[hint + 2].first)) {
        // Very likely alternative has occurred: ie age has moved one one step
        // Update the hint ready for next time and return the new state
        ++hint;
        result = story[hint].second;
    }
    else if(age >= 0.0 && age <= story.front().first) {
        // This is a sightly unusual case, but it is to avoid calling std::upper_bound below
        // with age at or before the first data point
        hint = 0;
        result = story.front().second;
    } else {

        // The obvious places to look for age did not work so do a full search

        // Upper bound returns an iterator to the first entry >= age, since age is a floating point
        // quantity we ignore the possibility that they are equal and assume the time point is > age.
        // (If there is some reason to expect exact equality - for example because molecule updates
        //  are based on the same time discretisation as the molecule story generation - and the
        //  resulting error is expected to be large then the case of time and age being exactly equal
        //  will require special handling.)
        // Based on the assumption, subtract 1 to get an interator to the last point in the story
        // before age is reached.
        Utils::FirstOfPairComparator<real, molecule_state> comp;
        const molecule_story::const_iterator it = std::upper_bound(story.begin(), story.end(), age, comp) - 1;

        // Try looking in the same place again next time
        hint = std::distance(story.begin(), it);
        result = it->second;
    }
    return result;
}

/*!
 * @param[in]   t       Time around which molecule formed
 * @param       rng     Function to generate uniform integers on [a, b] by the call rng(a, b)
 *
 * @return      Id for a molecule formed around the specified time
 *
 * @exception   std::logic_error    Empty database
 *
 * Changes to this method should generally be replicated in selectMoleculNearPosition
 */
Database::molecule_id Database::selectMoleculeNearTime(const Sweep::real t, int (*rng)(int, int)) const {
    // See if there are any molecules formed exactly at the specified time.
    // This is rather unlikely and so the range will generally contain 0 elements.
    Utils::FirstOfPairComparator<real, molecule_id> comp;
    std::pair<formation_lookup::const_iterator, formation_lookup::const_iterator> range
      = std::equal_range(mFormationTimes.begin(), mFormationTimes.end(), t, comp);

    // Store the range of molecules from which to select in this pair
    std::pair<formation_lookup::const_iterator, formation_lookup::const_iterator> rangeToUse = range;

    if(std::distance(range.first, range.second) == 0) {
    // Initial search did not find any molecules that formed exactly at t so
    // look at the closest molecules that formed at previous and later times.

        // Indicate if there were any molecules formed before t
        const bool havePrevious = (range.first != mFormationTimes.begin());

        // Indicate if there were any molecules formed after t
        const bool haveNext = (range.second != mFormationTimes.end());

        // Find a range of molecule ids from which to make a random selection
        if(havePrevious) {
            // Last time prior to t at which a molecule was formed
            const real tPrevious = (range.first - 1)->first;

            if(haveNext) {
                // First time after t at which a molecule was formed
                const real tNext = (range.second)->first;

                if((t - tPrevious) > (tNext - t)) {
                    // The molecules formed after t are closer than those formed before,
                    // so find all the molecules formed at tNext.
                    // Note that range.first already points to the first molecule
                    // formed at tNext.
                    range.second = std::upper_bound(range.first, mFormationTimes.end(), tNext, comp);
                }
                else {
                    // The molecules formed before t are closer than those formed after,
                    // so find all the molecules formed at tPrevious.
                    // Note that range.second already points to one past the last molecule
                    // formed at tPrevious.
                    range.first = std::lower_bound(mFormationTimes.begin(), range.second, tPrevious, comp);
                }
            }
            else {
                // Only have molecules formed prior to t
                std::cerr << "WARNING: " << t << " is after the last molecule formation time\n";

                // Have to look at the molecules formed prior to t
                // Note that range.second already points to one past the last molecule
                // formed at tPrevious.
                range.first = std::lower_bound(mFormationTimes.begin(), range.second, tPrevious, comp);
            }
        }
        else {
            // No particles formed before t
            if(haveNext) {
                std::cerr << "WARNING: " << t << " is before the first molecule formation time\n";

                // Look at the molecules formed after t since there are none formed any earlier
                range.second = std::upper_bound(range.first, mFormationTimes.end(), (range.first)->first, comp);

            }
            else {
                throw std::logic_error("Database::selectMoleculeNearTime called on empty database");
            }
        }
    }
    // Should now have a non-empty range of ids of molecules formed around t in [range.first, range.second)

    const formation_lookup::size_type numRelevantMolecules = std::distance(range.first, range.second);
    if(numRelevantMolecules < 10)
        std::cerr << "Only found " << numRelevantMolecules << " (Database::selectMoleculeNearTime)\n";

    //@todo The range should probably be extended to ensure minimum size

    // Get an iterator to one molecule chosen uniformly at random
    formation_lookup::const_iterator selectedMolecule(range.first);
    std::advance(selectedMolecule, rng(0, numRelevantMolecules - 1));
    return selectedMolecule->second;
}

/*!
 * @param[in]   x       Position around which molecule formed
 * @param       rng     Function to generate uniform integers on [a, b] by the call rng(a, b)
 *
 * @return      Id for a molecule formed around the specified position
 *
 * @exception   std::logic_error    Empty database
 *
 * Changes to this method should generally be replicated in selectMoleculNearTime
 */
Database::molecule_id Database::selectMoleculeNearPosition(const Sweep::real x, int (*rng)(int, int)) const {
    // See if there are any molecules formed exactly at the specified time.
    // This is rather unlikely and so the range will generally contain 0 elements.
    Utils::FirstOfPairComparator<real, molecule_id> comp;
    std::pair<formation_lookup::const_iterator, formation_lookup::const_iterator> range
      = std::equal_range(mFormationPositions.begin(), mFormationPositions.end(), x, comp);

    // Store the range of molecules from which to select in this pair
    std::pair<formation_lookup::const_iterator, formation_lookup::const_iterator> rangeToUse = range;

    if(std::distance(range.first, range.second) == 0) {
    // Initial search did not find any molecules that formed exactly at x so
    // look at the closest molecules that formed at previous and later positions.

        // Indicate if there were any molecules formed before x
        const bool havePrevious = (range.first != mFormationPositions.begin());

        // Indicate if there were any molecules formed after x
        const bool haveNext = (range.second != mFormationPositions.end());

        // Find a range of molecule ids from which to make a random selection
        if(havePrevious) {
            // Last time prior to x at which a molecule was formed
            const real xPrevious = (range.first - 1)->first;

            if(haveNext) {
                // First time after t at which a molecule was formed
                const real xNext = (range.second)->first;

                if((x - xPrevious) > (xNext - x)) {
                    // The molecules formed after x are closer than those formed before,
                    // so find all the molecules formed at xNext.
                    // Note that range.first already points to the first molecule
                    // formed at xNext.
                    range.second = std::upper_bound(range.first, mFormationPositions.end(), xNext, comp);
                }
                else {
                    // The molecules formed before x are closer than those formed after,
                    // so find all the molecules formed at xPrevious.
                    // Note that range.second already points to one past the last molecule
                    // formed at xPrevious.
                    range.first = std::lower_bound(mFormationPositions.begin(), range.second, xPrevious, comp);
                }
            }
            else {
                // Only have molecules formed prior to x
                //std::cerr << "WARNING: " << x << " is after the last molecule formation position\n";

                // Have to look at the molecules formed prior to x
                // Note that range.second already points to one past the last molecule
                // formed at xPrevious.
                range.first = std::lower_bound(mFormationPositions.begin(), range.second, xPrevious, comp);
            }
        }
        else {
            // No particles formed before x
            if(haveNext) {
                std::cerr << "WARNING: " << x << " is before the first molecule formation position\n";

                // Look at the molecules formed after t since there are none formed any earlier
                range.second = std::upper_bound(range.first, mFormationPositions.end(), (range.first)->first, comp);

            }
            else {
                throw std::logic_error("Database::selectMoleculeNearPosition called on empty database");
            }
        }
    }
    // Should now have a non-empty range of ids of molecules formed around t in [range.first, range.second)

    const formation_lookup::size_type numRelevantMolecules = std::distance(range.first, range.second);
    if(numRelevantMolecules < 10)
        std::cerr << "Only found " << numRelevantMolecules << " (Database::selectMoleculeNearPosition)\n";

    //@todo The range should probably be extended to ensure minimum size

    // Get an iterator to one molecule chosen uniformly at random
    formation_lookup::const_iterator selectedMolecule(range.first);
    std::advance(selectedMolecule, rng(0, numRelevantMolecules - 1));
    return selectedMolecule->second;
}
