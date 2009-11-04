/*!
 * \file   reset_chemistry.cpp
 * \author Robert I A Patterson
 *
 * \brief  Routines for solving with predictor corrector coupling
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
#include "reset_chemistry.h"

#include "gpc_mech.h"
#include "comostrings.h"
#include "mops_reactor.h"

#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>
#include <cstdlib>
#include <cassert>

using namespace Brush;

/*!
 * For any element of mInputChemistryData the range
 * [0, sNumNonSpeciesData) will contain all of the data that
 * is not species mass or mole fractions.
 *
 * Index 0 is the spatial position to which the
 * data applied, the other indices in the range should
 * have their own static constants.
 */
const size_t ResetChemistry::sNumNonSpeciesData = 4;

/*!
 * Position of temperature data in the elements of mInputChemistryData.
 * The data will always have units K.
 */
const size_t ResetChemistry::sTemperatureIndex = 1;

/*!
 * Position of density data in the elements of mInputChemistryData.
 * The data will always be have units \f$\mathrm{kg\,m^{-3}
 */
const size_t ResetChemistry::sDensityIndex = 2;

/*!
 * Position of velocity data in the elements of mInputChemistryData.
 * The data will always have units \f$\mathrm{m\,s^{-1}}\f$
 */
const size_t ResetChemistry::sVelocityIndex = 3;

/**
 * Construct an object from a data file and associate the concentration data
 * with the indices defined by the mechanism.
 *
 * The data file should contain one header row, which is case sensitive followed
 * by complete rows of numerical data, completely blank lines are permitted and
 * ignored.  The columns may be separated by spaces or tabs.
 *
 * Mandatory columns in the file are (Camflow format \see InputFileType)
 * - x Spatial position to which the row of data applies (\f$\mathrm{m}\f$)
 * - T Temperature (K)
 * - rho Mass density of the gas mixture (\f$\mathrm{kg\,m^{-3}}\f$)
 * - u Velocity of gas mixture (\f$\mathrm{m\,s^{-1}}\f$)
 * - A column of mass fraction data must be provided for each species in the
 *   mechanism and the column headings for these species must be identical to the
 *   strings specified as species names when the mechanism was constructed.
 *
 * Mandatory columns in the file are (Premix format \see InputFileType)
 * - x Spatial position to which the row of data applies (\f$\mathrm{cm}\f$)
 * - T Temperature (K)
 * - rho Mass density of the gas mixture (\f$\mathrm{g\,cm^{-3}}\f$)
 * - u Velocity of gas mixture (\f$\mathrm{cm\,s^{-1}}\f$)
 * - A column of mole fraction data must be provided for each species in the
 *   mechanism and the column headings for these species must be identical to the
 *   strings specified as species names when the mechanism was constructed.
 *
 *\param[in]    fname       Path of file from which to read the data
 *\param[in]    file_type   Style of data in file
 *\param[in]    mech        Mechanism defining species indexes
 *\param[in]    verbosity   Level of debugging information sent to standard out (higher is more)
 *
 *\exception    std::runtime_error  No data for a species present in the mechanism
 */
Brush::ResetChemistry::ResetChemistry(const string &fname, const InputFileType file_type, 
                                      const Sprog::Mechanism& mech, const int verbosity) {
    // Delimeters to use when splitting lines of data from file into the individual column entries
    const std::string delims(" ,\t");

    // Iterators to the species in the mechanism
    Sprog::Mechanism::const_sp_iterator spIt = mech.SpBegin();
    const Sprog::Mechanism::const_sp_iterator spItEnd = mech.SpEnd();

    // Vector of names of the species in the mechanism
    std::vector<std::string> speciesNames;
    // Also a vector that will hold the column numbers for these species in the data file
    std::vector<int> speciesFileIndices;

    // Position, temperature, velocity and density will both be needed as well as species names.
    // The order of these insertions is important for the interpretation of the data
    // in the apply method
    switch(file_type) {
        case Camflow:
            speciesNames.push_back("x");
            speciesNames.push_back("T");
            speciesNames.push_back("rho");
            speciesNames.push_back("u");
            mMassFractionData = true;
            break;
        case Premix:
            speciesNames.push_back("X[cm]");
            speciesNames.push_back("T[K]");
            speciesNames.push_back("RHO[g/cm3]");
            speciesNames.push_back("V[cm/s]");
            mMassFractionData = false;
            break;
    }
            

    // Chemical species names
    while(spIt != spItEnd) {
        // Dereference unicremented value
        speciesNames.push_back((*spIt++)->Name());
    }

    // Now we know what to look for, open the file
    std::ifstream dataFile(fname.c_str());
    if(dataFile.good()) {
        //Read in the first line which should contain the column names
        std::string lineText;
        std::getline(dataFile, lineText);

        // Split out the column names into a vector, one element for each entry in this title line
        std::vector<std::string> lineEntries;
        Strings::split(lineText, lineEntries, delims); 

        // Find the column in the file which corresponds to each species from the mechanism
        const std::vector<std::string>::const_iterator columnNamesBegin = lineEntries.begin();
        const std::vector<std::string>::const_iterator columnNamesEnd   = lineEntries.end();
        
        //BOOST_FOREACH(std::string speciesName, speciesNames)
        const std::vector<std::string>::const_iterator speciesEnd = speciesNames.end();
        for(std::vector<std::string>::const_iterator speciesIt = speciesNames.begin();
            speciesIt != speciesEnd;
            ++speciesIt) {

            const std::string &speciesName = *speciesIt;

            // Get an iterator to the occurrence of the name in the column headings vector
            const std::vector<std::string>::const_iterator foundIt = std::find(columnNamesBegin, columnNamesEnd, speciesName);

            // Check that the search was successful
            if(foundIt != columnNamesEnd) {
                // Use iterator arithmetic to find the index
                speciesFileIndices.push_back(foundIt - columnNamesBegin);

                if(verbosity > 0) {
                    std::cout << "Found column for " << speciesName << " in position " << speciesFileIndices.back() << '\n';
                }
            }
            else {
                std::string msg("no column of data for ");
                msg += speciesName;
                msg += " in " + fname;
                throw std::runtime_error(msg);
            }
        }

        assert(speciesFileIndices.size() == speciesNames.size());

        // Temporary to hold data structure as it is built up
        data_collection tempData;

        // Now loop through the lines of the file reading the necessary columns
        while(!dataFile.eof()) {
            // Read in the next line
            lineText.clear();
            std::getline(dataFile, lineText);

            // Split the line into the values it contains, one per column
            lineEntries.clear();
            Strings::split(lineText, lineEntries, delims);

            // Check if this row is empty and skip it if there is nothing to do
            if(lineEntries.empty())
                continue;

            // Empty row of data for this line in the file
            data_point dataRow;

            for(size_t i = 0; i < speciesNames.size(); ++i) {
                // Read the appropriate (mass or mole fraction) floating point number from text.
                std::string fracText = lineEntries[speciesFileIndices[i]];
                real frac = atof(fracText.c_str());

                if(verbosity > 1) {
                    std::cout << speciesNames[i] << ' ' << frac << ' ';
                }

                // Unit conversions, if file not in SI units
                if(file_type == Premix) {
                    if(i == sDensityIndex)
                        // convert g cm^-3 to kg m^-3
                        frac *= 1e3;

                    if(i == sVelocityIndex)
                        // convert cm s^-1 to m s^-1
                        frac *= 1e-2;
                }

                // and assume the storage order should be the same as the species name list order
                dataRow.push_back(frac);
            }
            if(verbosity > 1) {
                std::cout << std::endl;
            }

            // Now store the data in the class member
            tempData.push_back(dataRow);
        }

        // Data read successfully so store it in the class member and finish
        mInputChemistryData.swap(tempData);
        return;
    }
    else {
        std::string msg("failed to open ");
        msg += fname;
        msg += " in PredCorrSolver::readChemistry()\n";
        throw std::runtime_error(msg);
     }
}

/*!
 * Get first value (position) of first data row.
 *
 * This will crash horribly if there is no data.
 */
real Brush::ResetChemistry::startLocation() const {
    return mInputChemistryData.front().front();
}

/*!
 * Get first value (position) of last data row.
 *
 * This will crash horribly if there is no data.
 */
real Brush::ResetChemistry::endLocation() const {
    return mInputChemistryData.back().front();
}

/*!
 * Take the stored data and overwrite the chemistry
 * data in a reactor with values taken from this class
 * interpolated as necessary.  If values are required
 * outside the range of the data the closest point will
 * be used (constant extrapolation).
 *
 * Behaviour is not defined if this ResetChemistry was
 * not constructed with an identical mechanism to that
 * of the 1d reactor.
 *
 *\param[in]        x       Position to which to interpolate the data
 *\param[in,out]    reac    Reactor whose chemistry will be replaced
 */
void Brush::ResetChemistry::apply(const real x, Mops::Reactor &reac) const {
    // Update the chemistry in each sub-reactor
    data_point dummyDataPoint(mInputChemistryData.front().size(), 0);

    // Set the position on a vector of the same type as the data points
    // so that it can be passed to the comparison function in the
    // call to lower_bound
    dummyDataPoint.front() = x;

    data_collection::const_iterator itData =
            std::lower_bound(mInputChemistryData.begin(),
                             mInputChemistryData.end(),
                             dummyDataPoint,
                             DataPointPositionComparator());

    data_point dataToUse;
    if(itData == mInputChemistryData.end()) {
        // Cell is past the end of the data so use the last data point
        dataToUse = mInputChemistryData.back();
    }
    else if(itData == mInputChemistryData.begin()) {
        // Cell is before the start of the data so use the first data point
        dataToUse = mInputChemistryData.front();
    }
    else {
        // Interpolate the data
        dataToUse = interpolate(x, *(itData - 1), *itData);
    }

    // Build a chemical mixture object
    Sprog::Thermo::IdealGas chemMixture(reac.Mech()->Species());

    // Set the species data
    fvector speciesData(dataToUse.begin() + sNumNonSpeciesData, dataToUse.end());
    if(mMassFractionData) {
        chemMixture.SetMassFracs(speciesData);
    }
    else {
        chemMixture.SetFracs(speciesData);
    }


    // Set the bulk properties
    chemMixture.SetTemperature(dataToUse[sTemperatureIndex]);
    // Mass density cannot be set until after mass fractions are set,
    // because there is an internal conversion to molar density that
    // requires the mass fractions.
    chemMixture.SetMassDensity(dataToUse[sDensityIndex]);
    chemMixture.SetVelocity(dataToUse[sVelocityIndex]);
    
    if(reac.Mixture() == NULL)
        reac.Fill(*(new Mops::Mixture(reac.Mech()->ParticleMech())));
    reac.Mixture()->SetGasPhase(chemMixture);
}

/*!
 * Linearly interpolate between two data points to get estimated data
 * for x.  The first entry in a data point is assumed to be its position.
 * The order of the data_point arguments is not important; they are named
 * left and right for clarity and if used in this way negative intermediate
 * quantities will be avoided.
 *
 *\param[in]    x           Position for which to interpolate
 *\param[in]    leftData    First point for interpolation
 *\param[in]    rightData   Second point for interpolation
 *
 *\return       Interpolated data
 */
 ResetChemistry::data_point Brush::ResetChemistry::interpolate(const real x,
                                               const data_point &leftData,
                                               const data_point &rightData) const {
     // Calculate the distance between the two data points
     const real leftPos = leftData.front();
     const real rightPos = rightData.front();
     const real distance = rightPos - leftPos;

     // And the weights to use in the interpolation
     real leftWeight  = (rightPos - x) / distance;
     real rightWeight = (x - leftPos)  / distance;

     // Now interpolate term by term
     data_point result(leftData.size());
     for(size_t i = 0; i < leftData.size(); ++i) {
         result[i] = leftWeight * leftData[i] + rightWeight * rightData[i];
     }

     return result;
 }
