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
#include "swp_cell.h"

#include <fstream>
#include <vector>
#include <string>
#include <ostream>
#include <algorithm>
#include <stdexcept>
#include <cstdlib>
#include <cassert>
#include <cctype>
#include <iomanip>

#include <boost/algorithm/string.hpp>

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
const size_t ResetChemistry::sNumNonSpeciesData = 9;

/*!
 * Index of spatial position data in the elements of mInputChemistryData.
 * The data will always have units m.
 */
const size_t ResetChemistry::sPositionIndex = 0;

/*!
 * Index of temperature data in the elements of mInputChemistryData.
 * The data will always have units K.
 */
const size_t ResetChemistry::sTemperatureIndex = 1;

/*!
 * Index of density data in the elements of mInputChemistryData.
 * The data will always be have units \f$\mathrm{kg\,m^{-3}}\f$
 */
const size_t ResetChemistry::sDensityIndex = 2;

/*!
 * Index of velocity data in the elements of mInputChemistryData.
 * The data will always have units \f$\mathrm{m\,s^{-1}}\f$
 */
const size_t ResetChemistry::sVelocityIndex = 3;

/*!
 * Rate of PAH formation for use with PAH-PP model
 * OR
 * Factor used in ABF surface reaction rates for soot, but
 * could potentially be used to scale other rates
 */
const size_t ResetChemistry::sPAHFormationIndex = 4;

/*!
 * Index of
 * \f[ \frac{\partial T}{\partial x}, \f]
 * for use in flamelet calculations.
 */
const size_t ResetChemistry::sGradientTemperatureIndex = 5;

/*!
 * Index of mixture fraction diffusion coefficient (only for use
 * in flamelet calculations).
 * The data will have units
 * \f$\mathrm{m^2\,s^{-1}}\f$
 */
const size_t ResetChemistry::sMixFracDiffCoeffIndex = 6;

/*!
 * Index of spatial gradient (in physical space) of mixture
 * fraction, that is
 * \f[ \frac{\partial Z}{\partial x} \f],
 * for use in flamelet calculations.
 */
const size_t ResetChemistry::sGradientMixFracIndex = 7;

/*!
 * Index of Laplacian (in physical space) of mixture
 * fraction, that is
 * \f[ \frac{\partial^2 Z}{\partial x^2} \f],
 * for use in flamelet calculations.
 */
const size_t ResetChemistry::sLaplacianMixFracIndex = 8;

/*!
 * Factor used in ABF surface reaction rates for soot, but
 * could potentially be used to scale other rates.
 */


/**
 * Construct an object from a data file and associate the concentration data
 * with the indices defined by the mechanism.
 *
 * The data file should contain one header row, which is case sensitive followed
 * by complete rows of numerical data, completely blank lines are permitted and
 * ignored.  The columns may be separated by spaces or tabs.
 *
 * Mandatory columns in the file for the Camflow format (\ref InputFileType) are
 * - x Spatial position to which the row of data applies (\f$\mathrm{m}\f$)
 * - T Temperature (K)
 * - rho Mass density of the gas mixture (\f$\mathrm{kg\,m^{-3}}\f$)
 * - u Velocity of gas mixture (\f$\mathrm{m\,s^{-1}}\f$)
 * - A column of mass fraction data must be provided for each species in the
 *   mechanism and the column headings for these species must be identical to the
 *   strings specified as species names when the mechanism was constructed.
 *
 * In addition to the above mentioned columns the following are mandatory when
 * the Camflow Flamelet format (\ref InputFileType) is used.  The u column is
 * still required, but will generally contain zeros.  The wdotA4 is provided
 * for use with the PAH-PP soot particle model and arbitrary values may be used
 * with other soot models.
 *  - The x column heading is replaced by Z
 *  - wdotA4 The formation rate of pyrene (\f$\mathrm{mol\,m^{-3}\,s^{-1}}\f$)
 *  - D_z The mixture fraction diffusion coefficient  (\f$\mathrm{m^2\,s^{-1}}\f$)
 *  - GradZ Gradient in real space of the mixture fraction (\f$\mathrm{m^{-1}}\f$)
 *  - LaplZ Laplacian in real space of the mixture fraction (\f$\mathrm{m^{-2}}\f$)
 *  - GradT Gradient in real space of the temperature (\f$\mathrm{K\,m^{-1}}\f$)
 *
 * Mandatory columns in the file are (Premix format \ref InputFileType)
 * - x Spatial position to which the row of data applies (\f$\mathrm{cm}\f$)
 * - T Temperature (K)
 * - rho Mass density of the gas mixture (\f$\mathrm{g\,cm^{-3}}\f$)
 * - u Velocity of gas mixture (\f$\mathrm{cm\,s^{-1}}\f$)
 * - wdotA4 The formation rate of pyrene (\f$\mathrm{mol\,cm^{-3}\,s^{-1}}\f$)
 * - GradT Gradient in real space of the temperature (\f$\mathrm{K\,cm^{-1}}\f$)
 * - A column of mole fraction data must be provided for each species in the
 *   mechanism and the column headings for these species must be identical to the
 *   strings specified as species names when the mechanism was constructed.
 *
 * The PremixAlpha file format is the same as the premix format, except that
 * the requirement for a wdotA4 column is replaced by the requirement for an
 * Alpha column.
 *
 *\param[in]    fname       Path of file from which to read the data
 *\param[in]    file_type   Style of data in file
 *\param[in]    mech        Mechanism defining species indexes
 *\param[in]    verbosity   Level of debugging information sent to standard out (higher is more)
 *
 *\exception    std::runtime_error  No data for a species present in the mechanism
 */
Brush::ResetChemistry::ResetChemistry(const std::string &fname, const InputFileType file_type,
                                      const Sprog::Mechanism& mech, const int verbosity) {
    // Delimeters to use when splitting lines of data from file into the individual column entries
    const std::string delims(" ,\t");

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
            speciesNames.push_back("wdotA4");
            speciesNames.push_back("GradT");
            mMassFractionData = false;
            break;
        case PremixAlpha:
            speciesNames.push_back("X[cm]");
            speciesNames.push_back("T[K]");
            speciesNames.push_back("RHO[g/cm3]");
            speciesNames.push_back("V[cm/s]");
            speciesNames.push_back("Alpha");
            speciesNames.push_back("GradT[K/cm]");
            mMassFractionData = false;
            break;
        case CamflowFlamelet:
            speciesNames.push_back("Z");
            speciesNames.push_back("T");
            speciesNames.push_back("rho");
            speciesNames.push_back("u");
            speciesNames.push_back("wdotA4");
            speciesNames.push_back("GradT");
            speciesNames.push_back("D_z");
            speciesNames.push_back("GradZ");
            speciesNames.push_back("LaplZ");
            mMassFractionData = true;
            break;
    }
            


    // Iterators to the chemical species in the mechanism
    Sprog::Mechanism::const_sp_iterator spIt = mech.SpBegin();
    const Sprog::Mechanism::const_sp_iterator spItEnd = mech.SpEnd();
    while(spIt != spItEnd) {
        // Dereference unincremented value
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
        boost::algorithm::split(lineEntries, lineText, boost::algorithm::is_any_of(delims), boost::algorithm::token_compress_on);


        // Find the column in the file which corresponds to each species from the mechanism
        const std::vector<std::string>::const_iterator columnNamesBegin = lineEntries.begin();
        const std::vector<std::string>::const_iterator columnNamesEnd   = lineEntries.end();

        if(verbosity > 2) {
            std::cerr << "Read following columns of chemistry data from " << fname << '\n';
            std::vector<std::string>::const_iterator it = columnNamesBegin;
            while(it != columnNamesEnd-1){
                std::cerr << *it << ':';
                ++it;
            }
            std::cerr << "\nEnd of chemistry data columns" << std::endl;
        }
        
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

            //skip lines that do not contain anything (other than perhaps some form of line end character
            if(lineText.length() < 2)
                continue;

            // Split the line into the values it contains, one per column
            lineEntries.clear();
            boost::algorithm::split(lineEntries, lineText, boost::algorithm::is_any_of(delims), boost::algorithm::token_compress_on);

            // Check if this row is empty and skip it if there is nothing to do
            // Try to handle muddled line end sequences
            if(lineEntries.empty() || 
               ((lineEntries.size() == 1) &&
                ((lineEntries.at(0) == "\r") || (lineEntries.at(0) == "\n"))))
                continue;

            // Empty row of data for this line in the file, with one entry for
            // each item of non-species data.  The unknown number of species
            // entries will be added to the end of the vector.
            data_point dataRow(sNumNonSpeciesData);

            for(size_t i = 0; i < speciesNames.size(); ++i) {
                // Read the appropriate (mass or mole fraction) floating point number from text.
                //std::cerr << i << ' ' << speciesFileIndices[i] << ' ' << speciesNames[i] << ' ' << lineEntries[speciesFileIndices[i]].length() << std::endl;
                std::string fracText = lineEntries[speciesFileIndices[i]];
                real frac = atof(fracText.c_str());

                if(verbosity > 1) {
                    std::cout << speciesNames[i] << ' ' << frac << ' ';
                }

                // Unit conversions, if file not in SI units
                if((file_type == Premix) || (file_type == PremixAlpha)) {
                    if(i == sPositionIndex)
                        // convert cm to m
                        frac *= 1e-2;

                    if(i == sDensityIndex)
                        // convert g cm^-3 to kg m^-3
                        frac *= 1e3;

                    if(i == sVelocityIndex)
                        // convert cm s^-1 to m s^-1
                        frac *= 1e-2;

                    // The PAHFormation column is used for alpha values in the premix alpha input format
                    // and because alpha values are dimensionless, they should not be rescaled
                    if((i == sPAHFormationIndex) && (file_type != PremixAlpha))
                        // convert cm^-3 to m^-3
                        frac *= 1e6;

                    if(i == sGradientTemperatureIndex)
                        // convert K cm^-1 to K m^-1
                        frac *= 1e2;

                }

                // Space has already been allocated for the non-species data
                // The input file should not contain columns for PAH formation
                // and flamelet related quantities unless it is a flamelet
                // input file.
                if((i <= sVelocityIndex) ||
                    ((file_type == CamflowFlamelet) && (i < sNumNonSpeciesData)) ||
                    (((file_type == Premix) || (file_type == PremixAlpha)) && (i <= sGradientTemperatureIndex))) {
                    dataRow[i] = frac;
                }
                else {
                    // and assume the storage order should be the same as the species name list order
                    dataRow.push_back(frac);
                }
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
 * Method for use when coupling to flamelet codes.
 * 
 *@param[in]    x           Positions to which data apply (\f$\mathrm{m}\f$)
 *@param[in]    Temp        Temperature data (K)
 *@param[in]    rho         Mixture density data (\f$\mathrm{kg\,m^{-3}}\f$)
 *@param[in]    u           Bulk gas velocity data (\f$\mathrm{m\,s^{-1}}\f$)
 *@param[in]    PAH         PAH formation rate
 *@param[in]    D_Z         Mixture fraction diffusion coefficient
 *@param[in]    grad_Z      Gradient of mixture fraction in physical space
 *@param[in]    lapl_Z      Laplacian of mixture fraction in physical space
 *@param[in]    grad_T      Gradient of temperature
 *@param[in]    massFracs   Vector of vectors of mass fractions
 *
 * The first nine arguments must all have the same length, which must also be
 * the length of the elements of massFracs.  This is because they are the
 * columns of a rectangular 2d array of data, where one data set is one row.
 *
 * It might make sense to implement this constructor in terms of ranges rather than
 * containers.
 *
 *@exception     std::invalid_argument   Input vectors have differing lengths
 */
ResetChemistry::ResetChemistry(const fvector &x, const fvector &Temp,
                               const fvector &rho, const fvector &u,
                               const fvector &PAH, const fvector &D_Z,
                               const fvector &grad_Z, const fvector &lapl_Z,
                               const fvector &grad_T,
                               const std::vector<fvector> &massFracs) {
    //===== Check all input vectors have the same length ===========
    const size_t len = x.size();
    if(len != Temp.size()) {
        std::ostringstream msg;
        msg << "Length of temperature vector is " << Temp.size()
            << ", but it must match the x vector length of " << len
            << " (ResetChemistry::ResetChemistry)";
        throw std::invalid_argument(msg.str());
    }

    if(len != rho.size()) {
        std::ostringstream msg;
        msg << "Length of density vector is " << rho.size()
            << ", but it must match the x vector length of " << len
            << " (ResetChemistry::ResetChemistry)";
        throw std::invalid_argument(msg.str());
    }

    if(len != u.size()) {
        std::ostringstream msg;
        msg << "Length of velocity vector is " << u.size()
            << ", but it must match the x vector length of " << len
            << " (ResetChemistry::ResetChemistry)";
        throw std::invalid_argument(msg.str());
    }

    if(len != PAH.size()) {
        std::ostringstream msg;
        msg << "Length of PAH production rate vector is " << PAH.size()
            << ", but it must match the x vector length of " << len
            << " (ResetChemistry::ResetChemistry)";
        throw std::invalid_argument(msg.str());
    }

    if(len != D_Z.size()) {
        std::ostringstream msg;
        msg << "Length of mixture fraction diffusion coefficient vector is "
            << D_Z.size()
            << ", but it must match the x vector length of " << len
            << " (ResetChemistry::ResetChemistry)";
        throw std::invalid_argument(msg.str());
    }

    if(len != grad_Z.size()) {
        std::ostringstream msg;
        msg << "Length of mixture fraction gradient vector is " << grad_Z.size()
            << ", but it must match the x vector length of " << len
            << " (ResetChemistry::ResetChemistry)";
        throw std::invalid_argument(msg.str());
    }

    if(len != lapl_Z.size()) {
        std::ostringstream msg;
        msg << "Length of mixture fraction Laplacian vector is " << lapl_Z.size()
            << ", but it must match the x vector length of " << len
            << " (ResetChemistry::ResetChemistry)";
        throw std::invalid_argument(msg.str());
    }

    if(len != grad_T.size()) {
        std::ostringstream msg;
        msg << "Length of gradient of temperature vector is "
            << grad_T.size()
            << ", but it must match the x vector length of " << len
            << " (ResetChemistry::ResetChemistry)";
        throw std::invalid_argument(msg.str());
    }

    const std::vector<fvector>::const_iterator itEnd = massFracs.end();
    for(std::vector<fvector>::const_iterator it = massFracs.begin();
        it != itEnd; ++it) {
            if(len != it->size()) {
                std::ostringstream msg;
                msg << "Length of a mass fraction vector "
                    << std::distance(massFracs.begin(), it)
                    << " is " << it->size()
                    << ", but it must match the x vector length of " << len
                    << " (ResetChemistry::ResetChemistry)";
                throw std::invalid_argument(msg.str());
            }
    }

    //Input data vectors all have same length, if this point is reached
    //================================================================

    // Temporary to hold data structure as it is built up
    data_collection tempData;

    for(size_t i = 0; i != len; ++i) {
        // Build up row of data
        data_point dataRow(massFracs.size() + sNumNonSpeciesData);
        dataRow[sPositionIndex] = x[i];
        dataRow[sTemperatureIndex] = Temp[i];
        dataRow[sDensityIndex] = rho[i];
        dataRow[sVelocityIndex] = u[i];
        dataRow[sPAHFormationIndex] = PAH[i];
        dataRow[sMixFracDiffCoeffIndex] = D_Z[i];
        dataRow[sGradientMixFracIndex] = grad_Z[i];
        dataRow[sLaplacianMixFracIndex] = lapl_Z[i];
        dataRow[sGradientTemperatureIndex] = grad_T[i];

        // copy mass fraction data for position x[i]
        const std::vector<fvector>::const_iterator itFracEnd = massFracs.end();

        // non species data has been copied above
        data_point::iterator itTarget = dataRow.begin() + sNumNonSpeciesData;

        //std::cout << "Mass fracs from vector:";
        for(std::vector<fvector>::const_iterator itFrac = massFracs.begin();
            itFrac != itFracEnd; ++itFrac, ++itTarget) {
                *itTarget = itFrac->operator[](i);
                //std::cout << ' ' << *itTarget;
            }
        //std::cout << std::endl;

        // Now store the row that has just been built up
        tempData.push_back(dataRow);
    }

    // Data read successfully so store it in the class member and finish
    mInputChemistryData.swap(tempData);
    mMassFractionData = true;
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
void Brush::ResetChemistry::apply(const real x, Sweep::Cell &reac) const {
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
    Sprog::Thermo::IdealGas chemMixture(*(reac.GasPhase().Species()));

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

    // Other properties that have been addded to the mixture
    chemMixture.SetVelocity(dataToUse[sVelocityIndex]);
    chemMixture.SetPAHFormationRate(dataToUse[sPAHFormationIndex]);
    chemMixture.SetMixFracDiffCoeff(dataToUse[sMixFracDiffCoeffIndex]);
    chemMixture.SetGradientMixFrac(dataToUse[sGradientMixFracIndex]);
    chemMixture.SetLaplacianMixFrac(dataToUse[sLaplacianMixFracIndex]);
    chemMixture.SetGradientTemperature(dataToUse[sGradientTemperatureIndex]);
    // PAH formation is reused for ABF alpha since they can never appear together and
    // saves re-engineering this function (riap2 27 Oct 2011)
    chemMixture.SetAlpha(dataToUse[sPAHFormationIndex]);

    reac.SetGasPhase(chemMixture);

    // Uncomment this code to check molar concentrations
    /*unsigned int logIndices[14];
    for(unsigned int i = 0; i < reac.ParticleModel()->Species()->size(); ++i) {
        if(reac.ParticleModel()->Species()->at(i)->Name() == "H2")
            logIndices[0] = i;
        else if(reac.ParticleModel()->Species()->at(i)->Name() == "H")
            logIndices[1] = i;
        else if(reac.ParticleModel()->Species()->at(i)->Name() == "O")
            logIndices[2] = i;
        else if(reac.ParticleModel()->Species()->at(i)->Name() == "O2")
            logIndices[3] = i;
        else if(reac.ParticleModel()->Species()->at(i)->Name() == "OH")
            logIndices[4] = i;
        else if(reac.ParticleModel()->Species()->at(i)->Name() == "H2O")
            logIndices[5] = i;
        else if(reac.ParticleModel()->Species()->at(i)->Name() == "CH4")
            logIndices[6] = i;
        else if(reac.ParticleModel()->Species()->at(i)->Name() == "CO")
            logIndices[7] = i;
        else if(reac.ParticleModel()->Species()->at(i)->Name() == "CO2")
            logIndices[8] = i;
        else if(reac.ParticleModel()->Species()->at(i)->Name() == "C2H2")
            logIndices[9] = i;
        else if(reac.ParticleModel()->Species()->at(i)->Name() == "C2H4")
            logIndices[10] = i;
        else if(reac.ParticleModel()->Species()->at(i)->Name() == "AR")
            logIndices[11] = i;
        else if(reac.ParticleModel()->Species()->at(i)->Name() == "N2")
            logIndices[12] = i;
        else if(reac.ParticleModel()->Species()->at(i)->Name() == "A4")
            logIndices[13] = i;

    }

    // Output in PSDF_input.dat format
    std::cout << x << '\t' << std::scientific << std::setprecision(6)
              << reac.GasPhase().Temperature() << '\t';
    for(unsigned int j = 0; j != 14; ++j) {
        std::cout <<  reac.GasPhase().MolarConc(logIndices[j]) << '\t';
    }
    std::cout << reac.GasPhase().Pressure() << '\t'
              << reac.GasPhase().Density() << '\n';*/
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
