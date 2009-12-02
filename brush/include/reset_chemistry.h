/*!
 * \file   reset_chemistry.h
 * \author Robert I A Patterson
 *
 * \brief Interpolated chemical species amounts onto a 1d reactor
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
#ifndef BRUSH_RESET_CHEMISTRY_H
#define	BRUSH_RESET_CHEMISTRY_H

#include "brush_params.h"

#include <vector>

// Forward declaration
namespace Mops {
    class Reactor;
}

namespace Brush{

// Forward declarations
class Reactor1d;
    
//! Hold a spatially resolved chemistry profile and interpolate onto 1d reactors.
/*!
 * Time varying data is not currently supported, only spatial variation is possible.
 * Data for each interpolation node is held in an object of type data_point and
 * when an interpolation is requested these are searched using the functor defined
 * within this class to find the two nearest points for linear interpolation.  All
 * instances of data_point (which is a typedef for an std::vector of real numbers
 * should be of the same length), but there is no checking of this at present.
 *
 * The comparison functor performs real number ordering using only the first entry
 * of the data_point vectors.  The first few entries in these vectors contain bulk
 * property information (see the constructor documentation for details).  The remaining
 * enties are species mass fractions.
 */
class ResetChemistry {
public:
    //! Format of input file
    enum InputFileType {
        //! Mass fractions and SI units
        Camflow,

        //! Mole fractions and cgs units
        Premix,
    };

    //! Read in the data needed for the given mechanism
    ResetChemistry(const string &fname, const InputFileType file_type, const Sprog::Mechanism& mech, const int verbosity);

    //! Use data from vectors
    ResetChemistry(const fvector &x, const fvector &Temp, const fvector &rho,
                   const fvector &u, const std::vector<fvector> &massFracs);

    //! Overwrite chemistry information with that stored in this object
    void apply(const real x, Mops::Reactor &reac) const;
    
    //! Position of first data point
    real startLocation() const;

    //! Position of final data point
    real endLocation() const;

private:
    //! Type of one row of data
    typedef fvector data_point;

    //! Type of full set of data
    typedef std::vector<data_point> data_collection;

    //! One row per input point, x ordinate is first item in each inner vector
    data_collection mInputChemistryData;

    //! Functor defining spatial ordering of data points
    struct DataPointPositionComparator {
        /*!
         * The position to which a data point applies is assumed to be given
         * by the first element of the vector representing the data point.
         *
         *\param[in]    a   First data point to compare
         *\param[in]    b   Second data point to compare
         *
         *\return       True iff the position of a is strictly before the position of b
         */
        bool operator()(const data_point &a, const data_point &b) {
            return (a.front() < b.front());
        }
    };

    //! Linearly interpolate between the two data points
    data_point interpolate(const real x, const data_point &leftData,
                           const data_point &rightData) const;

    //! Number of data items that are not species mass fractions
    static const size_t sNumNonSpeciesData;

    //! Is the data stored in this instance stored as mass fractions or molefractions
    bool mMassFractionData;

    //! Index of temperature data in data_point
    static const size_t sPositionIndex;

    //! Index of temperature data in data_point
    static const size_t sTemperatureIndex;

    //! Index of density data in data_point
    static const size_t sDensityIndex;

    //! Index of velocity data in data_point
    static const size_t sVelocityIndex;

};

}

#endif	/* BRUSH_RESET_CHEMISTRY_H */

