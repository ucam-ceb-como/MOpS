/*!
 * \file   cam_sdr.h
 * \author L. R. McGlashan
 *
 * \brief Class for calculating the scalar dissipation rate.
 *
 *  Copyright (C) 2011 L. R. McGlashan.
 *

 Licence:
    This file is part of "camflow".

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

#ifndef _SDR_H
#define	_SDR_H

#include "array.h"
#include "cam_math.h"

namespace Camflow
{

/*! @brief SDR class.
 *
*/
class ScalarDissipationRate
{

        enum sdr
        {
            constant_fromStrainRate,
            profile_fromStrainRate,
            constant_fromCFD,
            profile_fromCFD
        };

        int sdrType_;

        //! The first index is 'time coordinate'.
        //! The second index is 'mixture fraction coordinate'.
        Array2D scalarDissipationRate_;
        doublereal scalarDissipationRateRef_;

    public:

        //! Default constructor.
        ScalarDissipationRate
        (
            const int n_TimePoints,
            const int n_MixtureFractionPoints
        );

        //! Destructor.
        ~ScalarDissipationRate();


        const doublereal calculate
        (
            const doublereal mixtureFraction,
            const doublereal strainRate
        );


        inline const doublereal getScalarDissipationRate(const int mixFracCoord)
        {
            return scalarDissipationRate_(0,mixFracCoord);
        };

        void setSdrType(const int sdrType) {sdrType_ = sdrType;};

        const doublereal getRefSDR() const
        {
            return scalarDissipationRateRef_;
        }

        /*inline const std::vector<doublereal> getScalarDissipationRate()
        {
            return scalarDissipationRate_(0);
        }*/

}; // End ScalarDissipationRate class declaration.

} // End Camflow namespace.

#endif	/* _SDR_H */
