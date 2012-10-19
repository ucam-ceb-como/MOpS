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
#include "linear_interpolator.hpp"

namespace Camflow
{

/*! @brief SDR class.
 *
*/
class ScalarDissipationRate
{

    enum sdr
    {
        NONE,
        notFromCFD,
        constant_fromCFD,
        profile_fromCFD
    };

    int sdrType_;

    const double& stoichZ_;
    double strainRate_;
    double stoichSDR_;

    const std::vector<double>& mixFracCoords_;

    //! Stoichiometric scalar dissipation rate that has a time history.
    std::vector<double> v_stoichSDR;

    //! Time profile of the scalar dissipation rates.
    std::vector<double> v_time;

    //! SDR profile with a time history.
    //! The second index is 'time coordinate'.
    //! The first index is 'mixture fraction coordinate'.
    std::vector< std::vector<double> > scalarDissipationRate_;

    Utils::LinearInterpolator<double, double>*  interpolator_;
    Utils::LinearInterpolator<double, double>*  interpolatorZeroTime_;
    Utils::LinearInterpolator<double, double>*  interpolatorNextTime_;

    void readStrainRate(const std::string& inputFileName);
    double calculate(const double& mixtureFraction) const;
    double scalarDissipationRate(const double& mixtureFraction) const;
    double strainRate(const double& mixtureFraction) const;

public:

    //! Default constructor.
    ScalarDissipationRate
    (
        const std::string& inputFileName,
        const double& stoichZ,
        const std::vector<double>& mixFracCoords,
        const int n_TimePoints
    );

    //! Destructor.
    ~ScalarDissipationRate();

    void setStrainRate(const double strainRate);

    void setSDRRate(const double sdr);

    void setExternalScalarDissipationRate
    (
        const std::vector<double>& time,
        const std::vector<double>& sdr
    );

    double operator()(const double& Z, const double& time) const;

    inline const double& getStoichSDR() const
    {
        return stoichSDR_;
    }

}; // End ScalarDissipationRate class declaration.

} // End Camflow namespace.

#endif	/* _SDR_H */
