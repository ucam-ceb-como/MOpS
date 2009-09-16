/*!
 * \file   linear_interpolator.hpp
 * \author Robert I A Patterson
 *
 * \brief  Utility method for linear interpolation
 *
 Copyright (C) 2009 Robert I A Patterson.

 Licence:
 
    This utility file is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This file is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have file a copy of the GNU General Public License
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

#ifndef UTILS_LINEAR_INTERPOLATOR
#define UTILS_LINEAR_INTERPOLATOR

#include <vector>
#include <utility>
#include <algorithm>
#include <cassert>

namespace Utils {

/*!
 * \brief Simple linear interpolation class
 *
 *@tparam	X 	Floating point type.
 *@tparam       Y	Type supporting multiplication by type X and addition of another instance of Y, effectively a vector field over X
 *
 * Data is stored internally as pair<X, Y> where X represents the position of
 * the data and Y is the data value.
 */
template<typename X, typename Y> class LinearInterpolator {
public:
    //! Intialise with the given data
    LinearInterpolator(const std::vector<std::pair<X, Y> > &data);

    //! Get an interpolated value
    const Y interpolate(const X& x) const;

private:
    //! Type of internal data structure
    typedef std::vector<std::pair<X, Y> > data_vector_type;

    //! The data to use for interpolation
    data_vector_type mData;

    //! Comparison functor for use in std algorithms
    class LessThan {
    public:
        /*!
         * Comparison is done on the value of the first element of the pairs
         *
         *\param[in]    lhs    First pair for comparison
         *\param[in]    rhs    Second pair for comparison
         *
         *\return       lhs.first() < rhs.first()
         */
        bool operator()(const std::pair<X, Y> &lhs, const std::pair<X, Y> &rhs) {
            return (lhs.first < rhs.first);
        }
    };
};


} //namespace Utils


/*!
 * Sort the provided data if necessary and use it to initialise the interpolator
 *
 *@tparam	X 	Floating point type.
 *@tparam       Y	Type supporting multiplication by type X and addition of another instance of Y, effectively a vector field over X
 *
 *\param[in]    data    Pairs of position value information 
 */
template<typename X, typename Y> Utils::LinearInterpolator<X, Y>::LinearInterpolator(const std::vector<std::pair<X, Y> > &data)
: mData(data)
{
    std::sort(mData.begin(), mData.end(), LessThan());
}

/*!
 * Calculate a value of Y linearly interpolate to x.  For points outside the range of the
 * data the nearest data point will be used (constant extrapolation).
 *
 *@tparam	X 	Floating point type.
 *@tparam       Y	Type supporting multiplication by type X and addition of another instance of Y, effectively a vector field over X
 *
 *\param[in]    x    Position for which to interpolate
 *
 *\return	Interpolated value
 */
template<typename X, typename Y> const Y Utils::LinearInterpolator<X, Y>::interpolate(const X &x) const {

    Y result;
    // Need a pair<X, Y> with x as its first element to use in the lookup, the contents of the second element
    // of the pair is ignored.  Note this is a call to the constructor of std::pair, not the evaluation of 
    // some other function.
    std::pair<X, Y> lookup(x, result);
    typename data_vector_type::const_iterator itLower = std::lower_bound(mData.begin(), mData.end(), lookup, LessThan());

    if(itLower == mData.end()) {
    // At or beyond the upper end of the data
        return mData.back().second;
    }
    else if(itLower == mData.begin()) {
    // At or before the lower end of the data
        return mData.front().second;
    }
    else {
    // Linear interpolation
        const X xLeft = itLower->first;
        const X xRight = (itLower - 1)->first;
        const X distance = xRight - xLeft;
        
        // Calculate the interpolation weights
        const X leftWeight = (xRight - x) / distance;
        const X rightWeight = (x - xLeft) / distance;

        // Do the interpolation
        result = leftWeight * itLower->second + rightWeight * (itLower - 1)->second;
    }
    return result;
}

#endif //UTILS_LINEAR_INTERPOLATOR
