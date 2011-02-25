/*!
 * \file   choose_index.hpp
 * \author Robert I A Patterson
 *
 * \brief  Utility method to choose an index based on a vector of weights
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
#ifndef UTILS_CHOOSE_INDEX
#define UTILS_CHOOSE_INDEX

#include <vector>
#include <numeric>
#include <algorithm>

/*!
 *@tparam		T		Real number type
 *@param[in]		weights		Vector of weights (all >= 0)
 *@param[in]		rng		Function that returns U[0,1] deviates
 *
 *@return		Index of an entry in weights with probability proportional to that entry
 *
 * The return value is \f$ i \in \left[0, \mathrm{weights.size\left(\right)}\right), \mathrm{P}\left(i = j\right) \propto \mathrm{weights\left[i\right]}\f$

 *
 */
template<typename T> int chooseIndex(const std::vector<T> &weights, T (*rng)()) {
    // partialSums[i] will contain the sum of the first i+1 elements of weights
    std::vector<T> partialSums(weights.size());
    std::partial_sum(weights.begin(), weights.end(), partialSums.begin());


    // Multiply a U[0,1] variable by the sum of the elements in weights
    T r = rng() * partialSums.back();

    // Perform an inverse transform to sample an index from weights with
    // probability of choosing i proportional to weights[i]
    return std::distance(partialSums.begin(), std::lower_bound(partialSums.begin(), partialSums.end(), r));
}

#endif //UTILS_CHOOSE_INDEX
