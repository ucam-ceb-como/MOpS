/*!
 * \file   first_of_pair_comparator.hpp
 * \author Robert I A Patterson
 *
 * \brief  Functor for comparing std::pair<T, X> based on the first element
 *
 Copyright (C) 2010 Robert I A Patterson.

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
#ifndef UTILS_FIRST_OF_PAIR_COMPARATOR
#define UTILS_FIRST_OF_PAIR_COMPARATOR

#include <utility>
#include <limits>

namespace Utils {
/*!
 *@tparam		X		Real number type
 *@tparam		T		Any type
 *
 * This class is intended to facilitate the use of standard algorithms with
 * sorted sequences of pair<T, X> where ordering is based only on the first
 * element of the pair.  It is further assumed here that the type, T, of the
 * first element is a floating point type.  The class is based on that
 * described in Item 23 of "Effective STL: 50 Specific Ways to Improve Your
 * Use of the Standard Template Library" by Scott Meyers.
 *
 * All the logic is in the operator()(const X, constX) method.  The remaining methods
 * provide interfaces simply unpack pairs that may be provided by the
 * standard algorithms.
 */
template<typename X, typename T> class FirstOfPairComparator {
public:
    //! The type of one value in the sequence
    typedef std::pair<X, T> pair_type;

    /*!
     * @param[in]   lhs     First  pair
     * @param[in]   rhs     Second pair
     *
     * @return      True iff the first element of the first pair sorts strictly
     *              before the first element of the second pair.
     */
    bool operator()(const pair_type &lhs, const pair_type &rhs) const {
        return operator()(lhs.first, rhs.first);
    }

    /*!
     * @param[in]   lhs     First  pair
     * @param[in]   k2      Second key value
     *
     * @return      True iff the first element of the first pair sorts strictly
     *              before the second key value.
     */
    bool operator()(const pair_type &lhs, const X &k2) const {
        return operator()(lhs.first, k2);
    }

    /*!
     * @param[in]   k1      First  key value
     * @param[in]   rhs     Second pair
     *
     * @return      True iff the first key value sorts strictly
     *              before the first element of the second pair.
     */
    bool operator()(const X &k1, const pair_type &rhs) const {
        return operator()(k1, rhs.first);
    }

    /*!
     * @param[in]   k1      First  key value
     * @param[in]   k2      Second key value
     *
     * @return      True iff the first key value is less than second by a relative margin of 10 epsilon
     *
     * The 10 epsilon tolerance is introduced to hide very small floating point precision differences,
     * however, it does mean that this method may return a value != (k1 < k2)
     */
    bool operator()(const X &k1, const X &k2) const {
        //std::cout << k1 << ',' << k2 << ' '
        //          << (k1 < k2 * (static_cast<X>(1.0) - 10 * std::numeric_limits<X>::epsilon()))
        //          << '\n';
        return k1 < k2 * (static_cast<X>(1.0) - 10 * std::numeric_limits<X>::epsilon());
    }
};

} //namespace Utils
#endif //UTILS_FIRST_OF_PAIR_COMPARATOR
