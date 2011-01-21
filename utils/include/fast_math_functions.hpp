/*!
 * \file   fast_math_functions.hpp
 * \author Laurence R McGlashan
 *
 * \brief  Set of maths functions for speeding up calculations.
 * pow(x,n) for when n is an integer.
 * halley_cbrt2d(x) for a faster cube root.
 *
 Copyright (C) 2011 Laurence R McGlashan.

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

#ifndef UTILS_FAST_MATH_FUNCTIONS
#define UTILS_FAST_MATH_FUNCTIONS

//! Namespace for fast math functions.
namespace fastMath {

    //! Faster pow functions for integer powers.
    inline double pow2(const double& x) {return x*x;}
    inline double pow3(const double& x) {return x*pow2(x);}
    inline double pow4(const double& x) {return pow2(x)*pow2(x);}

    //! Iterative cube root approximation using Halley's method (double)
    inline double cbrta_halleyd(const double a, const double R)
    {
	    const double a3 = a*a*a;
        const double b= a * (a3 + R + R) / (a3 + a3 + R);
	    return b;
    }

    //! Cube root approximation using bit hack for 64-bit float adapted from Kahan's cbrt
    inline double cbrt_5d(double d)
    {
	    const unsigned int B1 = 715094163;
	    double t = 0.0;
	    unsigned int* pt = (unsigned int*) &t;
	    unsigned int* px = (unsigned int*) &d;
	    pt[1]=px[1]/3+B1;
	    return t;
    }

    //! Cube root approximation using 2 iterations of Halley's method (double)
    inline double halley_cbrt2d(double d)
    {
	    double a = cbrt_5d(d);
	    a = cbrta_halleyd(a, d);
	    return cbrta_halleyd(a, d);
    }


} // End namespace fastMathFunctions

#endif // UTILS_FAST_MATH_FUNCTIONS
