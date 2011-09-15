/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Random number generators.  Uses the Mersenne Twister generator.
    Uses several algorithms from the Numerical Recipes in C++ book.

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
    Dr Markus Kraft
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

#ifndef SWEEP_RNG_H
#define SWEEP_RNG_H

#include "swp_params.h"

namespace {
// Natural log of the gamma function.
// From Numerical Recipes, chapter 6.1.
double gammln(double xx)
{
    double x, y, tmp, ser;
    const static double cof[6] = { 76.18009172947146,
                          -86.50532032941677,
                           24.01409824083091,
                          -1.231739572450155,
                           0.1208650973866179e-2,
                          -0.5395239384953e-5};
    int j;
    y = x = xx;
    tmp = x + 5.5;
    tmp -= (x+0.5)*log(tmp);
    ser = 1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}
}

namespace Utils
{
// Returns a Poisson deviate.
// From Numerical Recipes, chapter 7.3.
template <typename U> int ignpoi(double mu, U &rand_u01)
{
    double g, sq=0.0, lnmu=0.0;
    double em=0.0, t=0.0, y=0.0;

    const double PI = 3.141592653589793238462643383279502884197169399;

    if (mu<12.0) {
        // Direct method.
        g = exp(-mu);
        em = -1;
        t = 1.0;
        do {
            ++em;
            t *= rand_u01();
        } while (t > g);
    } else {
        // Rejection method.
        sq = sqrt(2.0*mu);
        lnmu = log(mu);
        g = (mu * lnmu) - gammln(mu+1.0);
        do {
            do {
                // y = deviate from a Lorentzian comparison
                // function.
                y = tan(PI * rand_u01());
                em = (sq*y) + mu;
            } while (em < 0.0);
            em=floor(em);
            t = 0.9 * (1.0+(y*y)) * exp(em*lnmu-gammln(em+1.0)-g);
        } while (rand_u01() > t);
    }
    return (int)em;
}
}

#endif
