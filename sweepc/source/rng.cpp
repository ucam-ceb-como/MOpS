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


#include "rng.h"

// Anonymous namespace for function only used in this file
namespace {
    // Natural log of the gamma function.
    // From Numerical Recipes, chapter 6.1.
    Sweep::real gammln(Sweep::real xx)
    {
        Sweep::real x, y, tmp, ser;
        static Sweep::real cof[6] = { 76.18009172947146,
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

    // Returns a Poisson deviate.
    // From Numerical Recipes, chapter 7.3.
    int Sweep::ignpoi(real mu, real(*rand_u01)())
    {
        static real oldmu=(-1.0), g, sq=0.0, lnmu=0.0;
        real em=0.0, t=0.0, y=0.0;

        if (mu<12.0) {
            // Direct method.
            if (mu != oldmu) {
                oldmu = mu;
                g = exp(-mu);
            }
            em = -1;
            t = 1.0;
            do {
                ++em;
                t *= (*rand_u01)();
            } while (t > g);
        } else {
            // Rejection method.
            if (mu!=oldmu) {
                oldmu = mu;
                sq = sqrt(2.0*mu);
                lnmu = log(mu);
                g = (mu * lnmu) - gammln(mu+1.0);
            }
            do {
                do {
                    // y = deviate from a Lorentzian comparison
                    // function.
                    y = tan(PI*(*rand_u01)());
                    em = (sq*y) + mu;
                } while (em < 0.0);
                em=floor(em);
                t = 0.9 * (1.0+(y*y)) * exp(em*lnmu-gammln(em+1.0)-g);
            } while ((*rand_u01)() > t);
        }
        return (int)em;
    }

    // Generates binomial deviates.
    // From Numerical Recipes, chapter 7.3.
    int Sweep::ignbin(
        int n,  // Number of trials.
        real pp, // Probability for accepting single trial.
        real(*rand_u01)() // Uniform [0,1] generator
        )
    {
        int j;
        static int nold=(-1);
        real am, em, g, angle, p, bnl, sq, t, y;
        static real pold=(-1.0), pc, plog, pclog, en, oldg;

        p = (pp <= 0.5 ? pp : 1.0-pp);

        // Calculate mean of deviate produced.
        am = n * p;

        if (n < 25) {
            // Use direct method for small n.
            bnl = 0.0;
            // At most 25 calls to (*rand_u01)().
            for (j=0; j<n; ++j) {
                if ((*rand_u01)() < p) ++bnl;
            }
        } else if (am < 1.0) {
            // Fewer than 1 event is expected and there are more than
            // 25 trials, then can assume a Poisson distribution.  Use
            // the direct Poisson method.
            g = exp(-am);
            t = 1.0;
            for (j=0; j<=n; ++j) {
                t *= (*rand_u01)();
                if (t < g) break;
            }
            bnl = (j<=n ? j : n);
        } else {
            // Use the rejection method.
            if (n != nold) {
                // if n has changed then calculate useful
                // quantities.
                en   = n;
                oldg = gammln(en + 1.0);
                nold = n;
            }
            if (p != pold) {
                // If p has changed then compute useful
                // quantities.
                pc    = 1.0 - p;
                plog  = log(p);
                pclog = log(pc);
                pold  = p;
            }
            // Now use rejection method with a Lorentzian
            // comparison function.
            sq = sqrt(2.0 * am * pc);
            do {
                do {
                    angle = PI * (*rand_u01)();
                    y = tan(angle);
                    em = sq * y + am;
                } while (em < 0.0 || em >= (en+1.0)); // Reject.
                em = floor(em); // Trick for integer-valued distribution.
                t = 1.2 * sq * (1.0+y*y) * exp(oldg - gammln(em+1.0) - 
                                               gammln(en-em+1.0) + 
                                               em*plog +(en-em)*pclog);
            } while ((*rand_u01)() > t); // Reject.  This occurs about 1.5 times per deviate.
            bnl = em;
        }
        if (p != pp) bnl = n - bnl; // Remember to undo symmetry transformation.
        return (int)bnl;
    }

    /*!
     * Return a real number from a normal (Gaussian) distribution with given
     * mean and standard deviation by polar form of Box-Muller transformation
     *
     *@param[in]    mean    Mean of Normal distribution
     *@param[in]    stddev  Standard deviation of Normal distribution
     *@param[in,out]    rand_u01    Pointer to function that generates U[0,1] deviates
     *
     *@return       Sample from N(mean, stddev)
     */
     double Sweep::randNorm(const double mean, const double stddev, real(*rand_u01)()) {
        double x, y, r;
        do
        {
                x = 2.0 * (*rand_u01)() - 1.0;
                y = 2.0 * (*rand_u01)() - 1.0;
                r = x * x + y * y;
        }
        while ( r >= 1.0 || r == 0.0 );
        double s = std::sqrt( -2.0 * std::log(r) / r );
        return mean + x * s * stddev;
    }
