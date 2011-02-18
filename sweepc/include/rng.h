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

namespace Sweep
{
    // Returns a Poisson deviate.
    // From Numerical Recipes, chapter 7.3.
    int ignpoi(real mu, real(*rand_u01)());

    // Generates binomial deviates.
    // From Numerical Recipes, chapter 7.3.
    int ignbin(
        int n,  // Number of trials.
        real pp, // Probability for accepting single trial.
        real(*rand_u01)() // Uniform [0,1] generator
        );

    //! Box-Muller generation of Normal sample
    double randNorm(const double mean, const double stddev, real(*rand_u01)());
};

#endif
