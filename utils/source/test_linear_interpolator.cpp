/*!
 * \file   test_linear_interpolator.cpp
 * \author Robert I A Patterson
 *
 * \brief  Test harness for linear interpolation utility class
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

#include "../include/linear_interpolator.hpp"

#include <iostream>
#include <utility>
#include <vector>

/*!
 * Perform a very simple test that checks constant extrapolation at both ends of
 * the data and interpolation to one point inside the data range 
 */
int main(int argc, char *argv[]) {
  std::cout << "Testing linear interpolator\n";

  std::vector<std::pair<double, float> > vector1;
  std::pair<double, float> pair1 = std::make_pair<double, float>(1.0, 10.0);
  vector1.push_back(pair1);
  pair1 = std::make_pair<double, float>(2.0, 30.0);
  vector1.push_back(pair1);
  
  Utils::LinearInterpolator<double, float> interp1(vector1);
  std::cout << "Interpolating to 0.5 should give 10.0, actual result: " << interp1.interpolate(0.5) << "\n";
  std::cout << "Interpolating to 1.5 should give 20.0, actual result: " << interp1.interpolate(1.5) << "\n";
  std::cout << "Interpolating to 5.5 should give 30.0, actual result: " << interp1.interpolate(5.5) << "\n";

  return 0;
}
