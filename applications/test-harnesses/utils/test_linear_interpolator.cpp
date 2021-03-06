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

#include "linear_interpolator.hpp"

#include <iostream>
#include <utility>
#include <vector>
#include <cmath>

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
  double result = interp1.interpolate(0.5);
  std::cout << "Interpolating to 0.5 should give 10.0, actual result: "
            << result << "\n";
  if(std::abs(result - 10.0) > 1e-8)
    return 1;

  result = interp1.interpolate(1.5);
  std::cout << "Interpolating to 1.5 should give 20.0, actual result: "
            << result << "\n";
  if(std::abs(result - 20.0) > 2e-8)
    return 2;

  result = interp1.interpolate(5.5);
  std::cout << "Interpolating to 5.5 should give 30.0, actual result: "
            << result << "\n";
  if(std::abs(result - 30.0) > 3e-8)
    return 3;

  std::cout << "Testing linear interpolator with vector inputs for position and data values\n";

  std::vector<double> vector3;
  std::vector<float> vector4;
  vector3.push_back(1.0);
  vector3.push_back(2.0);
  vector4.push_back(10.0);
  vector4.push_back(30.0);

  Utils::LinearInterpolator<double, float> interp2(vector3,vector4);

  result = interp2.interpolate(0.5);
  std::cout << "Interpolating to 0.5 should give 10.0, actual result: "
            << result << "\n";
  if(std::abs(result - 10.0) > 1e-8)
    return 4;

  result = interp2.interpolate(1.5);
  std::cout << "Interpolating to 1.5 should give 20.0, actual result: "
            << result << "\n";
  if(std::abs(result - 20.0) > 2e-8)
    return 5;

  result = interp2.interpolate(5.5);
  std::cout << "Interpolating to 5.5 should give 30.0, actual result: "
            << result << "\n";
  if(std::abs(result - 30.0) > 3e-8)
    return 6;

  // If we got this far then all the tests must have passed.
  return 0;
}
