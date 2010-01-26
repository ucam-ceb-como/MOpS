#!/bin/bash

#  Copyright (C) 2009 Robert I A Patterson.
#
#
# Licence:
#    This file is part of "brush".
#
#    brush is free software; you can redistribute it and/or
#    modify it under the terms of the GNU General Public License
#    as published by the Free Software Foundation; either version 2
#    of the License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#
#  Contact:
#    Prof Markus Kraft
#    Dept of Chemical Engineering
#    University of Cambridge
#    New Museums Site
#    Pembroke Street
#    Cambridge
#    CB2 3RA
#    UK
#
#    Email:       mk306@cam.ac.uk
#    Website:     http://como.cheng.cam.ac.uk

# Clean up any old files
rm regress2*.csv

echo "Test 2a: Diffusion jump process"
../bin/brush_d.x -v 2 -b ./regress2/brush2a.xml -c ./regress2/chem.inp -d ./regress2/chemsoln2a.dat -t ./regress2/therm.dat -s ./regress2/sweep2a.xml -g ./regress2/geometry.xml -a ./regress2/partsoln2a.xml
if(($? != 0))
  then
    echo "Simulation failed"
    echo "**************************"
    echo "****** TEST FAILURE ******"
    echo "**************************"
    exit $?
fi

count123=`grep "^0.2,1\.45" "regress2adiffusion123_psl.csv" | wc -l`
count124=`grep "^0.1,1\.35" "regress2adiffusion124_psl.csv" | wc -l`

if(($count123 != 112)) 
  then
    # Regression test has failed; print summary message and exit with non zero
    # value
    echo "Found $count123 particles, when 112 expected"
    echo "**************************"
    echo "****** TEST FAILURE ******"
    echo "**************************"
    exit $count123
fi

if(($count124 != 12)) 
  then
    # Regression test has failed; print summary message and exit with non zero
    # value
    echo "Found $count124 particles, when 12 expected"
    echo "**************************"
    echo "****** TEST FAILURE ******"
    echo "**************************"
    exit $i
fi

./regress2/regress2b.pl

rm regress2*.csv

# All tests passed if we get to here
echo "All tests passed"
exit 0

