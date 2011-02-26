#!/bin/bash

#  Copyright (C) 2010 Rebecca C Riehl.
#
#
# Licence:
#    This file is part of "mops".
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

cd cstrtest

#Choose the windows or linux names for the executable
uname -s | grep --ignore-case CYGWIN 
if(($?==0))
then
	program="../../../bin/release/mops.exe"
else 
	program="../../../bin/release/mops" 
fi

dos2unix chem.inp

$program -gpc -p -diag4 -rr mops.inx -s sweep.xml -c chem.inp -t therm.dat

if(($?==0)) 
then
  echo "Finished simulation"
else
  echo "Simulation failed"
  exit $?
fi
echo "========================"

dos2unix "J3-chem.csv"
dos2unix "J3-sensi.csv"
./cstrtest.pl
if(($?!=0)) 
  then
    cd ..
    exit $?
fi


rm -f J3*
# All tests passed
#echo "All tests passed"



exit 0

