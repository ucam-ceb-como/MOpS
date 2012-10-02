#!/bin/bash

#  Copyright (C) 2011 William J Menz.
#
#   This script tests the finite-rate sintering in the surfvol model
#   All three cases use the viscous flow model with varying sintering parameters:
#       1. finite sintering (26 nm particles)
#       2. instantaneous sintering (42 nm particles)
#       3. no sintering (0.49 nm particles)
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

#Path to executable should be supplied as first argument to
#this script.  Script will fail and return a non-zero value
#if no executable specified.
program=$1

if test -z "$program"
  then
    echo "No executable supplied to $0"
    exit 255
fi

# An optional second argument may specify the working directory
if test -n "$2"
  then
    cd "$2"
    echo "changed directory to $2"
fi

cd sinter1

"$program" -p -strang -rr mops-finite.inx -s sweep-finite.xml
r1=$?
"$program" -p -strang -rr mops-spherical.inx -s sweep-spherical.xml
r2=$?
"$program" -p -strang -rr mops-nosinter.inx -s sweep-nosinter.xml
r3=$?

if((r1==0 && r2==0 && r3==0))
then
  echo "Finished simulation"
else
  echo "Simulation failed"
  cd ..
  exit 2
fi
echo "========================"

dos2unix "silane-finite-psl(0.08s).csv"
dos2unix "silane-nosinter-psl(0.08s).csv"
dos2unix "silane-spherical-psl(0.08s).csv"
tclsh sintertest.tcl
postprocessResult=$?
if((postprocessResult!=0)) 
  then
    cd ..
    exit $postprocessResult
fi

rm -f silane*

cd ..
exit 0

