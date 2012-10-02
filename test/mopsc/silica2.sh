#!/bin/bash

#  Copyright (C) 2011 William J Menz.
#
#   This script tests the Silica model (developed by Shraddha and Markus
#   Sander) with the weighted transition kernl (w1 weightrule)
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

cd silica2
rm -f silica*

"$program" -p -strang -rr mops-fm.inx
fmResult=$?

if(($fmResult==0)) 
then
  echo "Finished simulation"
else
  echo "FM Simulation failed"
  exit $fmResult
fi
echo "========================"

"$program" -p -strang -rr mops-sf.inx
sfResult=$?

if(($sfResult==0)) 
then
  echo "Finished simulation"
else
  echo "SF Simulation failed"
  exit $sfResult
fi
echo "========================"

dos2unix "silica-fm-part.csv"
dos2unix "silica-sf-part.csv"
perl test-silica2.pl
postprocessResult=$?

if(($postprocessResult!=0)) 
  then
    exit $postprocessResult
fi

rm -f silica*

exit 0
