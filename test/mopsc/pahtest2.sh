#!/bin/bash

#  Copyright (C) 2012 Dongping Chen.
#
#
# Licence:
#    This file is part of "mops".
#
#    mops is free software; you can redistribute it and/or
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

cd pahtest2
dos2unix chem.inp
dos2unix sweep.xml
dos2unix mops.inx
dos2unix gasphase.inp
dos2unix therm.dat
dos2unix pahtest2.pl

"$program" -p --flamepp --ppah > /dev/null
#R CMD BATCH --no-save --no-restore stats.r Routput.txt &
#R2pid=$!
#wait $R2pid
perl ./pahtest2.pl 

#Capture the exit value
testresult=$?
if((testresult!=0)) 
  then
    exit $testresult
fi

rm -f pahtest2-bintree-serializer*
#rm -f stats.csv
#rm -f Routput.txt
# All tests passed
echo "All tests passed"
exit 0

