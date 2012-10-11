#!/usr/bin/perl

#  Copyright (C) 2012 Robert I A Patterson.
#
#
# Licence:
#    This file is part of "camflow".
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

use strict;
use warnings;

# Open the newly simulated solution file
my $profileFile;
open($profileFile, "<profile.dat") or die "ERR: failed to open profile.dat file: $!";

# Open the file containing the reference solution
my $originalProfileFile;
open($originalProfileFile, "<profileOriginal.dat") or die "ERR: failed to open profileOriginal.dat file: $!";

# Variable to hold the lines as they are read from the files
my @oldSolution;
my @newSolution;

# Loop through and keep the last line, which should be the converged solution
while(<$profileFile>) {
    @newSolution = split(" +");
}
#print "New solution @newSolution->[0][1] \n";

# Take the converged solution from the reference calculation too
while(<$originalProfileFile>) {
    @oldSolution = split(" +");
}
#print "Old solution @oldSolution\n";

#Compare temperature in old and new output files, skip the first line (line 0)
# because it contains the column headings.
for(my $i=1;$i<13;++$i) 
{
    #print  " @oldSolution->[$j][$i] and @newSolution->[$j][$i] \n";
    if (abs($oldSolution[$i] - $newSolution[$i]) > 1e-6 * abs($oldSolution[$i])) 
    {
      print "Old value was $oldSolution[$i], but new value is $newSolution[$i]\n";
      print "This was in column $i.\n"; 
      print "**************************\n";
      print "****** TEST FAILURE ******\n";
      print "**************************\n";
      exit 1;
    }
}
exit 0;
