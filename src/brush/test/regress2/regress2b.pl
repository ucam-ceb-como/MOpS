#!/usr/bin/perl

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

use strict;
use warnings;

print "Test 2b: Advection jump process (deprecated)\n";

# See if this is a windows system
my $windows = ($ENV{'OS'} =~ /windows.*/i);

# Choose the windows executable name if appropriate
my $program = "../bin/brush_d.x";
if($windows) {
    $program = "../bin/brush_d.exe";
}

# Arguments for simulation
my @simulationCommand = ($program,
                         "-g", "regress2/geometry.xml",
                         "-c", "regress2/chem.inp",
                         "-d", "regress2/chemsoln2b.dat",
                         "-s", "regress2/sweep2b.xml",
                         "-t", "regress2/therm.dat",
                         "-b", "regress2/brush2b.xml",
                         "-a", "regress2/partsoln2b.xml");

# Run the simulation and wait for it to finish
system(@simulationCommand) == 0 or die "ERR: simulation failed: $!";

# Collect all the moment data together
system("../bin/merge-partstats.sh regress2badvection") == 0 or die "ERR: failed to merge moment files: $!";

# Parse the moments file
my $momentFile;
open($momentFile, "<regress2badvectionMerged_partstats.csv") or die "ERR: failed to open merged moment file: $!";

my $m0 = 0;

while(<$momentFile>) {
  my @fields = split /,/;

  if(($fields[0] =~ /^0\.3/) && ((abs($fields[1] - 1.35) < 1e-6 ) || (abs($fields[1] - 1.25) < 1e-6) || (abs($fields[1] - 1.45) < 1e-6))) {
      $m0 += $fields[3];
      #print "$fields[3], ";
  }
}

print "$m0\n";
if($m0 < 6.1e-7) {
  print "Simulated M0 near x=1.35 was $m0, when at least 6.1e-7 cm^-3 expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 1;
}


#print "All tests passed\n";
exit 0;
