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

# The purpose of this test is to check that the shared pointer feature works correctly, i.e., when a doubling event occurs PAHs in a particle which are copied point to the same memory location as PAHs in the original particle
# Clean up any outputs from previous simulations
my @outputFiles = glob("pahtest4-test-sharedPointers*");
if($#outputFiles > 0) {
  print "Cleaning up old output files\n";
  system("rm " . '"' . join('" "', @outputFiles) . '"');
}

# Path of executable should be supplied as first argument to this script
my $program = $ARGV[0];

# Arguments for simulation
my @simulationCommand = ($program, "--flamepp", "-p",
                         "-g", "pahtest4/gasphase.inp",
                         "-c",  "pahtest4/chem.inp",
                         "-t",  "pahtest4/therm.dat",
                         "-s",  "pahtest4/sweep.xml",
                         "-r", "pahtest4/mops.inx");

# Run the simulation and wait for it to finish
system(@simulationCommand) == 0 or die "ERR: simulation failed: $!";

# Parse the moments file
my $momentFile;
open($momentFile, "<pahtest4-test-sharedPointers-part.csv") or die "ERR: failed to open moment file: $!";

my $m0 = 0;
my $m1 = 0;

while(<$momentFile>) {
  my @fields = split /,/;

  # Look for a line that begina with a number and has the first entry (the time)
  # equal (upto a small tolerance) to 0.030816
  if(($fields[0] =~ /^\d+/) && (abs($fields[1] - 0.030816) < 1e-6 )) {
      # Third field should be the zeroth moment
      $m0 = $fields[4];
      #print "4: $fields[4], ";

      $m1 = $fields[16];
      #print "16: $fields[16] \n";

      last;
  }
}

# Precalculated value: M0=2.19e+20, Fv=1.51e-6

# 20 repetitions
# mean values and 99% confidence interval widths
# m0 (4.23+-1.82)e+20 m^-3
# fv (1.59+-0.12)e-6

print "$m0, $m1\n";
if(abs($m0 - 4.23e+20) > 2.05e+20) {
  print "Simulated mean M0 was $m0, when  4.23e+20m^-3 expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 1;
}

if(abs($m1 - 1.59e-6) > 0.34e-6) {
  print "Simulated mean Fv was $m1, when 1.59e-6 expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 2;
}

#print "All tests passed\n";
system("rm pahtest4-test-sharedPointers*");
exit 0;
