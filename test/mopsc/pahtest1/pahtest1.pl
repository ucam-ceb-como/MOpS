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

# this test is designed for checking whether PAH-PP model work correctly by comparing to precalculated value
# Clean up any outputs from previous simulations
my @outputFiles = glob("pahtest1-test-pah-kmc-only*");
if($#outputFiles > 0) {
  print "Cleaning up old output files\n";
  system("rm -f" . '"' . join('" "', @outputFiles) . '"');
}

# Path of executable should be supplied as first argument to this script
my $program = $ARGV[0];

# Arguments for simulation
my @simulationCommand = ($program, "--flamepp", "-p",
                         "-g", "pahtest1/gasphase.inp",
                         "-c",  "pahtest1/chem.inp",
                         "-t",  "pahtest1/therm.dat",
                         "-s",  "pahtest1/sweep.xml",
                         "-r", "pahtest1/mops.inx");

# Run the simulation and wait for it to finish
system(@simulationCommand) == 0 or die "ERR: simulation failed: $!";

# Parse the moments file
my $momentFile;
open($momentFile, "<pahtest1-test-pah-kmc-only-part.csv") or die "ERR: failed to open moment file: $!";

my $m0 = 0;
my $m1 = 0;

while(<$momentFile>) {
  my @fields = split /,/;

  # Look for a line that begina with a number and has the first entry (the time)
  # equal (upto a small tolerance) to 0.006
  if(($fields[0] =~ /^\d+/) && (abs($fields[1] - 0.006) < 1e-6 )) {
      # Third field should be the zeroth moment
      $m0 = $fields[4];
      #print "4: $fields[4], ";

      $m1 = $fields[16];
      #print "16: $fields[16] \n";

      last;
  }
}

# Precalculated value: M0 = 2.87e+19, Fv = 5.74e-9

# 20 repetitions
# mean values and 99% confidence interval widths
# m0 (2.87+-0.097)e+19 m^-3
# fv (5.74+-0.125)e-9

print "$m0, $m1\n";
if(abs($m0 -  2.87e19) > 9.70e17) {
  print "Simulated mean M0 was $m0, when  2.87e19m^-3 expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 1;
}

if(abs($m1 - 5.74e-9) > 1.25e-10) {
  print "Simulated mean Fv was $m1, when 5.74e-9 expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 2;
}

#print "All tests passed\n";
system("rm -f pahtest1-test-pah-kmc-only*");
system("rm -f DIMER.csv");
system("rm -f MONOMER.csv");
exit 0;
