#!/usr/bin/perl

#  Copyright (C) 2011 Robert I A Patterson.
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

# Clean up any outputs from previous simulations
my @outputFiles = glob("moonmd1results*");
if($#outputFiles > 0) {
  system("rm @outputFiles");
}

# Path of executable should be supplied as first argument to this script
my $program = $ARGV[0];

print "Test moonmd1\n";

# Arguments for simulation
my @simulationCommand = ($program,
                         "moonmd1/chem.inp",
                         "moonmd1/therm.dat",
                         "moonmd1/brush.xml",
                         "moonmd1/sweep.xml",
                         "moonmd1/partsoln.xml",
                         "moonmd1/chemsoln.dat");

# Run the simulation and wait for it to finish
system(@simulationCommand) == 0 or die "ERR: simulation failed: $!";

# Parse the moments file
my $momentFile;
open($momentFile, "<moonmd1results-moments.csv") or die "ERR: failed to open moment file: $!";

my $m0 = 0;
my $m0sq = 0;
my $m1 = 0;
my $m1sq = 0;
my $avgcoag = 0;
my $maxcoag = 0;
my $count = 0;

while(<$momentFile>) {
  my @fields = split /,/;

  # Look for lines that begin with a number and have the second entry (the position)
  # equal (upto a small tolerance) to 0.71
  if(($fields[0] =~ /^\d+/) && (abs($fields[1] - 0.71) < 1e-4 )) {
      # Third field should be the zeroth moment
      $m0 += $fields[3];
      $m0sq += $fields[3] * $fields[3];
      #print "$fields[3], ";

      $m1 += $fields[11];
      $m1sq += $fields[11] * $fields[11];
      #print "$fields[11] \n";

      $avgcoag += $fields[15];
      $maxcoag = ($maxcoag > $fields[16]) ? $maxcoag : $fields[16];

      ++$count;
  }
}

# Get mean values
$m0 /= $count;
$m0sq = 1.96 * sqrt(($m0sq / $count - $m0 * $m0) / $count);
$m1 /= $count;
$m1sq = 1.96 * sqrt(($m1sq / $count - $m1 * $m1) / $count);
$avgcoag /= $count;

# Analytic solns, const coag and inception
# 0 initial and inflow conditions
# velocity 1
# Inception rate I = 5.31e9 m^-3 s^-1
# Coagulation kernel K = 5e-9 m^3 s^-1
# Velocity u = 2 m s^-1
# Mass of incepted particle m = 2 kg
# m0(x) = sqrt(2I/K)tanh(sqrt(IK/2)*x/u)
# m1(x) = Ixm/u
#
# Analytical solutions at x=0.71 are:
# m0(0.71) = 1.25e9
# m1(0.71) = 3.77e9
#
# Ten repetitions with git a8fcecd2f...
# gives the following mean and standard
# deviations for the results:
# m0: (1.26+-0.05)e9
# m1: (3.77+-0.12)e9
#
# avg and max coag are very rough checks that the code does not
# change unexpectedly.  Correct values are not known although maxcoag
# should be an integer.

printf "%G %G %G %G \n", $m0, $m0sq, $m1, $m1sq;
if(!(abs($m0 - 1.25e9) < 1e8)) {
  print "Simulated mean M0 was $m0, when 1.25e9 m^-3 expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 1;
}

if(!(abs($m1 - 3.77e9) < 2e8)) {
  print "Simulated mean M1 was $m1, when 3.77e9 kg m^-3 expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 2;
}

if(($avgcoag < 1e-4) || ($avgcoag > 0.1)) {
  print "Average number of coagulations per cell has changed\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 3;
}


if(($maxcoag != 1) && ($maxcoag != 2)) {
  print "Maximum number of coagulations per cell has changed\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 4;
}

# Clean outputs, there should always be some files to delete.
@outputFiles = glob("moonmd1results*");
#print "Files to remove: @outputFiles\n";
system("rm @outputFiles");

#print "All tests passed\n";
exit 0;
