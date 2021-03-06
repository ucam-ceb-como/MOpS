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

# Clean up any outputs from previous simulations
my @outputFiles = glob("regress1a*");
if($#outputFiles > 0) {
  system("rm -f @outputFiles");
}

# Path of executable should be supplied as first argument to this script
my $program = $ARGV[0];

print "Test 1a: Nucleation without transport\n";

# Arguments for simulation
my @simulationCommand = ($program,
                         "-g", "regress1/geom.xml",
                         "-c", "regress1/species.inp",
                         "-d", "regress1/profile.dat",
                         "-s", "regress1/sweep1a.xml",
                         "-t", "regress1/therm.dat",
                         "-a", "regress1/partsoln.xml",
                         "-b", "regress1/brush1a.xml");

# Run the simulation and wait for it to finish
system(@simulationCommand) == 0 or die "ERR: simulation failed: $!";

# Collect all the moment data together
system("../../applications/solvers/brush/bin/merge-partstats.sh regress1a-Nuc") == 0 or die "ERR: failed to merge moment files: $!";

# Parse the moments file
my $momentFile;
open($momentFile, "<regress1a-NucMerged_partstats.csv") or die "ERR: failed to open merged moment file: $!";

my $m0 = 0;
my $m0var = 0;
my $m1 = 0;
my $m1var = 0;
my $count = 0;

while(<$momentFile>) {
  my @fields = split /,/;

  # Look for lines that begin with a number and have the first entry (the time)
  # equal (upto a small tolerance) to 0.1
  if(($fields[0] =~ /^\d+/) && (abs($fields[0] - 0.1) < 1e-4 )) {
      # Third field should be the zeroth moment
      $m0 += $fields[3];
      $m0var += $fields[3] * $fields[3];
      #print "$fields[3], ";

      $m1 += $fields[11];
      $m1var += $fields[11] * $fields[11];
      #print "$fields[11] \n";
      ++$count;
  }
}

# Estimate confidence intervals
$m0var = 3.29*sqrt(($m0var - $m0 * $m0 / $count)/($count - 1)/$count);
$m1var = 3.29*sqrt(($m1var - $m1 * $m1 / $count)/($count - 1)/$count);

# Get mean values
$m0 /= $count;
$m1 /= $count;

# 20 simulations (3 separate cells in each simulation) giving a total
# of 60 samples give a sample mean and 99% confidence interval for the
# population mean of
#m0: (3.010 +- 0.038) e11 cm^-3
#m1: (1.919 +- 0.024) e-10 g cm^-3
# svn r821
# The above are consistent with calculations using 20 reps with 16384 max particles (riap 01 Aug 2011), which give
#m0: (2.998+-0.001)e17 m^-3
#m1: (1.9117+-0.006)e-7 kg m^-3

print "$m0, $m0var, $m1, $m1var\n";
if(!(abs($m0 - 2.998e17) < 4e15)) {
  print "Simulated mean M0 was $m0, when 2.998e17 m^-3 expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 1;
}

if(!(abs($m1 - 1.9117e-7) < 2.4e-9)) {
  print "Simulated mean M1 was $m1, when 1.9117e-7 kg m^-3 expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 2;
}

# Clean outputs, there should always be some files to delete.
@outputFiles = glob("regress1a*");
system("rm -f @outputFiles");

#print "All tests passed\n";
exit 0;
