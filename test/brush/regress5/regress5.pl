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
my @outputFiles = glob("regress5-addw*");
if($#outputFiles > 0) {
  system("rm", "-f", @outputFiles);
}

# Path of executable should be supplied as first argument to this script
my $program = $ARGV[0];

print "Test 5: Inception on inflow boundary with coagulation\n";

# Arguments for simulation
my @simulationCommand = ($program,
                         "-g", "regress5/geom.xml",
                         "-c", "regress5/species.inp",
                         "-d", "regress5/chemsoln.dat",
                         "-s", "regress5/sweep.xml",
                         "-t", "regress5/therm.dat",
                         "-a", "regress5/partsoln.xml",
                         "-b", "regress5/brush.xml");

# Run the simulation and wait for it to finish
system(@simulationCommand) == 0 or die "ERR: simulation failed: $!";

# Collect all the moment data together
system("../../applications/solvers/brush/bin/merge-partstats.sh regress5-addw") == 0 or die "ERR: failed to merge moment files: $!";

# Parse the moments file
my $momentFile;
open($momentFile, "<regress5-addwMerged_partstats.csv") or die "ERR: failed to open merged moment file: $!";

my $m0 = 0;
my $m0var = 0;
my $m2 = 0;
my $m2var = 0;
my $count = 0;

while(<$momentFile>) {
  my @fields = split /,/;

  # Look for lines that begin with a number and have the first entry (the time)
  # greater than 1.2 so that equilibrium should have been reached and select
  # the cell from 0.9 to 1.0 (centre at 0.95).
  if(($fields[0] =~ /^\d+/) && ($fields[0] > 1.2) && (abs($fields[1]-0.95) < 0.01)) {
      # Third field should be the zeroth moment
      $m0 += $fields[3];
      $m0var += $fields[3] * $fields[3];
      #print "$fields[3], ";

      $m2 += $fields[13];
      $m2var += $fields[13] * $fields[13];
      #print "$fields[13] \n";
      ++$count;
  }
}

# Estimate confidence intervals
$m0var = 3.29*sqrt(($m0var - $m0 * $m0 / $count)/($count - 1)/$count);
$m2var = 3.29*sqrt(($m2var - $m2 * $m2 / $count)/($count - 1)/$count);

# Get mean values
$m0 /= $count;
$m2 /= $count;

# Analytic solution is
#m0: 0.4904 m^-3
#m1: 1.000 kg m^-3 (not tested)
#m2: 4.158 kg^2 m^-3

# 10 simulations giving a total of 160 samples with 
# max 16384 computational particles
# give a sample mean and 99% confidence interval for the
# population mean of
#m0: (0.4845 +- 0.0015) m^-3
#m2: (4.124  +- 0.027) kg^2 m^-6
# The actual test settings (2 runs max 1024 particles, but slightly longer
# run times to give a total of 40 samples) should
# have confidence intervals sqrt(64) times larger.

print "$m0, $m0var, $m2, $m2var, ($count)\n";
if(abs($m0 - 0.4904) > 0.01) {
  print "Simulated mean M0 was $m0, when 0.4904 m^-3 expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 1;
}

if(abs($m2 - 4.158) > 0.20) {
  print "Simulated mean M2 was $m2, when 4.158 kg^-2 m^-6 expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 2;
}

# Clean outputs, there should always be some files to delete.
@outputFiles = glob("regress5-addw*");
system("rm", "-f", @outputFiles);

#print "All tests passed\n";
exit 0;
