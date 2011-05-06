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
my @outputFiles = glob("pahtest2output*");
if($#outputFiles > 0) {
  print "Cleaning up old output files\n";
  system("rm " . '"' . join('" "', @outputFiles) . '"');
}

# Path of executable should be supplied as first argument to this script
my $program = $ARGV[0];

# Arguments for simulation
my @simulationCommand = ($program, "-flamepp", "-p",
                         "-gp", "gasphase.inp",
                         "-c",  "chem.inp",
                         "-t",  "therm.dat",
                         "-s",  "sweep.xml",
                         "-rr", "mops.inx");


# Run the simulation and wait for it to finish
system(@simulationCommand) == 0 or die "ERR: simulation failed: $!";

# Parse the moments file
my $momentFile;
open($momentFile, "<pahtest2output-part.csv") or die "ERR: failed to open moment file: $!";

my $m0 = 0;
my $secondary_m0 = 0;
my $m1 = 0;

while(<$momentFile>) {
  my @fields = split /,/;

  # Look for a line that begins with a number and has the first entry (the time)
  # equal (upto a small tolerance) to 0.04
  if(($fields[0] =~ /^\d+/) && (abs($fields[1] - 0.02051) < 1e-6 )) {
      # Third field should be the zeroth moment
      $m0 = $fields[4];
      #print "4: $fields[4], ";

      $m1 = $fields[16];
      #print "16: $fields[16] \n";

  }
  elsif (($fields[0] =~ /^\d+/) && (abs($fields[1] - 0.0035) < 1e-6 )) {
      $secondary_m0 = $fields[34];
      #print "30: $fields[34] \n";
  }
}


print "$m0, $m1, $secondary_m0\n";

# Using 20 runs with 2048 main computational particles get
# mean m0 of 9.806e11, with 99% confidence interval for this mean of +-1.80e10
# and Fv = 3.122e-8 +- 2.3e-10.
# For secondary particles get mean m0 2.851e11 with a confidence interval
# of +-2.2e9
# svn r821
# Note that the noise in the test results will be considerably bigger than the 
# confidence intervals reported above, because fewer samples and fewer computational
# particles are used.
# The standard deviations of the distributions from which the regression test samples
# twice (in that it does two runs) are approximately
# m0: 3e10, Fv: 1e-9, secondary m0: 9e9
# These values were estimated from 10 runs using ordinary test settings and svn r821,
# they are standard deviations for one sample, not for the mean of a number of samples.
# All data quoted in this comment is in cgs units.

if(abs($m0 -  9.81e17) > 5.4e16) {
  print "Simulated mean M0 was $m0, when 9.81e17 m^-3 expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 1;
}

if(abs($m1 - 3.12e-8) > 1e-9) {
  print "Simulated mean Fv was $m1, when 3.12e-8 expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 2;
}

if(abs($secondary_m0 -  2.85e17) > 1.5e16) {
  print "Simulated mean secondary M0 was $secondary_m0, when 2.85e17 m^-3 expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 3;
}

#print "All tests passed\n";

# Clean up output files, since the test passed
@outputFiles = glob("pahtest2output*");
system("rm " . '"' . join('" "', @outputFiles) . '"');

exit 0;
