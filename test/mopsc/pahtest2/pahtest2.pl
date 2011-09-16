#!/usr/bin/perl

#  Copyright (C) 2011 Dongping Chen
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

use strict;
use warnings;
# this test is designed for testing doubling algorithm after coupled with KMC model
# Clean up any outputs from previous simulations
my @outputFiles = glob("pahtest2-test-doubling-algorithm*");
if($#outputFiles > 0) {
  print "Cleaning up old output files\n";
  system("rm " . '"' . join('" "', @outputFiles) . '"');
}

# Path of executable should be supplied as first argument to this script
my $program = $ARGV[0];

# Arguments for simulation
my @simulationCommand = ($program, "-flamepp", "-p",
                         "-gp", "pahtest2/gasphase.inp",
                         "-c",  "pahtest2/chem.inp",
                         "-t",  "pahtest2/therm.dat",
                         "-s",  "pahtest2/sweep.xml",
                         "-rr", "pahtest2/mops.inx");

# Run the simulation and wait for it to finish
system(@simulationCommand) == 0 or die "ERR: simulation failed: $!";

# Parse the moments file
my $momentFile;
open($momentFile, "<pahtest2-test-doubling-algorithm-part.csv") or die "ERR: failed to open moment file: $!";

my $sp = 0;
my $m0 = 0;
my $m1 = 0;

while(<$momentFile>) {
  my @fields = split /,/;

  # Look for a line that begina with a number and has the first entry (the time)
  # equal (upto a small tolerance) to 0.0058
  if(($fields[0] =~ /^\d+/) && (abs($fields[1] - 0.0058) < 1e-6 )) {
      # Third field should be the zeroth moment
      $sp = $fields[2];
      #print "2: $fields[2], ";

      $m0 = $fields[4];
      #print "4: $fields[4] \n";

      $m1 = $fields[20];
      last;
  }
}

# Precalsulated value: num of computational particles (sp) =855, M0=71.25e17+-1e15

# git 00706668ed.. gives the following values for the mean 
# and 99% confidence interval for the mean, using 20 repeitions
# boost 1.42.0
# sp = 854+-24.5
# m0 = (1.251+-0.036)e17 m^-3
# fv = (3.846+-0.076)e-7 kg m^-3
# with boost 1.47.0 git 58056fdbc (which should be the same as 00706668ed for
# PAH-PP purposes) gives the following mean values
# sp 861.7
# m0 1.262e17 m^-3
# fv 3.872e-7

print "$sp, $m0, $m1\n";
if(abs($sp -  855) > 20) {
  print "Simulated sp was $sp, when  855 expected\n";
  print "if pahtest1 passes and this test fails, it will indicate that the doubling algorithm works in a wrong way.";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 1;
}

if(abs($m0 - 1.25e17) > 4e15) {
  print "Simulated mean M0 was $m0, when  1.25e17m^-3 expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 2;
}

if(abs($m1 - 3.85e-7) > 1e-8) {
  print "Simulated mean M1 was $m1, when  3.85e-7 kg m^-3 expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 2;
}

#print "All tests passed\n";
system("rm pahtest2-test-doubling-algorithm*");
system("rm DIMER.csv");
system("rm MONOMER.csv");
exit 0;
