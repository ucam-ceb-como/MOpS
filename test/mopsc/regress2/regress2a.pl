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
my @outputFiles = glob("regression2a-nuc-coag-pyr*");
if($#outputFiles > 0) {
  system("rm", "-f", @outputFiles);
}

# Path of executable should be supplied as first argument to this script
my $program = $ARGV[0];

# Arguments for simulation
my @simulationCommand = ($program, "--flamepp", "-p",
                         "-g", "regress2/regress2.inp",
                         "-s", "regress2/regress2a.xml",
			"-c", "regress2/chem.inp",
			"-t", "regress2/therm.dat",
                         "-r", "regress2/regress2a.inx");

# Run the simulation and wait for it to finish
system(@simulationCommand) == 0 or die "ERR: simulation failed: $!";

# Parse the moments file
my $momentFile;
open($momentFile, "<regression2a-nuc-coag-pyr-part.csv") or die "ERR: failed to open moment file: $!";

my $m0 = 0;
my $m1 = 0;

while(<$momentFile>) {
  my @fields = split /,/;

  # Look for a line that begina with a number and has the first entry (the time)
  # equal (upto a small tolerance) to 0.04
  if(($fields[0] =~ /^\d+/) && (abs($fields[1] - 0.04) < 1e-6 )) {
      # Third field should be the zeroth moment
      $m0 = $fields[4];
      #print "4: $fields[4], ";

      $m1 = $fields[16];
      #print "16: $fields[16] \n";

      last;
  }
}

# 20 runs using 2048 computational particles gives the following estimates and
# 99% confidence intervals for the estimates
# m0 2.94+-0.06e10 cm^-3
# fv 1.26+-0.003e-11

# With git 0d294425... 200 repetitions with other settings unchanged gives
# m0 (2.936+-0.037)e17 m^-3
# fv (1.260+-0.007)e-8

print "$m0, $m1\n";
if(abs($m0 - 2.94e17) > 1.5e16) {
  print "Simulated mean M0 was $m0, when 2.94e17 m^-3 expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 1;
}

if(abs($m1 - 1.26e-8) > 4e-10) {
  print "Simulated mean Fv was $m1, when 1.26e-8 expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 2;
}

# Clean outputs, there should always be some files to delete.
@outputFiles = glob("regression2a-nuc-coag-pyr*");
system("rm", "-f", @outputFiles);

#print "All tests passed\n";
exit 0;
