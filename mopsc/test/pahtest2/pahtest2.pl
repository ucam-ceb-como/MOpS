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
#system("rm pahtest2-*");

# See if this is a windows system
my $windows = ($ENV{'OS'} =~ /windows.*/i);

# Choose the windows executable name if appropriate
my $program = "../../bin/mops_d.x";
if($windows) {
    $program = "../../bin/mops_d.exe";
}

# Arguments for simulation
my @simulationCommand = ($program, "-flamepp", "-p",);

# Run the simulation and wait for it to finish
system(@simulationCommand) == 0 or die "ERR: simulation failed: $!";

# Parse the moments file
my $momentFile;
open($momentFile, "<pahtest2-part.csv") or die "ERR: failed to open moment file: $!";

my $m0 = 0;
my $secondary_m0 = 0;
my $m1 = 0;

while(<$momentFile>) {
  my @fields = split /,/;

  # Look for a line that begina with a number and has the first entry (the time)
  # equal (upto a small tolerance) to 0.04
  if(($fields[0] =~ /^\d+/) && (abs($fields[1] - 0.02051) < 1e-6 )) {
      # Third field should be the zeroth moment
      $m0 = $fields[4];
      #print "4: $fields[4], ";

      $m1 = $fields[16];
      #print "16: $fields[16] \n";

  }
  elsif (($fields[0] =~ /^\d+/) && (abs($fields[1] - 0.0035) < 1e-6 )) {
      $secondary_m0 = $fields[26];
      #print "26: $fields[26] \n";
  }
}


print "$m0, $m1, $secondary_m0\n";

# Using 10 runs with 1024 main computational particles get
# mean m0 of 9.76e11, with 99% confidence interval for this mean of +-4.3e10
if(abs($m0 -  9.76e11) > 4.3e10) {
  print "Simulated mean M0 was $m0, when 9.76e11 cm^-3 expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 1;
}

# Using 10 runs with 1024 main computational particles get
# mean fv of 3.13e-8, with 99% confidence interval for this mean of +-5e-10
if(abs($m1 - 3.13e-8) > 5e-10) {
  print "Simulated mean M1 was $m1, when 3.13e-8 g cm^-3 expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 2;
}

# Using 10 runs with 1024 main computational particles get
# mean secondary m0 of 2.85e11, with 99% confidence interval 
# for this mean of +-8e9.  However, 2 runs with 512 particles
# gives a much larger noise.
if(abs($secondary_m0 -  2.85e11) > 1.5e10) {
  print "Simulated mean secondary M0 was $secondary_m0, when 2.85e11 cm^-3 expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 3;
}

#print "All tests passed\n";

exit 0;
