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

print "Test 3b: Advection, diffusion and thermophoresis";

# Clean up any outputs from previous simulations
my @outputFiles = glob("regress3a*");
if($#outputFiles > 0) {
  system("rm @outputFiles");
}

# Path of executable should be supplied as first argument to this script
my $program = $ARGV[0];


# Arguments for simulation
my @simulationCommand = ($program,
                         "-g", "regress3/geometry.xml",
                         "-c", "regress3/chem.inp",
                         "-d", "regress3/chemsoln3b.dat",
                         "-r", "regress3/tran.dat",
                         "-s", "regress3/sweep3b.xml",
                         "-t", "regress3/therm.dat",
                         "-b", "regress3/brush3b.xml",
                         "-a", "regress3/partsoln3.xml",
                         "--premix-chem");

# Run the simulation and wait for it to finish
system(@simulationCommand) == 0 or die "ERR: simulation failed: $!";

# Collect all the moment data together
system("../../applications/solvers/brush/bin/merge-partstats.sh regress3b-adv-diffn-thermoph") == 0 or die "ERR: failed to merge moment files: $!";

# Parse the moments file
my $momentFile;
open($momentFile, "<regress3b-adv-diffn-thermophMerged_partstats.csv") or die "ERR: failed to open merged moment file: $!";

my $m0a = 0;
my $counta = 0;

my $m0b = 0;
my $countb = 0;

#This test should exactly match the results of 3a, because it has
#adds a constant thermophoresis velocity, which is exactly offset by an increase
#in the bulk velocity.  As a result the random number sequences should be the
#same, although this might be affected by some floating point quirks resulting
#from the switch from Camflow chemistry input format (chemsoln3a.dat) to
#PREMIX format (chemsoln3b.dat).

while(<$momentFile>) {
  my @fields = split /,/;

  if(($fields[0] =~ /^0\.2$/) && (abs($fields[1] - 0.017) < 1e-6)) {
      $m0a += $fields[3];
      #print "t=0.2 x=0.017 $fields[3]\n";
      $counta++;
  }

  if(($fields[0] =~ /^0\.3$/) && (abs($fields[1] - 0.0245) < 1e-6)) {
      $m0b += $fields[3];
      #print "t=0.3 x=0.0245 $fields[3]\n";
      $countb++;
  }
}

# take the mean of the m0 values
$m0a = $m0a / $counta;
$m0b = $m0b / $countb;

print "$m0a\n";
if(abs($m0a - 1.7e2) > 2e1) {
  print "Simulated M0 at t=0.2 near x=0.017 was $m0a, when analytic solution is 1.7e2\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 1;
}

print "$m0b\n";
if(abs($m0b - 2.6e2) > 3e1) {
  print "Simulated M0 at t=0.3 near x=0.0245 was $m0b, when analytic solution is 2.6e2\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 2;
}


# Clean outputs, there should always be some files to delete.
@outputFiles = glob("regress3b*");
system("rm @outputFiles");


#print "All tests passed\n";
exit 0;
