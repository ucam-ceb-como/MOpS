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

print "Diffusion using Li Wang drag coefficients\n";

#This is diffusion and advection without any particle processes.
#Diffusion coefficients are 2.833e-5 for particles of size 10
# and 2.082e-5 for particles of size 15.

# Clean up any outputs from previous simulations
my @outputFiles = glob("regress3c*");
if($#outputFiles > 0) {
  system("rm @outputFiles");
}

# Path of executable should be supplied as first argument to this script
my $program = $ARGV[0];

# Arguments for simulation
my @simulationCommand = ($program,
                         "-g", "regress3/geometry.xml",
                         "-c", "regress3/chem.inp",
                         "-d", "regress3/chemsoln3c.dat",
                         "-s", "regress3/sweep3c.xml",
                         "-t", "regress3/therm.dat",
                         "-b", "regress3/brush3c.xml",
                         "-a", "regress3/partsoln3.xml");

# Run the simulation and wait for it to finish
system(@simulationCommand) == 0 or die "ERR: simulation failed: $!";

# Collect all the moment data together
system("../../applications/solvers/brush/bin/merge-partstats.sh regress3c-adv-diffn") == 0 or die "ERR: failed to merge moment files: $!";

# Parse the moments file
my $momentFile;
open($momentFile, "<regress3c-adv-diffnMerged_partstats.csv") or die "ERR: failed to open merged moment file: $!";

my $m0a = 0;
my $m0asq = 0;
my $counta = 0;

my $m2b = 0;
my $m2bsq = 0;
my $countb = 0;

while(<$momentFile>) {
  my @fields = split /,/;

  if(($fields[0] =~ /^2(\.0)?$/) && (abs($fields[1] - 0.017) < 1e-6)) {
      $m0a += $fields[3];
      $m0asq += $fields[3] * $fields[3];
      #print "t=2.0 x=0.017 $fields[3]\n";
      $counta++;
  }

  if(($fields[0] =~ /^3(\.0)?$/) && (abs($fields[1] - 0.0245) < 1e-6)) {
      $m2b += $fields[13];
      $m2bsq += $fields[13] * $fields[13];
      #print "t=3.0 x=0.0245 $fields[3]\n";
      $countb++;
  }
}

# take the mean of the m0 values
$m0a = $m0a / $counta;
$m0asq = sqrt(($m0asq / $counta - $m0a * $m0a) / $counta);
$m2b = $m2b / $countb;
$m2bsq = sqrt(($m2bsq / $countb - $m2b * $m2b) / $countb);

print "$m0a $m0asq\n";
print "$m2b $m2bsq\n";

if(abs($m0a - 198) > 1e1) {
  print "Simulated M0 at t=2.0 near x=0.017 was $m0a, when analytic solution is 198\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 1;
}

if(abs($m2b - 3.03e-47) > 1e-48) {
  print "Simulated M2 at t=3.0 near x=0.0245 was $m2b, when analytic solution is 3.03e-47\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 2;
}


# Clean outputs, there should always be some files to delete.
@outputFiles = glob("regress3c*");
system("rm @outputFiles");


#print "All tests passed\n";
exit 0;
