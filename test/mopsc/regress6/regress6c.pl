#!/usr/bin/perl

#  Copyright (C) 2009 Robert I A Patterson.
#     modified 2011 William J Menz
#
#   This test checks the weighted particle transition kernel (w3 rule)
#     based on mops.regress2c.
#
# Licence:
#    This file is part of "mopsc".
#
#    mopsc is free software; you can redistribute it and/or
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
my @outputFiles = glob("regression6c-nuc-coag-OH*");
if($#outputFiles > 0) {
  print "Cleaning up old output files\n";
  system("rm", "-f", @outputFiles);
}

# Path of executable should be supplied as first argument to this script
my $program = $ARGV[0];

# Arguments for simulation
my @simulationCommand = ($program, "--flamepp", "-p",
                         "-g", "regress6/regress6.inp",
                        "-c", "regress6/chem.inp",
                        "-t", "regress6/therm.dat",
                         "-s", "regress6/regress6c.xml",
                         "-r", "regress6/regress6c.inx");

# Run the simulation and wait for it to finish
system(@simulationCommand) == 0 or die "ERR: simulation failed: $!";

# Parse the moments file
my $momentFile;
open($momentFile, "<regression6c-nuc-coag-OH-part.csv") or die "ERR: failed to open moment file: $!";

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
# DSA results:
# L = 40, N = 2048: m0 5.51+-0.08e10 g cm^-3, fv 2.23+-0.04e-11
# Weighted particle (w3) results:
# L = 5, N = 512: m0 = 5.46+-0.35 e16 1/cm^3, fv = 2.19+-0.14 e-12 1/cm^6

print "$m0, $m1\n";
if(abs($m0 - 5.46e16) > 3.46e15) {
  print "Simulated mean M0 was $m0, when 5.46e16 m^-3 expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 1;
}

if(abs($m1 - 2.19e-11) > 1.41e-12) {
  print "Simulated mean Fv was $m1, when 2.19e-11 expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 2;
}

# Should always be some files to remove
@outputFiles = glob("regression6c-nuc-coag-OH*");
system("rm", "-f", @outputFiles);
exit 0;
