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

# See if this is a windows system
my $windows = ($ENV{'OS'} =~ /windows.*/i);

# Choose the windows executable name if appropriate
my $program = "../../bin/debug/brush";
if($windows) {
    $program = "../../bin/debug/brush.exe";
}

# Arguments for simulation
my @simulationCommand = ($program,
                         "-g", "regress4/geometry4.xml",
                         "-c", "regress4/chem.inp",
                         "-d", "regress4/chemsoln4.dat",
                         "-s", "regress4/sweep4a.xml",
                         "-t", "regress4/therm.dat",
                         "-b", "regress4/brush4a.xml",
                         "-a", "regress4/partsoln4.xml",
                         "-e", "1066");

# Run the simulation and wait for it to finish
system(@simulationCommand) == 0 or die "ERR: simulation failed: $!";

# Parse the moments file
my $momentFile;
open($momentFile, "<regress4a1189_partstats.csv") or die "ERR: failed to open moment file: $!";

my @m1samples;

# Read the line of headings
<$momentFile>;

# Now read the numeric data
while(<$momentFile>) {
  my @fields = split /,/;

  # only accept samples collected from t=0.12 onwards and at cell with centre 
  if(($fields[0] > 0.1199999) && (abs($fields[1] - 0.079) < 1e-6)) {
      push @m1samples, $fields[11];
  }
}

if (46 != @m1samples) {
  my $numSamples = @m1samples;
  die "Found $numSamples, but there should have been 46";
}

# calculate the mean and variance of the samples
my $sum = 0;
my $sumsq = 0;

my $val;
foreach $val (@m1samples) {
  $sum += $val;
  $sumsq += $val * $val;
}

my $mean = $sum / 46;
my $var  = $sumsq / 46 - $mean * $mean;

print "$mean, $var\n";
if(abs($mean - 8.374e-21) > 2e-22) {
  print "Simulated mean M1 at x=0.079 was $mean, when analytic solution is 8.374e-21 kg m^-3\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 1;
}

if(abs($var - 2.219e-43) > 3e-44) {
  print "Simulated M1 variance at x=0.079 was $var, when analytic solution is 2.219e-43 kg^2 m^-6\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 2;
}

system("rm regress4a*.csv");

#print "All tests passed\n";
exit 0;
