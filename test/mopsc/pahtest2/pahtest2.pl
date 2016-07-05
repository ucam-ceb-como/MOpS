#!/usr/bin/perl

#  Copyright (C) 2012 Dongping Chen.
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

# this test is designed for checking whether bintree_serializer works correctly with the PAHPrimary Class by comparing to precalculated value

# Parse the moments file
my $momentFile;
open($momentFile, "<pahtest2-bintree-serializer-part.csv") or die "ERR: failed to open moment file: $!";

my $m0 = 0;
my $m1 = 0;

while(<$momentFile>) {
  my @fields = split /,/;

  # Look for a line that begin with a number and has the first entry (the time)
  # equal (upto a small tolerance) to 0.006
  if(($fields[0] =~ /^\d+/) && (abs($fields[1] - 0.006) < 1e-6 )) {
      # Third field should be the zeroth moment
      $m0 = $fields[4];
      #print "4: $fields[4], ";

      $m1 = $fields[16];
      #print "16: $fields[16] \n";

      last;
  }
}

# Precalculated value: M0=2.26e18+-1e16, Fv=7.48e-9+-1e-7

# 20 repetitions
# mean values and 99% confidence interval widths
# m0 (1.41+-0.052)e19 m^-3
# fv (2.57+-0.0556)e-8 

print "$m0, $m1\n";
if(abs($m0 -  1.41e19) > 5.2e17) {
  print "Simulated mean M0 was $m0, when  1.41e19m^-3 expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 1;
}

if(abs($m1 - 2.57e-8) > 5.56e-10) {
  print "Simulated mean Fv was $m1, when 2.57e-8 expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 2;
}

my $postprocessPAHFile;
open($postprocessPAHFile, "<pahtest2-bintree-serializer-postprocess-PAH(0.006s).csv") or die "ERR: failed to open stats file: $postprocessPAHFile";

my $sum = 0;
my $counter = 1;
my $average = 0;

# Ignore lines which contain non-numeric values: the header (line 1) and the three lines containing the words 1runs, 2runs and 3runs (lines 2, 747 and 1489)
while (my $line = <$postprocessPAHFile>) {
  if (($counter == 1) || ($counter == 2) || ($counter == 747) || ($counter == 1489)) {
  } else {
      chomp $line;

      my @fields = split "," , $line;
      $sum += $fields[6];
  }
  $counter += 1;
}

# ($counter - 5) corresponds to the total number of PAHs
$average = $sum / ($counter - 5);

# Precalculated average over 3 runs = 455.6522316
# 20 runs: mean = 455.2484913; 99% confidence level = 4.826335734

if(abs($average -  455) > 5) {
  print "Average PAH mass (u) was $average, when 455 expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 3;
}

#print "All tests passed\n";
exit 0;
