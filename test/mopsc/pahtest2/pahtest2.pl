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

my $statsFile;
open($statsFile, "<stats.csv") or die "ERR: failed to open stats file: $statsFile";

my $mea = 0;
while(<$statsFile>) {
  my @fields = split /,/;

  # Look for a line that begin with a number and has the first entry (the time)
  # equal (upto a small tolerance) to 0.006
  if($fields[0] =~ /^\d+/) {
      $mea = $fields[0];
      print "$mea\n";
      last;
  }
}

# Precalculated value: 314.18
# 20 repetitions, stats=318.1

if(abs($mea -  318) > 10) {
  print "Stats was $mea, when 314 expected\n";
  print "**************************\n";
  print "****** TEST FAILURE ******\n";
  print "**************************\n";
  exit 3;
}

#print "All tests passed\n";
exit 0;
